"""A low-degree numerator with a larger denominator must not need symmetric bounds."""
import json
from math import comb
from pathlib import Path
import subprocess
import tempfile
import unittest

ROOT=Path(__file__).resolve().parents[1]


class AsymmetricReconstruction(unittest.TestCase):
    def test_large_denominator_and_corrupted_holdout(self):
        p=1000003
        with tempfile.TemporaryDirectory() as name:
            work=Path(name);family=work/'family.json';family.write_text(json.dumps(dict(prime=p,t0='0')))
            points=[]
            for x in range(17,42):
                value=(5+3*x)*pow(pow(1+x,10,p),-1,p)%p
                path=work/f'point{x}.json';path.write_text(json.dumps(dict(prime=p,t_point=str(x),DE_values={'M':{'M':value}})))
                points.append(str(path))
            manifest=work/'manifest.json';manifest.write_text(json.dumps(dict(family=str(family),points=points)))
            output=work/'answer.json';binary=ROOT/'build/reconstruct_partial'
            run=subprocess.run([str(binary),str(manifest),str(output)],capture_output=True,text=True)
            self.assertEqual(run.returncode,0,run.stderr)
            answer=json.loads(output.read_text())['de']['M']['rhs']['M']
            self.assertEqual(answer['numerator_tau'],[5,3])
            self.assertEqual(answer['denominator_tau'],[comb(10,k) for k in range(11)])
            bad=json.loads(Path(points[-1]).read_text());bad['DE_values']['M']['M']=(bad['DE_values']['M']['M']+1)%p
            Path(points[-1]).write_text(json.dumps(bad))
            run=subprocess.run([str(binary),str(manifest),str(output)],capture_output=True,text=True)
            self.assertEqual(run.returncode,2)
            self.assertEqual(json.loads(output.read_text())['unreconstructed'],[['M','M']])


if __name__=='__main__':unittest.main()
