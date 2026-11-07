"""The DE verifier requires a real completed run and nonempty comparison set."""
import json
from pathlib import Path
import subprocess
import tempfile
import unittest

ROOT=Path(__file__).resolve().parents[1]


class VerifyConnection(unittest.TestCase):
    def test_coefficients_zeros_and_certificate_inputs(self):
        with tempfile.TemporaryDirectory() as td:
            d=Path(td);m0,m1='J{1,0}','J{0,1}'
            answer=dict(status='complete',prime=101,t0='2/5',masters=[m0,m1],
                de={m0:dict(rhs={m0:dict(numerator_tau=[3,1],denominator_tau=[1])}),m1:dict(rhs={})})
            point=dict(prime=101,t_point=7,masters_at_point=[m0,m1],unclosed_derivatives=[],
                       DE_values={m0:{m0:50},m1:{}})
            manifest=dict(points=[str(d/'point.json')])
            def run():
                for name,value in [('answer',answer),('point',point),('manifest',manifest)]:
                    (d/(name+'.json')).write_text(json.dumps(value))
                return subprocess.run([str(ROOT/'build/verify_connection'),str(d/'answer.json'),
                    str(d/'manifest.json'),str(d/'check.json')],capture_output=True,text=True)
            self.assertEqual(run().returncode,0)
            self.assertEqual(json.loads((d/'check.json').read_text())['scalar_checks_including_zeros'],4)
            point['DE_values'][m1][m0]=1
            self.assertIn('coefficient mismatch',run().stderr)
            point['DE_values'][m1]={}
            del answer['status']
            self.assertNotEqual(run().returncode,0)
            answer['status']='complete';manifest['points']=[]
            self.assertIn('nonempty point manifest',run().stderr)
            manifest['points']=None
            self.assertIn('nonempty point manifest',run().stderr)


if __name__=='__main__':unittest.main()
