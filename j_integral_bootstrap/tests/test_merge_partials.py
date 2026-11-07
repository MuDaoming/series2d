import json,subprocess,sys,tempfile,unittest
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
class Merge(unittest.TestCase):
 def test_conflict_field_and_implicit_zero(self):
  with tempfile.TemporaryDirectory() as td:
   d=Path(td);(d/'points').mkdir()
   files={'family':dict(prime=101,t0='2/5'),'basis':['M0','M1'],'a':dict(prime=101,t0='2/5',de={'M0':{'rhs':{'M1':dict(numerator_tau=[7],denominator_tau=[1])}}})}
   for name,v in files.items():(d/(name+'.json')).write_text(json.dumps(v))
   (d/'points/point7.json').write_text(json.dumps(dict(prime=101,t_point='7',masters_at_point=files['basis'],unclosed_derivatives=[],DE_values={'M0':{'M1':7},'M1':{}})))
   cmd=[sys.executable,str(ROOT/'tests/merge_partials.py'),'--family',str(d/'family.json'),'--basis',str(d/'basis.json'),'--input',str(d/'a.json'),'--points',str(d/'points'),'--section','de','--output',str(d/'out.json')]
   self.assertEqual(subprocess.run(cmd,capture_output=True).returncode,0)
   self.assertEqual(json.loads((d/'out.json').read_text())['status'],'complete')
   bad=json.loads(json.dumps(files['a']));bad['de']['M0']['rhs']['M1']['numerator_tau']=[8];(d/'bad.json').write_text(json.dumps(bad))
   self.assertNotEqual(subprocess.run(cmd+['--input',str(d/'bad.json')],capture_output=True).returncode,0)
   bad['prime']=103;(d/'bad.json').write_text(json.dumps(bad))
   self.assertNotEqual(subprocess.run(cmd+['--input',str(d/'bad.json')],capture_output=True).returncode,0)
   bad=files['a'];bad['de']['M1']={'rhs':{'M0':dict(numerator_tau=[1],denominator_tau=[1])}};(d/'a.json').write_text(json.dumps(bad))
   self.assertNotEqual(subprocess.run(cmd,capture_output=True).returncode,0)
if __name__=='__main__':unittest.main()
