"""Native verifier must check identities, zeros, and corrupted coefficients."""
import json
from pathlib import Path
import subprocess
import tempfile
import unittest
ROOT=Path(__file__).resolve().parents[1]
class VerifyReductions(unittest.TestCase):
 def test_complete_and_corrupt(self):
  with tempfile.TemporaryDirectory() as td:
   d=Path(td);m='J{1}';target='J{2}';zero='Bz0{2}'
   answer=dict(status='complete',prime=101,t0='2/5',masters=[m],tested_objects=[m,target,zero],reductions={target:{m:dict(numerator_tau=[3,1],denominator_tau=[1])},zero:{}})
   point=dict(prime=101,t_point='7',masters_at_point=[m],unclosed_derivatives=[],reduction_values={m:{m:1},target:{m:50},zero:{}})
   # 2/5 = 61 mod 101; tau=47, tau+3=50.
   (d/'point.json').write_text(json.dumps(point));(d/'manifest.json').write_text(json.dumps(dict(points=[str(d/'point.json')])))
   def run():
    (d/'answer.json').write_text(json.dumps(answer))
    return subprocess.run([str(ROOT/'build/verify_reductions'),str(d/'answer.json'),str(d/'manifest.json'),str(d/'check.json')],capture_output=True,text=True)
   self.assertEqual(run().returncode,0)
   point['reduction_values'][zero]={m:1};(d/'point.json').write_text(json.dumps(point))
   self.assertNotEqual(run().returncode,0)
   point['reduction_values'][zero]={};(d/'point.json').write_text(json.dumps(point))
   del answer['reductions'][zero]
   self.assertIn('missing target',run().stderr)
 def test_rejects_empty_inputs_and_missing_status(self):
  with tempfile.TemporaryDirectory() as td:
   d=Path(td);m='J{1}'
   answer=dict(status='complete',prime=101,t0='0',masters=[m],tested_objects=[m],reductions={})
   point=dict(prime=101,t_point=2,masters_at_point=[m],unclosed_derivatives=[],reduction_values={m:{m:1}})
   manifest=dict(points=[str(d/'point.json')])
   def run():
    for name,value in [('answer',answer),('point',point),('manifest',manifest)]:
     (d/(name+'.json')).write_text(json.dumps(value))
    return subprocess.run([str(ROOT/'build/verify_reductions'),str(d/'answer.json'),str(d/'manifest.json'),str(d/'check.json')],capture_output=True,text=True)
   self.assertEqual(run().returncode,0)
   del answer['status'];self.assertNotEqual(run().returncode,0)
   answer['status']='complete';manifest['points']=[]
   self.assertIn('nonempty point manifest',run().stderr)
   manifest['points']=[str(d/'point.json')];answer['tested_objects']=[]
   self.assertIn('nonempty tested objects',run().stderr)
if __name__=='__main__':unittest.main()
