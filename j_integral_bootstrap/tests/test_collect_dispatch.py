"""Exercise the real collector: short jobs release slots without a batch wait."""
import json,subprocess,sys,tempfile,unittest
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
class Dispatcher(unittest.TestCase):
 def test_dynamic_dispatch_and_resume(self):
  with tempfile.TemporaryDirectory() as td:
   d=Path(td)
   for name,data in [('family',dict(prime=101)),('basis',['J']),('targets',[[1]])]:(d/(name+'.json')).write_text(json.dumps(data))
   binary=d/'mock.py';binary.write_text('''#!/usr/bin/env python3
import json,sys,time
from pathlib import Path
point=int(sys.argv[sys.argv.index('--point')+1]);start=time.time();time.sleep(.4 if point==1 else .02)
Path(sys.argv[2]).write_text(json.dumps(dict(prime=101,t_point=str(point),masters_at_point=['J'],unclosed_derivatives=[],DE_values={'J':{}},reduction_values={'J':{}},mock_start=start,mock_end=time.time())))
''');binary.chmod(0o755)
   cmd=[sys.executable,str(ROOT/'tests/collect_points.py'),'full',str(d/'work'),'--family',str(d/'family.json'),'--basis',str(d/'basis.json'),'--targets',str(d/'targets.json'),'--binary',str(binary),'--first','1','--count','3','--workers','2']
   subprocess.run(cmd,check=True,capture_output=True)
   one=json.loads((d/'work/points/point1.json').read_text());three=json.loads((d/'work/points/point3.json').read_text())
   self.assertLess(three['mock_start'],one['mock_end'])
   before={f.name:f.read_bytes() for f in (d/'work/points').glob('*.json')}
   subprocess.run(cmd,check=True,capture_output=True)
   self.assertEqual(before,{f.name:f.read_bytes() for f in (d/'work/points').glob('*.json')})
   self.assertEqual(len(list((d/'work/receipts').glob('point*.json'))),3)
if __name__=='__main__':unittest.main()
