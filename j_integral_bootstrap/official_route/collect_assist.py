"""Collect restricted R-only points, with atomic outputs and explicit inputs."""
import argparse,json,subprocess,concurrent.futures,time
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
p=argparse.ArgumentParser();p.add_argument('work');p.add_argument('--stop',type=int,default=225);p.add_argument('--workers',type=int,default=2);a=p.parse_args();O=Path(a.work).resolve();start=time.monotonic()
def one(t):
 out=O/'points'/f'point{t}.json'
 if out.exists():return t,True
 tmp=out.with_suffix('.tmp.json')
 cmd=[str(ROOT/'build/global_point_subset'),str(O/'family.json'),str(tmp),'--point',str(t),'--targets-file',str(O/'targets.json'),'--basis-file',str(O/'basis.json'),'--derivative-names-file',str(O/'no_derivatives.json'),'--integration','atomic','--values','1','--compact','1']
 with (O/f'point{t}.log').open('w') as log:subprocess.run(cmd,check=True,stdout=log,stderr=log)
 d=json.loads(tmp.read_text());assert not d['unclosed_derivatives'];assert d['masters_at_point']==json.loads((O/'basis.json').read_text());tmp.replace(out);return t,False
points=[t for t in range(17,a.stop) if not 151<=t<=158]
with concurrent.futures.ThreadPoolExecutor(max_workers=a.workers) as pool:
 for i,(t,cached) in enumerate(pool.map(one,points)):
  if i%20==0:print(json.dumps(dict(completed=i+1,t=t,cached=cached,seconds=time.monotonic()-start)),flush=True)
