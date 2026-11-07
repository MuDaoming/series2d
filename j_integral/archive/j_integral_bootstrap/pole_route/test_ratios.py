"""A high-degree shared numerator must be cancellable without high-order sampling."""
import json,subprocess,tempfile,sys
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'tests'))
from rational_learning import make_learner
p=1000003;center=200001;L=make_learner(p,center)
N=L.one;D=L.one
for i in range(50):N*=L.poly([300+i,1])
for i in range(45):D*=L.poly([600+i,1])
H=N*L.poly([7,3]);E=D*L.poly([13,1])
with tempfile.TemporaryDirectory() as td:
 d=Path(td);known=dict(prime=p,t0='2/5',reductions={'R':{'anchor':dict(numerator_tau=[int(v) for v in N],denominator_tau=[int(v) for v in D])}})
 (d/'known.json').write_text(json.dumps(known));points=[]
 for t in range(17,29):
  x=(t-center)%p;point=dict(prime=p,t_point=str(t),masters_at_point=['anchor','target'],reduction_values={'R':{'anchor':int(N(x))*pow(int(D(x)),-1,p)%p,'target':int(H(x))*pow(int(E(x)),-1,p)%p}})
  f=d/f'point{t}.json';f.write_text(json.dumps(point));points.append(str(f))
 m=dict(prime=p,t0='2/5',masters=['anchor','target'],holdout=3,section='reduction_values',points=points,known=[str(d/'known.json')],ratios=True)
 def run():
  (d/'manifest.json').write_text(json.dumps(m));subprocess.run([str(ROOT/'build/preconditioned'),str(d/'manifest.json'),str(d/'out.json')],check=True,capture_output=True);return json.loads((d/'out.json').read_text())
 out=run()['reductions']['R']['target'];assert not L.poly(out['numerator_tau'])*E-H*L.poly(out['denominator_tau'])
 v=json.loads(Path(points[-1]).read_text());v['reduction_values']['R']['target']=(v['reduction_values']['R']['target']+1)%p;Path(points[-1]).write_text(json.dumps(v));assert 'target' not in run()['reductions'].get('R',{})
 print('51-degree numerator learned from12 points by exact ratio cancellation; corrupted holdout rejected')
