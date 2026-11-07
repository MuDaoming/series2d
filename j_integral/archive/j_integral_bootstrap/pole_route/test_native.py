"""Meaningful regression for the native preconditioned rational learner."""
import json,sys,tempfile,subprocess
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'tests'));sys.path.insert(0,str(ROOT/'pole_route'))
from rational_learning import make_learner
from preconditioned_fit import make_preconditioned_fit
P=1000003;center=200001;L=make_learner(P,center);fit=make_preconditioned_fit(L,P,center,holdout=3)
def den(n):
 d=L.one
 for i in range(n):d*=L.poly([200+i,1])
 return d
cases=[('low_numerator_high_denominator',L.poly([7,3,5]),den(50),L.one,False),('partial_poles_missing_E',L.poly([6,4,3,2]),den(24),L.poly([17,8,1]),False),('corrupt_check_point',L.poly([6,4,3,2]),den(24),L.poly([17,8,1]),True)]
# Add an actual reconstructed DP coefficient, independently checked against the point files below.
real=ROOT/'results/dp_pole_route/native_de43.json'
if real.exists():
 d=json.load(open(real));added=False
 for row,vs in d['de'].items():
  for col,c in vs['rhs'].items():
   if len(c['numerator_tau'])<38:
    cases.append(('real_dp_entry',L.poly(c['numerator_tau']),L.poly(c['denominator_tau']),L.one,False));added=True;break
  if added:break
with tempfile.TemporaryDirectory(prefix='pole_regression_') as temp:
 t=Path(temp);known={'prime':P,'t0':'2/5','reductions':{}};ps=[];expected={}
 for name,H,D,E,bad in cases:known['reductions'][name]={'anchor':{'numerator_tau':[1],'denominator_tau':[int(c) for c in D]}}
 (t/'known.json').write_text(json.dumps(known))
 for k in range(41):
  time=17+k;x=(time-center)%P;v={}
  for name,H,D,E,bad in cases:
   value=int(H(x))*pow(int((D*E)(x)),-1,P)%P
   if bad and k==40:value=(value+1)%P
   v[name]={'target':value,'anchor':pow(int(D(x)),-1,P)}
  p=t/f'point{k}.json';p.write_text(json.dumps({'prime':P,'masters_at_point':['anchor','target'],'t_point':str(time),'reduction_values':v}));ps.append(str(p))
 manifest={'prime':P,'t0':'2/5','masters':['anchor','target'],'holdout':3,'section':'reduction_values','points':ps,'known':[str(t/'known.json')]};(t/'manifest.json').write_text(json.dumps(manifest))
 subprocess.run([str(ROOT/'build/preconditioned'),str(t/'manifest.json'),str(t/'out.json')],check=True,capture_output=True)
 out=json.load(open(t/'out.json'))['reductions']
 for name,H,D,E,bad in cases:
  values=[json.load(open(p))['reduction_values'][name]['target'] for p in ps]
  py=fit(D,list(range(17,58)),values);cpp=out.get(name,{}).get('target')
  if bad:assert py is None and cpp is None,(name,py,cpp)
  else:
   assert py is not None and cpp is not None,name
   assert list(py.num)==cpp['numerator_tau'] and list(py.den)==cpp['denominator_tau'],name
   expected=L.poly(cpp['numerator_tau'])*(D*E)-H*L.poly(cpp['denominator_tau']);assert not expected,name
  print(name,'PASS')
 # Input identities and overlapping observations are part of the reconstruction contract.
 for kind in ['wrong_prime','wrong_basis','duplicate_numeric_conflict','duplicate_rational_conflict','known_tau_mismatch','known_prime_mismatch']:
  m=json.loads(json.dumps(manifest))
  if kind=='wrong_prime':m['prime']=1000033
  elif kind=='wrong_basis':m['masters']=['target','anchor']
  elif kind in ['known_tau_mismatch','known_prime_mismatch']:
   d=json.loads(json.dumps(known));d['t0']='3/5' if kind=='known_tau_mismatch' else '2/5';d['prime']=1000033 if kind=='known_prime_mismatch' else P
   f=t/'mismatched_known.json';f.write_text(json.dumps(d));m['known']=[str(f)]
  elif kind=='duplicate_numeric_conflict':
   d=json.load(open(ps[0]));d['reduction_values'][cases[0][0]]['target']=(d['reduction_values'][cases[0][0]]['target']+1)%P
   f=t/'conflicting_point.json';f.write_text(json.dumps(d));m['points']=ps+[str(f)]
  else:
   d=json.loads(json.dumps(known));d['reductions'][cases[0][0]]['anchor']['numerator_tau']=[2]
   f=t/'conflicting_known.json';f.write_text(json.dumps(d));m['known']=manifest['known']+[str(f)]
  (t/'bad_manifest.json').write_text(json.dumps(m))
  result=subprocess.run([str(ROOT/'build/preconditioned'),str(t/'bad_manifest.json'),str(t/'bad_out.json')],capture_output=True)
  assert result.returncode!=0,kind
  print(kind,'PASS')
print('All tests use three internal checks; no final blind verification points are read.')
