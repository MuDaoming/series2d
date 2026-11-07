exec(open(__file__.replace('new_fixed.py','check_new_coefficients.py')).read().split('seen={};bad=[]')[0])
from reference_space import BootstrapSpec,RatMat,build_pl_ansatz
from native_fit import joint_fit
import time
rf=lambda v:RatFunc.make(v['numerator_tau'],v['denominator_tau'],p)
fixed={(i,j):rf(v) for i,j,v in c['fixed']}
for name in ['shared_de.json','shared_de_pending_v2.json','preconditioned_de.json','preconditioned_de_all.json','supplemental_v1.json']:
 path=ROOT/'results/dp_pole_route'/name
 if not path.exists():continue
 for m,row in json.loads(path.read_text()).get('de',{}).items():
  if m not in rp:continue
  for j,v in row['rhs'].items():
   if j in cp:fixed[rp[m],cp[j]]=rf(v)
for (i,j),r in fixed.items():c['denominators'][i][j]=r.den;c['numerator_caps'][i][j]=len(r.num)-1
spec=BootstrapSpec.from_inputs(family={'source':'new-known-C-diagnosis'},variable='tau',assignments={},masters=c['masters'],target=None,de={'A':RatMat({(i,j):rf(v) for i,j,v in c['A']},tuple(c['A_shape']),p),'B':RatMat({(i,j):rf(v) for i,j,v in c['B']},tuple(c['B_shape']),p)},denominators=c['denominators'],numerator_caps=c['numerator_caps'],prime=p,mode='de')
t=time.perf_counter();a=build_pl_ansatz(spec,backend='PL-Hasse',laurent_mode='structural-refined-all');print('BUILT',a.dimension,time.perf_counter()-t,flush=True)
fit=joint_fit(a,[],fixed,c['identity_points']);report={'known_functions':len(fixed),'dimension':a.dimension,'fit':{k:v for k,v in fit.items() if k!='parameters'}};(out/'new_fixed.json').write_text(json.dumps(report,indent=2));print(json.dumps(report),flush=True)
