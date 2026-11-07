"""Real solved DP target: source/target prefix learning without new IBP points."""
import sys,json,time,os
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'tests'))
from reference_space import BootstrapSpec,RatFunc,RatMat,build_pl_ansatz
from analytic_bootstrap.learning import learn_profiles
from analytic_bootstrap.learning_compare import FactorSpaceReference
from analytic_bootstrap.pl import plan_pl_orders
D=json.loads((ROOT/'results/dp_complete_run/connection.json').read_text());R=json.loads((ROOT/'results/dp_reduction_run/reduction.json').read_text());p=D['prime'];known=set(D['known_rows'])
choices=sorted((sum(len(v['numerator_tau']) for v in row.values()),n) for n,row in R['reductions'].items() if n.startswith('J') and n not in D['masters'] and set(row)<=known and 100<sum(len(v['numerator_tau']) for v in row.values())<200)
_,name=choices[0];truth=R['reductions'][name];deps=set(truth);front=list(deps)
while front:
 for j in D['de'][front.pop()]['rhs']:
  if j not in deps:deps.add(j);front.append(j)
sub=[m for m in D['masters'] if m in deps];assert set(sub)<=known
rf=lambda r:RatFunc.make(r['numerator_tau'],r['denominator_tau'],p)
zero=RatFunc.make([0],[1],p);B={(i,j):rf(D['de'][m]['rhs'][n]) for i,m in enumerate(sub) for j,n in enumerate(sub) if n in D['de'][m]['rhs']}
ep=[200001*200000%p,400001,1];shift=RatFunc.make([400001*3%p,6],ep,p)
for i in range(len(sub)):B[i,i]=B.get((i,i),zero)+shift
spec=BootstrapSpec.from_inputs(family={'source':'real-solved-DP-target'},variable='tau',assignments={},masters=sub,target=[name],de=RatMat(B,(len(sub),len(sub)),p),denominators=[[truth[j]['denominator_tau'] if j in truth else [1] for j in sub]],numerator_caps=[[len(truth[j]['numerator_tau'])-1 if j in truth else -1 for j in sub]],prime=p,mode='reduction')
t=time.perf_counter();plain=build_pl_ansatz(spec,backend='PL-Hasse',laurent_mode='structural-refined-all');plain_time=time.perf_counter()-t
t=time.perf_counter();profiles=learn_profiles(spec);learning=time.perf_counter()-t
t=time.perf_counter();learned=build_pl_ansatz(spec,backend='PL-Hasse',laurent_mode='structural-refined-all',order_profile=profiles['order'],target_prefix_profile=profiles['prefix']);production=time.perf_counter()-t
assert plain.dimension==learned.dimension
assert all(FactorSpaceReference(a,p).equal(b) for a,b in zip(plain.blocks,learned.blocks,strict=True))
reference=plan_pl_orders(spec);prefix=plan_pl_orders(spec,order_profile=profiles['order'],target_prefix_profile=profiles['prefix'])
report=dict(target=name,sub_rank=len(sub),dimension=plain.dimension,plain_seconds=plain_time,learning_seconds=learning,production_seconds=production,exact_factor_spaces_equal=True,reference_orders=[x['order'] for x in reference['factors']],prefix_orders=[x['order'] for x in prefix['factors']],profile_timing=profiles['timing'],cache=profiles['cache'])
(ROOT/'results/dp_official/reduction_prefix.json').write_text(json.dumps(report,indent=2));print(json.dumps(report),flush=True)
