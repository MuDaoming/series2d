import sys,json,os,time,argparse
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'tests'))
from reference_space import BootstrapSpec,RatFunc,RatMat,build_pl_ansatz
from analytic_bootstrap.upstream.fitting import ProductionFitter
from analytic_bootstrap.sampling import InconsistentSampleError
from analytic_bootstrap.poly import poly_mul
os.environ['ANALYTIC_BOOTSTRAP_UPSTREAM_KERNEL_DIR']=str(ROOT/'results/dp_official/upstream_native')
os.environ['ANALYTIC_BOOTSTRAP_FIT_THREADS']='1'
OUT=ROOT/'results/dp_diagnosis'
c=json.loads((OUT/'capsule.json').read_text());p=c['prime'];rf=lambda v:RatFunc.make(v['numerator_tau'],v['denominator_tau'],p)
A={(i,j):rf(v) for i,j,v in c['A']};B={(i,j):rf(v) for i,j,v in c['B']}
spec=lambda a,b:BootstrapSpec.from_inputs(family={'source':'diagnosis-capsule'},variable='tau',assignments={},masters=c['masters'],target=None,de={'A':RatMat(a,tuple(c['A_shape']),p),'B':RatMat(b,tuple(c['B_shape']),p)},denominators=c['denominators'],numerator_caps=c['numerator_caps'],prime=p,mode='de')
base=spec(A,B);factors=[list(fs.factor) for fs in base.envelope.factors];records=[]
for label,selected,exponent in [('endpoint+3',[f for f in factors if f in [[200000,1],[200001,1]]],3),('endpoint+9',[f for f in factors if f in [[200000,1],[200001,1]]],9),('other+3',[f for f in factors if f not in [[200000,1],[200001,1]]],3)]:
 h=[1]
 for f in selected:
  for _ in range(exponent):h=poly_mul(h,f,p)
 shift=RatFunc.make([(i*h[i])%p for i in range(1,len(h))],h,p);a=dict(A);b=dict(B);zero=RatFunc.make([0],[1],p)
 for i in range(c['A_shape'][0]):a[i,i]=a.get((i,i),zero)+shift
 for i in range(c['B_shape'][0]):b[i,i]=b.get((i,i),zero)+shift
 t=time.perf_counter();ansatz=build_pl_ansatz(spec(a,b),backend='PL-Hasse',laurent_mode='structural-refined-all');construct=time.perf_counter()-t
 fitter=ProductionFitter(ansatz,backend='stream');t=time.perf_counter();status='compatible';res={}
 try:
  for s in c['samples']:
   res=fitter.add(s)
   if res['parameters'] is not None:break
 except InconsistentSampleError:status='inconsistent'
 finally:fitter.close()
 rec={'label':label,'dimension':ansatz.dimension,'construction_seconds':construct,'fit_seconds':time.perf_counter()-t,'status':status,'rank':res.get('rank'),'points':len(res.get('sample_points',[]))};records.append(rec);(OUT/'normalization.json').write_text(json.dumps(records,indent=2));print(json.dumps(rec),flush=True)
