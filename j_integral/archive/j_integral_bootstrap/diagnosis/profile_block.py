import sys,json,shutil,argparse,time
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'tests'))
import bootstrap_global as bg
OUT=ROOT/'results/dp_diagnosis'; WORK=OUT/'snapshot';WORK.mkdir(exist_ok=True)
shutil.copy2(ROOT/'results/dp_complete_run/connection.json',WORK/'connection.json')
original=bg.build_pl_ansatz; fit=bg.joint_fit; state={}
class Done(BaseException):pass
def build(spec,**kw):
 state['spec']=spec;state['kw']=kw
 return original(spec,**kw)
def fitting(ansatz,samples,fixed,identity_points):
 records=[]
 for extra in [0,4,12]:
  t=time.perf_counter()
  a=ansatz if extra==0 else original(state['spec'],**state['kw'],extra_order=extra)
  construction=time.perf_counter()-t if extra else a.diagnostics['construction_seconds']
  r=fit(a,samples,fixed,identity_points)
  rec={'extra_order':extra,'construction_seconds':construction,'dimension':a.dimension,'ambient':a.ambient_parameters,'fit':{k:v for k,v in r.items() if k!='parameters'},'diagnostics':a.diagnostics,'factors':[{'factor':b.factor,'metadata':b.metadata} for b in a.blocks]}
  records.append(rec);(OUT/'profile_block.json').write_text(json.dumps(records,indent=2,default=str));print('PROFILE',extra,rec['dimension'],r['compatible'],construction,flush=True)
 raise Done()
bg.build_pl_ansatz=build;bg.joint_fit=fitting
args=argparse.Namespace(work=str(WORK),family=str(ROOT/'results/dp_full/input/family.json'),basis=str(ROOT/'results/dp_full/input/basis.json'),points=[str(ROOT/'results'/x/'points') for x in ['dp_full','dp_extra_de','dp_exceptional']],diagonals=[str(ROOT/'results/dp_cuts/reconstruction.json'),str(ROOT/'results/dp_exceptional/diagonal33.json')],partial=[str(ROOT/'results/dp_driven/partial_33.json'),str(ROOT/'results/dp_exceptional/partial33.json'),str(ROOT/'results/dp_extra_de/partial41.json')],backend='PL-structural',witness_budget=0,internal_holdouts=3,relax_all_poles=False,restart=False)
try:bg.main(args)
except Done:pass
