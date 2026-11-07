"""Exact supported-factor space comparison for classic and complete runtimes."""
import sys,os,json,time
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parent))
import test_preflight as regression
from analytic_bootstrap.pl import build_pl_ansatz
from analytic_bootstrap.learning_compare import FactorSpaceReference
from analytic_bootstrap.backends.pl.memory_transport import close_all,observe_calls
from preflight import preflight
ROOT=regression.ROOT;OUT=ROOT/'results/dp_official'
saved=[];original=regression.preflight

def capture(ansatz,samples,**kw):
    if not saved or ansatz is not saved[-1][0]:saved.append((ansatz,samples))
    return original(ansatz,samples,**kw)
regression.preflight=capture
start=time.perf_counter();regression.main();classic_total=time.perf_counter()-start
reports=[]
for label,(old,samples) in zip(['failed-six-by-66','solved-three-row-sector'],saved):
    close_all();os.environ['ANALYTIC_BOOTSTRAP_KERNEL_DIR']=str(ROOT/'build/production')
    start=time.perf_counter()
    with observe_calls() as calls:new=build_pl_ansatz(old.spec,backend='PL-Hasse',laurent_mode='structural-refined-all')
    seconds=time.perf_counter()-start
    record=dict(label=label,old_dimension=old.dimension,new_dimension=new.dimension,construction_seconds=seconds,worker_calls=calls)
    assert len(old.blocks)==len(new.blocks)
    start=time.perf_counter()
    equal=[FactorSpaceReference(a,old.spec.prime).equal(b) for a,b in zip(old.blocks,new.blocks,strict=True)]
    assert all(equal),[i for i,x in enumerate(equal) if not x]
    record.update(exact_factor_spaces_equal=all(equal),factors=len(equal),space_comparison_seconds=time.perf_counter()-start)
    record['preflight']=preflight(new,samples,max_points=len(samples),kernel_directory=ROOT/'build/upstream')
    print(json.dumps(record),flush=True);reports.append(record)
close_all()
(OUT/'runtime_comparison.json').write_text(json.dumps(dict(classic_regression_total_seconds=classic_total,reports=reports),indent=2))
