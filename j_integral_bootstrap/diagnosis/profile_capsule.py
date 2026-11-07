"""Replay the frozen 6x66 candidate, without reading live checkpoints or blind data."""
exec(open(__file__.replace('profile_capsule.py','normalize_capsule.py')).read().split('base=spec(A,B);')[0])
from native_fit import joint_fit
records=[]
fixed={(i,j):rf(v) for i,j,v in c['fixed']}
for extra in [0,4,12]:
 start=time.perf_counter();a=build_pl_ansatz(spec(A,B),backend='PL-Hasse',laurent_mode='structural-refined-all',extra_order=extra)
 construction=time.perf_counter()-start;fit=joint_fit(a,c['samples'],fixed,c['identity_points'])
 rec={'extra_order':extra,'dimension':a.dimension,'construction_seconds':construction,'fit':{k:v for k,v in fit.items() if k!='parameters'},'diagnostics':a.diagnostics,'factors':[{'factor':b.factor,'metadata':b.metadata} for b in a.blocks]};records.append(rec)
 (OUT/'profile_capsule.json').write_text(json.dumps(records,indent=2,default=str));print(extra,a.dimension,fit['compatible'],flush=True)
