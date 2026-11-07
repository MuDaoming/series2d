"""Reconstruct remaining numerical J reductions, require coverage, verify fresh data."""
import argparse,json,subprocess,sys,os
from pathlib import Path
HERE=Path(__file__).resolve().parent;ROOT=HERE.parent
sys.path.insert(0,str(ROOT/'tests'))
from reference_space import RatFunc
p=argparse.ArgumentParser();p.add_argument('--work',required=True);p.add_argument('--base',required=True);p.add_argument('--targets',required=True);p.add_argument('--verification',required=True);p.add_argument('--witness-only');a=p.parse_args()
O=Path(a.work).resolve();O.mkdir(exist_ok=True,parents=True);base=json.loads(Path(a.base).read_text());prime=base['prime'];masters=base['masters'];points=sorted((O/'samples/points').glob('point*.json'),key=lambda p:int(p.stem[5:]));assert len(points)==400,'sampling coverage is not400'
manifest=dict(prime=prime,t0=base['t0'],masters=masters,section='reduction_values',holdout=8,ratios=True,ratio_budget=12,max_denominator_degree=512,points=[str(p) for p in points],known=[str(Path(a.base).resolve())])
if a.witness_only:manifest['witness_only']=[dict(path=str(Path(a.witness_only).resolve()),section='de')]
(O/'reconstruct_manifest.json').write_text(json.dumps(manifest))
subprocess.run([str(ROOT/'build/preconditioned'),str(O/'reconstruct_manifest.json'),str(O/'reconstructed.json')],check=True)
new=json.loads((O/'reconstructed.json').read_text());known=base['reductions'];zero=dict(numerator_tau=[0],denominator_tau=[1])
def rf(v):return RatFunc.make(v['numerator_tau'],v['denominator_tau'],prime)
for row,vs in new['reductions'].items():
 dest=known.setdefault(row,{})
 for col,v in vs.items():
  if col in dest and not (rf(dest[col])-rf(v)).is_zero():raise ValueError('rational overlap differs')
  dest[col]=v
support={}
for path in points:
 d=json.loads(path.read_text());assert d['prime']==prime and d['masters_at_point']==masters and not d['unclosed_derivatives']
 for row,vs in d['reduction_values'].items():support.setdefault(row,set()).update(c for c,v in vs.items() if v)
complete=set(base['known_targets']);missing={}
for row,cols in support.items():
 absent=cols-set(known.get(row,{}))
 if absent:missing[row]=sorted(absent)
 else:complete.add(row)
targets=json.loads(Path(a.targets).read_text());remaining=[name for name in targets if name not in complete and name not in masters]
answer=dict(status='partial' if remaining else 'complete',prime=prime,t0=base['t0'],masters=masters,known_targets=sorted(complete),tested_objects=targets,reductions=known,pending_targets=remaining)
(O/'reduction.json').write_text(json.dumps(answer,separators=(',',':')))
(O/'coverage.json').write_text(json.dumps(dict(complete=len(set(targets)&complete),total=len(targets),missing_rows=remaining,missing_coefficients=missing),indent=2))
if remaining:raise SystemExit('remaining rows: '+str(len(remaining)))
subprocess.run([str(ROOT/'build/verify_reductions'),str(O/'reduction.json'),a.verification,str(O/'fresh_verification.json')],check=True)
print(json.dumps(dict(status='complete',targets=len(targets),fresh_verified=True)),flush=True)
