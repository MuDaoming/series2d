"""Reuse rational witnesses, reconstruct restricted R, and verify fresh points."""
import sys,json,argparse,subprocess,concurrent.futures
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parent))
from upstream_j import ROOT
from restriction_inputs import expand_name
from restricted_numeric import indices,label
sys.path.insert(0,str(ROOT/'tests'));from rational_learning import make_learner
p=argparse.ArgumentParser();p.add_argument('work');p.add_argument('--finish',action='store_true');a=p.parse_args();O=Path(a.work).resolve();mapping=json.loads((O/'mapping.json').read_text());mp=dict(zip(mapping['original_masters'],mapping['compressed_masters']));active=mapping['active_positions'];known={};R=ROOT/'results'
def compact(name):
 kind,nu=indices(name)
 if any(v and i not in active for i,v in enumerate(nu)):return None
 return label(kind,[nu[i] for i in active])
for f in ['dp_full/reductions_partial33.json','dp_extra_reductions/partial41.json','dp_pole_route/supplemental_v2.json','dp_pole_route/native_reductions45.json','dp_pole_route/native_reductions49.json','dp_reduction_final/reduction.json']:
 for row,vs in json.loads((R/f).read_text())['reductions'].items():
  name=compact(row)
  if name is None:continue
  for col,v in vs.items():
   if col in mp:
    dest=known.setdefault(name,{})
    if mp[col] in dest:assert dest[mp[col]]==v,(name,col)
    dest[mp[col]]=v
(O/'known.json').write_text(json.dumps(dict(prime=1000003,t0='2/5',masters=mapping['compressed_masters'],reductions=known)))
points=sorted((p for p in (O/'points').glob('point*.json') if '.tmp.' not in p.name),key=lambda p:int(p.stem[5:]));manifest=dict(prime=1000003,t0='2/5',masters=mapping['compressed_masters'],section='reduction_values',holdout=8,ratios=True,ratio_budget=8,max_denominator_degree=512,points=[str(p) for p in points],known=[str(O/'known.json')]);(O/'manifest.json').write_text(json.dumps(manifest))
with (O/'reconstruct.log').open('w') as log:subprocess.run([str(ROOT/'build/preconditioned'),str(O/'manifest.json'),str(O/'reconstructed.json')],check=True,stdout=log,stderr=log)
new=json.loads((O/'reconstructed.json').read_text())['reductions']
for n,vs in new.items():known.setdefault(n,{}).update(vs)
original=json.loads((O/'original_targets.json').read_text());targets=[compact(n) for n in original];data=[json.loads(p.read_text()) for p in points];complete={};missing=[]
for row in targets:
 columns=set().union(*(set(d['reduction_values'][row]) for d in data));missing.extend((row,c) for c in columns if c not in known.get(row,{}));complete[row]={c:known[row][c] for c in columns if c in known.get(row,{})}
report=dict(target_rows=len(targets),learning_points=len(points),missing=missing,coefficients=sum(map(len,complete.values())))
(O/'reconstruction_status.json').write_text(json.dumps(report,indent=2));print(json.dumps(report),flush=True)
if missing or not a.finish:sys.exit(0)
L=make_learner(1000003,200001);polys={r:{c:(L.poly(v['numerator_tau']),L.poly(v['denominator_tau'])) for c,v in vs.items()} for r,vs in complete.items()}
def verify(d):
 x=(int(d['t_point'])-200001)%1000003;checks=0
 for row in targets:
  values=d['reduction_values'][row]
  for col in mapping['compressed_masters']:
   if col in polys[row]:
    n,den=polys[row][col];assert int(den(x)) and int(n(x))==values.get(col,0)*int(den(x))%1000003,(row,col,d['t_point'])
   else:assert not values.get(col,0),(row,col,d['t_point'],'implicit-zero')
   checks+=1
 return checks
checks=sum(verify(d) for d in data);(O/'verification').mkdir(exist_ok=True)
def fresh(t):
 out=O/'verification'/f'point{t}.json'
 with (O/'verification'/f'point{t}.log').open('w') as log:subprocess.run([str(ROOT/'build/global_point_subset'),str(O/'family.json'),str(out),'--point',str(t),'--targets-file',str(O/'targets.json'),'--basis-file',str(O/'basis.json'),'--derivative-names-file',str(O/'no_derivatives.json'),'--integration','atomic','--values','1','--compact','1'],check=True,stdout=log,stderr=log)
 return json.loads(out.read_text())
with concurrent.futures.ThreadPoolExecutor(max_workers=2) as pool:freshpoints=list(pool.map(fresh,[1801,1802,1803]))
freshchecks=sum(verify(d) for d in freshpoints)
expanded={expand_name(r,mapping):{expand_name(c,mapping):v for c,v in vs.items()} for r,vs in complete.items()};out=dict(prime=1000003,t0='2/5',reductions=expanded,source='generic restricted physical subfamily; full17/18 compatibility and local independent points verified',internal_holdouts=8)
# Check reconstructed functions against every old full pilot, including zeros
# outside the restricted coordinate basis. No global fresh-validation R is read.
fullchecks=0
for t in range(17,50):
 full=json.loads((R/'dp_full/points'/f'point{t}.json').read_text());x=(t-200001)%1000003
 for original_row in original:
  row=compact(original_row);actual=full['reduction_values'][original_row]
  assert set(actual)<=set(mapping['original_masters'])
  for local_col,original_col in zip(mapping['compressed_masters'],mapping['original_masters']):
   if local_col in polys[row]:
    n,d=polys[row][local_col];assert int(d(x)) and int(n(x))==actual.get(original_col,0)*int(d(x))%1000003
   else:assert not actual.get(original_col,0)
  fullchecks+=len(full['masters_at_point'])
(O/'complete_expanded.json').write_text(json.dumps(out));report.update(full_reference_checks=fullchecks,learning_checks=checks,independent_points=[1801,1802,1803],independent_checks=freshchecks,failures=0,global_fresh_verification_read=False);(O/'checks.json').write_text(json.dumps(report,indent=2));print(json.dumps(report),flush=True)
