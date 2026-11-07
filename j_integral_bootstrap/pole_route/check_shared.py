import json,sys
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'tests'))
from rational_learning import make_learner
R=ROOT/'results';L=make_learner(1000003,200001)
for stem,sec,extra,outsec in [('shared_reductions','reduction_values','dp_extra_reductions','reductions'),('shared_de','DE_values','dp_extra_de','de')]:
 d=json.load(open(R/'dp_pole_route'/f'{stem}.json')); points={}
 for dr in ['dp_full',extra]:
  for f in (R/dr/'points').glob('point*.json'):
   v=json.load(open(f));points.setdefault(int(v['t_point']),{}).update(v[sec])
 total=heldout8=0;failed=[];high=[]
 for row,vs in d[outsec].items():
  if outsec=='de':vs=vs['rhs']
  ts=sorted(t for t in points if row in points[t]);train=[t for t in ts if t<=49]
  for col,c in vs.items():
   total+=1;vals=[points[t][row].get(col,0) for t in train]
   f=L.fit(L.poly(c['denominator_tau']),train,vals)
   if f is None:high.append([row,col]);continue
   if list(f.num)!=c['numerator_tau'] or list(f.den)!=c['denominator_tau']:failed.append([row,col]);continue
   for t in ts:
    x=(t-200001)%1000003
    if int(L.poly(list(f.num))(x))!=points[t][row].get(col,0)*int(L.poly(list(f.den))(x))%1000003:failed.append([row,col]);break
   else:
    if len(ts)==41:heldout8+=1
 report={'coefficients':total,'first33_plus_extra8':heldout8,'needs_more_than_first33':len(high),'failed':failed}
 (R/'dp_pole_route'/f'{stem}_check.json').write_text(json.dumps(report));print(stem,report)
