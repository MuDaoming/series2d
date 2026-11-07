"""Fit numerical DE entries against denominator witnesses from other entries.
Never reads dp_verification. Each fit reserves three points internally.
"""
import sys,json,time,collections,itertools
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'tests'))
from rational_learning import make_learner
P=1000003; CENTER=200001; R=ROOT/'results';OUT=R/'dp_pole_route';OUT.mkdir(exist_ok=True)
L=make_learner(P,CENTER); start=time.monotonic()
points={}
for directory in ['dp_full','dp_extra_reductions']:
 for path in sorted((R/directory/'points').glob('point*.json')):
  d=json.loads(path.read_text());t=int(d['t_point']);points.setdefault(t,{}).update(d['reduction_values'])
known={}
for fn in ['dp_full/reductions_partial33.json','dp_extra_reductions/partial41.json','dp_reduction_run/reduction.json','dp_pole_route/shared_reductions.json']:
 d=json.loads((R/fn).read_text())
 for row,v in d['reductions'].items():known.setdefault(row,{}).update(v)
rowdens=collections.defaultdict(collections.Counter);coldens=collections.defaultdict(collections.Counter)
for row,coeffs in known.items():
 for col,v in coeffs.items():
  den=tuple(v['denominator_tau']);rowdens[row][den]+=1;coldens[col][den]+=1
pending=json.loads((R/'dp_reduction_run/reduction.json').read_text())['pending_targets']
rows=pending
new={}; trials=0;counts=collections.Counter()
for ri,row in enumerate(rows):
 ts=sorted(t for t,v in points.items() if row in v)
 if len(ts)<33:continue
 cols=set().union(*(set(points[t][row]) for t in ts))
 candidates=list(rowdens[row].most_common(30))
 for col in sorted(cols):
  if col in known.get(row,{}):continue
  vals=[points[t][row].get(col,0) for t in ts]
  ds=[d for d,n in candidates]+[d for d,n in coldens[col].most_common(20)]
  # Shared row/column denominators are much tighter than a whole-block LCM.
  base=list(dict.fromkeys(ds));ds=base[:]
  for a,b in itertools.product(base[:5], [d for d,n in coldens[col].most_common(5)]):
   ds.append(tuple(int(x) for x in L.lcm(L.poly(list(a)),L.poly(list(b)))))
  expanded=[]
  for den in list(dict.fromkeys(ds)):
   dp=L.poly(list(den))
   for factor,power in dp.factor()[1]:
    for candidate in [dp//factor,dp*factor,dp*factor*factor]:
     expanded.append(tuple(int(x) for x in candidate))
  ds+=expanded
  for den in dict.fromkeys(ds):
   if len(den)>100:continue
   trials+=1;f=L.fit(L.poly(list(den)),ts,vals)
   if f is None:continue
   v={'numerator_tau':list(f.num),'denominator_tau':list(f.den)}
   new.setdefault(row,{'rhs':{}})['rhs'][col]=v;known.setdefault(row,{})[col]=v
   dd=tuple(f.den);rowdens[row][dd]+=1;coldens[col][dd]+=1;counts[len(ts)]+=1
   print(json.dumps({'row':row,'col':col,'degree':[len(f.num)-1,len(f.den)-1],'points':len(ts),'new':sum(counts.values()),'seconds':time.monotonic()-start}),flush=True)
   break
 if ri%10==0 or row in pending:
  (OUT/'shared_reductions_v2.json').write_text(json.dumps({'prime':P,'t0':'2/5','reductions':{r:v['rhs'] for r,v in new.items()},'new_count':sum(counts.values()),'counts_by_points':dict(counts),'trials':trials,'seconds':time.monotonic()-start}))
print('DONE',sum(counts.values()),trials,time.monotonic()-start,flush=True)
(OUT/'shared_reductions_v2.json').write_text(json.dumps({'prime':P,'t0':'2/5','reductions':{r:v['rhs'] for r,v in new.items()},'new_count':sum(counts.values()),'counts_by_points':dict(counts),'trials':trials,'seconds':time.monotonic()-start}))
