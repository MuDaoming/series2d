"""Check the historically blocked 6x47 C against old candidate envelopes."""
import sys,json,time
from pathlib import Path
from fractions import Fraction
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'tests'))
from reference_space import RatFunc
from rational_learning import make_learner
from analytic_bootstrap.poly import poly_divmod,poly_mul,poly_degree
from analytic_bootstrap.rational import factor_powers
r=ROOT/'results';out=r/'dp_diagnosis/restricted';out.mkdir(exist_ok=True)
d=json.loads((r/'dp_final/connection_frozen.json').read_text());mapping=json.loads((r/'dp_restricted/mapping.json').read_text());p=d['prime'];f=Fraction(d['t0']);center=f.numerator*pow(f.denominator,-1,p)%p;L=make_learner(p,center);one,ep=L.one,L.ep
rf=lambda v:RatFunc.make(v['numerator_tau'],v['denominator_tau'],p)
top=mapping['original_top'];topset=set(top)
graph={m:set(row['rhs'])-{m} for m,row in d['de'].items()};deps=set();todo=list(top)
while todo:
 m=todo.pop()
 for j in graph[m]:
  if j not in topset and j not in deps:deps.add(j);todo.append(j)
sub=[m for m in d['known_rows'] if m in deps]
D={(m,j):rf(v) for m,row in d['de'].items() if m in deps or m in topset for j,v in row['rhs'].items() if j in deps or j in topset}
# SCC partition of the verified connection.
from sector_de import components
masters=d['masters'];idx={m:i for i,m in enumerate(masters)};comp={masters[i]:k for k,g in enumerate(components({idx[m]:{idx[j] for j in v} for m,v in graph.items()})) for i in g}
learned={}
for file in ['dp_pole_route/merged_de49.json','dp_pole_route/refined_pending.json']:
 for m,row in json.loads((r/file).read_text())['de'].items():
  if m in topset:
   for j,v in row['rhs'].items():learned[m,j]=rf(v)
qs=one;bare={m:one for m in top};rowq={m:one for m in top};cols={j:one for j in sub};diag={j:one for j in sub}
for m in top:
 for j in top:
  if (m,j)in D:bare[m]=L.lcm(bare[m],L.poly(D[m,j].den));qs=L.lcm(qs,L.poly(D[m,j].den))
 rowq[m]=bare[m]
 for j in sub:
  if (m,j)in learned:rowq[m]=L.lcm(rowq[m],L.poly(learned[m,j].den));qs=L.lcm(qs,L.poly(learned[m,j].den))
for m in sub:
 for j in sub:
  if (m,j) in D:
   cols[j]=L.lcm(cols[j],L.poly(D[m,j].den))
   if comp[m]==comp[j]:diag[j]=L.lcm(diag[j],L.poly(D[m,j].den))
settings=[]
for scope in ['row-diagonal','row-column','block-column']:
 for power in [0,1,2]:
  failures=[];worst_growth=-1000
  for m in top:
   for j in sub:
    if (m,j) not in D:continue
    actual=D[m,j]
    if (m,j)in learned:den=L.poly(learned[m,j].den);cap=len(learned[m,j].num)-1
    else:
     qa=bare[m] if scope=='row-diagonal' else rowq[m] if scope=='row-column' else qs;qb=diag[j] if scope=='row-diagonal' else cols[j]
     den=L.lcm(L.lcm(qa,qb),ep) if power==0 else (qa*qb*ep)**power;den=L.monic(den);cap=int(den.degree())+3
    denarr=[int(v) for v in den];q,rem=poly_divmod(denarr,actual.den,p);growth=len(actual.num)-len(actual.den);worst_growth=max(worst_growth,growth)
    if any(rem):failures.append({'row':m,'col':j,'reason':'denominator','actual_degrees':[len(actual.num)-1,len(actual.den)-1],'candidate_den_degree':int(den.degree()),'actual_factors':factor_powers(actual.den,p),'candidate_factors':factor_powers(denarr,p)})
    elif poly_degree(poly_mul(actual.num,q,p))>cap:failures.append({'row':m,'col':j,'reason':'cap','actual_growth':growth,'cap':cap,'actual_in_candidate_num_degree':poly_degree(poly_mul(actual.num,q,p))})
  settings.append({'scope':scope,'power':power,'failures':failures,'worst_actual_growth':worst_growth})
report={'top':top,'sub':sub,'sub_rank':len(sub),'old_learned_couplings':sum((m,j)in learned for m in top for j in sub),'nonzero_C':sum((m,j)in D for m in top for j in sub),'settings':settings}
(out/'old_envelope.json').write_text(json.dumps(report,indent=2));print(json.dumps({'sub_rank':len(sub),'learned':report['old_learned_couplings'],'nonzero_C':report['nonzero_C'],'settings':[{k:v for k,v in s.items() if k!='failures'}|{'failures':len(s['failures'])} for s in settings]},indent=2))
