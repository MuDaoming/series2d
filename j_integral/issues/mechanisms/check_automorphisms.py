"""Exact U,F-preserving propagator permutations for controlled deformations."""
import json,itertools,sys
from pathlib import Path
import sympy as s
base=Path(__file__).resolve().parent;root=base.parents[2]
sys.path.insert(0,str(root/'j_integral/issues/sector_symmetries'))
from audit import polynomials
rows=[]
for p in sorted(base.glob('*/*/request.json')):
 if p.parent.parent.name=='branch_choice':continue
 f=json.load(open(p));params={s.Symbol(k):s.Rational(v) for k,v in f['point'].get('parameters',{}).items()}
 f['family']['gram']=[[str(s.sympify(v).subs(params)) for v in row] for row in f['family']['gram']]
 for prop in f['family']['propagators']:prop['mass2']=str(s.sympify(prop['mass2']).subs(params))
 x,U,F=polynomials(f);br=f['j']['branches'];groups=[[i for i in range(len(x)) if br[i]==b] for b in range(3)];orig=dict(s.Poly(F,*x).terms());autos=[]
 for bp in itertools.permutations(range(3)):
  if any(len(groups[b])!=len(groups[bp[b]]) for b in range(3)):continue
  for orders in itertools.product(*(itertools.permutations(groups[bp[b]]) for b in range(3))):
   mapping={i:j for b in range(3) for i,j in zip(groups[b],orders[b])};trans={}
   for ex,c in orig.items():
    q=[0]*len(x)
    for i,k in enumerate(ex):q[mapping[i]]=k
    trans[tuple(q)]=c
   if trans==orig:autos.append({'branch_map':list(bp),'propagator_map':{str(i+1):j+1 for i,j in mapping.items()}})
 assert autos
 rows.append({'case':str(p.parent.relative_to(base)),'automorphisms':autos});print(rows[-1]['case'],len(autos),flush=True)
(base/'automorphism_checks.json').write_text(json.dumps(rows,indent=2)+'\n')
