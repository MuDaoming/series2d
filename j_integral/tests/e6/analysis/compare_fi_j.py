"""Audit e6 counts and find an exact last-variable IBP modulo lower sectors."""
import json,gzip
from pathlib import Path
from flint import nmod_poly,nmod_mat
from analytic_bootstrap import RatFunc
ROOT=Path(__file__).resolve().parents[4];CASE=ROOT/'j_integral/tests/e6';OUT=CASE/'analysis';data=json.loads(gzip.decompress((CASE/'reference/physical_de.json.gz').read_bytes()));fi=json.loads((ROOT/'j_integral/archive/j_integral_bootstrap/reference/fi_master_counts/e6.json').read_text());p=data['prime'];zero=RatFunc.make([0],[1],p);one=RatFunc.make([1],[1],p)
def rf(v):return RatFunc.make(v['numerator_tau'],v['denominator_tau'],p)
def indices(n):return tuple(map(int,n[n.index('{')+1:-1].split(',')))
def sector(nu):return tuple(int(v>0) for v in nu)
from fractions import Fraction
family=json.loads((CASE/'input/family_request.json').read_text());blade=json.loads((ROOT/'j_integral/archive/j_integral_bootstrap/reference/fi_master_counts/e6_input.json').read_text());f=family['family'];E=len(f['gram']);momenta=['l1','l2']+['p'+str(i+1) for i in range(E)]
def linear(values):return '+'.join('('+str(Fraction(v))+')*'+name for v,name in zip(values,momenta) if Fraction(v)) or '0'
assert [ '('+linear(prop['loop']+prop['shift'])+')^2-('+str(prop['mass2'])+')' for prop in f['propagators'] ]==blade['propagators'][:7]
assert Fraction(family['point']['dimension'])==Fraction(blade['dimension'])
assert {name:Fraction(value) for name,value in blade['scalar_products']}=={'p'+str(i+1)+'*p'+str(j+1):Fraction(f['gram'][i][j]) for i in range(E) for j in range(i,E)}
from collections import Counter
F=Counter(sector(nu[:7]) for nu in fi['master_indices']);J=Counter(sector(indices(name)) for name in data['masters'] if name.startswith('J'))
report=dict(same_FI_family_dimension_and_kinematics_checked=True,FI_total=sum(F.values()),J_total=len(data['masters']),ordinary_J_total=sum(J.values()),boundary=[n for n in data['masters'] if n.startswith('B')],by_active_propagators=[dict(active=k,FI=sum(v for m,v in F.items() if sum(m)==k),ordinary_J=sum(v for m,v in J.items() if sum(m)==k)) for k in range(3,8)],sector_differences=[dict(mask=''.join(map(str,m)),FI=F[m],ordinary_J=J[m]) for m in sorted(set(F)|set(J),key=lambda m:(sum(m),m)) if F[m]!=J[m]])
(OUT/'counts.json').write_text(json.dumps(report,indent=2)+'\n')
mask=(0,1,1,0,1,0,1);top=[n for n in data['masters'] if n.startswith('J') and sector(indices(n))==mask];A=[[rf(data['de'][r]['rhs'][c]) if c in data['de'][r]['rhs'] else zero for c in top] for r in top];r=len(top)
from fractions import Fraction
q=Fraction(data['t0']);center=q.numerator*pow(q.denominator,-1,p)%p
t=RatFunc.make([center,1],[1],p);endpoint=t*(one-t)
def lcm(a,b):return (a//a.gcd(b))*b
common=nmod_poly([1],p)
for row in A:
 for value in row:common=lcm(common,nmod_poly(value.den,p))
attempts=[]
for den in [nmod_poly([1],p),common,common*common]:
 for degree in range(13):
  columns=[];generators=[]
  for j in range(r):columns.append([-one if i==j else zero for i in range(r)])
  for i in range(r):
   for k in range(degree+1):
    value=endpoint*RatFunc.make([0]*k+[1],[int(c) for c in den],p);g=[value if i==j else zero for j in range(r)];generators.append(g)
    columns.append([g[j].derivative()+sum((g[a]*A[a][j] for a in range(r)),zero) for j in range(r)])
  equations=[]
  for j in range(r):
   D=nmod_poly([1],p)
   for col in columns:D=lcm(D,nmod_poly(col[j].den,p))
   polys=[(D//nmod_poly(col[j].den,p))*nmod_poly(col[j].num,p) for col in columns]
   for k in range(max(v.degree() for v in polys)+1):equations.append([int(poly[k]) for poly in polys])
  null,n=nmod_mat(equations,p).nullspace();attempts.append(dict(denominator_degree=den.degree(),polynomial_degree=degree,unknowns=len(columns),nullity=n))
  picked=next((k for k in range(n) if any(int(null[j,k]) for j in range(r))),None)
  if picked is None:continue
  vector=[int(null[j,picked]) for j in range(len(columns))];scale=pow(next(x for x in vector[:r] if x),-1,p);vector=[v*scale%p for v in vector];g=[sum((generators[k][j]*vector[r+k] for k in range(len(generators))),zero) for j in range(r)]
  for j in range(r):assert (g[j].derivative()+sum((g[a]*A[a][j] for a in range(r)),zero)-vector[j]).is_zero()
  lower=sorted({m for n in top for m in data['de'][n]['rhs'] if m not in top});remainder={}
  for m in lower:
   v=sum((g[i]*rf(data['de'][name]['rhs'][m]) for i,name in enumerate(top) if m in data['de'][name]['rhs']),zero)
   if not v.is_zero():remainder[m]=dict(numerator_tau=v.num,denominator_tau=v.den)
  cert=dict(prime=p,t0=data['t0'],sector=''.join(map(str,mask)),top_masters=top,FI_sector_masters=[nu[:7] for nu in fi['master_indices'] if sector(nu[:7])==mask],constant_top_coefficients=dict(zip(top,vector[:r])),primitive={name:dict(numerator_tau=v.num,denominator_tau=v.den) for name,v in zip(top,g)},lower_terms_in_derivative=remainder,identity='d(sum_i g_i J_i)/dt = sum_i c_i J_i + sum_lower R_lower(t) J_lower(t); exact rational identity, endpoint terms NOT discarded',exact_top_identity_checked=True,attempts=attempts)
  # Validate every coefficient in the full79-dimensional connection, not just A_top.
  for name in data['masters']:
   lhs=sum((g[i]*rf(data['de'][row]['rhs'][name]) for i,row in enumerate(top) if name in data['de'][row]['rhs']),zero)
   if name in top:lhs=lhs+g[top.index(name)].derivative();rhs=RatFunc.make([vector[top.index(name)]],[1],p)
   else:rhs=rf(remainder[name]) if name in remainder else zero
   assert (lhs-rhs).is_zero(),name
  def evaluate(v,x):
   def poly(cs):
    answer=0
    for c in reversed(cs):answer=(answer*x+c)%p
    return answer
   return poly(v.num)*pow(poly(v.den),-1,p)%p
  points=[];checks=0
  for path in sorted((CASE/'runs/original_cli_001/numeric/points').glob('point*.json'))[:2]:
   native=json.loads(path.read_text());a=int(native['point']['t']);x=(a-center)%p
   for name in data['masters']:
    lhs=sum(evaluate(g[i],x)*native['derivatives'][row].get(name,0) for i,row in enumerate(top))%p
    if name in top:lhs=(lhs+evaluate(g[top.index(name)].derivative(),x))%p;rhs=vector[top.index(name)]
    else:rhs=evaluate(rf(remainder[name]),x) if name in remainder else 0
    assert lhs==rhs,(a,name);checks+=1
   points.append(a)
  cert.update(full79_coefficient_identity_checked=True,independent_cached_IBP_points=points,independent_scalar_checks=checks,endpoint_terms_kept=True)
  (OUT/'t_ibp_certificate.json').write_text(json.dumps(cert,indent=2)+'\n');print(json.dumps(cert),flush=True);raise SystemExit(0)
print('No certificate in tested finite ansatz',attempts)
