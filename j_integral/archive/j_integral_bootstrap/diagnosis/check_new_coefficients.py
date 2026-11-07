import sys,json
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'tests'))
from reference_space import RatFunc
from analytic_bootstrap.poly import poly_divmod,poly_degree,poly_mul
from analytic_bootstrap.rational import factor_powers
out=ROOT/'results/dp_diagnosis';c=json.loads((out/'capsule.json').read_text());p=c['prime'];nr,nc=c['A_shape'][0],c['B_shape'][0];subs=c['masters'][:nc];tops=c['masters'][nc:];rp={m:i for i,m in enumerate(tops)};cp={m:i for i,m in enumerate(subs)}
seen={};bad=[]
for name in ['shared_de.json','shared_de_pending_v2.json','preconditioned_de.json','preconditioned_de_all.json','supplemental_v1.json']:
 path=ROOT/'results/dp_pole_route'/name
 if not path.exists():continue
 for m,row in json.loads(path.read_text()).get('de',{}).items():
  if m not in rp:continue
  for j,v in row['rhs'].items():
   if j not in cp:continue
   i,k=rp[m],cp[j];r=RatFunc.make(v['numerator_tau'],v['denominator_tau'],p);seen[i,k]=r
   q,rem=poly_divmod(c['denominators'][i][k],r.den,p)
   if any(rem):bad.append({'source':name,'row':m,'col':j,'reason':'missing-denominator-factor','actual_den':r.den,'actual_factors':factor_powers(r.den,p),'candidate_factors':factor_powers(c['denominators'][i][k],p)})
   elif poly_degree(poly_mul(r.num,q,p))>c['numerator_caps'][i][k]:bad.append({'source':name,'row':m,'col':j,'reason':'cap','actual_degree':poly_degree(poly_mul(r.num,q,p)),'cap':c['numerator_caps'][i][k]})
result={'new_available':len(seen),'bad':bad};(out/'new_coefficients_envelope.json').write_text(json.dumps(result,indent=2));print(json.dumps(result))
