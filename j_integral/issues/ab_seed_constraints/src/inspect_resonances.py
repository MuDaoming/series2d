"""Inspect diagonal resonant primitive exponents; no target J reduction is run."""
import argparse,collections,json
from pathlib import Path
from analytic_bootstrap import Kd
from analytic_bootstrap.poly import poly_der,poly_divmod
from analytic_bootstrap.rational import factor_powers
parser=argparse.ArgumentParser();parser.add_argument('snapshots',nargs='+',type=Path);parser.add_argument('--scan',type=int,default=128);parser.add_argument('--output',type=Path,required=True);args=parser.parse_args()
all_results=[]
for path in args.snapshots:
 data=json.loads(path.read_text());p=data['prime'];finite=[];infinity=[];nonfuchsian=[];count_factors=0
 for i,row in enumerate(data['connection']):
  diagonal=row[i];num,den=diagonal['n'],diagonal['d']
  if not any(num):
   infinity.append({'component':i,'master':data['masters'][i],'primitive_degree':0})
   continue
  for f,power in factor_powers(den,p):
   count_factors+=1
   if power!=1:
    nonfuchsian.append({'component':i,'factor':list(f),'order':power});continue
   rest,rem=poly_divmod(den,list(f),p);assert not any(rem)
   kd=Kd(p,tuple(f));ratio=kd.mul(kd.reduce_poly(num),kd.inv(kd.mul(kd.reduce_poly(rest),kd.reduce_poly(poly_der(list(f),p)))))
   coefficients=kd.padded(ratio)
   if not any(coefficients[1:]) and 1<=coefficients[0]<=args.scan:
    n=coefficients[0];finite.append({'component':i,'master':data['masters'][i],'factor':list(f),'primitive_denominator_power':n,'resonant_integrand_pole_order':n+1,'numerator_basis_size':len(f)-1})
  assert len(num)<len(den),'diagonal is not strictly proper at infinity'
  lam=num[-1]*pow(den[-1],-1,p)%p if len(num)+1==len(den) else 0
  degree=(-lam)%p
  if degree<=args.scan:infinity.append({'component':i,'master':data['masters'][i],'primitive_degree':degree})
 expected=sum(x['numerator_basis_size'] for x in finite if x['primitive_denominator_power']<=16)+sum(x['primitive_degree']<=16 for x in infinity)
 assert expected==len(data['reference_seeds']),('default completion seed list mismatch',path,expected,len(data['reference_seeds']))
 result={'snapshot':str(path),'point':data['point'],'FBI_rank':len(data['masters']),'positive_finite_primitive_exponents_scanned_through':args.scan,'finite_factors_checked':count_factors,'finite_primitive_power_counts':dict(collections.Counter(x['primitive_denominator_power'] for x in finite for _ in range(x['numerator_basis_size']))),'infinity_primitive_degree_counts':dict(collections.Counter(x['primitive_degree'] for x in infinity)),'finite_nonfuchsian_diagonal_entries':nonfuchsian,'default_completion_seed_count_verified':expected,'finite_resonances':finite,'infinity_resonances':infinity}
 all_results.append(result)
 print(path,'finite powers',result['finite_primitive_power_counts'],'infinity degrees',result['infinity_primitive_degree_counts'],'default seed count',expected)
args.output.parent.mkdir(parents=True,exist_ok=True);args.output.write_text(json.dumps({'scan_bound':args.scan,'results':all_results},indent=2)+'\n')
