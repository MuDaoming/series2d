import json,sys
from pathlib import Path
from analytic_bootstrap import RatMat,RatFunc,Kd,build_local_expansion,stable_jets
from analytic_bootstrap.rational import factor_valuation,rational_local_coefficients
p=Path(sys.argv[1]);data=json.loads(p.read_text());summary=json.loads(Path(sys.argv[2]).read_text());prime=data['prime'];n=len(data['masters'])
def rf(x):return RatFunc.make(x['n'],x['d'],prime)
B=RatMat({(i,j):rf(x) for i,row in enumerate(data['connection']) for j,x in enumerate(row) if any(x['n'])},(n,n),prime)
name=summary['unrepresented_physical_rows'][0];row=[rf(x) for x in data['rows'][name]];bad=[]
for f in summary['finite_factors']:
 if len(f)!=2:continue
 k=max(factor_valuation(x.den,f,prime) for x in row)
 if not k:continue
 kd=Kd(prime,f);exp=build_local_expansion(B,tuple(f),kd,k+20);jet=stable_jets(exp,k,k+20,route='auto').jets
 coeff=[rational_local_coefficients(x,tuple(f),kd,-1) if factor_valuation(x.den,f,prime)>0 else {} for x in row];found=[]
 for power in range(-k,0):
  cv=[0]*(k*n)
  for j,c in enumerate(coeff):
   for s,value in c.items():
    d=power-s
    if 0<=d<k:cv[d*n+j]=(cv[d*n+j]+kd.padded(value)[0])%prime
  vs=[sum(a*b for a,b in zip(v,cv))%prime for v in jet.table()]
  if any(vs):found.append({'power':power,'value':next(v for v in vs if v)})
 if found:bad.append({'factor':f,'pole_order':k,'analytic_solution_count':jet.nrows,'forbidden':found})
print(json.dumps({'target':name,'violations':bad},indent=2))
