"""Vacuum m,0,m: independent Gamma evaluation and original J/FI quadratures."""
import json
from pathlib import Path
import mpmath as m
m.mp.dps=40
h=m.mpf('1.3');mass=m.mpf(1)
pref=-m.gamma(3-2*h)*mass**(2*h-3)
# Cheng-Wu alpha1+alpha3=1, alpha2 in [0,infinity).
cheng_wu=pref*m.beta(2-h,2-h)/(h-1)
tad=-m.gamma(1-h)*mass**(h-1)
ibp=(h-1)/(2*h-3)/mass*tad*tad
assert abs((cheng_wu-ibp)/ibp)<m.mpf('1e-35')

def J(t):
 a=t*(1-t)
 # z=1-u^10, collect powers before evaluating at endpoints.
 def integrand(u):
  q=u**10;z=1-q
  return 10*u**(10*(2-h)-1)*(a*q+z)**(-h)*(t*q+z)**(2*h-3)
 return pref*m.quad(integrand,[0,m.mpf('.5'),1])
# Independent Appell F1 using its convergent transformed series in mpmath.
rows=[]
for t in [m.mpf('.37'),m.mpf('.5')]:
 a=t*(1-t);x=1-1/a;y=-(1-t)/t
 # z -> 1-z transforms F1 arguments to x/(x-1),y/(y-1).
 val=pref*t**(2*h-3)*a**(-h)/(2-h)*(1-x)**(-h)*(1-y)**(-(3-2*h))*m.appellf1(2-h,h,3-2*h,3-h,x/(x-1),y/(y-1))
 direct=J(t);err=abs((val-direct)/direct);assert err<m.mpf('1e-28'),err
 rows.append(dict(t=str(t),J=m.nstr(direct,30),relative_error=m.nstr(err,6)))
report={'d':str(2*h),'m2':'1','Cheng_Wu_exact_Gamma':m.nstr(cheng_wu,35),'IBP_tadpole_product':m.nstr(ibp,35),'J_analytic_checks':rows,'FI_relation':'I111=(d-2)/(2*(d-3)*m2)*T^2'}
Path(__file__).with_name('checks.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report,indent=2))
