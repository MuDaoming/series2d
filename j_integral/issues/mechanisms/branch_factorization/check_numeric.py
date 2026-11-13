"""Independent Euler-integral verification of exact t separation."""
import json
from pathlib import Path
import mpmath as m
m.mp.dps=45
rows=[]
for A,B,C,h in [(1,2,1,m.mpf('1.3')),(2,3,1,m.mpf('1.8')),(2,3,2,m.mpf('2.3'))]:
 v=A+B+C;f=m.mpf(2);k=m.mpf('.7')
 kernel=lambda w:w**(2*h-A-C-1)*(1-w)**(A+C-h-1)*(f-(1-w)*k)**(2*h-v)
 const=m.quad(kernel,[0,m.mpf('.5'),1])
 for t in [m.mpf('.23'),m.mpf('.61')]:
  a=t*(1-t)
  def integrand(z):
   L=a+(1-a)*z
   return t**(A-1)*(1-t)**(C-1)*(1-z)**(A+C-h-1)*z**(2*h-A-C-1)*L**(v-3*h)*(L*f-a*(1-z)*k)**(2*h-v)
  direct=m.quad(integrand,[0,m.mpf('.5'),1]);closed=t**(h-C-1)*(1-t)**(h-A-1)*const;err=abs((direct-closed)/closed);assert err<m.mpf('1e-25'),err
  rows.append({'A':A,'B':B,'C':C,'d':str(2*h),'t':str(t),'relative_error':m.nstr(err,8)})
alpha=m.mpf('.8');beta=m.mpf('.6');u=m.mpf('.37')
lhs=m.quad(lambda t:(t-u)**(alpha-1)*(1-t)**(beta-1)*t**(-alpha-beta),[u,(u+1)/2,1])
rhs=m.beta(alpha,beta)*u**(-beta)*(1-u)**(alpha+beta-1)
euler_error=abs((lhs-rhs)/rhs);assert euler_error<m.mpf('1e-25'),euler_error
Path(__file__).with_name('numeric_checks.json').write_text(json.dumps({'precision':m.mp.dps,'checks':rows,'old_projection_Euler_integral_relative_error':m.nstr(euler_error,8)},indent=2)+'\n');print(json.dumps(rows,indent=2))
