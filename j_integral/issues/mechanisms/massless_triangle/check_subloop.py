"""A massless triangle with one null external leg is a divided difference of bubbles."""
import json
from pathlib import Path
import mpmath as m
m.mp.dps=45;h=m.mpf('2.3');A=-m.mpf('1.2');B=-m.mpf('2.7')
# alpha3=u, alpha1=(1-u)y, alpha2=(1-u)(1-y).
# u=w^10 smooths the integrable u=0 endpoint.
beta_integral=m.quad(lambda w:10*w**(10*(h-2)-1)*(1-w**10)**(h-2),[0,m.mpf('.5'),1])
y_integral=m.quad(lambda y:(-A*y-B*(1-y))**(h-3),[0,1])
triangle=-m.gamma(3-h)*beta_integral*y_integral
bubble=lambda s:m.gamma(2-h)*m.gamma(h-1)**2/m.gamma(2*h-2)*(-s)**(h-2)
rhs=-(2*h-3)/(h-2)*(bubble(A)-bubble(B))/(A-B)
err=abs((triangle-rhs)/triangle);assert err<m.mpf('1e-30'),err
result={'dimension':str(2*h),'A':str(A),'B':str(B),'triangle_direct':m.nstr(triangle,35),'divided_difference':m.nstr(rhs,35),'relative_error':m.nstr(err,8),'identity':'Triangle(0,A,B)=-2*(d-3)/(d-4)*(Bubble(A)-Bubble(B))/(A-B)'}
Path(__file__).with_name('subloop_checks.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))
