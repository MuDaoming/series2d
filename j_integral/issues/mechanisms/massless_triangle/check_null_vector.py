"""Null leg removes the quadratic term; the resulting parameter vector field moves t."""
import json
from pathlib import Path
import sympy as s
x,y,c,b,A,B,q=s.symbols('x y c b A B q')
V=lambda f:s.expand(B*s.diff(f,x)-A*s.diff(f,y))
F=c*(A*x+B*y)+q*x*y
assert s.factor(V(F)-q*(B*y-A*x))==0
t=(x+y)/(x+y+b)
assert s.factor(V(t)-b*(B-A)/(x+y+b)**2)==0
# General subloop polynomial when the two adjacent lines have a null difference.
C=s.symbols('C');Fbox=A*x+B*y+C+q*x*y
assert s.factor(V(Fbox)-q*(B*y-A*x))==0
r={'triangle_F':'c*(A*x+B*y)+p2*x*y','box_F':'A*x+B*y+C+p2*x*y','vector_field':'B*d/dx-A*d/dy','V_F':'p2*(B*y-A*x)','V_t':'b*(B-A)/(x+y+b)^2','exact_checks':True,'interpretation':'At p2=0 the vector field preserves F and yields a simpler subloop IBP, but it moves the fixed-t slice; its integrated identity therefore does not become a pointwise algebraic J relation.'}
Path(__file__).with_name('null_vector_checks.json').write_text(json.dumps(r,indent=2)+'\n');print(json.dumps(r,indent=2))
