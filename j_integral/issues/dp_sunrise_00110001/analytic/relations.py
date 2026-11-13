"""Exact symbolic z-total-derivative certificates, no interpolation."""
import json
from pathlib import Path
import sympy as s
z,t,h=s.symbols('z t h');a=t*(1-t);L=a+(1-a)*z;R=a*(1-z)**2+z
Az=(1-h)/(z-1)+(3-3*h)*s.diff(L,z)/L+(2*h-3)*s.diff(R,z)/R
At=(3-3*h)*s.diff(L,t)/L+(2*h-3)*s.diff(R,t)/R
K=(2*h-3)*z*L/R
G=(1-z)*((2*t-1)-(2*t-1)*(t*t-t+2)/(t*(t-1))*z)
A=-2*(h-1)*(2*t-1)/(t*(t-1));C=2*(2*t-1)/(t*(t-1))
assert s.cancel(At-A-C*K-s.diff(G,z)-G*Az)==0
H=(1-z)*(t*(t-1)**2-(t-1)*(t*t-t+1)*z)
A2=-(3*h-4)*(t-1);C2=2*(t-1)
K2=(2*h-3)*(1-z)*(1-t)*L/R
assert s.cancel(K2-A2-C2*K-s.diff(H,z)-H*Az)==0
# Endpoint primitives vanish at z=1 for Re h<2; at z=0 they equal G(0)B,H(0)B.
assert s.simplify(G.subs(z,0)-(2*t-1))==0
assert s.factor(H.subs(z,0)-t*(t-1)**2)==0
# Weighted boundary period / tadpole^2, using Gamma recursion.
beta_ratio=-(h-1)**2/(3-2*h)
assert s.factor(-beta_ratio/4+(h-1)**2/(4*(2*h-3)))==0
native=json.loads(Path(__file__).with_name('numeric_ibp.json').read_text());prime=native['point']['prime']
def mod(expr):
    value=s.cancel(expr.subs({t:17,h:s.Rational(13,7)}));n,d=s.fraction(value)
    return int(n)%prime*pow(int(d)%prime,-1,prime)%prime
assert {key:mod(value) for key,value in {'J{1,1,1}':A,'J{1,1,2}':C,'Bz0{1,1,1}':-(2*t-1)}.items()}==native['derivatives']['J{1,1,1}']
assert {key:mod(value) for key,value in {'J{1,1,1}':A2,'J{1,1,2}':C2,'Bz0{1,1,1}':-t*(t-1)**2}.items()}==native['reductions']['J{1,2,1}']
result={'exact_symbolic_certificates':True,'dimension':'d=2h','m2':'1 in certificates; dimensions restored in README','G':str(s.factor(G)),'H':str(s.factor(H)),'J111_derivative':{'J111':str(A),'J112':str(C),'B':str(-(2*t-1))},'J121_reduction':{'J111':str(A2),'J112':str(C2),'B':str(-t*(t-1)**2)},'FI_relation':'I112=(3d-8)/(8m2)*I111-(d-2)^2/(16(d-3)m2^2)*T^2','proof':'exact z-total-derivative identity plus t<->1-t symmetry and FI propagator permutation; all endpoint terms retained'}
Path(__file__).with_name('relations.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))
