"""General-s sunrise identity: identify the term that vanishes at s=m2."""
import json
from pathlib import Path
import sympy as s
z,t,h,q=s.symbols('z t h q');a=t*(1-t);L=a+(1-a)*z;R=L-q*a*z*(1-z)
Az=(1-h)/(z-1)+(3-3*h)*s.diff(L,z)/L+(2*h-3)*s.diff(R,z)/R
K=(2*h-3)*z*L/R;K2=(2*h-3)*(1-z)*(1-t)*L/R
H=(1-z)*(t*(t-1)**2-(t-1)*(t*t-t+1)*z)
D=1-a*(q-1);A=-(3*h-4)*(t-1);C=2*(t-1);boundary=-t*(t-1)**2
assert s.cancel(D*K2-A-C*K-s.diff(H,z)-H*Az)==0
native=json.load(open(Path(__file__).with_name('j_result.json')));p=native['point']['prime']
def mod(v):
 n,d=s.fraction(s.cancel(v.subs({t:17,h:s.Rational(13,7),q:s.Rational(3,2)})))
 return int(n)%p*pow(int(d)%p,-1,p)%p
expected={name:mod(v/D) for name,v in {'J{1,1,1}':A,'J{1,1,2}':C,'Bz0{1,1,1}':boundary}.items()}
assert expected==native['reductions']['J{1,2,1}']
report={'generic_dimension_symbolic_identity':True,'generic_s_symbolic_identity':True,'independent_numeric_IBP_match':True,'sigma':'s/m2','identity':'[1-t(1-t)(sigma-1)]m2 J121=-(3h-4)(t-1)J111+2(t-1)m2 J112-t(t-1)^2 B','integrated_identity':'I112=(3d-8)/(8m2)I111-(d-2)^2/(16(d-3)m2^2)T^2+(sigma-1)/2 * integral[t(1-t)J121 dt]','explanation':'The extra moment vanishes from the relation only at sigma=1; equal-mass permutations alone do not suffice away from this point.'}
Path(__file__).with_name('identity_checks.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report,indent=2))
