"""Threshold bubble DE, explicit J differential identity, and FI reduction."""
import json
from pathlib import Path
import sympy as s
import mpmath as mp
base=Path(__file__).resolve().parent;t,x,h,m=s.symbols('t x h m')
F=t*(m-x+x*t)
lt=(h-2)*s.diff(F,t)/F;lx=(h-2)*s.diff(F,x)/F
acoef=((h-2)*x+(h-1)*m)/(x*(x-m));g=m*t/(x*(m-x))
assert s.cancel(lx-acoef-s.diff(g,t)-g*lt)==0
assert s.cancel((2*h-3)*m-(x-m)*(x*lx-(h-2))-m*(1+t*lt))==0
# Numerical check of the full one-loop DE including the tadpole term.
mp.mp.dps=45;hh=mp.mpf('1.3');xx=-mp.mpf('.7');mm=mp.mpf('1.2')
def density(u):return 10*u**(10*(hh-1)-1)*(mm-xx+xx*u**10)**(hh-2)
b=mp.gamma(2-hh)*mp.quad(density,[0,mp.mpf('.5'),1])
bp=mp.gamma(2-hh)*mp.quad(lambda u:density(u)*(hh-2)*(u**10-1)/(mm-xx+xx*u**10),[0,mp.mpf('.5'),1]);T=-mp.gamma(1-hh)*mm**(hh-1)
lhs=xx*(xx-mm)*bp;rhs=((hh-2)*xx+(hh-1)*mm)*b-(hh-1)*T;err=abs((lhs-rhs)/rhs);assert err<mp.mpf('1e-35'),err
native=json.load(open(base/'j_result.json'));p=native['point']['prime'];tt=int(native['point']['t']);dd=26*pow(7,-1,p)%p;de=native['derivatives']['J{1,1,1,1}'];lower=native['reductions']['J{2,1,0,1}']
actual={k:(tt*v-(dd-4 if k=='J{1,1,1,1}' else 0))%p for k,v in de.items()};actual={k:v for k,v in actual.items() if v};assert actual==lower
report={'generic_d_m2_s_symbolic_bubble_identity':True,'bubble_DE_numeric_relative_error':mp.nstr(err,8),'J_identity':'d(t J1111)/dt=(d-3)J1111+J2101','J_identity_all_basis_coefficients_checked':True,'FI_identity':'I1111=-I1102/(d-3)+(d-2)/(2*(d-3)*m2)*I1011','endpoint':'[t J1111]_0^1=(d-2)/(2*m2)*T*Bubble_mm(Q2)','t_dependent_fractional_mass':'M2(t)=m2/(1-t)'}
if (base/'blade_checks.json').exists():
 blade=json.load(open(base/'blade_checks.json'));assert blade['identity_exactly_zero'];report['independent_Blade_full_reduction_identity']=True
(base/'checks.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report,indent=2))
