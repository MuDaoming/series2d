"""Analytic J formulas: exact diagonal-DE checks and convergent Euler integrals."""
import json
from pathlib import Path
from fractions import Fraction
import mpmath as mp
import sympy as sp
HERE=Path(__file__).resolve().parent;ISSUE=HERE.parent
mp.mp.dps=50
sP=mp.mpf(149)/60;sQ=mp.mpf(38113)/4620;sR=mp.mpf(1)

def phase(power):return mp.exp(-mp.j*mp.pi*power)
def sunrise(s,t,d):
    h=d/2
    return -mp.gamma(3-d)*phase(d-3)*s**(d-3)*mp.beta(d-2,h-1)*(t*(1-t))**(h-2)
def sunrise_z(s,t,d):
    h=d/2;a=t*(1-t)
    integrand=lambda z:z**(d-3)*(1-z)**(h-2)*(a+(1-a)*z)**(3-3*h)
    return -mp.gamma(3-d)*phase(d-3)*s**(d-3)*a**(d-3)*mp.quad(integrand,[0,mp.mpf('.5'),1])
def boundary_constant(d):return mp.gamma(4-d)*phase(d-4)*sQ**(d-4)*mp.beta(d-3,d-3)
def boundary(t,d):
    h=d/2;return boundary_constant(d)*t**(-h)*(1-t)**(h-3)
def boundary_y(t,d):
    h=d/2;U=t*(1-t);X2=1-t
    return mp.gamma(4-d)*X2*U**(4-3*h)*(U*X2*sQ)**(d-4)*phase(d-4)*mp.quad(lambda y:y**(d-4)*(1-y)**(d-4),[0,mp.mpf('.5'),1])
def H(y):return sP-(sP+sQ-sR)*y+sQ*y*y
def z_integrand(z,y,t,d,nu):
    h=d/2;v=sum(nu);n1,n2,n3,n4=nu
    X1=(1-z)*t;X2=(1-z)*(1-t);X3=z;U=X1*X2+X1*X3+X2*X3
    F=-X1*X2*X3*H(y)-U*X2*sQ*y*(1-y)
    pref=(-1)**v*mp.gamma(v-d)/mp.fprod(mp.gamma(n) for n in nu)
    return pref*(1-z)*X1**(n1-1)*X2**(n2+n3-1)*X3**(n4-1)*U**(v-3*h)*(-F)**(d-v)*phase(d-v)*y**(n3-1)*(1-y)**(n2-1)
def euler_result(y,t,d,nu):
    h=d/2;v=sum(nu);n1,n2,n3,n4=nu;q=sQ*y*(1-y)
    pref=(-1)**v*mp.gamma(v-d)*mp.gamma(h-n4)/(mp.gamma(n1)*mp.gamma(n2)*mp.gamma(n3)*mp.gamma(h))
    return pref*t**(n1+n4-h-1)*(1-t)**(h-n1-1)*y**(n3-1)*(1-y)**(n2-1)*q**(d-v)*phase(d-v)*mp.hyp2f1(v-d,n4,h,-t*H(y)/q)
def top_J(t,d,nu):return mp.quad(lambda y:euler_result(y,t,d,nu),[0,mp.mpf('.5'),1])
def relative(a,b):return abs(a-b)/max(abs(a),abs(b),mp.mpf('1e-40'))
def out(v):return dict(real=mp.nstr(mp.re(v),30),imaginary=mp.nstr(mp.im(v),30))

report=dict(precision=mp.mp.dps,normalization='physical FI/FBI/J including Gamma factors; (-s-i0)^a=e^(-i*pi*a)*s^a for positive s',checks=[])
# Verify the rational diagonal entries exactly at the saved finite-field setting.
reference=json.loads((ISSUE/'reference/j_de.json').read_text());p=reference['prime'];tau=sp.Symbol('tau');t=tau+sp.Rational(2,5);h=sp.Rational(13,7)
expected={'J{1,0,1,1}':(h-2)*(1-2*t)/(t*(1-t)),'J{1,1,0,1}':(h-2)*(1-2*t)/(t*(1-t)),'Bz0{1,1,1,1}':(h-3*t)/(t*(t-1))}
for name,value in expected.items():
    saved=reference['de'][name]['rhs'];assert set(saved)=={name}
    N=sum(c*tau**i for i,c in enumerate(saved[name]['numerator_tau']));D=sum(c*tau**i for i,c in enumerate(saved[name]['denominator_tau']));num,den=sp.fraction(sp.cancel(value));difference=sp.Poly(sp.expand(num*D-den*N),tau,modulus=p);assert difference.is_zero
report['exact_diagonal_DE_checks']=list(expected)
# Verify the variable substitution and F identity independently of any DE.
t,z,w,y,sP0,sQ0,sR0=sp.symbols('t z w y sP sQ sR');a=t*(1-t);L=a+(1-a)*z;X1=(1-z)*t;X2=(1-z)*(1-t);X3=z;U=X1*X2+X1*X3+X2*X3;HH=sP0-(sP0+sQ0-sR0)*y+sQ0*y*y;F=-X1*X2*X3*HH-U*X2*sQ0*y*(1-y)
assert sp.cancel(U-(1-z)*L)==0
assert sp.cancel(F+(1-z)**2*(1-t)*(L*sQ0*y*(1-y)+t*z*HH))==0
Z=a*w/(1-(1-a)*w);assert sp.cancel((z/L).subs(z,Z)-w)==0
report['exact_coordinate_and_F_checks']=True
hs=sp.Symbol('h');g=t*(1-t)*(2-hs-t)/((1-hs)*(2-hs));ab=(hs-3*t)/(t*(t-1))
assert sp.cancel(sp.diff(g,t)+g*ab-1)==0
report['exact_boundary_primitive_check']=True
for d in [mp.mpf(26)/7,mp.mpf(11)/2]:
    for t in [mp.mpf('.23'),mp.mpf('.61')]:
        for s in [sP,sR]:
            error=relative(sunrise(s,t,d),sunrise_z(s,t,d));assert error<mp.mpf('1e-30'),error
            report['checks'].append(dict(test='sunrise z integral',d=str(d),t=str(t),s=str(s),relative_error=mp.nstr(error,8)))
        error=relative(boundary(t,d),boundary_y(t,d));assert error<mp.mpf('1e-30'),error
        report['checks'].append(dict(test='boundary y integral',d=str(d),t=str(t),relative_error=mp.nstr(error,8)))
# d=11/2 gives a convergent z integral for all three top objects; general d follows by analytic continuation.
for nu in [(1,1,1,1),(1,1,1,2),(1,1,2,1)]:
    y=mp.mpf('.29');t=mp.mpf('.37');d=mp.mpf(11)/2
    direct=mp.quad(lambda z:z_integrand(z,y,t,d,nu),[0,mp.mpf('.5'),1]);closed=euler_result(y,t,d,nu);error=relative(direct,closed);assert error<mp.mpf('1e-28'),error
    report['checks'].append(dict(test='top z integral to Gauss 2F1',indices=nu,relative_error=mp.nstr(error,8)))
# Full lower FI normalization is also checkable in an ordinary convergence strip.
d=mp.mpf(26)/7;h=d/2
lhs=mp.beta(d-2,h-1)*mp.beta(h-1,h-1);rhs=mp.gamma(h-1)**3/mp.gamma(3*h-3);assert relative(lhs,rhs)<mp.mpf('1e-45')
report['sunrise_FI_Gamma_identity']=True
# There is no common convergence strip for this boundary t integral. Record analytic continuation explicitly.
report['boundary_period_prescription']='Beta-function analytic continuation in independent endpoint exponents; not an ordinary convergent integral'
report['boundary_period_beta_arguments']=['1-d/2','d/2-2'];report['boundary_period_beta_argument_sum']=-1;assert mp.rgamma(-1)==0
report['boundary_period_analytic_continuation']=0
report['boundary_regularization_samples']=[dict(epsilon=str(eps),beta_value=mp.nstr(mp.beta(1-h+eps,h-2+eps),18)) for eps in [mp.mpf('.01'),mp.mpf('.001'),mp.mpf('.0001')]]
report['example_lower_values_at_t_0_37']=dict(J1011=out(sunrise(sR,mp.mpf('.37'),d)),J1101=out(sunrise(sP,mp.mpf('.37'),d)),B1111=out(boundary(mp.mpf('.37'),d)))
(HERE/'checks.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report,indent=2),flush=True)
