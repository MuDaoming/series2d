"""Check explicit Lauricella series against the original physical z integral."""
import json
from pathlib import Path
import mpmath as m
import sympy as sp
m.mp.dps=45
N=320

def fd(a,bs,c,xs):
    # z -> 1-z in Euler's formula: all transformed arguments are in (0,1).
    ys=[x/(x-1) for x in xs]
    pref=m.fprod((1-x)**(-b) for x,b in zip(xs,bs))
    coeff=[m.mpf(1)]+[m.mpf(0)]*N
    for b,y in zip(bs,ys):
        seq=[m.mpf(1)]
        for k in range(1,N+1):seq.append(seq[-1]*(b+k-1)*y/k)
        coeff=[sum(coeff[j]*seq[k-j] for j in range(k+1)) for k in range(N+1)]
    weight=m.mpf(1);value=coeff[0]
    for k in range(1,N+1):
        weight*= (c-a+k-1)/(c+k-1)
        value+=weight*coeff[k]
    return pref*value

def analytic(t,h,nu):
    n1,n2,n3=nu;v=sum(nu);a=t*(1-t)
    xs=[1-1/a,-t/(1-t),-(1-t)/t]
    bs=[3*h-v,v-2*h,v-2*h]
    pref=(-1)**v*m.gamma(v-2*h)*m.gamma(n1+n2-h)/(m.gamma(n1)*m.gamma(n2)*m.gamma(v-h))
    return pref*t**(n1-1)*(1-t)**(n2-1)*a**(-h)*fd(n3,bs,v-h,xs)

def direct(t,h,nu):
    n1,n2,n3=nu;v=sum(nu)
    # 1-z=u^14 removes the weak endpoint singularity at d=26/7.
    def integrand(u):
        q=u**14;z=1-q;a=t*(1-t)
        return 14*t**(n1-1)*(1-t)**(n2-1)*u**(14*(n1+n2-h)-1)*z**(n3-1)*(a*q+z)**(v-3*h)*(a*q*q+z)**(2*h-v)
    return (-1)**v*m.gamma(v-2*h)/m.fprod(m.gamma(n) for n in nu)*m.quad(integrand,[0,m.mpf('.5'),1])

t,z=sp.symbols('t z');a=t*(1-t);x1=(1-z)*t;x2=(1-z)*(1-t);u=x1*x2+(x1+x2)*z
assert sp.factor(u-x1*x2*z-(x1+x2)*(x1+z)*(x2+z))==0
assert sp.cancel(u-(1-z)*(a+(1-a)*z))==0
hcheck=m.mpf('.7')
beta_period=-m.gamma(3-2*hcheck)*m.beta(1-hcheck,1-hcheck)
tadpole_period=(2*hcheck-2)*m.gamma(1-hcheck)**2
assert abs((beta_period-tadpole_period)/tadpole_period)<m.mpf('1e-40')
rows=[]
for h,t in [(m.mpf('1.3'),m.mpf('.37')),(m.mpf('1.3'),m.mpf('.5')),(m.mpf(13)/7,m.mpf('.37'))]:
    for nu in [(1,1,1),(1,1,2)]:
        a=analytic(t,h,nu);b=direct(t,h,nu);err=abs((a-b)/b)
        assert err<m.mpf('1e-27'),err
        rows.append(dict(t=str(t),d=str(2*h),indices=nu,analytic=m.nstr(a,35),direct=m.nstr(b,35),relative_error=m.nstr(err,6)))
result={'precision':m.mp.dps,'exact_U_F_factorization':True,'boundary_tadpole_identity_checked':True,'series_total_degree':N,'checks':rows,'note':'m2=s=1. Both d=2.6 and d=26/7 have convergent z integrals; 1-z=u^14 smooths the latter endpoint.'}
Path(__file__).with_name('checks.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))
