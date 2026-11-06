"""Shared finite-field pole witnesses and numerator interpolation.

These routines learn candidates from numerical data. Three internal checks are
kept per fit; final independent integral samples remain mandatory.
"""
from types import SimpleNamespace
from functools import lru_cache
from reference_space import RatFunc
from flint import nmod_poly


def make_learner(p, center):
    def poly(v):return nmod_poly(v,p)
    one=poly([1]);ep=poly([center,1])*poly([(center-1)%p,1])
    def monic(v):return v*pow(int(v.leading_coefficient()),-1,p)
    def lcm(a,b):return monic((a*b)//a.gcd(b))
    interpolation_bases={}
    interpolation_checks={}
    @lru_cache(maxsize=4096)
    def denominator_values(coefficients,xs):
        den=poly(list(coefficients))
        return tuple(int(den(x)) for x in xs)
    def known_denominator_fit(den,times,values,growth=None):
        degree=len(times)-4 if growth is None else min(int(den.degree())+growth,len(times)-4)
        count=degree+1
        if count<1 or count+3>len(times):return None
        xs=tuple((t-center)%p for t in times)
        key=xs[:count]
        if key not in interpolation_bases:
            modulus=one
            for x in key:modulus*=poly([-x%p,1])
            basis=[]
            for x in key:
                b=modulus//poly([-x%p,1]);basis.append(b*pow(int(b(x)),-1,p))
            interpolation_bases[key]=basis
        check_key=(key,xs[count:])
        if check_key not in interpolation_checks:
            interpolation_checks[check_key]=[
                (k,[int(b(xs[k])) for b in interpolation_bases[key]])
                for k in range(count,len(xs))]
        den_values=denominator_values(tuple(int(c) for c in den),xs)
        y=[values[k]*den_values[k]%p for k in range(count)]
        # A wrong pole hypothesis is usually rejected by one scalar residual;
        # constructing its full interpolation polynomial is unnecessary.
        for k,weights in interpolation_checks[check_key]:
            if sum(w*v for w,v in zip(weights,y))%p!=values[k]*den_values[k]%p:return None
        numerator=poly([0])
        for k,b in enumerate(interpolation_bases[key]):
            numerator+=b*y[k]
        if any(int(numerator(x))!=v*int(den(x))%p for x,v in zip(xs,values)):return None
        return RatFunc.make([int(c) for c in numerator],[int(c) for c in den],p)
    denominator_factors={}
    def prune_denominator(den,times,values,budget=128):
        den=monic(den);key=tuple(int(c) for c in den)
        if key not in denominator_factors:
            denominator_factors[key]=sorted(den.factor()[1],key=lambda item:-int(item[0].degree()))
        factors=denominator_factors[key];state=(0,)*len(factors)
        seen={state};frontier=[(den,state)];trials=0
        for depth in range(3):
            following=[]
            for current,state in frontier:
                for j,(factor,power) in enumerate(factors):
                    if state[j]>=power:continue
                    nxt=list(state);nxt[j]+=1;nxt=tuple(nxt)
                    if nxt in seen:continue
                    seen.add(nxt);candidate=current//factor
                    fitted=known_denominator_fit(candidate,times,values)
                    trials+=1
                    if fitted is not None:return fitted
                    if trials>=budget:return None
                    following.append((candidate,nxt))
            frontier=following
        return None
    def boost_denominator(den,times,values):
        den=monic(den);key=tuple(int(c) for c in den)
        if key not in denominator_factors:
            denominator_factors[key]=sorted(den.factor()[1],key=lambda item:-int(item[0].degree()))
        for factor,_ in denominator_factors[key]:
            candidate=den
            for extra in range(2):
                candidate*=factor
                fitted=known_denominator_fit(candidate,times,values)
                if fitted is None:fitted=prune_denominator(candidate,times,values,budget=16)
                if fitted is not None:return fitted
        return None
    return SimpleNamespace(poly=poly, one=one, ep=ep, monic=monic, lcm=lcm,
                           fit=known_denominator_fit, prune=prune_denominator,
                           boost=boost_denominator)
