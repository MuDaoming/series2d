"""Isolate a rejected known DE extension; test one scalar gauge factor at a time."""
import argparse
import json
from pathlib import Path
from reference_space import BootstrapSpec,RatFunc,SampleFitter,build_de_ansatz,poly_eval
from analytic_bootstrap import Fp,FpMatrix
from analytic_bootstrap.linalg import nullspace
from analytic_bootstrap import factor_powers
from analytic_bootstrap.poly import poly_der,poly_mul


def main(args):
    ref=json.loads(Path(args.reference).read_text());p=ref['prime'];zero=RatFunc.make([0],[1],p)
    top=args.target
    if top not in ref['masters']:raise ValueError('target must be a master')
    dep=set()
    def follow(m):
        for n in ref['de'][m]['rhs']:
            if n!=top and n not in dep:
                dep.add(n);follow(n)
    follow(top)
    sub=[n for n in ref['masters'] if n in dep]
    def val(m,n):
        v=ref['de'][m]['rhs'].get(n)
        return RatFunc.make(v['numerator_tau'],v['denominator_tau'],p) if v else zero
    truth=[val(top,n) for n in sub]
    all_factors={tuple(f) for r in truth for f,k in factor_powers(r.den,p)}
    reports=[]
    for factors in [[],*[[f] for f in sorted(all_factors)],sorted(all_factors)]:
        h=[1]
        for f in factors:h=poly_mul(h,list(f),p)
        shift=RatFunc.make(poly_der(h,p),h,p)
        spec=BootstrapSpec.from_inputs(family={},variable='tau',assignments={},masters=sub+[top],target=None,
            de={'B':[[val(m,n)+(shift if m==n else zero) for n in sub] for m in sub],
                'A':[[val(top,top)+shift]]},
            denominators=[[r.den for r in truth]],
            numerator_caps=[[-1 if r.is_zero() else len(r.num)-1 for r in truth]],prime=p,mode='de')
        ansatz=build_de_ansatz(spec);fitter=SampleFitter(ansatz)
        compatible=True;error=None
        try:
            for x in range(17,80):
                values=[[poly_eval(r.num,x,p)*pow(poly_eval(r.den,x,p),-1,p)%p for r in truth]]
                fitter.add(dict(prime=p,x=x,values=values))
                if fitter.result()['parameters'] is not None:
                    got=ansatz.assemble(fitter.result()['parameters'])[0]
                    if any(not (a-b).is_zero() for a,b in zip(got,truth)):
                        raise ValueError('unique interpolant differs from exact truth')
                    break
        except ValueError as e:compatible=False;error=str(e)
        report=dict(gauge_factors=factors,dimension=ansatz.dimension,compatible=compatible,
                    unique=fitter.result()['parameters'] is not None,error=error)
        reports.append(report);print(json.dumps(report),flush=True)
    result=dict(reference=args.reference,target=top,sub_masters=sub,
                diagonal=dict(numerator=val(top,top).num,denominator=val(top,top).den),tests=reports)
    # At a linear factor with Fuchsian B, test whether a double-pole coupling
    # forces T~k/u for a regular lower-sector initial vector S0.
    residues=[]
    for factor in sorted(all_factors):
        if len(factor)!=2:continue
        u=RatFunc.make(list(factor),[1],p);alpha=(-factor[0])%p
        try:
            R=[[(u*val(m,n)).evaluate(alpha) for n in sub] for m in sub]
            c=[(u*u*r).evaluate(alpha) for r in truth]
            a=(u*val(top,top)).evaluate(alpha)
        except ZeroDivisionError:continue
        ker=nullspace(FpMatrix(Fp(p),R)).table()
        pair=[sum(x*y for x,y in zip(c,s))%p for s in ker]
        pole=[-v*pow((1+a)%p,-1,p)%p for v in pair] if (1+a)%p else None
        residues.append(dict(factor=factor,residue_B=R,diagonal_residue=a,
                             double_pole_row=c,regular_initial_kernel=ker,
                             pairings=pair,forced_simple_pole=pole))
    result['residue_diagnostics']=residues
    Path(args.output).write_text(json.dumps(result,indent=2)+'\n')


if __name__=='__main__':
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('reference');ap.add_argument('target');ap.add_argument('output')
    main(ap.parse_args())
