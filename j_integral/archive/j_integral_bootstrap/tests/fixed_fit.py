"""Fix already reconstructed entries exactly before asking IBP for unknowns."""
from reference_space import SampleFitter
from analytic_bootstrap.sampling import evaluation_rows,InconsistentSampleError
from analytic_bootstrap.ansatz import _native_conditions_eliminate
from analytic_bootstrap.poly import poly_divmod,poly_mul,poly_degree


def fixed_entry_fitter(ansatz,entries,points):
    """entries maps (row,column) to a known rational function.

    Entry numerators have the exact input caps. cap+1 ordinary evaluations
    therefore impose a rational identity, rather than extra IBP samples.
    """
    fitter=SampleFitter(ansatz)
    if not entries:return fitter,0
    p=ansatz.spec.prime;nc=ansatz.shape[1];equations={};counts={k:0 for k in entries}
    for (a,j),r in entries.items():
        quotient,remainder=poly_divmod(ansatz.spec.envelope.denominators[a][j],r.den,p)
        if any(remainder) or poly_degree(poly_mul(r.num,quotient,p))>ansatz.spec.envelope.numerator_caps[a][j]:
            raise ValueError('known entry does not obey its denominator and numerator cap')
    for x in points:
        needed=[k for k in entries if counts[k]<=ansatz.spec.envelope.numerator_caps[k[0]][k[1]]]
        if not needed:break
        try:rows=evaluation_rows(ansatz,x)
        except ZeroDivisionError:continue
        for a,j in needed:
            r=entries[a,j]
            try:rhs=r.evaluate(x)
            except ZeroDivisionError:continue
            equation=dict(rows[a*nc+j])
            if rhs:equation[fitter.n]=rhs
            equations[len(equations)]=equation;counts[a,j]+=1
    if any(counts[k]<=ansatz.spec.envelope.numerator_caps[k[0]][k[1]] for k in entries):
        raise ValueError('not enough ordinary points to impose a known entry identity')
    pivots=_native_conditions_eliminate(p,fitter.n+1,equations)
    if pivots is None:raise RuntimeError('native fixed-entry equation solver unavailable')
    if fitter.n in pivots:raise InconsistentSampleError('known entries incompatible with analytic ansatz')
    fitter.pivots={c:({j:v for j,v in row.items() if j!=fitter.n},row.get(fitter.n,0))
                   for c,row in pivots.items()}
    # Check an explicit affine witness for every original equation.
    values=[0]*fitter.n
    for c,(row,rhs) in sorted(fitter.pivots.items(),reverse=True):
        values[c]=(rhs-sum(v*values[j] for j,v in row.items() if j!=c))%p
    if any((sum(v*values[j] for j,v in row.items() if j!=fitter.n)-row.get(fitter.n,0))%p
           for row in equations.values()):raise ArithmeticError('fixed-entry witness residual')
    return fitter,len(pivots)
