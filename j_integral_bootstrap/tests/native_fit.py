"""Joint ambient fit, following upstream sampling_batch's native sparse route.

Known rational entries are added as affine equations without expanding the
global-cap kernel at every evaluation point.
"""
from time import perf_counter
from pathlib import Path
import sys
from reference_space import SampleFitter
from analytic_bootstrap.sampling import ambient_evaluation_rows
from analytic_bootstrap.ansatz import _native_conditions_eliminate
from analytic_bootstrap.poly import poly_divmod,poly_mul,poly_degree


def joint_fit(ansatz,samples,fixed=None,identity_points=None):
    fixed=fixed or {}
    # Identity certificates need cap+1 distinct elements of Fp, not distinct
    # integer representatives (for example, 0 and p are the same node).
    identity_points=list(dict.fromkeys(int(x)%ansatz.spec.prime for x in (identity_points or [])))
    samples=list(samples)
    root=Path(__file__).resolve().parents[1]
    kernels=root/'build/upstream'
    precheck=None
    if kernels.is_dir() and ansatz.dimension>=512:
        if str(root) not in sys.path:sys.path.insert(0,str(root))
        from official_route.preflight import preflight
        precheck=preflight(ansatz,samples,max_points=32,kernel_directory=kernels)
        if precheck['status']=='incompatible':
            return dict(compatible=False,dimension=ansatz.dimension,rank=precheck.get('rank',0),
                        parameters=None,preflight=precheck,sample_count=precheck['used_samples'],
                        fixed_entries=len(fixed),total_seconds=precheck['seconds'])
    p=ansatz.spec.prime;n=ansatz.ambient_parameters;nc=ansatz.shape[1];started=perf_counter()
    equations={i:dict(row) for i,(_,row) in enumerate(sorted(ansatz.conditions.pivots.items()))}
    cap_rank=len(equations);target_masks=[0]*n;degrees=[0]*n
    def equation(row,rhs,output):
        row=dict(row)
        for c in row:target_masks[c]|=1<<(output//nc);degrees[c]+=1
        if rhs:row[n]=rhs
        equations[len(equations)]=row
    counts={k:0 for k in fixed}
    for (a,j),r in fixed.items():
        quotient,remainder=poly_divmod(ansatz.spec.envelope.denominators[a][j],r.den,p)
        if any(remainder) or poly_degree(poly_mul(r.num,quotient,p))>ansatz.spec.envelope.numerator_caps[a][j]:
            raise ValueError('known function violates its input envelope')
    for x in identity_points:
        needed=[k for k in fixed if counts[k]<=ansatz.spec.envelope.numerator_caps[k[0]][k[1]]]
        if not needed:break
        try:rows=ambient_evaluation_rows(ansatz,x)
        except ZeroDivisionError:continue
        for a,j in needed:
            try:rhs=fixed[a,j].evaluate(x)
            except ZeroDivisionError:continue
            equation(rows[a*nc+j],rhs,a*nc+j);counts[a,j]+=1
    if any(counts[k]<=ansatz.spec.envelope.numerator_caps[k[0]][k[1]] for k in fixed):
        raise ValueError('insufficient points for a fixed rational identity')
    identity_rows=len(equations)-cap_rank;checker=SampleFitter(ansatz);unique={}
    for sample in samples:
        x,ys=checker._sample(sample)
        if x in unique:
            if unique[x]!=ys:raise ValueError('conflicting sample')
            continue
        unique[x]=ys
        for output,(row,rhs) in enumerate(zip(ambient_evaluation_rows(ansatz,x),ys)):
            equation(row,rhs,output)
    order=sorted(range(n),key=lambda c:(target_masks[c].bit_count(),
        target_masks[c] if target_masks[c].bit_count()==1 else 0,degrees[c],c))
    inverse={c:j for j,c in enumerate(order)};inverse[n]=n
    row_order=sorted(range(cap_rank,len(equations)),key=lambda i:len(equations[i]))
    row_order+=sorted(range(cap_rank),key=lambda i:len(equations[i]))
    permuted={i:{inverse[c]:v for c,v in equations[r].items()} for i,r in enumerate(row_order)}
    assembled=perf_counter();pivots=_native_conditions_eliminate(p,n+1,permuted)
    if pivots is None:raise RuntimeError('native joint fitter unavailable')
    del permuted
    rank=len(pivots)-int(n in pivots)-cap_rank
    result=dict(compatible=n not in pivots,dimension=ansatz.dimension,rank=rank,
                parameters=None,sample_count=len(unique),fixed_entries=len(fixed),
                fixed_equation_count=identity_rows,equation_rows=len(equations),
                assembly_seconds=assembled-started,elimination_seconds=perf_counter()-assembled)
    if precheck is not None:result["preflight"]=precheck
    if n in pivots:return result
    permuted_z=[0]*n
    for c,row in sorted(pivots.items(),reverse=True):
        permuted_z[c]=(row.get(n,0)-sum(v*permuted_z[j] for j,v in row.items() if j!=c and j!=n))%p
    z=[permuted_z[inverse[c]] for c in range(n)]
    if any((sum(v*z[j] for j,v in row.items() if j!=n)-row.get(n,0))%p for row in equations.values()):
        raise ArithmeticError('joint fit original-equation witness failed')
    parameters=[z[c] for c in ansatz.conditions.free_columns]
    if ansatz.local_parameters(parameters)!=z:raise ArithmeticError('joint fit cap roundtrip failed')
    for x,ys in unique.items():
        if [v for row in ansatz.evaluate(x,parameters) for v in row]!=ys:raise ArithmeticError('joint fit independent evaluation failed')
    if rank==ansatz.dimension:result['parameters']=parameters
    result['total_seconds']=perf_counter()-started
    return result
