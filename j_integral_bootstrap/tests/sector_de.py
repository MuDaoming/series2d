"""Bottom-up unknown DE-block bootstrap, with oracle envelopes for validation.

Diagonal blocks come from reference cuts. Off-diagonal truth is used only for
envelopes and final comparison. Known lower DE is built from fitted blocks.
"""
import argparse
from fractions import Fraction
import json
from pathlib import Path
from time import perf_counter

from reference_space import (BootstrapSpec, RatFunc, RatMat, SampleFitter,
                             build_pl_ansatz, build_de_ansatz, poly_eval)
from point_oracle import NativeOracle


def components(graph):
    """Dependency-first SCCs; DE edges go from a row to its dependencies."""
    counter = 0
    stack, active, ids, low, result = [], set(), {}, {}, []
    def visit(i):
        nonlocal counter
        ids[i] = low[i] = counter
        counter += 1
        stack.append(i)
        active.add(i)
        for j in sorted(graph[i]):
            if j not in ids:
                visit(j)
                low[i] = min(low[i], low[j])
            elif j in active:
                low[i] = min(low[i], ids[j])
        if low[i] == ids[i]:
            block = []
            while True:
                j = stack.pop()
                active.remove(j)
                block.append(j)
                if j == i:
                    break
            result.append(sorted(block))
    for i in graph:
        if i not in ids:
            visit(i)
    return result


def main(args):
    ref = json.loads(Path(args.reference).read_text())
    p, masters = ref['prime'], ref['masters']
    index = {m:i for i,m in enumerate(masters)}
    zero = RatFunc.make([0], [1], p)
    matrix = []
    for m in masters:
        row = ref['de'][m]['rhs']
        matrix.append([RatFunc.make(row[n]['numerator_tau'], row[n]['denominator_tau'], p)
                       if n in row else zero for n in masters])
    graph = {i: {j for j,r in enumerate(row) if not r.is_zero() and i != j}
             for i,row in enumerate(matrix)}
    groups = components(graph)
    t0 = Fraction(ref['t0'])
    t0 = t0.numerator*pow(t0.denominator, -1, p) % p
    oracle = NativeOracle(args.binary, args.config, masters, p, t0, args.work) if args.config else None
    built, known, reports = {}, [], []
    begin = perf_counter()
    used_points = set()
    for block in groups:
        # Only the transitive dependency closure is needed in B; unrelated
        # previously solved sectors would enlarge the analytic module for no gain.
        dependencies = set()
        def follow(i):
            for j in graph[i]:
                if j not in block and j not in dependencies:
                    dependencies.add(j)
                    follow(j)
        for i in block:
            follow(i)
        sub = [i for i in known if i in dependencies]
        if any(i not in known for i in dependencies):
            raise AssertionError('SCC ordering is not dependency first')
        for i in block:
            for j in block:
                built[i,j] = matrix[i][j]  # known diagonal/cut input
        if not sub:
            reports.append(dict(top=[masters[i] for i in block], sub_rank=0,
                                diagonal_only=True))
            known.extend(block)
            continue
        truth = [[matrix[i][j] for j in sub] for i in block]
        spec = BootstrapSpec.from_inputs(
            family={'experiment': 'reference-envelope-sector-de'}, variable='tau', assignments={},
            masters=[masters[i] for i in sub+block], target=None,
            de={'B': [[built.get((i,j),zero) for j in sub] for i in sub],
                'A': [[built[i,j] for j in block] for i in block]},
            denominators=[[r.den for r in row] for row in truth],
            numerator_caps=[[-1 if r.is_zero() else len(r.num)-1 for r in row] for row in truth],
            prime=p, mode='de')
        if args.gauge:
            # A uniform scalar change S->h*S, T->h*T leaves C unchanged.
            # This is an explicit envelope-derived validation experiment,
            # not an automatic claim about analyticity or optimality.
            from analytic_bootstrap.poly import poly_mul,poly_gcd,poly_divmod,poly_der
            h=[1]
            if args.gauge_policy=='endpoints':
                h=poly_mul([t0,1],[(t0-1)%p,1],p)
            elif args.gauge_policy=='t1':
                h=[(t0-1)%p,1]
            else:
                for row in truth:
                    for r in row:
                        g=poly_gcd(h,r.den,p)
                        quotient,rem=poly_divmod(r.den,g,p)
                        if any(rem):raise AssertionError('LCM division')
                        h=poly_mul(h,quotient,p)
            shift=RatFunc.make(poly_der(h,p),h,p)*args.gauge
            for matrix_name in ('A','B'):
                mat=getattr(spec,matrix_name)
                for i in range(mat.shape[0]):
                    mat.entries[i,i]=mat.entries.get((i,i),zero)+shift
            # from_inputs compiled valuation metadata before the gauge change.
            from analytic_bootstrap.spec import compile_envelope
            compile_envelope(spec)
        start = perf_counter()
        ansatz = (build_de_ansatz(spec) if args.backend == 'jets'
                  else build_pl_ansatz(spec, backend=args.backend))
        seconds = perf_counter()-start
        fitter = SampleFitter(ansatz)
        def sample(x):
            if oracle:
                data = oracle.at(x)['DE_values']
                values = [[data[masters[i]].get(masters[j],0) for j in sub] for i in block]
            else:
                values = [[poly_eval(r.num,x,p)*pow(poly_eval(r.den,x,p),-1,p)%p
                           for r in row] for row in truth]
            return dict(prime=p,x=x,values=values)
        x=17
        print(json.dumps(dict(stage='fit-start',top=[masters[i] for i in block],
                              dimension=ansatz.dimension,gauge=args.gauge)),flush=True)
        while fitter.result()['parameters'] is None:
            try:
                fitter.add(sample(x))
                used_points.add(x)
            except ZeroDivisionError:
                pass
            x+=1
            if x>1000:
                raise RuntimeError('parameter fit not completed')
        fitted=ansatz.assemble(fitter.result()['parameters'])
        for a,i in enumerate(block):
            for b,j in enumerate(sub):
                if not (fitted[a][b]-matrix[i][j]).is_zero():
                    raise AssertionError('off-diagonal exact rational identity failed')
                built[i,j]=fitted[a][b]
        while len(fitter.holdout_points)<8:
            try:
                fitter.check_holdout(sample(x))
            except ZeroDivisionError:
                pass
            x+=1
        report=dict(top=[masters[i] for i in block], sub_rank=len(sub),
                    ambient=sum(c+1 for row in spec.envelope.numerator_caps for c in row),
                    constrained=ansatz.dimension, construction_seconds=seconds,
                    training_points=len(fitter.samples),holdout_points=8,
                    exact_rational_identity=True,status=ansatz.status)
        reports.append(report)
        known.extend(block)
        print(json.dumps(report),flush=True)
    for i in range(len(masters)):
        for j in range(len(masters)):
            if not (built.get((i,j),zero)-matrix[i][j]).is_zero():
                raise AssertionError('assembled global DE mismatch')
    result=dict(reference=args.reference,backend=args.backend,
                scope='all off-diagonal SCC blocks; oracle envelopes and diagonal blocks',
                sample_source='fresh native J-IBP' if oracle else 'reference evaluation',
                prime=p, master_count=len(masters),block_count=len(groups),blocks=reports,
                ambient=sum(r.get('ambient',0) for r in reports),
                constrained=sum(r.get('constrained',0) for r in reports),
                distinct_training_points=len(used_points),exact_global_DE_identity=True,
                total_seconds=perf_counter()-begin,
                native_calls=oracle.calls if oracle else 0,
                native_seconds=oracle.seconds if oracle else 0)
    result['gauge_power']=args.gauge
    result['gauge_policy']=args.gauge_policy
    result['masters']=masters
    result['t0']=ref['t0']
    result['de']={masters[i]: {'rhs': {
        masters[j]: {'numerator_tau':built[i,j].num,'denominator_tau':built[i,j].den}
        for j in range(len(masters)) if (i,j) in built and not built[i,j].is_zero()}}
        for i in range(len(masters))}
    Path(args.output).parent.mkdir(parents=True,exist_ok=True)
    Path(args.output).write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k not in ('blocks','de','masters')}),flush=True)


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('reference')
    parser.add_argument('output')
    parser.add_argument('--backend',choices=['jets','PL-Hasse','PL-dK'],default='PL-Hasse')
    parser.add_argument('--config')
    parser.add_argument('--binary',default='j_integral_global/build/global_point')
    parser.add_argument('--work',default='j_integral_bootstrap/results/native_points')
    parser.add_argument('--gauge',type=int,default=0)
    parser.add_argument('--gauge-policy',choices=['envelope','endpoints','t1'],default='endpoints')
    main(parser.parse_args())
