"""Validate upstream bootstrap against saved J connections, never a from-scratch solve.

Reference coefficients supply envelopes and a sample oracle. They do not enter
the ansatz builder except for the known connection (or known diagonal blocks).
All reconstructed coefficients are checked as exact rational functions.
"""
import argparse
import json
import os
import sys
from fractions import Fraction
from pathlib import Path
from time import perf_counter

ROOT = Path(__file__).resolve().parents[1]
# Full upstream runtime is validated here on its nonpersistent fused branch.
# Persistent workers currently fail on DP; explicit caller settings win.
production = ROOT/'build/production'
required = ('laurent_fused_kernel', 'laurent_higher_fused_kernel',
            'coordinate_kernel.so', 'learning_rank_kernel')
if 'ANALYTIC_BOOTSTRAP_KERNEL_DIR' not in os.environ:
    selected = production if all((production/name).is_file() for name in required) else ROOT/'build/kernels'
    os.environ['ANALYTIC_BOOTSTRAP_KERNEL_DIR'] = str(selected)
if Path(os.environ['ANALYTIC_BOOTSTRAP_KERNEL_DIR']).resolve() == production.resolve():
    os.environ.setdefault('ANALYTIC_BOOTSTRAP_PERSISTENT_PURE', '0')
    os.environ.setdefault('ANALYTIC_BOOTSTRAP_EARLY_SUPPORT', '0')
sys.path[:0] = [str(ROOT / '.deps'), str(ROOT / 'extern/analytic_bootstrap-master')]
from analytic_bootstrap import (BootstrapSpec, RatFunc, RatMat, SampleFitter,
                                build_reduction_ansatz, build_de_ansatz, build_pl_ansatz)
from analytic_bootstrap.poly import poly_eval
from point_oracle import NativeOracle


def run(args):
    ref = json.loads(Path(args.reference).read_text())
    p = ref['prime']
    masters = ref['masters']
    zero = RatFunc.make([0], [1], p)

    def rf(value):
        return RatFunc.make(value['numerator_tau'], value['denominator_tau'], p)

    de_source=ref
    if args.connection:
        de_source=json.loads(Path(args.connection).read_text())
        if (de_source['prime']!=p or de_source['masters']!=masters or de_source['t0']!=ref['t0']):
            raise ValueError('connection field, basis or expansion coordinate differs')
    connection = [[rf(de_source['de'][a]['rhs'][b]) if b in de_source['de'][a]['rhs'] else zero
                   for b in masters] for a in masters]
    if args.mode == 'reduction':
        names = [a for a in ref['reductions'] if a not in masters]
        if args.limit:
            names = names[:args.limit]
        truth = [[rf(ref['reductions'][a][b]) if b in ref['reductions'][a] else zero
                  for b in masters] for a in names]
        known = connection
        output_masters = masters
    else:
        n = args.sub_rank
        if not 0 < n < len(masters):
            raise ValueError('DE test requires 0 < sub-rank < master count')
        if any(not connection[i][j].is_zero() for i in range(n) for j in range(n, len(masters))):
            raise ValueError('chosen subsector block is not closed')
        names = masters[n:]
        truth = [row[:n] for row in connection[n:]]
        known = {'B': [row[:n] for row in connection[:n]],
                 'A': [row[n:] for row in connection[n:]]}
        output_masters = masters[:n]
    # Deliberately oracle envelopes: this is a controlled correctness experiment.
    denominators = [[r.den for r in row] for row in truth]
    caps = [[-1 if r.is_zero() else len(r.num)-1 for r in row] for row in truth]
    ambient = sum(c+1 for row in caps for c in row)
    spec = BootstrapSpec.from_inputs(
        family={'reference': str(Path(args.reference).resolve()), 'experiment': 'oracle-envelope'},
        variable='tau', assignments={}, target=names if args.mode == 'reduction' else None,
        masters=masters, de=known, denominators=denominators, numerator_caps=caps,
        prime=p, mode=args.mode)
    start = perf_counter()
    if args.backend == 'jets':
        ansatz = (build_reduction_ansatz if args.mode == 'reduction' else build_de_ansatz)(spec)
    else:
        ansatz = build_pl_ansatz(spec, backend=args.backend, extra_order=args.extra_order)
    construction = perf_counter()-start
    print(json.dumps({'stage': 'ansatz', 'ambient': ambient, 'dimension': ansatz.dimension,
                      'seconds': construction}), flush=True)
    fitter = SampleFitter(ansatz)
    oracle=None
    if args.config:
        center=Fraction(ref['t0'])
        center=center.numerator*pow(center.denominator,-1,p)%p
        oracle=NativeOracle(args.binary,args.config,masters,p,center,args.work)

    def sample(x):
        if oracle:
            data=oracle.at(x)
            section=data['reduction_values'] if args.mode=='reduction' else data['DE_values']
            return {'prime':p,'x':x,'values':[[section[a].get(b,0) for b in output_masters]
                                            for a in names]}
        return {'prime': p, 'x': x, 'values': [
            [poly_eval(r.num, x, p)*pow(poly_eval(r.den, x, p), -1, p) % p for r in row]
            for row in truth]}

    start = perf_counter()
    x = 17
    while fitter.result()['parameters'] is None:
        if x > 10000:
            raise RuntimeError('fit did not determine parameters')
        try:
            fitter.add(sample(x))
        except ZeroDivisionError:
            pass
        x += 1
    fitting = perf_counter()-start
    assembled = ansatz.assemble(fitter.result()['parameters'])
    checked = 0
    for row, expected in zip(assembled, truth):
        for got, wanted in zip(row, expected):
            if not (got-wanted).is_zero():
                raise AssertionError('exact rational identity failed')
            checked += 1
    while len(fitter.holdout_points) < 8:
        try:
            fitter.check_holdout(sample(x))
        except ZeroDivisionError:
            pass
        x += 1
    result = dict(reference=args.reference, mode=args.mode, backend=args.backend,
                  extra_order=args.extra_order, prime=p, variable='tau=t-t0', t0=ref.get('t0'),
                  envelope_source='complete reference; validation only',
                  target_count=len(names), coefficient_count=checked, master_count=len(masters),
                  output_master_count=len(output_masters), ambient_parameters=ambient,
                  constrained_parameters=ansatz.dimension,
                  parameter_reduction=ambient/ansatz.dimension if ansatz.dimension else None,
                  construction_seconds=construction, fitting_seconds=fitting,
                  training_points=len(fitter.samples), holdout_points=len(fitter.holdout_points),
                  exact_rational_identity=True, ansatz_status=ansatz.status,
                  sample_source='native J-IBP' if oracle else 'reference evaluation',
                  new_native_calls=oracle.calls if oracle else 0,
                  new_native_seconds=oracle.seconds if oracle else 0,
                  diagnostics=ansatz.diagnostics)
    result['masters']=masters
    if args.mode=='reduction':
        result['de']=de_source['de']
        result['DE_source']=args.connection or args.reference
        result['reductions']={a: {b: {'numerator_tau':r.num,'denominator_tau':r.den}
                                  for b,r in zip(masters,row) if not r.is_zero()}
                              for a,row in zip(names,assembled)}
        result['tested_objects']=names
        result['zero_objects']=ref.get('zero_objects',[])
    Path(args.output).parent.mkdir(parents=True, exist_ok=True)
    Path(args.output).write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k not in ('de','reductions','masters','tested_objects','zero_objects')}), flush=True)


if __name__ == '__main__':
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('reference')
    ap.add_argument('output')
    ap.add_argument('--mode', choices=['reduction', 'de'], default='reduction')
    ap.add_argument('--sub-rank', type=int, default=1)
    ap.add_argument('--backend', choices=['jets', 'PL-Hasse', 'PL-dK'], default='PL-Hasse')
    ap.add_argument('--extra-order', type=int, default=0)
    ap.add_argument('--limit', type=int, default=0)
    ap.add_argument('--config')
    ap.add_argument('--connection',help='use a separately constructed DE instead of reference DE')
    ap.add_argument('--binary',default='j_integral_global/build/global_point')
    ap.add_argument('--work',default='j_integral_bootstrap/results/native_points')
    run(ap.parse_args())
