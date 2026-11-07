#!/usr/bin/env python3
"""Run a family or a supplied FBI connection; output is reproducible JSON."""
from pathlib import Path
import argparse
import json
import sys
import time

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'src'))
from family import TwoLoopFamily
from series_engine import Model, generate_series
from native_engine import generate_native
from relations import search_relations, solve_rational_relations, search_master_de


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('input', type=Path)
    ap.add_argument('output', type=Path)
    ap.add_argument('--order', type=int)
    ap.add_argument('--train-order', type=int)
    ap.add_argument('--skip-flatness-check', action='store_true')
    ap.add_argument('--skip-de', action='store_true')
    ap.add_argument('--backend', choices=('native', 'reference'), default='native')
    ap.add_argument('--verify-all', action='store_true', help='Check every full polynomial integration certificate')
    ap.add_argument('--series-only', action='store_true', help='Generate coefficients without relation search')
    args = ap.parse_args()
    cfg = json.loads(args.input.read_text())
    order = args.order if args.order is not None else cfg['order']
    train = args.train_order if args.train_order is not None else cfg.get('train_order', order*2//3)
    args.output.mkdir(parents=True, exist_ok=True)
    before = time.perf_counter()
    if cfg.get('kind') == 'two_loop_family':
        fam = TwoLoopFamily(cfg, log=lambda text: print(text, flush=True))
        data = fam.connection_model()
    elif cfg.get('kind') == 'fbi_connection':
        data = cfg
    else:
        raise ValueError('kind must be two_loop_family or fbi_connection')
    (args.output/'fbi_model.json').write_text(json.dumps(data, indent=2))
    model = Model(data)
    if not args.skip_flatness_check:
        model.flatness_check()
    prep_seconds = time.perf_counter()-before
    print(f'ordinary t0={cfg["t0"]}; FBI rank={model.rank}; Q degree={model.reducer.q}', flush=True)
    result = (generate_native(model, order, args.verify_all, args.output)
              if args.backend == 'native' else
              generate_series(model, order, progress=lambda n: print(f'coefficient {n}/{order}', flush=True)))
    keys = result['period_keys']
    rows = {label: [[row.get(key, 0) for row in series] for key in keys]
            for label, series in result['coordinates'].items()}
    (args.output/'series.json').write_text(json.dumps({
        'prime':model.p,'t0':cfg['t0'],'local_variable':'tau=t-t0',
        'order':order,'period_keys':keys,'targets':rows}, indent=2))
    (args.output/'certificates.json').write_text(json.dumps(result['certificates'], indent=2))
    if args.series_only:
        (args.output/'summary.json').write_text(json.dumps({
            'name':cfg.get('name'), 'backend':args.backend, 'order':order,
            'fbi_rank':model.rank, 'period_coordinate_count':len(keys),
            'timings':{'preparation_seconds':prep_seconds,**result['timings']},
            'operations':result['operations'],
            'verification':'Native full polynomial certificates checked at n<5 and last order; all orders with --verify-all. Reference checks all.',
        }, indent=2))
        return 0
    before = time.perf_counter()
    search = search_relations(result['coordinates'], model.p, cfg.get('relation_degree', 3), train)
    solved = solve_rational_relations(list(result['coordinates']), search['relations'], model.p)
    de = {} if args.skip_de else search_master_de(result['coordinates'], solved['masters'], model.p,
                                                  cfg.get('de_degree', 6), train)
    summary = {'name':cfg.get('name'), 'input':str(args.input.resolve()), 'backend':args.backend,
               'prime':model.p,'t0':cfg['t0'],'order':order,'metadata':data.get('metadata',{}),
               'fbi_rank':model.rank,'period_coordinate_count':len(keys),
               'certificate_count':result['operations']['certificate_checks'],
               'operations':result['operations'],'timings':{
                   'preparation_seconds':prep_seconds,**result['timings'],
                   'search_seconds':time.perf_counter()-before},
               'search':search,'solution':solved,'de':de,
               'caveat':'Finite-order validated relations; period coordinates may be overcomplete. '
                        'regulated_zero endpoints require the stated analytic continuation assumptions.'}
    (args.output/'summary.json').write_text(json.dumps(summary, indent=2))
    print(json.dumps({k:summary[k] for k in ('fbi_rank','period_coordinate_count','certificate_count','operations','timings')}, indent=2))
    print(f'candidate J masters={solved["master_count"]}; rejected held-out relations={search["rejected_count"]}')
    return 0


if __name__ == '__main__':
    try:
        sys.exit(main())
    except Exception as exc:
        print(f'j_integral: {type(exc).__name__}: {exc}', file=sys.stderr)
        raise
