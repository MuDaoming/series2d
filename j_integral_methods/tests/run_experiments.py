#!/usr/bin/env python3
"""Run the C++ algorithms and collect correctness/timing reports; no reduction code here."""
import argparse
import json
from pathlib import Path
import shutil
import statistics
import subprocess
import time

ROOT = Path(__file__).resolve().parents[2]
HERE = ROOT / 'j_integral_methods'
BIN = HERE / 'build'
CASES = {
    'vac': ('vac.json', 100, 60, 6, 25),
    'db': ('db_topcut.json', 200, 120, 8, 30),
    'dp': ('dp_topcut.json', 200, 120, 8, 60),
}


def read(path):
    return json.loads(path.read_text())


def run(command, log, timeout=120):
    start = time.perf_counter()
    with log.open('w') as handle:
        subprocess.run([str(x) for x in command], cwd=ROOT, stdout=handle,
                       stderr=subprocess.STDOUT, check=True, timeout=timeout)
    return time.perf_counter() - start


def signature(data):
    return {'masters': sorted(data['masters']), 'reductions': data['reductions'],
            'DE': {key: value['rhs'] for key, value in data['de'].items()},
            'zeros': sorted(data['zero_objects'])}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--repeats', type=int, default=3)
    args = parser.parse_args()
    assert args.repeats > 0
    report = {'scope': 'original t,z coordinates; dot<=2; vacuum, DB cut, DP cut',
              'timing_protocol': 'sequential subprocesses; baseline once; sampled methods median of repeats; symbolic once; compilation excluded',
              'repeats': args.repeats, 'cases': {}}
    for case, (source, order, train, degree, de_degree) in CASES.items():
        work = HERE / 'results' / case
        work.mkdir(parents=True, exist_ok=True)
        family = ROOT / 'j_integral' / 'examples' / source
        print(f'{case}: Taylor baseline', flush=True)
        series, search = work / 'baseline_series.json', work / 'baseline_search.json'
        series_wall = run([BIN/'baseline_series', family, series, '--max-dot', 2,
                           '--boundaries', '--order', order, '--verify-all'], work/'baseline_series.log')
        search_wall = run([BIN/'baseline_search', series, search, '--degree', degree,
                           '--train', train, '--de-degree', de_degree], work/'baseline_search.log')
        baseline, base_series = read(search), read(series)
        result = {'input': str(family.relative_to(ROOT)),
                  'baseline': {'order': order, 'train': train,
                               'series_compute_seconds': base_series['timings']['process_compute_seconds'],
                               'series_wall_seconds': series_wall, 'search_wall_seconds': search_wall,
                               'total_wall_seconds': series_wall + search_wall,
                               'master_count': baseline['master_count'],
                               'max_de_degree': max(v['degree'] for v in baseline['de'].values())},
                  'methods': {}}
        for method in ['hermite', 'ansatz', 'symbolic']:
            print(f'{case}: {method}', flush=True)
            runs = []
            for repeat in range(1 if method == 'symbolic' else args.repeats):
                output = work / f'{method}_{repeat}.json'
                command = [BIN/('symbolic_j' if method == 'symbolic' else 'direct_j'), family, output]
                if method != 'symbolic':
                    command += ['--method', method]
                wall = run(command, work/f'{method}_{repeat}.log', timeout=90)
                data = read(output)
                assert signature(data) == signature(baseline), (case, method, 'baseline disagreement')
                runs.append({'wall_seconds': wall, 'compute_seconds': data['timings']['total_compute_seconds'],
                             'file': str(output.relative_to(HERE))})
            median = statistics.median(x['wall_seconds'] for x in runs)
            representative = min(runs, key=lambda x: abs(x['wall_seconds'] - median))
            canonical = work / f'{method}.json'
            shutil.copyfile(HERE / representative['file'], canonical)
            points = work / f'{method}_points.json'
            run([BIN/'verify_j_points', family, canonical, points], work/f'{method}_points.log')
            data, verification = read(canonical), read(points)
            assert verification['all_pass'] and len(verification['points']) == 8
            result['methods'][method] = {
                'median_wall_seconds': median,
                'median_compute_seconds': statistics.median(x['compute_seconds'] for x in runs),
                'runs': runs, 'master_count': data['master_count'], 'masters': data['masters'],
                'max_de_degree': max(v['degree'] for v in data['de'].values()),
                'DE_degrees': {k:v['degree'] for k,v in data['de'].items()},
                'sample_count': len(data['sample_points']), 'training_points': data.get('training_points', 0),
                'heldout_points': data.get('heldout_points', 0), 'Taylor_order': data['Taylor_order'],
                'diagnostics': data['diagnostics'], 'operations': data.get('operations', {}),
                'exactly_matches_Taylor_reduction_and_DE': True,
                'independent_point_checks': verification}
        # A different finite field and expansion origin; no Taylor series is generated.
        other = read(family)
        other.update(prime=1000033, t0='3/7')
        second_input = work / 'input_second.json'
        second_input.write_text(json.dumps(other, indent=2)+'\n')
        second_outputs = {}
        for method in ['hermite', 'ansatz']:
            output = work / f'{method}_second.json'
            run([BIN/'direct_j', second_input, output, '--method', method], work/f'{method}_second.log')
            data = read(output)
            assert set(data['masters']) == set(baseline['masters'])
            points = work / f'{method}_second_points.json'
            run([BIN/'verify_j_points', second_input, output, points], work/f'{method}_second_points.log')
            assert read(points)['all_pass']
            second_outputs[method] = data
        assert signature(second_outputs['hermite']) == signature(second_outputs['ansatz'])
        result['second_field'] = {'prime': 1000033, 't0': '3/7',
                                  'master_count': second_outputs['hermite']['master_count'],
                                  'max_de_degree': max(v['degree'] for v in second_outputs['hermite']['de'].values()),
                                  'backends_agree': True, 'independent_points_per_backend': 8}
        # Change the bounds of the dense IBP space, retaining all targets and their first derivatives.
        stability = work / 'ansatz_stability.json'
        run([BIN/'direct_j', family, stability, '--method', 'ansatz', '--degree-limit', 6,
             '--extra-poles', 1], work/'ansatz_stability.log')
        assert signature(read(stability)) == signature(baseline)
        result['ansatz_bound_stability'] = {'degree_limit': 6, 'extra_poles': 1, 'unchanged': True}
        report['cases'][case] = result
        (HERE/'results'/'summary.json').write_text(json.dumps(report, indent=2)+'\n')
        print(f'{case}: all methods and checks PASS', flush=True)
    print(json.dumps({k: {'baseline':v['baseline']['total_wall_seconds'],
                         **{m:d['median_wall_seconds'] for m,d in v['methods'].items()}}
                      for k,v in report['cases'].items()}, indent=2), flush=True)


if __name__ == '__main__':
    main()
