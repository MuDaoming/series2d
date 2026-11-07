"""FLINT polynomial-LRR backend for the known FBI model, not a supplied J DE."""
from pathlib import Path
import json
import subprocess
import tempfile
import time

ROOT = Path(__file__).resolve().parents[1]


def polynomial_coefficients(model, value):
    p = model.reducer.polynomial(value)
    return [int(p.get((j,), 0)) % model.p for j in range(p.degree()+1)]


def export_model(model, order, verify_all=False):
    poly = lambda x: polynomial_coefficients(model, x)
    jet = lambda x: [poly(p) for p in x]
    return {'prime': model.p, 'rank': model.rank, 'order': order,
            'verify_all': verify_all,
            'q': poly(model.reducer.Qrat),
            'B': [[poly(model.reducer.K.field.new(p)) for p in row]
                  for row in model.reducer.N],
            'qt': jet(model.qt), 'nt': [[jet(x) for x in row] for row in model.nt],
            'endpoints': [model.reducer.endpoints[e] for e in (0, 1)],
            'targets': [{'label': label, 'den': jet(d), 'numer': [jet(x) for x in b]}
                        for label, (d, b) in model.object_data.items()]}


def generate_native(model, order, verify_all=False, output_dir=None):
    if not 0 <= order+1 < model.p:
        raise ValueError('Require 0 <= order+1 < prime')
    subprocess.run(['make', '-s', '-C', str(ROOT/'native')], check=True)
    payload = export_model(model, order, verify_all)
    def run(directory):
        directory = Path(directory)
        src, dst = directory/'native_input.json', directory/'native_output.json'
        src.write_text(json.dumps(payload))
        before = time.perf_counter()
        subprocess.run([str(ROOT/'build/j_series'), str(src), str(dst)], check=True)
        result = json.loads(dst.read_text())
        elapsed = time.perf_counter()-before
        coords = {label: [{tuple(int(v) if i else v for i, v in enumerate(k.split(':'))): c
                            for k, c in row.items()} for row in series]
                  for label, series in result['coordinates'].items()}
        keys = sorted({k for series in coords.values() for row in series for k in row}, key=str)
        result.update(coordinates=coords, period_keys=keys, certificates={})
        result['timings']['native_process_seconds'] = elapsed
        return result
    if output_dir is not None:
        return run(output_dir)
    with tempfile.TemporaryDirectory(prefix='j-series-') as directory:
        return run(directory)
