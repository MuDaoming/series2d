"""Compile the upstream PL kernels against the same FLINT as our J worker."""
import argparse
import hashlib
import json
import os
import sys
import tempfile
from pathlib import Path
import subprocess

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('--prefix', default=os.environ.get('FLINT_PREFIX','/usr/local'), help='FLINT installation prefix')
parser.add_argument('--production', action='store_true', help='also build the complete upstream algebraic runtime')
args = parser.parse_args()
prefix=Path(args.prefix).resolve()
ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / 'extern/analytic_bootstrap-master'
OUT = ROOT / 'build/kernels'
OUT.mkdir(parents=True, exist_ok=True)
sources = {
    'hasse_kernel': 'native/native_hasse_kernel.c',
    'derivative_kernel': 'analytic_bootstrap/backends/pl/derivative_kernel.c',
    'support_kernel': 'native/support_kernel.c',
    'conditions_kernel': 'native/conditions_kernel.c',
    'laurent_kernel': 'analytic_bootstrap/backends/pl/laurent_derivative_kernel.c',
    'laurent_higher_kernel': 'analytic_bootstrap/backends/pl/laurent_higher_kernel.c',
}
records = []
for name, relative in sources.items():
    source = SOURCE / relative
    temporary = OUT / (name + '.tmp-' + str(os.getpid()))
    command = ['cc', '-O3', '-I'+str(prefix/'include'), str(source), '-o', str(temporary),
               '-L'+str(prefix/'lib'), '-Wl,-rpath,'+str(prefix/'lib'), '-lflint', '-lgmp']
    try:
        subprocess.run(command, check=True)
        os.replace(temporary, OUT/name)
    finally:
        temporary.unlink(missing_ok=True)
    records.append(dict(command=command, source_sha256=hashlib.sha256(source.read_bytes()).hexdigest(),
                        binary_sha256=hashlib.sha256((OUT/name).read_bytes()).hexdigest()))
manifest = OUT / ('build.json.tmp-' + str(os.getpid()))
manifest.write_text(json.dumps(records, indent=2)+'\n')
os.replace(manifest, OUT/'build.json')
print(OUT)

# ProductionFitter uses its own optional incremental native libraries. Build
# them with the unchanged upstream compiler entry rather than allowing its
# silent Python fallback on production ansatzes.
env = dict(os.environ)
env['PYTHONPATH'] = os.pathsep.join([str(SOURCE), str(ROOT/'.deps'),
                                      env.get('PYTHONPATH', '')])
upstream = ROOT/'build/upstream'
upstream.mkdir(parents=True, exist_ok=True)
with tempfile.TemporaryDirectory(prefix='upstream-build-', dir=ROOT/'build') as staged:
    staged = Path(staged)
    subprocess.run([sys.executable, '-m', 'analytic_bootstrap.upstream.build_native',
                    '--prefix', str(prefix), '--output', str(staged)],
                   env=env, check=True)
    for name in ('upstream_evaluate.so', 'upstream_fit.so', 'build.json'):
        os.replace(staged/name, upstream/name)
print(upstream)

if args.production:
    subprocess.run([sys.executable, str(ROOT/'official_route/build_production_runtime.py')], env=dict(os.environ,FLINT_PREFIX=str(prefix)), check=True)
