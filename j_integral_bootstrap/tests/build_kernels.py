"""Compile the upstream PL kernels against the same FLINT as our J worker."""
import hashlib
import json
from pathlib import Path
import subprocess

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
    command = ['cc', '-O3', '-I/usr/local/include', str(source), '-o', str(OUT/name),
               '-L/usr/local/lib', '-Wl,-rpath,/usr/local/lib', '-lflint', '-lgmp']
    subprocess.run(command, check=True)
    records.append(dict(command=command, source_sha256=hashlib.sha256(source.read_bytes()).hexdigest(),
                        binary_sha256=hashlib.sha256((OUT/name).read_bytes()).hexdigest()))
(OUT/'build.json').write_text(json.dumps(records, indent=2)+'\n')
print(OUT)
