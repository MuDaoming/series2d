#!/bin/bash
set -euo pipefail
ROOT=$(cd "$(dirname "$0")" && pwd)
cd "$ROOT"
PYTHON=${PYTHON:-python3}
FLINT_PREFIX=${FLINT_PREFIX:-/usr/local}
if [ ! -f "$FLINT_PREFIX/include/flint/flint.h" ]; then
  if [ -f /usr/include/flint/flint.h ] && [ "$FLINT_PREFIX" = /usr/local ]; then FLINT_PREFIX=/usr; else
    echo "FLINT headers not found; set FLINT_PREFIX to your installation prefix." >&2; exit 1
  fi
fi
export FLINT_PREFIX
export PYTHONPATH="$ROOT/vendor/analytic_bootstrap${PYTHONPATH:+:$PYTHONPATH}"
"$PYTHON" -c 'import sys,flint; assert sys.version_info >= (3,10); print(sys.version.split()[0],"python-flint",flint.__version__)'
for TOOL in make cc c++ pkg-config; do command -v "$TOOL" >/dev/null; done
mkdir -p build/pkgconfig
"$PYTHON" - <<'PY'
from pathlib import Path
import os,re
prefix=Path(os.environ['FLINT_PREFIX']);header=prefix/'include/flint/flint.h'
version=re.search(r'#define\s+FLINT_VERSION\s+"([^"]+)"',header.read_text()).group(1)
Path('build/pkgconfig/flint.pc').write_text(f'prefix={prefix}\nlibdir=${{prefix}}/lib\nincludedir=${{prefix}}/include\nName: FLINT\nDescription: Existing FLINT installation\nVersion: {version}\nLibs: -L${{libdir}} -lflint -lgmp -lmpfr\nCflags: -I${{includedir}}\n')
print('Native FLINT:',version,'prefix:',prefix)
PY
export PKG_CONFIG_PATH="$ROOT/build/pkgconfig${PKG_CONFIG_PATH:+:$PKG_CONFIG_PATH}"
make -C j_integral/src/IBP FLINT_PREFIX="$FLINT_PREFIX"
AB="$ROOT/vendor/analytic_bootstrap"
for OPTION in '' --hasse --derivative --support --support-legacy --conditions; do
  if [ -n "$OPTION" ]; then "$PYTHON" "$AB/tools/build_native_pole_kernel.py" "$OPTION"; else "$PYTHON" "$AB/tools/build_native_pole_kernel.py"; fi
done
for OPTION in '' --higher --fused --higher-fused; do
  if [ -n "$OPTION" ]; then "$PYTHON" "$AB/tools/build_laurent_kernel.py" "$OPTION"; else "$PYTHON" "$AB/tools/build_laurent_kernel.py"; fi
done
"$PYTHON" "$AB/tools/build_coordinate_kernel.py"
"$PYTHON" "$AB/tools/build_pure_runtime.py"
"$PYTHON" -m analytic_bootstrap.upstream.build_native --prefix "$FLINT_PREFIX" --output "$ROOT/build/upstream"
echo 'Build complete. Run: bash run_e6.sh RUN_NAME'
