#!/bin/bash
E6_BUNDLE_ROOT=$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)
if [ -x "$E6_BUNDLE_ROOT/.venv/bin/python" ]; then
  export PATH="$E6_BUNDLE_ROOT/.venv/bin:$PATH"
  export E6_PYTHON="$E6_BUNDLE_ROOT/.venv/bin/python"
else
  export E6_PYTHON="${PYTHON:-python3}"
fi
export PYTHONPATH="$E6_BUNDLE_ROOT/vendor/analytic_bootstrap${PYTHONPATH:+:$PYTHONPATH}"
export ANALYTIC_BOOTSTRAP_KERNEL_DIR="$E6_BUNDLE_ROOT/vendor/analytic_bootstrap/build/hjj_native"
export ANALYTIC_BOOTSTRAP_UPSTREAM_KERNEL_DIR="$E6_BUNDLE_ROOT/build/upstream"
export ANALYTIC_BOOTSTRAP_PERSISTENT_PURE=0
export ANALYTIC_BOOTSTRAP_EARLY_SUPPORT=0
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
export ANALYTIC_BOOTSTRAP_FIT_THREADS=1
