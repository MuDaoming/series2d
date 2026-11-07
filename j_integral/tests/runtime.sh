#!/bin/bash
# Select the unchanged AB source package and its already-built native runtime.
J_TEST_ROOT=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
J_AB_ROOT="$J_TEST_ROOT/../archive/j_integral_bootstrap"
export PYTHONPATH="$J_AB_ROOT/.deps:$J_AB_ROOT/extern/analytic_bootstrap-master${PYTHONPATH:+:$PYTHONPATH}"
export ANALYTIC_BOOTSTRAP_KERNEL_DIR="$J_AB_ROOT/build/production"
export ANALYTIC_BOOTSTRAP_UPSTREAM_KERNEL_DIR="$J_AB_ROOT/build/upstream"
# Use the package's previously validated nonpersistent native implementation.
export ANALYTIC_BOOTSTRAP_PERSISTENT_PURE=0
export ANALYTIC_BOOTSTRAP_EARLY_SUPPORT=0
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
export ANALYTIC_BOOTSTRAP_FIT_THREADS=1
