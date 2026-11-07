#!/bin/bash
#SBATCH --job-name=dp_J_ab_warm
#SBATCH --partition=cp6
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=65536M
#SBATCH --time=0
#SBATCH --output=j_integral_bootstrap/official_route/tianhe/logs/%j-warm.out
#SBATCH --error=j_integral_bootstrap/official_route/tianhe/logs/%j-warm.err
set -euo pipefail
cd "$SLURM_SUBMIT_DIR"
export AB_PYTHON=/fs2/software/python/3.10_manbaforge_4.14.0-0/bin/python3
export AB_FLINT_PREFIX=/fs2/home/mdm/.local
export AB_GCC_PREFIX=/fs2/software/gcc/11.3.0
export CXX="$AB_GCC_PREFIX/bin/g++"
export CC="$AB_GCC_PREFIX/bin/gcc"
# Upstream builders literally invoke cc; the shared GCC prefix has gcc only.
mkdir -p j_integral_bootstrap/build/toolchain/bin
ln -sfn "$CC" j_integral_bootstrap/build/toolchain/bin/cc
export PATH="$SLURM_SUBMIT_DIR/j_integral_bootstrap/build/toolchain/bin:$AB_GCC_PREFIX/bin:$PATH"
export LD_LIBRARY_PATH="$AB_GCC_PREFIX/lib64:${LD_LIBRARY_PATH:-}"
export LD_LIBRARY_PATH="$AB_FLINT_PREFIX/lib:${LD_LIBRARY_PATH:-}"
export PYTHONPATH="$SLURM_SUBMIT_DIR/j_integral_bootstrap/.deps:${PYTHONPATH:-}"
export ANALYTIC_BOOTSTRAP_PERSISTENT_PURE=0
export ANALYTIC_BOOTSTRAP_EARLY_SUPPORT=0
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1
export ANALYTIC_BOOTSTRAP_FIT_THREADS=1
: "${AB_PYTHON:?Set the verified private Python >=3.10 executable before submission}"
"$AB_PYTHON" -c 'import sys,flint; assert sys.version_info >= (3,10); print(sys.executable,flint.__version__)'
"$AB_PYTHON" j_integral_bootstrap/official_route/compiler_probe.py --prefix "$AB_FLINT_PREFIX"
"$AB_PYTHON" j_integral_bootstrap/official_route/build_backend.py --prefix "$AB_FLINT_PREFIX" --json-include j_integral_bootstrap/extern/include
"$AB_PYTHON" j_integral_bootstrap/tests/build_kernels.py --prefix "$AB_FLINT_PREFIX" --production
"$AB_PYTHON" j_integral_bootstrap/official_route/run_case.py --case j_integral_bootstrap/official_route/cases/dp_warm.json --stage production
