#!/bin/bash
#SBATCH --job-name=dp_J_ab_production
#SBATCH --partition=cp6
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=65536M
#SBATCH --time=0
#SBATCH --output=j_integral_bootstrap/official_route/tianhe/logs/%j-production.out
#SBATCH --error=j_integral_bootstrap/official_route/tianhe/logs/%j-production.err
set -euo pipefail
cd "$SLURM_SUBMIT_DIR"
export AB_PYTHON=/fs2/software/python/3.10_manbaforge_4.14.0-0/bin/python3
export AB_FLINT_PREFIX=/fs2/home/mdm/.local
export CXX=/usr/bin/g++
export CC=/usr/bin/gcc
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
"$AB_PYTHON" j_integral_bootstrap/official_route/run_case.py --case j_integral_bootstrap/official_route/cases/dp.json --stage production
