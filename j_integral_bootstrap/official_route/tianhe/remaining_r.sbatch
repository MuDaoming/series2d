#!/bin/bash
#SBATCH --job-name=dp_J_remaining_R
#SBATCH --partition=cp6
#SBATCH --nodes=1
#SBATCH --ntasks=48
#SBATCH --cpus-per-task=1
#SBATCH --mem=131072M
#SBATCH --time=0
#SBATCH --output=j_integral_bootstrap/official_route/tianhe/logs/%j-remaining-R.out
#SBATCH --error=j_integral_bootstrap/official_route/tianhe/logs/%j-remaining-R.err
set -euo pipefail
cd "$SLURM_SUBMIT_DIR"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
export LD_LIBRARY_PATH="/fs2/home/mdm/.local/lib:${LD_LIBRARY_PATH:-}"
PYTHON=/fs2/software/python/3.10_manbaforge_4.14.0-0/bin/python3
"$PYTHON" -c 'import os; assert len(os.sched_getaffinity(0))>=48'
"$PYTHON" j_integral_bootstrap/tests/collect_points.py subset j_integral_bootstrap/results/dp_remaining_global/samples --family j_integral_global/examples/dp_physical.json --basis j_integral_global/examples/dp496_basis.json --binary j_integral_bootstrap/build/global_point_subset --targets j_integral_bootstrap/results/dp_remaining_global/targets.json --derivatives j_integral_bootstrap/results/dp_remaining_global/no_derivatives.json --first 401 --count 400 --workers 48
