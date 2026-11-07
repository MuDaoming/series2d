#!/bin/bash
#SBATCH --job-name=dp_J_finish_R
#SBATCH --partition=cp6
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=65536M
#SBATCH --time=0
#SBATCH --output=j_integral_bootstrap/official_route/tianhe/logs/%j-finish-R.out
#SBATCH --error=j_integral_bootstrap/official_route/tianhe/logs/%j-finish-R.err
set -euo pipefail
cd "$SLURM_SUBMIT_DIR"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
export LD_LIBRARY_PATH="/fs2/home/mdm/.local/lib:/fs2/software/gcc/11.3.0/lib64:${LD_LIBRARY_PATH:-}"
export PYTHONPATH="$SLURM_SUBMIT_DIR/j_integral_bootstrap/.deps:${PYTHONPATH:-}"
CXX=/fs2/software/gcc/11.3.0/bin/g++
PYTHON=/fs2/software/python/3.10_manbaforge_4.14.0-0/bin/python3
for NAME in preconditioned verify_reductions; do
  if [ "$NAME" = preconditioned ]; then SOURCE=j_integral_bootstrap/pole_route/preconditioned.cpp; else SOURCE=j_integral_bootstrap/native/verify_reductions.cpp; fi
  "$CXX" -O3 -std=c++17 -I/fs2/home/mdm/.local/include -Ij_integral_bootstrap/extern/include "$SOURCE" -L/fs2/home/mdm/.local/lib -Wl,-rpath,/fs2/home/mdm/.local/lib -lflint -lgmp -o "j_integral_bootstrap/build/$NAME.tmp-$SLURM_JOB_ID"
  mv "j_integral_bootstrap/build/$NAME.tmp-$SLURM_JOB_ID" "j_integral_bootstrap/build/$NAME"
done
"$PYTHON" j_integral_bootstrap/official_route/finish_global_reductions.py --work j_integral_bootstrap/results/dp_remaining_global --base j_integral_bootstrap/results/dp_remaining_global/base_verified.json --targets j_integral_bootstrap/results/dp_remaining_global/all_targets.json --verification j_integral_bootstrap/results/dp_remaining_global/fresh_manifest.json --witness-only j_integral_bootstrap/results/dp_original_warm/anchor/connection.json
