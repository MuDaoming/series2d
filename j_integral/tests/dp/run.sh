#!/bin/bash
set -euo pipefail
CASE_DIR=$(cd "$(dirname "$0")" && pwd)
source "$CASE_DIR/../runtime.sh"
RUN_DIR="$CASE_DIR/runs/${1:-$(date +%Y%m%d_%H%M%S)}"
python "$CASE_DIR/../prepare_run.py" --case "$CASE_DIR" --run "$RUN_DIR"
python -m analytic_bootstrap.upstream learn \
  --config "$RUN_DIR/backend.json" --definition "$CASE_DIR/input/definition.json" \
  --setting "$CASE_DIR/input/setting.json" --output "$RUN_DIR/bundle.json" \
  --timings "$RUN_DIR/learn_timings.json" \
  >"$RUN_DIR/learn.stdout.log" 2>"$RUN_DIR/learn.stderr.log"
python -m analytic_bootstrap.upstream solve-de \
  --config "$RUN_DIR/backend.json" --bundle "$RUN_DIR/bundle.json" \
  --setting "$CASE_DIR/input/setting.json" --output "$RUN_DIR/solution.json" \
  --timings "$RUN_DIR/solve_timings.json" \
  >"$RUN_DIR/solve.stdout.log" 2>"$RUN_DIR/solve.stderr.log"
python "$CASE_DIR/../verify_result.py" --case "$CASE_DIR" --run "$RUN_DIR" >"$RUN_DIR/verification.log" 2>&1
cat "$RUN_DIR/verification.log"
python "$CASE_DIR/../export_physical.py" --case "$CASE_DIR" --run "$RUN_DIR" >"$RUN_DIR/export.log" 2>&1
