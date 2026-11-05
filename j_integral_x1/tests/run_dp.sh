#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")/../.."
make -C j_integral_x1/native all test
mkdir -p j_integral_x1/results/dp_topcut
for sample in primary second; do
  input=j_integral_x1/examples/dp_topcut.json
  order=200
  train=120
  suffix=
  if [[ "$sample" == second ]]; then
    input=j_integral_x1/examples/dp_topcut_second.json
    order=220
    train=140
    suffix=_second
  fi
  j_integral_x1/build/j_series "$input" \
    "j_integral_x1/results/dp_topcut/series${suffix}.json" \
    --max-dot 2 --boundaries --order "$order" --verify-all
  j_integral_x1/build/j_search \
    "j_integral_x1/results/dp_topcut/series${suffix}.json" \
    "j_integral_x1/results/dp_topcut/search${suffix}.json" \
    --degree 8 --train "$train" --de-degree 60
  j_integral_x1/build/verify_j_points "$input" \
    "j_integral_x1/results/dp_topcut/search${suffix}.json" \
    "j_integral_x1/results/dp_topcut/points${suffix}.json"
done
