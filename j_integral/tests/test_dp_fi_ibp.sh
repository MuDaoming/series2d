#!/usr/bin/env bash
set -euo pipefail
project_root=$(cd "$(dirname "$0")/.." && pwd)
test_output=$(mktemp -d /tmp/j-dp-fi-XXXXXX)
make -s -C "$project_root/native" all ../build/verify_j_points
"$project_root/build/j_series" "$project_root/examples/dp_topcut.json" "$test_output/series.json" --max-dot 2 --boundaries --order 120 --verify-all
"$project_root/build/j_search" "$test_output/series.json" "$test_output/j.json" --degree 8 --train 75 --de-degree 50 > "$test_output/search.log"
"$project_root/build/verify_j_points" "$project_root/examples/dp_topcut.json" "$test_output/j.json" "$test_output/points.json"
"$project_root/build/fi_ibp" "$project_root/examples/dp_fi_ibp_blade_basis.json" "$test_output/fi.json" "$test_output/j.json" "$project_root/tests/data/dp_blade_master_info.json"
"$project_root/build/compare_blade" "$test_output/fi.json" "$project_root/tests/data/dp_blade_dot2_cut_modp.json" "$test_output/comparison.json"
printf 'Regression artifacts: %s\n' "$test_output"
