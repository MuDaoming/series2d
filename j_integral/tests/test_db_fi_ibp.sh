#!/usr/bin/env bash
set -euo pipefail
project_root=$(cd "$(dirname "$0")/.." && pwd)
test_output=$(mktemp -d /tmp/j-db-fi-XXXXXX)
make -s -C "$project_root/native" all
run_case() {
  local name=$1 family=$2 config=$3 reference=$4
  "$project_root/build/j_series" "$project_root/examples/$family" "$test_output/$name-series.json" --max-dot 2 --boundaries --order 80 --verify-all
  "$project_root/build/j_search" "$test_output/$name-series.json" "$test_output/$name-j.json" --degree 6 --train 50 --de-degree 24 > "$test_output/$name-search.log"
  "$project_root/build/fi_ibp" "$project_root/examples/$config" "$test_output/$name-fi.json" "$test_output/$name-j.json"
  "$project_root/build/compare_blade" "$test_output/$name-fi.json" "$project_root/tests/data/$reference" "$test_output/$name-comparison.json"
}
run_case lower db_101111_cut.json db_101111_fi_ibp.json db_101111_blade_dot2_modp.json
run_case top db_topcut.json db_fi_ibp_topcut.json db_topcut_blade_dot2_modp.json
printf 'Regression artifacts: %s\n' "$test_output"
