#!/usr/bin/env bash
set -euo pipefail
project_root=$(cd "$(dirname "$0")/.." && pwd)
test_output=$(mktemp -d /tmp/j-vac-fi-XXXXXX)
make -s -C "$project_root/native" all
"$project_root/build/j_series" "$project_root/examples/vac.json" "$test_output/series.json" --max-dot 3 --boundaries --order 80 --verify-all
"$project_root/build/j_search" "$test_output/series.json" "$test_output/j.json" --degree 6 --train 50 --de-degree 12 > "$test_output/search.log"
"$project_root/build/fi_ibp" "$project_root/examples/vac_fi_ibp.json" "$test_output/fi.json" "$test_output/j.json"
"$project_root/build/compare_blade" "$test_output/fi.json" "$project_root/tests/data/vac_blade_dot3_modp.json" "$test_output/comparison.json"
printf 'Regression artifacts: %s\n' "$test_output"
