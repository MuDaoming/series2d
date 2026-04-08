#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT_DIR="$(cd "${SCRIPT_DIR}/../../.." && pwd)"
RUNNER_DIR="${ROOT_DIR}/expand/test/generator/fi_pipeline_runner"
RUNNER="${RUNNER_DIR}/fi_pipeline_runner"

echo "[1/4] Rebuilding fi_pipeline_runner..."
make -C "${RUNNER_DIR}" clean >/dev/null
make -C "${RUNNER_DIR}" >/dev/null

run_one_example() {
    local name="$1"
    local case_dir="${SCRIPT_DIR}/${name}"
    local output_file="${case_dir}/output"
    local expected_file="${case_dir}/series"

    echo "[2/4] Running example: ${name}"
    (
        cd "${case_dir}"
        "${RUNNER}" "S" "config" "target" "output"
    )

    echo "[3/4] Compare output vs expected for ${name}"
    if cmp -s "${output_file}" "${expected_file}"; then
        echo "  PASS: ${name}"
        rm -f "${output_file}"
        return 0
    else
        echo "  FAIL: ${name}"
        echo "  Keep output for inspection: ${output_file}"
        return 1
    fi
}

fail=0
run_one_example "vac" || fail=1
run_one_example "db" || fail=1

echo "[4/4] Summary"
if [[ "${fail}" -eq 0 ]]; then
    echo "All example checks passed."
    exit 0
else
    echo "Some example checks failed."
    exit 1
fi
