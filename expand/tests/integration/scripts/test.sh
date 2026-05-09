#!/usr/bin/env bash
set -euo pipefail
"$(dirname "$0")/run_cases.sh"
"$(dirname "$0")/check_expected.sh"
