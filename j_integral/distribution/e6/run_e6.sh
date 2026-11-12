#!/bin/bash
set -euo pipefail
ROOT=$(cd "$(dirname "$0")" && pwd)
bash "$ROOT/j_integral/tests/e6/run.sh" "${1:-$(date +%Y%m%d_%H%M%S)}"
