#!/bin/bash
set -euo pipefail
ROOT=$(cd "$(dirname "$0")" && pwd)
cd "$ROOT"
PYTHON=${PYTHON:-python3}
"$PYTHON" -m venv .venv
.venv/bin/python -m pip install -r requirements.txt
PYTHON="$ROOT/.venv/bin/python" bash build.sh
