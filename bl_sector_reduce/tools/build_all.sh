#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

make -C "$SCRIPT_DIR/bl_sector_reducer"
make -C "$SCRIPT_DIR/eval_reduction"
