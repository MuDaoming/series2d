#!/usr/bin/env bash
set -euo pipefail
DIR=$(cd "$(dirname "$0")" && pwd)
"$DIR/../solve_de" "$DIR" "$DIR/f0" 20 "$DIR/solution"
