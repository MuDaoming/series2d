#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "$0")" && pwd)"

fail=0
for d in "$ROOT"/*; do
  [[ -d "$d" ]] || continue
  [[ -f "$d/Makefile" ]] || continue
  name="$(basename "$d")"
  echo "[build] $name"
  if ! make -C "$d" -j; then
    fail=1
  fi
done

if [[ $fail -ne 0 ]]; then
  echo "build_all failed"
  exit 1
fi

echo "build_all passed"
