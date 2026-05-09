#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../../../.." && pwd)"
TOOLS="$ROOT/expand/tools"
BASE="$ROOT/expand/tests/integration"
IN="$BASE/inputs"
OUT="$BASE/results"
EXP="$BASE/expected"

families=(vac db dp_planar)
stages=(stage1_cz stage2_series2d_cachedeg0 stage3_series2d_targetdeg20 stage4_series)

build_tools() {
  for t in dump_cz dump_2dseries fi_pipeline_runner; do
    make -C "$TOOLS/$t" -j >/dev/null
  done
}

set_config_kv() {
  local file="$1" key="$2" val="$3"
  if grep -q "^${key}[[:space:]]*=" "$file"; then
    sed -i "s|^${key}[[:space:]]*=.*|${key} = ${val}|" "$file"
  else
    echo "${key} = ${val}" >> "$file"
  fi
}

normalize_config_keys() {
  local file="$1"
  if grep -q '^reduce_mode[[:space:]]*=' "$file"; then
    local v
    v=$(sed -n 's/^reduce_mode[[:space:]]*=[[:space:]]*//p' "$file" | head -n1)
    sed -i '/^reduce_mode[[:space:]]*=/d' "$file"
    set_config_kv "$file" "reduceMode" "$v"
  fi
}

prepare_family() {
  local fam="$1" work_root="$2"
  local src="$IN/$fam"
  local wf="$work_root/$fam"

  for s in "${stages[@]}"; do
    mkdir -p "$wf/$s" "$OUT/$s/$fam"
    cp -f "$src/S" "$wf/$s/S"
    cp -f "$src/config" "$wf/$s/config"
    cp -f "$src/target" "$wf/$s/target"
    normalize_config_keys "$wf/$s/config"
  done

  set_config_kv "$wf/stage2_series2d_cachedeg0/config" "deg" "0"
  set_config_kv "$wf/stage2_series2d_cachedeg0/config" "print2DMode" "cache"

  set_config_kv "$wf/stage3_series2d_targetdeg20/config" "deg" "20"
  set_config_kv "$wf/stage3_series2d_targetdeg20/config" "print2DMode" "target"

  set_config_kv "$wf/stage4_series/config" "deg" "20"
}

run_family() {
  local fam="$1" work_root="$2"
  local wf="$work_root/$fam"

  prepare_family "$fam" "$work_root"

  echo "[${fam}] stage1_cz"
  "$TOOLS/dump_cz/dump_cz" \
    "$wf/stage1_cz/S" "$wf/stage1_cz/config" "$OUT/stage1_cz/$fam/cz_topsector"

  echo "[${fam}] stage2_series2d_cachedeg0"
  "$TOOLS/dump_2dseries/dump_2dseries" \
    "$wf/stage2_series2d_cachedeg0/S" "$wf/stage2_series2d_cachedeg0/config" "$wf/stage2_series2d_cachedeg0/target" "$OUT/stage2_series2d_cachedeg0/$fam/cache_deg0.txt"

  echo "[${fam}] stage3_series2d_targetdeg20"
  "$TOOLS/dump_2dseries/dump_2dseries" \
    "$wf/stage3_series2d_targetdeg20/S" "$wf/stage3_series2d_targetdeg20/config" "$wf/stage3_series2d_targetdeg20/target" "$OUT/stage3_series2d_targetdeg20/$fam/target_deg20.txt"

  echo "[${fam}] stage4_series"
  "$TOOLS/fi_pipeline_runner/fi_pipeline_runner" \
    "$wf/stage4_series/S" "$wf/stage4_series/config" "$wf/stage4_series/target" "$OUT/stage4_series/$fam/series_deg20.txt"
}

run_all_cases() {
  local work_root
  work_root="$(mktemp -d /tmp/expand_integration_work.XXXXXX)"
  trap "rm -rf '$work_root'" EXIT

  mkdir -p "$OUT"
  build_tools
  for fam in "${families[@]}"; do
    run_family "$fam" "$work_root"
  done
  echo "done: outputs in $OUT"
}

compare_one() {
  local got="$1" exp="$2"
  diff -u "$exp" "$got" >/dev/null
}

compare_cz_after_convert() {
  local got="$1" exp="$2"
  local tg te
  tg="$(mktemp /tmp/cz_got.XXXXXX)"
  te="$(mktemp /tmp/cz_exp.XXXXXX)"
  trap "rm -f '$tg' '$te'" RETURN

  awk '/^\[after_convert:/{flag=1} flag{print}' "$got" > "$tg"
  awk '/^\[after_convert:/{flag=1} flag{print}' "$exp" > "$te"
  diff -u "$te" "$tg" >/dev/null
}

check_expected_all() {
  local fails=0
  for s in "${stages[@]}"; do
    for fam in "${families[@]}"; do
      local dgot="$OUT/$s/$fam"
      local dexp="$EXP/$s/$fam"
      if [[ ! -d "$dexp" ]]; then
        echo "[FAIL] missing expected dir: $dexp"
        fails=$((fails+1))
        continue
      fi
      for f in "$dexp"/*; do
        local bn
        bn="$(basename "$f")"
        if [[ ! -f "$dgot/$bn" ]]; then
          echo "[FAIL] missing result file: $dgot/$bn"
          fails=$((fails+1))
          continue
        fi
        if [[ "$s" == "stage1_cz" && "$bn" == "cz_topsector.txt" ]]; then
          if compare_cz_after_convert "$dgot/$bn" "$f"; then
            echo "[PASS] $s/$fam/$bn (after_convert)"
          else
            echo "[FAIL] $s/$fam/$bn (after_convert)"
            fails=$((fails+1))
          fi
          continue
        fi

        if [[ "$s" == "stage1_cz" && "$bn" == "cz_topsector_data.wl" ]]; then
          continue
        fi

        if compare_one "$dgot/$bn" "$f"; then
          echo "[PASS] $s/$fam/$bn"
        else
          echo "[FAIL] $s/$fam/$bn"
          fails=$((fails+1))
        fi
      done
    done
  done

  if [[ $fails -ne 0 ]]; then
    echo "integration check failed: $fails mismatches"
    return 1
  fi
  echo "integration check passed"
}
