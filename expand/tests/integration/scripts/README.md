# integration/scripts

- `common.sh`: shared functions (build tools, prepare staged configs, run cases, compare expected).
- `run_cases.sh`: run all integration cases and generate `../results`.
- `check_expected.sh`: compare `../results` against `../expected`.
- `test.sh`: run + check (CI-style entry).

## What Each Stage Checks

- `stage1_cz`:
  - compare `cz_topsector.txt` only on the `[after_convert: ...]` block (`diff`).
  - reason: `before_convert` term order is not stable across runs.

- `stage2_series2d_cachedeg0`:
  - compare `cache_deg0.txt` by `diff`.
  - compare `cache_deg0.txt.keys` by `diff`.

- `stage3_series2d_targetdeg20`:
  - compare `target_deg20.txt` by `diff`.

- `stage4_series`:
  - compare `series_deg20.txt` by `diff`.
  - checks FI/BFI/BBFI delta-series projection. For 2D expansion degree `deg`,
    each line contains the determined coefficients of `delta^0` through
    `delta^deg`.

Typical usage:

```bash
./test.sh
```
