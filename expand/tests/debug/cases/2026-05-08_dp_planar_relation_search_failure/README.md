# 2026-05-08 dp_planar relation search failure

## Trigger
Initial symptom was on `search` side: dp_planar could not find topsector relations.

Tried standard escalation first:
- enlarge integral set
- increase expansion degree
- increase polynomial ansatz complexity in search

All still failed, so suspicion moved to upstream expansion correctness.

## Investigation Timeline
1. Validate `expand` 2D series against 2-fold-integral IBP constraints.
- Substituting the produced 2D series into IBP did not give zero.
- This strongly suggested the 2D series was wrong, but did not yet isolate where.

2. Cross-check with `reconstruct` 2D results (top-down narrowing).
- First compare target `deg=20` 2D series: mismatch.
- Then compare `deg=0` cache/reduction values: still mismatch.
- Then compare C/z around conversion:
  - `before convert`: match
  - `after convert`: mismatch

3. Locate fault stage.
- From the above chain, error is localized to convert stage (GiNaC -> FlintMod), not later recurrence/integration.
- **During C/z construction from S, very large rational numbers can appear; GiNaC handles them correctly, but the old convert path mapped GiNaC numerics through native integer (`long long`) construction for finite-field values, so out-of-range values overflowed before mod-p reduction, producing wrong `after convert` C/z. This did not show up in vac/db because their intermediate coefficients stayed below that boundary, while dp_planar crossed it.**

4. Fix and rerun full chain.
- The fix was to route GiNaC integer/rational numerators and denominators through full-precision integer-string conversion into FLINT (`fmpz_set_str`-based construction), avoiding fixed-width narrowing before mod-p reduction.
- Re-ran C/z diagnostics, deg=0 cache comparison, deg=20 target comparison, and end-to-end checks.

## Outcome
After fix:
- C/z agreed with Mathematica reference.
- deg=0 cache-level values matched reconstruct reduction values.
- dp_planar target deg=20 2D series matched reconstruct DE-solved result.
- Search side recovered: updated dp_planar (and dp) data could find topsector relations.

## Note
This file is a compact incident record. Detailed command logs and temporary debug outputs were kept in live debugging areas during the fix process.
