# FI/BFI/BBFI Delta Series Implementation Plan

## Goal

Extend `expand` so the same pipeline can output one-dimensional delta series for:

- `FI{nu}`: two integrations
- `BFI[boundary]{nu}`: one boundary evaluation and one integration
- `BBFI[boundary,boundary]{nu}`: two boundary evaluations

For a two-dimensional expansion through total degree `deg`, all outputs should contain the determined coefficients from `delta^0` through `delta^deg`. Do not omit the leading zeros for FI or BFI, and do not pad unknown higher orders with zeros.

## Main Design

1. Add an `IntegralTag` layer.

   `nu` should remain only the propagator-index vector. The head (`FI`, `BFI`, `BBFI`) and boundary information (`XU`, `XD`, `YU`, `YD`) should live outside `nu`.

   This avoids the experimental synthetic-index problem where a tag like `{4,1,1,1}` can be accidentally sorted or reduced as if `4` were a dot.

2. Reuse the same two-dimensional integrand series.

   For each `nu`, compute:

   ```text
   G_nu(X,Y) = P(X,Y) * I_tilde_nu(X,Y)
   ```

   once, then project it differently for FI/BFI/BBFI.

3. Add a `DeltaProjector`.

   The projector owns the precomputed powers and integration weights:

   ```text
   (-a)^k, (1-a)^k, (-b)^k, (1-b)^k
   intX[k] = ((1-a)^(k+1) - (-a)^(k+1))/(k+1)
   intY[k] = ((1-b)^(k+1) - (-b)^(k+1))/(k+1)
   ```

   Projection is then one pass over the nonzero 2D coefficients.

## File-Level Changes

1. `expand/include/integral_tag.hpp`

   Add:

   - `IntegralHead`
   - `BoundaryAxis`
   - `BoundarySide`
   - `BoundaryTag`
   - `IntegralTag`
   - helpers for parsing/formatting if they are small and header-only

2. `expand/include/io.hpp` and `expand/src/io.tpp`

   Change `TargetConfig` to:

   ```cpp
   struct TargetConfig {
       std::vector<IntegralTag> targets;
   };
   ```

   Keep backward compatibility:

   ```text
   {1,1,1}
   ```

   parses as:

   ```text
   FI{1,1,1}
   ```

   New explicit formats:

   ```text
   FI{1,1,1}
   BFI[XU]{1,1,1}
   BBFI[XU,YD]{1,1,1}
   ```

3. `expand/include/integrand_expander.hpp` and `expand/src/integrand_expander.tpp`

   Add:

   ```cpp
   Series<ST> getIntegrand2DSeries(const std::vector<int>& nu) const;
   ```

4. `expand/include/delta_projector.hpp` and `expand/src/delta_projector.tpp`

   New module for the 2D-to-1D projection.

   It should validate:

   - `FI`: no boundaries
   - `BFI`: exactly one boundary
   - `BBFI`: exactly two boundaries, one `X*` and one `Y*`

   Invalid tags should throw clear errors.

5. `expand/src/fi_pipeline.tpp`

   Replace the target loop:

   ```cpp
   for (const auto& nu : targetCfg.nus)
   ```

   with:

   ```cpp
   for (const auto& tag : targetCfg.targets)
   ```

   Use a local cache keyed by `tag.nu`:

   ```cpp
   map<vector<int>, Series<FlintMod>> integrand2DCache;
   ```

   Then call:

   ```cpp
   projector.project(g2d, tag)
   ```

## Efficiency Notes

- Do not recompute `getFBISeries` or `Series::mulPoly` for multiple heads with the same `nu`.
- Precompute all boundary powers and integral weights needed for 2D coefficients through total degree `deg`. FI and BFI need powers up to `deg + 1` inside the integration weights; BBFI uses direct powers up to `deg`.
- Skip zero coefficients during projection.
- Keep projection `O(deg^2)`. The expensive part remains solving FBI and multiplying by `P`.

## Output Convention

Every output line is:

```text
{c0,c1,...,c_deg}
```

No implicit deletion of leading zeros:

- FI writes `c0=0` and `c1=0` explicitly.
- BFI writes `c0=0` explicitly when it vanishes.
- BBFI starts naturally at `delta^0`.
- Coefficients beyond `delta^deg` are not generated from a 2D expansion of degree `deg`, and are not padded with zeros.

This is important for later mixed FI/BFI/BBFI searches.

## Validation Plan

1. Backward compatibility:

   Run an existing FI-only target file. The values from `delta^2` onward should match the old output, with the old omitted zeros now present if the main runner is switched to the new convention.

2. Vac sanity checks:

   For `FI{1,1,1}`, verify:

   ```text
   d/delta FI111
   = (1-a) BFI[XU]111 + a BFI[XD]111
     + (1-b) BFI[YU]111 + b BFI[YD]111
   ```

   using the explicit full series convention.

3. BFI/BBFI relation checks:

   Reproduce the experimental lowest-degree search pattern:

   - `dFI111` relation at `m=0`
   - `dBFI` relations at `m=2`
   - `dBBFI` relations at `m=4`

4. Deg consistency:

   Check that all output series have exactly `deg + 1` entries.

## Open Implementation Points

- Decide whether `TargetConfig` keeps a deprecated `nus` field during migration. Prefer updating all in-repo callers to `targets` and keeping compatibility only in the parser syntax.
- Decide whether target output should include labels. Current series files are line-aligned with target order; keeping that format is simplest for `search`.
- Search-side support for head-aware objects is separate. This plan only updates `expand`; `search` should later receive a matching tag parser instead of synthetic leading indices.
