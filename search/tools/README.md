# search/tools

This directory contains standalone command-line tools built on top of `search` core code.

## Build

Each tool has its own `Makefile`.

```bash
cd search/tools/<tool_name>
make
```

## Common Input Files

### `G`
- Search target integral list.
- One integral tag per line, e.g. `FI{1,1,1,1}`, `BFI[XU]{1,1,1,1}`,
  or `BBFI[XU,YD]{1,1,1,1}`.
- A bare vector such as `{1,1,1,1}` is accepted as legacy input and means
  `FI{1,1,1,1}`.

### `config`
(Used by `poly_relation_searcher` and `search_de`)
- Contains finite-field and search settings.
- Typical keys include:
  - `p`
  - `deg`
  - `m`
  - `ncheck` (optional, default `1`)

### `series`
(Used by `poly_relation_searcher`)
- 1D expansion input for integrals in `G`.

### `poly_relation`
(Used by `integral_solver`)
- Output of `poly_relation_searcher`.
- Must contain header fields and `[rref]` block.

## Tools

### 1) `poly_relation_searcher`

Path: `search/tools/poly_relation_searcher`

Purpose:
- Build and solve polynomial-relation system from `G` and 1D delta expansion data.
- Output relation matrix / RREF for next-stage integral reduction solving.

Usage:

```bash
./poly_relation_searcher <config_path> <G_path> <series_list_path> <output_path>
```

Output:
- `<output_path>` with metadata header and `[rref]` block.

### 2) `integral_solver`

Path: `search/tools/integral_solver`

Purpose:
- Consume `G` and stage-I `poly_relation` result.
- Evaluate polynomial relations at a chosen delta value and build reductions
  among the integral tags in `G`.

Usage:

```bash
./integral_solver <G_path> <poly_relation_path> <delta_value> <output_path>
```

Output order in `output_path`:
1. `# p`, `# m`, `# delta`, `# |G|`
2. `# integral variables`, `# integral pivot columns`, `# integral free columns`
3. `#MIs`
4. `[reductions]`
5. `[relations]`
6. `[integral_rref]`

### 3) `search_de`

Path: `search/tools/search_de`

Purpose:
- For each master integral `M_i`, search a polynomial relation between
  `dM_i/delta` and the full master set.
- The derivative object is represented internally by shifting the same
  FI/BFI/BBFI tag by `nu -> nu + {shift,...,shift}`. Output maps it back to
  `d(M_i)`.
- The tool does not write intermediate target, series, or poly-relation files.

Usage:

```bash
./search_de <config_path> <G_path> <series_list_path> <masters_path> <output_path> [shift]
```

- `m` is read from `config_path` and used as the maximum degree to try.
- `search_de` tries polynomial degrees with an exponential schedule:
  `0, 1, 2, 4, 8, ...`, with a final `m` attempt if the schedule would
  otherwise skip the configured maximum. It stops at the first certified
  relation and does not search for the minimal working degree.
- `ncheck` is read from `config_path` with default `1`. Since `dM` is formed
  from `M` by differentiating the delta series, the effective derivative degree
  is `deg - 1`; `search_de` solves on `0..deg-1-ncheck` and requires the
  candidate relation to pass the held-out degrees before writing it as
  certified. The output also records whether the full train nullspace passed
  the held-out non-shrink check.
- Default `shift`: `100`.
- Output is written exactly to `output_path`.

### 4) `de_to_wl.wl`

Path: `search/tools/search_de/de_to_wl.wl`

Purpose:
- Convert `search_de` output into a Mathematica/Wolfram Language file.
- It reads `ADelta` from the relations, applies `delta -> 1 - x`, and
  reconstructs rational coefficients from `GF(p)`.

Usage:

```bash
wolframscript -file de_to_wl.wl <dmaster_relations> <output_wl>
```

Output:
- `<output_wl>` containing only:
  - `masters = {...};`
  - `AX = {...};`
