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
(Used by `poly_relation_searcher`)
- Contains finite-field and search settings.
- Typical keys include:
  - `p`
  - `m`

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
./poly_relation_searcher <G_path> <series_path> <config_path> <output_path>
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
./integral_solver <G_path> <poly_relation_path> <delta_value> [output_path]
```

- Default `output_path`: `integral_solution`

Output order in `integral_solution`:
1. `# p`, `# m`, `# delta`, `# |G|`
2. `# integral variables`, `# integral pivot columns`, `# integral free columns`
3. `#MIs`
4. `[reductions]`
5. `[relations]`
6. `[integral_rref]`
