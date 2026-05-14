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
- One `nu` vector per line, e.g. `{1,1,1,1}`.

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
(Used by `fi_solver`)
- Output of `poly_relation_searcher`.
- Must contain header fields and `[rref]` block.

## Tools

### 1) `poly_relation_searcher`

Path: `search/tools/poly_relation_searcher`

Purpose:
- Build and solve polynomial-relation system from `G` and 2D expansion data.
- Output relation matrix / RREF for next-stage FI solving.

Usage:

```bash
./poly_relation_searcher <G_path> <series_path> <config_path> <output_path>
```

Output:
- `<output_path>` with metadata header and `[rref]` block.

### 2) `fi_solver`

Path: `search/tools/fi_solver`

Purpose:
- Consume `G` and stage-I `poly_relation` result.
- Build FI-level reductions and FI relations.

Usage:

```bash
./fi_solver <G_path> <poly_relation_path> <delta_value> [output_path]
```

- Default `output_path`: `fi_solution`

Output order in `fi_solution`:
1. `# p`, `# m`, `# delta`, `# |G|`
2. `# FI variables`, `# FI pivot columns`, `# FI free columns`
3. `#MIs`
4. `[fi_reductions]`
5. `[fi_relations]`
6. `[fi_rref]`
