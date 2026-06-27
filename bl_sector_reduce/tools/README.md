# bl_sector_reduce/tools

This directory contains standalone command-line tools built on top of
`bl_sector_reduce` core code.

## Build

Each tool has its own `Makefile`.

```bash
cd bl_sector_reduce/tools/<tool_name>
make
```

Or build all tools:

```bash
bash bl_sector_reduce/tools/build_all.sh
```

## Common Input Files

### `config`

Used by `bl_sector_reducer`.

Required keys:

- `N = <int>`: propagator-index vector length.
- `deg = <int>`: available delta expansion degree `D`.
- `m = <int>`: user maximum polynomial degree `m_user`.
- `p = <uint64>`: finite-field prime.

Optional keys:

- `K_safety = <int>`: safety order, default `10`.
- `K_cert = <int>`: certification order, default `10`.

For a sector with `r` masters and trial degree `m`, the BL working order is:

```text
(r + 1) * (m + 1) + K_safety + K_cert
```

### `sector_series_list`

Used by `bl_sector_reducer`.

One sector boundary condition per line:

```text
sector={1,1,1} <series_path> <target_path> <master_path>
sector={1,1,0} <series_path> <target_path> <master_path>
```

Paths are resolved relative to the `sector_series_list` file.

- `sector={...}` is the FBI sector support vector.
- `series_path` is a 1D delta series file.
- `target_path` is line-aligned with `series_path`.
- `master_path` lists the masters for this sector.

The `series_path` / `target_path` format follows `expand`:

- `target_path`: one object label per line, such as `FI{1,1,1}`,
  `BFI[XU]{1,1,1}`, or `BBFI[XU,YD]{1,1,1}`.
- `series_path`: one `{c0,c1,...,cD}` coefficient list per target line.

The `master_path` format:

- one master object label per non-comment line;
- lines beginning with `#` are ignored.

### `object_list`

Optional input for `bl_sector_reducer`.

- One object label per line.
- If omitted, `bl_sector_reducer` processes all objects appearing in the loaded
  sector series data.

## Tools

### 1) `bl_sector_reducer`

Path: `bl_sector_reduce/tools/bl_sector_reducer`

Purpose:

- Run the sector-wise BL symbolic reduction pipeline.
- For each object and sector, build the sector contribution and reconstruct:

```text
P0(delta) * C_s(A) = sum_j Pj(delta) * M_j
```

Usage:

```bash
./bl_sector_reducer <config_path> <sector_series_list_path> <output_path> [object_list_path]
```

Inputs:

- `<config_path>`: finite-field, expansion degree, and BL degree/order settings.
- `<sector_series_list_path>`: sector data sources.
- `[object_list_path]`: optional subset of objects to reduce.

Output:

- `<output_path>` with sections:
  - `[sector_reductions]`: machine-readable symbolic reductions.
  - `[global_reductions]`: text summary.
- The header records metadata such as `p`, `D`, `m_user`, `K_safety`,
  `K_cert`, sector count, input master count, and reduction count.

In `[sector_reductions]`, each block contains:

```text
sector={...}
object=FI{...}
free_master=0
den={d0,d1,...}
term FI{...}={n0,n1,...}
term BFI[XU]{...}={m0,m1,...}
```

This means:

```text
object contribution at sector =
  (sum term_polynomial(delta) * master(delta)) / den(delta)
```

If subtracting all more-complex sector contributions leaves an exactly zero
series, the block remains successful but has no `den` or `term`:

```text
sector={...}
object=FI{...}
status=success
zero
```

### 2) `eval_reduction`

Path: `bl_sector_reduce/tools/eval_reduction`

Purpose:

- Evaluate one symbolic reduction from `bl_sector_reducer` at a finite-field
  value of `delta`.
- This is a verification/helper tool; it does not run reduction.

Usage:

```bash
./eval_reduction <reduction_path> <object_label> <delta_value> <nu_size>
```

Inputs:

- `<reduction_path>`: output of `bl_sector_reducer`.
- `<object_label>`: object to evaluate, for example `FI{1,1,2}`.
- `<delta_value>`: finite-field value for `delta`.
- `<nu_size>`: propagator-index vector length.

Output:

- One evaluated relation printed to stdout, for example:

```text
FI{1,1,2} = c1*FI{1,1,1} + c2*BFI[YD]{1,1,1} + ...
```

## Example: vac check

The repository contains a vac test input:

```text
bl_sector_reduce/runs/vac_bc1/sector_series_list
bl_sector_reduce/runs/vac_bc1/object_FI112
```

Run one-object symbolic reduction:

```bash
bl_sector_reduce/tools/bl_sector_reducer/bl_sector_reducer \
  search/runs/vac/bc1/config \
  bl_sector_reduce/runs/vac_bc1/sector_series_list \
  bl_sector_reduce/runs/vac_bc1/reduction_FI112 \
  bl_sector_reduce/runs/vac_bc1/object_FI112
```

Evaluate at `delta = 571`:

```bash
bl_sector_reduce/tools/eval_reduction/eval_reduction \
  bl_sector_reduce/runs/vac_bc1/reduction_FI112 \
  'FI{1,1,2}' \
  571 \
  3
```

This evaluated result matches:

```text
search/runs/vac/bc1/integral_solution_dot2deg300m10_delta571
```
