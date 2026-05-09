# expand/tools

This directory contains standalone command-line tools built on top of `expand` core code.

## Build

Each tool has its own `Makefile`.

```bash
cd expand/tools/<tool_name>
make
```

## Common Input Files

### `S`
- A Mathematica-style matrix expression: `{{...},{...},...}`.
- Variables are `X` and `Y`.
- Matrix must be square with size `(B + N) x (B + N)`.

### `config`
Required keys:
- `B = <int>`
- `N = <int>`
- `deg = <int>`
- `p = <uint64>`
- `a = <uint64>`
- `b = <uint64>`
- `d = <uint64>`
- `reduceMode = normal | maximal_cut`
- `bc = {u0,u1,...}`  (size must match number of masters)

Optional keys:
- `print2DMode = target | cache` (used by `dump_2dseries`, default `target`)
- `sector = {n1,n2,...,nN}` (used by `dump_cz`; if omitted, topsector `{1,...,1}`)

### `target`
- One integral index vector per line, e.g. `{1,1,1,1}`.
- Each vector length must be `N`.

## Tools

### 1) `fi_pipeline_runner`

Path: `expand/tools/fi_pipeline_runner`

Purpose:
- Compute final FI 1D series for all entries in `target`.

Usage:

```bash
./fi_pipeline_runner <S_path> <config_path> <target_path> <output_path>
```

Output:
- `<output_path>`: one FI 1D series per target line.

### 2) `dump_2dseries`

Path: `expand/tools/dump_2dseries`

Purpose:
- Dump FBI 2D series.
- Behavior controlled by `print2DMode` in `config`.

Usage:

```bash
./dump_2dseries <S_path> <config_path> <target_path> <output_path>
```

Output:
- If `print2DMode = target`:
  - `<output_path>` contains 2D series for each target integral.
- If `print2DMode = cache`:
  - `<output_path>` contains full cache dump.
  - `<output_path>.keys` contains extracted cache keys (`nu;delta`).

### 3) `dump_cz`

Path: `expand/tools/dump_cz`

Purpose:
- Dump `C/z` data before and after GiNaC -> Flint conversion.
- Sector is selected by `sector` in `config` (or topsector if omitted).

Usage:

```bash
./dump_cz <S_path> <config_path> <output_prefix>
```

Output:
- `<output_prefix>.txt`: human-readable report.
- `<output_prefix>_data.wl`: Mathematica data for comparison.

---


### 4) `dump_masters`

Path: `expand/tools/dump_masters`

Purpose:
- Detect and dump master sectors (`nu`) for a family.

Usage:

```bash
./dump_masters <S_path> <config_path> <output_path>
```

Output:
- `<output_path>`: one master `nu` per line, e.g. `{1,0,1,...}`.


### 4) `dump_masters`

Path: `expand/tools/dump_masters`

Purpose:
- Detect and dump master sectors (`nu`) for a family.

Usage:

```bash
./dump_masters <S_path> <config_path> <output_path>
```

Output:
- `<output_path>`: one master `nu` per line, e.g. `{1,0,1,...}`.
