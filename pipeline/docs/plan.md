# Pipeline Plan for Symmetry-Aware Reduction

This plan describes the orchestration layer needed around the existing
`expand`, `search`, and `bl_sector_reduce` modules. The goal is to organize
their inputs and outputs, not to rewrite their core algorithms.

## 1. New Pipeline Layer

Add a top-level directory:

```text
pipeline/
├── docs/
├── tools/
└── runs/
```

`pipeline/tools/` should contain small glue tools. They may be written in
Wolfram Language, Python, shell, or C++ depending on the input format, but they
should only generate metadata, targets, configs, run lists, and reports.

The existing module tools remain responsible for computation:

- `expand/tools/fi_pipeline_runner`
- `search/tools/poly_relation_searcher`
- `search/tools/integral_solver`
- `bl_sector_reduce/tools/bl_sector_reducer`

## 2. Data Products

The pipeline should produce the following intermediate files.

### 2.1 `sectormap`

Format: existing `expand` sector-map format.

```text
{{source sector} -> {target sector}, {source propagator -> target propagator, ...}}
```

Used by:

- `expand` config through `sectorMap = <path>`;
- target generation;
- representative-sector selection.

### 2.2 `symmetry_rejected`

Records maps found by the symmetry search but not used.

Each entry should include:

- source sector;
- target sector;
- propagator map;
- rejection reason.

Initial rejection reasons:

- `two_branch_only`;
- `cross_branch_propagator_map`;
- `not_representative_needed`;
- `verification_failed`.

### 2.3 `sector_table`

One row per valid 3-branch sector.

Fields:

- `sector`;
- `representative`;
- `is_representative`;
- `case`: case0/case1/case2/case3/case4, using the project convention;
- `has_cut`: true only for representative case1 cut sectors;
- `prop_count`;
- `branch_support`.

Used by all later target and run generation.

### 2.4 `noncase_routes`

Routes case2/3/4 sectors to lower case1 cut sectors.

Each row:

```text
source_sector -> {cut_sector_1, cut_sector_2, ...}
```

The target cut sectors are the nearest lower layer containing case1
representatives. If that layer has no valid cut sector, continue downward.

Used by:

- `expand` target generation;
- `search` local `G` generation.

### 2.5 `zero_routes`

Routes shifted-$\nu$ integrals that vanish on their own cut to lower cut
sectors.

Each row:

```text
cut_sector, nu -> {lower_cut_sector_1, lower_cut_sector_2, ...}
```

The probe should use all dot <= 2 integrals in the current representative
sector and degree 10 series. A `nu` is recorded only if all FBI components for
that `nu` are zero on the current cut.

Used by:

- `search` local `G` generation;
- optional `expand` target generation for lower cuts, if the required series
  are not already present.

### 2.6 Per-Cut `expand` Inputs

For each representative cut sector `s`:

```text
pipeline/runs/<family>/expand/<s>/
├── config
├── target
└── series
```

`config` should be compatible with `fi_pipeline_runner`.

`target` should include:

- integrals in representative sector `s`;
- integrals equivalent to `s` by `sectormap`;
- dot <= 2 integrals whose sector contains `s` or an equivalent sector as a
  subsector.

The `sectormap` and `noncase_routes` are needed to generate this target
correctly.

### 2.7 Per-Cut `search` Inputs and Outputs

For each representative cut sector `s`:

```text
pipeline/runs/<family>/search/<s>/
├── config
├── G
├── series_list
├── poly_relation
├── integral_solution
└── masters
```

`G` is the local search set:

$$
G_s
+G_{\mathrm{noncase1\ append}}
+G_{\mathrm{zero\ append}}.
$$

`series_list` points to the available per-cut series files. Matching must be by
full integral tag, not by `nu` only.

`masters` is extracted from `integral_solution`.

### 2.8 `bl_sector_reduce` Inputs and Outputs

For a family run:

```text
pipeline/runs/<family>/bl_sector_reduce/
├── config
├── sector_series_list
├── object_list_dot1
└── reduction_dot1
```

`sector_series_list` should point to each representative sector's:

- `series`;
- `target`;
- `masters`.

`object_list_dot1` contains all dot <= 1 test integrals.

## 3. Proposed Tools

### 3.1 `pipeline/tools/generate_symmetry.wl`

Purpose:

Generate usable sector symmetry data from the existing Mathematica symmetry
finder.

Inputs:

- family definition or `{U,F}` data;
- branch assignment of propagators;
- valid sector list;
- options:
  - keep only 3-branch sectors;
  - reject 2-branch-only maps;
  - reject maps converting propagators across branches.

Outputs:

- `sectormap`;
- `symmetry_rejected`;
- `symmetry_orbits`.

Notes:

- This can reuse and slightly modify `symmetry/wolfram/SectorSymmetry.wl`.
- Every accepted map should be verified before writing `sectormap`.

### 3.2 `pipeline/tools/build_sector_table`

Purpose:

Classify sectors and choose representative cut sectors.

Inputs:

- `S` matrix;
- base `expand` config;
- `sectormap`;
- branch assignment.

Outputs:

- `sector_table`;
- `cut_sectors`.

Required information:

- sector case;
- representative sector;
- whether the sector has an independent cut.

Implementation options:

- C++ tool using `expand` `Family`/`Sector` code;
- or WL tool if the family case data are easier to obtain there.

### 3.3 `pipeline/tools/build_noncase_routes`

Purpose:

Route case2/3/4 sectors to the nearest lower case1 cut sectors.

Inputs:

- `sector_table`;
- `sectormap`.

Outputs:

- `noncase_routes`.

Rules:

- Start from a non-case1 sector.
- Inspect lower subsectors by one fewer active propagator.
- Keep representative case1 cut sectors.
- If none exist, continue downward.
- Record all nearest valid cut sectors.

### 3.4 `pipeline/tools/build_expand_targets`

Purpose:

Generate per-cut `expand` configs and targets.

Inputs:

- base `expand` config;
- `S` path;
- `sector_table`;
- `sectormap`;
- `noncase_routes`;
- max dot for target generation, initially `2`.

Outputs for each cut sector:

- `expand/<sector>/config`;
- `expand/<sector>/target`;
- a run manifest line for `fi_pipeline_runner`.

Target rules:

- include current representative sector integrals;
- include equivalent-sector integrals by `sectormap`;
- include dot <= 2 integrals whose sector contains the current sector or an
  equivalent sector as a subsector.

### 3.5 `pipeline/tools/probe_zero_integrals`

Purpose:

Detect shifted-$\nu$ integrals whose series vanish on their own cut.

Inputs:

- `S` path;
- base or generated `expand` configs;
- `sector_table`;
- `sectormap`;
- max dot, initially `2`;
- probe degree, initially `10`.

Outputs:

- `zero_probe/<sector>/target`;
- `zero_probe/<sector>/series`;
- `zero_routes`;
- `zero_report`.

Computation:

- For each representative cut sector, generate all dot <= 2 integrals in that
  sector.
- Run `fi_pipeline_runner` with `deg = 10`.
- Group by `nu`; if every FBI component for that `nu` is zero, record it.
- Route the recorded `nu` to first-layer valid cut subsectors.

### 3.6 `pipeline/tools/build_search_inputs`

Purpose:

Generate local `search` inputs for each cut sector.

Inputs:

- `sector_table`;
- `noncase_routes`;
- `zero_routes`;
- per-cut `expand` targets and series;
- search parameter template: `p`, `deg`, `m`, `ncheck`.

Outputs for each cut sector:

- `search/<sector>/config`;
- `search/<sector>/G`;
- `search/<sector>/series_list`;
- run manifest lines for:
  - `poly_relation_searcher`;
  - `integral_solver`.

`G` must include:

- current representative-sector search objects;
- objects routed from non-case1 sectors;
- objects routed from zero-on-cut probes.

### 3.7 `pipeline/tools/extract_search_masters`

Purpose:

Normalize `search` output into per-sector master files.

Inputs:

- `search/<sector>/integral_solution`.

Outputs:

- `search/<sector>/masters`.

Rules:

- Preserve full integral tags, including FI/BFI/BBFI heads and boundary tags.
- Ignore comments and non-master sections.

### 3.8 `pipeline/tools/build_bl_inputs`

Purpose:

Build `bl_sector_reduce` inputs from per-cut expand/search outputs.

Inputs:

- `sector_table`;
- per-cut `expand/<sector>/target`;
- per-cut `expand/<sector>/series`;
- per-cut `search/<sector>/masters`;
- BL config template.

Outputs:

- `bl_sector_reduce/config`;
- `bl_sector_reduce/sector_series_list`;
- `bl_sector_reduce/object_list_dot1`.

Rules:

- `sector_series_list` contains representative cut sectors only.
- `object_list_dot1` contains all dot <= 1 test objects.

### 3.9 `pipeline/tools/run_manifest`

Purpose:

Provide a simple reproducible execution list.

Inputs:

- generated manifests from previous tools.

Outputs:

- `run_expand.sh`;
- `run_search.sh`;
- `run_bl_sector_reduce.sh`;
- optional `status_report`.

These scripts should call existing compiled tools and write outputs under
`pipeline/runs/<family>/`.

## 4. End-to-End Flow

For a family `<family>`:

1. Generate symmetry:
   - input: family polynomials/branches/sectors;
   - output: `sectormap`, `symmetry_rejected`, `symmetry_orbits`.

2. Build sector metadata:
   - input: `S`, base config, `sectormap`;
   - output: `sector_table`, `cut_sectors`.

3. Route non-case1 sectors:
   - input: `sector_table`;
   - output: `noncase_routes`.

4. Build and run zero probes:
   - input: `sector_table`, `S`, base config;
   - output: `zero_routes`, `zero_report`.

5. Build and run expand:
   - input: `sector_table`, `sectormap`, `noncase_routes`;
   - output: per-cut `target` and `series`.

6. Build and run search:
   - input: per-cut series, `noncase_routes`, `zero_routes`;
   - output: per-cut `masters`.

7. Build and run BL sector reduction:
   - input: per-cut series, targets, masters;
   - output: reductions for all dot <= 1 test objects.

8. Summarize:
   - list accepted/rejected symmetry maps;
   - list representative cuts;
   - list non-case1 routes;
   - list zero-on-cut routes;
   - list masters per sector;
   - list successful and failed dot <= 1 reductions.

## 5. Minimal First Implementation

The first useful version should implement only enough to reproduce the current
db/scarecrow workflow without manual file edits:

1. `generate_symmetry.wl`
2. `build_sector_table`
3. `build_noncase_routes`
4. `probe_zero_integrals`
5. `build_expand_targets`
6. `build_search_inputs`
7. `build_bl_inputs`

The actual execution can initially be done by generated shell manifests. A
single monolithic driver is not necessary for the first version.
