# Reduction From Cut Expansions

## Part 1. Problem Definition

### 1. Original Reduction Problem

For integrals $I_i(\delta)$, the target is

$$
I_i(\delta)=\sum_j R_{ij}(\delta) M_j(\delta),
\qquad
R_{ij}(\delta)\in \mathbb{Z}_p(\delta),
$$

where $M_j$ are master integrals. The problem is to determine both the masters
and the rational coefficients.

### 2. Cut Expansion Form

Each integral decomposes over cuts/boundary conditions:

$$
I_i(\delta)=\sum_{c\in \mathrm{cuts}} ir^{(c)} I_i^{(c)}(\delta).
$$

Here $ir^{(c)}$ is the irrational leading coefficient of cut $c$, and
$I_i^{(c)}$ is the normalized contribution after that leading coefficient is
fixed.

Substituting the cut decomposition into the original reduction relation gives

$$
I_i^{(c)}(\delta)
=
\sum_j R_{ij}(\delta) M_j^{(c)}(\delta),
\qquad \forall c.
$$

Thus the unknown irrational coefficients are irrelevant. The needed data are
finite Taylor expansions

$$
\{I_i^{(c),N}(\delta)\}_{i,c}.
$$

For sufficiently large $N$, these expansions should contain enough information
to reconstruct the finite-degree rational functions $R_{ij}(\delta)$.

### 3. Primitive Reconstruction

For a large enough finite integral set $G=\{I_i\}$, make the ansatz

$$
\sum_{i\in G} P_i(\delta) I_i^{(c),N}(\delta)=0,
\qquad \forall c,
$$

with sufficiently high but bounded-degree unknown polynomials $P_i(\delta)$.
Expanding in $\delta$ for all cuts gives a linear system over $\mathbb{Z}_p$.

Solving these systems gives polynomial relations among integrals in $G$.
Choosing an independent subset as masters then gives rational reductions of the
remaining integrals. This primitive algorithm is inefficient, but it shows why
finite cut expansions are sufficient in principle.

## Part 2. Current Algorithm

### 1. Expansion

The expansion step generates $I_i^{(c),N}(\delta)$. Details are in
`expand/docs/`.

The current workflow also uses symmetry. If two sectors are mapped to each
other without changing the branches, their FBI are the same, and the remaining
part of the integrand is also unchanged. Therefore the full objects are
identical. Such sectors form one orbit, and only one representative sector is
kept.

For each representative cut sector, the expansion target should include:

- integrals in the representative sector;
- integrals equivalent to the representative sector by `sectormap`;
- dot <= 2 integrals whose sector contains the current sector or an equivalent
  sector as a subsector.

These targets should be generated in cut mode. Target generation therefore
uses both the `sectormap` and the case2/3/4 routing table.

### 2. Master Search

The search is done cut by cut. For each cut $c$, we choose a local integral set
$G_c$ and search for relations using the series on that cut. The relation
search itself is documented in `search/docs/`.

The key point is the distinction between sector and cut:

- Sector: determined by the FBI/FI index vector $\nu$.
- Cut: an independent MFBI leading coefficient.

They are not one-to-one. Non-case1 sectors have no MFBI cut. Case1 sectors in
the same symmetry orbit correspond to the same object, so only the
representative sector gives a cut to be considered.

For a cut associated with sector $s$, integrals below $s$ vanish on that cut,
so the basic local choice is the integrals in the representative sector $s$.
This is the analogue of ordinary sector-cut searches.

Two additions are needed for completeness:

- Non-case1 sectors can still contribute masters. Their integrals are appended
  to the first-lower subsector with
  valid cuts.
- Some shifted-$\nu$ integrals vanish on the current cut but can still be
  masters globally. These are also passed to the first-lower subsectors with
  valid cuts.

Schematically,

$$
G_c
=
G_s
+G_{\mathrm{noncase1\ append}}
+G_{\mathrm{zero\ append}}.
$$

The purpose is to keep each local search large enough to find the full master
space relevant to later reductions.

### 3. BL Sector Reduction

After the masters and relations are found, `bl_sector_reduce` performs the
sector-by-sector reconstruction. Its details are in `bl_sector_reduce/docs/`.

If the master search is complete, this step should only apply the already
detected reduction structure.

## Part 3. Operation Flow

This layer should mostly organize existing `expand`, `search`, and
`bl_sector_reduce` code. It should not require rewriting their core algorithms,
though small tools may be added.

### 1. Generate Symmetry Data

Start from the Mathematica symmetry finder used before, with two filters:

- discard maps only between 2-branch sectors;
- discard maps that convert propagators across different branches.

The remaining maps define sector orbits and representative sectors. Record the
discarded maps as well, since they may be useful later.

### 2. Route Non-Case1 Sectors

For every case2/3/4 sector, find the nearest lower layer containing case1
subsectors. Record all such case1 subsectors.

This routing table is later used in two places:

- generating `expand` targets;
- enlarging the search set $G_c$ for master search.

### 3. Detect Zero-on-Cut Integrals

For each representative sector, generate all dot <= 2 integrals in that sector
and compute short series, e.g. degree 10.

If all FBI components for a given $\nu$ are zero on the current cut, record
that $\nu$. These zero-on-cut integrals should be passed to the first-layer
subsectors with valid cuts, because they may still be global masters.

### 4. Run Expansion

For each representative cut sector, build the `expand` target in cut mode.
The target contents are the ones specified in Part 2.

### 5. Search Masters

For each cut, build the local search set $G_c$ using the rules in Part 2:

$$
G_c
=
G_s
+G_{\mathrm{noncase1\ append}}
+G_{\mathrm{zero\ append}}.
$$

Then run the existing `search` workflow on the corresponding cut series.

### 6. Run BL Sector Reduction

Use the masters and relations found by `search` to run `bl_sector_reduce`.

For testing, reduce all dot <= 1 integrals. This should be computationally
acceptable and is broad enough to expose missing masters or missing routing
data.
