# Integral-tag polynomial relation search and reduction

## 1. Research Objective

`search` consumes one-dimensional delta expansions produced by `expand` and searches linear relations among a finite set of integral tags.

An integral tag is one of:

```text
FI{nu}
BFI[XU]{nu}
BFI[XD]{nu}
BFI[YU]{nu}
BFI[YD]{nu}
BBFI[XU,YU]{nu}
BBFI[XU,YD]{nu}
BBFI[XD,YU]{nu}
BBFI[XD,YD]{nu}
```

Here `nu` is the propagator-index vector. `FI`, `BFI`, and `BBFI` are part of the integral identity; two tags with the same `nu` but different heads or boundary tags are different search variables.

The input set is

$$
G = \{\alpha_1,\ldots,\alpha_{N_G}\},
$$

where each $\alpha$ is a full integral tag. For each $\alpha \in G$, consider its delta-dependent object

$$
I_\alpha(\delta).
$$

The goal has two stages:

1. Search polynomial-coefficient relations among the $I_\alpha(\delta)$.
2. Evaluate those relations at a chosen finite-field value of $\delta$ and solve reductions among the $I_\alpha(\delta_0)$.

## 2. Stage I: Polynomial Relations in Delta

For a fixed polynomial degree bound $m$, search relations of the form

$$
\sum_{\alpha \in G}
P_\alpha(\delta) I_\alpha(\delta)=0,
$$

where

$$
P_\alpha(\delta)=\sum_{k=0}^m c_{\alpha,k}\delta^k,
\qquad c_{\alpha,k}\in \mathbb{Z}_p.
$$

For each FBI boundary condition $b$, the known truncated series is

$$
I_\alpha^{(b)}(\delta)
=\sum_{r=0}^{d} a_{\alpha,b,r}\delta^r.
$$

Substituting the series into the ansatz gives, for every $b$ and every $0\le n\le d$,

$$
\sum_{\alpha\in G}
\sum_{k=0}^{\min(n,m)}
c_{\alpha,k}a_{\alpha,b,n-k}=0.
$$

This is a homogeneous linear system

$$
A^{(\delta)}x=0.
$$

Rows are indexed by $(b,n)$, columns by $(\alpha,k)$, and

$$
A^{(\delta)}_{(b,n),(\alpha,k)} =
\begin{cases}
a_{\alpha,b,n-k}, & n\ge k,\\
0, & n<k.
\end{cases}
$$

The number of rows is $(d+1)N_{\mathrm{bc}}$, where $N_{\mathrm{bc}}$ is carried in code as `numFBIMasters`, because these boundary conditions come from the FBI layer.

## 3. Nullspace Non-shrink Check

Finite-order data can create accidental relations. To detect this, the stage-I search can split the available orders into a training window and a held-out check window.

Given `ncheck`, define

$$
d_{\mathrm{train}}=d-\mathrm{ncheck}.
$$

The solver computes the nullspace using rows up to $d_{\mathrm{train}}$, then checks whether the same nullspace satisfies the held-out rows $d_{\mathrm{train}}+1,\ldots,d$.

If the nullspace shrinks on the held-out rows, the output is marked as no certified solution. Otherwise, the relation data are written as an RREF and can be consumed by the second-stage solver.

## 4. Stage II: Delta Evaluation and Integral Reductions

Given one stage-I nullspace basis vector, evaluate its polynomial coefficients at a chosen finite-field value $\delta_0$:

$$
r_\alpha =
\sum_{k=0}^{m} c_{\alpha,k}\delta_0^k.
$$

This produces a relation among integral tags:

$$
\sum_{\alpha\in G} r_\alpha I_\alpha(\delta_0)=0.
$$

All such relations are assembled into a second linear system

$$
A^{(I)}y=0,
$$

whose columns are full integral tags $\alpha$. Gaussian elimination then expresses more complex integral tags in terms of simpler tags. Free columns are reported as master integrals.

The output sections are:

```text
#MIs
[reductions]
[relations]
[integral_rref]
```

## 5. Integral Ordering

The reduction ordering is defined on full integral tags.

The head is compared first:

```text
FI > BFI > BBFI
```

Thus BFI is always simpler than FI, and BBFI is always simpler than BFI, regardless of `nu`.

If the head is the same, compare the `nu` vector:

1. Larger `props(nu)` is more complex.
2. If `props` is equal, larger `dots(nu)` is more complex.
3. If still tied, use reverse lexicographic comparison of `nu`.

Boundary tags do not represent physical complexity in this ordering. They are used only as deterministic tie-breaks.

For stage-I variables $(\alpha,k)$, compare $\alpha$ first. If the integral tag is identical, larger $k$ is more complex.

## 6. Target and Series Matching

The target file `G` contains one integral tag per line. Supported forms are:

```text
FI{1,1,1}
BFI[XU]{1,1,1}
BBFI[XU,YD]{1,1,1}
```

A bare vector is accepted as legacy input:

```text
{1,1,1}
```

and is interpreted as:

```text
FI{1,1,1}
```

Series files are line-aligned with their own target files. When a stage-I `G` is a subset of a larger series target, matching is done by full integral tag, not only by `nu`. This is essential because `FI{nu}`, `BFI[XU]{nu}`, and `BBFI[XU,YD]{nu}` are distinct objects.

## 7. Tool Flow

Stage I:

```bash
poly_relation_searcher <config> <G_file> <series_list> <output>
```

Stage II:

```bash
integral_solver <G_path> <poly_relation_path> <delta_value> [output_path]
```

`poly_relation_searcher` writes the stage-I RREF. `integral_solver` rebuilds the same ordered stage-I variables from `G` and `m`, evaluates coefficients at `delta_value`, solves the integral reduction system, and writes the master basis and reductions.
