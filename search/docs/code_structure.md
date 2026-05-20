# `search` Code Structure

## 1. Overview

`search` implements two-stage linear algebra for integral tags:

1. Find polynomial relations among delta series.
2. Evaluate those relations at a chosen delta value and solve reductions among integral tags.

An integral tag is represented by a head, optional boundary tags, and a `nu` vector:

```text
FI{nu}
BFI[XU]{nu}
BBFI[XU,YD]{nu}
```

The FBI terminology is retained only for the boundary-condition count `numFBIMasters`; it refers to the upstream FBI-series basis and is independent of whether the final tag is FI, BFI, or BBFI.

## 2. Files

```text
search/
├── include/
│   ├── relation_types.hpp
│   ├── io.hpp
│   ├── linear.hpp
│   ├── relation_matrix_builder.hpp
│   ├── relation_searcher.hpp
│   ├── coefficient_relation_expander.hpp
│   ├── integral_reduction_builder.hpp
│   ├── integral_reduction_searcher.hpp
│   ├── relation_formatter.hpp
│   └── search_pipeline.hpp
├── src/
│   ├── relation_types.tpp
│   ├── io.tpp
│   ├── linear.tpp
│   ├── relation_matrix_builder.tpp
│   ├── relation_searcher.tpp
│   ├── coefficient_relation_expander.tpp
│   ├── integral_reduction_builder.tpp
│   ├── integral_reduction_searcher.tpp
│   ├── relation_formatter.tpp
│   └── search_pipeline.tpp
└── tools/
    ├── poly_relation_searcher/
    └── integral_solver/
```

## 3. Core Types

### 3.1 Integral Tags

Defined in `include/relation_types.hpp`:

```cpp
enum class IntegralHead {
    FI,
    BFI,
    BBFI
};

enum class BoundaryAxis {
    X,
    Y
};

enum class BoundarySide {
    U,
    D
};

struct BoundaryTag {
    BoundaryAxis axis;
    BoundarySide side;
};

struct IntegralLabel {
    IntegralHead head = IntegralHead::FI;
    std::vector<BoundaryTag> boundaries;
    std::vector<int> nu;
};
```

`IntegralLabel` is the identity of a search object. Equality compares head, boundary tags, and `nu`. Matching by `nu` alone is not valid.

### 3.2 Series and Search Input

```cpp
struct SeriesLabel {
    IntegralLabel integral;
    int bcIndex = 0;
};

template<typename T>
struct SeriesSample {
    SeriesLabel label;
    std::vector<T> coeffs;
};

template<typename T>
struct SearchInput {
    int degreeD = 0;
    int maxDeltaDegreeM = 0;
    int numFBIMasters = 0;
    std::vector<IntegralLabel> targets;
    std::vector<SeriesSample<T>> samples;
};
```

`numFBIMasters` is the number of FBI boundary-condition series files. It is not renamed when FI is extended to BFI/BBFI.

### 3.3 Relation Variables and Results

```cpp
struct RelationVariable {
    IntegralLabel integral;
    int k = 0;
};

template<typename T>
struct RelationSearchResult {
    std::vector<RelationVariable> variables;
    std::vector<std::vector<T>> rrefMatrix;
    std::vector<int> pivotColumns;
    std::vector<int> freeColumns;
};

template<typename T>
struct CoefficientAssignment {
    std::vector<RelationVariable> variables;
    std::vector<T> values;
    int chosenFreeColumn = -1;
};
```

### 3.4 Integral Relations and Reductions

```cpp
template<typename T>
struct IntegralRelation {
    std::vector<IntegralLabel> integrals;
    std::vector<T> coeffs;
};

template<typename T>
struct IntegralReductionResult {
    std::vector<IntegralLabel> integrals;
    std::vector<std::vector<T>> rrefMatrix;
    std::vector<int> pivotColumns;
    std::vector<int> freeColumns;
};
```

These replace the old FI-only names. The second-stage solver works on general integral tags.

## 4. Ordering

The reduction order is implemented by:

```cpp
RelationVariableMoreComplexFirst
IntegralLabelMoreComplexFirst
IntegralLabelLess
```

Complexity order:

```text
FI > BFI > BBFI
```

Within the same head, compare `nu` by:

1. larger `props` first;
2. larger `dots` first;
3. reverse lexicographic `nu` tie-break.

Boundary tags are only deterministic tie-breaks. `IntegralLabelLess` is a stable map key order and is not the reduction complexity order.

For `RelationVariable`, compare the integral tag first, then larger delta power `k` first.

## 5. IO

Implemented in `include/io.hpp` and `src/io.tpp`.

### 5.1 Config

`parseSearchConfigFile` reads:

```text
N
deg
bc
ncheck
p
```

`bc` is used only to count `numFBIMasters`.

### 5.2 Targets

`parseSearchTargetFile` supports:

```text
FI{1,1,1}
BFI[XU]{1,1,1}
BBFI[XU,YD]{1,1,1}
{1,1,1}
```

Bare `{...}` is parsed as legacy `FI{...}`.

Validation rules:

```text
FI   : no boundary tags
BFI  : exactly one boundary tag
BBFI : exactly two boundary tags, one X* and one Y*
```

### 5.3 Series

`parseSeriesFile` reads line-aligned series for a given target list. The tool-level `poly_relation_searcher` can also read a large series target and filter it to a smaller `G`; this filtering uses full `IntegralLabel` identity.

## 6. Stage I Components

### 6.1 `RelationMatrixBuilder<T>`

Builds variables:

```text
(IntegralLabel, k), 0 <= k <= m
```

and constructs the matrix rows:

```text
row = (bcIndex, delta order n)
column = (integral tag alpha, polynomial degree k)
```

Matrix entry:

```text
sample(alpha, bcIndex).coeffs[n-k], if n >= k
0, otherwise
```

Samples are found by full integral tag.

### 6.2 `RelationSearcher<T>`

Uses `LinearSystem<T>` to compute the RREF of the stage-I matrix. The `poly_relation_searcher` tool adds the train/check split for nullspace non-shrink validation.

## 7. Stage II Components

### 7.1 `CoefficientRelationExpander<T>`

`expandAssignments` converts free columns of the stage-I RREF into concrete coefficient assignments.

`buildIntegralRelations` evaluates each assignment at a finite-field delta value:

```text
r_alpha = sum_k c_{alpha,k} delta^k
```

The result is an `IntegralRelation<T>`.

### 7.2 `IntegralReductionBuilder<T>`

Collects all integral tags from all `IntegralRelation<T>` objects, sorts them by `IntegralLabelMoreComplexFirst`, and builds the second-stage reduction matrix.

### 7.3 `IntegralReductionSearcher<T>`

Runs `LinearSystem<T>` on the second-stage matrix and returns an `IntegralReductionResult<T>`. Free columns are master integrals.

## 8. Formatting

Implemented in `relation_formatter`.

Stage-I output sections:

```text
[relations]
[rref]
```

Stage-II output sections:

```text
# integral variables
# integral pivot columns
# integral free columns
#MIs
[reductions]
[relations]
[integral_rref]
```

Integral variables are printed as full tags, for example:

```text
FI{1,1,1}
BFI[XU]{1,1,1}
BBFI[XU,YD]{1,1,1}
```

## 9. Tools

### 9.1 `poly_relation_searcher`

```bash
poly_relation_searcher <config> <G_file> <series_list> <output>
```

`series_list` contains one line per FBI boundary condition:

```text
series_path target_path
```

The target file for the series may contain more tags than `G`. Filtering uses full integral tags.

### 9.2 `integral_solver`

```bash
integral_solver <G_path> <poly_relation_path> <delta_value> [output_path]
```

It reconstructs the stage-I variable order from `G` and `m`, reads the stage-I RREF, evaluates at `delta_value`, and solves the integral reduction system.

Default output path:

```text
integral_solution
```

## 10. Current Validation Pattern

The current validation checks include:

1. FI-only vac, dp, and dp_planar runs reproduce the previous master sets and reductions after normalizing output names.
2. Mixed FI/BFI/BBFI dp and dp_planar runs produce master files grouped by `nu`.
3. For dp dot3 with `deg=1000`, `m=10`, and `delta=571`, the mixed search gives 43 master integrals.
4. For dp_planar dot2 with `deg=500`, `m=10`, and `delta=571`, the mixed search gives 32 master integrals.
