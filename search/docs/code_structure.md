# `search` 模块 `code_structure`

## 1. Overview

`search` 模块实现 `problem_solution.md` 中的两阶段关系搜索。

第一阶段对应 `problem_solution.md` 第 3 节：给定截断级数数据，构造

$$
A^{(\delta)} x = 0,
$$

其中未知量是

$$
x_{\vec{\nu},k} = c_k^{\vec{\nu}}.
$$

第二阶段对应 `problem_solution.md` 第 4 节：从第一阶段零空间生成具体系数赋值，构造 FI 关系

$$
\sum_{\vec{\nu} \in G} r_{\vec{\nu}} FI_{\vec{\nu}} = 0,
$$

再将这些关系组织成矩阵

$$
A^{(FI)} y = 0,
$$

并做第二次消元，输出主积分候选集合与 reduction 关系。

顶层入口函数是

```cpp
void runRelationSearchPipeline(const std::vector<std::string>& seriesPaths,
                               const std::string& configPath,
                               const std::string& targetPath,
                               const std::string& maxSearchDegPath,
                               const std::string& outputPath);
```

它按固定顺序执行：输入读取 -> 第一阶段矩阵构造与消元 -> 具体赋值展开 -> 第二阶段矩阵构造与消元 -> 结果格式化输出。

## 2. Code conventions

### 2.1 File layout

```text
search/
├── docs/
│   ├── problem_solution.md
│   ├── code_structure.md
│   └── code_structure.backup.md
├── include/
│   ├── linear.hpp
│   ├── io.hpp
│   ├── relation_types.hpp
│   ├── relation_matrix_builder.hpp
│   ├── relation_searcher.hpp
│   ├── coefficient_relation_expander.hpp
│   ├── fi_reduction_builder.hpp
│   ├── fi_reduction_searcher.hpp
│   ├── relation_formatter.hpp
│   └── search_pipeline.hpp
├── src/
│   ├── linear.tpp
│   ├── io.tpp
│   ├── relation_types.tpp
│   ├── relation_matrix_builder.tpp
│   ├── relation_searcher.tpp
│   ├── coefficient_relation_expander.tpp
│   ├── fi_reduction_builder.tpp
│   ├── fi_reduction_searcher.tpp
│   ├── relation_formatter.tpp
│   └── search_pipeline.tpp
└── test/
    ├── search_relations/
    └── FI_solve/
```

### 2.2 Symbol table

| Math symbol | Code name | Type / carrier |
|:---|:---|:---|
| $\vec{\nu}$ | `IntegralLabel::nu` | `std::vector<int>` |
| $(\vec{\nu}, b)$ | `SeriesLabel` | `IntegralLabel` + `int` |
| $FI_{\vec{\nu}}^{(b)}(\delta)$ truncated series | `SeriesSample<T>` | `coeffs[r] = a_r^{\vec{\nu},(b)}` |
| $d$ | `degreeD` | `int` |
| $m$ | `maxDeltaDegreeM` | `int` |
| $N_{\mathrm{bc}}$ | `numFBIMasters` | `int` |
| $c_k^{\vec{\nu}}$ | `RelationVariable{integral, k}` | first-stage variable |
| $A^{(\delta)}$ | `std::vector<std::vector<T>>` from `RelationMatrixBuilder` | first-stage matrix |
| $\ker A^{(\delta)}$ | `RelationSearchResult<T>` | first-stage elimination result |
| one concrete coefficient solution | `CoefficientAssignment<T>` | free-column assignment |
| $r_{\vec{\nu}} = \sum_{k=0}^{m} c_k^{\vec{\nu}}$ | `FIRelation<T>::coeffs` | second-stage relation coefficient |
| $A^{(FI)}$ | `std::vector<std::vector<T>>` from `FIReductionBuilder` | second-stage matrix |
| free FI variables | `FIReductionResult<T>::freeColumns` | master candidates |

### 2.3 Ordering convention

变量按 `problem_solution.md` 第 2.5 节的复杂度规则排序。

- 第一阶段排序器：`RelationVariableMoreComplexFirst`
- 第二阶段排序器：`IntegralLabelMoreComplexFirst`

排序规则必须保持：

1. `props` 大者在前；
2. 若 `props` 相同，`dots` 大者在前；
3. 第一阶段中若积分相同，较大的 $k$ 在前；
4. 最后按固定字典序比较 `nu`。

### 2.4 Template design

除顶层 pipeline 入口外，搜索模块的核心组件都以 `template<typename T>` 编写，当前实际标量类型为 `FlintMod`。因此所有算法层都只依赖如下抽象能力：

- `T(0)` 与 `T(1)`
- 加减乘
- 相等比较
- 流输出

## 3. Architecture

### 3.1 Class hierarchy

```text
SearchInput<T>
    ├── RelationMatrixBuilder<T>
    │       └── RelationSearcher<T>
    │               └── RelationSearchResult<T>
    │                       └── CoefficientRelationExpander<T>
    │                               ├── CoefficientAssignment<T>
    │                               └── FIRelation<T>
    └── RelationFormatter<T>

std::vector<FIRelation<T>>
    ├── FIReductionBuilder<T>
    │       └── FIReductionSearcher<T>
    │               └── FIReductionResult<T>
    └── RelationFormatter<T>
```

`LinearSystem<T>` 是两阶段共同依赖的底层消元组件，不承载业务语义，只负责对给定矩阵做消元并暴露 RREF、pivot 列和 free 列。

### 3.2 Data flow

```text
config / target / max degree / series files
    ↓
parseSearchConfigFile / parseSearchTargetFile / parseSeriesFile
    ↓
SearchInput<FlintMod>
    ↓
RelationMatrixBuilder<FlintMod>::buildVariables()
    ↓
RelationMatrixBuilder<FlintMod>::buildMatrix()
    ↓
A^(delta)
    ↓
LinearSystem<FlintMod>::eliminate()
    ↓
RelationSearchResult<FlintMod>
    ↓
CoefficientRelationExpander<FlintMod>::expandAssignments()
    ↓
CoefficientRelationExpander<FlintMod>::buildFIRelations()
    ↓
std::vector<FIRelation<FlintMod>>
    ↓
FIReductionBuilder<FlintMod>::buildIntegralVariables()
    ↓
FIReductionBuilder<FlintMod>::buildMatrix()
    ↓
A^(FI)
    ↓
LinearSystem<FlintMod>::eliminate()
    ↓
FIReductionResult<FlintMod>
    ↓
RelationFormatter<FlintMod>
```

### 3.3 Call graph

```text
runRelationSearchPipeline(...)
    -> parseSearchConfigFile(...)
    -> parseMaxSearchDegreeFile(...)
    -> parseSearchTargetFile(...)
    -> parseSeriesFile<FlintMod>(...) [for each boundary condition]
    -> RelationSearcher<FlintMod>::search()
         -> RelationMatrixBuilder<FlintMod>::buildVariables()
         -> RelationMatrixBuilder<FlintMod>::buildMatrix(...)
         -> LinearSystem<FlintMod>::eliminate()
    -> CoefficientRelationExpander<FlintMod>::expandAssignments(...)
    -> CoefficientRelationExpander<FlintMod>::buildFIRelations(...)
    -> FIReductionSearcher<FlintMod>::search()
         -> FIReductionBuilder<FlintMod>::buildIntegralVariables()
         -> FIReductionBuilder<FlintMod>::buildMatrix(...)
         -> LinearSystem<FlintMod>::eliminate()
    -> RelationFormatter<FlintMod>::write*(...)
```

## 4. Class-by-class specification

### 4.1 `IntegralLabel`, `SeriesLabel`, `SeriesSample<T>`, `SearchInput<T>`

**Files**: `search/include/relation_types.hpp`, `search/src/relation_types.tpp`

**Purpose**: 表示第一阶段输入数据与标签系统，对应 `problem_solution.md` 第 2.1–2.4 节。

```cpp
struct IntegralLabel {
    std::vector<int> nu;
};

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

- **Pre**:
  - 所有 `targets` 的 `nu` 维度一致。
  - `samples` 覆盖全部 `(target, bcIndex)` 组合。
  - 每条 `coeffs` 长度为 `degreeD + 1`。
- **Post**:
  - `SearchInput<T>` 足以唯一确定第一阶段矩阵构造。

### 4.2 `RelationVariable`, `RelationSearchResult<T>`, `CoefficientAssignment<T>`, `FIRelation<T>`, `FIReductionResult<T>`

**Files**: `search/include/relation_types.hpp`, `search/src/relation_types.tpp`

**Purpose**: 承载两阶段中间结果和最终结果。

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

template<typename T>
struct FIRelation {
    std::vector<IntegralLabel> integrals;
    std::vector<T> coeffs;
};

template<typename T>
struct FIReductionResult {
    std::vector<IntegralLabel> integrals;
    std::vector<std::vector<T>> rrefMatrix;
    std::vector<int> pivotColumns;
    std::vector<int> freeColumns;
};
```

- **Implements**:
  - `RelationVariable`: `problem_solution.md` 第 3.1 节
  - `RelationSearchResult<T>`: `problem_solution.md` 第 3.4 节
  - `CoefficientAssignment<T>` 与 `FIRelation<T>`: `problem_solution.md` 第 4.1–4.2 节
  - `FIReductionResult<T>`: `problem_solution.md` 第 4.4 节

### 4.3 `RelationVariableMoreComplexFirst` 与 `IntegralLabelMoreComplexFirst`

**Files**: `search/include/relation_types.hpp`, `search/src/relation_types.tpp`

**Purpose**: 实现 `problem_solution.md` 第 2.5 节的复杂度排序规则。

```cpp
struct RelationVariableMoreComplexFirst {
    bool operator()(const RelationVariable& lhs,
                    const RelationVariable& rhs) const;
};

struct IntegralLabelMoreComplexFirst {
    bool operator()(const IntegralLabel& lhs,
                    const IntegralLabel& rhs) const;
};
```

- **Pre**: `lhs`、`rhs` 的 `nu` 维度可比较。
- **Post**:
  - 排序结果稳定。
  - 第一阶段比较器额外使用 `k` 打破同一积分上的顺序。

### 4.4 `RelationMatrixBuilder<T>`

**File**: `search/include/relation_matrix_builder.hpp`, `search/src/relation_matrix_builder.tpp`

**Purpose**: 实现 `problem_solution.md` 第 3.2–3.3 节，把截断级数数据变成第一阶段线性系统。

```cpp
template<typename T>
class RelationMatrixBuilder {
public:
    explicit RelationMatrixBuilder(const SearchInput<T>& input);

    std::vector<RelationVariable> buildVariables() const;
    std::vector<std::vector<T>> buildMatrix(
        const std::vector<RelationVariable>& variables) const;

private:
    const SearchInput<T>& input_;

    const SeriesSample<T>& findSample(const std::vector<int>& nu, int bcIndex) const;
};
```

- **Pre**:
  - `input_.targets` 非空。
  - `input_.samples` 完整覆盖所有 `(target, bcIndex)`。
  - 每条 `coeffs` 长度为 `input_.degreeD + 1`。
- **Post**:
  - `buildVariables()` 返回全部 $(\vec{\nu},k)$，数量为 $|G|(m+1)$。
  - `buildMatrix()` 返回矩阵大小为 $(d+1)N_{\mathrm{bc}} \times |G|(m+1)$。
  - 矩阵元满足

    $$
    A^{(\delta)}_{(b,n),(\vec{\nu},k)} =
    \begin{cases}
    a_{n-k}^{\vec{\nu},(b)}, & n \ge k, \\
    0, & n < k.
    \end{cases}
    $$

**Critical method pseudocode**:

```text
buildVariables():
    vars = []
    for target in input_.targets:
        for k in [0, maxDeltaDegreeM]:
            vars.push_back((target, k))
    sort vars by RelationVariableMoreComplexFirst
    return vars

buildMatrix(variables):
    matrix = []
    for bcIndex in [0, numFBIMasters):
        for n in [0, degreeD]:
            row = zero vector of length variables.size()
            for col, var in enumerate(variables):
                if n < var.k:
                    continue
                sample = findSample(var.integral.nu, bcIndex)
                row[col] = sample.coeffs[n - var.k]
            matrix.push_back(row)
    return matrix
```

### 4.5 `RelationSearcher<T>`

**File**: `search/include/relation_searcher.hpp`, `search/src/relation_searcher.tpp`

**Purpose**: 封装第一阶段完整搜索，对应 `problem_solution.md` 第 3.4 节。

```cpp
template<typename T>
class RelationSearcher {
public:
    explicit RelationSearcher(const SearchInput<T>& input);

    RelationSearchResult<T> search() const;

private:
    const SearchInput<T>& input_;
};
```

- **Pre**: `RelationMatrixBuilder<T>` 的前置条件全部满足。
- **Post**:
  - `result.variables` 与第一阶段矩阵列严格对应。
  - `result.rrefMatrix`、`pivotColumns`、`freeColumns` 对应 $A^{(\delta)}$ 的 RREF 结构。

### 4.6 `CoefficientRelationExpander<T>`

**File**: `search/include/coefficient_relation_expander.hpp`, `search/src/coefficient_relation_expander.tpp`

**Purpose**: 实现 `problem_solution.md` 第 4.1–4.2 节，从零空间参数化结果构造具体关系，并在 $\delta = 1$ 处得到 FI 关系。

```cpp
template<typename T>
class CoefficientRelationExpander {
public:
    std::vector<CoefficientAssignment<T>> expandAssignments(
        const RelationSearchResult<T>& result) const;

    std::vector<FIRelation<T>> buildFIRelations(
        const std::vector<CoefficientAssignment<T>>& assignments) const;
};
```

- **Pre**:
  - `result.rrefMatrix` 与 `result.variables` 列数一致。
  - `pivotColumns` 与 `freeColumns` 来自合法 RREF。
- **Post**:
  - 每个自由列生成一组具体 `CoefficientAssignment<T>`。
  - 每组赋值满足“所选自由列取 1，其余自由列取 0”。
  - 每条 `FIRelation<T>` 的系数等于同一积分上所有 $k$ 系数之和。

**Critical method pseudocode**:

```text
expandAssignments(result):
    assignments = []
    for freeCol in result.freeColumns:
        assignment.values = zero vector
        assignment.values[freeCol] = 1
        for row from last pivot row downto first:
            pivotCol = result.pivotColumns[row]
            sum = 0
            for col in (pivotCol + 1) .. end:
                sum += result.rrefMatrix[row][col] * assignment.values[col]
            assignment.values[pivotCol] = -sum
        assignments.push_back(assignment)
    return assignments

buildFIRelations(assignments):
    relations = []
    for assignment in assignments:
        relation = empty coefficient map indexed by IntegralLabel
        for i in [0, assignment.variables.size()):
            relation[assignment.variables[i].integral] += assignment.values[i]
        relations.push_back(relation)
    return relations
```

### 4.7 `FIReductionBuilder<T>`

**File**: `search/include/fi_reduction_builder.hpp`, `search/src/fi_reduction_builder.tpp`

**Purpose**: 实现 `problem_solution.md` 第 4.3 节，构造第二阶段矩阵 $A^{(FI)}$。

```cpp
template<typename T>
class FIReductionBuilder {
public:
    explicit FIReductionBuilder(const std::vector<FIRelation<T>>& relations);

    std::vector<IntegralLabel> buildIntegralVariables() const;
    std::vector<std::vector<T>> buildMatrix(
        const std::vector<IntegralLabel>& integrals) const;

private:
    const std::vector<FIRelation<T>>& relations_;
};
```

- **Pre**: 每条 `FIRelation<T>` 中 `integrals.size() == coeffs.size()`。
- **Post**:
  - `buildIntegralVariables()` 收集第二阶段出现的全部积分并去重排序。
  - `buildMatrix()` 逐行写入每条 FI 关系的系数。

**Critical method pseudocode**:

```text
buildIntegralVariables():
    integrals = []
    for relation in relations_:
        for integral in relation.integrals:
            if integral not yet present:
                integrals.push_back(integral)
    sort integrals by IntegralLabelMoreComplexFirst
    return integrals

buildMatrix(integrals):
    matrix = []
    for relation in relations_:
        row = zero vector of length integrals.size()
        for col, integral in enumerate(integrals):
            if integral occurs in relation:
                row[col] = matching coefficient
        matrix.push_back(row)
    return matrix
```

### 4.8 `FIReductionSearcher<T>`

**File**: `search/include/fi_reduction_searcher.hpp`, `search/src/fi_reduction_searcher.tpp`

**Purpose**: 封装第二阶段完整消元，对应 `problem_solution.md` 第 4.4 节。

```cpp
template<typename T>
class FIReductionSearcher {
public:
    explicit FIReductionSearcher(const std::vector<FIRelation<T>>& relations);

    FIReductionResult<T> search() const;

private:
    const std::vector<FIRelation<T>>& relations_;
};
```

- **Pre**: `FIReductionBuilder<T>` 的前置条件满足。
- **Post**:
  - `result.integrals` 与第二阶段矩阵列严格对应。
  - `result.freeColumns` 给出当前关系集下的自由积分，即 master candidates。
  - `result.pivotColumns` 给出可由更简单 FI 表示的积分。

### 4.9 `RelationFormatter<T>`

**File**: `search/include/relation_formatter.hpp`, `search/src/relation_formatter.tpp`

**Purpose**: 将两阶段结果转成可读文本输出。

```cpp
template<typename T>
class RelationFormatter {
public:
    static void writeSummary(
        std::ostream& out,
        const RelationSearchResult<T>& result);

    static void writeRelations(
        std::ostream& out,
        const RelationSearchResult<T>& result);

    static void writeRREF(
        std::ostream& out,
        const RelationSearchResult<T>& result);

    static void writeAssignments(
        std::ostream& out,
        const std::vector<CoefficientAssignment<T>>& assignments);

    static void writeFIRelations(
        std::ostream& out,
        const std::vector<FIRelation<T>>& relations);

    static void writeFIReductionSummary(
        std::ostream& out,
        const FIReductionResult<T>& result);

    static void writeFIMasterBasis(
        std::ostream& out,
        const FIReductionResult<T>& result);

    static void writeFIReductions(
        std::ostream& out,
        const FIReductionResult<T>& result);

    static void writeFIRREF(
        std::ostream& out,
        const FIReductionResult<T>& result);
};
```

- **Pre**: 所有结果对象内部列顺序一致。
- **Post**:
  - 第一阶段输出包含摘要、relation、RREF、具体赋值。
  - 第二阶段输出包含 FI relations、FI reduction 摘要、reduction、FI RREF。

**Note**:

- `writeFIMasterBasis()` 已存在，但当前 `runRelationSearchPipeline(...)` 未调用它；保留该接口是为了显式输出主积分候选集合。

### 4.10 `io.hpp`

**File**: `search/include/io.hpp`, `search/src/io.tpp`

**Purpose**: 读取配置、目标积分、最大搜索次数和 series 文件，构造第一阶段输入。

```cpp
struct SearchConfig {
    int nuSize = 0;
    int degreeD = 0;
    int numFBIMasters = 0;
    mp_limb_t p = 0;
};

SearchConfig parseSearchConfigFile(const std::string& path);
int parseMaxSearchDegreeFile(const std::string& path);
std::vector<IntegralLabel> parseSearchTargetFile(const std::string& path, int expectedNuSize);

template<typename T>
std::vector<SeriesSample<T>> parseSeriesFile(
    const std::string& path,
    const std::vector<IntegralLabel>& targets,
    int degreeD,
    int bcIndex);
```

- **Pre**:
  - config 文件提供 `N`、`deg`、`bc`、`p`。
  - target 文件中的每行都能解析为一个 `nu`。
  - series 文件行数与 target 数目一致。
- **Post**:
  - 解析结果足以无损构造 `SearchInput<T>`。

### 4.11 `runRelationSearchPipeline(...)`

**File**: `search/include/search_pipeline.hpp`, `search/src/search_pipeline.tpp`

**Purpose**: 提供顶层固定工作流。

```cpp
void runRelationSearchPipeline(const std::vector<std::string>& seriesPaths,
                               const std::string& configPath,
                               const std::string& targetPath,
                               const std::string& maxSearchDegPath,
                               const std::string& outputPath);
```

- **Pre**:
  - `seriesPaths.size() == cfg.numFBIMasters`
  - `outputPath` 可写
- **Post**:
  - 输出文件按固定顺序包含两阶段结果
  - 顶层调用始终使用 `FlintMod` 作为有限域标量类型

**Execution order pseudocode**:

```text
read config / degree / targets
set modulus
parse all series files
assemble SearchInput<FlintMod>
run first-stage search
expand assignments
build FI relations
run second-stage reduction search
write summaries, relations, assignments, reductions, and RREF blocks
```

## 5. Verification strategy

### 5.1 Unit checks

- `relation_types.*`
  - 验证 `countProps()`、`countDots()` 与复杂度排序
  - 验证 `equalNu()` 与字符串化函数
- `io.*`
  - 验证合法输入可完整解析
  - 验证维度不匹配、行数不匹配、文件缺失时抛异常
- `relation_matrix_builder.*`
  - 对手工小例子检查矩阵尺寸和卷积条目
- `coefficient_relation_expander.*`
  - 对手工 RREF 检查回代
  - 检查同一积分不同 $k$ 的系数确实求和为 $r_{\vec{\nu}}$
- `fi_reduction_builder.*`
  - 检查积分去重与排序
  - 检查每条 FI relation 都被正确写入第二阶段矩阵

### 5.2 Stage-I checks

`test/search_relations/` 必须至少验证：

1. 变量总数等于 $|G|(m+1)$；
2. 矩阵行数等于 $(d+1)N_{\mathrm{bc}}$；
3. `freeColumns.size()` 与人工预期一致；
4. 至少一条显式 `c_k^{\vec{\nu}}` relation 与人工分析一致。

### 5.3 Stage-II checks

`test/FI_solve/` 必须至少验证：

1. 每个自由列都生成一条具体赋值；
2. `buildFIRelations()` 后的系数与 $\sum_k c_k^{\vec{\nu}}$ 一致；
3. 第二阶段排序与复杂度规则一致；
4. `writeFIReductions()` 输出的 pivot FI 只依赖排在其后的更简单 FI。

### 5.4 Regression outputs

顶层 pipeline 输出至少检查以下区段：

- `[relations]`
- `[assignments]`
- `[fi_relations]`
- `[fi_reductions]`
- `[fi_rref]`
- `[rref]`

### 5.5 Incremental verification order

推荐自底向上的验证顺序：

1. `relation_types.*`
2. `io.*`
3. `relation_matrix_builder.*`
4. `relation_searcher.*`
5. `coefficient_relation_expander.*`
6. `fi_reduction_builder.*`
7. `fi_reduction_searcher.*`
8. `relation_formatter.*`
9. `search_pipeline.*`

## 6. References

- `search/docs/problem_solution.md`
- `search/include/relation_types.hpp`
- `search/include/relation_matrix_builder.hpp`
- `search/include/relation_searcher.hpp`
- `search/include/coefficient_relation_expander.hpp`
- `search/include/fi_reduction_builder.hpp`
- `search/include/fi_reduction_searcher.hpp`
- `search/include/relation_formatter.hpp`
- `search/include/io.hpp`
- `search/include/search_pipeline.hpp`
