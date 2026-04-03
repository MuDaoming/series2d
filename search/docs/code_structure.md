# 代码结构与实现文档

## 1. 概述

本文档描述 `search` 模块的代码结构与实现分工。数学问题本身见
[`problem_and_workflow.md`](./problem_and_workflow.md)。

整个设计遵循两条基本原则：

1. `KISS`
   只保留当前问题必须的对象与接口，不为了“可能的远期情况”提前引入复杂层次。

2. `SOLID`
   每个类只负责一类明确工作；对象之间通过清晰的数据结构解耦，便于后续扩展到其它“搜索线性关系”问题。

`search` 模块的任务不是计算 $\mathrm{FI}_{\vec{\nu}}$，而是：

1. 读取已经得到的一维级数；
2. 按给定积分集合 $G$ 与次数上界 $m$ 构造线性系统；
3. 求解齐次方程；
4. 先得到 $c_k^{\vec{\nu}}$ 的关系；
5. 再将其整理成最终的 FI 约化关系。

## 2. 目录结构

第一版建议的目录结构如下：

```text
search/
├── docs/
│   ├── problem_and_workflow.md
│   └── code_structure.md
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
    │   ├── S
    │   ├── config
    │   ├── target
    │   ├── maxsearchdeg
    │   ├── test_search.cpp
    │   └── Makefile
    └── FI_solve/
```

其中：

- `linear.hpp/.tpp` 已从 `reconstruct` 复制过来，作为底层消元工具；
- 其余文件由 `search` 自己负责，专门处理关系搜索问题。
- `test/search_relations/` 用于第一阶段测试，即从 FI 级数搜索 $c_k^{\vec{\nu}}$ 的关系；
- `test/FI_solve/` 用于第二阶段测试，即从第一阶段结果出发生成最终 FI 约化关系。

## 3. 总体架构

### 3.1 模块层次

```text
第一阶段：
LinearSystem<T>
    ↓
RelationMatrixBuilder<T>
    ↓
RelationSearcher<T>

第二阶段：
第一阶段零空间
    ↓
CoefficientRelationExpander
    ↓
FIReductionBuilder
    ↓
FIReductionSearcher
    ↓
最终 FI 约化关系
```

### 3.2 数据流

```text
输入：
  - 多组边界条件下的 FI 一维级数
  - 待搜索的积分集合 G
  - 关系次数上界 m
        ↓
  [IO]
        ↓
  规范化为统一的级数数据结构
        ↓
  [RelationMatrixBuilder]
        ↓
  构造第一阶段齐次线性系统 A x = 0
        ↓
  [LinearSystem]
        ↓
  高斯消元，得到第一阶段零空间结构
        ↓
  对自由变量逐个赋值
        ↓
  令 δ = 1
        ↓
  得到一批 FI 之间的线性关系
        ↓
  再次构造线性系统
        ↓
  [LinearSystem]
        ↓
  输出最终 FI 约化关系
```

## 4. 核心数据对象

### 4.1 `SeriesSample`

表示一条输入级数，即某个边界条件下某个积分的一维展开。

```cpp
struct SearchLabel {
    std::vector<int> nu;
    int bcIndex;
};

template<typename T>
struct SeriesSample {
    SearchLabel label;
    std::vector<T> coeffs;
};
```

语义：

- `SearchLabel`：统一表示“哪一个积分、哪一个边界条件”；
- `label.nu`：积分指标 $\vec{\nu}$；
- `label.bcIndex`：边界条件编号；
- `coeffs[r]`：$\delta^r$ 的系数。

引入 `SearchLabel` 的原因是：

- `nu` 和 `bcIndex` 会在多个对象中反复出现；
- 今后如果加入其它类型的搜索任务，可以继续在这个标识对象中增加字段；
- 这样由“样本身份”决定的信息都集中在一个小对象中，避免重复定义与模块耦合。

### 4.2 `SearchInput`

表示一次完整搜索所需的全部输入。

```cpp
template<typename T>
struct SearchInput {
    int degreeD;
    int maxDeltaDegreeM;
    int numFBIMasters;
    std::vector<std::vector<int>> nus;
    std::vector<SeriesSample<T>> samples;
};
```

语义：

- `degreeD`：已知级数的最高阶 $d$；
- `maxDeltaDegreeM`：关系次数上界 $m$；
- `numFBIMasters`：FBI 主积分个数，也就是独立边界条件数；
- `nus`：待搜索的积分集合 $G$；
- `samples`：所有边界条件下的所有输入级数。

### 4.3 `RelationVariable`

表示一个未知量 $c_k^{\vec{\nu}}$。

```cpp
struct RelationVariable {
    SearchLabel label;
    int k;
};
```

语义：

- `label`：变量所属的样本标识；
- `k`：对应 $\delta^k$；

这个对象只表达“这是哪一个未知量”，不直接缓存 `props`、`dots` 之类可由 `label.nu` 计算出来的信息。

这样做的原因是：

- `RelationVariable` 只保留身份，不混入排序逻辑；
- 变量复杂度是比较规则，不是变量本体数据；
- 后续如果比较规则改变，不需要修改变量对象本身。

因此，与变量排序有关的信息应由单独的比较器负责，例如：

```cpp
struct RelationVariableLess {
    bool operator()(const RelationVariable& lhs,
                    const RelationVariable& rhs) const;
};
```

### 4.4 `RelationSearchResult`

表示一次消元后的结果。

```cpp
template<typename T>
struct RelationSearchResult {
    std::vector<RelationVariable> variables;
    std::vector<std::vector<T>> rrefMatrix;
    std::vector<int> pivotColumns;
    std::vector<int> freeColumns;
};
```

语义：

- `variables`：列与变量的对应关系；
- `rrefMatrix`：消元后的矩阵；
- `pivotColumns`：主元列；
- `freeColumns`：自由变量列。

第一版只需要这些信息，就足够生成显式关系。

### 4.5 `CoefficientAssignment`

表示对一组自由变量赋值之后得到的一组具体
$c_k^{\vec{\nu}}$ 解。

```cpp
template<typename T>
struct CoefficientAssignment {
    std::vector<RelationVariable> variables;
    std::vector<T> values;
    int chosenFreeColumn;
};
```

语义：

- `variables[i]` 与 `values[i]` 一一对应；
- `chosenFreeColumn` 记录这组解由哪一个自由变量取 1 得到。

这个对象是第一阶段与第二阶段之间的桥梁。

### 4.6 `FIRelation`

表示一条已经令 $\delta=1$ 之后得到的 FI 关系：

$$
\sum_{\vec{\nu}\in G} r_{\vec{\nu}} \mathrm{FI}_{\vec{\nu}} = 0
$$

```cpp
template<typename T>
struct FIRelation {
    std::vector<IntegralLabel> integrals;
    std::vector<T> coeffs;
};
```

语义：

- `integrals[i]` 与 `coeffs[i]` 一一对应；
- `coeffs[i]` 是对应 $\mathrm{FI}_{\vec{\nu}}$ 的系数。

### 4.7 `FIReductionResult`

表示第二阶段消元后的最终结果。

```cpp
template<typename T>
struct FIReductionResult {
    std::vector<IntegralLabel> integrals;
    std::vector<std::vector<T>> rrefMatrix;
    std::vector<int> pivotColumns;
    std::vector<int> freeColumns;
};
```

语义：

- 变量已经不再是 $c_k^{\vec{\nu}}$；
- 列对象直接对应原始的 $\mathrm{FI}_{\vec{\nu}}$；
- 目标是把复杂积分写成简单积分。

## 5. 核心类设计

### 5.1 `LinearSystem<T>`

文件：

- `search/include/linear.hpp`
- `search/src/linear.tpp`

来源：

- 直接复制自 `reconstruct/include/linear.hpp`

当前能力：

1. 保存矩阵；
2. 做高斯消元；
3. 记录主元列与自由变量列；
4. 输出变量表达式。

在 `search` 中的定位：

- 它只负责“给定矩阵后如何消元”；
- 不负责构造矩阵；
- 不负责决定变量顺序；
- 不负责把结果格式化为关系搜索的输出。

因此，`search` 不直接围绕 `LinearSystem` 写业务逻辑，而是把它当作底层工具。

### 5.2 `RelationMatrixBuilder<T>`

文件：

- `search/include/relation_matrix_builder.hpp`
- `search/src/relation_matrix_builder.tpp`

作用：

- 从输入级数构造齐次线性系统 $A x = 0$；
- 固定变量顺序；
- 建立矩阵列与变量 $c_k^{\vec{\nu}}$ 的对应关系。

建议接口：

```cpp
template<typename T>
class RelationMatrixBuilder {
public:
    RelationMatrixBuilder(const SearchInput<T>& input);

    std::vector<RelationVariable> buildVariables() const;
    std::vector<std::vector<T>> buildMatrix(
        const std::vector<RelationVariable>& variables) const;
};
```

实现要点：

1. 先生成所有变量 $c_k^{\vec{\nu}}$；
2. 用独立比较器按“复杂到简单”排序：
   - `props` 大的在前；
   - 若相同，`dots` 大的在前；
   - 若相同，`k` 大的在前；
   - 若仍相同，再用固定规则比较 `nu`；
3. 对每个边界条件 `b` 和每个阶数 `n=0,\dots,d` 生成一行；
4. 该行第 `(\vec{\nu},k)` 列的元素为
   $a_{n-k}^{\vec{\nu},(b)}$ 或 0。

### 5.3 `RelationSearcher<T>`

文件：

- `search/include/relation_searcher.hpp`
- `search/src/relation_searcher.tpp`

作用：

- 组织一次完整的搜索；
- 调用 `RelationMatrixBuilder` 构造矩阵；
- 调用 `LinearSystem` 消元；
- 整理成统一结果对象。

建议接口：

```cpp
template<typename T>
class RelationSearcher {
public:
    RelationSearcher(const SearchInput<T>& input);

    RelationSearchResult<T> search() const;
};
```

固定流程：

1. 生成变量列表；
2. 构造矩阵；
3. 建立 `LinearSystem<T>`；
4. 调用 `eliminate()`；
5. 读取 pivot / free variable 信息；
6. 返回 `RelationSearchResult<T>`。

### 5.4 `RelationFormatter<T>`

文件：

- `search/include/relation_formatter.hpp`
- `search/src/relation_formatter.tpp`

作用：

- 把消元结果写成适合阅读和后续检查的形式。

建议接口：

```cpp
template<typename T>
class RelationFormatter {
public:
    static void writeRelations(
        std::ostream& out,
        const RelationSearchResult<T>& result);

    static void writeRREF(
        std::ostream& out,
        const RelationSearchResult<T>& result);
};
```

输出形式分两种：

1. 显式关系  
   例如
   $$
   c_9^{\vec{\nu}_1} = \alpha_1 c_1^{\vec{\nu}_2} + \alpha_2 c_0^{\vec{\nu}_3}
   $$

2. 行最简形矩阵  
   用于进一步人工分析或调试。

第一版中，`writeRelations` 应优先实现。

但在两阶段结构中，`RelationFormatter<T>` 不再只负责第一阶段输出，还应支持：

1. 第一阶段输出  
   $c_k^{\vec{\nu}}$ 的 relation 和 RREF；

2. 第二阶段输出  
   最终的 FI reduction。

也就是说，这个类保留原名字，但职责会扩展为“统一格式化第一阶段与第二阶段结果”。

### 5.5 `CoefficientRelationExpander<T>`

文件：

- `search/include/coefficient_relation_expander.hpp`
- `search/src/coefficient_relation_expander.tpp`

作用：

- 从第一阶段零空间结果生成一组具体的 $c_k^{\vec{\nu}}$ 解；
- 再令 $\delta=1$；
- 输出 FI 关系。

建议接口：

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

实现要点：

1. 对每个自由列单独生成一组解；
2. 该自由列取 1，其余自由列取 0；
3. 用第一阶段 RREF 回代得到所有 pivot 变量；
4. 对每个固定的 $\vec{\nu}$，把
   $$
   \sum_{k=0}^{m} c_k^{\vec{\nu}}
   $$
   作为 $\mathrm{FI}_{\vec{\nu}}$ 的系数；
5. 这样就得到一批 `FIRelation<T>`。

### 5.6 `FIReductionBuilder<T>`

文件：

- `search/include/fi_reduction_builder.hpp`
- `search/src/fi_reduction_builder.tpp`

作用：

- 把一批 `FIRelation<T>` 重新组织成第二阶段的线性系统；
- 固定 FI 的变量顺序。

建议接口：

```cpp
template<typename T>
class FIReductionBuilder {
public:
    explicit FIReductionBuilder(const std::vector<FIRelation<T>>& relations);

    std::vector<IntegralLabel> buildIntegralVariables() const;
    std::vector<std::vector<T>> buildMatrix(
        const std::vector<IntegralLabel>& integrals) const;
};
```

实现要点：

1. 变量直接是 $\mathrm{FI}_{\vec{\nu}}$；
2. 变量顺序仍按“复杂到简单”；
3. 每条 `FIRelation<T>` 对应矩阵的一行。

### 5.7 `FIReductionSearcher<T>`

文件：

- `search/include/fi_reduction_searcher.hpp`
- `search/src/fi_reduction_searcher.tpp`

作用：

- 调用 `FIReductionBuilder<T>` 构造第二阶段矩阵；
- 再调用 `LinearSystem<T>` 做消元；
- 输出最终的 FI reduction。

建议接口：

```cpp
template<typename T>
class FIReductionSearcher {
public:
    explicit FIReductionSearcher(const std::vector<FIRelation<T>>& relations);

    FIReductionResult<T> search() const;
};
```

固定流程：

1. 生成 FI 变量列表；
2. 构造第二阶段矩阵；
3. 调用 `LinearSystem<T>` 消元；
4. 返回 `FIReductionResult<T>`。

### 5.8 `IO`

文件：

- `search/include/io.hpp`
- `search/src/io.tpp`

作用：

- 读取搜索任务输入；
- 将外部文件整理成 `SearchInput<T>`。

第一版不需要设计复杂格式，接口只要能完成：

1. 读取待搜索的 $\vec{\nu}$ 集合；
2. 读取边界条件个数；
3. 读取每组边界条件下每个积分的一维级数；
4. 组装成 `SearchInput<T>`。

### 5.9 `runRelationSearchPipeline(...)`

文件：

- `search/include/search_pipeline.hpp`
- `search/src/search_pipeline.tpp`

作用：

- 提供顶层入口；
- 从输入路径读数据；
- 顺序执行两阶段搜索；
- 输出结果。

建议接口：

```cpp
template<typename T>
void runRelationSearchPipeline(
    const std::string& inputPath,
    const std::string& outputPath);
```

两阶段固定流程应为：

1. 读取 series 与第一阶段搜索参数；
2. 调用 `RelationSearcher<T>`；
3. 调用 `CoefficientRelationExpander<T>`；
4. 调用 `FIReductionSearcher<T>`；
5. 用 `RelationFormatter<T>` 输出第一阶段与第二阶段结果。

## 6. 测试与 Debug 循环

`search` 的开发与排错必须分成两个阶段进行，不能一开始就把所有问题混在一起。

### 6.1 第一阶段：`test/search_relations/`

这一组测试只负责验证：

1. 输入的一维级数是否正确读入；
2. 第一阶段矩阵是否正确构造；
3. 变量顺序是否符合“复杂到简单”；
4. 第一阶段零空间是否正确；
5. 输出的 $c_k^{\vec{\nu}}$ 关系是否合理。

这一阶段的调试输出应优先检查：

- 读入的 `series`
- 变量列表
- 第一阶段矩阵的行列数
- pivot 列与 free 列
- 第一阶段 relation 输出

如果第一阶段结果不对，不应继续做第二阶段。

### 6.2 第二阶段：`test/FI_solve/`

这一组测试只负责验证：

1. 是否能从第一阶段零空间正确生成具体的 $c_k^{\vec{\nu}}$ 解；
2. 在每组具体解中令 $\delta=1$ 后，是否得到正确的 FI 关系；
3. 第二阶段矩阵是否正确构造；
4. 最终是否能把复杂 FI 写成简单 FI 的线性组合。

这一阶段的调试输出应优先检查：

- 每个自由变量对应生成的具体 $c_k^{\vec{\nu}}$
- 每条由 $\delta=1$ 得到的 FI 关系
- 第二阶段矩阵
- 最终 reduction 输出

### 6.3 推荐开发顺序

推荐的实现与调试顺序是：

1. 先让 `search_relations` 稳定通过；
2. 再实现 `FI_solve`；
3. 最后再把两阶段合并进统一的 pipeline。

这样做的原因是：

- 第一阶段错了，第二阶段一定错；
- 第二阶段的输入完全依赖第一阶段的零空间结果；
- 分阶段测试能显著降低 debug 难度。

第一版允许把所有输入都组织在一个目录下，只要路径语义清楚即可。

## 7. 实现顺序

建议按下面顺序实现：

1. 保留并复用 `linear.hpp/.tpp`；
2. 定义 `SeriesSample`、`SearchInput`、`RelationVariable`、`RelationSearchResult`；
3. 实现 `RelationMatrixBuilder`；
4. 实现 `RelationSearcher`；
5. 实现 `RelationFormatter`；
6. 最后补 `IO` 与 `search_pipeline`。

原因：

- 真正的核心是“矩阵构造是否正确”；
- 只要矩阵构造正确，底层消元已经现成；
- 输出与输入都放在后面做，不会阻塞核心验证。

## 8. 与现有模块的关系

### 8.1 与 `expand` 的关系

`expand` 负责生成 $\mathrm{FI}_{\vec{\nu}}$ 的一维级数。

`search` 只消费这些结果，不参与 FI 级数本身的计算。

### 8.2 与 `reconstruct` 的关系

`search` 当前只复用 `reconstruct` 中的 `linear`：

- 复制 `linear.hpp`
- 复制 `linear.tpp`

其余 `reconstruct` 模块与当前关系搜索问题无直接耦合。

## 9. 第一版输出要求

第一版只要求支持下面两类输出：

1. 显式关系输出  
   目标是让结果尽量呈现为“复杂变量由简单变量表示”。

2. 行最简形矩阵输出  
   便于人工检查秩、自由变量与零空间结构。

不在第一版中处理：

- 输入格式自动推断；
- 多种输出风格切换；
- 复杂的关系筛选策略。
