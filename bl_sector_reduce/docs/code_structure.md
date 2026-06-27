# BL Sector Reduction：代码结构与实现文档

## 1. 概述

本文档把 [`problem_solution.md`](./problem_solution.md) 中的 sector-wise BL 符号约化问题映射到代码结构。`bl_sector_reduce` 只实现该文档定义的路线：

1. 读取对象在各 FBI sector 边界条件下的一维 $\delta$ 级数。
2. 读取每个 sector 的 master 列表。
3. 根据传播子支持集构造 sector tree，并从根到叶处理。
4. 对每个对象构造当前 sector contribution $C_s(A)$。
5. 用 BL / approximant basis 搜索

   $$
   P_0 C_s(A)=\sum_j P_j M_j.
   $$

6. 输出符号多项式/有理函数约化。

`expand` 和 `search` 只作为文件风格、有限域类型和标签文本格式参考；本模块不复用 `search` 的全局关系搜索、integral solver 或 dense reduction 逻辑。

## 2. 文件组织

```text
bl_sector_reduce/
├── docs/
│   ├── problem_solution.md
│   └── code_structure.md
├── include/
│   ├── bl_config.hpp
│   ├── label.hpp
│   ├── polynomial_1d.hpp
│   ├── sector_tree.hpp
│   ├── series_store.hpp
│   ├── master_data.hpp
│   ├── contribution.hpp
│   ├── approximant_basis.hpp
│   ├── reducer.hpp
│   ├── io.hpp
│   └── formatter.hpp
├── src/
│   ├── bl_config.tpp
│   ├── label.tpp
│   ├── polynomial_1d.tpp
│   ├── sector_tree.tpp
│   ├── series_store.tpp
│   ├── master_data.tpp
│   ├── contribution.tpp
│   ├── approximant_basis.tpp
│   ├── reducer.tpp
│   ├── io.tpp
│   └── formatter.tpp
└── tools/
    ├── bl_sector_reducer/
    │   ├── Makefile
    │   └── bl_sector_reducer.cpp
    └── build_all.sh
```

所有核心实现采用头文件包含 `.tpp` 的模板风格，和 `expand`、`search` 保持一致。

## 3. 核心数据模型

### 3.1 `ObjectLabel`

**文件**：`label.hpp`

`ObjectLabel` 是对象文本标签。它保留已有积分标签的解析能力，但 `bl_sector_reduce` 只要求 label 能稳定比较、打印并提供传播子指数向量用于 sector 支持集。

```cpp
enum class ObjectHead { FI, BFI, BBFI };
enum class BoundaryAxis { X, Y };
enum class BoundarySide { U, D };

struct BoundaryTag {
    BoundaryAxis axis;
    BoundarySide side;
};

struct ObjectLabel {
    ObjectHead head = ObjectHead::FI;
    std::vector<BoundaryTag> boundaries;
    std::vector<int> nu;
};
```

公共函数：

```cpp
ObjectLabel parseObjectLabel(const std::string& text, int expectedNuSize);
std::string objectLabelToString(const ObjectLabel& label);
bool equalObjectLabel(const ObjectLabel& lhs, const ObjectLabel& rhs);
std::vector<int> sectorSupport(const ObjectLabel& label);
```

### 3.2 `SectorId` 与 `SectorTree`

**文件**：`sector_tree.hpp`

`SectorId` 是传播子支持集的二进制向量：

```cpp
struct SectorId {
    std::vector<int> bits;
};
```

`SectorTree` 根据所有出现的 sector 构造父子关系。父节点定义为当前 sector 的直接上游 sector：在已知 sector 集合中，支持集严格包含当前 sector 且支持集大小最小者。若存在多个候选，按字典序选择一个确定父节点。Phase A 假设处理结构为树；若输入 sector 集不满足唯一父节点，代码用该确定性规则构造处理树并记录 warning。

接口：

```cpp
class SectorTree {
public:
    explicit SectorTree(const std::vector<SectorId>& sectors);
    int rootIndex() const;
    const std::vector<int>& processingOrder() const; // root to leaves
    int parentOf(int sectorIndex) const;
    std::vector<int> ancestorsOf(int sectorIndex) const;
    int indexOf(const SectorId& sector) const;
    const SectorId& sectorAt(int idx) const;
};
```

### 3.3 `Polynomial1D<T>`

**文件**：`polynomial_1d.hpp`

一维多项式

$$
P(\delta)=\sum_i p_i\delta^i
$$

用连续向量存储。

```cpp
template<typename T>
class Polynomial1D {
public:
    Polynomial1D();
    explicit Polynomial1D(int degree);
    explicit Polynomial1D(std::vector<T> coeffs);

    int degree() const;
    bool isZero() const;
    const std::vector<T>& coeffs() const;
    T coeff(int i) const;
    void setCoeff(int i, const T& value);
    void trim();

    T eval(const T& x) const;
};
```

### 3.4 `SeriesStore<T>`

**文件**：`series_store.hpp`

保存对象在每个 sector 边界条件下的级数。

```cpp
template<typename T>
class SeriesStore {
public:
    void addSeries(const SectorId& sector,
                   const ObjectLabel& object,
                   std::vector<T> coeffs);

    const std::vector<T>& getSeries(const SectorId& sector,
                                    const ObjectLabel& object) const;

    bool hasSeries(const SectorId& sector,
                   const ObjectLabel& object) const;

    int degree() const;
    std::vector<ObjectLabel> objects() const;
    std::vector<SectorId> sectors() const;
};
```

### 3.5 `MasterData`

**文件**：`master_data.hpp`

保存每个 sector 的 master 列表。

```cpp
class MasterData {
public:
    void setMasters(const SectorId& sector, std::vector<ObjectLabel> masters);
    const std::vector<ObjectLabel>& mastersFor(const SectorId& sector) const;
    bool isMaster(const SectorId& sector, const ObjectLabel& object) const;
};
```

## 4. Sector Contribution

### 4.1 `ReductionTerm<T>` 与 `SectorReduction<T>`

**文件**：`contribution.hpp`

一个 sector contribution 的符号约化为

$$
C_s(A)=\sum_\mu \frac{N_\mu}{D}M_\mu.
$$

代码表示：

```cpp
template<typename T>
struct ReductionTerm {
    ObjectLabel master;
    Polynomial1D<T> numerator;
};

template<typename T>
struct SectorReduction {
    SectorId sector;
    ObjectLabel object;
    Polynomial1D<T> denominator;
    std::vector<ReductionTerm<T>> terms;
    bool isFreeMaster = false;
};
```

### 4.2 `ContributionBuilder<T>`

**文件**：`contribution.hpp`

构造第 $s$ 个 sector 上的

$$
C_s(A)^{(s)}
=A^{(s)}-\sum_{a\in\operatorname{Ancestors}(s)}
\operatorname{EvalAtSector}(\operatorname{Red}_a(A),s).
$$

接口：

```cpp
template<typename T>
class ContributionBuilder {
public:
    ContributionBuilder(const SeriesStore<T>& series,
                        const SectorTree& tree);

    std::vector<T> buildContribution(
        const ObjectLabel& object,
        int sectorIndex,
        const std::vector<SectorReduction<T>>& knownReductions) const;
};
```

关键实现：

1. 取 `series.getSeries(s, object)` 作为起点。
2. 对每个祖先 reduction，把每个 master 的 sector-`s` 级数乘 numerator，多项式卷积后求和。
3. 若有非平凡 denominator，则用一维形式级数除法把 rational contribution 展开到目标阶数。
4. 从起点级数中减去该祖先贡献。

## 5. BL / Approximant Basis

### 5.1 `ApproximantRequest<T>`

**文件**：`approximant_basis.hpp`

```cpp
template<typename T>
struct ApproximantRequest {
    std::vector<T> target;                 // C_s(A)^s
    std::vector<std::vector<T>> masters;   // M_j^s
    int maxDegree = 0;                     // m
    int workOrder = 0;                     // K_work
};
```

### 5.2 `ApproximantResult<T>`

```cpp
template<typename T>
struct ApproximantResult {
    bool success = false;
    std::vector<Polynomial1D<T>> polynomials; // P0, P1, ..., Pr
};
```

### 5.3 `ApproximantBasisSolver<T>`

`ApproximantBasisSolver` 实现第 4-5 节的 BL / approximant basis 搜索。

接口：

```cpp
template<typename T>
class ApproximantBasisSolver {
public:
    ApproximantResult<T> solve(const ApproximantRequest<T>& request) const;
};
```

要求：

1. 构造行向量级数

   $$
   F=(C_s(A)^s,-M_1^s,\ldots,-M_r^s).
   $$

2. 在 order `workOrder` 下计算 approximant basis。
3. 从 basis 中选取满足 `deg(P_i)<=maxDegree` 且 `P0!=0` 的向量。
4. 返回对应 `P0,...,Pr`。

### 5.4 关键算法

实现采用逐阶 discrepancy 更新的 order-basis 算法：

1. 初始化 basis 为单位矩阵，每一行是一个多项式向量。
2. 对 $k=0,\ldots,K-1$：
   1. 计算每个 basis 行在 order $k$ 的 discrepancy。
   2. 若所有 discrepancy 为零，进入下一阶。
   3. 选择 discrepancy 非零且当前 degree 最小的 pivot 行。
   4. 对其他 discrepancy 非零的行执行

      $$
      b_i \leftarrow b_i - \frac{d_i}{d_p} b_p.
      $$

   5. pivot 行乘以 $\delta$。
3. 处理完 $K$ 阶后，basis 行生成 order-$K$ approximant module。
4. 在 basis 行中寻找满足 degree bound 和 $P_0\ne0$ 的候选，并直接验证

   $$
   F P =0 \bmod \delta^K.
   $$

该实现只解决本文档的一行向量 approximant 问题，不实现通用矩阵 approximant API。

## 6. Reducer Pipeline

### 6.1 配置

**文件**：`bl_config.hpp`

```cpp
struct BLSectorConfig {
    int nuSize = 0;
    int degreeD = 0;
    int maxDegree = 0;      // m_user
    int safetyOrder = 10;   // K_safety
    int certOrder = 10;     // K_cert
    mp_limb_t prime = 0;
};
```

### 6.2 `BLSectorReducer<T>`

**文件**：`reducer.hpp`

```cpp
template<typename T>
class BLSectorReducer {
public:
    BLSectorReducer(const BLSectorConfig& config,
                    const SectorTree& tree,
                    const SeriesStore<T>& series,
                    const MasterData& masters);

    std::vector<SectorReduction<T>> reduceAll(
        const std::vector<ObjectLabel>& objects) const;
};
```

流程对应 `problem_solution.md` 第 8 节：

1. 对每个对象按 root-to-leaf 遍历 sector。
2. 构造 `C_s(A)^s`。
3. 若 contribution 为零，跳过。
4. 若对象是 sector master，记录 `isFreeMaster=true`。
5. 计算 `m_supported`。
6. 按指数 degree schedule 尝试 `m`。
7. 成功则记录 `SectorReduction`。
8. 失败则抛出包含 object、sector、`r`、`D`、`m_max`、`K_safety`、`K_cert` 的错误。

## 7. IO

**文件**：`io.hpp`

### 7.1 Config

配置文件为 key-value：

```text
N = 3
deg = 300
m = 10
p = 2305843009213693951
K_safety = 10
K_cert = 10
```

### 7.2 Sector series list

工具读取一个 sector series list 文件：

```text
sector={1,1,1} series_path target_path master_path
sector={1,1,0} series_path target_path master_path
```

每行定义一个 sector 边界条件的数据来源。`series_path` 与 `target_path` 采用 `expand` 当前输出格式：target 每行一个对象标签，series 每行一个 `{c0,...,cD}`，二者行对齐。`master_path` 每行一个 master 标签，允许 `#` 注释。

### 7.3 Object list

object list 每行一个对象标签。若未显式传入 object list，工具默认使用所有 series target 中出现的对象。

## 8. Formatting

**文件**：`formatter.hpp`

输出结构：

```text
# p = ...
# D = ...
# m_user = ...
# K_safety = ...
# K_cert = ...

[masters]
sector={...}
FI{...}

[sector_reductions]
sector={...}
object = FI{...}
den = {d0,d1,...}
term FI{...} = {n0,n1,...}

[global_reductions]
FI{...} = ...
```

`[sector_reductions]` 是机器可读主输出；`[global_reductions]` 可先输出同样信息的文本摘要。Phase C 必须保证多项式系数完整输出。

## 9. Tools

### 9.1 `bl_sector_reducer`

```bash
bl_sector_reducer <config_path> <sector_series_list_path> <output_path> [object_list_path]
```

职责：

1. 解析 config。
2. 读取所有 sector 的 series、target、master。
3. 构造 `SectorTree`。
4. 运行 `BLSectorReducer<FlintMod>`。
5. 写出符号结果。

## 10. 验证策略

### 10.1 单元检查

1. `Polynomial1D`：乘法、除法展开、求值。
2. `SectorTree`：根、父节点、祖先链顺序。
3. `ContributionBuilder`：两层 sector 的祖先贡献剥离。
4. `ApproximantBasisSolver`：已知小例子，例如 $A=(1+\delta)M$ 应返回 $P_0=1$、$P_1=1+\delta$ 或等价关系。

### 10.2 vac 回归

vac 的数据可以作为退化 sector tree 测试。用同一个 sector 的 series 和 master list 运行通用 pipeline，选取一个非 master 对象，例如 `FI{1,1,2}`：

1. 得到符号多项式约化。
2. 对 $\delta=571$ 求值。
3. 与 `search/runs/vac/bc1/integral_solution_dot2deg300m10_delta571` 中同一对象的结果比较。

该测试不改变通用算法；vac 只是输入 sector 集较简单。

## 11. 实现顺序

1. `label`、`polynomial_1d`、`sector_tree`。
2. `series_store`、`master_data`、`io`。
3. `approximant_basis`。
4. `contribution`。
5. `reducer` 和 `formatter`。
6. `bl_sector_reducer` tool。
7. 编译和 vac 对照测试。
