# 代码结构与实现文档

## 1. 概述

本文档详细描述项目的代码结构和实现细节。完整的数学理论背景请参考 [`problem_solution.md`](./problem_solution.md)。

### 1.1 项目目标

在质数域 $\mathbb{Z}_p$ 下进行端到端计算：
- 重定义后 FBI（$\widetilde{I} = U^{\text{pow}_U} \cdot I$）的二维幂级数展开（逐阶递推）
- 统一二维被积函数级数构造（$G_{\nu}=P \cdot \widetilde{I}$，无需 $U^\gamma$ 级数）
- FI、BFI、BBFI 的一维 $\delta$ 级数输出（对同一个二维级数采用不同投影规则）

### 1.2 核心方法

将IBP约化和微分方程统一为**线性递推关系**（LRR）形式，直接计算级数系数，而非先计算符号有理函数再做级数运算。约化和微分方程中通过 `Redefinition` 结构体处理 $U$ 幂差因子。

## 2. 代码规范

### 2.1 文件组织

```
expand/
├── docs/                          # 文档
│   ├── problem_solution.md        # 问题定义与计算流程
│   ├── code_structure.md          # 代码架构（本文档）
│   └── todo.md                    # 待办事项
├── include/                       # 头文件（类声明）
│   ├── ff_type.hpp                # 有限域类型 FlintMod
│   ├── rational.hpp               # 多项式 Polynomial 和有理函数 Rational
│   ├── series.hpp                 # 二维幂级数 Series
│   ├── sector.hpp                 # Sector类（含辅助函数）
│   ├── family.hpp                 # Family类
│   ├── series_solver.hpp          # SeriesSolver类 + Redefinition结构体
│   ├── converter.hpp              # GiNaC到FlintMod类型转换
│   ├── integral_tag.hpp           # FI/BFI/BBFI target tag
│   ├── integrand_expander.hpp     # 统一二维被积函数级数构造
│   ├── delta_projector.hpp        # 二维级数到FI/BFI/BBFI一维delta级数
│   ├── io.hpp                     # S/config/target 输入解析
│   └── fi_pipeline.hpp            # 顶层端到端流程入口
├── src/                           # 实现文件（模板实现）
│   └── *.tpp                      # 对应头文件的模板实现
└── test/                          # 测试代码
    └── onedimseries/              # fi_pipeline_runner构建与入口
```

### 2.2 命名规范

- **文件名**：小写字母加下划线（snake_case），头文件用 `.hpp`，模板实现用 `.tpp`
- **类名**：大驼峰（PascalCase），如 `SeriesSolver`、`FlintMod`
- **成员变量**：小驼峰加下划线后缀，如 `numProps_`、`targetDeg_`
- **成员函数**：小驼峰（camelCase），如 `getCase()`、`solveLRRAtDeg()`
- **模板参数**：大写字母缩写，如 `RT`（有理函数类型）、`PT`（多项式类型）、`ST`（标量类型）

### 2.3 模板设计

| 模板参数 | 符号计算模式 | 数值计算模式 | 用途 |
|---------|-------------|-------------|------|
| `RT` | `GiNaC::ex` | `Rational<FlintMod>` | 有理函数类型 |
| `PT` | `GiNaC::ex` | `Polynomial<FlintMod>` | 多项式类型 |
| `ST` | `GiNaC::ex` | `FlintMod` | 标量类型 |

## 3. 代码架构

### 3.1 类层次结构

```
FlintMod          ← 有限域元素（底层数值类型）
    ↓
Polynomial<T>     ← 二变量多项式
    ↓
Rational<T>       ← 有理函数（分子/分母多项式）
    ↓
Series<T>         ← 二维幂级数
    
Sector<RT, PT>    ← 单个sector的处理（S矩阵、RREF、C和z）
    ↓
Family<RT, PT, ST> ← 管理所有sector，提供全局信息
    ↓
Redefinition<PT, ST> ← FBI重定义参数（pow_U、幂差计算）
    ↓
SeriesSolver<RT, PT, ST> ← 级数求解器（重定义后的约化+微分方程）
    ↓
IntegrandExpander<RT, PT, ST> ← 构造 Gν = P·Ĩ 二维级数
    ↓
IntegralTag / TargetConfig ← 描述 FI/BFI/BBFI 目标对象
    ↓
DeltaProjector<ST> ← 投影得到 FI/BFI/BBFI 一维 delta 级数
    ↓
runFI1DSeriesPipeline(...) ← 顶层流程
```

### 3.2 数据流

```
输入：topS矩阵 → Family构造 → Sector枚举/Case分类
         ↓
    转换为FlintMod类型
         ↓
    IntegrandExpander 构造 → 生成 shiftedU
         ↓
    Redefinition 构造（L, D_in, shiftedU）
         ↓
    SeriesSolver.setRedefinition()
      → 预计算 dRdX·U^L, dRdY·U^L, masterPowU_
         ↓
    IntegrandExpander.getIntegrand2DSeries(nu)
      → SeriesSolver.getFBISeries(nu, delta_in, targetDeg)
      → 按需递归（重定义后微分方程+约化）
         ↓
    Gν = P·Ĩ（多项式×级数，无需U^gamma）
         ↓
    DeltaProjector.project(Gν, tag) → FI/BFI/BBFI 一维 delta 级数输出
```

### 3.3 调用关系

```
IntegrandExpander::getIntegrand2DSeries(nu_target)
    └── SeriesSolver::getFBISeries(nu_target, delta_in, needDeg=targetDeg)
            ├── cache命中且阶数充足: 直接返回
            └── cache缺失或阶数不足:
                └── for deg = cachedDeg+1 .. needDeg:
                    ├── 若 (nu,delta) 是 master key:
                    │   └── solveMasterAtDeg(masterIdx, deg)
                    │       ├── deg=0: 从 masterBoundary_ 写入 (0,0)
                    │       └── deg>0:
                    │           ├── solveMasterCoeffX(masterIdx,p,q)
                    │           │   └── getFBISeries(nu+e_i+e_j, delta+1, p+q-1)
                    │           └── solveMasterCoeffY(masterIdx,q)
                    │               └── getFBISeries(nu+e_i+e_j, delta+1, q-1)
                    └── 否则:
                        └── reduceFBIAtDeg(result, nu, delta, deg)
                            ├── reduceCase0AtDeg
                            │   ├── !corner → case0IBPAtDeg
                            │   ├── corner 且 delta==targetDelta:
                            │   │   ├── 必要时 getFBISeries(masterNu, masterDelta, deg) 兜底补阶
                            │   │   └── 从 master cache 拷贝该 deg 系数到 result
                            │   └── corner 且 delta!=targetDelta:
                            │       ├── case0DimShiftUpAtDeg 或
                            │       └── case0DimShiftDownAtDeg
                            ├── reduceCase1AtDeg
                            ├── reduceCase2AtDeg
                            └── reduceCase3AtDeg
                                └── solveLRRAtDeg(g, D, polys, seriesPtrs, deg)
```

## 4. 核心类详解

### 4.1 FlintMod 类

**文件**：`include/ff_type.hpp`，`src/ff_type.tpp`

**功能**：封装FLINT库的有限域元素类型，提供 $\mathbb{Z}_p$ 中的算术运算。

```cpp
class FlintMod {
    static nmod_t mod_ctx;
    mp_limb_t value;
public:
    static void set_modulus(mp_limb_t p);
    template<typename IntType> FlintMod(IntType val);
    mp_limb_t get_value() const;
    // 算术运算（+, -, *, /）
};
```

### 4.2 Polynomial 类模板

**文件**：`include/rational.hpp`，`src/rational.tpp`

**功能**：表示二变量多项式 $P(X, Y) = \sum_{i,j} c_{ij} X^i Y^j$（稀疏存储）。

```cpp
template<typename T>
class Polynomial {
    std::unordered_map<Power, T, PowerHash> poly_;
    int deg_;
    size_t numOfMono_;
public:
    void addMonomial(const T& coeff, const Power& power);
    T getCoeff(int x, int y) const;
    Polynomial<T> derivativeX() const;
    Polynomial<T> derivativeY() const;
    auto begin() const; auto end() const;
};
```

### 4.3 Series 类模板

**文件**：`include/series.hpp`，`src/series.tpp`

**功能**：表示二维幂级数 $f(X, Y) = \sum_{i+j \leq d} c_{ij} X^i Y^j$。

```cpp
template<typename T>
class Series {
    std::vector<T> coefficients_;
    int deg_;
public:
    T getCoeff(int i, int j) const;
    void setCoeff(int i, int j, const T& coeff);
    static void mulPoly(Series<T>& result, const Series<T>& series, const Polynomial<T>& poly);
    static void divPoly(Series<T>& result, const Series<T>& series, const Polynomial<T>& poly);
};
```

### 4.4 Sector 类模板

**文件**：`include/sector.hpp`，`src/sector.tpp`

**功能**：处理单个sector的子矩阵 $S_{\text{sub}}$，执行RREF行化简，计算C和z系数，判定Case类型（0/1/2/3）。

提供分子/分母形式的系数（`denoCandZ_`, `numeC_`, `numeZ_`, `denoInvS_`, `numeInvS_`）供约化使用。

### 4.5 Family 类模板

**文件**：`include/family.hpp`，`src/family.tpp`

**功能**：管理整个FBI族，枚举所有有效sector，提供全局信息（Case分类、主积分列表、$\partial R / \partial X$、$\partial R / \partial Y$）。

### 4.6 Redefinition 结构体

**文件**：`include/series_solver.hpp`（与 SeriesSolver 同文件）

**功能**：封装FBI重定义参数，提供幂差计算。

```cpp
template<typename PT, typename ST>
struct Redefinition {
    ST D_in;           // 主积分 delta = L * D_Feynman / 2
    int L;             // 圈数
    PT shiftedU;       // 平移后的 U(X,Y) 多项式
    PT dUdX;           // ∂U/∂X
    PT dUdY;           // ∂U/∂Y

    Redefinition(int L_, const ST& D_in_, const PT& U_);

    bool sameParity(const ST& D) const;
    int etaMinus(const ST& D) const;   // {0, 2}
    int etaPlus(const ST& D) const;    // {-2, 0}

    /// 计算幂差 Δp_U = (ν_t_tot - ν_s_tot) - (L+1)/2 · (D̄_t - D̄_s)
    int deltaPowU(int nuTotT, const ST& DT, int nuTotS, const ST& DS) const;

    /// pow_U 标量值（微分方程 dlog 系数）
    ST powUScalar(int nuTot, const ST& D) const;
};
```

### 4.7 SeriesSolver 类模板

**文件**：`include/series_solver.hpp`，`src/series_solver.tpp`

**功能**：级数求解器，处理重定义后的微分方程和约化。

```cpp
template<typename RT, typename PT, typename ST>
class SeriesSolver {
    Family<RT, PT, ST>& family_;
    int targetDeg_, numMaster_, numProps_, numBranch_;
    
    // 主积分信息
    std::vector<std::vector<int>> masterNus_;
    std::vector<ST> masterDeltas_, masterBoundary_;
    
    // FBI级数缓存
    std::map<CacheKey, Series<ST>> cache_;
    std::map<CacheKey, int> cacheCurrentDeg_;
    
    // 重定义相关
    const Redefinition<PT, ST>* redef_;
    std::vector<std::vector<PT>> dRdXModified_;  // dRdX · U^L
    std::vector<std::vector<PT>> dRdYModified_;  // dRdY · U^L
    std::vector<ST> masterPowU_;                  // 每个主积分的 pow_U

public:
    void setRedefinition(const Redefinition<PT, ST>* redef);  // 必须在getFBISeries前调用
    const Series<ST>& getFBISeries(const std::vector<int>& nu, const ST& delta, int needDeg);
    
    static void solveLRRAtDeg(Series<ST>& g, const PT& D,
                               const std::vector<PT>& polys,
                               const std::vector<const Series<ST>*>& series, int deg);
    static PT multiplyPolys(const PT& a, const PT& b);
    static PT powPolyExpand(const PT& base, int exp);

private:
    void solveMasterAtDeg(int masterIdx, int deg);            // 按需补单个master某一阶
    void solveMasterCoeffX(int masterIdx, int p, int q);  // 重定义后DE
    void solveMasterCoeffY(int masterIdx, int q);          // 重定义后DE
    void applyRatioFactors(PT& D, std::vector<PT>& polys,
                           const std::vector<int>& deltaPs) const;
    // Case 0-3 约化（均调用 applyRatioFactors）
};
```

#### 实现要点

1. **`setRedefinition()`**：预计算 `dRdXModified_[i][j] = dRdX[i][j] · U^L`（避免运行时重复乘）、`masterPowU_[k]`。

2. **按需补阶主循环（`getFBISeries`）**：`getFBISeries(nu,delta,needDeg)` 会把该 key 从已缓存阶数逐层补到 `needDeg`。若 key 是 master（`delta==masterDelta`），每层调用 `solveMasterAtDeg`；否则每层调用 `reduceFBIAtDeg`。

3. **主积分补阶（`solveMasterAtDeg`）**：`deg=0` 直接读取 `masterBoundary_`；`deg>0` 时逐个 `(p,q)` 调用 `solveMasterCoeffX/Y` 计算该层系数。

4. **微分方程（`solveMasterCoeffX/Y`）**：使用 `dRdXModified_`/`dRdYModified_`（即 `dRdX·U^L`/`dRdY·U^L`），并按当前 $\nu$ 所在 sector 的活跃传播子索引取子块求和（等价于使用该 sector 子 family 的 $R$ 子矩阵）；额外计算 `dlogCoeff`（$\text{pow}_U \cdot [\partial U \cdot \widetilde{I}]$）和 `lhsCorrection`（$U$ 非常数项修正），递推公式参见 [`problem_solution.md`](./problem_solution.md) §6.4。

5. **约化函数（`case0IBP/DimShift/Case1-3`）**：构造 LRR 后调用 `applyRatioFactors(D, polys, deltaPs)`，该函数将 $U^{m}$ 乘到 $D$、$U^{\Delta p_i + m}$ 乘到各 $N_i$，统一为非负幂次的多项式。

### 4.8 IntegralTag 与 TargetConfig

**文件**：`include/integral_tag.hpp`，`src/io.tpp`

**功能**：描述 target 文件中的积分对象。每个 target 由 head、边界信息和传播子指数 $\vec{\nu}$ 组成。

```cpp
enum class IntegralHead {
    FI,      // two integrations
    BFI,     // one boundary, one integration
    BBFI     // two boundaries, no integration
};

enum class BoundaryAxis { X, Y };
enum class BoundarySide { U, D };

struct BoundaryTag {
    BoundaryAxis axis;
    BoundarySide side;
};

struct IntegralTag {
    IntegralHead head;
    std::vector<BoundaryTag> boundaries;  // FI:0, BFI:1, BBFI:2
    std::vector<int> nu;
};

struct TargetConfig {
    std::vector<IntegralTag> targets;
};
```

Target 文件格式支持显式 head：

```text
FI{1,1,1}
BFI[XU]{1,1,1}
BFI[XD]{1,1,1}
BFI[YU]{1,1,1}
BFI[YD]{1,1,1}
BBFI[XU,YU]{1,1,1}
BBFI[XU,YD]{1,1,1}
BBFI[XD,YU]{1,1,1}
BBFI[XD,YD]{1,1,1}
```

为兼容旧输入，只有 `{...}` 的行解析为 `FI{...}`。`nu` 始终只表示传播子指数，head 和 boundary 不编码进 `nu`。

### 4.9 IntegrandExpander 类模板

**文件**：`include/integrand_expander.hpp`，`src/integrand_expander.tpp`

**功能**：构造统一二维级数 $G_{\nu}(X,Y) = P(X,Y) \cdot \widetilde{I}_\nu^{D_{in}}(X,Y)$。FI、BFI、BBFI 对同一个 $\nu$ 共用该二维级数。

```cpp
template<typename RT, typename PT, typename ST>
class IntegrandExpander {
    SeriesSolver<RT, PT, ST>& solver_;
    int numLoops_, targetDeg_;
    ST feynmanD_, fbiDelta_, shiftA_, shiftB_;
    PT shiftedU_;

public:
    IntegrandExpander(SeriesSolver<RT, PT, ST>& solver,
                      int numLoops, int targetDeg,
                      const ST& feynmanD, const ST& shiftA, const ST& shiftB);
    
    /// Gν = J·W · Ĩ（多项式×级数，无需U^gamma）
    Series<ST> getIntegrand2DSeries(const std::vector<int>& nu) const;
    const PT& getShiftedU() const;

private:
    PT buildFIPolynomial(const std::vector<int>& nu) const;  // J·W（不含U）
    PT buildShiftedU() const;
    PT applyShift(const PT& xrYrPoly) const;
};
```

#### 实现要点

- `buildFIPolynomial(nu)`：在 $(X_r, Y_r)$ 上构造 $J \cdot X_0^{e_X} \cdot Y_0^{e_Y} \cdot Z_0^{e_Z}$，然后 `applyShift` 到平移坐标。不包含 $U^{\nu_{tot}}$（已被吸收进 $\widetilde{I}$）。
- `getIntegrand2DSeries(nu)`：调用 `solver_.getFBISeries(nu, fbiDelta_, targetDeg_)` 获取 $\widetilde{I}$，然后用 `Series::mulPoly` 乘以多项式 $P$。
- 顶层 pipeline 按 `nu` 缓存该二维级数；同一个 `nu` 的 FI、四个 BFI、四个 BBFI 不重复触发 `getFBISeries` 和 `mulPoly`。

### 4.10 DeltaProjector 类模板

**文件**：`include/delta_projector.hpp`，`src/delta_projector.tpp`

**功能**：将二维级数 $G_{\nu}(X,Y)$ 投影为一维 $\delta$ 级数。投影规则由 `IntegralTag` 决定。

```cpp
template<typename ST>
struct DeltaProjectionConfig {
    ST shiftA;
    ST shiftB;
    int degree;
};

template<typename ST>
class DeltaProjector {
public:
    explicit DeltaProjector(const DeltaProjectionConfig<ST>& config);
    std::vector<ST> project(const Series<ST>& series, const IntegralTag& tag) const;

private:
    std::vector<ST> intX_, intY_;
    std::vector<ST> xUpper_, xLower_, yUpper_, yLower_;
};
```

#### 实现要点

- 构造时预计算所有 $0 \leq k \leq \deg$ 的边界幂和积分权重：
  - $(-a)^k$、$(1-a)^k$、$(-b)^k$、$(1-b)^k$
  - $\frac{(1-a)^{k+1}-(-a)^{k+1}}{k+1}$
  - $\frac{(1-b)^{k+1}-(-b)^{k+1}}{k+1}$
- `project` 只遍历二维系数一次，根据 head 选择权重并写入对应的 $\delta$ 阶：
  - FI 写到 $p+q+2$
  - BFI 写到 $p+q+1$
  - BBFI 写到 $p+q$
- 若二维级数计算到总度数 `deg`，投影后的一维级数只保留 `{c0,...,c_deg}`。FI 的前两个零和 BFI 的前一个零显式保留；超过 `deg` 的项不由当前二维展开确定，因此不输出。

### 4.11 IO 与 Pipeline

**IO**：`include/io.hpp` — 解析 S、config、target 三类输入文件。target 解析为 `IntegralTag` 列表，并兼容旧的 `{nu}` 语法。

**Pipeline**：`include/fi_pipeline.hpp`，`src/fi_pipeline.tpp` — 顶层入口 `runFI1DSeriesPipeline()`，串联所有步骤：

```
Family 构造
  → IntegrandExpander
  → Redefinition
  → solver.setRedefinition()
  → 按 target 顺序处理 IntegralTag
      → 按 nu 查询/构造 Gν 二维级数 cache
      → DeltaProjector.project(Gν, tag)
      → 写出 {c0,...,c_deg}
```

## 5. 参考文献

[1] Li-Hong Huang, Rui-Jun Huang, and Yan-Qing Ma, "Tame multi-leg Feynman integrals beyond one loop," arXiv:2412.21053v1 [hep-ph], 2024.
