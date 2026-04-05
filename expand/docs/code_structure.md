# 代码结构与实现文档

## 1. 概述

本文档详细描述项目的代码结构和实现细节。完整的数学理论背景请参考 [`problem_and_workflow.md`](./problem_and_workflow.md)。

### 1.1 项目目标

在质数域 $\mathbb{Z}_p$ 下进行端到端计算：
- 重定义后 FBI（$\widetilde{I} = U^{\text{pow}_U} \cdot I$）的二维幂级数展开（逐阶递推）
- FI 被积函数二维级数构造（$P \cdot \widetilde{I}$，无需 $U^\gamma$ 级数）
- FI 的一维级数输出（对二维级数积分）

### 1.2 核心方法

将IBP约化和微分方程统一为**线性递推关系**（LRR）形式，直接计算级数系数，而非先计算符号有理函数再做级数运算。约化和微分方程中通过 `Redefinition` 结构体处理 $U$ 幂差因子。

## 2. 代码规范

### 2.1 文件组织

```
expand/
├── docs/                          # 文档
│   ├── problem_and_workflow.md    # 问题定义与计算流程
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
│   ├── integrand_expander.hpp     # FI被积函数二维级数构造
│   ├── series_integrator.hpp      # 二维级数到一维级数积分
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
IntegrandExpander<RT, PT, ST> ← 构造 FI = P·Ĩ 二维级数
    ↓
SeriesIntegrator<ST> ← 积分得到 FI 一维级数
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
    SeriesSolver.solve() → 逐阶递推（重定义后微分方程+约化）
         ↓
    IntegrandExpander.getFI2DSeries(nu)
      → P·Ĩ（多项式×级数，无需U^gamma）
         ↓
    SeriesIntegrator.integrate() → 一维级数输出
```

### 3.3 调用关系

```
SeriesSolver::solve()
    │
    ├── solveAtDeg(deg)
    │       ├── solveMasterCoeffX(k, p, q)   [重定义后DE]
    │       │     ├── getFBISeries(nu+e_i+e_j, delta+1) → 触发约化
    │       │     ├── dlogCoeff: pow_U * [dUdX · Ĩ]
    │       │     └── lhsCorrection: U非常数项修正
    │       └── solveMasterCoeffY(k, q)       [重定义后DE]
    │
    └── getFBISeries(nu, delta) → reduceFBIAtDeg()
            ├── case0IBPAtDeg()        + applyRatioFactors()
            ├── case0DimShiftUp/Down() + applyRatioFactors()
            ├── reduceCase1AtDeg()     + applyRatioFactors()
            ├── reduceCase2AtDeg()     + applyRatioFactors()
            └── reduceCase3AtDeg()     [Δp=0, 无需ratio]
                    └── solveLRRAtDeg()
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
    int targetDeg_, currentDeg_, numMaster_, numProps_, numBranch_;
    
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
    void setRedefinition(const Redefinition<PT, ST>* redef);  // 必须在solve()前调用
    void solve();
    const Series<ST>& getFBISeries(const std::vector<int>& nu, const ST& delta);
    
    static void solveLRRAtDeg(Series<ST>& g, const PT& D,
                               const std::vector<PT>& polys,
                               const std::vector<const Series<ST>*>& series, int deg);
    static PT multiplyPolys(const PT& a, const PT& b);
    static PT powPolyExpand(const PT& base, int exp);

private:
    void solveMasterCoeffX(int masterIdx, int p, int q);  // 重定义后DE
    void solveMasterCoeffY(int masterIdx, int q);          // 重定义后DE
    void applyRatioFactors(PT& D, std::vector<PT>& polys,
                           const std::vector<int>& deltaPs) const;
    // Case 0-3 约化（均调用 applyRatioFactors）
};
```

#### 实现要点

1. **`setRedefinition()`**：预计算 `dRdXModified_[i][j] = dRdX[i][j] · U^L`（避免运行时重复乘）、`masterPowU_[k]`。

2. **微分方程（`solveMasterCoeffX/Y`）**：使用 `dRdXModified_` 代替原始 `dRdX`，额外计算 `dlogCoeff`（$\text{pow}_U \cdot [\partial U \cdot \widetilde{I}]$）和 `lhsCorrection`（$U$ 非常数项修正），递推公式参见 [`problem_and_workflow.md`](./problem_and_workflow.md) §6.4。

3. **约化函数（`case0IBP/DimShift/Case1-3`）**：构造 LRR 后调用 `applyRatioFactors(D, polys, deltaPs)`，该函数将 $U^{m}$ 乘到 $D$、$U^{\Delta p_i + m}$ 乘到各 $N_i$，统一为非负幂次的多项式。

### 4.8 IntegrandExpander 类模板

**文件**：`include/integrand_expander.hpp`，`src/integrand_expander.tpp`

**功能**：构造 $\text{FI}(X,Y) = P(X,Y) \cdot \widetilde{I}_\nu^{D_{in}}(X,Y)$ 的二维级数。

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
    
    /// FI = J·W · Ĩ（多项式×级数，无需U^gamma）
    Series<ST> getFI2DSeries(const std::vector<int>& nu) const;
    const PT& getShiftedU() const;

private:
    PT buildFIPolynomial(const std::vector<int>& nu) const;  // J·W（不含U）
    PT buildShiftedU() const;
    PT applyShift(const PT& xrYrPoly) const;
};
```

#### 实现要点

- `buildFIPolynomial(nu)`：在 $(X_r, Y_r)$ 上构造 $J \cdot X_0^{e_X} \cdot Y_0^{e_Y} \cdot Z_0^{e_Z}$，然后 `applyShift` 到平移坐标。不包含 $U^{\nu_{tot}}$（已被吸收进 $\widetilde{I}$）。
- `getFI2DSeries(nu)`：调用 `solver_.getFBISeries(nu, fbiDelta_)` 获取 $\widetilde{I}$，然后用 `Series::mulPoly` 乘以多项式 $P$。

### 4.9 SeriesIntegrator 类模板

**文件**：`include/series_integrator.hpp`，`src/series_integrator.tpp`

**功能**：将二维级数积分为一维级数（积分域由平移参数 $a, b$ 定义）。

### 4.10 IO 与 Pipeline

**IO**：`include/io.hpp` — 解析 S、config、target 三类输入文件。

**Pipeline**：`include/fi_pipeline.hpp`，`src/fi_pipeline.tpp` — 顶层入口 `runFI1DSeriesPipeline()`，串联所有步骤：Family构造 → IntegrandExpander → Redefinition → SeriesSolver.solve() → getFI2DSeries → integrate → 输出。

## 5. 参考文献

[1] Li-Hong Huang, Rui-Jun Huang, and Yan-Qing Ma, "Tame multi-leg Feynman integrals beyond one loop," arXiv:2412.21053v1 [hep-ph], 2024.
