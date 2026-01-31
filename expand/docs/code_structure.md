# 代码结构与实现文档

## 1. 概述

本文档详细描述项目的代码结构和实现细节。完整的数学理论背景请参考 [`problem_and_workflow.md`](./problem_and_workflow.md)。

### 1.1 项目目标

在质数域 $\mathbb{Z}_p$ 下计算二圈Feynman积分（FBI）的二维幂级数展开，使用**逐阶递推**的方法。

### 1.2 核心方法

将IBP约化和微分方程统一为**线性递推关系**（LRR）形式，直接计算级数系数，而非先计算符号有理函数再做级数运算。

## 2. 代码规范

### 2.1 文件组织

```
expand/
├── .clinerules/                   # Cline规则文件
│   ├── project-overview.md        # 项目概述
│   ├── coding-rules.md            # 代码规范
│   └── doc-sync-rules.md          # 文档同步规则
├── docs/                          # 文档
│   ├── problem_and_workflow.md    # 问题定义与计算流程
│   ├── code_structure.md          # 代码架构（本文档）
│   └── todo.md                    # 待办事项
├── include/                       # 头文件（类声明）
│   ├── ff_type.hpp                 # 有限域类型 FlintMod
│   ├── rational.hpp               # 多项式 Polynomial 和有理函数 Rational
│   ├── series.hpp                 # 二维幂级数 Series
│   ├── sector.hpp                 # Sector类（含辅助函数）
│   ├── family.hpp                 # Family类
│   ├── series_solver.hpp          # SeriesSolver类
│   └── converter.hpp              # GiNaC到FlintMod类型转换
├── src/                           # 实现文件（模板实现）
│   ├── ff_type.tpp
│   ├── rational.tpp
│   ├── series.tpp
│   ├── sector.tpp
│   ├── family.tpp
│   ├── series_solver.tpp
│   └── converter.tpp
└── test/                          # 测试代码
```

### 2.2 命名规范

- **文件名**：小写字母加下划线（snake_case），头文件用 `.hpp`，模板实现用 `.tpp`
- **类名**：大驼峰（PascalCase），如 `SeriesSolver`、`FlintMod`
- **成员变量**：小驼峰加下划线后缀，如 `numProps_`、`targetDeg_`
- **成员函数**：小驼峰（camelCase），如 `getCase()`、`solveLRRAtDeg()`
- **模板参数**：大写字母缩写，如 `RT`（有理函数类型）、`PT`（多项式类型）、`ST`（标量类型）
- **局部变量**：小驼峰，如 `nuSum`、`maxIndex`

### 2.3 模板设计

项目采用双/三参数模板设计，支持符号计算和数值计算两种模式：

| 模板参数 | 符号计算模式 | 数值计算模式 | 用途 |
|---------|-------------|-------------|------|
| `RT` | `GiNaC::ex` | `Rational<FlintMod>` | 有理函数类型 |
| `PT` | `GiNaC::ex` | `Polynomial<FlintMod>` | 多项式类型 |
| `ST` | `GiNaC::ex` | `FlintMod` | 标量类型（级数系数、delta等） |

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
SeriesSolver<RT, PT, ST> ← 级数求解器（核心算法）
```

### 3.2 数据流

```
输入：topS矩阵（GiNaC符号表达式）
         ↓
    [Family构造]
         ↓
    枚举所有有效Sector
    对每个Sector: RREF → 计算C,z → Case分类
         ↓
    [转换为FlintMod类型]
         ↓
    Family<Rational<FlintMod>, Polynomial<FlintMod>, FlintMod>
         ↓
    [SeriesSolver构造]
         ↓
    设置主积分边界条件
         ↓
    [solve() 逐阶递推]
         ↓
    for deg = 0 to targetDeg:
        主积分: 微分方程递推
        非主积分: IBP约化递推（按需触发）
         ↓
输出：所有主积分的级数展开
```

### 3.3 调用关系

```
SeriesSolver::solve()
    │
    ├── solveAtDeg(deg)
    │       │
    │       ├── solveMasterCoeffX(masterIdx, p, q)  [p > 0]
    │       │       │
    │       │       └── getFBISeries(nu+e_i+e_j, delta+1) → 触发约化
    │       │
    │       └── solveMasterCoeffY(masterIdx, q)      [p = 0, q > 0]
    │               │
    │               └── getFBISeries(nu+e_i+e_j, delta+1) → 触发约化
    │
    └── getFBISeries(nu, delta)  [递归调用]
            │
            └── reduceFBIAtDeg(result, nu, delta, deg)
                    │
                    ├── reduceCase0AtDeg() → case0IBPAtDeg() / case0DimShift...()
                    ├── reduceCase1AtDeg()
                    ├── reduceCase2AtDeg()
                    └── reduceCase3AtDeg()
                            │
                            └── solveLRRAtDeg() → polySeriesCoeff()
```

## 4. 核心类详解

### 4.1 FlintMod 类

**文件**：`include/ff_type.hpp`，`src/ff_type.tpp`

**功能**：封装FLINT库的有限域元素类型，提供 $\mathbb{Z}_p$ 中的算术运算。

```cpp
class FlintMod {
private:
    static nmod_t mod_ctx;           // FLINT模数上下文（静态，全局共享）
    static bool ctx_initialized;     // 上下文是否已初始化（静态）
    mp_limb_t value;                 // 元素值（0 到 p-1）

public:
    // ==================== 核心成员函数 ====================
    static void set_modulus(mp_limb_t p);  // 设置全局模数 p（只需调用一次）
    
    template<typename IntType>
    FlintMod(IntType val);           // 从整数构造（支持所有整数类型，自动取模）
    
    mp_limb_t get_value() const;     // 获取原始数值
    
    // 算术运算（模 p）
    FlintMod operator+(const FlintMod& other) const;
    // ...
};
```

### 4.2 Polynomial 类模板

**文件**：`include/rational.hpp`，`src/rational.tpp`

**功能**：表示二变量多项式 $P(X, Y) = \sum_{i,j} c_{ij} X^i Y^j$。

```cpp
// 辅助结构：单项式幂次
struct Power {
    int x_power;    // X的指数
    int y_power;    // Y的指数
};

template<typename T>
class Polynomial {
private:
    std::unordered_map<Power, T, PowerHash> poly_;  // 幂次 → 系数的映射（稀疏存储）
    int deg_;                                        // 多项式度（-1表示零多项式）
    size_t numOfMono_;                               // 非零单项式数量

public:
    // ==================== 核心成员函数 ====================
    void addMonomial(const T& coeff, const Power& power);  // 添加单项式（自动合并同类项）
    T getCoeff(int x, int y) const;                        // 获取指定幂次的系数
    int getDegree() const;                                 // 获取多项式度
    
    Polynomial<T> derivativeX() const;   // 对X求偏导
    Polynomial<T> derivativeY() const;   // 对Y求偏导
    
    // 迭代器（遍历单项式）
    auto begin() const { return poly_.begin(); }
    auto end() const { return poly_.end(); }
};
```

### 4.3 Series 类模板

**文件**：`include/series.hpp`，`src/series.tpp`

**功能**：表示二维幂级数 $f(X, Y) = \sum_{i+j \leq d} c_{ij} X^i Y^j$。

```cpp
template<typename T>
class Series {
private:
    std::vector<T> coefficients_;  // 系数数组，按特定顺序存储
    int deg_;                      // 级数的最高度数
    
    // 索引转换：(i,j) → 一维索引
    // 存储顺序: (0,0), (1,0), (0,1), (2,0), (1,1), (0,2), ...
    inline int getIndex(int i, int j) const {
        int deg = i + j;
        return (deg * (deg + 1)) / 2 + j;
    }

public:
    // ==================== 核心成员函数 ====================
    Series(int deg = 0);                          // 构造指定度数的零级数
    T getCoeff(int i, int j) const;               // 获取系数（越界返回0）
    void setCoeff(int i, int j, const T& coeff);  // 设置系数
    int getDeg() const;                           // 获取度数
    
    // 静态方法（避免临时对象）
    static void mulPoly(Series<T>& result, const Series<T>& series, const Polynomial<T>& poly);
    static void divPoly(Series<T>& result, const Series<T>& series, const Polynomial<T>& poly);
};
```

#### 实现要点：级数除法

级数除以多项式 $g = f / P$ 的递推公式（参见 [`problem_and_workflow.md`](./problem_and_workflow.md) 5.4节）：

$$
g_{pq} = \frac{1}{P_{00}} \left( f_{pq} - \sum_{(a,b) \neq (0,0)} P_{ab} \cdot g_{p-a,q-b} \right)
$$

### 4.4 Sector 类模板

**文件**：`include/sector.hpp`，`src/sector.tpp`

**功能**：处理单个sector的子矩阵 $S_{\text{sub}}$，执行RREF行化简，计算C和z系数。

```cpp
template<typename RT, typename PT>  // RT: 有理函数类型, PT: 多项式类型
class Sector {
private:
    // ==================== 来自构造参数 ====================
    std::vector<std::vector<PT>> S_;   // 输入的(B+N)×(B+N)子矩阵
    int numProps_;                      // 该sector中非零传播子的数量N
    int numBranch_;                     // 分支数B

    // ==================== RREF结果 ====================
    std::vector<std::vector<RT>> reducedS_;      // RREF后的矩阵
    std::vector<std::vector<RT>> rowOperation_;  // 行变换矩阵（dimNull=0时为S^{-1}）
    int dimNull_;                                 // 零空间维度
    int z0_;                                      // 0或1

    // ==================== C和z系数 ====================
    std::vector<RT> candz_;     // (C_1,...,C_B, z_1,...,z_N)
    RT C_;                      // C = sum(C_b)

    // ==================== 分子/分母形式（用于LRR） ====================
    PT denoCandZ_;                       // candz_所有元素分母的lcm
    std::vector<PT> numeCandZ_;          // numeCandZ_[i] = candz_[i] * denoCandZ_
    PT numeC_;                           // numeC_ = C_ * denoCandZ_
    std::vector<PT> denoRowOperation_;   // 每行的分母lcm
    std::vector<std::vector<PT>> numeRowOperation_;  // 分子形式

public:
    // ==================== 核心成员函数 ====================
    Sector(const std::vector<std::vector<PT>>& S, int numProps, int numBranch);
    
    int getCase() const;        // 返回Case类型（0/1/2/3）
    int getDimNull() const;
    int getZ0() const;
    
    // 获取分子/分母形式（用于约化）
    const PT& getDenoCandZ() const;
    PT getNumeZ(int i) const;   // z_i * denoCandZ
    const PT& getNumeC() const;
    const std::vector<PT>& getDenoInvS() const;
    const std::vector<std::vector<PT>>& getNumeInvS() const;

private:
    void rowReduce();     // 执行RREF行化简
    void solveCandZ();    // 求解C和z系数
};
```

#### 实现要点：RREF与Case判定

Sector类的核心功能是对输入的子矩阵 $S$ 执行RREF（简化行阶梯形）化简：
1. 通过RREF判断零空间维度 $\dim(\text{Null}(S))$
2. 若 $\dim(\text{Null}(S)) = 0$，求解 $S \cdot x = (1,\ldots,1,0,\ldots,0)^T$ 得到 $(C, z)$
3. 若 $\dim(\text{Null}(S)) > 0$，取零空间基向量作为 $(C, z)$
4. 根据 $\dim(\text{Null}(S))$ 和 $C = \sum C_b$ 是否为零，确定Case类型（0/1/2/3）

### 4.5 Family 类模板

**文件**：`include/family.hpp`，`src/family.tpp`

**功能**：管理整个FBI族，枚举所有有效sector，提供全局信息。

```cpp
template<typename RT, typename PT, typename ST>
class Family {
private:
    // ==================== 来自构造参数 ====================
    std::vector<std::vector<PT>> topS_;  // 顶层S矩阵 (N+B)×(N+B)
    int numProps_;                        // 传播子数N
    int numBranch_;                       // 分支数B

    // ==================== 分支信息 ====================
    std::vector<int> branchIndices_;      // 每个传播子对应的分支索引

    // ==================== Sector信息 ====================
    std::vector<Sector<RT, PT>> sectors_;  // 所有有效sector
    std::vector<int> sectorIdxs_;          // 每个sector对应的二进制索引
    std::vector<int> cases_;               // cases_[idx] = Case类型（无效为-1）
    std::vector<int> masterIdxs_;          // 主积分索引列表

    // ==================== 主积分信息 ====================
    std::vector<ST> masterDeltas_;         // 每个主积分的目标delta

    // ==================== R矩阵导数 ====================
    std::vector<std::vector<PT>> dRdX_;    // dR/dX
    std::vector<std::vector<PT>> dRdY_;    // dR/dY

public:
    // ==================== 核心成员函数 ====================
    template<typename Symbol>
    Family(const std::vector<std::vector<PT>>& topS, int numProps, int numBranch,
           const Symbol& X, const Symbol& Y);
    
    int getCase(std::vector<int> nu) const;
    const Sector<RT, PT>* getSector(std::vector<int> nu) const;
    
    const std::vector<std::vector<PT>>& getDRdX() const;
    const std::vector<std::vector<PT>>& getDRdY() const;
    const std::vector<int>& getMasterIdxs() const;
    
    void setMasterDelta(ST Delta);

private:
    void findSectors();   // 枚举所有有效sector
};
```

#### 实现要点：Sector枚举与构造

Family类根据输入的顶层矩阵 `topS` 构造所有有效的Sector：
1. 枚举所有 $2^N$ 个可能的sector（用二进制索引表示）
2. 筛选覆盖所有 $B$ 个分支的有效sector
3. 对每个有效sector，从 `topS` 提取对应的子矩阵，创建 `Sector` 对象
4. 记录每个sector的Case类型，识别所有主积分（Case 0的角积分）
5. 预计算 $\partial R / \partial X$ 和 $\partial R / \partial Y$ 用于微分方程

### 4.6 SeriesSolver 类模板

**文件**：`include/series_solver.hpp`，`src/series_solver.tpp`

**功能**：级数求解器，统一处理微分方程求解和IBP约化。

```cpp
template<typename RT, typename PT, typename ST>
class SeriesSolver {
private:
    // ==================== 来自构造参数 ====================
    Family<RT, PT, ST>& family_;    // Family对象引用
    int targetDeg_;                  // 目标级数度数
    int currentDeg_;                 // 当前已求解的度数
    int numMaster_;                  // 主积分数量
    int numProps_;                   // 传播子数量
    int numBranch_;                  // 分支数量
    
    // ==================== 主积分信息 ====================
    std::vector<std::vector<int>> masterNus_;     // 每个主积分的nu
    std::vector<ST> masterDeltas_;                // 每个主积分的delta
    std::vector<ST> masterBoundary_;              // 边界条件（零阶系数）
    
    // ==================== FBI级数缓存 ====================
    struct CacheKey {
        std::vector<int> nu;
        ST delta;
    };
    mutable std::map<CacheKey, Series<ST>> cache_;
    mutable std::map<CacheKey, int> cacheCurrentDeg_;

public:
    // ==================== 核心成员函数 ====================
    SeriesSolver(Family<RT, PT, ST>& family, int targetDeg);
    
    void setMasterBoundary(int masterIdx, const ST& value);
    void setAllMasterBoundary(const ST& value);
    void solve();   // 主求解函数
    
    const Series<ST>& getFBISeries(const std::vector<int>& nu, const ST& delta);
    const Series<ST>& getMasterSeries(int masterIdx) const;
    
    // LRR求解器（静态方法）
    static void solveLRRAtDeg(Series<ST>& g, const PT& D,
                               const std::vector<PT>& polys,
                               const std::vector<const Series<ST>*>& series, int deg);

private:
    void solveAtDeg(int deg);
    void solveMasterCoeffX(int masterIdx, int p, int q);
    void solveMasterCoeffY(int masterIdx, int q);
    
    // Case 0-3 约化
    void reduceCase0AtDeg(Series<ST>& result, const std::vector<int>& nu, const ST& delta, int deg);
    void reduceCase1AtDeg(Series<ST>& result, const std::vector<int>& nu, const ST& delta, int deg);
    void reduceCase2AtDeg(Series<ST>& result, const std::vector<int>& nu, const ST& delta, int deg);
    void reduceCase3AtDeg(Series<ST>& result, const std::vector<int>& nu, const ST& delta, int deg);
};
```

#### 实现要点：逐阶递推算法

SeriesSolver实现了 [`problem_and_workflow.md`](./problem_and_workflow.md) 第5节描述的逐阶递推算法：

1. **递推逻辑**：假设已知主积分的 $\deg \leq N$ 项和其它积分的 $\deg < N$ 项，则
   - 通过约化计算其它积分的 $\deg = N$ 项
   - 通过微分方程计算主积分的 $\deg = N+1$ 项

2. **统一的LRR形式**：所有约化公式和微分方程都转换为标准的线性递推关系：
   $$D(X,Y) \cdot g(X,Y) = \sum_i N_i(X,Y) \cdot f_i(X,Y)$$
   然后调用 `solveLRRAtDeg()` 统一求解第 $\deg$ 阶的系数。

3. **按需约化**：微分方程右边需要的FBI级数通过 `getFBISeries()` 递归获取，自动触发相应的约化计算。

## 5. 使用示例

```cpp
#include <ginac/ginac.h>
#include "family.hpp"
#include "series_solver.hpp"
#include "converter.hpp"

int main() {
    FlintMod::set_modulus(1000000007);
    
    GiNaC::symbol X("X"), Y("Y");
    // ... 构造topS矩阵 ...
    
    Family<GiNaC::ex, GiNaC::ex, GiNaC::ex> ginacFamily(topS, numProps, numBranch, X, Y);
    auto family = convertFamily(ginacFamily, X, Y);
    family.setMasterDelta(FlintMod(4));
    
    SeriesSolver<Rational<FlintMod>, Polynomial<FlintMod>, FlintMod> solver(family, 100);
    solver.setAllMasterBoundary(FlintMod(1));
    solver.solve();
    
    return 0;
}
```

## 6. 参考文献

[1] Li-Hong Huang, Rui-Jun Huang, and Yan-Qing Ma, "Tame multi-leg Feynman integrals beyond one loop," arXiv:2412.21053v1 [hep-ph], 2024.
