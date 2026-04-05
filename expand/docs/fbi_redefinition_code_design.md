# FBI 重定义：代码结构与实现文档

## 1. 概述

本文档描述 FBI 重定义功能的代码实现。数学公式请参考 [`fbi_redefinition_design.md`](./fbi_redefinition_design.md)。

### 1.1 改动范围

```
expand/
├── include/
│   ├── rational.hpp          ← 重构：Polynomial → MonomialSum + 新 Polynomial
│   ├── redefinition.hpp      ← 新增
│   ├── series_solver.hpp     ← 修改：添加 redef_ 成员
│   ├── integrand_expander.hpp← 修改：添加重定义路径
│   └── series.hpp            ← 无变化（接口名不变）
├── src/
│   ├── rational.tpp          ← 对应重构
│   ├── redefinition.tpp      ← 新增
│   ├── series_solver.tpp     ← 修改：约化函数 + 微分方程
│   └── integrand_expander.tpp← 修改：新 getFI2DSeries 路径
```

### 1.2 向后兼容

所有改动通过 `const Redefinition<ST,PT>* redef_ = nullptr` 控制。`redef_ == nullptr` 时走原有路径，行为完全不变。

## 2. 多项式类型升级

### 2.1 MonomialSum 类模板

将当前 `Polynomial<T>` 的全部实现**改名**为 `MonomialSum<T>`。接口不变：

```cpp
template<typename T>
class MonomialSum {
private:
    std::unordered_map<Power, T, PowerHash> poly_;
    int deg_;
    size_t numOfMono_;

public:
    MonomialSum();
    void addMonomial(const T& coeff, const Power& power);
    T getCoeff(int x, int y) const;
    int getDegree() const;
    size_t getNumOfMonomials() const;

    MonomialSum<T> derivativeX() const;
    MonomialSum<T> derivativeY() const;

    MonomialSum<T> operator+(const MonomialSum<T>& other) const;
    MonomialSum<T>& operator+=(const MonomialSum<T>& other);
    MonomialSum<T> operator*(const T& scalar) const;
    MonomialSum<T>& operator*=(const T& scalar);

    auto begin() const -> decltype(poly_.begin());
    auto end() const -> decltype(poly_.end());
};
```

### 2.2 新 Polynomial 类模板

新的 `Polynomial<T>` 表示若干 `MonomialSum<T>` 的乘积：

$$
P = \prod_{k=0}^{n-1} F_k
$$

其中每个 $F_k$ 是一个 `MonomialSum<T>`。

```cpp
template<typename T>
class Polynomial {
private:
    std::vector<MonomialSum<T>> factors_;  // 乘积因子列表

public:
    // ==================== 构造 ====================
    Polynomial();                                    // 零多项式（factors_ 为空）
    Polynomial(const MonomialSum<T>& ms);            // 单因子
    static Polynomial<T> one();                      // 常数 1

    // ==================== 乘法（拼接 factors，不展开） ====================
    Polynomial<T> operator*(const Polynomial<T>& other) const;
    Polynomial<T>& operator*=(const Polynomial<T>& other);

    // ==================== 标量乘法 ====================
    Polynomial<T> operator*(const T& scalar) const;
    Polynomial<T>& operator*=(const T& scalar);

    // ==================== 加法（展开后相加） ====================
    Polynomial<T> operator+(const Polynomial<T>& other) const;
    Polynomial<T>& operator+=(const Polynomial<T>& other);

    // ==================== 展开为单个 MonomialSum ====================
    MonomialSum<T> expand() const;

    // ==================== 向后兼容接口（触发展开） ====================
    T getCoeff(int x, int y) const { return expand().getCoeff(x, y); }
    int getDegree() const { return expand().getDegree(); }
    size_t getNumOfMonomials() const { return expand().getNumOfMonomials(); }
    bool isEmpty() const;
    Polynomial<T> derivativeX() const;
    Polynomial<T> derivativeY() const;

    auto begin() const -> decltype(expand().begin());
    auto end() const -> decltype(expand().end());

    // ==================== 因子访问 ====================
    int numFactors() const { return factors_.size(); }
    const std::vector<MonomialSum<T>>& getFactors() const { return factors_; }
};
```

### 2.3 关键语义

| 操作 | 行为 | 展开？ |
|:---|:---|:---:|
| `Polynomial * Polynomial` | 拼接两个 factors_ 列表 | 否 |
| `Polynomial * scalar` | 在第一个 factor 上乘标量 | 否 |
| `Polynomial + Polynomial` | 各自 expand()，再相加，结果是单因子 | 是 |
| `derivativeX/Y` | 先 expand()，再求导 | 是 |
| `begin()/end()/getCoeff()` | 先 expand() | 是 |

### 2.4 expand() 实现

```cpp
template<typename T>
MonomialSum<T> Polynomial<T>::expand() const {
    if (factors_.empty()) return MonomialSum<T>();  // 零多项式
    MonomialSum<T> result = factors_[0];
    for (size_t i = 1; i < factors_.size(); ++i) {
        MonomialSum<T> temp;
        for (const auto& [pa, ca] : result) {
            for (const auto& [pb, cb] : factors_[i]) {
                temp.addMonomial(ca * cb,
                    Power(pa.x_power + pb.x_power, pa.y_power + pb.y_power));
            }
        }
        result = std::move(temp);
    }
    return result;
}
```

### 2.5 Rational 的变化

`Rational<T>` 的分子和分母类型从旧 `Polynomial<T>` 变为 `MonomialSum<T>`（保持展开态）：

```cpp
template<typename T>
struct Rational {
    MonomialSum<T> numerator;
    MonomialSum<T> denominator;
    // ... 其余不变
};
```

> **原因**：`Rational` 的加法需要交叉乘，只有展开态才能做。且 Sector 中预计算的 `denoCandZ_`、`numeC_` 等都是展开态，不需要因子化。

### 2.6 对 Series 的影响

`Series::mulPoly` 的签名改为接受新 `Polynomial<T>`。内部实现：

```cpp
template<typename T>
void Series<T>::mulPoly(Series<T>& result, const Series<T>& series,
                         const Polynomial<T>& poly) {
    if (poly.numFactors() <= 1) {
        // 单因子：走原有展开态路径
        mulMonomialSum(result, series, poly.numFactors() == 0 ?
                       MonomialSum<T>() : poly.getFactors()[0]);
    } else {
        // 多因子：链式乘法
        // series * factor[n-1] * factor[n-2] * ... * factor[0]
        Series<T> temp = series;  // 拷贝
        for (int k = poly.numFactors() - 1; k >= 1; --k) {
            Series<T> next(temp.getDeg());
            mulMonomialSum(next, temp, poly.getFactors()[k]);
            temp = std::move(next);
        }
        mulMonomialSum(result, temp, poly.getFactors()[0]);
    }
}
```

其中 `mulMonomialSum` 是原来的 `mulPoly` 实现（改名）。

### 2.7 对 solveLRRAtDeg 的影响

`solveLRRAtDeg` 中调用 `polySeriesCoeff(poly, series, p, q)` 和 `sumPolySeriesCoeff`。当 `poly` 是多因子时：

**方案**：在 `solveLRRAtDeg` 的开头，对每个 `(polys[i], series[i])` 预处理一次：
- 如果 `polys[i]` 只有 1 个 factor → 不变
- 如果 `polys[i]` 有多个 factor → 将 series[i] 预乘以 `factor[1] * factor[2] * ...`，只保留 `factor[0]` 作为 poly

这样主循环中的 `polySeriesCoeff` 只处理单因子 MonomialSum，性能不退化。

```cpp
template<typename RT, typename PT, typename ST>
void SeriesSolver<RT, PT, ST>::solveLRRAtDeg(
    Series<ST>& g, const PT& D,
    const std::vector<PT>& polys,
    const std::vector<const Series<ST>*>& series,
    int deg)
{
    // 预处理：将多因子 poly 的额外因子预乘到 series 上
    std::vector<Series<ST>> tempSeriesStorage;
    std::vector<MonomialSum<ST>> effectivePolys;
    std::vector<const Series<ST>*> effectiveSeries;

    for (size_t i = 0; i < polys.size(); ++i) {
        if (polys[i].numFactors() <= 1) {
            effectivePolys.push_back(
                polys[i].numFactors() == 0 ? MonomialSum<ST>() : polys[i].getFactors()[0]);
            effectiveSeries.push_back(series[i]);
        } else {
            // 链式预乘：series * factor[n-1] * ... * factor[1]
            Series<ST> temp = *series[i];
            for (int k = polys[i].numFactors() - 1; k >= 1; --k) {
                Series<ST> next(temp.getDeg());
                Series<ST>::mulMonomialSum(next, temp, polys[i].getFactors()[k]);
                temp = std::move(next);
            }
            tempSeriesStorage.push_back(std::move(temp));
            effectivePolys.push_back(polys[i].getFactors()[0]);
            effectiveSeries.push_back(&tempSeriesStorage.back());
        }
    }

    // D 也可能是多因子，需要展开（D 只有一个，展开代价可接受）
    MonomialSum<ST> Dexpanded = D.expand();

    // 以下与原有逻辑完全相同，只是用 effectivePolys 和 effectiveSeries
    ST D00 = Dexpanded.getCoeff(0, 0);
    // ... 原有 LRR 求解逻辑 ...
}
```

## 3. Redefinition 数据结构

**文件**：`include/redefinition.hpp`，`src/redefinition.tpp`

```cpp
enum class RedefScheme { U, UJXYZ, UIrrat, None };

template<typename ST, typename PT>
struct Redefinition {
    ST D_in;                        // 主积分 delta
    int L;                          // 圈数
    RedefScheme scheme;             // 当前方案
    MonomialSum<ST> shiftedU;       // 平移后的 U(X,Y)（展开态）
    MonomialSum<ST> dUdX;           // ∂U/∂X（展开态）
    MonomialSum<ST> dUdY;           // ∂U/∂Y（展开态）
    std::vector<int> branchOfProp;  // 每个传播子的分支索引

    // ==================== 核心方法 ====================

    /// 判断 D 与 D_in 是否同奇偶
    bool sameParity(const ST& D) const;

    /// 计算 η_-(D) ∈ {0, 2}
    int etaMinus(const ST& D) const;

    /// 计算 η_+(D) ∈ {-2, 0}
    int etaPlus(const ST& D) const;

    /// 计算幂差 Δp_U（整数）
    /// 公式：Δp = (ν_t_tot - ν_s_tot) - (L+1)/2 * (D̄_t - D̄_s)
    int deltaPowU(const std::vector<int>& nu_t, const ST& D_t,
                  const std::vector<int>& nu_s, const ST& D_s) const;

    /// pow_U 在 Z_p 中的标量值（微分方程 dlogR 系数）
    /// 公式：pow_U = nu_tot - (L+1)*D̄/(2)，在 Z_p 中计算
    ST powUScalar(const std::vector<int>& nu, const ST& D) const;

    // ==================== 辅助 ====================

    /// 从 Z_p 提取小整数偏移量
    static int toSignedInt(mp_limb_t val, mp_limb_t p);
};
```

### 3.1 sameParity 实现

```cpp
template<typename ST, typename PT>
bool Redefinition<ST, PT>::sameParity(const ST& D) const {
    mp_limb_t p = /* 从 FlintMod 获取模数 */;
    int offset = toSignedInt((D - D_in).get_value(), p);
    return (offset % 2 == 0);
}
```

### 3.2 deltaPowU 实现

```cpp
template<typename ST, typename PT>
int Redefinition<ST, PT>::deltaPowU(
    const std::vector<int>& nu_t, const ST& D_t,
    const std::vector<int>& nu_s, const ST& D_s) const
{
    int nuTotT = 0, nuTotS = 0;
    for (int v : nu_t) nuTotT += v;
    for (int v : nu_s) nuTotS += v;

    // D̄_t - D̄_s 的整数值
    bool parityT = sameParity(D_t);
    bool parityS = sameParity(D_s);
    mp_limb_t p = /* modulus */;
    int offsetT = toSignedInt((D_t - D_in).get_value(), p);
    int offsetS = toSignedInt((D_s - D_in).get_value(), p);
    int DbarT = offsetT + (parityT ? 0 : 1);  // D̄_t - D_in 的偏移
    int DbarS = offsetS + (parityS ? 0 : 1);  // D̄_s - D_in 的偏移
    int DbarDiff = DbarT - DbarS;  // 总是偶数

    return (nuTotT - nuTotS) - (L + 1) * DbarDiff / 2;
}
```

## 4. SeriesSolver 的改动

### 4.1 新增成员

```cpp
template<typename RT, typename PT, typename ST>
class SeriesSolver {
    // ... 原有成员 ...

    // 新增
    const Redefinition<ST, MonomialSum<ST>>* redef_ = nullptr;

public:
    void setRedefinition(const Redefinition<ST, MonomialSum<ST>>* redef) {
        redef_ = redef;
    }
};
```

### 4.2 约化函数的修改模板

所有 6 个约化函数（`case0IBPAtDeg`、`case0DimShiftUpAtDeg`、`case0DimShiftDownAtDeg`、`reduceCase1AtDeg`、`reduceCase2AtDeg`、`reduceCase3AtDeg`）遵循同一修改模板。

以 `reduceCase1AtDeg` 为例。当前结构：

```cpp
// 当前
PT D = denoC * factor;
std::vector<PT> polys;
std::vector<const Series<ST>*> seriesPtrs;

for (每个源项 alpha) {
    polys.push_back(sector->getNumeZ(iCur) * ST(-1));
    seriesPtrs.push_back(&getFBISeries(nuMinusEi, delta - ST(1)));
}

solveLRRAtDeg(result, D, polys, seriesPtrs, deg);
```

修改后：

```cpp
// 新增 ratio 因子逻辑
PT D = denoC * factor;
std::vector<PT> polys;
std::vector<const Series<ST>*> seriesPtrs;

int dpMin = 0;  // 收集最小 Δp

if (redef_) {
    // 第一遍：收集所有 Δp，找最小值
    for (每个源项 alpha) {
        int dp = redef_->deltaPowU(nu, delta, nuMinusEi, delta - ST(1));
        dpMin = std::min(dpMin, dp);
    }
    // 如果有负 Δp，将 U^{-dpMin} 乘到 D 上（通分）
    if (dpMin < 0) {
        MonomialSum<ST> Upow = powMonomialSum(redef_->shiftedU, -dpMin);
        D = D * PT(Upow);
    }
}

for (每个源项 alpha) {
    PT Ni = PT(sector->getNumeZ(iCur)) * ST(-1);

    if (redef_) {
        int dp = redef_->deltaPowU(nu, delta, nuMinusEi, delta - ST(1));
        int adjustedPow = dp - dpMin;  // >= 0
        // 将 U^adjustedPow 追加为 Ni 的因子
        for (int k = 0; k < adjustedPow; ++k) {
            Ni *= PT(redef_->shiftedU);  // 拼接因子，不展开
        }
    }

    polys.push_back(Ni);
    seriesPtrs.push_back(&getFBISeries(nuMinusEi, delta - ST(1)));
}

solveLRRAtDeg(result, D, polys, seriesPtrs, deg);
```

> **关键**：`Ni *= PT(redef_->shiftedU)` 只拼接 factors_ 列表，不做多项式展开。展开由 `solveLRRAtDeg` 内部的预处理步骤完成（链式 mulPoly）。

### 4.3 各约化函数的 Δp 公式

直接使用 [`fbi_redefinition_design.md`](./fbi_redefinition_design.md) 第3节的公式。每个函数需要的参数：

| 函数 | 目标 $(\nu_t, D_t)$ | 源 $(\nu_s, D_s)$ |
|:---|:---|:---|
| `case0IBPAtDeg` | `(nuPlus, delta)` | `(nu, delta-1)` 和 `(nuMinusEi, delta-1)` |
| `case0DimShiftUpAtDeg` | `(nu, delta)` | `(nu, delta+1)` 和 `(nuMinusEi, delta)` |
| `case0DimShiftDownAtDeg` | `(nu, delta)` | `(nu, delta-1)` 和 `(nuMinusEi, delta-1)` |
| `reduceCase1AtDeg` | `(nu, delta)` | `(nuMinusEi, delta-1)` |
| `reduceCase2AtDeg` | `(nu, delta)` | `(nuMinusEi, delta)` |
| `reduceCase3AtDeg` | `(nu, delta)` | `(nuShifted, delta)` |

### 4.4 微分方程的修改

`solveMasterCoeffX` 的修改。当前实现：

```cpp
// 当前
ST rhsCoeff = ST(0);
for (i, j) {
    rhsCoeff -= factor_ij * polySeriesCoeff(dRdX[i][j], source, p-1, q) / ST(2);
}
f_{p,q} = rhsCoeff / ST(p);
```

修改后（`redef_ != nullptr` 时）：

```cpp
// 新实现
// 1. 微分方程源项的 Δp
int dpDE = redef_->deltaPowU(nu, delta, nuShifted_ij, delta + ST(1));
int adjustedPow = dpDE + 1;  // +1 因为两侧乘了 U

// 2. 处理可能的负 adjustedPow
int dpMin = adjustedPow;  // 所有 (i,j) 的 adjustedPow 相同
// 如果 dpMin < 0，构造 U^{-dpMin} 作为额外因子

// 3. rhsCoeff：每个源项的 series 预乘 U^{adjustedPow - dpMin}
ST rhsCoeff = ST(0);
for (i, j) {
    // 获取源 series
    const Series<ST>& source = getFBISeries(nuShifted_ij, delta + ST(1));
    // 预乘 U^{adjustedPow - dpMin} 次
    Series<ST> modifiedSource = mulPolyPower(source, redef_->shiftedU, adjustedPow - dpMin);
    rhsCoeff -= factor_ij * polySeriesCoeff(dRdX[i][j], modifiedSource, p-1, q) / ST(2);
}
// 如果 dpMin < 0，rhsCoeff 中已经包含了公分母

// 4. dlogCoeff：自耦合项
ST powU = redef_->powUScalar(nu, delta);
ST dlogCoeff = powU * polySeriesCoeff(redef_->dUdX, masterSeries, p-1, q);

// 5. lhsCorrection：U 的非常数项修正
ST lhsCorrection = ST(0);
for ((a,b) in shiftedU, (a,b) != (0,0)) {
    if (p >= a && q >= b) {
        lhsCorrection += U_ab * ST(p - a) * f_{p-a, q-b};
    }
}

// 如果 dpMin < 0，dlogCoeff 和 lhsCorrection 也要乘 U^{-dpMin}
// 实现上：将 masterSeries 也预乘 U^{-dpMin}

// 6. 最终
f_{p,q} = (rhsCoeff + dlogCoeff - lhsCorrection) / (U_00 * ST(p));
// 如果 dpMin < 0，分母还要乘 U_00^{-dpMin} 的修正
```

> **注意**：微分方程中所有 $(i,j)$ 的 $\Delta p_{ij}$ 相同（都是 $-2 - \frac{L+1}{2}\eta_+(D)$），所以 `adjustedPow` 只需算一次。

**Y 方向**（`solveMasterCoeffY`）完全对称，将 `dRdX` → `dRdY`，`dUdX` → `dUdY`，`p` → `q`，`(p-1,q)` → `(0,q-1)`。

### 4.5 mulPolyPower 辅助

```cpp
/// 对 series 链式乘以 MonomialSum poly 共 power 次
/// power 必须 >= 0
template<typename ST>
static Series<ST> mulPolyPower(const Series<ST>& series,
                                const MonomialSum<ST>& poly, int power) {
    if (power == 0) return series;
    Series<ST> result(series.getDeg());
    Series<ST>::mulMonomialSum(result, series, poly);
    for (int k = 1; k < power; ++k) {
        Series<ST> tmp(series.getDeg());
        Series<ST>::mulMonomialSum(tmp, result, poly);
        result = std::move(tmp);
    }
    return result;
}
```

## 5. IntegrandExpander 的改动

### 5.1 新增成员

```cpp
template<typename RT, typename PT, typename ST>
class IntegrandExpander {
    // ... 原有成员 ...

    // 新增
    const Redefinition<ST, MonomialSum<ST>>* redef_ = nullptr;

    // 预计算的平移因子多项式（MonomialSum，各只有几个单项式）
    MonomialSum<ST> shiftedX0_;   // X0 在平移坐标下
    MonomialSum<ST> shiftedY0_;   // Y0 在平移坐标下
    MonomialSum<ST> shiftedZ0_;   // Z0 在平移坐标下
    MonomialSum<ST> shiftedJ_;    // J = 1 - X0 在平移坐标下
};
```

### 5.2 新的 getFI2DSeries

```cpp
template<typename RT, typename PT, typename ST>
Series<ST> IntegrandExpander<RT, PT, ST>::getFI2DSeries(
    const std::vector<int>& nu) const
{
    if (redef_ == nullptr) {
        // 原有路径，完全不变
        return getFI2DSeriesOriginal(nu);
    }

    // ===== 重定义路径（U-only, L=2）=====

    // 1. 获取重定义后的 FBI 级数
    const Series<ST>& fbi = solver_.getFBISeries(nu, fbiDelta_);

    // 2. 计算分支幂次
    const auto& branchIndices = solver_.getFamily().getBranchIndices();
    std::vector<int> nuBranch(3, 0);
    for (int i = 0; i < (int)nu.size(); ++i) {
        nuBranch[branchIndices[i]] += nu[i];
    }
    int eX = nuBranch[0] - 1;
    int eY = nuBranch[1] - 1;
    int eZ = nuBranch[2] - 1;

    // 3. 链式 mulPoly：fbi × J × X0^eX × Y0^eY × Z0^eZ
    Series<ST> result = mulPolyPower(fbi, shiftedJ_, 1);
    result = mulPolyPower(result, shiftedX0_, eX);
    result = mulPolyPower(result, shiftedY0_, eY);
    result = mulPolyPower(result, shiftedZ0_, eZ);

    return result;
    // 不需要 U^nuTot，不需要 U^gamma，不需要 multiplySeries
}
```

### 5.3 平移因子的构造

在构造函数中一次性计算（使用 `applyShift`）：

```
X0_raw = Xr                    → shiftedX0_ = applyShift(X0_raw)
Y0_raw = (1-Xr)*Yr             → shiftedY0_ = applyShift(Y0_raw)
Z0_raw = 1 - Xr - (1-Xr)*Yr   → shiftedZ0_ = applyShift(Z0_raw)
J_raw  = 1 - Xr                → shiftedJ_  = applyShift(J_raw)
```

每个平移后的 MonomialSum 只有约 2-4 个单项式。

## 6. 管线连接

在 `fi_pipeline.tpp`（或 runner 入口）中：

```cpp
// 1. 构造 Redefinition
Redefinition<FlintMod, MonomialSum<FlintMod>> redef;
redef.D_in = fbiDelta;  // = L * feynmanD / 2
redef.L = numLoops;
redef.scheme = RedefScheme::U;
redef.shiftedU = expander.getShiftedU();   // 已有的平移 U
redef.dUdX = redef.shiftedU.derivativeX();
redef.dUdY = redef.shiftedU.derivativeY();
redef.branchOfProp = family.getBranchIndices();

// 2. 传入
solver.setRedefinition(&redef);
expander.setRedefinition(&redef);

// 3. 正常 solve + 输出
solver.solve();
```

## 7. 实施顺序

| 步骤 | 内容 | 验证 |
|:---:|:---|:---|
| 1 | `MonomialSum` / 新 `Polynomial` / `Rational` 重构 | 原有测试全部通过 |
| 2 | `solveLRRAtDeg` 适配新 `Polynomial`（预处理多因子） | 原有测试全部通过 |
| 3 | 新建 `redefinition.hpp/tpp` | 单元测试 |
| 4 | `SeriesSolver` 添加 `redef_`，修改 6 个约化函数 | 对比 redef=None |
| 5 | `SeriesSolver` 修改微分方程 | 对比 redef=None |
| 6 | `IntegrandExpander` 新路径 | 对比 redef=None |
| 7 | 管线连接 | 端到端回归 |

每一步保持向后兼容，可逐步验证。
