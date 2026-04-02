# FI 二维展开与积分模块实现说明

## 1. 范围

本文档只服务于当前代码实现，目标只有两个：

1. 生成 FI 被积函数的二维级数
2. 将二维级数积分成一维级数

当前只需要实现两个新类：

- `IntegrandExpander`
- `SeriesIntegrator`

## 2. 当前问题中的固定定义

### 2.1 变量替换

原始三变量满足：

$$
X + Y + Z = 1
$$

先做二维化替换：

$$
X \to X,\qquad
Y \to (1-X)Y,\qquad
Z \to 1 - X - (1-X)Y
$$

再做局部平移：

$$
X \to X + a,\qquad
Y \to Y + b
$$

实现时，所有要进入二维级数展开的多项式，都必须按这个顺序处理：

1. 三变量到二维变量
2. 二维变量平移

### 2.2 `U`

原始

$$
U = XY + YZ + ZX
$$

实现中使用的 `U(X,Y)` 定义为：

$$
U(X,Y)
=
\left.
\left.
(XY+YZ+ZX)
\right|_{X\to X,\;Y\to(1-X)Y,\;Z\to1-X-(1-X)Y}
\right|_{X\to X+a,\;Y\to Y+b}
$$

### 2.3 `gamma`

固定定义：

$$
\gamma = -\frac{(L+1)D}{2}
$$

当前两圈问题 `L = 2`，所以

$$
\gamma = -\frac{3D}{2}
$$

这里的 `D` 是费曼积分参数。FBI 里使用的是独立参数 `\Delta`，两者关系为：

$$
\Delta = \frac{L D}{2}
$$

因此当前两圈 `L=2` 时有 `\Delta = D`，但在实现和文档中仍保留两个参数名，不混用。

实现要求：

- `IntegrandExpander` 同时保存 `feynmanD_` 和 `fbiDelta_`
- 其中 `fbiDelta_` 在构造时由 `feynmanD_` 和 `numLoops_` 计算
- `getFI2DSeries` 内部调用 FBI 时使用 `fbiDelta_`

### 2.4 多项式因子

当前 FI 被积函数中需要处理的多项式因子只有：

1. Jacobian：
   $$
   1 - X
   $$

2. 参数幂次：
   $$
   X^{\nu_X-1}Y^{\nu_Y-1}Z^{\nu_Z-1}
   $$

这些因子也都要先做：

1. 三变量到二维变量替换
2. 局部平移

最终把它们合并成一个二维多项式：

$$
P(X,Y)
$$

## 3. 开发规格

### 3.1 新增文件

```text
expand/include/integrand_expander.hpp
expand/src/integrand_expander.tpp
expand/include/series_integrator.hpp
expand/src/series_integrator.tpp
```

### 3.2 IntegrandExpander（类声明 + 注释）

```cpp
template<typename RT, typename PT, typename ST>
class IntegrandExpander {
private:
    // dependency
    SeriesSolver<RT, PT, ST>& solver_;

    // config
    int numLoops_;          // L
    int targetDeg_;         // max total degree for 2D series
    ST feynmanD_;           // D in Feynman integral
    ST fbiDelta_;           // Delta = L*D/2 (parameter for FBI solver)
    ST shiftA_;             // X -> X + a
    ST shiftB_;             // Y -> Y + b

    // derived constants
    ST gamma_;              // gamma = -(L+1)D/2
    PT shiftedU_;           // U after replacement and shift

    // cache of F = U^gamma (2D series)
    mutable Series<ST> uPowerSeriesCache_;
    mutable bool uPowerSeriesCached_;

public:
    // Constructor initializes all constants and validates data.
    IntegrandExpander(SeriesSolver<RT, PT, ST>& solver,
                      int numLoops,
                      int targetDeg,
                      const ST& feynmanD,
                      const ST& shiftA,
                      const ST& shiftB);

    // Main API:
    // input:
    //   nu = outer FI exponents (used by both P(X,Y) and I_nu^Delta)
    // output:
    //   2D series of FI integrand at degree <= targetDeg_
    //
    // exact steps:
    //   1) poly = buildFIPolynomial(nu)
    //   2) fbi  = solver_.getFBISeries(nu, fbiDelta_)
    //   3) tmp  = poly * fbi
    //   4) uPow = getUPowerSeries()
    //   5) ret  = uPow * tmp
    //   6) return ret
    Series<ST> getFI2DSeries(const std::vector<int>& nu) const;

    // debug/read-only helpers
    const PT& getShiftedU() const;
    ST getGamma() const;
    ST getFeynmanD() const;
    ST getFBIDelta() const;
    const Series<ST>& getUPowerSeries() const;

    // clear only U-power series cache
    void clearCache();

private:
    // scalar parameter computation
    ST computeFBIDelta() const;   // Delta = L*D/2
    ST computeGamma() const;      // gamma = -(L+1)D/2

    // polynomial construction
    PT buildShiftedU() const;
    PT buildFIPolynomial(const std::vector<int>& nu) const;
    PT applyReplacementAndShift(const PT& xyzPoly) const;

    // build U^gamma via PDE recurrence
    Series<ST> expandUPower() const;
};
```

### 3.3 IntegrandExpander：复杂实现逻辑

#### A. 构造函数必须做的事

1. store all inputs
2. `fbiDelta_ = computeFBIDelta()`
3. `gamma_ = computeGamma()`
4. `shiftedU_ = buildShiftedU()`
5. initialize cache flag (`uPowerSeriesCached_ = false`)

#### A2. `getFI2DSeries(nu)` 的固定执行顺序

这个函数不要做额外推断，按固定顺序执行：

1. `poly = buildFIPolynomial(nu)`
2. `fbi = solver_.getFBISeries(nu, fbiDelta_)`
3. `tmp = poly * fbi`
4. `uPow = getUPowerSeries()`
5. `ret = uPow * tmp`
6. 返回 `ret`

其中：

- `poly` 是 FI 的多项式因子（不含 `U^gamma`）
- `fbi` 是 `I_nu^Delta` 的二维级数
- `uPow` 是 `U^gamma` 的二维级数缓存
- 函数内不重新计算 `Delta`、`gamma`、`shiftedU_`

#### B. `buildFIPolynomial(nu)` 具体逻辑

目标：构造 `P(X,Y)`，只包含多项式因子，不包含 `U^gamma`。

步骤：

1. 构造三变量表达式的 Jacobian 因子 `J_xyz = 1 - X`
2. 从 `nu` 读取 `(nuX, nuY, nuZ)`，构造
   `W_xyz = X^(nuX-1) * Y^(nuY-1) * Z^(nuZ-1)`
3. 先合并三变量多项式：`P_xyz = J_xyz * W_xyz`
4. 调用 `applyReplacementAndShift(P_xyz)`：
   - `X -> X`
   - `Y -> (1-X)Y`
   - `Z -> 1-X-(1-X)Y`
   - 再做 `X -> X+a, Y -> Y+b`
5. 返回二维多项式 `P(X,Y)`

#### C. `expandUPower()` 的 PDE 递推（关键）

目标：构造 `F = U^gamma` 的二维级数，满足 `F(0,0)=1`。

方程：

- `U * dF/dX = gamma * (dU/dX) * F`
- `U * dF/dY = gamma * (dU/dY) * F`

递推方式：

1. 初始化 `F(0,0) = 1`
2. 总度 `deg = 1..targetDeg_`
3. 对 `p = 1..deg` (`q = deg-p`)：
   - 用 X 方程解 `F(p,q)`
4. 对 `p = 0` (`q = deg`)：
   - 用 Y 方程解 `F(0,q)`

实现提示：

- 系数求解形式与现有 `SeriesSolver` 的递推模式一致
- 可复用“多项式与级数卷积取系数”的工具函数
- 只在第一次请求 `getUPowerSeries()` 时计算，之后走缓存

### 3.4 SeriesIntegrator（类声明 + 注释）

```cpp
template<typename ST>
struct IntegrationConfig {
    ST shiftA_;
    ST shiftB_;
    int degree_;   // max degree of 1D output
};

template<typename ST>
class SeriesIntegrator {
private:
    IntegrationConfig<ST> config_;

public:
    explicit SeriesIntegrator(const IntegrationConfig<ST>& config);

    // input: 2D series
    // output: 1D coefficient array coeff[d], d=0..degree_
    std::vector<ST> integrate(const Series<ST>& series) const;

private:
    // monomial integral weight for X^p Y^q under current config
    ST monomialWeight(int xPower, int yPower) const;
};
```

SeriesIntegrator 的实现公式和权重计算，直接参考 `reconstruct/src/integratedseries.tpp`：

- `performIntegration(...)`：二维级数到一维系数的累加框架
- 权重计算代码段：`weight_n`、`weight_m`、`total_weight`
- `integrate(...)` 的行为与上述逻辑一一对应

### 3.5 实现顺序

1. 写 `integrand_expander.hpp/.tpp`
2. 构造函数 + `computeFBIDelta` + `computeGamma`
3. `applyReplacementAndShift`
4. `buildShiftedU`
5. `buildFIPolynomial`
6. `expandUPower`（PDE 递推）
7. `getUPowerSeries` + `clearCache`
8. `getFI2DSeries`
9. 写 `series_integrator.hpp/.tpp`
10. `integrate` + `monomialWeight`
