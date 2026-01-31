# SeriesSolver 设计文档

## 1. 概述

SeriesSolver 是级数求解器，统一处理 **微分方程求解** 和 **IBP约化**，直接计算FBI的二维幂级数展开。

### 1.1 核心思想

将所有问题统一为 **Linear Recurrence Relation (LRR)** 形式：

$$D(X,Y) \cdot g(X,Y) = \sum_i N_i(X,Y) \cdot f_i(X,Y)$$

其中：
- $D(X,Y)$ 是公分母多项式
- $N_i(X,Y)$ 是分子多项式
- $g(X,Y)$ 是待求级数
- $f_i(X,Y)$ 是已知级数

### 1.2 与旧方法的区别

| 方面 | 旧方法 (FBIReducer) | 新方法 (SeriesSolver) |
|-----|-------------------|---------------------|
| 约化结果 | 有理函数系数 | 级数系数 |
| 缓存内容 | (ν,Δ) → 有理函数向量 | (ν,Δ) → Series |
| 运算方式 | 先符号约化，后级数乘除 | 直接LRR递推级数 |
| 复杂度 | 数值插值重构有理函数开销大 | 直接数值递推更高效 |

## 2. 数学原理

### 2.1 问题统一形式

任意约化关系或微分方程都可以写成：

$$g = \sum_i \frac{N_i}{D_i} f_i$$

**标准化**：乘以最小公倍多项式 $D = \text{lcm}(D_1, D_2, \ldots)$，得到：

$$D \cdot g = \sum_i \tilde{N}_i \cdot f_i$$

其中 $\tilde{N}_i = N_i \cdot D / D_i$ 是多项式。

### 2.2 微分方程

对于主积分 $I_{\vec{\nu}}^{\Delta}$，微分方程为：

$$\frac{\partial I_{\vec{\nu}}^{\Delta}}{\partial X} = \sum_{i,j} \left(-\frac{1}{2}\right) \frac{\partial R_{ij}}{\partial X} \cdot \text{factor}_{ij} \cdot I_{\vec{\nu}+\vec{e}_i+\vec{e}_j}^{\Delta+1}$$

其中：
$$\text{factor}_{ij} = \begin{cases}
\nu_i \nu_j & \text{if } i \neq j \\
\nu_i(\nu_i + 1) & \text{if } i = j
\end{cases}$$

右边的 $I_{\vec{\nu}+\vec{e}_i+\vec{e}_j}^{\Delta+1}$ 需要通过约化得到。

### 2.3 IBP约化（四种Case）

根据 Sector 的 Case 类型，约化公式不同：

**Case 0** ($\dim(\text{Null}(S))=0$, $C\neq 0$)：

1. IBP约化：$\nu_{\max} I_{\vec{\nu}}^{\Delta} = \sum_{j} (S^{-1})_{\text{row}, j} \cdot (\text{rhs})_j$
2. 维度迁移：$I_{\vec{\nu}}^{\Delta} = \frac{C}{(2\Delta - \nu - B) z_0} I_{\vec{\nu}}^{\Delta-1} - \frac{1}{(2\Delta - \nu - B) z_0} \sum_{\alpha} z_\alpha I_{\vec{\nu}-\vec{e}_\alpha}^{\Delta-1}$

**Case 1** ($\dim(\text{Null}(S))=0$, $C=0$)：
$$(2\Delta - \nu - B) I_{\vec{\nu}}^{\Delta} = -\sum_{\alpha=1}^{N} z_\alpha I_{\vec{\nu}-\vec{e}_\alpha}^{\Delta-1}$$

**Case 2** ($\dim(\text{Null}(S))>0$, $C\neq 0$)：
$$C \cdot I_{\vec{\nu}}^{\Delta-1} = \sum_{\alpha=1}^{N} z_\alpha I_{\vec{\nu}-\vec{e}_\alpha}^{\Delta-1}$$

**Case 3** ($\dim(\text{Null}(S))>0$, $C=0$)：
$$I_{\vec{\nu}}^{\Delta} = -\sum_{\alpha \neq \beta} \frac{z_\alpha}{z_\beta} I_{\vec{\nu}+\vec{e}_\beta-\vec{e}_\alpha}^{\Delta}$$

### 2.4 LRR递推公式

给定 $D \cdot g = \sum_i N_i \cdot f_i$，对于总度数 $p+q=N$ 的系数 $g_{pq}$：

$$g_{pq} = \frac{1}{D_{00}} \left( \sum_i [N_i \cdot f_i]_{pq} - \sum_{\substack{a+b \leq d \\ (a,b) \neq (0,0)}} D_{ab} \cdot g_{p-a,q-b} \right)$$

其中：
- $[N_i \cdot f_i]_{pq}$ 是多项式与级数卷积的第$(p,q)$项
- $g_{p-a,q-b}$ 是已知的（$p-a+q-b < N$）
- 要求 $D_{00} \neq 0$

## 3. 类设计

### 3.1 SeriesSolver类

```cpp
template<typename T>
class SeriesSolver {
public:
    // 构造函数
    SeriesSolver(Family<T>& family, int targetDeg);
    
    // 设置主积分边界条件（零阶系数）
    void setMasterBoundary(int masterIdx, T value);
    
    // 获取FBI的级数展开
    const Series<T>& getFBISeries(const std::vector<int>& nu, T delta);
    
    // 求解到指定度数
    void solve();
    
    // 获取当前已求解的度数
    int getCurrentDeg() const;

private:
    Family<T>& family_;
    int targetDeg_;
    int currentDeg_;
    
    // dR/dX 和 dR/dY 的多项式矩阵（N×N）
    std::vector<std::vector<Polynomial<T>>> dRdX_;
    std::vector<std::vector<Polynomial<T>>> dRdY_;
    
    // 主积分边界条件
    std::vector<T> masterBoundary_;
    
    // 主积分级数缓存
    std::vector<Series<T>> masterSeries_;
    
    // 一般FBI级数缓存：(nu, delta) → Series
    std::map<std::pair<std::vector<int>, T>, Series<T>> cache_;
    
    // 初始化dR/dX和dR/dY
    void initDerivatives();
    
    // LRR求解器
    void solveLRR(Series<T>& g, 
                  const Polynomial<T>& D,
                  const std::vector<std::pair<Polynomial<T>, const Series<T>*>>& terms,
                  int deg);
    
    // 微分方程递推（主积分）
    void solveMasterDE(int masterIdx, int deg);
    
    // 约化递推
    Series<T>& reduceAndSolve(const std::vector<int>& nu, T delta);
    
    // 四种Case的约化
    void reduceCase0(Series<T>& result, const std::vector<int>& nu, T delta);
    void reduceCase1(Series<T>& result, const std::vector<int>& nu, T delta);
    void reduceCase2(Series<T>& result, const std::vector<int>& nu, T delta);
    void reduceCase3(Series<T>& result, const std::vector<int>& nu, T delta);
};
```

### 3.2 依赖类

**Polynomial<T>**：二变量多项式

```cpp
template<typename T>
class Polynomial {
public:
    Polynomial();
    
    // 获取/设置系数
    T getCoeff(int i, int j) const;
    void setCoeff(int i, int j, const T& coeff);
    
    // 多项式运算
    Polynomial<T> operator+(const Polynomial<T>& other) const;
    Polynomial<T> operator*(const Polynomial<T>& other) const;
    Polynomial<T> operator*(const T& scalar) const;
    
    // 求导
    Polynomial<T> derivativeX() const;
    Polynomial<T> derivativeY() const;
    
    // 获取度数
    int getDeg() const;

private:
    std::unordered_map<std::pair<int,int>, T, PairHash> coeffs_;
    int deg_;
};
```

**Series<T>**：二维幂级数（复用backup实现）

```cpp
template<typename T>
class Series {
public:
    Series(int deg);
    
    T getCoeff(int i, int j) const;
    void setCoeff(int i, int j, const T& coeff);
    
    int getDeg() const;
    
    // 级数加法
    Series<T>& operator+=(const Series<T>& other);
    
    // 静态方法：避免临时对象
    static void mulPoly(Series<T>& result, const Series<T>& s, const Polynomial<T>& p);
    static void addScaled(Series<T>& result, const Series<T>& s, const T& scale);

private:
    std::vector<T> coefficients_;
    int deg_;
    
    int getIndex(int i, int j) const;
};
```

## 4. 算法流程

### 4.1 主循环

```cpp
void SeriesSolver<T>::solve() {
    // 初始化主积分零阶系数
    for (int k = 0; k < numMaster; k++) {
        masterSeries_[k].setCoeff(0, 0, masterBoundary_[k]);
    }
    
    // 逐度数递推
    for (int deg = 1; deg <= targetDeg_; deg++) {
        // 对每个主积分，用微分方程递推
        for (int k = 0; k < numMaster; k++) {
            solveMasterDE(k, deg);
        }
        currentDeg_ = deg;
    }
}
```

### 4.2 主积分微分方程递推

```cpp
void SeriesSolver<T>::solveMasterDE(int masterIdx, int deg) {
    // 微分方程: Q * df/dX = P * f + Q * g
    // 其中 g = sum_{i,j} (-1/2) * dR_{ij}/dX * factor_{ij} * I_{nu+e_i+e_j}^{delta+1}
    
    // 1. 计算右边的 g（递归调用getFBISeries获取I_{nu+e_i+e_j}）
    
    // 2. 转换为LRR形式: Q * df/dX = P * f + Q * g
    //    即: Q * (m * f_{m,n}) = [P * f + Q * g]_{m-1,n}
    
    // 3. 对于每个 (p,q) 满足 p+q=deg:
    //    如果 p > 0，使用X方向方程
    //    如果 p = 0, q > 0，使用Y方向方程
    
    // 4. 递推计算 f_{p,q}
}
```

### 4.3 约化与级数求解

```cpp
Series<T>& SeriesSolver<T>::reduceAndSolve(const std::vector<int>& nu, T delta) {
    auto key = std::make_pair(nu, delta);
    
    // 检查缓存
    if (cache_.count(key)) {
        return cache_[key];
    }
    
    // 创建新级数
    cache_[key] = Series<T>(targetDeg_);
    Series<T>& result = cache_[key];
    
    // 根据Case类型约化
    int caseType = family_.getCase(nu);
    switch (caseType) {
        case 0: reduceCase0(result, nu, delta); break;
        case 1: reduceCase1(result, nu, delta); break;
        case 2: reduceCase2(result, nu, delta); break;
        case 3: reduceCase3(result, nu, delta); break;
    }
    
    return result;
}
```

### 4.4 LRR求解器

```cpp
void SeriesSolver<T>::solveLRR(
    Series<T>& g, 
    const Polynomial<T>& D,
    const std::vector<std::pair<Polynomial<T>, const Series<T>*>>& terms,
    int deg
) {
    T D00 = D.getCoeff(0, 0);
    
    for (int p = 0; p <= deg; p++) {
        int q = deg - p;
        
        // 计算 sum_i [N_i * f_i]_{p,q}
        T rhs = T(0);
        for (const auto& [Ni, fi] : terms) {
            // 卷积 [Ni * fi]_{p,q}
            for (auto& [pow, coeff] : Ni) {
                int a = pow.first, b = pow.second;
                if (p >= a && q >= b) {
                    rhs += coeff * fi->getCoeff(p - a, q - b);
                }
            }
        }
        
        // 减去 sum_{(a,b)!=(0,0)} D_{a,b} * g_{p-a,q-b}
        for (auto& [pow, coeff] : D) {
            int a = pow.first, b = pow.second;
            if ((a != 0 || b != 0) && p >= a && q >= b) {
                rhs -= coeff * g.getCoeff(p - a, q - b);
            }
        }
        
        // g_{p,q} = rhs / D00
        g.setCoeff(p, q, rhs / D00);
    }
}
```

## 5. 复杂度分析

设：
- 目标度数 $N$
- 多项式度数 $d$
- 主积分数量 $M$
- 非主积分约化调用次数 $K$

**单次LRR求解**：
- 每个度数有 $N+1$ 个系数
- 每个系数需要 $O(d^2)$ 次运算
- 总计 $O(N \cdot d^2)$

**整体复杂度**：
- 主积分递推：$O(M \cdot N \cdot d^2)$
- 约化调用：$O(K \cdot N \cdot d^2)$
- 总计 $O((M + K) \cdot N \cdot d^2)$

对于 $N=1000$, $d=100$, $M=20$, $K=100$：约 $1.2 \times 10^9$ 次运算。

## 6. 实现计划

### 6.1 第一阶段：基础设施

- [ ] 实现 `Polynomial<T>` 类
- [ ] 复用/修改 `Series<T>` 类
- [ ] 测试多项式和级数运算

### 6.2 第二阶段：LRR求解器

- [ ] 实现 `solveLRR()` 方法
- [ ] 测试单个LRR求解

### 6.3 第三阶段：约化

- [ ] 实现 Case 0 约化
- [ ] 实现 Case 1 约化
- [ ] 实现 Case 2 约化
- [ ] 实现 Case 3 约化
- [ ] 测试约化正确性

### 6.4 第四阶段：微分方程

- [ ] 实现 `initDerivatives()` 计算 dR/dX, dR/dY
- [ ] 实现 `solveMasterDE()` 方法
- [ ] 测试主积分递推

### 6.5 第五阶段：集成

- [ ] 实现完整的 `solve()` 流程
- [ ] 实现 `getFBISeries()` 接口
- [ ] 端到端测试
- [ ] 性能优化

## 7. 参考

- `backup/include/fbi_reducer.hpp` - 旧版约化器
- `backup/include/de_builder.hpp` - 旧版微分方程构建器
- `backup/include/series.hpp` - 级数类
- `backup/docs/problem_and_workflow.md` - 数学原理详解
