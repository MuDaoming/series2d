---
id: PS003
title: MFBI 的微分方程建立与求解
status: active
created: 2026-03-31
updated: 2026-04-01
depends:
  - B001
  - PS002
parts:
  - PS001
---

# MFBI 的微分方程建立与求解

## 1. 目标

在 [PS002](./PS002_fbi_reduction.md) 完成 FBI 约化后，所有积分都可以表示为主积分（MFBI）的线性组合。本文档说明如何建立并求解 MFBI 关于分支参数 $(X, Y)$ 的微分方程，从而计算 MFBI 的二维幂级数展开。

## 2. FBI 对分支参数的导数

### 2.1 导数公式

由 FBI 定义（见 [B001](./../background/B001_feynman_integral_and_fbi.md)）：

$$
I_{\vec{\nu}}^{\Delta}(X,Y) = \frac{(-1)^{\nu}\Gamma(\nu - \Delta)}{\prod_\alpha \Gamma(\nu_\alpha)} \int [dy] \frac{\prod_\alpha y_\alpha^{\nu_\alpha - 1}}{(\frac{1}{2} y^T R y)^{\nu - \Delta}}
$$

对分支参数 $X$ 求导，利用

$$
\frac{\partial}{\partial X} \left( \frac{1}{(\frac{1}{2} y^T R y)^{\nu - \Delta}} \right) = -(\nu - \Delta) \cdot \frac{1}{(\frac{1}{2} y^T R y)^{\nu - \Delta + 1}} \cdot \frac{1}{2} \sum_{i,j} \frac{\partial R_{ij}}{\partial X} y_i y_j
$$

整理 $\Gamma$ 函数比值后，得到 FBI 对分支参数的导数公式：

$$
\frac{\partial I_{\vec{\nu}}^{\Delta}}{\partial X} = -\frac{1}{2} \sum_{i,j} \frac{\partial R_{ij}}{\partial X} \cdot \text{factor}_{ij}(\vec{\nu}) \cdot I_{\vec{\nu}+\vec{e}_i+\vec{e}_j}^{\Delta+1}
$$

其中

$$
\text{factor}_{ij}(\vec{\nu}) = \begin{cases}
\nu_i \nu_j & i \neq j \\
\nu_i(\nu_i + 1) & i = j
\end{cases}
$$

对 $Y$ 的导数公式完全类似，只需将 $\partial R_{ij}/\partial X$ 替换为 $\partial R_{ij}/\partial Y$。

### 2.2 微分方程系统

对于 MFBI 向量 $\mathbf{f} = (f_1, \dots, f_M)^T$，利用上式和 FBI 约化（见 [PS002](./PS002_fbi_reduction.md)），得到一阶线性偏微分方程组：

$$
\frac{\partial \mathbf{f}}{\partial X} = A_X(X,Y) \cdot \mathbf{f}, \qquad \frac{\partial \mathbf{f}}{\partial Y} = A_Y(X,Y) \cdot \mathbf{f}
$$

其中 $A_X, A_Y$ 是 $M \times M$ 的有理函数矩阵，其元素是分支参数 $X, Y$ 的有理函数。

## 3. 常点求解二维微分方程

### 3.1 标准微分方程

对于第 $i$ 个MFBI $f_i$，在三角化基下，方程简化为：

$$
\begin{cases}
\displaystyle\frac{\partial f_i}{\partial X} = R_1(X, Y) \cdot f_i + g_1(X, Y) \\[10pt]
\displaystyle\frac{\partial f_i}{\partial Y} = R_2(X, Y) \cdot f_i + g_2(X, Y)
\end{cases}
$$

其中：
- $R_1, R_2$ 是有理函数系数
- $g_1, g_2$ 是耦合项，依赖于已求解的 $f_1, \ldots, f_{i-1}$（因此可以依次求解）

### 3.2 级数递推公式

给定初值条件 $f_i(0, 0) = f_i^{(0)}$，设：

$$
f_i(X, Y) = \sum_{m+n \leq d} c_{mn} X^m Y^n
$$

以及 $R_1 = P_1/Q_1$，$R_2 = P_2/Q_2$。

将微分方程改写为：
$$
Q_1 \frac{\partial f}{\partial X} = P_1 \cdot f + Q_1 \cdot g_1, \quad Q_2 \frac{\partial f}{\partial Y} = P_2 \cdot f + Q_2 \cdot g_2
$$

**递推计算**（按总度数 $m+n$ 递增）：

- **当 $m > 0$ 时**，使用 $X$ 方向方程：

$$
c_{mn} = \frac{1}{m \cdot Q_1^{(0,0)}} \left[ [Q_1 g_1]_{m-1,n} + \sum_{i,j} P_1^{(m-1-i, n-j)} c_{ij} - \sum_{i,j} (i+1) Q_1^{(m-1-i, n-j)} c_{i+1,j} \right]
$$

- **当 $m = 0, n > 0$ 时**，使用 $Y$ 方向方程：

$$
c_{0n} = \frac{1}{n \cdot Q_2^{(0,0)}} \left[ [Q_2 g_2]_{0,n-1} + \sum_{i,j} P_2^{(0-i, n-1-j)} c_{ij} - \sum_{i,j} (j+1) Q_2^{(0-i, n-1-j)} c_{i,j+1} \right]
$$

通过这个递推关系，从初值 $c_{00} = f_i^{(0)}$ 开始，逐阶计算所有系数。

### 3.3 求解结果

得到所有MFBI的级数展开：

$$
f_k(X, Y) = \sum_{m+n \leq d} c_{mn}^{(k)} X^m Y^n, \quad k = 1, \ldots, M
$$

### 3.4 级数与有理函数的乘法

给定级数 $f(X, Y) = \sum c_{mn} X^m Y^n$ 和有理函数 $R(X,Y) = P(X,Y)/Q(X,Y)$，计算：

$$
h(X, Y) = f(X, Y) \cdot R(X, Y)
$$

**算法：**

1. 计算 $\tilde{f} = f \cdot P$（级数与多项式乘法，卷积）
2. 从 $Q \cdot h = \tilde{f}$ 解出 $h$ 的系数（按总度数递增）：

$$
h_{mn} = \frac{1}{Q^{(0,0)}} \left[ \tilde{f}_{mn} - \sum_{(i,j) \neq (0,0)} Q^{(i,j)} h_{m-i, n-j} \right]
$$

### 3.5 最终FBI级数

利用约化结果和MFBI级数：

$$
I_{\vec{\nu}}^{\Delta}(X, Y) = \sum_{k=1}^{M} r_k^{(\vec{\nu}, \Delta)}(X, Y) \cdot f_k(X, Y)
$$

对每项执行级数与有理函数乘法，然后求和，得到：

$$
I_{\vec{\nu}}^{\Delta}(X, Y) = \sum_{m+n \leq d} I_{mn}^{(\vec{\nu}, \Delta)} X^m Y^n
$$

## 4. 约化-微分方程结合的常点逐阶递推

### 4.1 算法逻辑

假设当前已知：
- 主积分的 $\deg \leq N$ 的项
- 其它所需积分的 $\deg < N$ 的项

则递推步骤为：

1. **约化步**：通过约化公式（第3节）计算其它所需积分的 $\deg = N$ 的项
2. **微分步**：利用微分方程（第4节）计算主积分导数的 $\deg = N$ 的项，即主积分 $\deg = N+1$ 的项

递推后得到：
- 主积分的 $\deg \leq N+1$ 的项
- 其它所需积分的 $\deg < N+1$ 的项

因此只需给定主积分的边界条件（$\deg = 0$ 的项）即可启动递推。

### 4.2 核心思想

将IBP约化关系和微分方程统一为**线性递推关系**（LRR），逐阶计算级数系数。

**关键特点**：
- 所有中间量都是**二维幂级数**而非有理函数
- IBP和微分方程的系数是**多项式**（关于 $X, Y$）
- 每一阶的计算都是纯数值运算（在有限域 $\mathbb{Z}_p$ 中）

### 4.3 线性递推关系（LRR）的统一形式

所有约化公式和微分方程都可以写成LRR形式：

$$
D(X,Y) \cdot g(X,Y) = \sum_i N_i(X,Y) \cdot f_i(X,Y)
$$

其中：
- $D(X,Y)$ 是公分母多项式
- $N_i(X,Y)$ 是分子多项式
- $g(X,Y)$ 是待求级数
- $f_i(X,Y)$ 是已知级数

### 4.4 LRR递推公式

对于总度数 $p+q = d$ 的系数 $g_{pq}$：

$$
g_{pq} = \frac{1}{D_{00}} \left( \sum_i [N_i \cdot f_i]_{pq} - \sum_{\substack{(a,b) \neq (0,0) \\ a+b \leq \deg(D)}} D_{ab} \cdot g_{p-a,q-b} \right)
$$

其中：
- $[N_i \cdot f_i]_{pq}$ 是多项式与级数卷积的第 $(p,q)$ 项：$[N_i \cdot f_i]_{pq} = \sum_{a,b} (N_i)_{ab} \cdot (f_i)_{p-a,q-b}$
- $g_{p-a,q-b}$ 是已知的（$p-a+q-b < p+q$）
- 要求 $D_{00} \neq 0$（常数项非零）
