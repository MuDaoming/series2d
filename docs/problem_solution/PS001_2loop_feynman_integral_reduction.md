---
id: PS001
title: 二圈费曼积分的约化：问题定义
status: active
created: 2026-03-31
updated: 2026-03-31
---

# 二圈费曼积分的约化

## 1. 问题背景

两圈费曼积分在完成 FBI 表示后，具有以下特殊性质：

1. **Branch 结构**：两圈费曼积分只有三个 Branch（$B = 3$），考虑 Branch 的 $\delta$ 函数约束，本质上是二重积分。

2. **被积函数结构**：被积函数只包含以下三类因子：
   - $X$ 和 $Y$ 的有理函数 $Q(X, Y)$
   - 类单圈的 FBI（固定分支积分）
   - 关于 $U$ 的非整数幂次

这两点决定了在选取合适的方式之后，积分的复杂度不会太高。

## 2. 问题表述

我们需要处理的问题是：**利用 FBI 新表示完成对两圈费曼积分的约化**。

经过 FBI 表示和变量重标度后（详见 [B001](./../background/B001_feynman_integral_and_fbi.md)），两圈费曼积分可以写成如下形式：

$$
\text{FI}_{\vec{\nu}}^{D} = \int_0^1 dX \int_0^1 dY \; U(X,Y)^{-\gamma} \cdot \sum_i Q_i(X,Y) \cdot I_{\vec{\nu}_i}^{\Delta_i}(X,Y)
$$

其中：
- $\vec{\nu}$ 是费曼积分的指标向量
- $D$ 是时空维度
- $U(X,Y)$ 是 Symanzik 多项式
- $\gamma$ 是与维度相关的非整数参数
- $Q_i(X,Y)$ 是 $X, Y$ 的有理函数
- $I_{\vec{\nu}_i}^{\Delta_i}(X,Y)$ 是 FBI

**核心问题**：如何找到上述费曼积分之间的约化关系？

## 3. 方法概述：级数展开搜寻关系

我们采用**级数展开**的方法来搜寻费曼积分之间的关系。具体步骤如下：

### 3.1 变积分上限

引入参数 $\delta \in [0,1]$，将积分区域从 $[0,1]^2$ 拓展为：

$$
\text{FI}_{\vec{\nu}}^{D}(a,b;\delta) = \int_{a-a\delta}^{a+(1-a)\delta} dX \int_{b-b\delta}^{b+(1-b)\delta} dY \; U(X,Y)^{-\gamma} \cdot \sum_i Q_i(X,Y) \cdot I_{\vec{\nu}_i}^{\Delta_i}(X,Y)
$$

其中 $(a,b) \in (0,1)^2$ 是展开点。

- 当 $\delta = 1$ 时，积分区域回到 $[0,1]^2$，即原始费曼积分 $\text{FI}_{\vec{\nu}}^{D}(a,b;1) = \text{FI}_{\vec{\nu}}^{D}$
- 当 $\delta = 0$ 时，积分区域收缩为点 $(a,b)$，此时 $\text{FI}_{\vec{\nu}}^{D}(a,b;0)$ 只包含被积函数在点 $(a,b)$ 的值

**本质**：这是在点 $(a,b)$ 处的级数展开。

### 3.2 计算被积函数的二维展开

计算 $\text{FI}_{\vec{\nu}}^{D}(a,b;\delta)$ 的展开系数 $c_n(a,b)$ 的关键是**计算被积函数的二维展开**，具体步骤为：

**第一步：将一般指标的 FBI 约化到 MFBI**（详见 [PS002](./PS002_fbi_reduction.md)）

首先利用 IBP 恒等式和维度迁移公式，将被积函数中的任意 FBI 表示为 MFBI 的线性组合。约化后（暂不考虑 $U^{-\gamma}$ 部分），被积函数变为：

$$
\sum_k R_k(X,Y) \cdot I_{\vec{\nu}_k}^{\Delta_k}(X,Y)
$$

其中 $R_k(X,Y)$ 是 $X, Y$ 的有理函数，$I_{\vec{\nu}_k}^{\Delta_k}(X,Y)$ 是 MFBI。

**第二步：建立并求解 MFBI 的微分方程**（详见 [PS003](./PS003_mfbi_differential_equation.md)）

MFBI 作为 $X, Y$ 的函数，满足一阶线性偏微分方程组：

$$
\frac{\partial \mathbf{f}}{\partial X} = A_X(X,Y) \cdot \mathbf{f}, \qquad \frac{\partial \mathbf{f}}{\partial Y} = A_Y(X,Y) \cdot \mathbf{f}
$$

给定边界条件后，可通过递推求解 MFBI 的二维级数展开。

**第三步：乘以有理式，再处理 $U^{-\gamma}$**

乘以有理函数后仍为二维级数，按卷积公式计算系数：

$$
I_{\vec{\nu}_k}^{\Delta_k}(X,Y) = \sum_{p,q} f_{pq}^{(k)} X^p Y^q
$$

对于 $U(X,Y)^{-\gamma}$ 的非整数幂次，可将其吸收到 FBI 的定义中，改变微分方程的形式，使得展开计算仍然可行。

### 3.3 寻找积分关系

将 $\text{FI}_{\vec{\nu}}^{D}(a,b;\delta)$ 在 $\delta = 0$ 处展开为幂级数：

$$
\text{FI}_{\vec{\nu}}^{D}(a,b;\delta) = \sum_{n=0}^{N} c_n(a,b) \, \delta^n + \mathcal{O}(\delta^{N+1})
$$

假设不同 FBI 之间存在如下线性关系：

$$
\sum_k \left( c_{\vec{\nu}_k,0} + c_{\vec{\nu}_k,1} \delta + \cdots + c_{\vec{\nu}_k,d} \delta^d \right) \cdot I_{\vec{\nu}_k}^{\Delta_k}(a,b;\delta) = 0
$$

其中 $c_{\vec{\nu},n}$ 是待定的系数。通过比较 $\delta$ 各幂次的系数，可以建立费曼积分之间的约化关系。

## 4. 适用范围

**目前**：本项目处理不含 ISP（rank = 0）且 Branch 数目为 3 的两圈费曼积分。

**原因**：目前尚不清楚如何在统一框架下处理不同 Branch 数目的积分，以及如何处理含 ISP（rank > 0）的情形。这些问题留待后续研究。
