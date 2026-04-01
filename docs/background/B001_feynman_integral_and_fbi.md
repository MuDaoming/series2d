---
id: B001
title: Feynman Integral and FBI
status: active
created: 2026-03-31
updated: 2026-03-31
---

# 费曼积分与固定分支积分 (FBI)

---

## 第一部分：费曼积分

### 1.1 费曼积分家族

一个费曼积分家族由一组圈动量 $l_1, \dots, l_L$（$L$ 为圈数）和一组外动量 $p_1, \dots, p_E$ 定义。给定一组传播子（propagators）的集合：

$$D_1, D_2, \dots, D_N$$

其中每个传播子的形式为：

$$D_\alpha = q_\alpha^2 - m_\alpha^2$$

这里 $q_\alpha$ 是圈动量和外动量的线性组合，$m_\alpha$ 是对应的质量。

一个费曼积分家族中的通用积分定义为：

$$I_{\nu_1, \nu_2, \dots, \nu_N} = \int \prod_{i=1}^{L} \frac{d^D l_i}{i\pi^{D/2}} \frac{1}{D_1^{\nu_1} D_2^{\nu_2} \cdots D_N^{\nu_N}}$$

其中 $\nu_\alpha$ 是各传播子的幂次（指标），$D$ 是时空维数（在维数正规化中通常取 $D = d_0 - 2\epsilon$）。

**家族中的积分**由指标向量 $\vec{\nu} = (\nu_1, \nu_2, \dots, \nu_N)$ 唯一标识。

通常，$\nu_1,\cdots,\nu_n\geq 0$ 为**传播子分母**，$\nu_{n+1},\cdots,\nu_{N}\leq 0$ 为**分子**（或 irreducible scalar product，ISP）。

### 1.2 Sector（扇区）

#### 1.2.1 Sector 的定义

在费曼积分理论中，**Sector（扇区）**是由当前积分的传播子分母的集合定义的。具体地，给定一个指标向量 $\vec{\nu}$，Sector 由以下二进制向量表示：

$$
\vec{s}(\vec{\nu}) = (s_1, \ldots, s_N), \quad s_i = \begin{cases} 1 & \nu_i > 0 \\ 0 & \nu_i \leq 0 \end{cases}
$$

#### 1.2.2 Sub-sector 和 Top-sector

- **Sub-sector（子扇区）**：如果 $\vec{s}_{sub}$ 满足对所有 $i$ 都有 $(s_{sub})_i \leq s_i$，则 $\vec{s}_{sub}$ 是 $\vec{s}$ 的 sub-sector。一般来说，sub-sector 是指那些非平凡的子扇区，即子扇区和扇区不相等。

- **Top-sector（顶级扇区）**：如果 $\vec{s}$ 是 $\vec{s}_{top}$ 的 sub-sector，则 $\vec{s}_{top}$ 是 $\vec{s}$ 的 top-sector。一个家族的 top-sector 是指由所有传播子分母定义的 sector，即 $\vec{s} = (1,1,\ldots,1)$ 的 sector。

### 1.3 Dot、Rank 和 Degree

对于给定的指标向量 $\vec{\nu}$，定义：

- **Dot**：传播子幂次超过 1 的部分之和
$$\text{dot} = \sum_i \max(\nu_i - 1, 0) = \sum_{\nu_i \geq 1} (\nu_i - 1)$$

- **Rank**：传播子负幂次的绝对值之和
$$\text{rank} = \sum_i \max(-\nu_i, 0) = \sum_{\nu_i < 0} (-\nu_i)$$

- **Degree**：dot 与 rank 之和
$$\text{degree} = \text{dot} + \text{rank}$$

直观理解：dot 表示分母的总幂次，rank 表示分子或者 ISP 的总幂次。

### 1.4 Corner 积分

在一个给定的 Sector 中，如果所有传播子的指标均满足 $\nu_\alpha \in \{0, 1\}$，则称该积分为**Corner 积分（角积分）**。

Corner 积分是该 Sector 中最简形式的积分：dot = 0，rank = 0，degree = 0。

---

## 第二部分：FBI（固定分支积分）

FBI（Fixed-Branch Integral）是论文 (arXiv: 2412.21053) 提出的一种费曼积分新表示。该表示的核心思想是将多圈费曼积分分解为**分支参数 (branch parameters)** 上的积分，其被积函数是结构上类似于单圈积分的 FBI。

### 2.1 Branch（分支）的定义

**Branch（分支）**的划分规则是：将传播子中的外动量全部删除（即只保留圈动量部分），如果两个传播子删除外动量后具有相同的圈动量结构，则它们属于同一个 Branch。

- **单圈**（$L=1$）：只有 1 个 Branch（$B=1$），因为所有传播子都只含同一个圈动量 $l_1$
- **两圈**（$L=2$）：至多 3 个 Branch（$B \leq 3$），对应圈动量的三种独立组合 $l_1$, $l_2$, $l_1 \pm l_2$
- **三圈**（$L=3$）：至多 6 个 Branch（$B \leq 6$）

每个分支有一个**分支参数** $X_b$，满足积分测度 $[dX] = \prod_{b=1}^{B} dX_b \, \delta(1 - \sum_{b=1}^{B} X_b)$。

在每个分支内部，有费曼参数 $y_{(b,i)}$，满足 $[dy] = \prod_{\alpha=1}^{N} dy_\alpha \prod_{b=1}^{B} \delta(1 - \sum_{i=1}^{n_b} y_{(b,i)})$。

### 2.2 FBI 的定义

出发点是费曼参数化。对于一个具有 $B$ 个分支、第 $b$ 个分支有 $n_b$ 条传播子的费曼图，总传播子数为 $N = \sum_{b=1}^{B} n_b$。费曼积分可以写成：

$$M = \int [dX] \, \hat{M}(X)$$

其中被积函数 $\hat{M}(X)$ 可以展开为一组**固定分支积分 (FBI)** 的线性组合：

$$\hat{M}(X) = U^{-\frac{(L+1)D}{2}} \sum_{\Delta, \vec{\nu}'} K_{\vec{\nu}'}^{\Delta}(X) \, I_{\vec{\nu}'}^{\Delta}(X)$$

其中 FBI 定义为：

$$\boxed{I_{\vec{\nu}}^{\Delta}(X) = \frac{(-1)^{\nu}\Gamma(\nu - \Delta)}{\prod_{\alpha=1}^{N}\Gamma(\nu_\alpha)} \int [dy] \frac{\prod_{\alpha=1}^{N} y_\alpha^{\nu_\alpha - 1}}{(\frac{1}{2} y^T \cdot R \cdot y - i0^+)^{\nu - \Delta}}}$$

这里：
- $\Delta$ 是与时空维数相关的参数
- $\vec{\nu} = (\nu_1, \dots, \nu_N)$ 是指标向量
- $R$ 是一个对称矩阵，其元素是分支参数 $X$ 的齐次多项式（次数为 $L+1$）
- $U$ 和 $F$ 是 Symanzik 多项式，$F = \frac{1}{2} y^T \cdot R \cdot y$

**关键特性**：FBI 在结构上类似于**单圈费曼积分的费曼参数表示**，因此可以用类似单圈积分的方法处理。

### 2.3 S 矩阵的定义

对于一个 FBI 家族，论文定义了一个对称矩阵 $S$，其构造规则如下（推广自单圈记号 [Ref. 40]）：

设总共有 $B$ 个分支，第 $b$ 个分支有 $n_b$ 条传播子，$N = \sum_{b} n_b$。$S$ 是一个 $(B+N) \times (B+N)$ 的对称矩阵，索引从 $1$ 到 $B+N$，其中前 $B$ 个索引对应分支，后 $N$ 个索引对应传播子。具体定义为：

1. **$S_{\alpha,\beta} = R_{\alpha-B, \beta-B}$**，当 $\alpha > B$ 且 $\beta > B$ 时（即两个索引都对应传播子时，取 $R$ 矩阵的对应元素）

2. **$S_{\alpha,\beta} = 1$**，当 $\alpha \leq B$（对应分支）且传播子 $\beta - B$ 属于分支 $\alpha$ 时；或对称地，$\beta \leq B$ 且传播子 $\alpha - B$ 属于分支 $\beta$ 时

3. **$S_{\alpha,\beta} = 0$**，其他情况（包括分支-分支之间的矩阵元）

**例子**：若 $B = 3$，$(n_1, n_2, n_3) = (2, 1, 1)$，则 $N = 4$，$S$ 是 $7 \times 7$ 矩阵：

$$S = \begin{pmatrix} & & & 1 & 1 & 0 & 0 \\ & 0_{3\times 3} & & 0 & 0 & 1 & 0 \\ & & & 0 & 0 & 0 & 1 \\ 1 & 0 & 0 & & & & \\ 1 & 0 & 0 & & R & & \\ 0 & 1 & 0 & & & & \\ 0 & 0 & 1 & & & & \end{pmatrix}$$

其中左上角 $3 \times 3$ 块为零矩阵（分支-分支），右上和左下块编码分支-传播子的归属关系，右下 $4 \times 4$ 块为 $R$ 矩阵。

### 2.4 Master FBI（MFBI）

在 FBI 约化中，所有积分最终都会约化到各 sector 的 corner 积分。**Master FBI（MFBI）**是指在 top-sector 中属于特定 Case（Case 0）的 corner 积分，它们构成约化的基。所有其他 FBI 都可以表示为 MFBI 的线性组合。

---

## 总结

| 概念 | 含义 |
|------|------|
| **费曼积分家族** | 由一组传播子定义的积分集合，由指标 $\vec{\nu}$ 参数化 |
| **Sector** | 由非零传播子的集合定义，用二进制向量标识 |
| **Sub-sector** | 从当前 sector 删除某些传播子后形成的扇区 |
| **Top-sector** | 由所有传播子分母定义的扇区 |
| **Dot** | $\sum_{\nu_i \geq 1} (\nu_i - 1)$ |
| **Rank** | $\sum_{\nu_i < 0} (-\nu_i)$ |
| **Degree** | dot + rank |
| **Corner 积分** | 所有传播子指标 $\nu_\alpha \in \{0, 1\}$ 的积分 |
| **Branch（分支）** | 删除外动量后圈动量结构相同的传播子归为同一分支 |
| **FBI** | 固定分支参数 $X$ 后的积分，结构类似单圈费曼参数积分 |
| **$S$ 矩阵** | $(B+N)\times(B+N)$ 对称矩阵，编码分支归属和 $R$ 矩阵信息 |
| **MFBI** | Top-sector 中 Case 0 的 corner 积分，构成约化的基 |