---
id: PS004
title: 张量分解方法：将带 rank 积分表示为 FBI
status: active
created: 2026-04-01
updated: 2026-04-01
depends:
  - B001
---

# 费曼积分分子约化：张量分解方法

---

## 1. 动机

在计算费曼积分时，我们经常遇到分子中含有圈动量的情况：

$$\int \prod_{i=1}^{L} \frac{d^D l_i}{i\pi^{D/2}} \frac{\text{Num}(l_i \cdot p_j,\; l_i \cdot l_j,\; \dots)}{D_1^{\nu_1} \cdots D_N^{\nu_N}} \tag{1}$$

其中分子 $\text{Num}$ 是圈动量 $l_i$ 与外动量 $p_j$ 的标量积的多项式。目标是**将分子中的圈动量消除**，最终把带分子的积分约化为不含分子的标量积分（即 FBI）的线性组合。

本文档记录一种基于**辅助动量 $\vec{P}$ 的张量分解方法**，其核心思路是：

> 引入辅助动量 $\vec{P}$，将分子中的圈动量替换为对 $\vec{P}$ 的微分算子，从而先完成无分子的圈动量积分，再通过对 $\vec{P}$ 求导恢复分子的贡献。

---

## 2. 从费曼积分到 FBI 表示中的分子问题

### 2.1 费曼参数化与合并传播子

由论文 arXiv:2412.21053 的 Eq. (2)，引入分支参数 $X_b$ 和费曼参数 $y_{(b,i)}$，将传播子合并：

$$\boxed{\frac{1}{D_1^{\nu_1} \cdots D_N^{\nu_N}} = \frac{\Gamma(\nu)}{\prod_{\alpha=1}^{N}\Gamma(\nu_\alpha)} \int_0^\infty [dX][dy] \frac{\prod_{b} X_b^{\nu_b-1} \prod_{\alpha} y_\alpha^{\nu_\alpha-1}}{\left(\sum_{b=1}^{B}\sum_{i=1}^{n_b} X_b y_{(b,i)} D_{(b,i)}\right)^{\nu}}} \tag{2}$$

其中 $\nu = \sum_\alpha \nu_\alpha$，$[dX]$ 和 $[dy]$ 分别含有 $\delta$ 函数约束。

将此代入 (1)，得到：

$$\boxed{M = \frac{\Gamma(\nu)}{\prod \Gamma(\nu_\alpha)} \int [d^D l] \cdot \text{Num}(\vec{l}) \int [dX][dy] \frac{\prod X_b^{\nu_b-1} \prod y_\alpha^{\nu_\alpha-1}}{\left(\sum_b \sum_i X_b y_{(b,i)} D_{(b,i)}\right)^{\nu}}} \tag{3}$$

### 2.2 展开为关于圈动量的二次型

由论文 Eq. (5)，合并后的传播子展开为圈动量的二次型：

$$\sum_{b=1}^{B}\sum_{i=1}^{n_b} X_b y_{(b,i)} D_{(b,i)} = \sum_{i,j=1}^{L} A_{ij}\, l_i \cdot l_j + 2\sum_{i=1}^{L} B_i \cdot l_i + C \tag{4}$$

其中：
- $A_{ij}$ 是 $L \times L$ 对称矩阵，元素依赖于分支参数 $X$（不依赖 $y$）
- $B_i^\mu$ 是 $D$ 维矢量，依赖于 $X$、$y$ 和外动量
- $C$ 是标量，依赖于 $X$、$y$、外动量和质量

于是 (3) 变为：

$$
\boxed{
\begin{aligned}
M = & \frac{\Gamma(\nu)}{\prod \Gamma(\nu_\alpha)} \int [d^D l] \cdot \text{Num}(\vec{l}) \times \\
& \int [dX][dy] \frac{\prod X_b^{\nu_b-1} \prod y_\alpha^{\nu_\alpha-1}}{\left(\sum_{i,j} A_{ij}\, l_i \cdot l_j + 2\sum_i B_i \cdot l_i + C\right)^{\nu}}
\end{aligned}
}\tag{5} $$

---

## 3. 圈动量平移：消除线性项

### 3.1 完成配方

对分母中关于 $\vec{l}$ 的二次型做配方。令

$$\vec{l} \;\to\; \vec{l} - A^{-1}\vec{B}$$

即 $l_i^\mu \to l_i^\mu - (A^{-1})_{ij} B_j^\mu$。平移后：

$$\sum_{i,j} A_{ij}\, l_i \cdot l_j + 2\sum_i B_i \cdot l_i + C \;\longrightarrow\; \vec{l}^{\,T} A\, \vec{l} + C - \vec{B}^{\,T} A^{-1} \vec{B}$$

其中 $\vec{l}^{\,T} A\, \vec{l} \equiv \sum_{i,j} A_{ij}\, l_i \cdot l_j$，而标量部分恰好是：

$$\boxed{C - \vec{B}^{\,T} A^{-1} \vec{B} = -\frac{F}{U}} \tag{6}$$

这里 $U = \det(A)$，$F = \frac{1}{2} y^T R\, y$ 是 Symanzik 多项式。

### 3.2 分子的变换

平移后，分子也随之变换：

$$\text{Num}(\vec{l}) \;\longrightarrow\; \text{Num}(\vec{l} - A^{-1}\vec{B})$$

于是积分变为：

$$\boxed{M = \frac{\Gamma(\nu)}{\prod \Gamma(\nu_\alpha)} \int [dX][dy] \prod X_b^{\nu_b-1} \prod y_\alpha^{\nu_\alpha-1} \int [d^D l]\; \frac{\text{Num}(\vec{l} - A^{-1}\vec{B})}{\left(\vec{l}^{\,T} A\, \vec{l} - F/U\right)^{\nu}}} \tag{7}$$

---

## 4. 分子展开：多项式分解

### 4.1 展开为圈动量的单项式

将 $\text{Num}(\vec{l} - A^{-1}\vec{B})$ 按 $\vec{l}$ 展开为多项式。由于 $\text{Num}$ 原本是 $l_i \cdot p_j$、$l_i \cdot l_j$ 等标量积的多项式，展开后的每一项是圈动量洛伦兹分量的乘积：

$$\text{Num}(\vec{l} - A^{-1}\vec{B}) = \sum_{\text{tensor structure}} \mathcal{C}^{(r)}_{\mu_1 \cdots \mu_r}(A^{-1}\vec{B},\, p,\, m) \;\cdot\; l_{a_1}^{\mu_1} \cdots l_{a_r}^{\mu_r} \tag{8}$$

其中：

- 系数 $\mathcal{C}$ 依赖于 $A^{-1}\vec{B}$（含 $X$、$y$ 和外动量）以及外部运动学参数
- 由于 $\text{Num}$ 是标量，所有洛伦兹指标最终都是缩并的


### 4.2 积分的分离

将展开代入 (7)，得到：

$$\boxed{M = \sum_{\text{tensor structure}} \int [dX][dy]\; \Pi_{X,y}\cdot \mathcal{C}^{(r)}_{\mu_1 \cdots \mu_r}(X,y,p_i^{\mu},m_j) \cdot \frac{\Gamma(\nu)}{\prod \Gamma(\nu_\alpha)} \underbrace{\int [d^D l]\; \frac{l_{a_1}^{\mu_1} \cdots l_{a_r}^{\mu_r}}{\left(\vec{l}^{\,T} A\, \vec{l} - F/U\right)^{\nu}}}_{\text{张量型圈动量积分}}} \tag{9}$$

其中 $\Pi_{X,y}$ 表示省略的 $\Gamma$ 函数和参数幂次因子。

**核心问题**：如何计算分母仅含 $\vec{l}^{\,T} A\, \vec{l} - F/U$ 的**张量型**圈动量积分？

---

## 5. 辅助动量方法：张量分解的核心

### 5.1 引入辅助动量 $\vec{P}$

对每个圈动量 $l_i$（$i = 1, \dots, L$），引入一个辅助 $D$ 维矢量 $P_i^\mu$。定义**生成函数**：

$$\boxed{G^D_\nu(\vec{P}) \equiv \int \prod_{i=1}^{L} \frac{d^D l_i}{i\pi^{D/2}} \frac{1}{\left(\vec{l}^{\,T} A\, \vec{l} - F/U + \vec{P} \cdot \vec{l}\right)^{\nu}}} \tag{10}$$

其中 $\vec{P} \cdot \vec{l} \equiv \sum_{i=1}^{L} P_i^\mu\, l_{i,\mu}$。

### 5.2 分子 → 微分算子

关键观察：对 $G^D_\nu(\vec{P})$ 关于 $P_i^\mu$ 求导，每次导数会从分母中"拉下"一个 $l_i^\mu$ 因子（同时使分母幂次增加 1）：

$$\frac{\partial}{\partial P_a^{\mu}} \frac{1}{(\cdots + \vec{P} \cdot \vec{l})^{\nu}} = -\nu \cdot \frac{l_a^{\mu}}{(\cdots + \vec{P} \cdot \vec{l})^{\nu+1}}$$

一般地，$r$ 阶导数给出：

$$\frac{\partial^r G^D_\nu(\vec{P})}{\partial P_{a_1}^{\mu_1} \cdots \partial P_{a_r}^{\mu_r}} \bigg|_{\vec{P}=0} = \frac{(-1)^r\, \Gamma(\nu+r)}{\Gamma(\nu)} \int [d^D l]\; \frac{l_{a_1}^{\mu_1} \cdots l_{a_r}^{\mu_r}}{\left(\vec{l}^{\,T} A\, \vec{l} - F/U\right)^{\nu+r}} \tag{11}$$

由 (11) 令 $\nu \to \nu - r$，可得张量积分与生成函数的关系：

$$\boxed{\int [d^D l]\; \frac{l_{a_1}^{\mu_1} \cdots l_{a_r}^{\mu_r}}{\left(\vec{l}^{\,T} A\, \vec{l} - F/U\right)^{\nu}} = \frac{(-1)^r\, \Gamma(\nu-r)}{\Gamma(\nu)} \frac{\partial^r G^D_{\nu-r}(\vec{P})}{\partial P_{a_1}^{\mu_1} \cdots \partial P_{a_r}^{\mu_r}} \bigg|_{\vec{P}=0}} \tag{12}$$

> **注意**：要保持分母幂次为原始的 $\nu$，生成函数的幂次参数需要降为 $\nu - r$（而非 $\nu$），因为 $r$ 阶导数会使分母幂次增加 $r$。

### 5.3 先积圈动量，后求导

方法的核心流程是：

1. **先计算无分子的标量积分** $G^D_{\nu-r}(\vec{P})$
2. **再对结果求 $\vec{P}$ 的导数**，恢复分子的贡献
3. **最后令 $\vec{P} = 0$**

这样就把"带分子的张量积分"转化为"无分子的标量积分 + 微分运算"。

---

## 6. 标量生成函数的计算

### 6.1 对 $G^D_\nu(\vec{P})$ 再次配方

$G^D_\nu(\vec{P})$ 的被积函数分母中，$\vec{l}$ 出现在二次型和线性项中：

$$\vec{l}^{\,T} A\, \vec{l} + \vec{P} \cdot \vec{l} - F/U$$

对 $\vec{l}$ 再做一次配方，令 $\vec{l} \to \vec{l} - \frac{1}{2}A^{-1}\vec{P}$：

$$\vec{l}^{\,T} A\, \vec{l} + \vec{P} \cdot \vec{l} - F/U \;\longrightarrow\; \vec{l}^{\,T} A\, \vec{l} - \frac{1}{U}\left(F + \frac{1}{4}\vec{P}^{\,T} A^{\text{adj}} \vec{P}\right) \tag{13}$$

其中利用了 $A^{-1} = A^{\text{adj}}/U$（$U = \det A$）。

### 6.2 标准多圈高斯积分

平移后的积分是标准的多圈高斯型，利用附录 A 的公式可直接求值，得到：

$$\boxed{G^D_\nu(\vec{P}) = \frac{(-1)^\nu\, \Gamma(\nu - LD/2)}{\Gamma(\nu)} \cdot \frac{U^{\nu - (L+1)D/2}}{\left(F + \frac{1}{4}\vec{P}^{\,T} A^{\text{adj}} \vec{P}\right)^{\nu - LD/2}}} \tag{14}$$

> 注意：$G^D_\nu(\vec{P})$ 现在是 $\vec{P}$ 的**显式已知函数**，仅涉及 Symanzik 多项式 $F$（在分母）、$U$（在分子）和伴随矩阵 $A^{\text{adj}}$，可以直接对其求导。

### 6.3 代回 $M$ 的表达式

将 (12) 和 (14) 代回 (9)。对于分子展开中张量阶数为 $r$ 的项，需要使用 $G^D_{\nu-r}(\vec{P})$。由 (14) 令 $\nu \to \nu - r$：

$$G^D_{\nu-r}(\vec{P}) = \frac{(-1)^{\nu-r}\, \Gamma(\nu - r - LD/2)}{\Gamma(\nu-r)} \cdot \frac{U^{\nu - r - (L+1)D/2}}{\left(F + \frac{1}{4}\vec{P}^{\,T} A^{\text{adj}} \vec{P}\right)^{\nu - r - LD/2}} \tag{15}$$

代入后得到：

$$\boxed{\begin{aligned}
M =  & \sum_{\text{tensor structure}} \int [dX] \prod X_b^{\nu_b-1}\, U^{\nu-r-(L+1)D/2}\frac{(-1)^{\nu}\, \Gamma(\nu - r - LD/2)}{\prod \Gamma(\nu_\alpha)} \int [dy]\, \\
& \prod y_\alpha^{\nu_\alpha-1}\, \mathcal{C}^{(r)}_{\mu_1\cdots\mu_r} \cdot \frac{\partial^r}{\partial P_{a_1}^{\mu_1} \cdots \partial P_{a_r}^{\mu_r}} \frac{1}{\left(F + \tfrac{1}{4}\vec{P}^{\,T} A^{\text{adj}} \vec{P}\right)^{\nu-r-LD/2}} \bigg|_{\vec{P}=0}
\end{aligned}} \tag{16}$$

对 $\vec{P}$ 的求导由 Wick 缩并公式给出（见第 7 节）。

---

## 7. 求导与最终结果

### 7.1 对 $\vec{P}$ 求导的结构

由 (14)，$G^D_\nu(\vec{P})$ 关于 $\vec{P}$ 的依赖仅通过 $\vec{P}^{\,T} A^{\text{adj}} \vec{P} = \sum_{ij} (A^{\text{adj}})_{ij}\, P_i \cdot P_j$ 出现。因此：

- **奇数阶导数**在 $\vec{P} = 0$ 处为零（与洛伦兹对称性一致：奇数个圈动量的积分为零）
- **偶数阶导数**在 $\vec{P} = 0$ 处非零，结果由 $A^{\text{adj}}$ 的矩阵元和度规张量 $g^{\mu\nu}$ 的组合给出

### 7.2 一般 Wick 缩并公式（$r$ 为偶数，$r = 0, 2, 4, \ldots$）

对 $(F + \frac{1}{4}\vec{P}^{\,T} A^{\text{adj}} \vec{P})^{-\alpha}$ 求 $r$ 阶导数并令 $\vec{P}=0$，结果是所有可能的 **Wick 缩并**的求和（类似于自由场论中 Wick 定理）：

$$\boxed{
\begin{aligned}
\frac{\partial^{r}}{\partial P_{a_1}^{\mu_1} \cdots \partial P_{a_{r}}^{\mu_{r}}} \frac{1}{(F + \frac{1}{4}\vec{P}^{\,T} A^{\text{adj}} \vec{P})^{\alpha}} \bigg|_{\vec{P}=0} = & \frac{(-1)^{r/2}\Gamma(\alpha+r/2)}{\Gamma(\alpha)} \cdot \frac{1}{F^{\alpha+r/2}} \times \\
& \frac{1}{2^{r/2}}\sum_{\text{pairings}} \prod_{\text{pairs } (i,j)} (A^{\text{adj}})_{a_i a_j}\, g^{\mu_i \mu_j}
\end{aligned}
} \tag{17}$$

每次配对（Wick 缩并）消耗两个指标，并使 $F$ 在分母中的幂次增加 1（$\Delta$ 增加 1）。共有 $(r-1)!! = r!/(2^{r/2}(r/2)!)$ 种配对方式。

### 7.3 代入 Wick 缩并后的 $M$

将 (17) 代入 (16)（其中 $\alpha = \nu - r - LD/2$，$r$ 为偶数），得到最终结果：

$$\boxed{\begin{aligned}
M =  & \sum_{\text{tensor structure}} \int [dX] \prod X_b^{\nu_b-1}\, U^{\nu-r-(L+1)D/2} \cdot \frac{(-1)^{\nu-r/2}\, \Gamma(\nu - r/2 - LD/2)}{\prod \Gamma(\nu_\alpha)} \times \\
& \int [dy]\, \prod y_\alpha^{\nu_\alpha-1}\, \frac{1}{F^{\nu-r/2-LD/2}} \cdot \frac{1}{2^{r/2}} \sum_{\text{pairings}} \prod_{\text{pairs } (i,j)} (A^{\text{adj}})_{a_i a_j}\, g^{\mu_i \mu_j} \cdot \mathcal{C}^{(r)}_{\mu_1\cdots\mu_r}
\end{aligned}} \tag{18}$$
注意, $\mathcal{C}^{(r)}_{\mu_1\cdots\mu_r}$ 里面还有 $X$ 和 $y$；$r$ 随张量结构变化。

---

## 8. 完整流程总结

$$\boxed{\text{带分子的费曼积分} \;\xrightarrow{\text{张量分解}}\; \text{不同 } \Delta \text{ 的 FBI 的线性组合}}$$

具体步骤：

| 步骤 | 操作 | 公式 |
|------|------|------|
| 1 | 费曼参数化，合并传播子 | **(3)** $M = \cdots \int [dX][dy] \frac{\cdots}{(\sum X_b y_{b,i} D_{b,i})^\nu}$ |
| 2 | 展开为关于 $\vec{l}$ 的二次型 | **(5)** 分母 $= (l^T A l + 2B \cdot l + C)^\nu$ |
| 3 | 平移 $\vec{l} \to \vec{l} - A^{-1}\vec{B}$ | **(7)** 分母 $= (l^T A l - F/U)^\nu$；分子变为 $\text{Num}(l - A^{-1}B)$ |
| 4 | 展开分子为 $\vec{l}$ 的多项式 | **(9)** $M = \sum_{\text{tensor structure}} \mathcal{C}^{(r)} \cdot \int [d^D l] \frac{l_{a_1}^{\mu_1} \cdots l_{a_r}^{\mu_r}}{(\cdots)^\nu}$ |
| 5 | 张量积分 → $G^D_{\nu-r}(\vec{P})$ 的导数 | **(12)** 张量积分 $= \frac{(-1)^r \Gamma(\nu-r)}{\Gamma(\nu)} \partial^r_P G^D_{\nu-r}(\vec{P})\big\|_{P=0}$ |
| 6 | 计算 $G^D_{\nu-r}(\vec{P})$ | **(14)** $G^D_\nu(\vec{P}) \propto U^{\nu-(L+1)D/2} / (F + \frac{1}{4} P^T A^{\text{adj}} P)^{\nu-LD/2}$ |
| 7 | 对 $\vec{P}$ 求导（Wick 缩并）并令 $\vec{P} = 0$ | **(17)** → $(A^{\text{adj}})_{ab}\, g^{\mu\nu}$，$F$ 幂次增加 |
| 8 | 组合所有项 | **(18)** $M = \sum_{\text{tensor structure}} \cdots U^{\nu-r-(L+1)D/2} \cdot \frac{1}{F^{\nu-r/2-LD/2}} \cdots$ |

---

## 附录 A：标准多圈高斯积分公式

$L$-圈标量积分（$A$ 为 $L \times L$ 正定对称矩阵，$\Delta > 0$）：

$$\boxed{\int \prod_{i=1}^{L} \frac{d^D l_i}{i\pi^{D/2}} \frac{1}{(\vec{l}^{\,T} A\, \vec{l} - \Delta)^n} = \frac{(-1)^n}{U^{D/2}} \cdot \frac{\Gamma(n - LD/2)}{\Gamma(n)} \cdot \Delta^{LD/2 - n}} \tag{A.1}$$

其中 $U = \det(A)$。

## 附录 B：Wick 缩并公式

对偶数阶张量积分：

$$\boxed{
\begin{aligned}
\int \prod \frac{d^D l_i}{i\pi^{D/2}} \frac{l_{a_1}^{\mu_1} \cdots l_{a_{r}}^{\mu_{r}}}{(\vec{l}^{\,T} A\, \vec{l} - \Delta)^n} = &\frac{(-1)^{n-r/2}}{U^{D/2}} \cdot \frac{\Gamma(n-r/2-LD/2)}{\Gamma(n)} \cdot \frac{\Delta^{LD/2-n+r/2}}{2^{r/2}}\\
&\times \sum_{\text{pairings}} \prod_{\text{pairs}(i,j)} (A^{\text{adj}})_{a_i a_j}\, g^{\mu_i \mu_j}
\end{aligned}
} \tag{B.1}$$

其中 $U = \det(A)$，$r$ 为偶数，求和遍历 $r$ 个指标的所有 $(r-1)!!$ 种配对方式。

> **注意**：此处 $(A^{\text{adj}})_{a_i a_j}$ 等价于 $U \cdot (A^{-1})_{a_i a_j}$，每对缩并贡献一个 $A^{\text{adj}}$ 因子的同时，对应的 $U$ 因子在与前面的 $1/U^{D/2}$ 合并后构成最终结果。
