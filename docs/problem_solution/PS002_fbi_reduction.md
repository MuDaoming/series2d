---
id: PS002
title: FBI 约化
status: active
created: 2026-03-31
updated: 2026-03-31
depends:
  - B001
parts:
  - PS001
---

# FBI 约化

## 1. 背景

FBI 约化是利用递推关系将一般指标的 FBI 表示为少量主积分（MFBI）的线性组合。

在 FBI 表示中，给定指标向量 $\vec{\nu}$ 的积分记为 $I_{\vec{\nu}}^{\Delta}(X,Y)$，其中 $\Delta$ 是维度参数。

约化的核心工具是 IBP（分部积分）恒等式和维度迁移公式，它们给出了不同 FBI 之间的线性约束关系。

## 2. Dot 约化（Rank = 0 的情况）

首先考虑 $\text{rank} = 0$ 的情况，此时所有指标 $\nu_\alpha \geq 0$。

#### 核心公式

**IBP 递推关系**（论文 Eq. 13）：

$$
S \cdot \begin{pmatrix} t_1^\Delta \\ \vdots \\ t_B^\Delta \\ \nu_1 I_{\vec{\nu}+\vec{e}_1}^{\Delta} \\ \vdots \\ \nu_N I_{\vec{\nu}+\vec{e}_N}^{\Delta} \end{pmatrix} = \begin{pmatrix} -I_{\vec{\nu}}^{\Delta-1} \\ \vdots \\ -I_{\vec{\nu}}^{\Delta-1} \\ I_{\vec{\nu}-\vec{e}_1}^{\Delta-1} \\ \vdots \\ I_{\vec{\nu}-\vec{e}_N}^{\Delta-1} \end{pmatrix}
$$

**维度迁移公式**（论文 Eq. 14/15）：

$$
C \, I_{\vec{\nu}}^{\Delta-1} = (2\Delta - \nu - B) z_0 \, I_{\vec{\nu}}^{\Delta} + \sum_{\alpha=1}^{N} z_\alpha \, I_{\vec{\nu}-\vec{e}_\alpha}^{\Delta-1}
$$

其中 $(C_1, \dots, C_B, z_1, \dots, z_N)^T$ 满足：

$$
S \cdot (C_1, \dots, C_B, z_1, \dots, z_N)^T = (z_0, \dots, z_0, 0, \dots, 0)^T
$$

$C = \sum_{b=1}^B C_b$；若 $\det(S) \neq 0$ 则 $z_0 = 1$，否则 $z_0 = 0$。

**前置处理**：指标为 0 的传播子必须先删除，构造局部 $S$ 矩阵后再进行约化。

#### 四种 Case 的处理

根据 $\det(S)$ 和 $C$ 的值，分四种情况处理：

**Case 1：$\det(S) \neq 0$ 且 $C \neq 0$**

该 Sector 包含唯一的 MFBI（Corner 积分）。

- **Case 1a（已是 Corner，$\nu_\alpha = 1$）**：用 Eq. 14 调整维数
  - 降维（$\Delta > \Delta_{\text{target}}$）：
    $$I_{\vec{\nu}}^{\Delta} = \frac{1}{2\Delta - \nu - B} \left[ C\, I_{\vec{\nu}}^{\Delta-1} - \sum_{\alpha} z_\alpha\, I_{\vec{\nu}-\vec{e}_\alpha}^{\Delta-1} \right]$$
  - 升维（$\Delta < \Delta_{\text{target}}$）：
    $$I_{\vec{\nu}}^{\Delta} = \frac{1}{C} \left[ (2(\Delta+1) - \nu - B)\, I_{\vec{\nu}}^{\Delta+1} + \sum_{\alpha} z_\alpha\, I_{\vec{\nu}-\vec{e}_\alpha}^{\Delta} \right]$$

- **Case 1b（非 Corner，含 dot）**：用 Eq. 13 降 dot

  选取 $\nu_\alpha > 1$ 的传播子 $\alpha$，令 $\vec{\nu} = \vec{\mu} - \vec{e}_\alpha$，解出：
  $$I_{\vec{\mu}}^\Delta = \frac{1}{\mu_\alpha - 1} \left[ S^{-1} \cdot \text{RHS} \right]_{B+\alpha}$$

  结果含 $\Delta - 1$ 的 FBI，递归处理（最终进入 Case 1a）。

**Case 2：$\det(S) \neq 0$ 且 $C = 0$**

使用 Eq. 14（$z_0 = 1$，$C = 0$）：
$$I_{\vec{\nu}}^{\Delta} = -\frac{1}{2\Delta - \nu - B} \sum_{\alpha=1}^{N} z_\alpha \, I_{\vec{\nu}-\vec{e}_\alpha}^{\Delta-1}$$

该 Sector 无 MFBI，所有积分约化到 sub-sector。

**Case 3：$\det(S) = 0$ 且 $C \neq 0$**

使用 Eq. 14（$z_0 = 0$）：
$$I_{\vec{\mu}}^\Delta = \frac{1}{C} \sum_{\alpha=1}^N z_\alpha \, I_{\vec{\mu}-\vec{e}_\alpha}^\Delta$$

RHS 全是更少 dot 的同维数积分。该 Sector 无 MFBI。

**Case 4：$\det(S) = 0$ 且 $C = 0$**

使用 Eq. 14 做指标转移。选 $|z_\beta|$ 最大的 $\beta$：
$$I_{\vec{\mu}}^\Delta = -\sum_{\alpha \neq \beta} \frac{z_\alpha}{z_\beta} \, I_{\vec{\mu}+\vec{e}_\beta-\vec{e}_\alpha}^\Delta$$

将其他传播子的 dot 转移到 $\beta$，最终落入 sub-sector。该 Sector 无 MFBI。

---

## 3. Rank 约化

尽管我们暂时不知道如何处理含 Rank 的费曼积分（作为外层积分），但 FBI 的 IBP 和维度递推关系原则上也能处理含 Rank 的 FBI。

#### 分块矩阵与区域划分

将传播子分为两个区域：

$$\vec{\nu} = (\underbrace{\nu_1, \dots, \nu_{N_{\text{dot}}}}_{\text{dot 区}}, \underbrace{\nu_{N_{\text{dot}}+1}, \dots, \nu_N}_{\text{rank 区}})$$

相应地将 $S$ 矩阵分块：

$$
\begin{pmatrix} S_{\text{dot}} & T \\ T^T & R_{\text{rank}} \end{pmatrix}
$$

#### 特殊情况：$B_r \neq 0$

若存在某些 Branch 完全由 Rank 区传播子构成（记数量为 $B_r$），则可从 $T_1$ 块直接降 Rank：

$$T_1 \cdot (\nu_{\text{rank}} \cdot I_{\vec{\nu}+\vec{e}_{\text{rank}}}^{\Delta}) = (-I_{\vec{\nu}}^{\Delta-1}, \dots, -I_{\vec{\nu}}^{\Delta-1})^T$$

右端是 $(d, r)$ 积分，左端是 $(d, r-1)$ 积分，直接降 Rank。

#### 一般情况：$B_r = 0$

对 $S_{\text{dot}}$ 求解方程：

$$S_{\text{dot}} \cdot (C_1, \dots, C_B, z_1, \dots, z_{N_{\text{dot}}})^T = (z_0, \dots, z_0, 0, \dots, 0)^T$$

得到核心替代方程：

$$C_{\text{dot}} I_{\vec{\nu}}^{\Delta-1} = (2\Delta - \nu - B)z_0 I_{\vec{\nu}}^{\Delta} + \sum_{i=1}^{N_{\text{dot}}} z_i I_{\vec{\nu}-\vec{e}_i}^{\Delta-1} - \vec{Z}_{\text{dot}}^T \, T \, (\nu_{\text{rank}} \cdot I_{\vec{\nu}+\vec{e}_{\text{rank}}}^{\Delta})$$

**Case 1：$\det(S_{\text{dot}}) \neq 0$ 且 $C_{\text{dot}} \neq 0$**

由方程 (I) 解出左向量，代入方程 (II)，得到：

$$\text{Rank 目标} = f(\text{Rank}=d,r, \text{Dot}=d-1,r, \text{Rank}=d,r-1)$$

即用 $(d,r)$、$(d-1,r)$ 和 $(d,r-1)$ 的积分表示 $(d,r+1)$。

**Case 2：$\det(S_{\text{dot}}) \neq 0$ 且 $C_{\text{dot}} = 0$**

$$-(2\Delta - \nu - B) I_{\vec{\nu}}^{\Delta} = \sum_{i=1}^{N_{\text{dot}}} z_i I_{\vec{\nu}-\vec{e}_i}^{\Delta-1} - \vec{Z}_{\text{dot}}^T \, T \, (\nu_{\text{rank}} \cdot I_{\vec{\nu}+\vec{e}_{\text{rank}}}^{\Delta})$$

右端是 $(d-1,r)$ 和 $(d,r-1)$，直接表示 $(d,r)$。

**Case 3：$\det(S_{\text{dot}}) = 0$ 且 $C_{\text{dot}} \neq 0$**

$$C_{\text{dot}} I_{\vec{\nu}}^{\Delta-1} = \sum_{i=1}^{N_{\text{dot}}} z_i I_{\vec{\nu}-\vec{e}_i}^{\Delta-1} - \vec{Z}_{\text{dot}}^T \, T \, (\nu_{\text{rank}} \cdot I_{\vec{\nu}+\vec{e}_{\text{rank}}}^{\Delta})$$

同样是 $(d-1,r)$ 和 $(d,r-1)$ 表示 $(d,r)$。

**Case 4：$\det(S_{\text{dot}}) = 0$ 且 $C_{\text{dot}} = 0$**

不使用核心替代方程，对 Eq. 13 中同一 Branch 内的 dot 传播子行和 rank 传播子行做差：

$$I_{\vec{\nu}-\vec{e}_j}^{\Delta-1} = I_{\vec{\nu}-\vec{e}_i}^{\Delta-1} + \sum_k (R_{jk} - R_{ik}) \nu_k I_{\vec{\nu}+\vec{e}_k}^{\Delta}$$

右端包含 $(d-1,r)$、$(d+1,r)$ 和 $(d,r-1)$。其中 $(d-1,r)$ 和 $(d,r-1)$ 降低了复杂度，$(d+1,r)$ 将复杂度从 rank 转移到 dot。

## 4. 速查表

### Dot 约化

| Case | 条件 | 方程 | 约化效果 | MFBI |
|------|------|------|----------|------|
| 1 | $\det \neq 0, C \neq 0$ | Eq. 13 + Eq. 14 | 降 dot/调维 | 有 |
| 2 | $\det \neq 0, C = 0$ | Eq. 14 | 降 dot + 维数 | 无 |
| 3 | $\det = 0, C \neq 0$ | Eq. 14 | 降 dot | 无 |
| 4 | $\det = 0, C = 0$ | Eq. 14 | 指标转移 | 无 |

### Rank 约化

| Case | 条件 | 约化目标 | 约化结果 |
|------|------|----------|----------|
| 1 | $\det \neq 0, C \neq 0$ | $(d,r+1)$ | $(d,r)$、$(d-1,r)$、$(d,r-1)$ |
| 2 | $\det \neq 0, C = 0$ | $(d,r)$ | $(d-1,r)$、$(d,r-1)$ |
| 3 | $\det = 0, C \neq 0$ | $(d,r)$ | $(d-1,r)$、$(d,r-1)$ |
| 4 | $\det = 0, C = 0$ | $(d,r+1)$ | $(d-1,r)$、$(d+1,r)$、$(d,r-1)$ |
