# FBI 重定义：问题定义与计算流程

## 1. 目标

### 1.1 问题

当前 `IntegrandExpander::getFI2DSeries` 构造 FI 被积函数时，需要：

1. 将 $J \cdot X_0^{e_X} \cdot Y_0^{e_Y} \cdot Z_0^{e_Z} \cdot U^{\nu_{tot}}$ 全部展开为一个大多项式（单项式数量随 $\nu_{tot}$ 爆炸式增长）
2. 用 PDE 递推计算 $U^\gamma$ 的二维级数（$O(\text{deg}^2)$）
3. 做级数×级数卷积 `multiplySeries(uPowerSeries, tmp)`（$O(\text{deg}^4)$）

在 $\text{deg}$ 较大时，步骤2和3是主要瓶颈。

### 1.2 方案

引入**重定义因子** $R$，将上述因子（全部或部分）吸收到 FBI 的定义中：

$$
\widetilde{I}_\nu^D = R(\nu, D) \cdot I_\nu^D
$$

所有递推关系（约化 + 微分方程）改为关于 $\widetilde{I}$ 的形式。代价是每条递推公式多乘一个**多项式幂差因子**，但这远比末端的大卷积便宜。

### 1.3 第一阶段目标

只吸收 $U$ 相关因子（方案 U-only），验证正确性与提速效果。不改动 $X_0 / Y_0 / Z_0$ 相关部分。

## 2. 重定义因子的定义

### 2.1 奇偶对齐函数 $T(D)$

定义**对齐函数** $T$，将递推维度 $D$ 对齐到与输入维度 $D_{in}$ 相同的奇偶性：

$$
T(D) = \begin{cases}
D & \text{if } (D - D_{in}) \text{ 是偶数} \\
D + 1 & \text{if } (D - D_{in}) \text{ 是奇数}
\end{cases}
$$

记 $\bar{D} = T(D)$。

### 2.2 辅助量 $\eta$

定义两个辅助量（整数值，只取 $\{0, 2\}$ 或 $\{-2, 0\}$）：

$$
\eta_-(D) = \bar{D} - \overline{D-1}, \qquad \eta_+(D) = \bar{D} - \overline{D+1}
$$

| $(D - D_{in})$ 的奇偶 | $\eta_-(D)$ | $\eta_+(D)$ |
|:---:|:---:|:---:|
| 偶 | 0 | $-2$ |
| 奇 | 2 | 0 |

### 2.3 重定义因子（U-only 方案）

$$
R(\nu, D) = U(X,Y)^{\text{pow}_U(\nu, \bar{D})}
$$

$$
\text{pow}_U(\nu, \bar{D}) = \nu_{tot} - \frac{(L+1) \bar{D}}{2}
$$

其中 $\nu_{tot} = \sum_\alpha \nu_\alpha$，$L$ 是圈数。

> **注意**：$\text{pow}_U$ 本身可能不是整数（当 $(L+1)\bar{D}$ 是奇数时是半整数）。但**所有需要用到的幂差 $\Delta p_U$ 都是整数**（见第3节），所以不影响实现。

### 2.4 幂差因子

任何递推公式中，源项 $I_{\nu_s}^{D_s}$ 替换为 $\widetilde{I}_{\nu_s}^{D_s}$ 时，系数需要乘以：

$$
U^{\Delta p_U}, \quad \Delta p_U = \text{pow}_U(\nu_t, \bar{D}_t) - \text{pow}_U(\nu_s, \bar{D}_s)
$$

其中 $(\nu_t, D_t)$ 是等式左端的目标积分，$(\nu_s, D_s)$ 是当前源项。

展开：

$$
\Delta p_U = (\nu_{t,tot} - \nu_{s,tot}) - \frac{(L+1)}{2}(\bar{D}_t - \bar{D}_s)
$$

由于 $\bar{D}_t - \bar{D}_s$ 总是偶数（0 或 $\pm 2$），所以 $\Delta p_U$ 总是整数。

## 3. 每种 Case 的幂差公式

以下给出 U-only 方案下每种 Case 每个源项的 $\Delta p_U$。**这些是直接编码的整数公式**。

### 3.1 Case 0 情形1：IBP（非角积分）

目标：$(\nu+e_i, D)$。

| 源项 | $\nu_{t,tot} - \nu_{s,tot}$ | $\bar{D}_t - \bar{D}_s$ | $\Delta p_U$ |
|:---|:---:|:---:|:---:|
| $I_\nu^{D-1}$ | $+1$ | $\eta_-(D)$ | $1 - \frac{L+1}{2}\eta_-(D)$ |
| $I_{\nu-e_\alpha}^{D-1}$ | $+2$ | $\eta_-(D)$ | $2 - \frac{L+1}{2}\eta_-(D)$ |

### 3.2 Case 0 情形2a：向上维度迁移

目标：$(\nu, D)$。

| 源项 | $\nu_{t,tot} - \nu_{s,tot}$ | $\bar{D}_t - \bar{D}_s$ | $\Delta p_U$ |
|:---|:---:|:---:|:---:|
| $I_\nu^{D+1}$ | $0$ | $\eta_+(D)$ | $-\frac{L+1}{2}\eta_+(D)$ |
| $I_{\nu-e_\alpha}^D$ | $+1$ | $0$ | $1$ |

### 3.3 Case 0 情形2b：向下维度迁移

目标：$(\nu, D)$。

| 源项 | $\nu_{t,tot} - \nu_{s,tot}$ | $\bar{D}_t - \bar{D}_s$ | $\Delta p_U$ |
|:---|:---:|:---:|:---:|
| $I_\nu^{D-1}$ | $0$ | $\eta_-(D)$ | $-\frac{L+1}{2}\eta_-(D)$ |
| $I_{\nu-e_\alpha}^{D-1}$ | $+1$ | $\eta_-(D)$ | $1 - \frac{L+1}{2}\eta_-(D)$ |

### 3.4 Case 1（$\dim\text{Null}=0, C=0$）

目标：$(\nu, D)$。

| 源项 | $\Delta p_U$ |
|:---|:---:|
| $I_{\nu-e_\alpha}^{D-1}$ | $1 - \frac{L+1}{2}\eta_-(D)$ |

### 3.5 Case 2（$\dim\text{Null}>0, C\neq 0$）

目标：$(\nu, D)$。

| 源项 | $\Delta p_U$ |
|:---|:---:|
| $I_{\nu-e_\alpha}^D$ | $1$ |

### 3.6 Case 3（$\dim\text{Null}>0, C=0$）

目标：$(\nu, D)$。

| 源项 | $\Delta p_U$ |
|:---|:---:|
| $I_{\nu+e_\beta-e_\alpha}^D$ | $0$ |

### 3.7 二圈（$L=2$）的数值速查

| 情形 | $\Delta p_U$（偶偏移） | $\Delta p_U$（奇偏移） |
|:---|:---:|:---:|
| Case0-IBP, $I_\nu^{D-1}$ | $1$ | $-2$ |
| Case0-IBP, $I_{\nu-e_\alpha}^{D-1}$ | $2$ | $-1$ |
| Case0-up, $I_\nu^{D+1}$ | $3$ | $0$ |
| Case0-up/Case2, $I_{\nu-e_\alpha}^D$ | $1$ | $1$ |
| Case0-down, $I_\nu^{D-1}$ | $0$ | $-3$ |
| Case0-down/Case1, $I_{\nu-e_\alpha}^{D-1}$ | $1$ | $-2$ |
| Case3, $I_{\nu+e_\beta-e_\alpha}^D$ | $0$ | $0$ |

## 4. 重定义后的微分方程

### 4.1 原始微分方程

$$
\partial_X I_\nu^D = \sum_{i,j} \left(-\frac{1}{2}\right) \frac{\partial R_{ij}}{\partial X} \cdot \text{factor}_{ij} \cdot I_{\nu+e_i+e_j}^{D+1}
$$

### 4.2 重定义后

$$
\partial_X \widetilde{I}_\nu^D = \text{pow}_U \cdot \frac{\partial_X U}{U} \cdot \widetilde{I}_\nu^D + \sum_{i,j} \left(-\frac{1}{2}\right) \frac{\partial R_{ij}}{\partial X} \cdot \text{factor}_{ij} \cdot U^{\Delta p_{ij}} \cdot \widetilde{I}_{\nu+e_i+e_j}^{D+1}
$$

其中 $\Delta p_{ij}$ 的目标是 $(\nu, D)$，源是 $(\nu+e_i+e_j, D+1)$：

$$
\Delta p_{ij} = -2 - \frac{L+1}{2}\eta_+(D)
$$

> 二圈：偶偏移 $\Delta p = 1$，奇偏移 $\Delta p = -2$。

### 4.3 乘以 $U$ 消除分母

等式两侧乘以 $U$：

$$
U \cdot \partial_X \widetilde{I}_\nu^D = \text{pow}_U \cdot (\partial_X U) \cdot \widetilde{I}_\nu^D + \sum_{i,j} \left(-\frac{1}{2}\right) \frac{\partial R_{ij}}{\partial X} \cdot \text{factor}_{ij} \cdot U^{\Delta p_{ij} + 1} \cdot \widetilde{I}_{\nu+e_i+e_j}^{D+1}
$$

所有项现在都是**多项式×级数**，没有分母。

### 4.4 逐系数提取

设 $\widetilde{I}_\nu^D = \sum_{p,q} f_{pq} X^p Y^q$。

**X 方向**（$p > 0$）：提取 $(p-1, q)$ 系数：

$$
f_{p,q} = \frac{1}{U_{00} \cdot p} \left( \text{rhsCoeff} + \text{dlogCoeff} - \text{lhsCorrection} \right)
$$

- **rhsCoeff** = $\sum_{i,j} (-\frac{1}{2}) \cdot \text{factor}_{ij} \cdot [(\partial_X R_{ij}) \cdot h_{ij}]_{p-1,q}$，其中 $h_{ij} = U^{\Delta p_{ij}+1} \cdot \widetilde{I}_{\nu+e_i+e_j}^{D+1}$
- **dlogCoeff** = $\text{pow}_U \cdot [(\partial_X U) \cdot \tilde{f}]_{p-1,q}$
- **lhsCorrection** = $\sum_{(a,b) \neq (0,0)} U_{ab} \cdot (p-a) \cdot f_{p-a, q-b}$

**Y 方向**（$p = 0, q > 0$）完全对称。

### 4.5 `pow_U` 的数值

对于主积分（角积分，$D = D_{in}$）：

$$
\text{pow}_U = \nu_{tot} - \frac{(L+1) D_{in}}{2}
$$

这是一个 $\mathbb{Z}_p$ 中的常量。

## 5. IntegrandExpander 的变化

### 5.1 U-only 方案下的效果

FI 被积函数 = $P_0 \cdot U^{\nu_{tot}} \cdot U^\gamma \cdot I_\nu^\Delta$，其中 $P_0 = J \cdot X_0^{e_X} \cdot Y_0^{e_Y} \cdot Z_0^{e_Z}$。

重定义后剩余 $U$ 幂次为 $\nu_{tot} + \gamma - \text{pow}_U = \frac{(L+1)D_F}{2} \cdot \frac{L-2}{2}$。

对于 $L = 2$：**剩余幂次为零**。

> **结论**：二圈 U-only 方案下，$\text{FI} = P_0(X,Y) \cdot \widetilde{I}_\nu^{D_{in}}(X,Y)$，完全消除了 $U^\gamma$ 级数和 $U^{\nu_{tot}}$ 多项式。

### 5.2 新的计算流程

1. 获取重定义后的 FBI 级数 $\widetilde{I}_\nu^{D_{in}}$
2. 逐因子乘以 $J$、$X_0^{e_X}$、$Y_0^{e_Y}$、$Z_0^{e_Z}$（链式 `mulPoly`）
3. 不需要 $U^{\nu_{tot}}$，不需要 $U^\gamma$，不需要 `multiplySeries`

## 6. 参考文献

[1] Li-Hong Huang, Rui-Jun Huang, and Yan-Qing Ma, "Tame multi-leg Feynman integrals beyond one loop," arXiv:2412.21053v1 [hep-ph], 2024.
