# 二圈Feynman积分的计算：问题定义与工作流程

## 1. 研究目标

本项目的最终目标是**在质数域 $\mathbb{Z}_p$ 下，端到端计算二圈Feynman积分（FI）的一维级数展开**。

二圈Feynman积分具有如下形式：

$$
\text{FI} = \int_{X_1+X_2+X_3=1} Q(X_1, X_2, X_3) \cdot I_{\vec{\nu}}^{\Delta}(X_1, X_2, X_3) \, [dX]
$$

其中：
- $Q(X_1, X_2, X_3)$ 是有理函数
- $I_{\vec{\nu}}^{\Delta}(X_1, X_2, X_3)$ 是固定分支积分（Fixed-Branch Integral, FBI）
- 积分在单纯形 $\{X_1, X_2, X_3 \geq 0, X_1+X_2+X_3=1\}$ 上进行

经过变量重标度后，Feynman积分变为：

$$
\text{FI} = \int_0^1 \int_0^1 Q(X, Y) \cdot I_{\vec{\nu}}^{\Delta}(X, Y) \, dX \, dY
$$

**计算策略**：先求重定义后FBI的二维幂级数，再构造FI被积函数二维级数，最后积分得到FI一维级数。

## 2. 数学结构

### 2.1 FBI的定义

固定分支积分（FBI）由参考文献[1]的方程(11)定义。FBI是关于分支参数 $(X, Y)$ 的函数，由传播子指数 $\vec{\nu} = (\nu_1, \ldots, \nu_N)$ 和维度参数 $\Delta$ 标记。

### 2.2 分支结构

传播子按照**分支**（branch）分组，两圈积分分支数 $B = 3$。每个传播子 $\alpha$ 对应一个分支索引 $\text{branch}(\alpha) \in \{1, 2, 3\}$。

分支索引由 $S$ 矩阵的结构决定：传播子 $\alpha$ 属于分支 $b$ 当且仅当 $S_{\alpha, b} = 1$（参见2.3节矩阵 $S$ 的定义）。

### 2.3 矩阵 $S$ 的定义

对于给定的指数向量 $\vec{\nu}$，定义 $(N+B) \times (N+B)$ 对称矩阵 $S$（参考文献[1]方程12）：

$$
S = \begin{pmatrix}
0_{B \times B} & \mathbf{1}_{B \times N}^T \\
\mathbf{1}_{N \times B} & R(X, Y)
\end{pmatrix}
$$

其中：
- $0_{B \times B}$ 是 $B \times B$ 的零矩阵
- $\mathbf{1}_{B \times N}$ 是分支-传播子关联矩阵，$\mathbf{1}_{b, \alpha} = 1$ 当且仅当传播子 $\alpha$ 属于分支 $b$
- $R(X, Y)$ 是 $N \times N$ 的多项式矩阵，元素 $R_{ij}(X, Y)$ 是关于 $(X, Y)$ 的多项式

矩阵 $S$ 是关于 $(X, Y)$ 的多项式矩阵，其结构编码了积分的拓扑信息。

### 2.4 Sector的定义与分类

一个**sector**由非零传播子的集合定义。对于给定的 $\vec{\nu}$，sector 由以下二进制向量表示：

$$
\text{secvec} = (s_1, \ldots, s_N), \quad s_i = \begin{cases} 1 & \nu_i > 0 \\ 0 & \nu_i = 0 \end{cases}
$$

Sector 需要满足**有效性条件**：所有 $B$ 个分支都必须被覆盖（即每个分支至少有一个非零传播子）。

对于有效的 sector，根据其子矩阵 $S_{\text{sub}}$ 的性质分为四类：

| Case | $\dim(\text{Null}(S_{\text{sub}}))$ | $C = \sum_b C_b$ | 物理意义 |
|------|-----------------------------------|------------------|---------|
| 0 | $=0$ | $\neq 0$ | 包含主积分（MFBI） |
| 1 | $=0$ | $=0$ | 无主积分 |
| 2 | $>0$ | $\neq 0$ | 无主积分 |
| 3 | $>0$ | $=0$ | 无主积分 |

**主积分（MFBI）**：Case 0 的角积分（所有 $\nu_i \in \{0, 1\}$）称为主积分（Master FBI），构成约化的基。所有其他FBI都可以表示为MFBI的线性组合（系数为有理函数）。

### 2.5 C 和 z 系数

对于每个 sector，通过求解线性方程组 $S_{\text{sub}} \cdot (C_1, \ldots, C_B, z_1, \ldots, z_N)^T = (1, \ldots, 1, 0, \ldots, 0)^T$ 得到：

- **$C_b$**：与分支 $b$ 关联的系数（$b = 1, \ldots, B$）
- **$z_\alpha$**：与传播子 $\alpha$ 关联的系数（$\alpha = 1, \ldots, N$）
- **$C = \sum_b C_b$**：总 $C$ 系数

当 $\dim(\text{Null}(S_{\text{sub}})) > 0$ 时，$(C_b, z_\alpha)$ 取自零空间的一个基向量。

## 3. 约化

### 3.1 两个核心关系

#### IBP递推关系（方程13）

Integration-by-Parts（IBP）关系给出了不同FBI之间的线性约束：

$$
S \cdot \begin{pmatrix} t_1 \\ \vdots \\ t_B \\ \nu_1 I_{\vec{\nu}+\vec{e}_1}^{\Delta} \\ \vdots \\ \nu_N I_{\vec{\nu}+\vec{e}_N}^{\Delta} \end{pmatrix} = \begin{pmatrix} -I_{\vec{\nu}}^{\Delta-1} \\ \vdots \\ -I_{\vec{\nu}}^{\Delta-1} \\ I_{\vec{\nu}-\vec{e}_1}^{\Delta-1} \\ \vdots \\ I_{\vec{\nu}-\vec{e}_N}^{\Delta-1} \end{pmatrix}
$$

其中 $t_b$ 是辅助变量，$\vec{e}_\alpha$ 是第 $\alpha$ 个单位向量。

当 $\dim(\text{Null}(S)) = 0$ 时，可以通过 $S^{-1}$ 求解任意 $I_{\vec{\nu}+\vec{e}_\alpha}^{\Delta}$。

#### 维度迁移关系（方程14）

维度迁移公式允许在不同维度 $\Delta$ 之间转换：

$$
C \cdot I_{\vec{\nu}}^{\Delta-1} = (2\Delta - \nu - B) z_0 \cdot I_{\vec{\nu}}^{\Delta} + \sum_{\alpha=1}^{N} z_\alpha I_{\vec{\nu}-\vec{e}_\alpha}^{\Delta-1}
$$

其中：
- $\nu = \sum_\alpha \nu_\alpha$ 是传播子指数之和
- $z_0 = 1$ 当 $\dim(\text{Null}(S)) = 0$，否则 $z_0 = 0$

### 3.2 四种Case的约化算法

#### Case 0（$\dim(\text{Null}(S)) = 0$，$C \neq 0$）

**情形1：非角积分**（存在 $\nu_i > 1$）

首先使用IBP关系。选取最大指数 $\nu_i$，解出 $I^{\Delta}_{\vec{\nu}}$：

$$
I^{\Delta}_{\vec{\nu}+\vec{e}_i} = \frac{1}{\nu_i} \left[ -\sum_{b=1}^{B} (S^{-1})_{B+i, b} \cdot I_{\vec{\nu}}^{\Delta-1} + \sum_{\alpha=1}^{N} (S^{-1})_{B+i, B+\alpha} \cdot I_{\vec{\nu}-\vec{e}_\alpha}^{\Delta-1} \right]
$$

**情形2：角积分（主积分）**

若维度已等于目标维度 $\Delta_{\text{target}}$，则该积分即为主积分，无需约化。

若 $\Delta \neq \Delta_{\text{target}}$，使用维度迁移：

- **向上迁移**（$\Delta < \Delta_{\text{target}}$）：

$$
I_{\vec{\nu}}^{\Delta} = \frac{1}{C}\left[(2\Delta + 2 - \nu - B) z_0 \cdot I_{\vec{\nu}}^{\Delta+1} + \sum_{\alpha=1}^{N} z_\alpha I_{\vec{\nu}-\vec{e}_\alpha}^{\Delta}\right]
$$

- **向下迁移**（$\Delta > \Delta_{\text{target}}$）：

$$
I_{\vec{\nu}}^{\Delta} = \frac{1}{(2\Delta - \nu - B) z_0}\left[C \cdot I_{\vec{\nu}}^{\Delta-1} - \sum_{\alpha=1}^{N} z_\alpha I_{\vec{\nu}-\vec{e}_\alpha}^{\Delta-1}\right]
$$

#### Case 1（$\dim(\text{Null}(S)) = 0$，$C = 0$）

$$
I_{\vec{\nu}}^{\Delta} = \frac{-1}{(2\Delta - \nu - B) z_0} \sum_{\alpha=1}^{N} z_\alpha \cdot I_{\vec{\nu}-\vec{e}_\alpha}^{\Delta-1}
$$

#### Case 2（$\dim(\text{Null}(S)) > 0$，$C \neq 0$）

$$
I_{\vec{\nu}}^{\Delta} = \frac{1}{C} \sum_{\alpha=1}^{N} z_\alpha \cdot I_{\vec{\nu}-\vec{e}_\alpha}^{\Delta}
$$

#### Case 3（$\dim(\text{Null}(S)) > 0$，$C = 0$）

选取任意非零 $z_\beta \neq 0$：

$$
I_{\vec{\nu}}^{\Delta} = \frac{-1}{z_\beta} \sum_{\alpha \neq \beta} z_\alpha \cdot I_{\vec{\nu}+\vec{e}_\beta-\vec{e}_\alpha}^{\Delta}
$$

## 4. 微分方程

MFBI作为 $(X, Y)$ 的函数满足耦合的一阶线性偏微分方程系统：

$$
\frac{\partial I_{\vec{\nu}}^{\Delta}}{\partial X} = \sum_{i,j} \left(-\frac{1}{2}\right) \frac{\partial R_{ij}}{\partial X} \cdot \text{factor}_{ij} \cdot I_{\vec{\nu}+\vec{e}_i+\vec{e}_j}^{\Delta+1}
$$

其中：
$$
\text{factor}_{ij} = \begin{cases}
\nu_i \nu_j & i \neq j \\
\nu_i(\nu_i + 1) & i = j
\end{cases}
$$

右边的 $I_{\vec{\nu}+\vec{e}_i+\vec{e}_j}^{\Delta+1}$ 需要通过约化（第3节）得到主积分的线性组合。

## 5. 计算方法：逐阶递推

### 5.1 算法逻辑

假设当前已知：
- 主积分的 $\deg \leq N$ 的项
- 其它所需积分的 $\deg < N$ 的项

则递推步骤为：

1. **约化步**：通过约化公式（第3节）计算其它所需积分的 $\deg = N$ 的项
2. **微分步**：利用微分方程（第4节）计算主积分导数的 $\deg = N$ 的项，即主积分 $\deg = N+1$ 的项

### 5.2 核心思想

将IBP约化关系和微分方程统一为**线性递推关系**（LRR），逐阶计算级数系数。

**关键特点**：
- 所有中间量都是**二维幂级数**而非有理函数
- IBP和微分方程的系数是**多项式**（关于 $X, Y$）
- 每一阶的计算都是纯数值运算（在有限域 $\mathbb{Z}_p$ 中）

### 5.3 线性递推关系（LRR）的统一形式

所有约化公式和微分方程都可以写成LRR形式：

$$
D(X,Y) \cdot g(X,Y) = \sum_i N_i(X,Y) \cdot f_i(X,Y)
$$

其中：
- $D(X,Y)$ 是公分母多项式
- $N_i(X,Y)$ 是分子多项式
- $g(X,Y)$ 是待求级数
- $f_i(X,Y)$ 是已知级数

### 5.4 LRR递推公式

对于总度数 $p+q = d$ 的系数 $g_{pq}$：

$$
g_{pq} = \frac{1}{D_{00}} \left( \sum_i [N_i \cdot f_i]_{pq} - \sum_{\substack{(a,b) \neq (0,0) \\ a+b \leq \deg(D)}} D_{ab} \cdot g_{p-a,q-b} \right)
$$

其中：
- $[N_i \cdot f_i]_{pq}$ 是多项式与级数卷积的第 $(p,q)$ 项：$[N_i \cdot f_i]_{pq} = \sum_{a,b} (N_i)_{ab} \cdot (f_i)_{p-a,q-b}$
- $g_{p-a,q-b}$ 是已知的（$p-a+q-b < p+q$）
- 要求 $D_{00} \neq 0$（常数项非零）

## 6. FBI 重定义

### 6.1 动机

原始的 FI 被积函数包含 $U^{\gamma} \cdot U^{\nu_{tot}} \cdot I$，其中 $U^{\gamma}$ 需要通过 PDE 递推展开为二维级数，然后与 FBI 级数做级数×级数卷积（$O(\text{deg}^4)$）。这在高阶时是主要瓶颈。

通过引入重定义因子 $R$，将 $U$ 的幂次吸收进 FBI 的定义中，可以完全消除这个瓶颈。

### 6.2 重定义因子

定义重定义后的 FBI：

$$
\widetilde{I}_\nu^D = U^{\text{pow}_U(\nu, \bar{D})} \cdot I_\nu^D
$$

其中：

$$
\text{pow}_U(\nu, \bar{D}) = \nu_{tot} - \frac{(L+1) \bar{D}}{2}
$$

$\bar{D} = T(D)$ 是对齐函数，将 $D$ 对齐到与 $D_{in}$ 相同的奇偶性：

$$
T(D) = \begin{cases}
D & \text{if } (D - D_{in}) \text{ 是偶数} \\
D + 1 & \text{if } (D - D_{in}) \text{ 是奇数}
\end{cases}
$$

### 6.3 约化关系的变化

所有约化关系中，源项 $I_{\nu_s}^{D_s}$ 替换为 $\widetilde{I}_{\nu_s}^{D_s}$ 后，系数需要乘以幂差因子：

$$
U^{\Delta p_U}, \quad \Delta p_U = \text{pow}_U(\nu_t, \bar{D}_t) - \text{pow}_U(\nu_s, \bar{D}_s)
$$

$\Delta p_U$ 总是整数。为消除负幂次，等式两侧乘以 $U^m$（$m = \max(0, -\min(\Delta p_i))$），得到：

$$
D \cdot U^m \cdot \widetilde{I}_T = \sum_i (N_i \cdot U^{\Delta p_i + m}) \cdot \widetilde{I}_{S_i}
$$

仍为 LRR 形式，只是 $D$ 和 $N_i$ 变成了更高次的多项式。

### 6.4 微分方程的变化

重定义后的微分方程（乘以 $U$ 消除分母后）：

$$
U \cdot \partial_X \widetilde{I} = \text{pow}_U \cdot (\partial_X U) \cdot \widetilde{I} + \sum_{i,j} \left(-\frac{1}{2}\right) (\partial_X R_{ij} \cdot U^{\Delta p+1}) \cdot \text{factor}_{ij} \cdot \widetilde{I}_{\nu+e_i+e_j}^{D+1}
$$

逐系数提取的递推公式：

$$
\widetilde{I}_{p,q} = \frac{1}{U_{00} \cdot p} \left( \text{rhsCoeff} + \text{dlogCoeff} - \text{lhsCorrection} \right)
$$

- **rhsCoeff**：源项的贡献（使用 $\partial_X R_{ij} \cdot U^L$ 代替原始 $\partial_X R_{ij}$）
- **dlogCoeff** = $\text{pow}_U \cdot [(\partial_X U) \cdot \widetilde{I}]_{p-1,q}$
- **lhsCorrection** = $\sum_{(a,b) \neq (0,0)} U_{ab} \cdot (p-a) \cdot \widetilde{I}_{p-a, q-b}$

### 6.5 二圈的简化

对于 $L = 2$，FI 被积函数中 $U$ 的剩余幂次为：

$$
\nu_{tot} + \gamma - \text{pow}_U = \frac{(L+1)D_F}{2} \cdot \frac{L-2}{2} = 0
$$

**结论**：二圈时，$\text{FI} = J \cdot W \cdot \widetilde{I}_\nu^{D_{in}}$，完全消除了 $U^{\gamma}$ 级数和 $U^{\nu_{tot}}$ 多项式，不再需要级数×级数卷积。

## 7. 从重定义FBI到FI一维展开

### 7.1 变量替换

原始三变量满足 $X + Y + Z = 1$。先做二维化替换：

$$
X \to X, \quad Y \to (1-X)Y, \quad Z \to 1 - X - (1-X)Y
$$

再做局部平移：$X \to X + a, \quad Y \to Y + b$。

所有多项式（$U$、Jacobian、$X_0/Y_0/Z_0$ 幂次）都按此顺序处理。

### 7.2 FI 被积函数的构造

对于目标 $\vec{\nu}$，FI 被积函数的二维级数为：

$$
\text{FI}(X,Y) = P(X,Y) \cdot \widetilde{I}_\nu^{D_{in}}(X,Y)
$$

其中多项式因子 $P = J \cdot W$：
- $J = 1 - X$（Jacobian，变量替换后平移）
- $W = X_0^{e_X} \cdot Y_0^{e_Y} \cdot Z_0^{e_Z}$（参数幂次，$e_i = \nu_i - 1$）

计算步骤：
1. 构造平移后的多项式 $P(X,Y)$
2. 获取重定义后的 FBI 级数 $\widetilde{I}_\nu^{D_{in}}$
3. 用 `mulPoly` 计算 $P \cdot \widetilde{I}$（多项式×级数，$O(\text{deg}^2 \cdot \text{多项式项数})$）

### 7.3 二维到一维积分

对二维级数 $\sum_{p+q \leq d} c_{pq} X^p Y^q$ 积分得到一维级数 $\sum_d s_d \epsilon^d$：

$$
s_d = \sum_{p+q=d} c_{pq} \cdot w(p, q)
$$

其中 $w(p, q)$ 是单项式 $X^p Y^q$ 在积分域上的权重（由平移参数 $a, b$ 决定）。

## 8. 完整计算流程

```
输入（四个路径参数）：
  - S_path: 已平移后的 topS(X,Y)（维度应为 B+N）
  - config_path: 包含 N,B,deg,p,d,a,b,bc
  - target_path: 目标 nu 列表（每行一个 {nu_1,...,nu_N}）
  - output_path: 输出文件路径

步骤0：Family/Sector 初始化
  0.1 从 topS 构造分支索引 branchIndices
  0.2 计算 dR/dX 和 dR/dY
  0.3 枚举有效 sector 并分类 Case（0/1/2/3）
  0.4 识别主积分集合（Case 0 的角积分）

步骤1：Redefinition 初始化
  1.1 构造平移后的 U(X,Y) 多项式
  1.2 创建 Redefinition 结构体（L, D_in, shiftedU, dUdX, dUdY）
  1.3 预计算 dRdX·U^L 和 dRdY·U^L（用于重定义后微分方程）
  1.4 预计算每个主积分的 pow_U 标量值

步骤2：重定义FBI 二维级数求解（SeriesSolver）
  2.1 根据 bc 设置主积分零阶边界
  2.2 逐阶递推：
      - 微分方程求解主积分（使用重定义后的公式）
      - 按需触发约化（约化关系中乘以 U^{Δp} 幂差因子）

步骤3：FI 被积函数构造（IntegrandExpander）
  3.1 构造 FI 多项式 P = J·W（不含 U 因子）
  3.2 用 mulPoly 计算 P · Ĩ 得到 FI 二维级数

步骤4：二维到一维积分（SeriesIntegrator）
  4.1 对二维系数积分得到一维系数 {c0,...,c_deg}
  4.2 按 target 的 nu 顺序写入 output

输出：
  - 单文件文本
  - 每行一个 nu 对应的一维级数：{c0,...,c_deg}
```

## 9. 参考文献

[1] Li-Hong Huang, Rui-Jun Huang, and Yan-Qing Ma, "Tame multi-leg Feynman integrals beyond one loop," arXiv:2412.21053v1 [hep-ph], 2024.
