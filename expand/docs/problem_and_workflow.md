# 二圈Feynman积分的计算：问题定义与工作流程

## 1. 研究目标

本项目的最终目标是**在质数域 $\mathbb{Z}_p$ 下，计算二圈Feynman积分（FI）的二维幂级数展开**。

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

**计算策略**：通过级数展开方法，将FBI表示为关于 $(X, Y)$ 的二维幂级数，然后进行数值积分。

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
I^{\Delta}_{\vec{\nu}+\vec{e}_i} = \frac{1}{\nu_i}\sum_{j}(S^{-1})_{B+i, j} \cdot \begin{pmatrix} -I_{\vec{\nu}}^{\Delta-1} \\ \vdots \\ -I_{\vec{\nu}}^{\Delta-1} \\ I_{\vec{\nu}-\vec{e}_1}^{\Delta-1} \\ \vdots \\ I_{\vec{\nu}-\vec{e}_N}^{\Delta-1} \end{pmatrix}_j
$$

展开写成：

$$
I^{\Delta}_{\vec{\nu}+\vec{e}_i} = \frac{1}{\nu_i} \left[ -\sum_{b=1}^{B} (S^{-1})_{B+i, b} \cdot I_{\vec{\nu}}^{\Delta-1} + \sum_{\alpha=1}^{N} (S^{-1})_{B+i, B+\alpha} \cdot I_{\vec{\nu}-\vec{e}_\alpha}^{\Delta-1} \right]
$$

**情形2：角积分（主积分）**

若维度已等于目标维度 $\Delta_{\text{target}}$，则该积分即为主积分，无需约化。

若 $\Delta \neq \Delta_{\text{target}}$，使用维度迁移：

- **向上迁移**（$\Delta < \Delta_{\text{target}}$，将 $\Delta$ 替换为 $\Delta+1$）：

$$
I_{\vec{\nu}}^{\Delta} = \frac{1}{C}\left[(2\Delta + 2 - \nu - B) z_0 \cdot I_{\vec{\nu}}^{\Delta+1} + \sum_{\alpha=1}^{N} z_\alpha I_{\vec{\nu}-\vec{e}_\alpha}^{\Delta}\right]
$$

- **向下迁移**（$\Delta > \Delta_{\text{target}}$）：

$$
I_{\vec{\nu}}^{\Delta} = \frac{1}{(2\Delta - \nu - B) z_0}\left[C \cdot I_{\vec{\nu}}^{\Delta-1} - \sum_{\alpha=1}^{N} z_\alpha I_{\vec{\nu}-\vec{e}_\alpha}^{\Delta-1}\right]
$$

#### Case 1（$\dim(\text{Null}(S)) = 0$，$C = 0$）

由维度迁移公式（取 $C = 0$）：

$$
I_{\vec{\nu}}^{\Delta} = \frac{-1}{(2\Delta - \nu - B) z_0} \sum_{\alpha=1}^{N} z_\alpha \cdot I_{\vec{\nu}-\vec{e}_\alpha}^{\Delta-1}
$$

#### Case 2（$\dim(\text{Null}(S)) > 0$，$C \neq 0$）

由维度迁移公式（取 $z_0 = 0$）：

$$
I_{\vec{\nu}}^{\Delta} = \frac{1}{C} \sum_{\alpha=1}^{N} z_\alpha \cdot I_{\vec{\nu}-\vec{e}_\alpha}^{\Delta}
$$

注意：这里右边的维度参数仍为 $\Delta$（不变）。

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

$$
\frac{\partial I_{\vec{\nu}}^{\Delta}}{\partial Y} = \sum_{i,j} \left(-\frac{1}{2}\right) \frac{\partial R_{ij}}{\partial Y} \cdot \text{factor}_{ij} \cdot I_{\vec{\nu}+\vec{e}_i+\vec{e}_j}^{\Delta+1}
$$

其中：
$$
\text{factor}_{ij} = \begin{cases}
\nu_i \nu_j & i \neq j \\
\nu_i(\nu_i + 1) & i = j
\end{cases}
$$

右边的 $I_{\vec{\nu}+\vec{e}_i+\vec{e}_j}^{\Delta+1}$ 需要通过约化（第3节）得到主积分的线性组合。

对于主积分 $I_{\vec{\nu}}^{\Delta}$，设其级数为 $f(X,Y) = \sum_{p,q} f_{pq} X^p Y^q$，则微分方程转化为：

**X方向**（$p > 0$）：
$$
f_{p,q} = \frac{1}{p}[\text{rhs}_X]_{p-1,q}
$$

其中 $\text{rhs}_X = \sum_{i,j} (-\frac{1}{2}) \cdot \frac{\partial R_{ij}}{\partial X} \cdot \text{factor}_{ij} \cdot I_{\vec{\nu}+\vec{e}_i+\vec{e}_j}^{\Delta+1}$

**Y方向**（$p = 0, q > 0$）：
$$
f_{0,q} = \frac{1}{q}[\text{rhs}_Y]_{0,q-1}
$$

其中 $\text{rhs}_Y = \sum_{i,j} (-\frac{1}{2}) \cdot \frac{\partial R_{ij}}{\partial Y} \cdot \text{factor}_{ij} \cdot I_{\vec{\nu}+\vec{e}_i+\vec{e}_j}^{\Delta+1}$

## 5. 计算方法：逐阶递推

### 5.1 算法逻辑

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

## 6. 完整计算流程

```
输入：
  - 顶层矩阵 topS(X,Y)（(N+B)×(N+B) 多项式矩阵）
  - 分支数 B = 3，传播子数 N
  - 目标级数度数 d
  - 质数 p（有限域模数）
  - 主积分边界条件（零阶系数）

步骤0：Family/Sector初始化
  0.1 从 topS 构造分支索引 branchIndices
  0.2 计算 dR/dX 和 dR/dY 矩阵
  0.3 枚举所有有效 sector（2^N 个候选，保留覆盖所有分支的）
  0.4 对每个有效 sector：
      - 构造子矩阵 S_sub(X,Y)
      - 执行 RREF 行化简，得到 invS（若满秩）
      - 求解 C 和 z 系数
      - 分类 Case（0/1/2/3）
  0.5 识别所有主积分（Case 0 的角积分）

步骤1：初始化级数求解器
  1.1 设置主积分边界条件 f_{00} = 初始值
  1.2 初始化级数缓存

步骤2：逐阶递推
  for deg = 1 to d:
      for p = 0 to deg:
          q = deg - p
          
          // 2.1 约化递推
          计算微分方程右边所需FBI的deg=deg系数（按需触发）
          
          // 2.2 主积分递推（微分方程）
          for each master FBI:
              if p > 0:
                  使用X方向微分方程计算 f_{p,q}
              else if q > 0:
                  使用Y方向微分方程计算 f_{0,q}

步骤3：输出结果
  - 所有主积分的级数展开
  - 可选：任意FBI的级数展开（通过约化）

输出：
  - FBI级数展开 {f_{pq}} 对于所有 (p,q) 满足 p+q ≤ d
  - 所有系数都是有限域 Z_p 中的元素
```

## 7. 关键实现要点

### 7.1 多项式系数的存储

IBP和微分方程中的多项式系数需要预先计算并存储：
- 使用**哈希表**存储稀疏多项式（只存储非零单项式）
- 多项式结构：`{(x_power, y_power) -> coefficient}`

### 7.2 分子/分母形式的存储

为避免有理函数运算的复杂性，将 $C$, $z$, $S^{-1}$ 等有理函数系数分解为分子和公分母形式：
- `denoCandZ`：$C$ 和 $z$ 系数的公分母多项式
- `numeCandZ[i]`：$C_i$ 或 $z_i$ 乘以公分母后的分子多项式
- `denoInvS[row]`：$S^{-1}$ 第 row 行的公分母
- `numeInvS[row][col]`：$S^{-1}$ 第 (row,col) 元素乘以行公分母后的分子

### 7.3 级数缓存

使用 `map<(nu, delta), Series>` 缓存已计算的FBI级数：
- 避免重复计算
- 支持增量计算（当需要更高度数时继续递推）

### 7.4 数值稳定性

在有限域 $\mathbb{Z}_p$ 中计算时：
- 确保LRR公式中 $D_{00} \neq 0$（选择合适的质数 $p$）
- 使用模逆运算进行除法：$a / b \equiv a \cdot b^{-1} \pmod{p}$

## 8. 参考文献

[1] Li-Hong Huang, Rui-Jun Huang, and Yan-Qing Ma, "Tame multi-leg Feynman integrals beyond one loop," arXiv:2412.21053v1 [hep-ph], 2024.
