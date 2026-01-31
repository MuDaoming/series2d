# 二圈Feynman积分的计算：问题定义与流程

## 1. 研究目标

本项目的最终目标是**在质数域 $\mathbb{Z}_p$ 下，计算二圈Feynman积分（FI） 渐进展开**。

二圈Feynman积分具有如下形式：

$$
\text{FI} = \int_{X_1+X_2+X_3=1} Q(X_1, X_2, X_3) \cdot I_{\vec{\nu}}^{\Delta}(X_1, X_2, X_3) \, [dX]
$$

其中：
- $Q(X_1, X_2, X_3)$ 是有理函数
- $I_{\vec{\nu}}^{\Delta}(X_1, X_2, X_3)$ 是固定分支积分（Fixed-Branch Integral, FBI）
- 积分在单纯形 $\{X_1, X_2, X_3 \geq 0, X_1+X_2+X_3=1\}$ 上进行

**第一步：处理狄拉克函数。** 原始FI包含delta函数约束 $\delta(1-X_1-X_2-X_3)$，这使得积分实际上在二维单纯形上进行。通过消去delta函数，我们将三重积分化为二重积分：

$$
\text{FI}_{\vec{\nu}}^{\Delta} = \int_0^1 dX_1 \int_0^{1-X_1} dX_2 \, U^{\nu-\frac{L+1}{2}D} X_1^{\nu_{b_1}} X_2^{\nu_{b_2}} (1-X_1-X_2)^{\nu_{b_3}} \sum_{\Delta',\vec{\nu}'} Q(X_1,X_2,1-X_1-X_2) \text{FBI}_{\vec{\nu}'}^{\Delta'}(X_1,X_2,1-X_1-X_2)
$$

其中 $U=X_1X_2+X_2X_3+X_3X_1$ 是二圈动能项。

进一步通过变量重标度将积分域变为矩形 $[0,1]^2$。引入新变量 $(X, Y)$：

$$
X_1 = X, \quad X_2 = Y(1-X), \quad X_3 = (1-Y)(1-X)
$$

雅可比因子为 $1-X$。变换后Feynman积分变为：

$$
\text{FI} = \int_0^1 \int_0^1 \tilde{Q}(X, Y) \cdot I_{\vec{\nu}}^{\Delta}(X, Y) \, dX \, dY
$$

其中 $\tilde{Q}(X, Y) = Q(X, Y(1-X), (1-Y)(1-X)) \cdot (1-X)$ 包含了雅可比因子。

**第二步：引入渐进展开参数。** 为了得到FBI在相空间特定点 $(a, b)$ 附近的渐进展开，我们引入小参数 $\delta \in (0,1)$，将积分限制在点 $(a,b)$ 的邻域内：

$$
\text{FI}_{\vec{\nu}}^{\Delta} \approx \int_{a-a\delta}^{a+(1-a)\delta} dX \int_{b-b\delta}^{b+(1-b)\delta} dY \, Q_{\vec{\nu}'}(X, Y) \text{FBI}_{\vec{\nu}'}^{\Delta'}(X, Y)
$$

当 $\delta \to 0$ 时，FBI在 $(a,b)$ 附近的Taylor级数展开主导积分行为。这允许我们通过级数展开精确计算渐进贡献，最后取 $\delta = 1$ 恢复完整积分。

**计算策略：** 通过以下步骤实现：
1. **变量重标度**：将三重积分（含delta函数）变换为 $[0,1]^2$ 上的二重积分
2. **FBI约化**：将FBI表示为主积分（MFBI）的线性组合
3. **微分方程求解**：求解MFBI的二维级数展开
4. **级数运算**：计算FBI的级数展开
5. **数值积分**：对FBI级数进行数值积分，得到FI的值

本文档完备地定义问题的数学结构，并说明完整的计算流程。

## 2. 数学结构

### 2.1 FBI的定义

固定分支积分（FBI）由参考文献[1]的方程(11)定义。对于给定的分支参数 $\mathbf{X} = (X_1, X_2, \ldots, X_B)$（$B$ 为分支数），FBI定义为：

$$
I_{\vec{\nu}}^{\Delta}(\mathbf{X}) = \frac{(-1)^{\nu} \Gamma(\nu - \Delta)}{\prod_{\alpha=1}^{N} \Gamma(\nu_\alpha)} \int [dy] \prod_{\alpha=1}^{N} y_\alpha^{\nu_\alpha - 1} \left( \frac{1}{2} \mathbf{y}^T \cdot R(\mathbf{X}) \cdot \mathbf{y} - i0^+ \right)^{\nu - \Delta}
$$

其中：
- $\vec{\nu} = (\nu_1, \nu_2, \ldots, \nu_N)$ 是传播子指数向量，$N$ 是传播子数量
- $\nu = \sum_{\alpha=1}^{N} \nu_\alpha$ 是总指数
- $\Delta$ 是与时空维度相关的参数
- $R(\mathbf{X})$ 是关于分支参数的 $L+1$ 次齐次多项式矩阵（$L$ 为圈数）
- 积分测度：$[dy] = \prod_{\alpha=1}^{N} dy_\alpha \prod_{b=1}^{B} \delta\left(1 - \sum_{i=1}^{n_b} y_{(b,i)}\right)$

在FBI内部，Feynman参数 $\mathbf{y}$ 已经通过delta函数约束，FBI被视为分支参数 $\mathbf{X}$ 的函数。经过重标度后，FBI表示为 $(X, Y)$ 的函数，矩阵 $R(\mathbf{X})$ 相应地表示为 $R(X, Y)$。

### 2.2 分支结构

传播子按照**分支**（branch）分组。属于第 $b$ 个分支的传播子记为 $D_{(b,1)}, D_{(b,2)}, \ldots, D_{(b,n_b)}$，其中 $n_b$ 是第 $b$ 个分支的传播子数量，满足：

$$
\sum_{b=1}^{B} n_b = N
$$

每个传播子 $\alpha$ 对应一个分支索引 $\text{branch}(\alpha) \in \{1, 2, \ldots, B\}$。

### 2.3 矩阵 $S$ 的定义

对于给定的指数向量 $\vec{\nu}$，定义 $(N+B) \times (N+B)$ 对称矩阵 $S$（参考文献[1]方程12）：

1. 若 $\alpha > B$ 且 $\beta > B$，则 $S_{\alpha,\beta} = R_{\alpha-B, \beta-B}$
2. 若 $\alpha \leq B$ 且传播子 $\beta-B$ 属于分支 $\alpha$，则 $S_{\alpha,\beta} = 1$（对称位置同理）
3. 其他情况 $S_{\alpha,\beta} = 0$

例如，当 $B=3$，$(n_1, n_2, n_3) = (2, 1, 1)$ 时，矩阵 $S$ 的形式为：

$$
S = \begin{pmatrix}
0_{3 \times 3} & \begin{matrix} 1 & 1 & 0 & 0 \\ 0 & 0 & 1 & 0 \\ 0 & 0 & 0 & 1 \end{matrix} \\
\begin{matrix} 1 & 0 & 0 \\ 1 & 0 & 0 \\ 0 & 1 & 0 \\ 0 & 0 & 1 \end{matrix} & R
\end{pmatrix}
$$

这里 $R$ 是 $N \times N$ 的矩阵。

### 2.4 Sector的定义与分类

一个**sector**由非零传播子的集合定义，可用二进制向量 $\text{secvec} = (s_1, s_2, \ldots, s_N)$ 表示，其中 $s_i = 1$ 当且仅当 $\nu_i > 0$。

对每个sector，通过求解线性方程组：

$$
S \cdot (C_1, \ldots, C_B, z_1, \ldots, z_N)^T = (z_0, \ldots, z_0, 0, \ldots, 0)^T
$$

得到常数 $(C_b, z_\alpha, z_0)$，其中 $C = \sum_{b=1}^{B} C_b$。

根据 $\dim(\text{Null}(S))$ 和 $C$ 的取值，sector分为四类（参考文献[1]）：

| Case | $\dim(\text{Null}(S))$ | $C$ | 是否有MFBI |
|------|----------------------|-----|-----------|
| 0（论文Case 1） | $=0$ | $\neq 0$ | 是（唯一） |
| 1（论文Case 2） | $=0$ | $=0$ | 否 |
| 2（论文Case 3） | $>0$ | $\neq 0$ | 否 |
| 3（论文Case 4） | $>0$ | $=0$ | 否 |

**主积分（MFBI）：** Case 0的角积分（所有 $\nu_i=1$）称为主积分（Master FBI, MFBI），构成约化的基。

### 2.5 FBI Family

**FBI Family** 是指对于给定的Feynman图拓扑结构，所有可能的传播子指数组合 $\{\vec{\nu}\}$ 所对应的FBI集合。Family包含：

- 顶层矩阵 $R(\mathbf{X})$（或重标度后的 $R(X,Y)$）
- 所有可能的sector（$2^N$ 个）
- 每个sector的分类（Case 0/1/2/3）
- 所有主积分MFBI的标识

Family的作用是提供全局信息，用于FBI约化和微分方程的构建。

## 3. 计算流程

本项目的完整计算流程分为四个核心步骤。

### 步骤1：FBI约化

**目标：** 将任意FBI $I_{\vec{\nu}}^{\Delta}(X, Y)$ 表示为MFBI的线性组合。

#### 3.1 约化原理

基于分部积分（IBP）恒等式（参考文献[1]方程13和14），我们有两个基本关系：

**递推关系式**（方程13）：
$$
S \cdot (t_1, \ldots, t_B, \nu_1 I_{\vec{\nu}+\vec{e}_1}^{\Delta}, \ldots, \nu_N I_{\vec{\nu}+\vec{e}_N}^{\Delta})^T = (-I_{\vec{\nu}}^{\Delta-1}, \ldots, -I_{\vec{\nu}}^{\Delta-1}, I_{\vec{\nu}-\vec{e}_1}^{\Delta-1}, \ldots, I_{\vec{\nu}-\vec{e}_N}^{\Delta-1})^T
$$

**维度迁移关系式**（方程14）：
$$
C \cdot I_{\vec{\nu}}^{\Delta-1} = (2\Delta - \nu - B) z_0 \cdot I_{\vec{\nu}}^{\Delta} + \sum_{\alpha=1}^{N} z_\alpha I_{\vec{\nu}-\vec{e}_\alpha}^{\Delta-1}
$$

#### 4.2 约化策略（按Case分类）

##### Case 0：$\dim(\text{Null}(S)) = 0$，$C \neq 0$

这是包含MFBI的情形。约化分两步：

1. **降低指数**：使用递推关系（方程13）和 $S^{-1}$，将 $I_{\vec{\nu}}^{\Delta}$ 约化到角积分和子sector
   
   $$
   \nu_{\max} I_{\vec{\nu}}^{\Delta} = \sum_{j} (S^{-1})_{\text{row}, j} \cdot (\text{rhs})_j
   $$
   
   递归应用直至所有 $\nu_i = 1$（角积分）。

2. **维度迁移**：使用方程14调整角积分到目标维度 $\Delta_0$
   
   $$
   I_{\vec{\nu}}^{\Delta} = \frac{C}{(2\Delta - \nu - B) z_0} I_{\vec{\nu}}^{\Delta-1} - \frac{1}{(2\Delta - \nu - B) z_0} \sum_{\alpha} z_\alpha I_{\vec{\nu}-\vec{e}_\alpha}^{\Delta-1}
   $$

##### Case 1：$\dim(\text{Null}(S)) = 0$，$C = 0$

$$
(2\Delta - \nu - B) I_{\vec{\nu}}^{\Delta} = -\sum_{\alpha=1}^{N} z_\alpha I_{\vec{\nu}-\vec{e}_\alpha}^{\Delta-1}
$$

降低指数和维度，最终约化到子sector。

##### Case 2：$\dim(\text{Null}(S)) > 0$，$C \neq 0$

取 $z_0 = 0$：
$$
C \cdot I_{\vec{\nu}}^{\Delta-1} = \sum_{\alpha=1}^{N} z_\alpha I_{\vec{\nu}-\vec{e}_\alpha}^{\Delta-1}
$$

降低指数，约化到子sector。

##### Case 3：$\dim(\text{Null}(S)) > 0$，$C = 0$

取 $z_0 = 0$，选择最大的 $|z_\beta|$：
$$
I_{\vec{\nu}}^{\Delta} = -\sum_{\alpha \neq \beta} \frac{z_\alpha}{z_\beta} I_{\vec{\nu}+\vec{e}_\beta-\vec{e}_\alpha}^{\Delta}
$$

调整指数组合，约化到子sector。

#### 4.3 约化结果

约化的最终结果是有理函数系数 $r_k^{(\vec{\nu}, \Delta)}(X, Y)$，使得：

$$
I_{\vec{\nu}}^{\Delta}(X, Y) = \sum_{k=1}^{M} r_k^{(\vec{\nu}, \Delta)}(X, Y) \cdot f_k(X, Y)
$$

其中 $f_k(X, Y) = I_{\text{MFBI}_k}^{\Delta_0}(X, Y)$ 是 $M$ 个MFBI，$r_k$ 是有理函数。

### 步骤2：微分方程求解

**目标：** 求解所有MFBI的级数展开 $f_k(X, Y) = \sum c_{mn}^{(k)} X^m Y^n$。

#### 4.4 微分方程系统的构建

在求解MFBI的级数展开之前，需要首先构建微分方程系统的系数矩阵 $A_X$ 和 $A_Y$。

**构建步骤：**

1. **计算FBI对分支参数的导数**

对于FBI $I_{\vec{\nu}}^{\Delta}(X,Y)$，其对 $X$ 和 $Y$ 的导数可由以下公式给出：

$$
\frac{\partial I_{\vec{\nu}}^{\Delta}}{\partial X} = \sum_{i=1}^{N} \sum_{j=1}^{N} \left(-\frac{1}{2}\right) \frac{\partial R_{ij}}{\partial X} \cdot \text{factor}_{ij} \cdot I_{\vec{\nu}+\vec{e}_i+\vec{e}_j}^{\Delta+1}
$$

其中系数：
$$
\text{factor}_{ij} = \begin{cases}
\nu_i \nu_j & \text{if } i \neq j \\
\nu_i(\nu_i + 1) & \text{if } i = j
\end{cases}
$$

对 $Y$ 的导数公式类似，只需将 $\partial R/\partial X$ 替换为 $\partial R/\partial Y$。

2. **对每个MFBI计算导数**

对于第 $k$ 个MFBI $f_k(X,Y) = I_{\vec{\nu}_k}^{\Delta_0}(X,Y)$，使用上述公式计算：
- $\partial f_k / \partial X$ 得到FBI的线性组合
- $\partial f_k / \partial Y$ 得到FBI的线性组合

3. **FBI约化到MFBI**

导数中出现的FBI都需要约化为MFBI的线性组合：
$$
I_{\vec{\nu}}^{\Delta}(X,Y) = \sum_{j=1}^{M} r_j^{(\vec{\nu},\Delta)}(X,Y) \cdot f_j(X,Y)
$$

4. **组合得到微分方程系数**

通过组合导数公式中的系数和约化系数，得到：
$$
\frac{\partial f_k}{\partial X} = \sum_{j=1}^{M} [A_X]_{kj} \cdot f_j, \quad \frac{\partial f_k}{\partial Y} = \sum_{j=1}^{M} [A_Y]_{kj} \cdot f_j
$$

**最终的微分方程系统：**

MFBI作为 $(X, Y)$ 的函数满足耦合的一阶线性微分方程系统：

$$
\begin{cases}
\displaystyle\frac{\partial \mathbf{f}}{\partial X} = A_X(X, Y) \cdot \mathbf{f} \\[10pt]
\displaystyle\frac{\partial \mathbf{f}}{\partial Y} = A_Y(X, Y) \cdot \mathbf{f}
\end{cases}
$$

其中 $\mathbf{f} = (f_1, \ldots, f_M)^T$，$A_X$ 和 $A_Y$ 是 $M \times M$ 有理函数矩阵：

$$
[A_X]_{ij} = \frac{P_{ij}^{(X)}(X, Y)}{Q_{ij}^{(X)}(X, Y)}, \quad [A_Y]_{ij} = \frac{P_{ij}^{(Y)}(X, Y)}{Q_{ij}^{(Y)}(X, Y)}
$$

通过选择合适的基（三角化），可以将系统解耦为标准形式。

#### 4.5 标准微分方程

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

#### 4.6 级数递推公式

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

#### 4.7 求解结果

得到所有MFBI的级数展开：

$$
f_k(X, Y) = \sum_{m+n \leq d} c_{mn}^{(k)} X^m Y^n, \quad k = 1, \ldots, M
$$

### 步骤3：FBI级数计算

**目标：** 结合约化系数和MFBI级数，得到所有FBI的级数展开。

#### 4.8 级数与有理函数的乘法

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

#### 4.9 最终FBI级数

利用步骤1的约化结果和步骤2的MFBI级数：

$$
I_{\vec{\nu}}^{\Delta}(X, Y) = \sum_{k=1}^{M} r_k^{(\vec{\nu}, \Delta)}(X, Y) \cdot f_k(X, Y)
$$

对每项执行级数与有理函数乘法，然后求和，得到：

$$
I_{\vec{\nu}}^{\Delta}(X, Y) = \sum_{m+n \leq d} I_{mn}^{(\vec{\nu}, \Delta)} X^m Y^n
$$

### 步骤4：Feynman积分的数值计算

**目标：** 计算Feynman积分的数值。

利用步骤3得到的FBI级数展开和有理函数 $\tilde{Q}(X, Y)$，Feynman积分可表示为：

$$
\text{FI} = \int_0^1 \int_0^1 \tilde{Q}(X, Y) \cdot \sum_{m,n} I_{mn}^{(\vec{\nu}, \Delta)} X^m Y^n \, dX \, dY
$$

利用幂次积分的解析公式：

$$
\int_0^1 \int_0^1 X^m Y^n \, dX \, dY = \frac{1}{(m+1)(n+1)}
$$

以及 $\tilde{Q}(X, Y)$ 的多项式展开，可以直接计算出FI的值。

## 5. 数值计算：质数域

### 5.1 有限域运算

所有运算在质数域 $\mathbb{Z}_p$ 中进行（$p$ 为大质数）：

- 多项式系数、有理函数系数、级数系数都是 $\mathbb{Z}_p$ 中的元素
- 避免浮点数误差，保证符号计算的精确性
- 模逆运算：$a^{-1} \equiv a^{p-2} \pmod{p}$（Fermat小定理）

### 5.2 有理数重建（可选）

通过中国剩余定理（CRT），可以从多个质数 $p_1, p_2, \ldots, p_k$ 下的结果重建有理数系数：

1. 对每个 $p_i$ 独立计算，得到 $c \bmod p_i$
2. 使用CRT重建整数 $c \in \mathbb{Z}$
3. 通过连分数等方法识别有理数 $c = n/d$

## 6. 完整流程总结

给定二圈Feynman积分，完整的计算流程为：

```
输入：
  - 顶层矩阵 R(X,Y)（重标度后）
  - 分支数 B = 3，传播子数 N
  - 有理函数 Q(X,Y)
  - 目标FBI (ν, Δ)
  - 级数度数 d
  - 质数 p

步骤0：Family初始化
  0.1 构造分支索引
  0.2 枚举所有sector（2^N个）
  0.3 对每个sector：
      - 构造矩阵 S
      - 计算 dim(Null(S))、C、z值
      - 分类Case（0/1/2/3）
  0.4 识别所有MFBI（Case 0角积分）

步骤1：FBI约化
  1.1 对目标FBI (ν, Δ)：
      - 根据Case类型应用约化策略
      - 递归约化到MFBI
      - 得到有理函数系数 r_k(X,Y)

步骤2：微分方程求解
  2.1 构造MFBI的微分方程系数矩阵 A_X, A_Y
  2.2 设置初值条件 f(0,0)
  2.3 对每个MFBI：
      - 使用递推公式逐阶计算系数
      - 得到级数展开 f_k(X,Y)

步骤3：FBI级数计算
  3.1 级数与有理函数乘法：r_k(X,Y) * f_k(X,Y)
  3.2 对所有k求和
  3.3 得到FBI级数展开 I(X,Y)

步骤4：Feynman积分计算
  4.1 计算 Q(X,Y) * I(X,Y) 的级数展开
  4.2 利用幂次积分公式进行数值积分
  4.3 得到FI的值

输出：
  - Feynman积分FI的数值（在Zp中）
```

这个流程在质数域 $\mathbb{Z}_p$ 中完全确定且高效，避免了数值误差，适合高精度符号计算。

## 7. 参考文献

[1] Li-Hong Huang, Rui-Jun Huang, and Yan-Qing Ma, "Tame multi-leg Feynman integrals beyond one loop," arXiv:2412.21053v1 [hep-ph], 2024.
