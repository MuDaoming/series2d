# FBI Sector Symmetry：问题定义与计算流程

## 1. 研究目标

本模块在固定 FBI family 中寻找不同 sector 之间由 Feynman 参数重标记产生的离散对称关系。

FBI 的分母二次型为

$$
F_s(y;X,Y)=\frac12 y^T R_s(X,Y)y,
$$

其中：

- $s$ 是 sector；
- $R_s$ 是 family 矩阵 $R$ 在 sector 活跃传播子上的主子矩阵；
- $y_\alpha$ 是 branch 内 Feynman 参数；
- 当前工程中的 $X,Y$ 是经过单纯形到单位正方形变换并局部平移后的坐标。

目标是寻找两个 sector $s,t$ 之间所有满足条件的变换 $(\sigma,\tau)$：

- $\sigma$ 是活跃传播子以及对应 $y$ 参数的置换；
- $\tau$ 是 branch 的置换；
- $\sigma$ 与 $\tau$ 保持传播子到 branch 的归属；
- 二次型在变换后完全相同。

由此得到 FBI 函数关系

$$
I_{s,\nu}^{\Delta}(X,Y)
=
I_{t,\sigma_*\nu}^{\Delta}\bigl(g_\tau(X,Y)\bigr),
$$

其中 $g_\tau$ 是 branch 置换在当前二维坐标上诱导的有理变换。

本模块只寻找参数置换对称性，不寻找 IBP、维度迁移或一般参数变换产生的线性关系。

## 2. 数学对象

### 2.1 Family、branch 与 sector

设 family 有 $N$ 条传播子和 $B$ 个 branch。传播子 $\alpha$ 的 branch 记为

$$
\operatorname{br}(\alpha)\in\{0,\ldots,B-1\}.
$$

sector 是二进制向量

$$
s=(s_1,\ldots,s_N),\qquad s_\alpha\in\{0,1\}.
$$

活跃传播子集合为

$$
A_s=\{\alpha\mid s_\alpha=1\}.
$$

本模块沿用 `expand` 的有效性条件：每个 branch 至少包含一个活跃传播子，

$$
\forall b,\quad
A_s\cap\operatorname{br}^{-1}(b)\ne\varnothing.
$$

sector 子矩阵定义为

$$
R_s(X,Y)=R(X,Y)[A_s,A_s].
$$

### 2.2 当前二维 branch 坐标

对三 branch family，原始 branch 参数满足

$$
X_1+X_2+X_3=1.
$$

`expand` 使用如下正方形参数化和平移：

$$
\begin{aligned}
X_1 &= X+a,\\
X_2 &= (1-X-a)(Y+b),\\
X_3 &= (1-X-a)(1-Y-b).
\end{aligned}
$$

记该映射为

$$
\Phi_{a,b}(X,Y)=(X_1,X_2,X_3).
$$

其逆映射为

$$
\Phi_{a,b}^{-1}(X_1,X_2,X_3)
=
\left(
X_1-a,\;
\frac{X_2}{1-X_1}-b
\right).
$$

### 2.3 Branch 置换

branch 置换表示为

$$
\tau=(\tau(0),\ldots,\tau(B-1)).
$$

本文采用以下固定约定：置换后的第 $b$ 个 branch 参数取自原来的第 $\tau(b)$ 个 branch，

$$
X'_b=X_{\tau(b)}.
$$

对 $B=3$，它在当前二维坐标中诱导出

$$
g_\tau
=
\Phi_{a,b}^{-1}\circ\tau\circ\Phi_{a,b}.
$$

即

$$
\begin{aligned}
X' &= X_{\tau(0)}-a,\\
Y' &= \frac{X_{\tau(1)}}{1-X_{\tau(0)}}-b.
\end{aligned}
$$

$g_\tau(X,Y)$ 一般是有理变换。实现必须对代换结果做精确有理化简，不能把 $\tau$ 简化为 $X,Y$ 的直接交换。

## 3. 相容的传播子置换

固定 $\tau$ 后，传播子置换 $\sigma:A_s\to A_t$ 必须满足

$$
\operatorname{br}(\sigma(\alpha))
=
\tau^{-1}\bigl(\operatorname{br}(\alpha)\bigr),
$$

等价地，若把规范表示中的 branch $b$ 看作接收位置，则该位置只能放入原 branch $\tau(b)$ 的传播子。

因此，只有 branch 活跃传播子数满足

$$
|A_s\cap B_{\tau(b)}|
=
|A_t\cap B_b|
$$

时，给定 $\tau$ 才可能产生 $s\to t$ 的关系。

对活跃指数的推送定义为

$$
(\sigma_*\nu)_{\sigma(\alpha)}=\nu_\alpha.
$$

不活跃传播子的指数保持为零。

## 4. FBI Sector 等价条件

设 $P_\sigma$ 是 $\sigma$ 在活跃传播子空间上的置换矩阵。sector $s,t$ 在 $(\sigma,\tau)$ 下等价，当且仅当

$$
\boxed{
P_\sigma^T R_s(X,Y)P_\sigma
=
R_t\bigl(g_\tau(X,Y)\bigr)
}
$$

并且 $\sigma$ 与 $\tau$ 满足第 3 节的 branch 相容条件。

证明来自 FBI 积分中的变量替换：

$$
y=P_\sigma y'.
$$

二次型变为

$$
\frac12y^TR_sy
=
\frac12y'^TP_\sigma^TR_sP_\sigma y'.
$$

传播子置换保持乘积测度和单项式权重，branch 相容条件保证每个

$$
\delta\left(1-\sum_{\alpha\in B_b}y_\alpha\right)
$$

被映射到另一个完整 branch 的约束。于是得到

$$
I_{s,\nu}^{\Delta}(X,Y)
=
I_{t,\sigma_*\nu}^{\Delta}\bigl(g_\tau(X,Y)\bigr).
$$

## 5. Pak 型规范化

### 5.1 规范对象

对 sector $s$ 和固定 $\tau$，先对矩阵元素做坐标代换：

$$
R_s^\tau(X,Y)=R_s\bigl(g_\tau(X,Y)\bigr).
$$

随后在所有与 $\tau$ 相容的传播子排列下，对称地重排行列。

对一个排列后的 $n\times n$ 对称矩阵，采用下三角序列

$$
\mathcal S(R)
=
(R_{11},R_{21},R_{22},R_{31},R_{32},R_{33},\ldots,R_{nn})
$$

作为比较对象。每个矩阵元素先经过精确 `normal` 和展开，再序列化为确定字符串。

sector 的 canonical key 定义为所有允许 $(\sigma,\tau)$ 表示中字典序最小的序列：

$$
\operatorname{Can}(s)
=
\min_{\tau\in S_B}
\min_{\sigma\text{ compatible with }\tau}
\mathcal S\left(P_\sigma^TR_s^\tau P_\sigma\right).
$$

### 5.2 逐位置前缀剪枝

直接枚举 $\sigma$ 需要检查大量传播子排列。Pak 型算法逐个确定规范矩阵中的传播子位置。

假设前 $k-1$ 个位置已经确定。为第 $k$ 个位置尝试所有 branch 相容且未使用的传播子。加入该传播子后，规范序列新增前缀块

$$
(R_{k1},\ldots,R_{k,k-1},R_{kk}).
$$

算法只保留新增完整前缀字典序最小的候选；并列候选全部保留并进入下一位置。

被删除候选已经在最终序列的最早不同位置上严格大于保留候选，后续元素无法改变这个顺序。因此该剪枝与枚举所有相容 $\sigma$ 后取全局最小值完全等价。

### 5.3 显式枚举 $\tau$

branch 数通常远小于传播子数。第一版显式枚举全部

$$
\tau\in S_B
$$

并对每个 $\tau$ 独立运行第 5.2 节的 Pak 规范化。最后在所有 $\tau$ 的结果中选择全局最小 canonical key。

## 6. Orbit 与代表 Sector

若

$$
\operatorname{Can}(s)=\operatorname{Can}(t),
$$

则 $s,t$ 属于同一个 symmetry orbit。

规范化同时保存：

$$
T_s=(\sigma_s,\tau_s):s\to\operatorname{Can}(s),
$$

以及

$$
T_t=(\sigma_t,\tau_t):t\to\operatorname{Can}(t).
$$

由二者复合可构造 $s\to t$ 的候选关系。由于 branch 坐标变换的方向和传播子索引约定容易产生逆置换错误，实现不依赖形式推断作为最终判据，而是根据两边到共同规范矩阵的有序传播子列表直接恢复映射，并执行第 4 节的精确矩阵验证。

每个 orbit 选择一个确定代表 sector。第一版采用 sector 二进制索引最大的 sector 作为代表，即按 `expand::Family::idxFromSecvec` 的整数值选择最大者。

最终只需保存每个非代表 sector 到代表的关系；orbit 内任意两 sector 的关系可由这些映射复合得到。

## 7. 精确验证

canonical key 相同只用于发现候选。每个输出关系必须重新验证：

1. source 和 target 活跃传播子数相同；
2. $\sigma$ 是 source 活跃集合到 target 活跃集合的双射；
3. $\sigma$ 与 $\tau$ 的 branch 归属相容；
4. 对每个矩阵元素，精确检查

   $$
   \operatorname{normal}\left(
   [P_\sigma^TR_sP_\sigma]_{ij}
   -
   [R_t(g_\tau)]_{ij}
   \right)=0.
   $$

验证失败的候选不得进入最终结果。

## 8. 完整计算流程

输入：

- 当前 `expand` 使用的 $(B+N)\times(B+N)$ 矩阵 $S(X,Y)$；
- $B,N$；
- 局部平移参数 $a,b$。

步骤：

1. 从 $S$ 的 branch-incidence block 提取每个传播子的 branch。
2. 从 $S$ 的右下角提取 $R(X,Y)$。
3. 枚举所有覆盖全部 branch 的有效 sector。
4. 对每个 sector $s$：
   1. 提取 $R_s(X,Y)$；
   2. 枚举所有 branch 置换 $\tau$；
   3. 计算 $g_\tau(X,Y)$；
   4. 构造 $R_s(g_\tau(X,Y))$；
   5. 使用 Pak 前缀剪枝寻找与 $\tau$ 相容的最小传播子排列；
   6. 在全部 $\tau$ 结果中选择 sector canonical key；
   7. 保存所有达到同一最小 key 的 canonicalizing 变换。
5. 按 canonical key 将 sector 分成 symmetry orbit。
6. 为每个 orbit 选择代表 sector。
7. 从共同 canonical ordering 恢复每个 sector 到代表的 $(\sigma,\tau)$。
8. 对每个关系执行第 7 节的精确验证。
9. 输出 orbit、代表 sector、$\sigma$、$\tau$ 和 $g_\tau(X,Y)$。

## 9. 适用范围与限制

第一版：

- 直接使用当前 $R(X,Y)$，不恢复 $R(X_1,X_2,X_3)$；
- 支持 `expand` 当前三 branch 坐标映射；
- 使用精确 GiNaC 有理运算；
- 寻找全部参数置换 symmetry；
- 不处理 IBP、维度迁移、一般线性参数变换或仅在特殊 $X,Y$ 点成立的偶然关系。

## 10. 参考资料

1. A. Pak, “The Toolbox of modern multi-loop calculations: novel analytic and semi-analytic techniques”, arXiv:1111.0868.
2. `expand/docs/2412.21053v1.pdf`，FBI 定义及 branch 参数结构。
3. `expand/docs/problem_solution.md`，当前项目的 sector、$S$ 矩阵和二维坐标约定。

## 11. Phase A 审查点

进入架构和实现前，领域专家需要确认：

- 第 2.3 节的 $\tau$ 方向约定；
- 第 3 节的 $\sigma$ 与 $\tau$ 相容关系；
- 第 4 节的 FBI 等价公式；
- 第 6 节代表 sector 的选择规则。

本次实现按上述约定执行。
