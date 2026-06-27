# BL Sector Reduction：问题定义与工作流程

## 1. 研究目标

`bl_sector_reduce` 的目标是在有限域 $\mathbb{F}_p$ 上，利用一维 $\delta$ 级数直接重构符号积分约化关系。

给定一个待约化对象 $A$，本模块希望输出关于 $\delta$ 的有理函数系数表达：

$$
A(\delta)
=
\sum_{\mu} R_{A,\mu}(\delta) M_\mu(\delta),
\qquad
R_{A,\mu}(\delta)\in\mathbb{F}_p(\delta),
$$

或等价的公共分母形式：

$$
P_A(\delta) A(\delta)
=
\sum_{\mu} P_{A,\mu}(\delta) M_\mu(\delta),
\qquad
P_A,P_{A,\mu}\in\mathbb{F}_p[\delta].
$$

这里 $M_\mu$ 是已知主积分对象。主输出是符号多项式/有理函数约化；若后续需要某个数值点 $\delta_0$ 的约化，只需对这些多项式系数赋值。

本模块采用 sector-by-sector 的 Beckermann-Labahn（BL）/ approximant basis 重构路线。与全局关系搜索不同，它不一次性寻找所有积分之间的任意多项式关系，而是在每个 FBI sector 上直接重构当前 sector 的贡献到该 sector 的主积分。

## 2. 基本输入与对象

### 2.1 有限域与截断级数

所有计算在素域 $\mathbb{F}_p$ 中进行。输入给出若干对象在不同 FBI sector 边界条件下的一维截断级数：

$$
A^{(s)}(\delta)
=
\sum_{n=0}^{D} a_{A,s,n}\delta^n
\quad \bmod \delta^{D+1},
$$

其中：

- $A$ 是一个原始待约化对象或主积分对象；
- $s$ 是 FBI sector；
- $D$ 是当前展开阶数；
- $a_{A,s,n}\in\mathbb{F}_p$。

外部生成这些级数的方式不是本模块的一部分。当前工程中，这些级数通常来自 `expand`。对象名可以沿用已有的积分标签格式，例如 `FI{nu}`、`BFI[XU]{nu}` 或 `BBFI[XU,YD]{nu}`，但本模块的数学问题只要求它们是可索引的一维级数对象。

### 2.2 FBI sector 与边界条件

FBI sector 由传播子支持集定义。给定一个传播子指数向量 $\nu=(\nu_1,\ldots,\nu_N)$，其 sector 是

$$
\operatorname{sec}(\nu)
=
(s_1,\ldots,s_N),
\qquad
s_i=
\begin{cases}
1, & \nu_i>0,\\
0, & \nu_i=0.
\end{cases}
$$

本模块采用以下约定：

- 一个 FBI sector 对应一个边界条件；
- sector 按包含关系形成一棵处理树；
- 根节点是全局 top sector；
- 若 $t$ 是 $s$ 的父节点或祖先节点，则 $s$ 是 $t$ 的 subsector。

处理顺序由这棵 sector 树硬确定：从根节点向叶节点处理。sector 划分和处理顺序不是用户输入的自由参数。

### 2.3 原始对象、主对象与 sector contribution

输入对象分为两类：

1. 原始待约化对象 $A$；
2. 每个 sector 的主对象集合 $\mathcal{M}(s)$。

主对象集合 $\mathcal{M}(s)$ 由外部步骤给出，本模块不负责寻找 master。工程上可以读取类似 `search/runs/.../masters` 的文件，但具体格式由 Phase B 设计。

核心被约化的量不是完整的 $A$，而是 $A$ 在某个 sector $s$ 上的新贡献，记为

$$
C_s(A).
$$

这个记号的意义是：在 sector $s$ 的边界条件下，先从 $A^{(s)}$ 中减去所有祖先 sector 已经确定的贡献，剩下的部分就是当前 sector 需要解释的贡献。

因此 $A$ 可以仍然用原始标签记录，例如 `FI{nu}`；而实际送入 BL 的对象是带 sector 语义的

$$
C_s(A).
$$

这样避免把对象定义成任意线性组合，同时又能表达 subsector 上的已剥离贡献。

## 3. Sector tree 上的贡献剥离

### 3.1 根节点

设 $t$ 是 top sector。由于没有祖先节点，根节点贡献为

$$
C_t(A)^{(t)} = A^{(t)}.
$$

在 top sector 上，目标是寻找多项式

$$
P_{A,t},\quad P_{A,t,\mu}\in\mathbb{F}_p[\delta],
\qquad \mu\in\mathcal{M}(t),
$$

使得

$$
P_{A,t}(\delta) C_t(A)^{(t)}(\delta)
=
\sum_{\mu\in\mathcal{M}(t)}
P_{A,t,\mu}(\delta) M_\mu^{(t)}(\delta).
$$

### 3.2 子节点

设 $s$ 是某个非根 sector，祖先链为

$$
t=s_0 \to s_1 \to \cdots \to s_k=s.
$$

假设所有祖先 $s_0,\ldots,s_{k-1}$ 的贡献都已经被重构。把这些已知贡献在 sector $s$ 的边界条件下求值并从 $A^{(s)}$ 中减去，定义当前 sector contribution：

$$
C_s(A)^{(s)}
=
A^{(s)}
-
\sum_{i=0}^{k-1}
\operatorname{EvalAtSector}
\left(
\operatorname{Red}_{s_i}(A), s
\right).
$$

其中 $\operatorname{Red}_{s_i}(A)$ 表示祖先 sector $s_i$ 已经得到的符号约化项，$\operatorname{EvalAtSector}(\cdot,s)$ 表示把其中出现的对象都替换为它们在 sector $s$ 边界条件下的级数。

然后在 sector $s$ 上重构

$$
P_{A,s}(\delta) C_s(A)^{(s)}(\delta)
=
\sum_{\mu\in\mathcal{M}(s)}
P_{A,s,\mu}(\delta) M_\mu^{(s)}(\delta).
$$

### 3.3 更深层节点

上述定义递归适用于任意深度。若处理到 $s_{11}$，则其 contribution 是

$$
C_{s_{11}}(A)^{(s_{11})}
=
A^{(s_{11})}
-
\operatorname{EvalAtSector}(\operatorname{Red}_{t}(A),s_{11})
-
\operatorname{EvalAtSector}(\operatorname{Red}_{s_1}(A),s_{11})
-\cdots .
$$

也就是说，在每个 sector 边界条件下，当前对象和所有祖先 sector 的贡献都可能非零；本模块逐层剥离祖先贡献，然后只重构当前 sector 的新贡献。

本文档使用“sector contribution”或“subsector contribution”描述这个量，不使用“残差”作为核心术语。

## 4. BL 重构问题

### 4.1 固定 sector 的问题形式

固定一个 sector $s$ 和一个原始对象 $A$。设当前 sector 有 $r$ 个 master：

$$
\mathcal{M}(s)=\{M_1,\ldots,M_r\}.
$$

在该 sector 上，要寻找多项式向量

$$
P(\delta)
=
\left(P_0(\delta),P_1(\delta),\ldots,P_r(\delta)\right)^T
$$

满足

$$
P_0(\delta) C_s(A)^{(s)}(\delta)
-
\sum_{j=1}^{r}P_j(\delta)M_j^{(s)}(\delta)
=0
\quad \bmod \delta^{K}.
$$

若 $P_0\ne 0$，则得到 sector $s$ 的符号贡献约化：

$$
C_s(A)
=
\sum_{j=1}^{r}
\frac{P_j(\delta)}{P_0(\delta)}M_j.
$$

如果 $P_0=0$，则该关系只涉及 master 之间的依赖，不能用于约化 $C_s(A)$。

### 4.2 Approximant 表述

定义一维向量级数

$$
F_{A,s}(\delta)
=
\left(
C_s(A)^{(s)}(\delta),
-M_1^{(s)}(\delta),
\ldots,
-M_r^{(s)}(\delta)
\right).
$$

则重构条件为

$$
F_{A,s}(\delta)P(\delta)
=0
\quad \bmod \delta^K.
$$

这就是一个向量 Hermite-Pade / approximant basis 问题。BL 算法用于在给定 degree bound 和工作阶数 $K$ 下寻找满足条件的低复杂度多项式向量。

## 5. Degree 与使用阶数

### 5.1 统一 degree bound

第一版采用统一 degree bound。给定非负整数 $m$，要求

$$
\deg P_i \le m,
\qquad i=0,1,\ldots,r.
$$

未知系数数目为

$$
U(r,m)=(r+1)(m+1).
$$

本文档不使用 $U-1$ 作为工作阶数。虽然齐次关系整体缩放不重要，但使用 $U-1$ 个方程会导致任意输入都必然存在非零候选，不能作为有效搜索判据。因此最小工作阶数取为

$$
K_{\min}(r,m)=U(r,m).
$$

### 5.2 保险阶数与认证阶数

为了降低有限阶退化和伪关系风险，实际 BL 工作阶数取为

$$
K_{\mathrm{work}}(r,m)
=
U(r,m)+K_{\mathrm{safety}}+K_{\mathrm{cert}}.
$$

其中：

- $K_{\mathrm{safety}}$ 是保险阶数；
- $K_{\mathrm{cert}}$ 是认证阶数；
- 默认值为

  $$
  K_{\mathrm{safety}}=10,\qquad K_{\mathrm{cert}}=10.
  $$

这两个量是配置参数，不应在实现中写成 magic number。

给定展开阶数 $D$，可用系数数量是 $D+1$。degree $m$ 可被当前数据支持，当且仅当

$$
D+1 \ge K_{\mathrm{work}}(r,m).
$$

等价地，当前展开阶数支持的最大统一 degree 是满足

$$
(r+1)(m+1)+K_{\mathrm{safety}}+K_{\mathrm{cert}}
\le D+1
$$

的最大整数 $m$。

### 5.3 给定 degree 的接受判据

给定 $m$，使用阶数

$$
0,1,\ldots,K_{\mathrm{work}}(r,m)-1
$$

直接参与 BL 计算。若 BL 返回多项式向量

$$
P=(P_0,\ldots,P_r)^T,
$$

则接受该 sector contribution 的约化，当且仅当：

1. $P$ 不全为零；
2. $\deg P_i\le m$ 对所有 $i$ 成立；
3. $P_0\ne 0$；
4. $F_{A,s}(\delta)P(\delta)=0 \bmod \delta^{K_{\mathrm{work}}(r,m)}$。

因为 $K_{\mathrm{safety}}$ 和 $K_{\mathrm{cert}}$ 已经并入 BL 工作阶数，本模块默认不再设置独立的额外检查流程。

### 5.4 degree 搜索

用户输入最高允许 degree：

$$
m_{\mathrm{user}}.
$$

算法先根据 $D,r,K_{\mathrm{safety}},K_{\mathrm{cert}}$ 计算当前数据支持的最高 degree：

$$
m_{\mathrm{supported}}
=
\max\left\{
m\ge0:
(r+1)(m+1)+K_{\mathrm{safety}}+K_{\mathrm{cert}}
\le D+1
\right\}.
$$

实际搜索上限为

$$
m_{\max}=\min(m_{\mathrm{user}},m_{\mathrm{supported}}).
$$

若 $m_{\mathrm{supported}}<0$，则当前展开阶数不足以搜索该 sector 的任何 $m\ge0$ 关系。

degree 搜索使用指数序列：

$$
0,1,2,4,8,\ldots
$$

直到不超过 $m_{\max}$。若指数序列跳过 $m_{\max}$，最后补充一次 $m_{\max}$。

找到第一个通过第 5.3 节判据的 degree 后立即接受，不要求寻找最小 degree。

若搜索到 $m_{\max}$ 仍未找到可接受关系，则报告：

- 当前对象和 sector 未约化；
- 已搜索到 degree $\le m_{\max}$；
- 当前展开阶数 $D$ 在配置的 $K_{\mathrm{safety}},K_{\mathrm{cert}}$ 下最多支持 $m_{\mathrm{supported}}$；
- 若 $m_{\max}=m_{\mathrm{supported}}$，需要更高展开阶数或降低安全/认证余量；
- 若 $m_{\max}=m_{\mathrm{user}}<m_{\mathrm{supported}}$，需要提高用户 degree 上限。

## 6. 完整 sector-wise 算法

### 6.1 输入

`bl_sector_reduce` 的数学输入是：

1. 有限域特征 $p$；
2. 展开阶数 $D$；
3. 原始待约化对象集合 $\mathcal{A}$；
4. 每个对象和每个 sector 边界条件下的一维级数；
5. 每个 sector 的 master 集合 $\mathcal{M}(s)$；
6. 用户最高 degree $m_{\mathrm{user}}$；
7. 保险阶数 $K_{\mathrm{safety}}$，默认 10；
8. 认证阶数 $K_{\mathrm{cert}}$，默认 10。

sector tree、sector 划分和处理顺序由传播子支持集关系在代码中确定，不作为自由输入。

### 6.2 单个对象的处理

对每个原始对象 $A\in\mathcal{A}$，沿 sector tree 从根到叶处理。

在 sector $s$：

1. 读取 $A^{(s)}$。
2. 把所有祖先 sector 已确定的贡献在 sector $s$ 的边界条件下求值。
3. 构造当前 contribution：

   $$
   C_s(A)^{(s)}
   =
   A^{(s)}
   -
   \sum_{a\in\operatorname{Ancestors}(s)}
   \operatorname{EvalAtSector}(\operatorname{Red}_a(A),s).
   $$

4. 若 $C_s(A)^{(s)}$ 在已知阶数内为零，则该 sector 对 $A$ 没有新贡献。
5. 若 $A$ 在 sector $s$ 被声明为 master，则记录它为自由 master，不对其当前贡献做约化。
6. 否则，对 $C_s(A)$ 和 $\mathcal{M}(s)$ 执行第 4-5 节的 BL 重构。
7. 若成功，记录 $\operatorname{Red}_s(A)$。
8. 若失败，报告该 sector contribution 未约化，并继续或停止由 Phase B 的错误策略决定。

### 6.3 全局输出

对每个原始对象 $A$，最终符号约化由所有 sector contribution 相加得到：

$$
A
=
\sum_s C_s(A)
=
\sum_s
\sum_{\mu\in\mathcal{M}(s)}
R_{A,s,\mu}(\delta) M_\mu.
$$

其中

$$
R_{A,s,\mu}(\delta)
=
\frac{P_{A,s,\mu}(\delta)}{P_{A,s}(\delta)}.
$$

输出可以保存为有理函数形式，也可以保存为公共分母多项式形式。Phase B 需要指定具体文件格式，但数学内容必须保留每个 sector contribution 的分母和分子多项式。

## 7. 与现有模块的接口关系

### 7.1 与级数生成模块的关系

本模块不生成 $A^{(s)}(\delta)$。它只要求输入能够提供：

```text
object label
sector label
coefficient list {a_0,...,a_D}
```

当前工程中，这些数据通常由 `expand` 生成。`FI/BFI/BBFI` 等标签可以作为对象 label 的一种具体格式，但不属于本模块数学定义的必要组成部分。

### 7.2 与 master 生成步骤的关系

本模块不寻找每个 sector 的 master。它只读取：

```text
sector -> list of master object labels
```

这些 master 可以由当前 maximal-cut search 流程得到，也可以由其他工具生成。

### 7.3 与 `search` 的关系

`search` 的关系搜索路线是：

1. 全局搜索多项式关系；
2. 对关系求值或再做线性代数；
3. 得到 reduction。

`bl_sector_reduce` 的路线是：

1. 在每个 sector 上构造当前 contribution；
2. 用 BL 直接重构到该 sector master；
3. 沿 sector tree 逐层剥离并合成符号约化。

`search` 可以继续作为小规模验证工具或 master 生成工具；`bl_sector_reduce` 负责新的 sector-wise BL 符号约化路线。

## 8. 完整计算流程

```text
输入：
  p
  D
  object series for each sector boundary condition
  master list for each sector
  m_user
  K_safety
  K_cert

准备：
  1. 解析对象标签。
  2. 根据传播子支持集构造 sector tree。
  3. 按根到叶得到 sector 处理顺序。
  4. 为每个 sector 读取 master 列表。

对每个对象 A：
  for sector s in root-to-leaf order:
      build C_s(A)^s by subtracting ancestor contributions
      if C_s(A)^s is zero:
          continue
      if A is a master in sector s:
          record free master contribution
          continue
      r = number of masters in M(s)
      compute m_supported from D, r, K_safety, K_cert
      m_max = min(m_user, m_supported)
      for m in exponential schedule up to m_max:
          K_work = (r+1)(m+1) + K_safety + K_cert
          run BL with orders 0..K_work-1
          if accepted P with P0 != 0:
              record Red_s(A)
              break
      if no accepted P:
          report unresolved sector contribution

输出：
  symbolic sector contributions Red_s(A)
  global symbolic reductions A = sum_s Red_s(A)
  unresolved contributions, if any
```

## 9. 失败条件

一个 sector contribution 可能失败，主要原因包括：

1. 当前展开阶数 $D$ 不足以支持 $m=0$；
2. 用户 degree 上限 $m_{\mathrm{user}}$ 太低；
3. $K_{\mathrm{safety}}$ 或 $K_{\mathrm{cert}}$ 设置过大，导致可搜索 degree 太低；
4. 当前 sector 的 master 集合不完整；
5. 输入级数或 sector tree 归属不一致；
6. 当前有限域素数导致退化，需要换素数验证。

失败时必须报告 sector、对象、master 数 $r$、已搜索 degree 上限、$D$、$K_{\mathrm{safety}}$ 和 $K_{\mathrm{cert}}$。

## 10. 后续 Phase B 需要设计的内容

Phase A 已经固定数学问题和算法逻辑。Phase B 需要把以下内容转成代码结构：

1. object label、sector label、master list 的具体文件格式；
2. sector tree 的构造和遍历接口；
3. sector contribution 的缓存表示；
4. 祖先贡献在当前 sector 边界条件下的求值机制；
5. BL / approximant basis 调用接口；
6. 多项式向量、公共分母和有理函数输出格式；
7. 失败报告和可恢复策略。

## 11. 参考资料

1. `expand/docs/problem_solution.md`：一维 $\delta$ 级数生成背景。
2. `search/docs/problem_solution.md`：全局多项式关系搜索原型。
3. `search/docs/ideal_absorbing_poly_search.md`：关系搜索中的 approximant 与吸收思想。
4. Beckermann-Labahn 类型 approximant basis 算法。

## 12. Phase A 审查点

本文档将 `bl_sector_reduce` 定义为一个自洽的 sector-wise 符号约化问题：

- 输入是一维 $\delta$ 级数、sector tree 和每个 sector 的 master 列表；
- 核心对象是 $C_s(A)$，即原始对象 $A$ 在 sector $s$ 上的新贡献；
- 每个 sector 只使用自己的边界条件；
- degree 搜索采用统一 bound $m$ 和指数 schedule；
- 给定 $m$ 时 BL 工作阶数为

  $$
  (r+1)(m+1)+K_{\mathrm{safety}}+K_{\mathrm{cert}},
  $$

  默认 $K_{\mathrm{safety}}=10$、$K_{\mathrm{cert}}=10$；
- 输出是符号多项式/有理函数约化。

进入 Phase B 前，需要由领域专家确认上述数学定义与 sector contribution 剥离逻辑正确。
