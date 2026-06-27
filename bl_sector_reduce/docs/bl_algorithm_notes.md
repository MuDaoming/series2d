# BL / Approximant Basis 处理固定 sector 约化问题的说明

本文档只讨论固定一个 sector、固定一个待约化对象时的问题：

$$
P_0(\delta) I(\delta)
=
\sum_{i=1}^{r} P_i(\delta) M_i(\delta).
$$

这里 $I$ 可以理解为当前 sector contribution 的一维级数，$M_i$ 是该 sector 的 master 级数。本文档不讨论 sector 间 contribution 剥离。

## 1. 问题形式

输入是一组截断级数：

$$
I(\delta),\quad M_1(\delta),\ldots,M_r(\delta)
\quad \bmod \delta^{D+1}.
$$

目标是找多项式：

$$
P_0(\delta),P_1(\delta),\ldots,P_r(\delta)\in\mathbb{F}_p[\delta],
$$

满足

$$
P_0(\delta) I(\delta)
-
\sum_{i=1}^{r}P_i(\delta)M_i(\delta)
=0
\quad \bmod \delta^K.
$$

若 $P_0\ne 0$，则可以写成符号约化：

$$
I(\delta)
=
\sum_{i=1}^{r}\frac{P_i(\delta)}{P_0(\delta)}M_i(\delta).
$$

如果 $P_0=0$，得到的是 master 之间的关系，不能直接约化 $I$。

## 2. 转成 approximant 问题

定义一个长度为 $r+1$ 的向量级数：

$$
F(\delta)=
\left(
I(\delta),
-M_1(\delta),
\ldots,
-M_r(\delta)
\right).
$$

定义多项式列向量：

$$
P(\delta)=
\left(
P_0(\delta),
P_1(\delta),
\ldots,
P_r(\delta)
\right)^T.
$$

那么原问题就是：

$$
F(\delta)P(\delta)=0\quad \bmod \delta^K.
$$

这就是一个 order-$K$ 的向量 approximant 问题。所有满足上式的多项式向量 $P$ 构成一个 $\mathbb{F}_p[\delta]$-module。

BL 算法的任务不是只找一个关系，而是构造这个 approximant module 的一组好基，例如 minimal approximant basis。得到这个 basis 后，再从 basis 中选出满足 degree bound 且 $P_0\ne 0$ 的向量。

## 3. Degree bound 与工作阶数

当前采用统一 degree bound：

$$
\deg P_i\le m,\qquad i=0,\ldots,r.
$$

未知数数量为：

$$
U=(r+1)(m+1).
$$

为了让关系可判定，工作阶数至少要覆盖未知数数量。当前设计使用：

$$
K = U + K_{\mathrm{safety}} + K_{\mathrm{cert}}.
$$

其中：

- $K_{\mathrm{safety}}$ 用来避免退化情况下刚好阶数不够；
- $K_{\mathrm{cert}}$ 作为额外验证阶数；
- 实现上没有单独再跑一个后处理 check，而是把这两部分都并入 BL 使用的阶数。

给定展开阶数 $D$，必须满足：

$$
K\le D+1.
$$

因此在 $r,D,K_{\mathrm{safety}},K_{\mathrm{cert}}$ 固定时，支持的最大统一 degree 是：

$$
m_{\max}
=
\left\lfloor
\frac{D+1-K_{\mathrm{safety}}-K_{\mathrm{cert}}}{r+1}
\right\rfloor-1.
$$

## 4. 标准 BL / order basis 算法怎么处理

### 4.1 要计算的对象

固定

$$
n=r+1.
$$

我们要计算的是 order-$K$ approximant module：

$$
\mathcal{A}_K(F)
=
\{P\in\mathbb{F}_p[\delta]^n: F(\delta)P(\delta)=0\bmod \delta^K\}.
$$

BL 算法通常返回一个 $n\times n$ 多项式矩阵 $B(\delta)$，它的列向量构成 $\mathcal{A}_K(F)$ 的一组基：

$$
F(\delta)B(\delta)=0\bmod\delta^K.
$$

也就是说，每一列 $B_j$ 都是一个关系：

$$
F B_j=0\bmod\delta^K.
$$

所有 order-$K$ 关系都可以写成

$$
P(\delta)=B(\delta)Q(\delta),
\qquad Q\in\mathbb{F}_p[\delta]^n.
$$

BL 的关键不是只让 $B$ 满足 $FB=0$，还要让 $B$ 在某个 shift 下 minimal。minimality 保证低 degree 关系可以从基中系统地读出，而不是藏在很复杂的列组合里。

### 4.2 Shifted degree

给定 shift

$$
\mathbf{s}=(s_0,\ldots,s_{n-1}),
$$

一个多项式向量

$$
P=(P_0,\ldots,P_{n-1})^T
$$

的 shifted degree 定义为

$$
\deg_{\mathbf{s}} P
=
\max_i(\deg P_i+s_i).
$$

当前统一 degree bound

$$
\deg P_i\le m
$$

最简单可以先用零 shift：

$$
\mathbf{s}=(0,\ldots,0),
$$

然后在输出候选时检查每个坐标 degree 是否都不超过 $m$。更一般地，也可以把不同坐标的 degree 约束编码进 shift，但当前还没有走到这一步。

### 4.3 BL 的分治递归

BL / order basis 的核心递归如下。这里写成列基形式。

输入：

```text
F(delta): 1 x n 的截断级数向量
K:        order
s:        shift
```

输出：

```text
B(delta): n x n 多项式矩阵
          F B = 0 mod delta^K
          B 是 s-minimal order basis
```

递归步骤：

1. 若 $K=1$，处理 base case。
2. 否则取

   $$
   K_1=\lfloor K/2\rfloor,\qquad K_2=K-K_1.
   $$

3. 先计算前半阶的 order basis：

   $$
   B_1=\operatorname{BL}(F\bmod\delta^{K_1}, K_1,\mathbf{s}).
   $$

   它满足：

   $$
   F B_1=0\bmod\delta^{K_1}.
   $$

4. 计算 residual：

   $$
   G(\delta)
   =
   \delta^{-K_1}F(\delta)B_1(\delta)
   \bmod\delta^{K_2}.
   $$

   因为 $F B_1$ 已经在前 $K_1$ 阶为零，所以可以除掉 $\delta^{K_1}$。$G$ 表示剩下后半阶还没有消掉的误差。

5. 更新 shift。设 $B_1$ 的各列 shifted degree 为

   $$
   \mathbf{s}'=\operatorname{cdeg}_{\mathbf{s}}(B_1).
   $$

   这里 $\operatorname{cdeg}_{\mathbf{s}}$ 表示每一列的 shifted column degree。这个更新很重要：后半段不是重新在原 shift 下 minimal，而是在前半段已经得到的基的 degree 基础上继续 minimal。

6. 对 residual 继续计算后半阶 order basis：

   $$
   B_2=\operatorname{BL}(G,K_2,\mathbf{s}').
   $$

   它满足：

   $$
   G B_2=0\bmod\delta^{K_2}.
   $$

7. 返回：

   $$
   B=B_1B_2.
   $$

验证：

$$
F B
=
F B_1B_2
=
\delta^{K_1}G B_2
=0\bmod\delta^{K_1+K_2}
=0\bmod\delta^K.
$$

这就是 BL 分治的核心。复杂度优势来自：不是逐阶做 $K$ 次大矩阵更新，而是把 order 分成两半，递归计算两个较小 order basis，并用多项式矩阵乘法/乘级数处理 residual。

### 4.4 Base case 怎么处理

当 $K=1$ 时，只需要消掉常数项：

$$
F_0 P_0 = 0,
$$

其中 $F_0=F(0)$ 是一个 $1\times n$ 常数向量，$P_0$ 是常数列向量。

若 $F_0=0$，则任意常数向量都是 order-1 approximant，可以返回单位矩阵：

$$
B=I_n.
$$

若 $F_0\ne0$，需要构造一个完整 rank 的 polynomial basis。对一行 $F_0$，可以选一个 pivot 位置 $q$，满足 $F_{0,q}\ne0$。则：

- 对每个 $j\ne q$，常数向量

  $$
  e_j-\frac{F_{0,j}}{F_{0,q}}e_q
  $$

  都满足 $F_0P=0$。

- 还需要一个方向补齐 rank。可以用

  $$
  \delta e_q
  $$

  因为它的常数项为零，所以也满足 order-1 条件。

这些列组成一个 order-1 basis。实际 BL 会根据 shift 选择 pivot，使返回的 basis 满足 shifted minimality。

### 4.5 从 basis 里读出约化关系

BL 返回 $B$ 后，候选关系来自 $B$ 的列。对每一列 $P$，检查：

$$
F P=0\bmod\delta^K,
$$

$$
P_0\ne0,
$$

$$
\deg P_i\le m,\quad i=0,\ldots,r.
$$

如果有这样的列，就得到：

$$
P_0 I=\sum_{i=1}^r P_iM_i.
$$

若使用的是正确的 shifted minimal basis，并且 degree bound/shift 设计匹配，那么存在的低 degree 关系应当能通过 basis 的低 shifted-degree 列被读出。

## 5. 当前算法实际怎么处理

当前实现位于：

- `include/approximant_basis.hpp`
- `src/approximant_basis.tpp`

当前输入是：

```cpp
ApproximantRequest<T> {
    target,      // I 的系数
    masters,     // M_i 的系数
    maxDegree,   // m
    workOrder    // K
}
```

代码内部先构造：

```text
f = (target, -masters[0], ..., -masters[r-1])
```

然后初始化一个单位矩阵形式的多项式 basis。代码把每个 basis vector 存成一行；数学上可把它理解成要找的多项式向量 $P$：

```text
basis = identity matrix
```

接着从 $k=0$ 到 $K-1$ 逐阶处理 discrepancy。

对第 $i$ 个 basis vector $b_i$，定义：

```text
disc_i(k) = coefficient of delta^k in F * b_i
```

用公式写就是：

$$
d_i(k)
=
[\delta^k]\left(F(\delta)b_i(\delta)\right).
$$

当前算法维护一个直观目标：处理完 $k-1$ 阶后，basis vectors 尽量满足

$$
F b_i=0\bmod\delta^k.
$$

第 $k$ 阶步骤如下：

1. 计算所有 $d_i(k)$。
2. 若所有 $d_i(k)=0$，则第 $k$ 阶已经满足，直接进入 $k+1$。
3. 否则在所有 $d_i(k)\ne0$ 的行中，选择普通 row degree 最小的一行作为 pivot，记为 $p$。
4. 对每个非 pivot 行 $i$，若 $d_i(k)\ne0$，执行：

   $$
   b_i \leftarrow b_i-\frac{d_i(k)}{d_p(k)}b_p.
   $$

   这样第 $i$ 行在第 $k$ 阶的 discrepancy 被消掉，因为：

   $$
   d_i(k)-\frac{d_i(k)}{d_p(k)}d_p(k)=0.
   $$

5. pivot 行乘以 $\delta$：

   $$
   b_p\leftarrow \delta b_p.
   $$

   乘以 $\delta$ 后，$F b_p$ 的低阶项整体后移，因此它在当前已处理阶不会破坏已满足的条件。

循环到 $K-1$ 后，代码逐行检查：

```text
row[0] != 0
degree(row_j) <= m
F * row = 0 mod delta^K
```

找到第一条满足条件的行就返回这条关系。

## 6. 当前算法和 BL 的具体差别

当前算法不是上面第 4 节描述的完整 Beckermann-Labahn 分治 minimal approximant basis。具体差别是：

1. **当前算法是逐阶 discrepancy 更新，不是 BL 分治递归**

   当前复杂度随 $K$ 逐阶循环；没有做

   $$
   K\to K_1,K_2,
   \qquad
   G=\delta^{-K_1}FB_1
   $$

   这种 residual 分治。

2. **当前算法没有维护 shifted minimality**

   标准 BL 每一步都带 shift，并且第二半递归用

   $$
   \mathbf{s}'=\operatorname{cdeg}_{\mathbf{s}}(B_1)
   $$

   更新 shift。当前代码只用普通 row degree 选 pivot，没有 shifted degree，也没有 Popov/minimal form 保证。

3. **当前算法的 basis 形态不保证适合读出所有低 degree 关系**

   标准 BL 返回 minimal basis 后，可以根据 shifted degree 读候选。当前算法只是得到一些满足 order 条件的 vectors，然后直接逐行检查。它可能在很多例子上能找到关系，但没有完整 BL 的“低 degree 关系一定以可读形式出现”的保证。

4. **当前算法没有显式 base-case kernel basis**

   标准 BL 的 $K=1$ 会从 $F(0)$ 的 kernel 构造 order-1 basis，再递归组合。当前算法从单位 basis 出发，通过逐阶 discrepancy 消元隐式地产生关系。

5. **当前算法目前没有不同坐标的 degree/shift 设计**

   所有 $P_i$ 用统一 degree bound $m$。标准 BL 可以通过 shift 更自然地表达不同坐标的 degree 权重。

所以现在应当把当前实现理解为：

```text
一个逐阶 approximant relation 搜索算法，
输入 F,K,m，
尝试构造并读取满足 F P = 0 mod delta^K 的低 degree 向量。
```

而不是严格意义上已经完整实现了 Beckermann-Labahn minimal approximant basis。

## 7. 对当前 db 失败的含义

在 `{0,0,1,1,1,1}` sector 上，当前数据是：

```text
r = 14
D = 1000
K_safety = 10
K_cert = 10
```

因此统一 degree 能支持到：

```text
m_max = floor((1001 - 20) / 15) - 1 = 64
```

我临时试过 `m=64`，当前实现仍没有找到关系。

这个结果目前只能说明：

```text
在当前 contribution 构造、当前 masters 输入、当前 approximant 实现下，
m <= 64 没有成功读出关系。
```

它还不能严格说明数学上不存在这样的关系，因为当前 approximant basis 实现还没有达到完整 BL / shifted minimal approximant basis 的标准。

## 8. 下一步需要讨论的问题

接下来需要决定的是：

1. 是否先把固定 sector 的 BL 算法严格实现为 shifted minimal approximant basis；
2. 是否需要给 $P_0,P_i$ 使用不同 degree/shift，而不是统一 $m$；
3. 是否要先用一个小规模可控例子验证：已知存在关系时，当前 basis 提取是否一定能读出；
4. db 的 `{0,0,1,1,1,1}` 失败，到底优先排查算法实现，还是优先排查 masters / contribution 输入是否匹配。
