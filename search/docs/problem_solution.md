# FI 关于 $\delta$ 的多项式关系搜索与 FI 约化

## 1. Research objective

设 $G$ 是一组待处理的 Feynman 积分标签集合，每个元素记为 $\vec{\nu}$。对每个 $\vec{\nu} \in G$，考虑其一参数形变积分

$$
FI_{\vec{\nu}}(\delta).
$$

已知的是这些积分在若干组独立边界条件下关于 $\delta$ 的一维截断级数展开。目标不是重新计算这些展开，而是利用这些展开搜索积分之间的线性关系，并最终得到原始积分

$$
FI_{\vec{\nu}} := FI_{\vec{\nu}}(1)
$$

之间的约化公式。

整个研究目标分成两个顺序相连的问题。

### 1.1 Problem I: search polynomial relations among $FI_{\vec{\nu}}(\delta)$

搜索非平凡关系

$$
\sum_{\vec{\nu} \in G}
P_{\vec{\nu}}(\delta) \, FI_{\vec{\nu}}(\delta) = 0,
$$

其中

$$
P_{\vec{\nu}}(\delta) = \sum_{k=0}^{m} c_k^{\vec{\nu}} \delta^k,
\qquad c_k^{\vec{\nu}} \in Z_p.
$$

这里 $m$ 是预先给定的多项式次数上界，未知量是全部 $c_k^{\vec{\nu}}$。

### 1.2 Problem II: derive undeformed FI relations and solve reductions

从 Problem I 得到的系数空间零空间出发，构造具体的多项式关系；然后在这些关系中取 $\delta = 1$，得到原始积分 $FI_{\vec{\nu}}$ 之间的线性关系；最后对这些关系再做一次消元，得到把复杂积分表示为简单积分、进一步表示为主积分的约化式。

目标形式是

$$
FI_{\vec{\nu}_{\mathrm{comp}}}
=
\sum_{\vec{\mu} \in M}
R_{\vec{\nu}_{\mathrm{comp}},\vec{\mu}} \, FI_{\vec{\mu}},
$$

其中 $M$ 可以是预先给定的主积分集合，也可以是第二次消元后剩下的自由积分集合。

## 2. Mathematical structures

### 2.1 Integral set

记待研究积分集合为

$$
G = \{\vec{\nu}^{(1)}, \dots, \vec{\nu}^{(N_G)}\}.
$$

这里 $N_G = |G|$。

### 2.2 Boundary conditions

设共有 $N_{\mathrm{bc}}$ 组独立边界条件，编号为

$$
b = 1,2,\dots,N_{\mathrm{bc}}.
$$

在目标应用中，通常希望

$$
N_{\mathrm{bc}} = N_{\mathrm{master}},
$$

并把这组边界条件理解为主积分空间中的一组基。但对本文档的数学形式化而言，只要求这些边界条件线性独立即可。

### 2.3 Truncated series data

对每个 $\vec{\nu} \in G$ 以及每个边界条件 $b$，已知截断展开

$$
FI_{\vec{\nu}}^{(b)}(\delta)
=
\sum_{r=0}^{d} a_r^{\vec{\nu},(b)} \delta^r,
\qquad a_r^{\vec{\nu},(b)} \in Z_p.
$$

其中 $d$ 是已知展开的最高阶。

### 2.4 Polynomial degree bound

给定非负整数 $m$，并只搜索次数不超过 $m$ 的多项式系数：

$$
P_{\vec{\nu}}(\delta) = \sum_{k=0}^{m} c_k^{\vec{\nu}} \delta^k.
$$

通常假设

$$
d > m,
$$

以保证已知级数包含足够信息来约束所搜索的关系。

### 2.5 Complexity ordering on integral labels

为了把最终结果整理成“复杂积分由简单积分表示”的形式，需要在积分标签上固定一个复杂度顺序。

对每个 $\vec{\nu}$ 定义

$$
\mathrm{props}(\vec{\nu}) := \#\{i \mid \nu_i \neq 0\},
$$

以及

$$
\mathrm{dots}(\vec{\nu}) := \sum_i \nu_i - \mathrm{props}(\vec{\nu}).
$$

其中：

- $\mathrm{props}(\vec{\nu})$ 表示非零 propagator 幂次数目；
- $\mathrm{dots}(\vec{\nu})$ 表示总幂次超出 $1$ 的部分。

规定 $\vec{\mu}$ 比 $\vec{\nu}$ 更简单，当且仅当依次满足如下比较规则：

1. $\mathrm{props}(\vec{\mu}) < \mathrm{props}(\vec{\nu})$；
2. 若 $\mathrm{props}$ 相同，则 $\mathrm{dots}(\vec{\mu}) < \mathrm{dots}(\vec{\nu})$；
3. 若两者仍相同，则用固定的标签比较规则打破平局。

在 Problem I 中，变量是 $(\vec{\nu}, k)$，因此在积分标签相同的情况下，再把较大的 $k$ 视为更复杂，从而在消元中优先消去高次 $\delta$ 系数。

## 3. Problem I: polynomial relation search for $FI_{\vec{\nu}}(\delta)$

### 3.1 Unknown variables

对每个 $\vec{\nu} \in G$ 与每个 $0 \le k \le m$，定义未知量

$$
x_{\vec{\nu},k} := c_k^{\vec{\nu}}.
$$

因此未知量总数为

$$
N_{\mathrm{var}} = |G|(m+1).
$$

### 3.2 Polynomial ansatz

目标关系写为

$$
\sum_{\vec{\nu} \in G}
\left(\sum_{k=0}^{m} c_k^{\vec{\nu}} \delta^k\right)
FI_{\vec{\nu}}(\delta) = 0.
$$

对固定的边界条件 $b$，代入截断展开

$$
FI_{\vec{\nu}}^{(b)}(\delta)
=
\sum_{r=0}^{d} a_r^{\vec{\nu},(b)} \delta^r,
$$

得到

$$
\sum_{\vec{\nu} \in G}
\left(\sum_{k=0}^{m} c_k^{\vec{\nu}} \delta^k\right)
\left(\sum_{r=0}^{d} a_r^{\vec{\nu},(b)} \delta^r\right).
$$

展开后，$\delta^n$ 的系数为

$$
\sum_{\vec{\nu} \in G}
\sum_{k=0}^{\min(n,m)}
c_k^{\vec{\nu}} a_{n-k}^{\vec{\nu},(b)}.
$$

由于目标关系必须对每一组边界条件同时成立，因此对所有 $b = 1,\dots,N_{\mathrm{bc}}$ 和 $n = 0,1,\dots,d$ 都要求

$$
\sum_{\vec{\nu} \in G}
\sum_{k=0}^{\min(n,m)}
c_k^{\vec{\nu}} a_{n-k}^{\vec{\nu},(b)} = 0.
$$

### 3.3 Linear system for Problem I

把全部方程堆叠起来，得到齐次线性系统

$$
A^{(\delta)} x = 0.
$$

行指标记为 $(b,n)$，列指标记为 $(\vec{\nu},k)$，则矩阵元为

$$
A^{(\delta)}_{(b,n),(\vec{\nu},k)} =
\begin{cases}
a_{n-k}^{\vec{\nu},(b)}, & n \ge k, \\
0, & n < k.
\end{cases}
$$

矩阵规模为

$$
\#\mathrm{rows} = (d+1)N_{\mathrm{bc}},
\qquad
\#\mathrm{cols} = |G|(m+1).
$$

### 3.4 Output of Problem I

Problem I 的直接输出是零空间

$$
\ker A^{(\delta)}.
$$

等价地，也可以输出以下任一形式：

1. 零空间的一组基向量；
2. $A^{(\delta)}$ 的行最简形以及主元列、自由列信息；
3. 把主变量写成自由变量线性组合的显式表达。

这一阶段的结果仍然是关于系数 $c_k^{\vec{\nu}}$ 的关系，而不是最终的 FI 约化关系。

## 4. Problem II: from $FI_{\vec{\nu}}(\delta)$ relations to $FI_{\vec{\nu}} = FI_{\vec{\nu}}(1)$ reductions

### 4.1 From nullspace basis to concrete polynomial relations

取 $\ker A^{(\delta)}$ 的一组基。为了从参数化的零空间结果得到具体关系，对每一个自由变量分别执行：

1. 令该自由变量取值为 $1$；
2. 令其余自由变量取值为 $0$；
3. 用行最简形回代出全部主变量。

于是得到一组具体系数

$$
\{c_k^{\vec{\nu}}\}_{\vec{\nu} \in G,\, 0 \le k \le m},
$$

并从而得到一条具体的多项式关系

$$
\sum_{\vec{\nu} \in G}
\left(\sum_{k=0}^{m} c_k^{\vec{\nu}} \delta^k\right)
FI_{\vec{\nu}}(\delta) = 0.
$$

### 4.2 Specialization at $\delta = 1$

把上述关系在 $\delta = 1$ 处求值，定义

$$
r_{\vec{\nu}} := \sum_{k=0}^{m} c_k^{\vec{\nu}}.
$$

于是得到原始积分之间的线性关系

$$
\sum_{\vec{\nu} \in G} r_{\vec{\nu}} FI_{\vec{\nu}} = 0.
$$

这一步把系数空间问题转成了原始积分空间问题，是第二次消元的输入。

### 4.3 Linear system for FI relations

把所有由 4.1 和 4.2 得到的 FI 关系堆叠起来，得到第二个齐次线性系统

$$
A^{(FI)} y = 0,
$$

其中向量 $y$ 的分量直接由 $FI_{\vec{\nu}}$ 标记。

每一行对应一条具体 FI 关系；在 $\vec{\nu}$ 这一列上的矩阵元就是该关系中 $FI_{\vec{\nu}}$ 的系数 $r_{\vec{\nu}}$。

### 4.4 Reduction target and masters

对 $A^{(FI)}$ 按积分复杂度顺序做消元后，希望读出的结果形如

$$
FI_{\vec{\nu}_{\mathrm{comp}}}
=
\sum_{\vec{\mu}\text{ simpler}}
R_{\vec{\nu}_{\mathrm{comp}},\vec{\mu}} FI_{\vec{\mu}}.
$$

如果关系系统足够完备，则最后保留下来的自由积分集合可视为主积分集合 $M$，此时上式可写为

$$
FI_{\vec{\nu}_{\mathrm{comp}}}
=
\sum_{\vec{\mu} \in M}
R_{\vec{\nu}_{\mathrm{comp}},\vec{\mu}} FI_{\vec{\mu}}.
$$

因此，求解将 $FI$ 表示为主积分，可以被理解为第二次消元的最终目标；而如果在具体数据中尚不能唯一识别外部给定的主积分集合，那么至少也能得到复杂积分表示为最终自由积分的 reduction 形式。

## 5. Implementation-specific extensions

### 5.1 Problem decomposition

整个任务分为两个顺序阶段：

1. 从截断级数数据中构造并求解关于 $c_k^{\vec{\nu}}$ 的齐次线性系统；
2. 从第一阶段得到的具体关系中导出原始积分之间的关系，并再次消元得到 reduction。

在第二阶段中，可以再细分为三个子步骤：

1. 从零空间得到具体多项式关系；
2. 令 $\delta = 1$ 得到 FI 关系；
3. 再次消元得到 reduction。

### 5.2 Default interpretation of masters

本文档默认采用下面的解释：

- 若先验指定主积分集合，则第二阶段的目标是把其余积分表示成该集合的线性组合；
- 若未先验指定主积分集合，则第二阶段消元后的自由积分集合就是当前数据下自然得到的 master-like set。

### 5.3 Data sufficiency caveat

本文档只形式化：如果已有这些截断级数数据，应如何构造两个线性系统并读取最终关系。

本文档不证明：

1. 给定的 $d$ 和 $m$ 一定足够恢复全部目标关系；
2. 给定的边界条件数一定足够区分所有待搜索关系；
3. 单个有限域 $Z_p$ 上的结果一定足以恢复最终有理系数结构。

这些问题属于后续分析或实现验证问题，不属于当前 phase-a 主体。

## 6. Complete computation flow

### 6.1 Stage I: search polynomial relations for $FI_{\vec{\nu}}(\delta)$

1. 固定积分集合 $G$。
2. 固定独立边界条件 $b = 1,\dots,N_{\mathrm{bc}}$。
3. 读入全部截断级数

   $$
   FI_{\vec{\nu}}^{(b)}(\delta) = \sum_{r=0}^{d} a_r^{\vec{\nu},(b)} \delta^r.
   $$

4. 固定多项式次数上界 $m$。
5. 引入全部未知量 $c_k^{\vec{\nu}}$。
6. 依据第 3 节中的卷积公式构造矩阵 $A^{(\delta)}$。
7. 求解 $A^{(\delta)}x = 0$ 的零空间。
8. 输出零空间基，或者输出其 RREF 与自由变量信息。

### 6.2 Stage II: derive and solve undeformed FI relations

9. 对每一个自由变量，构造一组具体系数赋值：该自由变量取 $1$，其余自由变量取 $0$。
10. 用回代求得全部 $c_k^{\vec{\nu}}$。
11. 对每个 $\vec{\nu}$ 计算

    $$
    r_{\vec{\nu}} = \sum_{k=0}^{m} c_k^{\vec{\nu}}.
    $$

12. 构造具体 FI 关系

    $$
    \sum_{\vec{\nu} \in G} r_{\vec{\nu}} FI_{\vec{\nu}} = 0.
    $$

13. 把所有此类关系堆叠成第二个矩阵 $A^{(FI)}$。
14. 按复杂到简单的顺序排列积分变量。
15. 对 $A^{(FI)}$ 做消元。
16. 读出复杂积分关于简单积分或主积分的线性表示。

## 7. References

- `search/docs/problem_and_workflow.md`

## 8. Review request

请只从数学正确性角度审查本文档，重点检查以下几点：

1. 第 3 节中由截断级数导出 $A^{(\delta)}$ 的公式是否符合预期；
2. 第 4.1 节中由零空间自由变量逐个取值生成具体关系的方法是否就是所需构造；
3. 第 4.2 节中在该阶段令 $\delta = 1$ 是否在目标应用中是正确且充分的；
4. 第 4.4 节中把第二次消元后的自由集合解释为主积分集合的做法是否符合预期含义。

在这四点得到确认之前，不建议进入 phase b。
