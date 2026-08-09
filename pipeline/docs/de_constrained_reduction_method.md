# 利用微分方程约束积分约化结构：方法总览

## 1. 文档目的

本文总结一种用于减少 BL 约化所需展开阶数的方法。核心不是用微分方程
代替约化，而是先用主积分微分方程确定约化系数允许属于什么结构空间，
再用较短的积分展开只确定该空间中的坐标。

本文给出：

1. 问题的统一数学定义；
2. 有限极点、无穷远和 spanning cut 下的核心算法；
3. 已经由实际约化验证的典型结果；
4. 当前方法已经解决和仍未解决的边界。

把同一方法用于建立微分方程的非对角块，只是一个应用：此时待约对象是
$M_i'$，约化结果正是未知 DE 块。它不是本文方法的定义本身。

## 2. 问题定义

在基础域 $F$ 上，设某个 cut 的 $N$ 个主积分组成列向量

$$
\boldsymbol M(\delta)
=
(M_1(\delta),\ldots,M_N(\delta))^{\mathsf T},
$$

并满足微分方程

$$
\frac{d\boldsymbol M}{d\delta}
=A(\delta)\boldsymbol M(\delta),
\qquad
A(\delta)\in\operatorname{Mat}_N(F(\delta)).
$$

待约积分 $I_a$ 的约化为

$$
I_a(\delta)
=\boldsymbol R_a(\delta)\boldsymbol M(\delta),
\qquad
\boldsymbol R_a(\delta)\in F(\delta)^{1\times N}.
$$

对 $\boldsymbol R_a$ 做标准 Apart 后，固定表示为

$$
\boxed{
\boldsymbol R_a(\delta)
=\boldsymbol P_a(\delta)
+\sum_{d\in\mathcal D}
 \sum_{r=1}^{e_d}
 \frac{\boldsymbol n_{a,d,r}(\delta)}{d(\delta)^r}
}
$$

其中：

- $d(\delta)\in F[\delta]$ 首一且不可约；
- $\deg_\delta \boldsymbol n_{a,d,r}<\deg d$；
- $\boldsymbol P_a(\delta)$ 是固定的正则多项式部分；
- $e_d$ 是该约化中因子 $d$ 的最高极点阶数。

若完全不知道结构，就要逐个确定所有 $\boldsymbol n_{a,d,r}$ 和
$\boldsymbol P_a$ 的系数。展开阶数的瓶颈，本质上就是这些待定系数太多。

本文要解决的问题是：在给定候选分母 $d$ 和最高重数 $e_d$ 后，能否只由
DE 预先计算联合分子

$$
(\boldsymbol n_{a,d,e_d},\ldots,\boldsymbol n_{a,d,1})
$$

允许属于的公共有限维空间，并让不同 $I_a$ 只确定其基坐标。

## 3. 核心数学逻辑

### 3.1 DE 给出主积分的局部 jet 空间

先考虑一次因子，并把其零点平移到 $x=0$。主积分在选定的整数 Taylor
区域写成

$$
\boldsymbol M(x)
=\sum_{q=0}^{\infty}\boldsymbol M^{(q)}x^q,
\qquad
\boldsymbol M^{(q)}=B_q\boldsymbol t.
$$

$\boldsymbol t\in F^{r_T}$ 是该区域剩余的自由参数，矩阵 $B_q$ 由完整
Frobenius/Taylor-jet 递推得到。因此 DE 的作用被编码为

$$
(\boldsymbol M^{(0)},\boldsymbol M^{(1)},\ldots)
=(B_0\boldsymbol t,B_1\boldsymbol t,\ldots).
$$

这里必须使用完整 Taylor-jet，而不能只看 $A_{-1}$ 的留数。即使 DE
表面上含 $x^{-2}$、$x^{-3}$，或者存在共振，更高阶递推也可能继续消去
候选自由参数。

当前有限点测试使用以下工作假设：物理解的整数区域没有负整数幂和不可消去
的 log。应联合求解这些消失条件，不能仅按候选 $\mu$ 逐个删除模式。

### 3.2 同一分母的所有极点层必须联合处理

设 $I_a$ 在 $x=0$ 的极点部分为

$$
I_a^{(x)}(x)
=\sum_{r=1}^{d_I}
\frac{\boldsymbol c_{a,r}\boldsymbol M(x)}{x^r},
\qquad
\boldsymbol c_{a,r}\in F^{1\times N}.
$$

若 $I_a$ 在该整数区域也没有负整数幂，则代入主积分展开后，所有负幂必须
相消。$x^{-s}$ 给出

$$
\boxed{
\sum_{r=s}^{d_I}\boldsymbol c_{a,r}B_{r-s}=0,
\qquad
s=1,\ldots,d_I.
}
$$

前几层为

$$
\begin{aligned}
\boldsymbol c_{a,d_I}B_0&=0,\\
\boldsymbol c_{a,d_I-1}B_0
+\boldsymbol c_{a,d_I}B_1&=0,\\
\boldsymbol c_{a,d_I-2}B_0
+\boldsymbol c_{a,d_I-1}B_1
+\boldsymbol c_{a,d_I}B_2&=0.
\end{aligned}
$$

所以 $x^{-1},x^{-2},\ldots$ 的分子一般不能各自独立选基。正确对象是
联合向量

$$
\boldsymbol C_a
=
(\boldsymbol c_{a,d_I},\ldots,\boldsymbol c_{a,1}).
$$

定义块 Toeplitz 矩阵

$$
\mathcal T_{d_I}(B)
=
\begin{pmatrix}
B_0&B_1&\cdots&B_{d_I-1}\\
0&B_0&\cdots&B_{d_I-2}\\
\vdots&\ddots&\ddots&\vdots\\
0&\cdots&0&B_0
\end{pmatrix}.
$$

全部条件等价于

$$
\boxed{
\boldsymbol C_a\mathcal T_{d_I}(B)=0.
}
$$

因此允许的约化结构空间为

$$
\boxed{
\mathcal S_{d_I}
=\operatorname{LeftNull}_F\mathcal T_{d_I}(B),
\qquad
q_{d_I}=\dim_F\mathcal S_{d_I}.
}
$$

若 $\{\boldsymbol V_\beta\}_{\beta=1}^{q_{d_I}}$ 是其基，则任何待约积分
都只能写成

$$
\boldsymbol C_a
=\sum_{\beta=1}^{q_{d_I}}
\lambda_{a\beta}\boldsymbol V_\beta.
$$

原本的 $Nd_I$ 个分子系数被减少为 $q_{d_I}$ 个标量
$\lambda_{a\beta}$。公共向量 $\boldsymbol V_\beta$ 由 DE 决定；展开只用来
确定具体积分的坐标 $\lambda_{a\beta}$。

### 3.3 一般不可约因子

对 $m=\deg d>1$ 的不可约因子，使用局部剩余类域

$$
K_d=F[\delta]/(d),
\qquad t=d(\delta).
$$

在 $K_d[[t]]$ 中构造 $\delta(t)$，把 DE 化为

$$
\frac{d\boldsymbol M}{dt}
=B(t)\boldsymbol M.
$$

标准 Apart 分子 $\boldsymbol n_{a,d,r}(\delta)$ 在 $\delta(t)$ 下也会产生
$t$ 展开。因此先用一个可逆的块上三角映射，把 proper numerators 变成真正
位于 $t^{-r}$ 前的局部向量，再应用完全相同的 Taylor-jet 和 Toeplitz
约束，最后映回标准 Apart 分子。

这一步不会改变全局 Apart 的正则部分。局部极点项产生的 $t^0,t^1,\ldots$
只是同一个有理项确定的局部尾部，不是额外移动出来的自由正则项。

局部结构适合在 $K_d$ 上描述，实际待定系数应回到 $F$ 上计数。若某结构在
$K_d$ 上维数为 $D_K$，则限制标量后

$$
D_F=mD_K.
$$

实际从展开中要确定的是 $D_F$ 个 $F$ 系数。

### 3.4 无穷远和正则多项式

令

$$
z=\frac1\delta.
$$

$\boldsymbol P_a(\delta)$ 变成 $z=0$ 处的负幂部分。无穷远不能简单要求
$I_a$ 的所有负幂都消失，因为目标积分自身可能发散。

记：

- $D_{\rm Reduce}$ 为约化正则多项式的最高次数；
- $D_M$ 为主积分整数区域的最高极点阶数；
- $D_I$ 为目标积分自身允许的最高整数极点阶数。

只有比 $I_a$ 自身更发散的层必须相消，层数为

$$
\boxed{
N_{\rm cancel}
=D_{\rm Reduce}+D_M-D_I.
}
$$

把这些最高 Laurent 层写成主积分 jet 与多项式系数的卷积，再求相应线性
系统的核，便得到正则多项式的联合结构。若 $D_I$ 尚不可靠，就应暂时只给
$\boldsymbol P_a$ 一个次数上界，不能强行使用无穷远约束。

## 4. 实际约化算法

### 4.1 单个 cut

给定当前 cut 的主积分展开、DE、候选分母与重数上界：

1. 分解候选有限分母，并把结果统一成 canonical Apart 口径；
2. 对每个 $d$ 和试探最高重数 $e_d$，由完整 Taylor-jet 构造
   $\mathcal S_{d,e_d}$；
3. 用结构基坐标代替所有无结构 proper-numerator 系数；
4. 若无穷远行为可靠，则同时构造正则部分的联合结构；否则保留给定次数的
   普通多项式系数；
5. 将整个结构化 ansatz 代入短展开，求有限域线性方程；
6. 至少保留一段没有参加求解的高阶展开，检查残差严格为零；
7. 若无解或留出验证失败，只增加相关分母、重数或正则次数，再重新求解。

一次成功重构至少应报告：无结构系数数、结构坐标数、线性系统秩、实际使用
阶数和独立留出区间。

### 4.2 分母和重数不是由局部 DE 自动全部决定的

DE 的严格结论是：给定 $d$ 和 $e_d$ 后，分子必须属于什么空间。当前算法
并未证明只凭当前对角 DE 就能先验确定约化会出现哪些 $d$ 以及最高到几重。

实际候选信息来自：

1. 当前 cut 的 DE 因子；
2. 更高 sector 已知约化所继承的因子；
3. 同 cut、同复杂度或更简单积分给出的保守历史上界；
4. 从低重数开始、由无解或留出失败触发的逐因子升重。

历史结果只能给候选集合与上界，不能读取当前目标积分的精确非零 support，
否则会把事后信息误当成盲重构的压缩。

### 4.3 spanning cut

对目标积分按 sector 从高到低推进。到达 cut $c$ 时，先减去已完成的高层
贡献，得到 residual

$$
I_c^{\rm res}
=I_c-sum_{h>c}\boldsymbol R_h\boldsymbol M_h\big|_c.
$$

然后只把 $I_c^{\rm res}$ 约化到当前 sector 主积分。此 residual 除了当前
DE 的奇点，还会继承高层系数的分母；所以当前层分母不必是当前对角 DE
分母的子集。

若完整块上三角 DE 已知，高层极点关系可在更低 cut 唯一延拓为

$$
\frac1{d^r}
\left(
\boldsymbol v^{\rm high}\boldsymbol M_{\rm high}
+\boldsymbol v^{\rm low}\boldsymbol M_{\rm low}
\right),
$$

从而直接固定继承分母在低层主积分前的分子。若非对角 DE 尚未知，这部分
扩张数据一般不能仅由两个对角块推出；必须由展开确定一次，或者暂时使用
peer/holdout 基并保留退回完整 proper numerator 的路径。

### 4.4 用同一方法约束 DE 本身

若要建立块上三角 DE 的未知非对角块

$$
\frac d{d\delta}
\begin{pmatrix}\boldsymbol H\\\boldsymbol M\end{pmatrix}
=
\begin{pmatrix}A_H&C\\0&A_M\end{pmatrix}
\begin{pmatrix}\boldsymbol H\\\boldsymbol M\end{pmatrix},
$$

则对每个高层主积分定义

$$
J_i
=H_i'-(A_H)_{i*}\boldsymbol H
=C_{i*}\boldsymbol M.
$$

这正是一次普通约化，因而可使用前面的局部结构算法约束 $C_{i*}$。但是
$C$ 是两个微分模之间的扩张数据，存在

$$
C\longmapsto C+G'+GA_M-A_HG
$$

的基变换自由度；两个对角块一般不能唯一决定它。DE 约束可以显著减少
未知数，但新出现的扩张方向仍需由展开确定。

## 5. 已验证的典型结果

### 5.1 Box–triangle 的有限一次因子

在 $6\times6$ 和 $11\times11$ 闭合对角 DE 上，有限线性极点得到以下
联合结构维数。表中数值就是公共基已知后，每个积分仍需确定的标量坐标数。

| DE | 极点 | $d_I=1$ | $d_I=2$ | $d_I=3$ | $d_I=4$ |
|---|---|---:|---:|---:|---:|
| $6\times6$ | $\delta=1$ | 5 | 10 | 15 | 20 |
| $6\times6$ | $\delta=-25/12$ | 3 | 6 | 9 | 12 |
| $11\times11$ | $\delta=-73/75$ | 1 | 2 | 3 | 4 |
| $11\times11$ | $\delta=1$ | 7 | 14 | 21 | 28 |
| $11\times11$ | $\delta=-25/12$ | 5 | 10 | 15 | 20 |

其余三个 $11\times11$ 有限线性极点在 $d_I=1,2,3,4$ 时也分别为
$1,2,3,4$。已有 $42$ 个 $6\times6$ 样本和 $57$ 个 $11\times11$
样本全部属于对应的 DE 结构，包括实际出现的 $d_I=4$ 样本。

数值恰好线性增长不表示不同极点幂次可以独立处理；实际基中已经观察到
幂次之间的非平凡耦合。

### 5.2 高次不可约因子

一般不可约因子的 $K_d$ 局部算法与回到 $F$ 的算法给出一致结果，全部满足
$D_F=(\deg d)D_K$。DB 的 $102$ 个实际样本和 DP topsector 的 $216$ 个
实际样本全部通过。

一个重要的非平凡结果是：DP 有 $6$ 个因子在 DE 中以 $d^{-2}$ 出现，完整
局部递推后允许结构在 $K_d$ 上均只有一维，已有约化全部落在该一维空间中。

### 5.3 无穷远

在 box–triangle 的 $6\times6$ 例子中，主积分整数区域从 $z^{-2}$ 开始并
含两个自由参数。对从 $z^{-2}$ 开始或没有整数区域的目标积分，正次数正则
系数由

$$
6D\longrightarrow4D;
$$

对从 $z^{-4}$ 开始的 `FI`，则为

$$
6D\longrightarrow4D+4.
$$

把每个目标积分自身的允许发散度计入后，$72$ 个含正次正则多项式的约化
全部通过。这同时说明：若错误假设目标积分无负幂，实际具有负幂的样本会明确
失败；因此有限点“$M$ 与 $I$ 均无负整数幂”的测试成功不是自动得到的空结论。

### 5.4 DP 的结构化 spanning-cut 约化

对 `FI31111111`：

- topsector `11111111` 有 $43$ 个主积分、$23$ 个有限分母；有限极点从
  $989$ 个无结构系数降到 $231$ 个 DE 坐标，加上 $258$ 个正则系数后共
  $489$ 个待定量；使用前 $550$ 阶唯一求解，并验证到 $\delta^{1000}$；
- 单零 cut `11101111` 在 topsector 减除后有 $27$ 个主积分。第一版由
  $51$ 个当前 DE 极点坐标、$29$ 个继承极点坐标和 $162$ 个正则坐标组成，
  共 $242$ 个待定量；前 $300$ 阶唯一求解，验证到 $\delta^{2000}$；
- 双零 cut `01101111` 在完整 $115\times115$ 块上三角 DE 已知后，
  $23$ 个继承因子全部稳定且唯一延拓，固定了 $2100$ 个坐标；最终只剩
  $62$ 个 child 对角 DE 坐标和 $84$ 个正则坐标，即 $146$ 个待定量，
  用 $146$ 个标量方程唯一确定并在全部可用阶数上零失配。

这些结果验证了完整逻辑：局部 DE 负责限制当前层分子；完整三角 DE 负责
把高层极点结构延拓到低层；最后只有结构坐标和暂未约束的正则部分需要展开。

### 5.5 用于非对角 DE 的盲测试

在 `01101111` cut 上，只读取父层 DE 和前 $1500$ 个展开系数，未读取旧的
目标非对角块，成功重建了 $43\times14$ 的 top-to-child 块。最后 $13$ 个
复杂行平均需要 $648.54$ 个盲待定系数，最复杂行需要

$$
1209=1078+75+56
$$

个，分别来自新扩张 image、target 齐次结构和正则部分。结果在独立的
$1500$ 至 $2799$ 阶以及另一份 $6000$ 阶数据上均零失配，并与旧精确矩阵
逐项相同。

这一结果说明方法确实也能约束 DE 本身，但也显示其边界：若某行产生许多
新的扩张方向，DE 约束后的盲未知数仍可能很大。事后只看到 $165$ 个非零
系数不能用于盲算，因为事前并不知道它们的 support。

## 6. 当前结论与边界

已经确定的是：

1. 给定分母和最高重数后，DE 能系统、有限维地确定所有极点层的联合分子
   空间；
2. 线性因子和一般不可约因子使用同一个局部 Taylor-jet 原理；
3. 有限极点的实际约化在 DB 与 DP 测试中全部落入预测空间；
4. 无穷远在计入目标积分自身行为后也可按同样的“过度发散项相消”原则处理；
5. spanning cut 中必须同时处理当前 DE 极点和高层继承极点；完整三角 DE
   可以把后者的低层分子固定下来。

仍需继续解决的是：

1. 仅由 DE 自动、可靠地预测约化实际分母集合和最高重数；
2. 在未知 $D_I$ 时如何稳健利用无穷远约束；
3. 如何在不知道精确 support 的前提下进一步压缩新的跨 sector 扩张方向；
4. 将目前已验证的脚本统一成可从 topsector 自动推进到底层的生产工具。

## 7. 详细结果索引

本文只保留方法主线和典型数据。详细推导、完整表格与可复现文件见：

- [有限线性极点的 Frobenius–Toeplitz 推导](../runs/db/de_reduction_pole_constraints/docs/frobenius_structured_pole_constraints.md)
- [有限线性极点的完整结构维数表](../runs/db/de_reduction_pole_constraints/docs/finite_linear_pole_structure_dimensions.md)
- [一般不可约因子的理论、DB/DP 表格与代码入口](../runs/db/de_reduction_pole_constraints/docs/high_degree_irreducible_factor_constraints.md)
- [Box–triangle 的总测试总结](../runs/db/de_reduction_pole_constraints/docs/box_triangle_structure_summary.md)
- [无穷远 $6\times6$ 的完整分类和计数](../runs/db/de_reduction_pole_constraints/docs/infinity_integer_region_structure_dimensions.md)
- [无穷远展开中空与 $6\times6$、$11\times11$ 的比较](../runs/db/de_reduction_pole_constraints/docs/infinity_series_gap_analysis.md)
- [DP spanning-cut 的方法与候选上界策略](../runs/dp/de_structured_spanning_reduce_FI31111111/docs/workflow.md)
- [FI31111111 的 topsector 与单零 cut 结果](../runs/dp/de_structured_spanning_reduce_FI31111111/docs/progress.md)
- [继承 top 极点在 subsector 的关系与完整 DE 对照](../runs/dp/de_structured_spanning_reduce_FI31111111/docs/inherited_top_pole_relations.md)
- [九个 DE 块及非对角块结构分析](../runs/dp/de_structured_spanning_reduce_FI31111111/doublezero_01101111/docs/nine_de_blocks_and_structured_cross_reconstruction.md)
- [01101111 的 1500 阶盲重构](../runs/dp/de_structured_spanning_reduce_FI31111111/doublezero_01101111/docs/blind1500_reconstruction_of_01101111_de.md)
- [01101111 完整约化、146 参数结果与独立验证](../runs/dp/de_structured_spanning_reduce_FI31111111/doublezero_01101111/docs/progress.md)
