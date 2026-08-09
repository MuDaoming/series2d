# 微分方程给出的约化极点结构

文件 `de_reduction_pole_constraints.wl` 在给定线性微分方程

$$
\frac{d\boldsymbol M}{dx}=A(x)\boldsymbol M
$$

及一次极点位置 $x=x_0$ 后，计算约化极点阶数
$d_I=1,\ldots,d_{\max}$ 所允许的联合分子结构。

当前采用的物理解假设是：整数区域中的负整数幂分支已由边界条件去掉，
并且不保留含对数的整数区域分支。因此程序求的是解析 Taylor 解空间

$$
\boldsymbol M(x_0+t)
=\sum_{q\geq 0}\boldsymbol M^{(q)}t^q .
$$

## 调用

```wl
Get["pipeline/tools/de_reduction_pole_constraints/de_reduction_pole_constraints.wl"];

result =
  DEReductionPoleConstraints`ComputePoleConstraintStructures[
    A,
    x,
    x0,
    3,
    "Modulus" -> 2305843009213693951
  ];
```

这里 `A` 是方阵，`x` 是微分变量，`x0` 是极点，`3` 表示同时计算
$d_I=1,2,3$。

关键输出为：

- `result["DEPoleOrder"]`：$A$ 在该点的极点阶数 $d_{\rm DE}$；
- `result["StabilizedQ"]`：保留的 Taylor jet 空间是否已经稳定；
- `result["StabilityHistory"]`：随方程展开阶数增加的稳定性记录；
- `result["TaylorJetBasisRows"]`：允许的 Taylor jet 基；
- `result["Structures"][dI]`：给定 $d_I$ 的联合分子结构。

对固定的 $d_I$，联合分子按

$$
\boldsymbol c
=
\left(
\boldsymbol c_{-d_I},
\boldsymbol c_{-(d_I-1)},
\ldots,
\boldsymbol c_{-1}
\right)
$$

排列。结构中最重要的量是：

- `AllowedJointNumeratorBasis`：所有允许的联合分子向量的一组基；
- `AllowedDimension`：该结构所需的独立标量系数数目；
- `ToeplitzMatrix`：联合分子必须湮灭的 Taylor-jet 矩阵；
- `ToeplitzRank`：独立约束数；
- `EqualsLayerwiseD1ProductQ`：该结构是否只是 $d_I=1$ 结构在每一层的独立复制。

检验一个已有联合分子：

```wl
structure = result["Structures"][2];

DEReductionPoleConstraints`PoleTupleContainedQ[
  Join[cMinus2, cMinus1],
  structure
]
```

返回 `True` 表示已有约化的极点部分满足微分方程给出的全部 Taylor
区域约束。

## 数学口径

若

$$
I_{\rm pole}(t)
=
\sum_{r=1}^{d_I}
\frac{\boldsymbol c_{-r}\boldsymbol M(t)}{t^r},
$$

则 $I$ 在该区域没有负幂要求

$$
\sum_{r=s}^{d_I}
\boldsymbol c_{-r}\boldsymbol M^{(r-s)}
=0,
\qquad
s=1,\ldots,d_I.
$$

程序先由微分方程构造所有允许的
$(\boldsymbol M^{(0)},\ldots,\boldsymbol M^{(d_I-1)})$，再把上述条件
组成块 Toeplitz 线性系统。它一次求解整个
$(\boldsymbol c_{-d_I},\ldots,\boldsymbol c_{-1})$，不会把不同极点
幂次误当作互相独立的 $d_I=1$ 问题。

计算可在有理数域进行（`"Modulus" -> 0`），也可在有限域
$\mathbb F_p$ 上进行。项目中的大规模检验使用后者。
