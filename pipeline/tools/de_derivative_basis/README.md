# DE 导数基约化

这个工具检验已知约化

$$
F_a=\boldsymbol R_a(x)\boldsymbol M
$$

能否改写为

$$
F_a=\sum_{n=0}^{r}\boldsymbol p_{a,n}(x)\boldsymbol M^{(n)},
\qquad \deg\boldsymbol p_{a,n}\leq s.
$$

程序不显式构造可能很大的高阶有理矩阵，而是在有限域普通点计算 DE
矩阵的 Taylor jet，再递推主积分各阶导数。对多项式坐标建立有限域
线性系统，并用没有参与求解的额外点验证。

主要接口为 AnalyzeDerivativeBasisGrid。输入 DE 矩阵、变量、目标约化
行向量的 Association、导数阶列表和多项式次数列表。

结果中：

- UnknownCount 是未经消除冗余的 $N(r+1)(s+1)$；
- DesignRank 是除去当前有限 ansatz 内 syzygy 后的函数空间维数；
- KernelDimension 是两者之差；
- PassedCount 和 FailedCount 记录已有约化是否可由该 ansatz 表示；
- Fits 给出每个目标的一组多项式导数坐标。

Passed 表示训练点精确求解且全部独立 holdout 点通过。它是大有限域上的
高可信恒等式检验，但不是符号恒等式证明。
