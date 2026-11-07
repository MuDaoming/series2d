# 固定 X1、积分 X2：定义和测试方法

## 1. 目标与共同输入

原测试将 $X_1/(X_1+X_2)$ 留作参数；本测试将 $X_1$ 本身留作参数，比较两种部分积分家族的主积分数及微分方程复杂度。原始二圈费曼积分 $I_{\boldsymbol\nu}^d$、固定分支积分 $\mathcal I_{\boldsymbol\nu}^{\Delta}$、维度 $d$、传播子指标 $\boldsymbol\nu$、branch 划分和 Gamma 归一化完全沿用 [共同定义](../../j_integral/docs/problem_solution.md) 第 2—3 节。

令 $v=\sum_\alpha\nu_\alpha$，$\sigma_b=\sum_{\alpha\in\mathcal B_b}\nu_\alpha$，$W_{\boldsymbol\nu}=\prod_bX_b^{\sigma_b-1}$。共同的 FBI 表示为

$$
I_{\boldsymbol\nu}^d=\int[dX]\,W_{\boldsymbol\nu}(X)
U(X)^{v-3d/2}\mathcal I_{\boldsymbol\nu}^d(X),
\qquad [dX]=dX_1dX_2dX_3\,\delta_{\rm D}(1-X_1-X_2-X_3).
$$

$U$ 是第一 Symanzik 多项式；$\delta_{\rm D}$ 是 Dirac 函数。这里不引入辅助变形参数。

## 2. 先写三角形上的部分积分，再归一化

直接固定 $X_1=x$、令 $s=X_2$，则 $X_3=1-x-s$。部分积分是

$$
J_{\boldsymbol\nu}^{(X_1)}(x)
=\int_0^{1-x}ds\,
W_{\boldsymbol\nu}(x,s,1-x-s)
U(x,s,1-x-s)^{v-3d/2}
\mathcal I_{\boldsymbol\nu}^d(x,s,1-x-s).
$$

按照两个变量均归一化到 $[0,1]$ 的要求，令 $s=(1-x)y$，于是

$$
X_1=x,\quad X_2=(1-x)y,\quad X_3=(1-x)(1-y),
\qquad ds=(1-x)dy.
$$

定义完整被积函数

$$
G_{\boldsymbol\nu}^{(X_1)}(x,y)
=(1-x)W_{\boldsymbol\nu}(X)
U(X)^{v-3d/2}\mathcal I_{\boldsymbol\nu}^d(X).
$$

因此

$$
\boxed{J_{\boldsymbol\nu}^{(X_1)}(x)=\int_0^1G_{\boldsymbol\nu}^{(X_1)}(x,y)dy,
\qquad I_{\boldsymbol\nu}^d=\int_0^1J_{\boldsymbol\nu}^{(X_1)}(x)dx.}
$$

标准三个 branch 的

$$
U=x(1-x)+(1-x)^2y(1-y).
$$

Jacobian 为 $1-x$，不是原参数化的 $1-z$；这会同时影响目标的权重及 FBI 连接的规范重定义。

## 3. 边界对象及其含义

新内积分的端点是 $y=0$ 对应 $X_2=0$、$y=1$ 对应 $X_3=0$。定义

$$
B_{\boldsymbol\nu,\epsilon}^{(X_1)}(x)
=G_{\boldsymbol\nu}^{(X_1)}(x,\epsilon),\qquad\epsilon=0,1,
$$

其中仅在端点正则时使用直接赋值；否则需相应的正规化定义。在一般 $x$ 处两个端点均有 $U=x(1-x)\ne0$，没有原坐标 $z=1$ 时的 $U=0$ 退化。但 $U\ne0$ 本身不足以保证 FBI 正则；必须检查当前 cut FBI 连接以及目标的端点极点。

本轮输入将两个端点均设为 `regular`；程序只有在共同分母在端点非零时才能执行这一模式，否则会报错，不将端点自动置零。两个边界对象都纳入候选关系。零边界与可相互约化的非零边界单独统计。

## 4. 通用程序与分目录原则

`examples/*.json` 给出 `branch_parameters` 和 `jacobian`，使用引擎变量名 `t,z` 表示 $x,y$。通用引擎检查三个 branch 参数之和为 1，并检查 Jacobian 与坐标行列式相等或相反。积分域确实覆盖目标单纯形、Jacobian 在实区域内取正，则由本节明确的坐标映射保证。

引擎从传播子、质量、Gram 矩阵构造 FBI 模型，不加入 DP 专用公式。原 `j_integral` 没有指定坐标的输入仍使用原映射。新目录的 Makefile 编译共享源码到本目录 `build/`，结果只写入本目录 `results/`，不复制或分叉一套约化算法。

在一般点 $x=x_0$ 展开，先由已知 FBI 方程生成系数，再完成固定 $x_0$ 的 y 积分，随后搜索 $J^{(X_1)}$ 的关系和 DE。展开方法见 [系数生成推导](../../j_integral/docs/delta_vs_t_expansion.md)；此处仅更换输入映射和 Jacobian。

## 5. 比较与验证

先以等质量二圈真空检查新 Jacobian 和已知完整被积函数的连接。随后使用原 DP topsector cut 相同的传播子、运动学、维度和 dot 集合：全部 $\nu_i\ge1$、$\sum_i(\nu_i-1)\le2$，共 45 个 J，加每个目标的两个边界对象，共 135 个候选。

第一组使用 $p=1000003$、$x_0=2/5$、200 阶，第二组改变有限域与展开点。记录训练阶数、留出阶数、关系次数、主对象及 DE 是否闭合。对于闭合的候选结果，再在未参与 Taylor 训练的新 x 点直接检查完整积分证书和 DE。

两种 $J$ 的定义不同，不能逐项比较其 Taylor 系数或把它们的基视为同一组函数；比较的是相同原始家族下部分积分对象的维数和 DE 复杂度。主积分数是在当前 dot、关系次数和验证范围内得到的结果，DE 闭合与多组取样用于检验稳定性。
