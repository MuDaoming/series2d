# 部分分支积分 J(t) 的约化：问题定义与工作流程

## 1. 研究目标

本项目的最终目标是**利用固定分支积分（FBI）表示，实现原始二圈费曼积分的约化**。

三个 branch 的二圈积分在 FBI 表示下只剩两个独立分支参数。本模块先对其中一个参数完成积分，得到一族单变量函数 $J_{\boldsymbol\nu}(t)$，研究它们之间的有理函数约化关系及关于 $t$ 的微分方程。原始积分由最后一次 $t$ 积分恢复。

**计算策略**：在一般常点 $t=t_0$ 生成 $J_{\boldsymbol\nu}(t)$ 的展开，再通过展开系数搜索积分关系。生成展开时使用已知的 FBI 约化和求导公式，不以 $J$ 的约化或微分方程为前提。

生产实现采用 C++/FLINT 的多项式系数递推；Python 有理函数原型只保留作低阶交叉检查。单变量函数的系数数目为 $O(M)$，当前实现仍需处理随阶数增长的多项式，不能称为 $O(M)$ 标量算法；第 9 节给出具体分析。

## 2. 原始费曼积分与分支参数

### 2.1 费曼积分家族

考虑两个圈动量 $\ell_1,\ell_2$ 和 $E$ 个独立外动量 $p_1,\ldots,p_E$。一个家族包含 $N$ 条传播子：

$$
D_\alpha=q_\alpha^2-m_\alpha^2+i0,
\qquad
q_\alpha=c_{\alpha1}\ell_1+c_{\alpha2}\ell_2
+\sum_{a=1}^{E}e_{\alpha a}p_a,
\qquad \alpha=1,\ldots,N.
$$

其中 $c_{\alpha i}$ 和 $e_{\alpha a}$ 是给定的动量组合系数，$m_\alpha$ 是质量，$i0$ 是费曼处方。外部运动学由质量平方及标量积 $p_a\cdot p_b$ 确定。

在时空维度 $d$ 下，积分定义为

$$
\operatorname{FI}_{\boldsymbol\nu}^{d}
=\int\frac{d^d\ell_1}{i\pi^{d/2}}
      \frac{d^d\ell_2}{i\pi^{d/2}}
\prod_{\alpha\in\mathcal A_{\boldsymbol\nu}}D_\alpha^{-\nu_\alpha},
\qquad
\boldsymbol\nu=(\nu_1,\ldots,\nu_N).
$$

本阶段处理无分子的标量积分：$\nu_\alpha$ 为非负整数，活跃传播子集合为

$$
\mathcal A_{\boldsymbol\nu}=\{\alpha:\nu_\alpha>0\},
\qquad
\nu_{\mathrm{tot}}=\sum_{\alpha=1}^{N}\nu_\alpha.
$$

$\nu_\alpha=0$ 表示删除该传播子，不在参数积分中使用 $\Gamma(0)$。维度 $d$ 先按维数正规化解释，进行代数计算时可在一般有理数或有限域值上取样。

### 2.2 Branch 划分

删除传播子的外动量部分后，具有相同圈动量二次结构的传播子属于同一 branch。本阶段考虑三个非空 branch：

$$
\mathcal B_1,\mathcal B_2,\mathcal B_3,
\qquad
\mathcal B_1\sqcup\mathcal B_2\sqcup\mathcal B_3=\{1,\ldots,N\}.
$$

对目标积分定义

$$
n_b=|\mathcal B_b\cap\mathcal A_{\boldsymbol\nu}|,
\qquad
\nu_b=\sum_{\alpha\in\mathcal B_b}\nu_\alpha,
\qquad b=1,2,3.
$$

$n_b$ 是当前 sector 在该 branch 的活跃传播子数，$\nu_b$ 是指标之和，两者不同。本阶段要求每个 branch 都被覆盖，即 $n_b\ge1$。

### 2.3 参数测度

对每个 branch 引入参数 $X_b\ge0$，对每条活跃传播子引入 branch 内参数 $y_\alpha\ge0$，令通常的费曼参数为 $x_\alpha=X_b y_\alpha$。约束为

$$
X_1+X_2+X_3=1,
\qquad
\sum_{\alpha\in\mathcal B_b\cap\mathcal A_{\boldsymbol\nu}}y_\alpha=1.
$$

相应测度定义为

$$
[dX]=\prod_{b=1}^{3}dX_b\,
\delta_{\mathrm D}(1-X_1-X_2-X_3),
$$

$$
[dy]_{\boldsymbol\nu}
=\prod_{\alpha\in\mathcal A_{\boldsymbol\nu}}dy_\alpha
\prod_{b=1}^{3}\delta_{\mathrm D}
\left(1-\sum_{\alpha\in\mathcal B_b\cap\mathcal A_{\boldsymbol\nu}}y_\alpha\right).
$$

$\delta_{\mathrm D}$ 表示狄拉克函数；本模块不引入旧方法中控制积分域的辅助参数 $\delta$。

## 3. FBI 表示

### 3.1 U、F 和 R

对合并后的传播子作圈动量配方：

$$
\sum_\alpha X_{\operatorname{br}(\alpha)}y_\alpha
(q_\alpha^2-m_\alpha^2)
=\sum_{i,j=1}^{2}A_{ij}\ell_i\cdot\ell_j
+2\sum_{i=1}^{2}\ell_i\cdot B_i+C_{\mathrm{kin}}.
$$

其中 $\operatorname{br}(\alpha)$ 是传播子的 branch 编号；$A$ 是 $2\times2$ 对称矩阵，$B_i$ 是外动量的线性组合，$C_{\mathrm{kin}}$ 是标量项。利用 branch 内约束，$A$ 只依赖于 $X_1,X_2,X_3$。

定义

$$
U(X)=\det A(X),
\qquad
F(y;X)=\sum_{i,j=1}^{2}\operatorname{adj}(A)_{ij}\,B_i\cdot B_j-U C_{\mathrm{kin}}.
$$

$\operatorname{adj}(A)$ 是伴随矩阵。在 branch 约束面上，$F$ 可写成关于 $y$ 的二次型：

$$
F(y;X)=\frac12y^{\mathsf T}R(X)y.
$$

$R$ 是对称矩阵，其行列对应活跃传播子；矩阵表达式在约束面外不必唯一。对于标准三个 branch 对应 $\ell_1,\ell_2,\ell_1+\ell_2$ 的选择，

$$
U=X_1X_2+X_2X_3+X_3X_1.
$$

### 3.2 固定分支积分的定义

以 $\Delta$ 表示 FBI 的维度参数。固定 $X=(X_1,X_2,X_3)$ 后，定义

$$
\mathcal I_{\boldsymbol\nu}^{\Delta}(X)
=\frac{(-1)^{\nu_{\mathrm{tot}}}\Gamma(\nu_{\mathrm{tot}}-\Delta)}
{\prod_{\alpha\in\mathcal A_{\boldsymbol\nu}}\Gamma(\nu_\alpha)}
\int[dy]_{\boldsymbol\nu}
\frac{\prod_{\alpha\in\mathcal A_{\boldsymbol\nu}}y_\alpha^{\nu_\alpha-1}}
{(F(y;X)-i0)^{\nu_{\mathrm{tot}}-\Delta}}.
$$

$\mathcal I$ 表示 FBI，$\operatorname{FI}$ 表示原始费曼积分。二圈标量积分的物理 FBI 参数为 $\Delta=d$，但 FBI 递推中也使用 $\Delta=d+s$，其中 $s$ 为整数维度偏移。

完成两个圈动量积分后，原始积分为

$$
\boxed{
\operatorname{FI}_{\boldsymbol\nu}^{d}
=\int[dX]\,
W_{\boldsymbol\nu}(X)\,
U(X)^{\nu_{\mathrm{tot}}-3d/2}\,
\mathcal I_{\boldsymbol\nu}^{d}(X)
},
\qquad
W_{\boldsymbol\nu}(X)=\prod_{b=1}^{3}X_b^{\nu_b-1}.
$$

这里的 Gamma 函数和符号已包含在 FBI 定义中，不能为不同 $\boldsymbol\nu$ 任意删除不同的归一化因子。

### 3.3 Sector、S 矩阵与 MFBI

sector 定义为二进制向量

$$
s_\alpha=\begin{cases}1,&\nu_\alpha>0,\\0,&\nu_\alpha=0.\end{cases}
$$

corner 指标为该二进制向量本身。设当前 sector 有 $n$ 条活跃传播子，定义关联矩阵 $E_{b\alpha}=1$ 当传播子 $\alpha$ 属于 branch $b$，否则为零。则

$$
S=\begin{pmatrix}0_{3\times3}&E\\E^{\mathsf T}&R\end{pmatrix}
$$

是 $(3+n)\times(3+n)$ 矩阵。

令辅助系数 $\kappa_b,\zeta_\alpha,\zeta_0$ 满足

$$
S\begin{pmatrix}\kappa_1\\\kappa_2\\\kappa_3\\\zeta_1\\\vdots\\\zeta_n\end{pmatrix}
=\begin{pmatrix}\zeta_0\\\zeta_0\\\zeta_0\\0\\\vdots\\0\end{pmatrix},
\qquad C=\kappa_1+\kappa_2+\kappa_3.
$$

$\zeta_\alpha$ 是递推辅助系数，与积分变量 $z$ 不同。非奇异时取 $\zeta_0=1$；奇异时从零空间选择 $\zeta_0=0$ 的解。

沿用 `expand` 的编号：

| Case | S | C | FBI 层的处理 |
|---|---|---|---|
| 0 | 非奇异 | 非零 | corner 是一个 MFBI；其余指标降 dot、对齐维度 |
| 1 | 非奇异 | 零 | 降指标并迁移维度，约到子 sector |
| 2 | 奇异 | 非零 | 降指标，约到子 sector |
| 3 | 奇异 | 零 | 用非零 $\zeta_\alpha$ 转移指标，约到子 sector |

设保留的 MFBI corner 指标为 $\boldsymbol\mu_1,\ldots,\boldsymbol\mu_r$，$r$ 是 FBI 基的维数。它不等于待搜索的 $J$ masters 数目。

## 4. 从原始积分到 J(t)

### 4.1 变量变换

令

$$
X_1=(1-z)t,\qquad X_2=(1-z)(1-t),\qquad X_3=z,
\qquad 0\le t,z\le1.
$$

该变换自动满足分支参数约束。消去狄拉克函数后的 Jacobian 为

$$
[dX]=(1-z)\,dt\,dz.
$$

以下 $U(t,z)$、$W_{\boldsymbol\nu}(t,z)$ 和 $\mathcal I_{\boldsymbol\nu}^{\Delta}(t,z)$ 均表示将上述变换代入后的函数。

### 4.2 被积函数与部分积分

定义完整被积函数

$$
G_{\boldsymbol\nu}(t,z)
=(1-z)W_{\boldsymbol\nu}(t,z)
U(t,z)^{\nu_{\mathrm{tot}}-3d/2}
\mathcal I_{\boldsymbol\nu}^{d}(t,z).
$$

定义研究对象

$$
\boxed{J_{\boldsymbol\nu}(t)=\int_0^1G_{\boldsymbol\nu}(t,z)\,dz},
\qquad
\boxed{\operatorname{FI}_{\boldsymbol\nu}^{d}=\int_0^1J_{\boldsymbol\nu}(t)\,dt}.
$$

$z$ 的积分从一开始就是完整的 $[0,1]$。$J$ 不是在一条移动边界上取值的旧 BFI，也不是通过截断积分域定义的变形积分。

### 4.3 约化目标

设 $\mathcal T$ 是待处理的指标集合。希望从中选出主函数 $J_{\boldsymbol\eta_1},\ldots,J_{\boldsymbol\eta_m}$，使

$$
J_{\boldsymbol\nu}(t)
=\sum_{a=1}^{m}\mathscr R_{\boldsymbol\nu,a}(t)
J_{\boldsymbol\eta_a}(t),\qquad
\mathscr R_{\boldsymbol\nu,a}(t)\in\mathbb K(t).
$$

$\mathbb K$ 是固定运动学和维度后的系数域。有限域测试时取 $\mathbb K=\mathbb F_p$，其中 $p$ 是大素数。

约化之后，进一步研究主函数列向量 $\boldsymbol J$ 的方程

$$
\frac{d\boldsymbol J}{dt}=A_J(t)\boldsymbol J.
$$

$A_J$ 是要从展开中求出的结果，不是展开输入。

## 5. 用于展开的已知 FBI 关系

### 5.1 FBI 约化与求导

FBI 的已知关系包含 IBP 和维度迁移：

$$
S\begin{pmatrix}T_1^\Delta\\T_2^\Delta\\T_3^\Delta\\
\nu_1\mathcal I_{\boldsymbol\nu+\boldsymbol e_1}^{\Delta}\\\vdots\\
\nu_n\mathcal I_{\boldsymbol\nu+\boldsymbol e_n}^{\Delta}\end{pmatrix}
=\begin{pmatrix}-\mathcal I_{\boldsymbol\nu}^{\Delta-1}\\
-\mathcal I_{\boldsymbol\nu}^{\Delta-1}\\-\mathcal I_{\boldsymbol\nu}^{\Delta-1}\\
\mathcal I_{\boldsymbol\nu-\boldsymbol e_1}^{\Delta-1}\\\vdots\\
\mathcal I_{\boldsymbol\nu-\boldsymbol e_n}^{\Delta-1}\end{pmatrix},
$$

$$
C\mathcal I_{\boldsymbol\nu}^{\Delta-1}
=(2\Delta-\nu_{\mathrm{tot}}-3)\zeta_0\mathcal I_{\boldsymbol\nu}^{\Delta}
+\sum_\alpha\zeta_\alpha\mathcal I_{\boldsymbol\nu-\boldsymbol e_\alpha}^{\Delta-1}.
$$

$T_b^\Delta$ 是 IBP 辅助量，$\boldsymbol e_\alpha$ 是传播子指标空间的单位向量。对分支坐标 $u=t$ 或 $u=z$，有

$$
\partial_u\mathcal I_{\boldsymbol\nu}^{\Delta}
=-\frac12\sum_{\alpha,\beta}\partial_uR_{\alpha\beta}
\,h_{\alpha\beta}(\boldsymbol\nu)
\mathcal I_{\boldsymbol\nu+\boldsymbol e_\alpha+\boldsymbol e_\beta}^{\Delta+1},
$$

$$
h_{\alpha\beta}(\boldsymbol\nu)
=\begin{cases}\nu_\alpha\nu_\beta,&\alpha\ne\beta,\\
\nu_\alpha(\nu_\alpha+1),&\alpha=\beta.\end{cases}
$$

求和只涉及当前 sector 的活跃传播子。上述公式用于 FBI 层，不要求 $J$ 已经约化。

### 5.2 FBI 基与完整被积函数的表示

定义裸 MFBI 列向量

$$
\boldsymbol i(t,z)
=\bigl(\mathcal I_{\boldsymbol\mu_1}^{d},\ldots,
\mathcal I_{\boldsymbol\mu_r}^{d}\bigr)^{\mathsf T}.
$$

通过已知 FBI 递推得到有理行向量 $\boldsymbol\rho_{\boldsymbol\nu}$ 及有理矩阵 $C_t,C_z$：

$$
\mathcal I_{\boldsymbol\nu}^{d}=\boldsymbol\rho_{\boldsymbol\nu}\boldsymbol i,
\qquad \partial_t\boldsymbol i=C_t\boldsymbol i,
\qquad \partial_z\boldsymbol i=C_z\boldsymbol i.
$$

为吸收 $U$ 的非整数幂并包含分支权重，定义

$$
P_j(t,z)=(1-z)W_{\boldsymbol\mu_j}(t,z)
U(t,z)^{|\boldsymbol\mu_j|},\qquad
P=\operatorname{diag}(P_1,\ldots,P_r),
$$

$$
\boldsymbol f(t,z)=U(t,z)^{-3d/2}P(t,z)\boldsymbol i(t,z).
$$

$|\boldsymbol\mu_j|$ 表示指标之和。$f_j$ 正是该 MFBI corner 对应的完整被积函数。于是

$$
G_{\boldsymbol\nu}=\boldsymbol r_{\boldsymbol\nu}\boldsymbol f,
\qquad
\boldsymbol r_{\boldsymbol\nu}
=(1-z)W_{\boldsymbol\nu}U^{\nu_{\mathrm{tot}}}
\boldsymbol\rho_{\boldsymbol\nu}P^{-1},
$$

$$
\partial_u\boldsymbol f=A_u\boldsymbol f,
\qquad
A_u=PC_uP^{-1}
+\operatorname{diag}\left(
\frac{\partial_uP_j}{P_j}-\frac{3d}{2}\frac{\partial_uU}{U}
\right),\quad u=t,z.
$$

$\boldsymbol r_{\boldsymbol\nu},A_t,A_z$ 都是有理函数。这些 $A_t,A_z$ 描述被积函数的 FBI 基，与未知的 $A_J$ 不同。

## 6. 常点逐阶展开

### 6.1 展开点和初始函数

选择一般常点 $t_0$，定义局部坐标 $\tau=t-t_0$。要求积分族在该点解析，并可在一个收敛域内将 $t$ 微分移入 $z$ 积分，再按维度解析延拓。

标准三个 branch 的 $U$ 为

$$
U(t,z)=(1-z)\bigl[z+(1-z)t(1-t)\bigr].
$$

$t=0,1$ 会在 $z=0$ 引入额外退化，故本阶段选择避开这类退化和运动学奇点的常点。

定义初始函数列向量

$$
\boldsymbol f_0(z)=\boldsymbol f(t_0,z),\qquad
\frac{d\boldsymbol f_0}{dz}=A_0(z)\boldsymbol f_0(z),
\qquad A_0(z)=A_z(t_0,z).
$$

$\boldsymbol f_0(z)$ 仍是关于 $z$ 的函数。它不能在积分前被替换成一组任意数值常数。

### 6.2 FBI 展开的固定宽度递推

定义矩阵级数

$$
\boldsymbol f(t_0+\tau,z)=\Phi(\tau,z)\boldsymbol f_0(z),
\qquad \Phi(\tau,z)=\sum_{n=0}^{\infty}\Phi_n(z)\tau^n,
\qquad \Phi_0=1_r.
$$

$1_r$ 是 $r\times r$ 单位矩阵。通分已知 FBI $t$ 连接，写为

$$
Q(t,z)\partial_t\boldsymbol f=\mathcal N(t,z)\boldsymbol f,
\qquad \mathcal N=Q A_t.
$$

$Q$ 为标量多项式，$\mathcal N$ 为多项式矩阵。平移后有有限展开

$$
Q(t_0+\tau,z)=\sum_{j=0}^{d_Q}Q_j(z)\tau^j,
\qquad
\mathcal N(t_0+\tau,z)=\sum_{j=0}^{d_{\mathcal N}}\mathcal N_j(z)\tau^j.
$$

比较 $\tau^n$ 得

$$
\boxed{
\Phi_{n+1}
=\frac{1}{(n+1)Q_0}
\left[
\sum_{j=0}^{\min(d_{\mathcal N},n)}\mathcal N_j\Phi_{n-j}
-\sum_{j=1}^{\min(d_Q,n)}Q_j(n-j+1)\Phi_{n-j+1}
\right]}.
$$

$Q_0(z)$ 必须不是零多项式。$d_Q,d_{\mathcal N}$ 由输入决定，不随展开阶数增长。

### 6.3 目标被积函数的系数

对每个目标，取公分母 $d_{\boldsymbol\nu}(t,z)$，使

$$
\boldsymbol r_{\boldsymbol\nu}
=\frac{\boldsymbol b_{\boldsymbol\nu}}{d_{\boldsymbol\nu}}.
$$

$\boldsymbol b_{\boldsymbol\nu}$ 为多项式行向量。分别将分子、分母平移并按 $\tau$ 写成有限多项式，其系数记为 $\boldsymbol b_{\boldsymbol\nu,j}(z)$、$d_{\boldsymbol\nu,j}(z)$。

定义有理行向量系数 $\boldsymbol v_{\boldsymbol\nu,n}(z)$：

$$
G_{\boldsymbol\nu}(t_0+\tau,z)
=\sum_{n=0}^{\infty}\boldsymbol v_{\boldsymbol\nu,n}(z)
\boldsymbol f_0(z)\tau^n.
$$

由 $d_{\boldsymbol\nu}\boldsymbol v_{\boldsymbol\nu}
=\boldsymbol b_{\boldsymbol\nu}\Phi$ 得

$$
\boxed{
\boldsymbol v_{\boldsymbol\nu,n}
=\frac{1}{d_{\boldsymbol\nu,0}}
\left[
\sum_{j=0}^{\min(\deg_\tau\boldsymbol b_{\boldsymbol\nu},n)}
\boldsymbol b_{\boldsymbol\nu,j}\Phi_{n-j}
-\sum_{j=1}^{\min(\deg_\tau d_{\boldsymbol\nu},n)}
d_{\boldsymbol\nu,j}\boldsymbol v_{\boldsymbol\nu,n-j}
\right]}.
$$

最后，真正需要的 $J$ 展开系数为

$$
\boxed{
J_{\boldsymbol\nu}(t_0+\tau)=\sum_{n=0}^{\infty}c_{\boldsymbol\nu,n}\tau^n,
\qquad
c_{\boldsymbol\nu,n}=\int_0^1
\boldsymbol v_{\boldsymbol\nu,n}(z)\boldsymbol f_0(z)\,dz}.
$$

$\Phi_n$ 和 $\boldsymbol v_{\boldsymbol\nu,n}$ 是本试验算法的中间对象；$c_{\boldsymbol\nu,n}$ 才是标量展开系数。

### 6.4 C++ 中的多项式递推

第 6.2、6.3 节给出数学递推，生产程序不直接保存这些有理函数。选择第 7 节的固定极点多项式 $q(z)$，取固定整数 $e\ge1$，使 $Q_0$ 和所有目标的 $d_{\boldsymbol\nu,0}$ 都整除 $H=q^e$。定义

$$
F_n=H^n\Phi_n,\qquad
V_{\boldsymbol\nu,n}=H^{n+1}\boldsymbol v_{\boldsymbol\nu,n}.
$$

这里 $F_n$ 是多项式矩阵，$V_{\boldsymbol\nu,n}$ 是多项式行。其递推为

$$
nF_n=\sum_{j=0}^{\min(d_{\mathcal N},n-1)}
\mathcal N_jH^j\frac{H}{Q_0}F_{n-1-j}
-\sum_{j=1}^{\min(d_Q,n)}
(n-j)Q_jH^{j-1}\frac{H}{Q_0}F_{n-j},
$$

$$
V_{\boldsymbol\nu,n}=
\sum_{j=0}^{\min(\deg\boldsymbol b_{\boldsymbol\nu},n)}
\boldsymbol b_{\boldsymbol\nu,j}H^j
\frac{H}{d_{\boldsymbol\nu,0}}F_{n-j}
-\sum_{j=1}^{\min(\deg d_{\boldsymbol\nu},n)}
d_{\boldsymbol\nu,j}H^{j-1}
\frac{H}{d_{\boldsymbol\nu,0}}V_{\boldsymbol\nu,n-j}.
$$

所有递推系数都是预先算好的固定多项式。逐阶循环只做有限域多项式运算和标量 $1/n$ 乘法，不做符号有理函数通分、约分，不做级数乘级数卷积。只需保留固定宽度的历史窗口；积分输入为 $V_{\boldsymbol\nu,n}/q^{e(n+1)}$。

## 7. 完整 z 积分的系数处理

### 7.1 总导数恒等式

对于关于 $z$ 的有理行向量 $\boldsymbol g$，定义

$$
\nabla_z\boldsymbol g=\boldsymbol g'+\boldsymbol g A_0.
$$

由 $\boldsymbol f_0'=A_0\boldsymbol f_0$，有

$$
\partial_z(\boldsymbol g\boldsymbol f_0)
=(\nabla_z\boldsymbol g)\boldsymbol f_0.
$$

将每个系数行写为

$$
\boldsymbol v_{\boldsymbol\nu,n}
=\boldsymbol h_{\boldsymbol\nu,n}+\nabla_z\boldsymbol g_{\boldsymbol\nu,n},
$$

其中 $\boldsymbol h$ 是保留的积分余项，$\boldsymbol g$ 是有理证书。于是

$$
c_{\boldsymbol\nu,n}
=\int_0^1\boldsymbol h_{\boldsymbol\nu,n}\boldsymbol f_0\,dz
+\left[\boldsymbol g_{\boldsymbol\nu,n}\boldsymbol f_0\right]_{z=0}^{z=1}.
$$

程序可在有理函数域逐项检验 $\boldsymbol v-\boldsymbol h-\nabla_z\boldsymbol g=0$。该检验验证代数恒等式；端点项仍须按积分的边界行为处理。

### 7.2 有限极点的逐层消去

设 $q(z)$ 为覆盖所有系数极点的首一无重因子多项式，$s=\deg q$。当前系数积分算法要求

$$
A_0(z)=\frac{B(z)}{q(z)},\qquad \deg B<s.
$$

$B$ 是多项式矩阵，该条件对应有限点的 Fuchs 型连接及无穷远无多项式部分。写当前待处理行向量为 $\boldsymbol v=\boldsymbol P/q^k$。

当 $k\ge2$ 时，选择

$$
\boldsymbol g=\frac{\boldsymbol a}{q^{k-1}},\qquad \deg\boldsymbol a<s.
$$

则

$$
\nabla_z\boldsymbol g
=\frac{q\boldsymbol a'+\boldsymbol a B-(k-1)q'\boldsymbol a}{q^k}.
$$

令 $\boldsymbol a$ 满足

$$
\boldsymbol a\bigl(B-(k-1)q'1_r\bigr)
\equiv\boldsymbol P\pmod q,
$$

即可将最高极点层降到 $k-1$。这是至多 $rs$ 个未知数的有限域线性系统。

若线性系统无解，选择像空间部分进行消去，并将不能消去的余项保留为固定点积分对象；不能将余项设为零。

### 7.3 无穷远部分的消去

降到 $\boldsymbol P/q$ 后，若 $\deg\boldsymbol P=d_P\ge s-1$，取

$$
\boldsymbol g=\boldsymbol a z^m,\qquad m=d_P-s+1.
$$

记 $B_{s-1}$ 为 $B$ 的 $z^{s-1}$ 系数矩阵，则 $q\nabla_z\boldsymbol g$ 的最高次系数为

$$
\boldsymbol a(m1_r+B_{s-1}).
$$

解该线性系统以消去最高次分子。奇异时同样保留不能消去的方向，继续处理更低次数。

### 7.4 固定点积分与端点常数

统一定义积分对象

$$
\Pi_{kij}=\int_0^1\frac{z^j}{q(z)^k}f_{0,i}(z)\,dz,
$$

其中 $k$ 为分母幂次，$i$ 是 FBI 基分量，$j$ 是分子幂次。正则端点还可能产生

$$
E_{ei}=f_{0,i}(e),\qquad e\in\{0,1\}.
$$

把所有实际保留的 $\Pi_{kij}$ 和 $E_{ei}$ 统一记为 $C_\lambda$。每阶系数表示为

$$
c_{\boldsymbol\nu,n}=\sum_\lambda a_{\boldsymbol\nu,\lambda,n}C_\lambda.
$$

$C_\lambda$ 固定在同一个 $t_0$，不随 $n$ 改变。它们是系数积分的表示对象，不是已经确定的 $J$ masters；其集合也可能包含进一步的依赖。

正则端点显式加入 $\boldsymbol g(e)\boldsymbol f_0(e)$。只有在适当解析正则化域内证明该边界项为零后，才可通过解析延拓使用零边界条件。不能仅根据连接有极点就自动丢弃边界项。

### 7.5 C++ 的极点分层实现

对多项式分子作唯一的 $q$ 进展开

$$
V=\sum_{j\ge0}v_j(z)q(z)^j,\qquad \deg v_j<\deg q.
$$

程序用分治多项式除法转换这些系数。第 7.2 节一个极点层的消元只会改变当前和下一层，不必每次重写整个高次多项式。相同极点层的逆映射、余项映射和下一层更新矩阵在所有目标、所有阶数之间复用。无穷远部分使用固定宽度的移位更新。

**核方向产生的附加关系也必须加入。**若有限极点或无穷远消元矩阵奇异，不能只保留其像空间之外的余项。对核向量构造对应的 $g$，最高层虽抵消，$\nabla_z g$ 的低层部分仍给出

$$
\int_0^1(\nabla_z g)f_0\,dz-[gf_0]_0^1=0.
$$

程序将这些精确恒等式转换为固定点 period/端点坐标之间的关系，行化简后统一作用于全部展开系数。输出的 `period_completion` 记录核方向、生成函数及原始关系。当前默认检查多项式次数 0–16 和极点层 2–17 的核，覆盖已测试算例的 n=2 无穷远共振；这一有限上限不能自动解释为所有输入的完备性证明。

遗漏这一步曾使 DB 和 DP 的 J 基分别多保留一个对象；只验证高阶系数和 DE 闭合无法发现这种冗余。当前回归包含专门的核方向关系测试及独立 Blade FI 对照。

`--verify-all` 对每个系数重建多项式总导数证书并精确检查；默认检查每个目标的前五阶和最后一阶。检查数写入结果，不能把抽查数称作全部系数验证数。

## 8. 关系搜索、sector 与 J 的 DE

### 8.1 多项式关系搜索

给定有限目标集合 $\mathcal T$ 及最高多项式次数 $b$，搜索

$$
\sum_{\boldsymbol\nu\in\mathcal T}P_{\boldsymbol\nu}(t)J_{\boldsymbol\nu}(t)=0,
\qquad
P_{\boldsymbol\nu}(t_0+\tau)=\sum_{k=0}^{b}p_{\boldsymbol\nu,k}\tau^k.
$$

代入第 7.4 节的表示，对每个固定点对象 $C_\lambda$ 和展开阶数 $n$，要求

$$
\sum_{\boldsymbol\nu\in\mathcal T}\sum_{k=0}^{\min(b,n)}
p_{\boldsymbol\nu,k}a_{\boldsymbol\nu,\lambda,n-k}=0.
$$

这给出有限域齐次线性系统。矩阵行是 $(\lambda,n)$，列是 $(\boldsymbol\nu,k)$。训练阶数用于求候选关系，未参与求解的更高阶用于验证。

在所有 $C_\lambda$ 坐标上成立是积分关系成立的充分条件。若这些固定点对象含未消除的依赖，该表示可能漏掉关系，因此有限阶搜索得到的自由对象先称为候选 masters。

### 8.2 按 sector/cut 处理

可保留原项目按支持集和 cut 分层的组织方式。给定 cut 支持集 $c$，FBI 递推中只保留支持集包含 $c$ 的源项。topsector cut 是 $c=(1,\ldots,1)$ 的情形。

cut 模式生成的是 $J$ 在所选 FBI 投影上的数据，不等于完整 $J$。比较新旧 masters 数时，必须保持 cut、运动学、目标指标范围和允许对象类别一致。

### 8.3 搜索主函数的 DE

若得到候选主函数集合，对其展开求导：

$$
J_{\boldsymbol\eta_a}'(t_0+\tau)
=\sum_{n=0}^{\infty}(n+1)c_{\boldsymbol\eta_a,n+1}\tau^n.
$$

随后搜索该导数到候选主函数集合的有理关系。找不到关系可能表示目标集合尚未闭合、多项式次数不足、展开阶数不足或系数积分表示尚未消除足够依赖；不能直接宣布不存在 DE。

### 8.4 边界对象：哪些 BFI/BBFI 进入这里

对于正则 $z$ 端点，定义

$$
B^{(z)}_{a,\boldsymbol\nu}(t)=G_{\boldsymbol\nu}(t,a),\qquad a=0,1.
$$

这些是完成 $z$ 积分时出现的 BFI 型对象，与 $J(t)$ 具有同一个自变量。其展开也由已知 FBI 递推获得：在第 6.4 节的多项式行中取 $z=a$，再除以 $H(a)^{n+1}$，坐标落在相同的端点常数 $E_{ai}$ 上。程序名称为 `Bz0{...}`、`Bz1{...}`。不能把它们设成与 $J$ 无关的随机投影，也不能丢掉非零端点常数。

原二维积分的另外两条边固定 $t=0,1$，仍保留 $z$；双边界 BBFI 则对应两个坐标同时到端点。它们不是一般 $t$ 上的另一组任意函数，不能把旧变域 $F(\delta)$ 的四个 BFI、四个 BBFI 标签直接复制为 $J(t)$ 的独立 masters。

奇异端点必须先说明极限或正则化。当前 `regulated_zero` 表示在适当解析正则化域令该端点的总导数边界项为零，再解析延拓；相应直接端点对象在这一约定下为零。它不计算非零渐近首项或有限部分。如需那类边界对象，必须另行定义并实现。

DB 当前测试的 $z=0$ 是正则端点，显式保留；$z=1$ 使用上述零边界约定。最早 DB 原型将两端都设为零，其结果已经被本轮测试取代。

## 9. 复杂度与当前方法需要验证的问题

设 $M$ 为展开截断阶数，$N$ 始终表示传播子数。固定 FBI 基维数、目标集合和第 6 节多项式次数后：

- 第 6.2、6.3 节每增加一阶，只执行固定数量的 $\mathbb K(z)$ 上有理函数运算，因此是 $O(M)$ 次有理函数域运算。
- 本算法的中间对象 $\Phi_n(z)$、$\boldsymbol v_{\boldsymbol\nu,n}(z)$ 仍依赖 $z$，其分子分母次数可能随 $n$ 增长；一次域运算不能当成一次有限域标量运算。
- 第 7 节每阶还需处理相应的极点层与积分余项。其成本必须单独计入，不能因最后只有一个 $t$ 变量而忽略。

生产 C++ 实现通过第 6.4、7.5 节消除了逐阶符号有理函数计算，但 $F_n,V_n$ 的多项式次数仍可按 $n$ 增长，极点层数也按 $n$ 增长。固定输入下，递推的多项式系数总工作量为二次规模；系数积分还包括多项式进制转换的成本。因此尚未实现 $O(M)$ 标量算法，不能据此承诺百万阶在一秒内完成。实际耗时见 [测试记录](testing.md)。

## 10. 完整计算流程

```text
输入：
  两圈传播子及三个 branch 的划分
  质量平方和外动量 Gram 矩阵
  时空维度 d、有限域素数 p
  常点 t0、目标指标集合、可选 cut
  展开阶数 M、关系次数上限、训练/留出阶数
  各端点的处理条件

步骤 0：构造 FBI 家族
  从传播子构造 U、F、R 和各 sector 的 S。
  用已知 FBI 关系确定 MFBI、目标约化行和 FBI 连接。

步骤 1：定义完整被积函数
  代入 X1=(1-z)t、X2=(1-z)(1-t)、X3=z。
  纳入 Jacobian、分支权重、U 的幂及 FBI 归一化。

步骤 2：常点递推
  令 tau=t-t0，初始化 Phi0=1。
  用第 6 节递推计算 Phi_n 和目标系数行 v_n。

步骤 3：完整系数积分
  在固定 t0 使用 FBI z 连接，构造 h_n 和有理证书 g_n。
  检验 v_n=h_n+g_n'+g_n A0，加入正确的端点项。
  得到固定点积分对象前的有限域系数。

步骤 4：搜索 J 关系
  用训练前缀求多项式关系，在留出阶数逐项验证。
  增加目标范围、次数和展开阶数，检查候选 master 数是否稳定。

步骤 5：搜索 DE 与独立验证
  对候选 masters 的展开求导，搜索导数关系。
  检查闭合性，变更 prime/t0，并与独立积分或已知关系比较。

输出：
  J 的常点展开坐标、系数积分证书、候选约化及 DE
  cut/目标范围、训练与验证阶数、失败项和各阶段耗时
```

## 11. 当前适用范围与验证要求

物理前端面向三个 branch、非负传播子指标的二圈标量家族。拓扑、质量、运动学、展开点和目标均由输入给出。一般 FBI 矩阵连接也可作为独立输入。

当前系数积分后端要求第 7.2 节的连接形式；不满足时应明确报告，需要另行处理连接或改进算法。奇异端点非零广义边界项也未包含在当前后端中。

验证先做 vac 的完整对象，再做 db 的 topsector cut。每个算例至少检查：FBI 连接相容性、系数积分证书、关系留出阶数、候选 masters 对目标集合和次数的稳定性，以及独立正确性依据。低阶通过不能代替 master 完备性检查。

## 12. 符号索引

| 符号 | 定义与含义 |
|---|---|
| $d,\Delta$ | 原始积分的时空维度、FBI 递推的维度参数；物理 FBI 取 $\Delta=d$ |
| $N,E$ | 传播子总数、选用的独立外动量数 |
| $\ell_i,p_a,q_\alpha,D_\alpha,m_\alpha$ | 圈动量、外动量、传播子动量、传播子分母、质量，见第 2.1 节 |
| $\boldsymbol\nu,\nu_{\mathrm{tot}},\nu_b$ | 传播子指标向量、指标总和、branch 内指标之和 |
| $\mathcal A_{\boldsymbol\nu},\mathcal B_b,n_b$ | 活跃传播子集合、第 $b$ 个 branch、该 branch 的活跃传播子数 |
| $X_b,y_\alpha,x_\alpha$ | 分支参数、分支内参数、通常的费曼参数；$x_\alpha=X_{\operatorname{br}(\alpha)}y_\alpha$ |
| $[dX],[dy]_{\boldsymbol\nu}$ | 包含各自狄拉克约束的参数测度，见第 2.3 节 |
| $A,B_i,C_{\mathrm{kin}},U,F,R$ | 圈动量配方矩阵及其线性、常数项，两个参数多项式及 FBI 二次型矩阵，见第 3.1 节 |
| $\operatorname{FI}_{\boldsymbol\nu}^{d},\mathcal I_{\boldsymbol\nu}^{\Delta}$ | 原始费曼积分、固定分支积分；两者的测度和归一化分别见第 2、3 节 |
| $W_{\boldsymbol\nu}$ | 分支权重 $\prod_b X_b^{\nu_b-1}$ |
| $S,C,\kappa_b,\zeta_\alpha,\zeta_0$ | FBI 递推矩阵及辅助系数，见第 3.3 节 |
| $t,z$ | 变换后的两个分支坐标；$z=X_3$ 已在 $J$ 中完整积分，$t$ 保留 |
| $G_{\boldsymbol\nu},J_{\boldsymbol\nu}$ | 完整的 $dt\,dz$ 被积函数、完成 $z$ 积分后的待约化函数 |
| $\boldsymbol\mu_j,r,\boldsymbol i$ | MFBI corner 指标、MFBI 数目、裸 MFBI 列向量 |
| $P,\boldsymbol f,\boldsymbol r_{\boldsymbol\nu}$ | 对角归一化矩阵、含完整权重的 FBI 列向量、目标到该向量的有理行 |
| $C_t,C_z,A_t,A_z$ | 裸 FBI 连接和归一化后的 FBI 连接；均属于已知 FBI 层 |
| $t_0,\tau,\boldsymbol f_0,A_0$ | 一般展开常点、$t-t_0$、固定在 $t_0$ 的 FBI 函数列及其 $z$ 连接 |
| $Q,\mathcal N,\Phi_n,\boldsymbol v_{\boldsymbol\nu,n}$ | FBI 固定宽度递推的通分多项式、分子矩阵、FBI 展开矩阵系数、目标展开行系数 |
| $c_{\boldsymbol\nu,n}$ | 已完成 $z$ 积分的 $J$ 的第 $n$ 阶标量系数 |
| $q,B,\nabla_z,\boldsymbol h,\boldsymbol g$ | 固定点极点多项式、连接分子矩阵、总导数算子、积分余项、有理积分证书 |
| $\Pi_{kij},E_{ei},C_\lambda$ | 固定点积分、正则端点值，以及二者的统一编号；它们不等于 $J$ masters |
| $\mathcal T,\boldsymbol\eta_a,m,A_J$ | 搜索目标集合、候选主函数指标、主函数数目、最终要得到的 $J$ 微分方程矩阵 |
| $\mathbb K,p,M,b$ | 系数域、有限域素数、展开截断阶数、关系多项式次数上限 |

局部复用的字母均以所在小节为准：例如第 7 节的 $B(z)$ 是固定点连接的分子矩阵，区别于第 3 节的外动量向量 $B_i$；第 3.3 节的 $E_{b\alpha}$ 是关联矩阵元素，区别于外动量数 $E$ 和端点值 $E_{ei}$。

## 13. 与原 expand 路线的对应及当前缺口

原 `expand` 路线从已知 FBI 的二元多项式微分关系直接生成二维数值系数。按旧方法定义的辅助参数和积分区域，对这些系数进行积分或边界投影，得到用于关系搜索的一维数据。它不需要事先约化外层 FI，也不需要先得到外层 FI 的 DE。

本模块保持“已知 FBI 关系 → 展开 → 搜索外层对象关系”的顺序，改变的是外层对象：这里从定义上完整积掉 $z$，搜索 $J(t)$ 的关系，不再搜索一般变形对象 $F(\delta)$ 的关系。

最早的 Python 原型只展开 $t$，逐阶构造关于 $z$ 的符号有理函数再积分，造成了不可接受的开销。C++ 已改为多项式分子递推和缓存的极点分层积分，前处理、展开、关系搜索均可直接运行 C++。仍未解决的是严格的线性标量成本；目前的高阶多项式和极点层总量仍增长，不能通过先输入待求的 $J$ DE 来回避这个问题。

生产代码已按第 6.4、7.5 节改为 C++ 多项式递推和极点分层消元；原型只用于回归。实际测试结果、未闭合项及耗时记录在 [测试记录](testing.md)，不将有限阶候选数当作完备的主积分数。

## 14. 参考文献与项目内依据

1. L.-H. Huang, R.-J. Huang, Y.-Q. Ma, *Tame multi-leg Feynman integrals beyond one loop*, arXiv:2412.21053：FBI 表示、S 矩阵及递推关系。
2. `docs/background/B001_feynman_integral_and_fbi.md`：积分、branch 与 FBI 的符号约定。
3. `expand/docs/problem_solution.md`：已知 FBI 关系和直接逐阶递推。
4. `search/docs/problem_solution.md`：由多 projection 级数搜索多项式关系。
5. `original_workflow_review.md`：原有各层方法与实现的对应记录。
