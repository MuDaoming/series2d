# 旧 δ 展开与当前 J(t) 展开：定义、推导与实际算法

本文回答一个具体问题：**没有预先求出待约化对象的 DE 时，旧 `expand` 怎样得到 $F(\delta)$ 的系数，当前 `j_integral` 又怎样得到 $J(t)$ 的系数？**

两者都从已知 FBI 约化和求导公式出发，都不需要先得到外层积分的约化或 DE。但是两者的系数生成算法确实不同：旧方法对二维局部展开作有限的单项式积分；当前方法对每个 $t$ 系数完成整个 $z$ 区间的积分，程序为此增加了固定 $t_0$ 的积分约化步骤。下面逐步说明，不能把这个差别简化成“二维变成一维，删去一个展开指标”。

## 1. 共同起点：原始 FI 和 FBI

考虑无分子的二圈标量积分

$$
I_{\boldsymbol\nu}^{d}
=\int\frac{d^d\ell_1}{i\pi^{d/2}}
      \frac{d^d\ell_2}{i\pi^{d/2}}
\prod_{\alpha\in\mathcal A}D_\alpha^{-\nu_\alpha},
\qquad
D_\alpha=q_\alpha^2-m_\alpha^2+i0.
$$

这里 $\ell_1,\ell_2$ 是圈动量，$q_\alpha$ 是圈动量和外动量的指定线性组合，$m_\alpha$ 是质量，$d$ 是时空维度，$\nu_\alpha$ 是非负整数，$\mathcal A=\{\alpha:\nu_\alpha>0\}$ 是活跃传播子集合。令 $v=\sum_\alpha\nu_\alpha$。

将传播子按圈动量结构分成三个 branch $\mathcal B_b$，$b=1,2,3$，并令 $\sigma_b=\sum_{\alpha\in\mathcal B_b}\nu_\alpha$。本文考虑三个 branch 都非空的情形。通常的费曼参数写为

$$
x_\alpha=X_b y_\alpha\quad(\alpha\in\mathcal B_b),
\qquad \sum_bX_b=1,
\qquad \sum_{\alpha\in\mathcal A\cap\mathcal B_b}y_\alpha=1.
$$

定义测度和 branch 权重

$$
[dX]=\prod_{b=1}^3dX_b\,\delta_{\rm D}(1-X_1-X_2-X_3),
\qquad W_{\boldsymbol\nu}(X)=\prod_bX_b^{\sigma_b-1},
$$

$$
[dy]_{\boldsymbol\nu}=\prod_{\alpha\in\mathcal A}dy_\alpha
\prod_b\delta_{\rm D}\left(1-\sum_{\alpha\in\mathcal A\cap\mathcal B_b}y_\alpha\right).
$$

$\delta_{\rm D}$ 是 Dirac 函数，与后面的辅助变形参数 $\delta$ 不同。令 $U$、$F$ 为原始积分的第一、第二 Symanzik 多项式代入 $x_\alpha=X_b y_\alpha$ 后的表达式。$U$ 只依赖 $X$；在 branch 约束面上写成 $F(y;X)=\frac12y^{\mathsf T}R(X)y$，这定义了所用对称矩阵 $R$。标准三个 branch 的 $U=X_1X_2+X_2X_3+X_3X_1$。

固定分支积分定义为

$$
\mathcal I_{\boldsymbol\nu}^{\Delta}(X)
=\frac{(-1)^v\Gamma(v-\Delta)}{\prod_{\alpha\in\mathcal A}\Gamma(\nu_\alpha)}
\int[dy]_{\boldsymbol\nu}
\frac{\prod_{\alpha\in\mathcal A}y_\alpha^{\nu_\alpha-1}}
{[F(y;X)-i0]^{v-\Delta}}.
$$

$\Delta$ 是 FBI 的维度参数；物理二圈积分对应 $\Delta=d$。于是

$$
\boxed{I_{\boldsymbol\nu}^{d}
=\int[dX]W_{\boldsymbol\nu}(X)
U(X)^{v-3d/2}\mathcal I_{\boldsymbol\nu}^{d}(X).}
$$

Gamma 归一化已包含在 $\mathcal I$ 中。不同指标之间不能任意删除不同的归一化因子。更完整的传播子、$U,F,R$ 和 FBI 关系定义见 [problem_solution.md](problem_solution.md) 第 2—3 节。

## 2. 两种方法已经知道什么？

共同的输入是 **FBI 层**的约化、维度递推和参数求导公式。例如，对任意 branch 坐标 $s$，直接对 FBI 定义求导得到

$$
\partial_s\mathcal I_{\boldsymbol\nu}^{\Delta}
=-\frac12\sum_{\alpha,\beta}
(\partial_s R_{\alpha\beta})
\nu_\alpha(\nu_\beta+\delta_{\alpha\beta})
\mathcal I_{\boldsymbol\nu+e_\alpha+e_\beta}^{\Delta+1}.
$$

$e_\alpha$ 是第 $\alpha$ 个单位指标向量，$\delta_{\alpha\beta}$ 是 Kronecker 符号。右侧的高 dot、移维 FBI 再由已知 FBI 关系处理。这里没有输入 $F(\delta)$ 的 DE，也没有输入 $J(t)$ 的 DE。

将所需 FBI 选成列向量 $\mathbf i$ 后，这些已知关系也可整理成

$$
\partial_s\mathbf i=C_s\mathbf i,
\qquad \mathcal I_{\boldsymbol\nu}^{d}=\rho_{\boldsymbol\nu}\mathbf i,
$$

其中 $C_s$ 是有理函数矩阵，$\rho_{\boldsymbol\nu}$ 是行向量；必要的下 sector 按当前完整家族或 cut 一并处理。这个 FBI 方程与以后搜索出来的 $J$ 方程是不同层次的对象。

## 3. 旧 expand：先二维局部展开，再作 δ 投影

### 3.1 坐标和完整被积函数

旧方形参数化取

$$
X_1=u,\quad X_2=(1-u)v_2,\quad X_3=(1-u)(1-v_2),
\qquad 0\le u,v_2\le1.
$$

这里 $v_2$ 是坐标，与指标总和 $v$ 区分。Jacobian 为 $1-u$。定义完整外层被积函数

$$
H_{\boldsymbol\nu}(u,v_2)
=(1-u)W_{\boldsymbol\nu}(X)U(X)^{v-3d/2}
\mathcal I_{\boldsymbol\nu}^{d}(X).
$$

选一般常点 $(a,b)$，置 $\xi=u-a$、$\eta=v_2-b$，展开

$$
H_{\boldsymbol\nu}(a+\xi,b+\eta)
=\sum_{p,q\ge0}h_{\boldsymbol\nu,pq}\xi^p\eta^q.
$$

### 3.2 二维系数怎么来？

`expand` 使用已知 FBI 关系直接递推二维系数。它吸收共同的 $U$ 幂并处理移维标签，再将关系变成固定多项式系数的递推，而非每阶操作越来越复杂的符号有理函数。

例如，某条已经清分母的关系为

$$
D(\xi,\eta)g(\xi,\eta)=\sum_i N_i(\xi,\eta)f_i(\xi,\eta),
\qquad D_{00}\ne0.
$$

$D,N_i$ 是固定多项式，$g,f_i$ 是待递推 FBI 或其重定义对象。比较系数得到

$$
g_{pq}=\frac1{D_{00}}\left(
\sum_i[N_if_i]_{pq}
-\sum_{(j,k)\ne(0,0)}D_{jk}g_{p-j,q-k}\right).
$$

负指标系数记为零。求导关系中 $[\partial_\xi f]_{pq}=(p+1)f_{p+1,q}$，另一方向类似。因此从已知低阶数据和固定点 MFBI 初值递推下一层，而不需要外层积分的 DE。MFBI 指 master FBI，与原始 FI 或 $J$ 的 master 不同。

对固定家族、固定目标集合，多项式移位的数量有固定上界，每个格点做固定量的标量运算。总次数不超过 $M$ 的格点数为 $(M+1)(M+2)/2$，这是旧方法 $O(M^2)$ 展开的来源。最后乘 branch 权重和 Jacobian，得到 $h_{\boldsymbol\nu,pq}$。

### 3.3 引入 δ：缩放两个积分区间

定义

$$
F_{\boldsymbol\nu}(\delta)
=\int_{-a\delta}^{(1-a)\delta}d\xi
 \int_{-b\delta}^{(1-b)\delta}d\eta\,
H_{\boldsymbol\nu}(a+\xi,b+\eta).
$$

因此 $F_{\boldsymbol\nu}(1)=I_{\boldsymbol\nu}^{d}$。这是对象定义；不要求在 $\delta=0$ 的 Taylor 级数直接收敛到 1，必要时需解析延拓。

定义常数矩

$$
m_p(a)=\frac{(1-a)^{p+1}-(-a)^{p+1}}{p+1}.
$$

单项式的积分恰好是

$$
\int_{-a\delta}^{(1-a)\delta}\xi^p d\xi
=m_p(a)\delta^{p+1}.
$$

故

$$
\boxed{[\delta^n]F_{\boldsymbol\nu}
=\sum_{p+q=n-2}h_{\boldsymbol\nu,pq}m_p(a)m_q(b).}
$$

**每个 δ 系数只取二维系数表的一条有限对角线。不存在“先把某个完整区间的积分约化掉”这一步。** 两个积分所需的信息已经由 $\delta$ 次数分级了。

### 3.4 原来的 BFI、BBFI

沿用同一个完整被积函数 $H$，将一个或两个变量固定在移动端点。以 $\xi=s_a\delta$ 为例，$s_a=-a$ 或 $1-a$，则

$$
\int_{-b\delta}^{(1-b)\delta}
H(a+s_a\delta,b+\eta)d\eta
=\sum_{p,q}h_{pq}s_a^pm_q(b)\delta^{p+q+1}.
$$

两个变量都固定时，系数为 $h_{pq}s_a^ps_b^q\delta^{p+q}$。这就是旧投影中 FI、BFI、BBFI 分别增加 2、1、0 个 $\delta$ 次数的原因；它们复用同一个二维局部展开。

## 4. 当前 J：只展开 t，但 z 已经完整积分

### 4.1 定义和展开点

新坐标为

$$
X_1=(1-z)t,\quad X_2=(1-z)(1-t),\quad X_3=z,
\qquad 0\le t,z\le1.
$$

Jacobian 为 $1-z$。定义

$$
G_{\boldsymbol\nu}(t,z)
=(1-z)W_{\boldsymbol\nu}(X)U(X)^{v-3d/2}
\mathcal I_{\boldsymbol\nu}^{d}(X),
$$

$$
\boxed{J_{\boldsymbol\nu}(t)=\int_0^1G_{\boldsymbol\nu}(t,z)dz,
\qquad I_{\boldsymbol\nu}^{d}=\int_0^1J_{\boldsymbol\nu}(t)dt.}
$$

选一般常点 $t_0$，令 $\tau=t-t_0$。通常不选 $t=0,1$，因为它们对应 branch 退化边界。要求所用连接和目标分母在 $t_0$ 处可以建立所需的局部展开；$z$ 端点的正规化另行处理。

$$
G_{\boldsymbol\nu}(t_0+\tau,z)=\sum_{n\ge0}g_{\boldsymbol\nu,n}(z)\tau^n,
\qquad
J_{\boldsymbol\nu}(t_0+\tau)=\sum_{n\ge0}c_{\boldsymbol\nu,n}\tau^n,
$$

$$
\boxed{c_{\boldsymbol\nu,n}=\int_0^1g_{\boldsymbol\nu,n}(z)dz.}
$$

交换展开、微分和积分先在共同的收敛及正规化区域解释，再作解析延拓。这里的积分区间始终是 $[0,1]$，没有随 $\tau$ 缩小。

### 4.2 已知 FBI 方程怎样产生 g 的系数？

为了包括共同的非整数 $U$ 幂，将 FBI 向量重定义为

$$
\mathbf f=U^{-3d/2}P\mathbf i,
$$

其中 $P$ 为可逆对角有理矩阵。例如第 $j$ 个 FBI corner 指标为 $\boldsymbol\mu_j$ 时，取

$$
P_{jj}=(1-z)W_{\boldsymbol\mu_j}U^{|\boldsymbol\mu_j|}.
$$

于是目标和已知连接写为

$$
G_{\boldsymbol\nu}=r_{\boldsymbol\nu}\mathbf f,
\qquad
\partial_t\mathbf f=A_t\mathbf f,
\qquad \partial_z\mathbf f=A_z\mathbf f,
$$

$$
r_{\boldsymbol\nu}=(1-z)W_{\boldsymbol\nu}U^v\rho_{\boldsymbol\nu}P^{-1},
$$

$$
A_s=P C_sP^{-1}
+\operatorname{diag}\left(\frac{\partial_sP_{jj}}{P_{jj}}
-\frac{3d}{2}\frac{\partial_sU}{U}\right),\qquad s=t,z.
$$

所有 $r,A$ 都是有理函数。若实现选用了等价的重定义，以下推导不变。

令 $\mathbf f_0(z)=\mathbf f(t_0,z)$。定义归一化的 $t$ 方向传输矩阵

$$
\mathbf f(t_0+\tau,z)=\Phi(\tau,z)\mathbf f_0(z),
\qquad \Phi=\sum_{n\ge0}\Phi_n(z)\tau^n,
\qquad \Phi_0=\mathbf1.
$$

将已知连接清分母，写成

$$
Q(\tau,z)\partial_\tau\mathbf f=\mathcal N(\tau,z)\mathbf f,
\quad Q=\sum_jQ_j(z)\tau^j,
\quad\mathcal N=\sum_j\mathcal N_j(z)\tau^j.
$$

这里 $Q$ 是标量多项式，$\mathcal N$ 是多项式矩阵，两者 $\tau$ 次数固定。比较系数，对 $n\ge1$ 得

$$
\boxed{nQ_0\Phi_n
=\sum_{j=0}^{n-1}\mathcal N_j\Phi_{n-1-j}
-\sum_{j=1}^{n}Q_j(n-j)\Phi_{n-j}.}
$$

超出多项式次数的项为零。再将目标行向量写为 $r_{\boldsymbol\nu}=b_{\boldsymbol\nu}/d_{\boldsymbol\nu}$，分子、分母也按 $\tau$ 展开。令

$$
r_{\boldsymbol\nu}\Phi=\sum_{n\ge0}v_{\boldsymbol\nu,n}(z)\tau^n.
$$

则

$$
\boxed{d_{\boldsymbol\nu,0}v_{\boldsymbol\nu,n}
=\sum_{j=0}^n b_{\boldsymbol\nu,j}\Phi_{n-j}
-\sum_{j=1}^n d_{\boldsymbol\nu,j}v_{\boldsymbol\nu,n-j}.}
$$

这两条递推只用了已知 FBI 输入。它们生成的是

$$
g_{\boldsymbol\nu,n}(z)=v_{\boldsymbol\nu,n}(z)\mathbf f_0(z),
$$

**还不是已经积分好的标量 $c_{\boldsymbol\nu,n}$。** $v_{\boldsymbol\nu,n}$ 中的有理函数来自保留了整个 $z$ 变量；旧二维标量递推没有这个保留变量。

### 4.3 C++ 为什么不直接递推符号有理函数？

选择包含所需 $z$ 极点的平方自由多项式 $q(z)$，再取 $H=q^e$，使 $Q_0$ 和所有目标的 $d_{\boldsymbol\nu,0}$ 都整除 $H$。定义多项式矩阵和行向量

$$
F_n=H^n\Phi_n,
\qquad V_{\boldsymbol\nu,n}=H^{n+1}v_{\boldsymbol\nu,n}.
$$

上一节递推等价于

$$
nF_n=\sum_{j=0}^{n-1}\mathcal N_jH^j\frac{H}{Q_0}F_{n-1-j}
-\sum_{j=1}^n Q_j(n-j)H^{j-1}\frac{H}{Q_0}F_{n-j},
$$

$$
V_{\boldsymbol\nu,n}
=\sum_{j=0}^n b_{\boldsymbol\nu,j}H^j\frac{H}{d_{\boldsymbol\nu,0}}F_{n-j}
-\sum_{j=1}^n d_{\boldsymbol\nu,j}H^{j-1}
\frac{H}{d_{\boldsymbol\nu,0}}V_{\boldsymbol\nu,n-j}.
$$

固定多项式乘子提前计算，主循环只做多项式加乘、标量除法，并用固定长度的历史窗口。当前 C++ 正是这一路线，避免了 Python 原型逐阶进行有理函数约分的开销。但 $F_n,V_n$ 的 $z$ 次数仍会随 $n$ 增长。

## 5. 当前多出来的步骤：固定 t₀ 的完整 z 积分

### 5.1 积分系数需要额外的约化

在固定点上，已知 FBI 方程给出

$$
\frac{d\mathbf f_0}{dz}=A_0(z)\mathbf f_0,
\qquad A_0=A_z(t_0,z).
$$

对有理行向量 $a(z)$ 定义

$$
\nabla a=a'+aA_0.
$$

于是

$$
\frac{d}{dz}(a\mathbf f_0)=(\nabla a)\mathbf f_0.
$$

若将某个系数的有理行向量分解成

$$
v_{\boldsymbol\nu,n}=h_{\boldsymbol\nu,n}+\nabla a_{\boldsymbol\nu,n},
$$

就得到

$$
\boxed{c_{\boldsymbol\nu,n}
=\int_0^1h_{\boldsymbol\nu,n}(z)\mathbf f_0(z)dz
+[a_{\boldsymbol\nu,n}(z)\mathbf f_0(z)]_0^1.}
$$

这就是当前提到的 **$z$ 方向 IBP**：它发生在生成 $t$ 展开系数的过程中，用来完成每阶剩下的完整区间积分。这里既没有预先求出 $J(t)$ 的约化，也没有输入 $J(t)$ 的 DE；但确实引入了旧 δ 投影不需要的固定点积分约化。

### 5.2 程序实际保留什么？

在当前实现所处理的正规奇点表示中，写 $A_0=B(z)/q(z)$。$B$ 是多项式矩阵。系数最初形如 $V_n/q^{e(n+1)}$。程序通过上述全导数恒等式降低极点阶数和无穷远多项式次数，剩下若干固定积分及端点值，例如

$$
\Pi_{kij}=\int_0^1\frac{z^j}{q(z)^k}f_{0,i}(z)dz,
\qquad E_{\epsilon i}=f_{0,i}(\epsilon),\quad\epsilon=0,1.
$$

$i$ 标记 FBI 分量，$j$ 标记分子次数，$k$ 标记极点阶数。保留下来的独立坐标统记为 $C_\lambda$，于是输出为

$$
c_{\boldsymbol\nu,n}=\sum_\lambda a_{\boldsymbol\nu,\lambda,n}C_\lambda.
$$

程序并不数值积分得到这些 $C_\lambda$，也不把它们任意赋成随机数；它计算的是对所有这些坐标成立的系数关系。有限域取样用于有理代数运算，不是在有限域中直接计算物理积分或 Gamma 函数。

例如，设 $a=u/q^{k-1}$，则

$$
\nabla a=
\frac{qu'+uB-(k-1)q'u}{q^k}.
$$

因此降极点的最高层线性映射来自

$$
u\bigl(B-(k-1)q'\mathbf1\bigr)\pmod q.
$$

这里公式中的 $u$ 是行多项式，与第 3 节旧坐标 $u$ 无关。最高层映射若有核，不能只把不能消去的余项当成独立 master；核方向的全导数还可能产生低层积分之间的关系。当前 `complete_relations` 会生成这类关系并再次消元。这一步曾使 DB、DP 的过完备坐标进一步减少，因此它是正确性所需的步骤。

当前默认的核关系补充范围是无穷远次数 $0\ldots16$、有限极点阶数 $2\ldots17$。这是实际实现的有限搜索范围，不是任意家族完备性的证明；扩大样本、改变有限域、独立验证仍然必要。

正则端点项必须保留；使用正规化后端点为零的模式时，需要相应的解析依据。不能因为端点奇异或某个因子趋零就直接删去 $a\mathbf f_0$ 的边界项。端点函数 $B_{\boldsymbol\nu,\epsilon}(t)=G_{\boldsymbol\nu}(t,\epsilon)$ 与固定点端点坐标 $E_{\epsilon i}$ 也不是同一个对象：前者是一族待搜索关系的函数，后者是这些函数展开所用的固定初值坐标。

## 6. 一个可手算的例子：看到额外积分发生在哪里

用标量函数 $G(t,z)=1/(1-tz)$ 演示算法结构；它只是教学例子，不代表具体的费曼积分家族。取 $0<t_0<1$，记 $D=1-t_0z$，则

$$
G(t_0+\tau,z)=\sum_{n\ge0}\frac{z^n}{D^{n+1}}\tau^n,
\qquad c_n=\int_0^1\frac{z^n}{D^{n+1}}dz.
$$

已知的“被积函数层”求导关系是 $\partial_tG=zG/(1-tz)$、$\partial_zG=tG/(1-tz)$。先由它们得到上面的局部系数，仍须积分才能得到 $c_n$。

对 $n\ge1$，利用

$$
\frac{d}{dz}\frac{z^n}{D^n}=\frac{nz^{n-1}}{D^{n+1}},
\qquad
\frac{t_0z^n}{D^{n+1}}
=\frac{z^{n-1}}{D^{n+1}}-\frac{z^{n-1}}{D^n},
$$

积分后得到

$$
\boxed{t_0c_n+c_{n-1}=\frac{(1-t_0)^{-n}}{n}.}
$$

初始积分 $c_0=-\log(1-t_0)/t_0$，右侧来自端点。这个例子中，固定点积分的 IBP 恰好直接给出了常数成本的标量递推；不需要先输入 $J$ 的 DE。一般矩阵家族能否同样得到这样高效的闭合系数递推，是需要证明和实现的额外问题，不能从“只有一个展开变量”直接推出。

## 7. 为什么不能把旧算法直接删去一个指标？

为排除两种坐标不同带来的干扰，在新坐标中也作二维局部展开。取 $z_0$ 为常点，写

$$
G(t_0+\tau,z_0+\zeta)=\sum_{p,q\ge0}a_{pq}\tau^p\zeta^q.
$$

如果仍把两个积分区间都按 δ 缩放，那么

$$
[\delta^n]F(\delta)=\sum_{p+q=n-2}a_{pq}m_p(t_0)m_q(z_0),
$$

仍然只用有限对角线。但是当前要的是

$$
[\tau^p]J(t_0+\tau)
=\int_0^1[\tau^p]G(t_0+\tau,z)dz.
$$

若局部 $z$ 级数可在整个区间逐项积分，它等于

$$
\sum_{q=0}^{\infty}a_{pq}m_q(z_0).
$$

**固定一个 $p$ 后，仍要所有 $q$，而非仅仅 $q\le M-p$。** 若局部级数不覆盖整个区间，还需要解析延拓。有限二维截断不能直接给出精确的 $J$ 系数。

也可定义一个同时包含两种操作的对象

$$
K(\tau,\delta)=\int_{-z_0\delta}^{(1-z_0)\delta}
G(t_0+\tau,z_0+\zeta)d\zeta
=\sum_{p,q}a_{pq}\tau^pm_q(z_0)\delta^{q+1}.
$$

当前对象是 $J(t_0+\tau)=K(\tau,1)$；将 $\delta$ 设成 1 需要完成整个 $z$ 积分的信息，不能以截取某一有限次数代替。当前固定点 IBP 正是处理这部分信息的一种办法。

## 8. 得到系数后才搜索 J 关系与 DE

对候选对象（包括需要的端点函数）搜索

$$
\sum_\nu P_\nu(t)J_\nu(t)=0.
$$

将多项式写为 $P_\nu(t_0+\tau)=\sum_{k=0}^{b}p_{\nu k}\tau^k$，代入固定积分坐标后，逐个 $\lambda,n$ 要求

$$
\sum_\nu\sum_{k=0}^{\min(b,n)}
p_{\nu k}a_{\nu,\lambda,n-k}=0.
$$

这是有限域线性代数。解出关系、选择 $J$ 的主积分后，利用 $J'$ 的系数 $(n+1)c_{n+1}$ 搜索它们的 DE，并以未参与求解的高阶系数等检验。

若固定积分坐标还有未发现的关系，逐坐标要求为零会过强，可能漏掉真实的 $J$ 关系，导致 master 计数偏大。因此不能只凭少量 Taylor 阶数或一个候选零空间宣称完备。

这里和之后用 $J$ 的约化、DE 对最后一次 $t$ 积分做 IBP、得到原始 FI 约化，是两个阶段。后一个阶段不属于本文的系数生成输入。

## 9. 复杂度：哪些结论已经成立？

固定家族、目标集合和有限域，令 $M$ 为最高展开阶数。

| 项目 | 旧 δ 方法 | 当前 J 方法 |
|---|---|---|
| 实际递推对象 | 二维局部标量系数 | 一维序列，但元素是 $z$ 多项式矩阵、行向量 |
| 第 $M$ 阶系数需要什么 | 有限二维对角线 | 第 $M$ 个 $t$ 系数在完整 $z$ 区间的积分 |
| 完成外层积分的手段 | 单项式矩的有限投影 | 固定 $t_0$ 的积分约化及边界项 |
| 预先需要待约化对象的 DE 吗 | 不需要 | 不需要 |
| 当前能否称为 $O(M)$ 标量展开 | 二维算法为 $O(M^2)$ | 不能 |

旧算法的固定多项式系数递推有 $O(M^2)$ 个标量格点；δ 投影也只访问这些格点。当前虽然只有 $O(M)$ 个 $t$ 系数，但在一般情况下第 $n$ 阶中间多项式的次数、分母极点阶数为 $O(n)$。固定窗口并不意味着每一步只有固定个标量操作。仅逐阶处理这些多项式便通常已有二次量级工作，完整积分约化和表示转换还会增加开销。

因此，“$J$ 只有一个自变量”保证的是外部系数数量线性，**没有证明目前这条系数生成路线的总成本线性**。教学例子说明线性标量递推可能存在；一般家族中如何在不预先完成 $J$ 约化的前提下得到同样有效的递推，需要单独解决。不能把生成完 $J$ 的 DE 后容易实现的线性展开，拿来作为当前前处理已经线性的论据。

## 10. 与实际代码的对应

| 步骤 | 实现位置 |
|---|---|
| 旧 FBI 二维系数递推 | [series_solver.tpp](../../expand/src/series_solver.tpp)：`getFBISeries`、`solveMasterCoeffX/Y` |
| 旧权重、Jacobian 和完整被积函数 | [integrand_expander.tpp](../../expand/src/integrand_expander.tpp) |
| 旧 FI/BFI/BBFI 的 δ 投影 | [delta_projector.tpp](../../expand/src/delta_projector.tpp)：`project` |
| 当前 FBI 连接和局部多项式输入准备 | [prepare.hpp](../native/prepare.hpp) |
| 当前 $F_n,V_n$ 多项式递推 | [series.cpp](../native/series.cpp)：主循环 |
| 当前固定点 $z$ 积分及核关系 | [series.cpp](../native/series.cpp)：`Reducer::reduce`、`complete_relations` |
| 由展开搜索 J 约化与 DE | [search.cpp](../native/search.cpp) |

本文描述的是当前 C++ 生产路线。旧 Python 有理函数实现只用于低阶交叉检查；`reconstruct` 是另一条验证路线，不是上述旧 δ 展开的前提。
