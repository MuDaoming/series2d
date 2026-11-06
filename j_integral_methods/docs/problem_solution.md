# 一重 J 积分的直接约化与 DE：对象、四条路线和实现

## 1. 问题与完整对象定义

研究目标是先处理二圈 FBI 表示中的一重积分 J，再由 J 的约化与 DE 处理最后一重积分，得到原始费曼积分 FI 的约化。本试验始终沿用原 t,z 参数化。

原始无分子二圈积分为

$$
I_{\boldsymbol\nu}^{d}=\int\frac{d^d\ell_1}{i\pi^{d/2}}
\frac{d^d\ell_2}{i\pi^{d/2}}
\prod_{\alpha\in\mathcal A}D_\alpha^{-\nu_\alpha},
\qquad D_\alpha=q_\alpha^2-m_\alpha^2+i0.
$$

$d$ 是时空维度，$\ell_1,\ell_2$ 是圈动量，$q_\alpha$ 是给定的圈动量与外动量线性组合，$m_\alpha$ 是质量。指标 $\nu_\alpha\ge0$，$\mathcal A=\{\alpha:\nu_\alpha>0\}$。将传播子按圈动量二次结构分为三个非空 branch $\mathcal B_b$，定义

$$
v=\sum_\alpha\nu_\alpha,\qquad
\sigma_b=\sum_{\alpha\in\mathcal B_b}\nu_\alpha,\qquad
W_{\boldsymbol\nu}(X)=\prod_{b=1}^{3}X_b^{\sigma_b-1}.
$$

在费曼参数 $x_\alpha=X_by_\alpha$ 下，$\sum_bX_b=1$，每个 branch 内 $\sum y_\alpha=1$。令 $U,F$ 为相应的 Symanzik 多项式，$[dy]_{\boldsymbol\nu}$ 是三个 branch 内带这些 Dirac 约束的参数测度。固定分支积分定义为

$$
\mathcal I_{\boldsymbol\nu}^{\Delta}(X)
=\frac{(-1)^v\Gamma(v-\Delta)}{\prod_{\alpha\in\mathcal A}\Gamma(\nu_\alpha)}
\int[dy]_{\boldsymbol\nu}
\frac{\prod_{\alpha\in\mathcal A}y_\alpha^{\nu_\alpha-1}}
{[F(y;X)-i0]^{v-\Delta}}.
$$

$\Delta$ 是 FBI 的维度参数；物理二圈积分取 $\Delta=d$。原始 FI 的 FBI 表示及 branch 内测度的完整定义见 [共同定义](../../j_integral/docs/problem_solution.md) 第 2—3 节。

取

$$
X_1=(1-z)t,\quad X_2=(1-z)(1-t),\quad X_3=z,
\qquad t,z\in[0,1].
$$

定义包含归一化和 Jacobian 的完整被积函数

$$
G_{\boldsymbol\nu}(t,z)=(1-z)W_{\boldsymbol\nu}(X)
U(X)^{v-3d/2}\mathcal I_{\boldsymbol\nu}^{d}(X).
$$

于是

$$
\boxed{J_{\boldsymbol\nu}(t)=\int_0^1G_{\boldsymbol\nu}(t,z)dz,
\qquad I_{\boldsymbol\nu}^{d}=\int_0^1J_{\boldsymbol\nu}(t)dt.}
$$

本试验的直接目标是找到 J 之间的有理函数约化，并构造所选主对象向量 $\mathbf J$ 的一阶系统 $\partial_t\mathbf J=A_J(t)\mathbf J$。不将 $A_J$ 当成已知输入。

## 2. 所有直接方法共同使用的已知 FBI 输入

从已知 FBI 约化、移维和求导公式出发，选定并重定义一个 MFBI 列向量 $\mathbf f(t,z)$，将 Jacobian、branch 权重及共同非整数 $U$ 幂一并纳入规范选择。MFBI 是 master FBI，与 J 主积分不同。得到

$$
G_{\boldsymbol\nu}=r_{\boldsymbol\nu}\mathbf f,
\qquad \partial_z\mathbf f=A_z\mathbf f,
\qquad \partial_t\mathbf f=A_t\mathbf f.
$$

$r_{\boldsymbol\nu}$ 是有理行向量，$A_z,A_t$ 是有理矩阵，向量长度为 $r$。具体规范重定义见 [展开推导](../../j_integral/docs/delta_vs_t_expansion.md) 第 4 节。这些式子属于 FBI 层，直接程序没有读入旧 J 约化或旧 J DE。

直接对目标求一次导数：

$$
\partial_tG_{\boldsymbol\nu}
=s_{\boldsymbol\nu}\mathbf f,
\qquad s_{\boldsymbol\nu}=\partial_t r_{\boldsymbol\nu}+r_{\boldsymbol\nu}A_t.
$$

所以若能完成有理行向量的 z 积分约化，就可以同时处理 J 与 J'；不必生成第 2 阶至第 M 阶的 t 导数或展开。

## 3. 核心恒等式：带端点的 z 方向 IBP

对有理行向量 $a(t,z)$ 定义

$$
\nabla_z a=\partial_z a+aA_z.
$$

于是

$$
\partial_z(a\mathbf f)=(\nabla_z a)\mathbf f,
$$

$$
\int_0^1(\nabla_z a)\mathbf f\,dz
=[a\mathbf f]_{z=0}^{z=1}.
$$

若目标行向量满足

$$
r_{\nu}-\sum_jc_j(t)r_{\mu_j}=\nabla_z a,
$$

就有

$$
J_\nu=\sum_jc_j(t)J_{\mu_j}+[a\mathbf f]_0^1.
$$

这是直接约化的基础。将左侧 $r_\nu$ 换成 $s_\nu$，就是直接构造 DE 的同一个问题。

为了正确保留边界，计算空间同时包含完整积分和正则端点的 $f_i(t,\epsilon)$，$\epsilon=0,1$。对这些对象取上述全导数关系生成的商空间。这可理解为带边界的积分约化；不能只对有理函数模全导数后将端点丢弃。

本轮输入沿用原测试约定：z=0 正则显式保留，z=1 为已声明的解析零边界。这个约定是物理输入；程序验证有理恒等式不等于自动证明该解析边界条件。正规化、积分与微分先在适用的收敛区域解释，再解析延拓。

## 4. 路线 A：现有 Taylor 展开搜索，作为基准

在一般点 $t_0$ 展开所有 J，逐阶完成 z 积分，再搜索多项式关系和 DE。旧程序并不需要预先得到 J DE；它已能利用 FBI 关系产生系数。

本次不修改这个算法，仅将其作为正确性和时间基准。它为获得 J 约化与 DE 生成了大量高阶系数，其分母极点阶数和中间多项式次数随展开阶数增长。完整推导见 [delta_vs_t_expansion](../../j_integral/docs/delta_vs_t_expansion.md)。

## 5. 路线 B：取样点的结构化降极点，然后重建

### 5.1 在 t=a 直接处理 J 与 J'

这里 a 是有限域中的一般取样点。专门化后令

$$
\mathbf f_a(z)=\mathbf f(a,z),\qquad
\partial_z\mathbf f_a=\frac{B_a(z)}{q_a(z)}\mathbf f_a.
$$

$q_a$ 是包含目标、导数目标及连接所需极点的平方自由多项式；$B_a$ 为多项式矩阵。当前共享降极点实现要求有限极点为单极点、无穷远没有多项式连接部分。这个限制沿用现有引擎，不能因此宣称已支持所有可能的 FBI 规范。

结构化降极点将每个有理行向量化为

$$
r_\nu(a,z)=h_\nu(z)+\nabla_z g_\nu(z).
$$

保留余项的积分及端点值为有限坐标 $C_\lambda(a)$：

$$
J_\nu(a)=\sum_\lambda v_{\nu\lambda}(a)C_\lambda(a).
$$

每次降极点都有精确有理恒等式检查。奇异降极点映射的核所产生的低层关系也纳入消元，避免先前 DB、DP 的过完备计数。对 $s_\nu(a,z)$ 同样处理，得到 J' 的坐标。

不求出 $C_\lambda(a)$ 的物理数值；它们是形式积分坐标。有限域计算只涉及有理代数。

### 5.2 先选基，再直接求约化和 DE 的值

将低 dot 对象和简单边界优先入基，在线性独立的 J/边界候选中选出主对象。在每个新点重新建立这些主对象的坐标矩阵，解出

$$
J_\nu(a)=\sum_jc_{\nu j}(a)J_{\mu_j}(a),
\qquad
J'_{\mu_i}(a)=\sum_j A_{J,ij}(a)J_{\mu_j}(a).
$$

如果所选基不能张成目标或一阶导数，程序不会将其默认为闭合。偶然奇异点会被跳过并记录，超过设定次数则报错。

这里的坐标基 $C_\lambda(a)$ 可以随取样点变化；每个点先在同一个局部坐标空间求完线性代数，重建的是最终 c 和 A_J，而不是强行重建不同点之间未对齐的余项坐标。

### 5.3 重建有理系数

设某个最终系数是 $c(t)=n(t)/d(t)$，对一组互异取样点 a_i 得到 c(a_i)。程序以 $\tau=t-t_0$ 为输出坐标，先插值得到多项式 p，令

$$
M(\tau)=\prod_i(\tau-(a_i-t_0)),
\qquad nd^{-1}\equiv p\pmod M.
$$

利用多项式扩展欧几里得算法重建 n,d。对于训练点数 T，当前尝试 $\deg n,\deg d\le\lfloor(T-1)/2\rfloor$；若所有系数尚未通过留出的 8 点，就继续增加样本。没有把 Taylor 关系次数上限预先当成重建结果。

之后还在另外 8 个未参与重建的有理点直接检查积分证书。重建与有限个新点验证是强的计算证据，但不自动成为一般 t 的符号恒等式证明；路线 D 提供了独立的符号 t 检查。

这条路线没有 Taylor 系数，输出中的 `Taylor_order` 为 0。取样点数由最终系数的有理函数次数驱动，不由人为要求的展开到 M 阶驱动。

## 6. 路线 C：取样点直接设总导数 ansatz

它与路线 B 使用同样的重建框架，但不逐级降极点，而是直接在线性空间中消去 IBP。

在固定点 a 写 $A_z=B/q$，选择共同分母 $q^K$，覆盖全部目标和一阶导数。以

$$
\mathcal P_{ij}=\int_0^1\frac{z^j}{q^K}f_{a,i}(z)dz,
\qquad 0\le j\le N,
$$

及正则端点值为原始坐标。对每个分量 i 和适用的整数 j，设

$$
g=\frac{z^j}{q^{K-1}}e_i^{\mathsf T}.
$$

$e_i$ 是分量单位列向量。清分母后的行向量为

$$
q^K\nabla_z g
=q\partial_z(z^je_i^{\mathsf T})
-(K-1)q'z^je_i^{\mathsf T}
+z^je_i^{\mathsf T}B.
$$

积分后减去 $[g\mathbf f_a]_0^1$，得到精确的 IBP 线性关系。程序对这些关系做一次行消元，再将全部目标和一阶导数投影到商空间。未知量的系数与积分变量 z 无关，属于当前有限域。

默认选择 $N=K\deg q+L-1$，L=16，并扩大到足以覆盖目标的最高分子次数。测试另用 L=6、K 再增加 1，检查基与全部约化、DE 是否变化。有限 ansatz 的稳定性不构成对任意家族完备性的证明。

路线 C 的好处是构造直观、无需高阶极点的逐级算法；代价是矩阵维数随 K、q 次数及 FBI 秩增长。路线 B 利用了极点结构，通常更适合扩大规模。

## 7. 路线 D：保留符号 t，直接在 F_p(t) 上做 IBP

将路线 C 的系数域由固定点的 $\mathbb F_p$ 换成 $\mathbb F_p(t)$，并在整个过程中保留符号 t：

$$
A_z(t,z)=B(t,z)/q(t,z).
$$

此时 q 可以取成关于 z 的首一多项式，其系数属于 $\mathbb F_p(t)$。同样设 $g=z^j e_i^{\mathsf T}/q^{K-1}$，清分母、比较 z 的幂次，并在 $\mathbb F_p(t)$ 上直接消元。最后对全部 J 和 J' 解出有理函数系数。

这条路线既不展开 t，也不在 t 上取样。所得约化与 DE 来自符号 t 的精确 IBP 线性代数，因此与路线 B、C 的插值重建构成独立检查。在给定有限域运动学和边界约定下，它们是该符号有理函数空间内的恒等关系；仍不能把运动学有限域取样当成一般运动学符号证明。

当前符号程序的 ansatz 默认 L=4。它会检查每个解在所建商空间内的余项为零；输出保留基、约化和 DE，但暂不逐条输出展开后的总导数生成函数 g。若需要独立审阅每条全 t 的 g，可进一步加入证书导出。

其主要成本是消元中有理函数的分子分母增长。与纯有限域线性代数相比，即使矩阵不大，也可能在 t 有理函数运算上花较多时间。因此是否优于重建，需要实际比较。

## 8. 如何解释复杂度和主积分计数

直接方法的主要规模是 FBI 秩 r、z 极点多项式次数、目标 dot 导致的极点阶数 K、有限 IBP 空间维数，以及最终关于 t 的有理函数次数。它不再需要一个高阶 Taylor 截断 M。

本轮 DP、DB cut 的 FBI 秩均为 1；真空测试的保留 FBI 秩也是 1。矩阵版本的方向及非零端点另由一个二分量连接的精确 IBP 测试检验。不能据此宣称大型多 MFBI 的完整未 cut 家族已达到同样性能。

数主对象时总是包括所需边界函数。DP 的八个是七个 J 加一个边界；它不是原始 FI 的十五个 master。不同路线沿用同一对象排序规则，以便直接比较全部有理约化与 DE。

DE 闭合、提高 ansatz 范围、两有限域和独立点验证共同检验当前基的稳定性。本轮目标是 dot≤2 的集合；没有将它写成任意高 dot 完备性的普遍证明。

## 9. 程序位置

- [direct_j.cpp](../native/direct_j.cpp)：取样点、两种 z-IBP 后端、选基、直接一阶导数约化及有理重建。
- [symbolic_j.cpp](../native/symbolic_j.cpp)：符号 t 的直接总导数 ansatz。
- [tests.cpp](../tests/tests.cpp)：有理重建、留出点拒绝、基求解、矩阵连接和两个非零端点。
- [run_experiments.py](../tests/run_experiments.py)：串行启动 C++，与展开法比较全部系数，记录计时并运行第二有限域、独立点和 ansatz 范围检查。

传播子和运动学来自原 `j_integral/examples/*.json`。共享 FBI 构造器、结构化降极点器及独立点检查器从原目录编译使用，不复制一套拓扑专用程序，不修改旧结果。
