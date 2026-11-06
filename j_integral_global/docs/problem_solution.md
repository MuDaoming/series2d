# 原 t,z 参数化下，DP 第一重 J 的全局约化

## 1. 目标与对象定义

本轮处理第一重积分 J 的约化和关于 t 的 DE，保留低 sector 耦合。暂不进行最后一次 t 积分，也不启动原始 FI 的完整 Blade 约化。

原始无分子二圈积分为

$$
I_{\boldsymbol\nu}^{d}=\int\frac{d^d\ell_1}{i\pi^{d/2}}
\frac{d^d\ell_2}{i\pi^{d/2}}
\prod_{\alpha\in\mathcal A}D_\alpha^{-\nu_\alpha},
\qquad D_\alpha=q_\alpha^2-m_\alpha^2+i0.
$$

$d$ 是时空维度，$\ell_i$ 是圈动量，$q_\alpha$ 是指定的圈动量与外动量组合，$m_\alpha$ 是质量，$\nu_\alpha\ge0$，$\mathcal A=\{\alpha:\nu_\alpha>0\}$。传播子按圈动量二次结构划分为三个 branch $\mathcal B_b$。定义

$$
v=\sum_\alpha\nu_\alpha,\qquad
\sigma_b=\sum_{\alpha\in\mathcal B_b}\nu_\alpha,\qquad
W_{\boldsymbol\nu}(X)=\prod_{b=1}^3X_b^{\sigma_b-1}.
$$

令通常的费曼参数 $x_\alpha=X_by_\alpha$，其中 $\sum_bX_b=1$，每个 branch 内 $\sum y_\alpha=1$。$U,F$ 为此替换后的 Symanzik 多项式，$[dy]_{\boldsymbol\nu}$ 表示 branch 内带这些 Dirac 约束的测度。FBI 定义为

$$
\mathcal I_{\boldsymbol\nu}^{\Delta}(X)
=\frac{(-1)^v\Gamma(v-\Delta)}{\prod_{\alpha\in\mathcal A}\Gamma(\nu_\alpha)}
\int[dy]_{\boldsymbol\nu}
\frac{\prod_{\alpha\in\mathcal A}y_\alpha^{\nu_\alpha-1}}
{[F(y;X)-i0]^{v-\Delta}}.
$$

二圈物理积分取 $\Delta=d$。测度和 FBI 矩阵 R 的完整定义见 [共同定义](../../j_integral/docs/problem_solution.md)。本轮始终采用

$$
X_1=(1-z)t,\qquad X_2=(1-z)(1-t),\qquad X_3=z,
\qquad t,z\in[0,1].
$$

包含 Jacobian 和全部归一化的对象为

$$
G_{\boldsymbol\nu}(t,z)=(1-z)W_{\boldsymbol\nu}(X)
U^{v-3d/2}\mathcal I_{\boldsymbol\nu}^d(X),
$$

$$
\boxed{J_{\boldsymbol\nu}(t)=\int_0^1G_{\boldsymbol\nu}(t,z)dz,
\qquad I_{\boldsymbol\nu}^d=\int_0^1J_{\boldsymbol\nu}(t)dt.}
$$

程序从传播子、Gram 矩阵和质量构造 FBI，不读入旧 J 约化或旧 J DE。

## 2. “全局”的范围

DP 三个 branch 的传播子数为 $(3,3,2)$。三个 branch 均有活跃传播子的 sector 共

$$
(2^3-1)(2^3-1)(2^2-1)=147.
$$

本轮普通 J 的目标为这些 sector 中的全部 dot≤2 对象：活跃指标至少为 1，且 $\sum_{\alpha\in\mathcal A}(\nu_\alpha-1)\le2$，合计 2952 个 J。

全局 FBI 求导、IBP 和 J 约化必须保留低 sector 项。不能分别求 maximal cut 后将其直接拼成全局约化。本轮处理的是整个耦合系统。

若整个 branch 缺失，上述普通三 branch J 定义不再直接适用。原始 FI 可以是非零因子化积分，但需要外层退化边界或分布表示。本轮不把“147 个普通 J sector 的约化”称作“所有原始 FI sector 的最终约化”。

## 3. 已知 FBI 系统与必要的恒等对象合并

选取规范化 MFBI 向量 $\mathbf f$，将共同非整数 U 幂、branch 权重与 Jacobian 纳入规范。MFBI 是 master FBI，与 J 主对象不同。得到

$$
G_\nu=r_\nu\mathbf f,\qquad
\partial_z\mathbf f=A_z\mathbf f,\qquad
\partial_t\mathbf f=A_t\mathbf f.
$$

原按 sector 选角点的方法保留了 135 个 MFBI，但不同 sector 的角点 FBI 可能完全相同。全局时必须先合并这些对象。只按 sector 保留代表会产生冗余表示，甚至使两种求导方式的形式系数不一致，差别实际上落在尚未合并的零关系上。

本实现对相同 branch 活跃数的 corner，枚举 branch 内置换，消去每个 branch 的一个 y。将约束面上的二次型写成

$$
F=c+L^{\mathsf T}u+\frac12u^{\mathsf T}Hu.
$$

$u$ 是剩余的 branch 内独立参数。精确比较 c、L、H 关于 t,z 的有理函数，而不是只在一个 t 点比较。相同 branch 数保证测度、Gamma 因子及规范权重也相同。这与原 symmetry/sector-map 方法的目的相同，当前实现使用约束二次型进行核验。

DP 合并了 20 个重复角点，MFBI 数为 115。对整个系统另行检查

$$
\partial_zA_t-\partial_tA_z+A_tA_z-A_zA_t=0.
$$

第一有限域、t=17 的全部 6670 个上三角矩阵元检查通过。此处相容性是必要的输入检查，不将其冒充 J 的 DE。

## 4. 固定 t 的 FBI 和一次导数

为了避免构造全局符号 t 有理函数，先记录 R、U 和 branch 参数的一次 t 导数，再在一般点 $t=a\in\mathbb F_p$ 专门化，然后在 $\mathbb F_p(z)$ 上构造 FBI 关系。

目标导数满足

$$
\partial_tG_\nu=(\partial_t r_\nu+r_\nu A_t)\mathbf f.
$$

不能对已经取值的系数直接求导而得到零。实现有两种可比较的导数路径：原 FBI 求导恒等式提高指标并移维；以及用双数环 $\mathbb F_p(z)[\tau]/(\tau^2)$ 对 FBI 约化系数自动求一次导数，$\tau=t-a$。后者没有生成 J 的高阶 Taylor 级数。

双数有理函数统一编码成分母与 $\tau$ 无关的形式。若 n、d 的零、一阶项为 $n_0,n_1,d_0,d_1$，则

$$
\frac{n_0+\tau n_1}{d_0+\tau d_1}
=\frac{n_0d_0+\tau(n_1d_0-n_0d_1)}{d_0^2}
\pmod{\tau^2}.
$$

这避免了双数分母表示造成的不必要增长。两种导数方式在合并正确的 FBI 恒等对象后精确一致，测试包括低 sector 和 dot3。

## 5. 利用 sector 三角结构完成完整 z-IBP

按高 sector 到低 sector 排列 MFBI，A_z 为上三角。对有理行向量 g 定义

$$
\nabla_zg=\partial_zg+gA_z,
\qquad \int_0^1(\nabla_zg)\mathbf f\,dz=[g\mathbf f]_0^1.
$$

若当前目标行为 v，在第 i 个分量求标量 Hermite 分解

$$
v_i=h_i+g_i'+g_iA_{z,ii},
$$

则必须向所有低 sector 转移

$$
v_j\longleftarrow v_j-g_iA_{z,ij},\qquad j>i,
$$

并保留 $[g_if_i]_0^1$。逐分量执行后，得到整个耦合系统的精确积分证书，而不是对角块 cut 的约化。

全局公共分母次数曾达到 361；直接套用“MFBI 秩×公共分母次数”的稠密降极点矩阵会很大。当前按不可约因子分别降极点，再用三角结构传递低层贡献。

对于因子 f 上的 k 阶极点，取 $g_i=a/f^{k-1}$。最高层消元映射是模 f 的

$$
\bar a\bigl(\operatorname{res}_f A_{z,ii}-(k-1)f'\bigr).
$$

非共振时在相应有限域扩张上求逆；共振时保留障碍项，并用核方向的完整总导数生成额外关系。无穷远也处理 $z^{-1}$ 层，即常数原函数；不能只消去多项式部分。每个标量步骤检查原始有理恒等式。

核方向不能按独立 cut 的零关系直接删除。其低 sector 贡献同样经完整三角系统传播后，才加入全局关系空间。

## 6. 全局端点：先区分一般 FBI 解和物理 FBI

全局 A_z 的部分非对角元在 z=0 有一阶或二阶极点，因此不能将所有端点分量视作互相独立的正则数。

实现先构造上三角系统在 z=0 的 Laurent-log 局部解，表达式形如

$$
f_i(z)=\sum_{n,l}c_{i,n,l}z^n(\log z)^l.
$$

系数 $c_{i,n,l}$ 是有限组局部初值的线性组合。对原函数与 FBI 的乘积，端点的正则化值取 $z^0(\log z)^0$ 系数。该定义在真实端点极限存在时等于普通极限。日志与负幂不是事先随意删除，而是先由耦合 ODE 表达，再施加有依据的物理条件。

第二个端点 z=1 沿用已声明的解析零条件。它不是通过有限域数值大小推断的，也不是“见到极点就置零”。

仅把一般局部解用于全局 J，得到 596 维的较大闭合表示。物理 FBI 还满足以下额外端点条件，必须纳入。

## 7. DP 的物理 z=0 条件及其依据

### 7.1 解析性与消失阶数

本例不是凭经验假设 FBI 解析。DP 给定运动学在 z=0 时有

$$
F(y;t,0)=U(t,0)\left[t\mathcal F_1+(1-t)\mathcal F_2\right],
\qquad U(t,0)=t(1-t),
$$

$$
\mathcal F_1=1-y_2y_3\ge\frac34,
\qquad
\mathcal F_2=y_4^2-\frac14y_5y_6\ge-\frac1{16}.
$$

这里原传播子编号为 1—8，branch 内 y 非负且和为 1；删除传播子等价于相应 y=0，故同一界对全部子 sector 有效。因此

$$
F(y;t,0)\ge U(t,0)\frac{13t-1}{16}>0,
\qquad \frac1{13}<t<1.
$$

所以 FBI 的紧致 y 积分在这一开区域内对 z 解析；所得恒等关系再解析延拓。Gamma 因子不改变 z 的解析性。对 corner MFBI，第三个 branch 活跃数为 $n_{3,i}$，完整规范化函数满足

$$
\boxed{f_i(t,z)=z^{n_{3,i}-1}\times\text{在 z=0 解析的函数}.}
$$

因此其负幂、log 系数，以及低于 $n_{3,i}-1$ 的 Taylor 系数都为零。程序通过 `z0_fbi_analytic` 显式启用这些条件；不将这个 DP 的证明无条件推广到别的家族，也禁止把它套在固定传播子 cut 的一般 FBI 解上。

### 7.2 相同物理端点的合并

当 $n_{3,i}=1$ 时，z=0 的 F 不再依赖第三个 branch 具体选了哪条传播子。对于相同 branch 活跃数、且端点约束二次型完全相同的两个 corner，Gamma、测度及规范权重相同，端点函数必须相等。程序仍用保留符号 t 的精确约束二次型证明这个等式，而不是从某个数值点猜测。

这些条件把一般局部初值的 115 个方向限制为 15 个物理边界方向。在第一有限域，J 候选表示由 596 降为 496；DE 保持闭合。

原 maximal-cut 中独立保留的边界对象，不一定是全局物理 FBI 的独立端点。原 cut 八维结果不能直接作为全局 topsector 的计数。

## 8. 有限域选基、约化和有理重构

在固定 t 点，完整 z-IBP 后的余项积分和物理端点是形式坐标；没有计算物理积分数值。以低 dot 边界优先、随后低 dot J 的方式选基，解出所有目标和主对象导数的系数。

本轮 DP 第一组得到 481 个 J 加 15 个边界函数，合计 496 个候选主对象。对这些对象的导数闭合是必要检查；有限域点上闭合和计数稳定仍不能代替完整 t 重构。

对多个 t 点的约化及 DE 系数进行有理重构，训练点与 8 个留出点分开。之后再在新的有理点重建完整未 cut FBI 系统，直接检验系数。原 cut 专用的单 MFBI 检查器不能冒充这一多 MFBI 验证。

数值 worker 一次只计算一个 t 点，共享这一点的昂贵 FBI 准备与完整耦合消元。不同点由外部动态池调度，使用独立输出、文件锁、退出码及闭合检查；输入和执行文件冻结，已验证点可以恢复使用。

## 9. 当前实现边界

- 普通 J 目标为三个 branch 都活跃的无 rank 积分，包含全部低 sector 耦合。
- 物理端点条件已有本 DP 运动学的解析依据；推广家族需重新证明。
- 本轮候选积分为 dot≤2，核关系及端点局部展开有显式有限范围。DE 闭合与稳定性检验不能当作任意家族任意 dot 的完备性定理。
- 496 是当前候选基大小；尚未完成 DP 的全部 t 有理重构时不称为最终主积分数或完整已验收的全局约化。
- 更大的生产取样拟采用 Tianhe 并行；当前 SSH 入口不可达，未提交远程作业。具体状态见 [testing](testing.md)。
