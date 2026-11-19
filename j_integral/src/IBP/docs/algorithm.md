# 当前J(t)层数值IBP怎样做

本文说明现有程序，不是新的AB播种实现。完整输入格式、归一化及运行方法见 [接口README](../README.md)。

z方向局部Hermite的完整数学步骤、终止条件和结构复杂度见 [z_hermite_algorithm.md](z_hermite_algorithm.md)。

## 1. 从FI到待约化对象

两个圈动量为 $\ell_1,\ell_2$，传播子

$$
D_i=(c_{i1}\ell_1+c_{i2}\ell_2+q_i)^2-m_i^2+i0,
\qquad
I_\nu^d=\int\frac{d^d\ell_1}{i\pi^{d/2}}
\frac{d^d\ell_2}{i\pi^{d/2}}\prod_{i:\nu_i>0}D_i^{-\nu_i}.
$$

$q_i$ 是输入外动量的线性组合；现有FBI接口只接受非负指标，不处理不可约分子。圈动量二次结构相同的传播子组成branch $\mathcal B_b$，$b=1,2,3$。活跃传播子的参数写为 $x_i=X_by_i$，各branch中 $\sum y_i=1$，$\sigma_b=\sum_{i\in\mathcal B_b}\nu_i$，$v=\sum_i\nu_i$。

配方 $\sum_iX_by_iD_i=\ell^{\mathsf T}A\ell+2\ell\cdot B+C$，定义 $U=\det A$、$F=B^{\mathsf T}\operatorname{adj}(A)B-UC$。FBI定义为

$$
\mathcal I_\nu^d(X)=\frac{(-1)^v\Gamma(v-d)}{\prod_{i:\nu_i>0}\Gamma(\nu_i)}
\int[dy]_\nu\frac{\prod_{i:\nu_i>0} y_i^{\nu_i-1}}{(F-i0)^{v-d}}.
$$

$[dy]_\nu$ 包含各活跃branch的simplex约束。采用

$$
X_1=(1-z)t,\quad X_2=(1-z)(1-t),\quad X_3=z,
$$

$$
G_\nu(t,z)=(1-z)\prod_bX_b^{\sigma_b-1}
U^{v-3d/2}\mathcal I_\nu^d(X),\qquad
J_\nu(t)=\int_0^1G_\nu(t,z)\,dz,
\quad I_\nu=\int_0^1J_\nu(t)\,dt.
$$

J的普通目标保留三个branch。缺branch的完整FI可能是非零的低sector或子圈乘积，不能一律扔掉。边界对象的定义亦不等于普通缺branch的J。

## 2. 单次计算需要什么

输入是完整FI家族、branch选择、端点约定、目标指标，以及数值 $t=a,d,p$。可提供固定J/边界主对象列表；不提供时在候选目标中选基。只做约化设 `derivatives=[]`，需要DE时额外请求主对象的一阶导数。

**不需要预先知道J的约化或J的DE。** 当前在一个 $t$ 点完成有限域约化；后续沿 $t$ 的有理重构、AB学习及最后的 $t$ 积分是其它步骤。

## 3. FBI内部：已知递推，不是原始动量空间的大范围播种

程序从家族生成FBI各sector的约束矩阵，使用FBI已知的升降指标及求导公式处理目标，保留低sector。小矩阵消元用于这些闭式递推的系数和主对象构造。

以包括完整权重的FBI主对象组成列向量 $H(t,z)$。程序同时得到

$$
G_\nu(t,z)=r_\nu(t,z)H(t,z),\qquad
\partial_zH=A_z(t,z)H,
$$

以及需要时的 $\partial_tH=A_tH$。这里 $r_\nu$ 是有理行向量，$A_z$ 已按sector排成三角形式。已知的是FBI层求导公式，不是事先完成了J层约化。

固定 $t=a$ 后，$r_\nu,A_z$ 属于 $\mathbb F_p(z)$，不会为了每个展开阶保留不断变复杂的 $t$ 有理函数。实现入口为 [fbi.hpp](../kernel/fbi.hpp) 的 `reduce`、`model`。

## 4. z积分：保留Hermite优势

对任意有理行 $g(z)$，定义

$$
\nabla_zg=g'+gA_z,
\qquad \partial_z(gH)=(\nabla_zg)H.
$$

Hermite逐个三角分量，将待积行约化为

$$
r=h+\nabla_zg.
$$

有限因子 $f$ 的高阶极点通过 $\operatorname{res}_f(A_{ii})-(k-1)f'$ 的局部除法降阶；非共振时不产生大范围全局未知量。无穷远的多项式层也按标量递推消去。相应系数为零的共振方向保留，不能强除；随后单独补充关系。

对每次标量消元，代码直接检查有理恒等式 $r-h-g'-gA_{ii}=0$。三角非对角项产生的下sector源项继续约化。当前 `atomic` 模式按分量、不可约因子、极点阶和分子次数缓存基本项的Hermite像，复用相同局部计算。

实现见 [z_integrator.hpp](../kernel/z_integrator.hpp) 的 `ScalarHermite`、`TriangularIntegrator`。这是当前主要优势，不能用一个巨大的统一分母Laporta ansatz替换后仍声称保留了它。

## 5. 端点与补充关系

采用声明的 $z=1$ 正规化零边界；$z=0$ 由三角FBI连接的Laurent-log局部解计算有限部分。对全导数，端点项保留：

$$
\int_0^1(\nabla_zg)H\,dz
=-\operatorname{FP}_{z=0}(gH).
$$

所以Hermite原函数产生的 $z=0$ 边界必须进入约化，不能擅自设为零。`fbi_analytic` 模式再加入由FBI原始积分定义给出的解析端点阶数条件、无log条件和已知角点等价关系。

`complete` 对标量Hermite剩下的共振原函数产生全导数关系，`Quotient` 用稀疏有限域消元存储独立关系。这不是对所有高dot/rank对象统一撒种子。有限的 `completion`、`germ_order` 是显式截断参数，截断稳定或数值闭合不自动成为任意指标完备性定理。

## 6. 目标约化及可选DE

将全部目标和边界变为上述商空间中的稀疏坐标，选出独立J/边界对象，解出目标的主对象系数。固定基请求会检验其独立性及是否张成请求目标。

若要求导数，使用FBI的 $t$ 求导公式或first-jet，将导数变成新的有理被积行，再通过同一个z约化器。由此建立J的DE。只求约化时不重复计算这一步；这也是规模测试统一使用 `derivatives=[]` 的原因。

程序调用和分阶段诊断见 [point.cpp](../kernel/point.cpp)、[main.cpp](../main.cpp)。最后一重 $t$ 积分的IBP和FI约化不属于本模块。

## 7. 怎样解释“方程数、变量数”

- FBI矩阵次数：局部闭式递推所需的小矩阵操作，不是完整IBP调用次数。
- z补充方程：提交给 `Quotient` 的非零行；变量为剩余极点积分坐标和端点常数，独立行数为该消元器的秩。
- 目标表示矩阵：目标/边界的坐标，用来选基和求系数，不额外算成IBP方程。
- 输出主对象数：上述约化后留下的J/边界对象数，不是中间变量数。

e6和DP的具体统计与比较限制见 [规模总结](../../../tests/IBP_SIZE_SUMMARY.md)。

## 与拟议AB改进的接口

新思路是对有限奇点的整数支要求从 $\mu=0$ 开始，不允许负整数幂；非整数支保留，暂不约束无穷远。应按AB构造允许分子组成的模，用其生成元参数化候选被积行或关系，减少未知参数，再考虑如何接入上述Hermite流程。

这不是把单独FBI分量的整数幂删掉。正式后端尚未接入，独立e6实验已得到可保持约化系数的模参数化和Hermite像空间方案，见 [实验结果](../../../issues/ab_seed_constraints/README.md)。矩阵组合之间可能发生解析极点相消；需要由整个组合的模条件判断。现有FBI的 $z$ 连接可作为这一步AB的输入，无需先有J的DE。新约束对物理密度、全导数原函数和边界项的具体适用方式需要分别证明和验证。
