# 用解析性 Bootstrap 处理第一重积分 J(t)

## 1. 从原始费曼积分到待约化对象

考虑无不可约分子的二圈积分

$$
I_{\boldsymbol\nu}^{d}
=\int\frac{d^d\ell_1}{i\pi^{d/2}}
       \frac{d^d\ell_2}{i\pi^{d/2}}
       \prod_{\alpha\in\mathcal A}D_\alpha^{-\nu_\alpha},
\qquad D_\alpha=q_\alpha^2-m_\alpha^2+i0.
$$

$d$ 是维度，$\ell_1,\ell_2$ 是圈动量，$q_\alpha$ 是由输入指定的圈动量与
外动量组合，$m_\alpha$ 是质量，$\nu_\alpha\geq0$，
$\mathcal A=\{\alpha:\nu_\alpha>0\}$。所有 branch 均活跃时，令

$$
v=\sum_\alpha\nu_\alpha,\qquad
\sigma_b=\sum_{\alpha\in\mathcal B_b}\nu_\alpha,
\qquad b=1,2,3.
$$

$\mathcal B_b$ 是按圈动量二次结构划分的三个传播子 branch。
将通常的费曼参数写成 $x_\alpha=X_by_\alpha$，其中
$\sum_bX_b=1$，每个 branch 内 $\sum_{\alpha\in\mathcal B_b}y_\alpha=1$。
以这些约束下的测度记为 $[dy]_{\boldsymbol\nu}$，定义 FBI

$$
\mathcal I_{\boldsymbol\nu}^{\Delta}(X)
=\frac{(-1)^v\Gamma(v-\Delta)}
        {\prod_{\alpha\in\mathcal A}\Gamma(\nu_\alpha)}
\int[dy]_{\boldsymbol\nu}
\frac{\prod_{\alpha\in\mathcal A}y_\alpha^{\nu_\alpha-1}}
     {[F(y;X)-i0]^{v-\Delta}}.
$$

$U,F$ 是上述替换后的第一、第二 Symanzik 多项式。二圈物理积分取
$\Delta=d$。沿用原来的坐标

$$
X_1=(1-z)t,\qquad X_2=(1-z)(1-t),\qquad X_3=z,
\qquad (t,z)\in[0,1]^2.
$$

令 $W_{\boldsymbol\nu}=\prod_bX_b^{\sigma_b-1}$，则包含全部权重的对象为

$$
G_{\boldsymbol\nu}(t,z)
=(1-z)W_{\boldsymbol\nu}(X)U^{v-3d/2}
  \mathcal I_{\boldsymbol\nu}^{d}(X),
$$

$$
J_{\boldsymbol\nu}(t)=\int_0^1G_{\boldsymbol\nu}(t,z)\,dz,
\qquad I_{\boldsymbol\nu}^{d}=\int_0^1J_{\boldsymbol\nu}(t)\,dt.
$$

端点按原程序的正则化约定处理。branch 缺失的原始 FI 需要额外退化边界处理，
不能把普通三 branch J 的结果称作所有 FI 的最终约化。
共同定义与有限域数值约化见
[全局 J 定义与算法](../../j_integral_global/docs/problem_solution.md)。

## 2. 新方法作用在哪一层

令 $M(t)$ 是独立的 J 主函数与必要边界函数组成的列向量。
所求约化和连接为

$$
J_\nu(t)=r_\nu(t)M(t),\qquad M'(t)=\mathsf D(t)M(t).
$$

固定维度、运动学与质数 $p$ 后，$r_\nu$、$\mathsf D$ 的元素属于
$\mathbb F_p(t)$。这次不再以长 Taylor 展开为主要重建手段。
现有 C++ 的一重积分有限域 IBP 提供某个 $t$ 点的有理约化系数值；
解析性 Bootstrap 先缩小这些有理函数的候选空间，再取样确定剩余常数。

论文与上游程序：https://arxiv.org/abs/2610.01154 ，
源码身份见 [upstream.md](upstream.md)。

## 3. 已知连接时，怎样约化新的 J

给定连接 $\mathsf D$、每个系数的分母包络 $d_j(t)$ 和分子次数上界 $b_j$，
原始候选为

$$
r_j(t)=\frac{\sum_{k=0}^{b_j}a_{jk}t^k}{d_j(t)}.
$$

上界对应允许的分母，实际系数可以约分。恒零分量用 $b_j=-1$ 表示，
必须有独立依据；不能因为少量点为零就排除。

对连接反复求导得到行

$$
\rho_i^{[0]}=e_i,\qquad
\rho_i^{[m+1]}=(\rho_i^{[m]})'+\rho_i^{[m]}\mathsf D.
$$

$e_i$ 是第 $i$ 个单位行向量。在满足正规解析部分的有限性条件时，
导数行给出允许的极点消去组合。上游 PL-Hasse 实现逐个不可约分母因子处理
有限阶导数源，消去辅助极点态，得到允许的分子模；随后与全部分子次数约束求交。
其输出可写作

$$
r(t)=\sum_{a=1}^{K}\kappa_a\,v_a(t),\qquad \kappa_a\in\mathbb F_p.
$$

$v_a$ 是已构造的有理行基，$K$ 是约束后的维数。
在普通点 $t=t_s$ 输入数值 IBP 行 $r(t_s)$，对 $\kappa_a$ 做线性拟合。
只有拟合秩达到 $K$ 才得到当前候选空间内的唯一结果；另外保留未参与拟合的点验证。

程序使用 $\tau=t-t_0$ 保存多项式，与旧结果的坐标一致。
因为 $d/d\tau=d/dt$，连接不额外乘 Jacobian；调用 C++ 时必须换回绝对坐标 $t$。

## 4. 尚无完整连接时，怎样逐层构造 DE

主函数按依赖关系排序，写成

$$
S'=B(t)S,\qquad T'=A(t)T+C(t)S.
$$

$S$ 是已处理的低层主函数列，$T$ 是当前层列，$B$ 是已经完成的低层连接，
$A$ 是当前层的对角连接，可由相应 cut 独立求得，$C$ 是尚未知的耦合块。
Bootstrap 的输入是 $A,B$ 与 $C$ 的分母及次数包络，**不输入 $C$ 的系数**。

在局部有限商模中，$A$ 的极点行构成源模，$B$ 给出低层允许分子模及其商。
未知 $C$ 必须使源模的关系在低层商模中成立，因此其不同行之间共享约束；
不能把每个未知条目分别做一遍普通插值。
上游 `BootstrapSpec(mode='de')` 与 `build_pl_ansatz` 实现这一构造。

确定当前层 $C$ 后，将已经解出的连接加入 $B$，再做下一层。
只保留当前层真正依赖的低层闭包，避免把无关 sector 放入每次模计算。

当前验证程序从参考 DE 的非零图提取强连通块并按依赖排序。
这验证的是逐层构造机制；从头计算 DP 时，应从 sector、边界及已知结构建立
依赖关系，并检查其闭包，不能把参考非零图当作已免费获得的物理输入。

## 5. t=1 处归一化：实际遇到的问题

44 维小型全局例子的标量块 $T=J_{\{1,0,1,1,1,0\}}$，在
$u=t-1$ 处具有

$$
B=\frac{R}{u}+O(1),\qquad
A=\frac{a}{u}+O(1),\qquad
C=\frac{c_{-2}}{u^2}+O(u^{-1}).
$$

如果正规低层解的初值 $s_0\in\ker R$ 满足 $c_{-2}s_0\ne0$，代入
$T=k/u+\cdots$ 得

$$
-k=ak+c_{-2}s_0,\qquad
k=-\frac{c_{-2}s_0}{1+a}.
$$

实际有限域残数中，这个配对非零。也就是说，低层的正规初值不一定能提升成
原归一化下的正规 $T$。直接套用这一归一化的候选模会排除真实的耦合系数。
这一局部证据并不等于证明特定物理边界值一定出现该极点。

对这个例子，使用统一标量规范变换

$$
\widetilde S=hS,\qquad \widetilde T=hT,\qquad h=t-1
$$

得到

$$
\widetilde B=B+\frac{h'}h\mathbf1,\qquad
\widetilde A=A+\frac{h'}h\mathbf1,\qquad
\widetilde C=C.
$$

未知耦合块保持不变，可以在归一化后的已知连接上构造其约束，再直接还原原对象。
此变换使 44 维例子的全部块通过；推广至其他 family 前仍需判断所需端点幂次。
不能把本例的一次幂当成所有二圈积分的定理。

将所有未知分母因子都乘入 $h$ 的对照也能通过，但本例中约束完全消失，
3509 个参数仍是 3509 个。这不是有效简化，所以保留为诊断对照。

## 6. 计算量与正确性标准

主要成本是按因子构造导数模、块扩张、次数约束求交，以及对 $K$ 个剩余参数拟合。
成本依赖主函数数目、稀疏耦合、极点阶、因子次数、源阶和 $K$；
不能仅因自变量是一个就称整个 Bootstrap 为线性复杂度。
收益应比较相同包络下的参数数目、取样点数，以及包含包络获取成本的总耗时。

当前测试区分：

1. 上游算术与接口测试。
2. 给定参考包络及对角块的受控重建实验。
3. 真正从头学习包络、确定 cut 对角块、构造未知全局 DE 的计算。

目前完成的是前两类。验证同时使用独立 C++ 数值取样、8 个留出点，
以及与旧结果的精确有理函数恒等式比较。上游有限导数源或有限窗口的状态
保留在输出中，经验通过不提升成特征零完备性证明。

## 7. 文件与运行

- `tests/build_kernels.py`：用本机 `/usr/local` 的 FLINT 编译原上游 C 内核。
- `tests/reference_space.py`：已知连接下的约化空间构造、拟合与精确核对。
- `tests/sector_de.py`：依赖闭包内逐块构造未知 DE 耦合，保存完整连接。
- `tests/point_oracle.py`：调用现有 C++ 数值 J-IBP，校验质数、坐标、主函数顺序与闭合。
- `tests/diagnose_block.py`：逐因子规范变换及残数诊断。

数据与结论见 [测试记录](results.md)。现阶段的包络来自参考结果，不应把这些
脚本的完成状态解读为已经完成 496 维 DP 全局系统的从头求解。
