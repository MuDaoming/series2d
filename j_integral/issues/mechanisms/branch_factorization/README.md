# 两个 branch 的内部一圈多项式为零：t 可以完全分离，J 与 FI 的约化等价

这是对 e6 全部 sector 的一个统一解析解释，比“换 branch 后计数下降”更强。

## 定义

两圈三个 branch 的圈动量依次为 $L_A=\ell_1,L_B=\ell_2,L_C=\ell_1+\ell_2$，各传播子为 $D_{bi}=(L_b+q_{bi})^2-m_{bi}^2+i0$。令

$$
\alpha_{bi}=X_by_{bi},\quad \sum_i y_{bi}=1,\quad
N_b=\sum_i\nu_{bi},\quad v=N_A+N_B+N_C,\quad h=d/2.
$$

定义 branch 内部的一圈参数多项式及平均平移：

$$
f_b(y)=\sum_i y_{bi}m_{bi}^2
-\sum_{i<j}y_{bi}y_{bj}(q_{bi}-q_{bj})^2,\qquad
P_b(y)=\sum_i y_{bi}q_{bi}.
$$

将圈动量平移后，完整 Symanzik 多项式为

$$
U=X_AX_B+X_BX_C+X_CX_A,
$$

$$
F=U(X_Af_A+X_Bf_B+X_Cf_C)-X_AX_BX_CK^2,
\qquad K=P_C-P_A-P_B.
$$

这是通常的完成平方恒等式，不假设图有传播子置换对称。

## 关键条件和分离公式

假设作为 y 的多项式恒等成立

$$
f_A=f_C=0.
$$

例如单个无质量传播子满足；两条无质量传播子之间的外动量无质量时也满足。取 **积掉 B** 的坐标

$$
X_A=(1-z)t,\quad X_C=(1-z)(1-t),\quad X_B=z.
$$

令 $a=t(1-t),L=a+(1-a)z$，则

$$
U=(1-z)L,\quad F=z(1-z)[Lf_B-a(1-z)K^2].
$$

作 projective 变换

$$
w=\frac zL,\quad D=1-(1-a)w,\quad
z=\frac{aw}{D},\quad1-z=\frac{1-w}{D},\quad dz=\frac a{D^2}dw.
$$

于是

$$
F=\frac{a^2w(1-w)}{D^3}[f_B-(1-w)K^2].
$$

物理 J 定义（含 FI/FBI 的原始 Gamma 归一化）为

$$
J_\nu(t)=\frac{(-1)^v\Gamma(v-2h)}{\prod_{bi}\Gamma(\nu_{bi})}
\int dy\prod_{bi}y_{bi}^{\nu_{bi}-1}
\int_0^1 dz\,(1-z)X_A^{N_A-1}X_C^{N_C-1}X_B^{N_B-1}
U^{v-3h}(F-i0)^{2h-v}.
$$

这里 $dy$ 包含每个 branch 的 simplex delta。变换后 D 的所有幂严格抵消，t 仅剩

$$
\boxed{J_\nu(t)=C_\nu\,t^{h-N_C-1}(1-t)^{h-N_A-1},}
$$

其中与 t 无关的 $C_\nu$ 明确为

$$
C_\nu=\frac{(-1)^v\Gamma(v-2h)}{\prod_{bi}\Gamma(\nu_{bi})}
\int dy\prod y_{bi}^{\nu_{bi}-1}\int_0^1dw\,
w^{2h-N_A-N_C-1}(1-w)^{N_A+N_C-h-1}
[f_B-(1-w)K^2-i0]^{2h-v}.
$$

再积 t，因此

$$
\boxed{J_\nu(t)=
\frac{t^{h-N_C-1}(1-t)^{h-N_A-1}}
{\mathrm B(h-N_C,h-N_A)}\,I_\nu.}
$$

公式先在适当收敛域成立，再按物理 i0/维数正规化解析延拓；不在 Beta 极点处单独截取某个因子。它适用于本问题无 rank、正整数有效指标、三个 branch 均保留的对象。

## 为什么这能说明主积分数

设 $P_\nu(t)$ 为上式 I 前的因子。任何 FI 约化

$$
I_\nu=\sum_\mu c_{\nu\mu}(d,\mathrm{kin})I_\mu
$$

立即给出固定 t 下的关系

$$
J_\nu=\sum_\mu c_{\nu\mu}\frac{P_\nu(t)}{P_\mu(t)}J_\mu.
$$

由于指标为整数，所有 Beta 比值可用 Gamma 递推化为 d 的有理函数，而 t 的幂次差也是整数。所以此系数仍在原来的有理函数域中。反向同理：除去明确的 P，J 的约化就是 FI 常数的约化。

因此在这些条件及相同非零子 sector 范围内，**这个坐标下的 J 与 FI 并非只碰巧数目相同，而是约化问题直接等价。**e6 中缺少某个 branch 的 sector 都含无标度 massless 单圈因子而为零，因此不会遗漏非零两-branch FI。

## 原来积掉 C 时，为什么会有额外 J

仍然假设 $f_A=f_C=0$，但采用旧坐标 $X_A=(1-z)t,X_B=(1-z)(1-t),X_C=z$。用同一个 w 变换，J 变成

$$
J_\nu^{\mathrm{old}}(t)=\mathcal N_\nu\,
 t^{N_A+N_C-h-1}(1-t)^{h-N_A-1}
\int dy\prod y^{\nu-1}\int_0^1dw\,
w^{N_C-1}(1-w)^{h-N_C-1}
[f_B-twK^2-i0]^{2h-v}.
$$

这里 $\mathcal N_\nu=(-1)^v\Gamma(v-2h)/\prod\Gamma(\nu_{bi})$。**t 留在核 $f_B-twK^2$ 内**，一般不能提出；不是上述纯幂函数乘以常数。原来 J 多出来的是这部分仍未完成的积分变换所保留的函数自由度。

可以直接看完整积分如何消掉它。令 $u=tw$，再交换 t、u 积分顺序，t 部分恰好是

$$
\int_u^1dt\,(t-u)^{\alpha-1}(1-t)^{\beta-1}t^{-\alpha-\beta}
=\mathrm B(\alpha,\beta)u^{-\beta}(1-u)^{\alpha+\beta-1},
$$

其中 $\alpha=h-N_C,\beta=h-N_A$。最后将 u 换成 $1-w$，就恢复前面的常数积分 C。也就是说，**完整 FI 可以解析完成这一重 Beta 积分，而固定 t 的旧 J 不能；改变切片则使这一步预先完成。**这给出了具体的积分机制，不只是一条未知系数的 t 全导数。

## e6 全家族适用，不只六传播子

e6 的 A 为 D1、D2，都是无质量，平移差 p3 满足 $p_3^2=0$，故 $f_A=0$。C 为 D6、D7，都是无质量，平移差 p4 满足 $p_4^2=0$，故 $f_C=0$。这对子 sector 的任意非空子集也成立。中间 B 为 D3、D4、D5。

所以积掉中间 branch 后，整个 e6 家族所有 J 的 t 依赖都精确分离，并非某个 sector 的偶然对称。没有传播子交换自对称的那些 sector 同样适用。

## 全局计数和实际验证

重新生成所有保留三个 branch、dot<=2 的 **1111 个目标**：

| 坐标/输入 | J 系统维数 | 边界函数 | DE |
| --- | ---: | ---: | --- |
| 原来积掉 C | 79 | 1 | 一般耦合 |
| 改为积掉 B | 61 | 0 | 全部为上面幂函数的对角 DE |

原 Blade 的 62 个 FI 代表中有一条显式重复：

$$
I_{\{1,0,1,0,0,1,0\}}=I_{\{0,1,0,0,1,0,1\}}.
$$

左边是外动量 p1 的 massless sunrise，右边平移圈动量后是外动量 -p2 的同一积分，且 $p_1^2=p_2^2=1$。因此合并这条重复后该清单为 61 个。该重复不是数值秩猜测，而是标量积分只依赖外动量平方的精确关系。

进一步直接把这 **61 个 Blade FI 代表对应的 J** 作为指定基底输入。程序确认它们独立并覆盖全部 1111 个目标，DE 闭合；在另一质数 $1000033$、t=19 复核。每个 DE 对角元均精确等于

$$
\frac{h-N_C-1}{t}-\frac{h-N_A-1}{1-t}.
$$

这里的“61”以当前 Blade 清单、显式重复合并和上述数值 IBP 完备范围为证据；一般的 J/FI 等价性来自解析公式，而非这些数值测试。

## 进一步验证：把 J 约化直接变成 FI 约化

`transfer_fi.py` 用 P 的有理比值把全部 1111 个 J 约化转成 FI 约化。同一质数、同一指定基底下，从 t=17 和 t=23 的两个独立 J 点得到的 **全部 FI 系数逐项相同**。

又针对六传播子 sector，用独立 Blade maximal-cut IBP 约化四个目标（三个 single-dot、一个 double-dot），与上述 FI 结果的顶层部分逐项相同。只比较了四个 cut 目标，不将其冒充全部 1111 个全局 FI 关系均与 Blade 独立核验。

FI 模 p 结果保存于 `e6_full/fi_dot2_modp.json.gz`，证据汇总见 `transfer_checks.json`。

## 文件与复现

- `check_theorem.py`：一般参数变换的符号检查、e6 两个零 branch 的检查，以及全部 DE 对角元核对。
- `check_numeric.py`：三个 branch 大小组合、两个 t 点，独立数值积分核验分离公式。
- `e6_full/request_dot2.json`/`j_dot2.json`：新的完整目标集合和自动基底结果。
- `request_original_dot2.json`/`j_original_dot2.json`：同一目标集合的旧坐标对照。
- `request_blade_basis*.json`/`j_blade_basis*.json`：指定 61 个 FI 代表对应的 J 基底，两个质数及两个同质数 t 点。
- `transfer_fi.py`：消掉解析 P 比值并核验 FI 系数的 t 独立性、已有 Blade 对照。

`e6_full/request.json`/`j_result.json` 是初次只用旧 78 个目标的诊断，得到 60 维，**不是本节主计数依据**。补齐所有 sector 的 dot<=2 后为 61；DE 闭合不足以保证初始少量目标已找全主积分，这个诊断特意保留。
