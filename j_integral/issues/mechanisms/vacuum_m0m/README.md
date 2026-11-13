# 一般动量 IBP 与单尺度齐次性：顶层 FI 全部消掉

## 家族与实际计数

$$
D_1=L_1^2-m^2+i0,\quad D_2=L_2^2+i0,\quad
D_3=(L_1+L_2)^2-m^2+i0.
$$

定义 $I_{abc}=\int d^dL_1/(i\pi^{d/2})\,d^dL_2/(i\pi^{d/2})D_1^{-a}D_2^{-b}D_3^{-c}$。

它对应 DP 中质量 $(m,0,m)$、外动量平方为零的四个等价三传播子 sector：`10001010`、`10000101`、`01001010`、`00100101`。其完整 FI 仅依赖外动量平方，等于零动量的真空积分；参数多项式也直接给出同一结果。

独立 Blade maximal cut 返回顶层 FI **0 个**，独立 J IBP 返回顶层普通 J **1 个**，另有一个边界函数，DE 闭合。这里不是把跨 sector 主积分代表放在别处造成的假象：下面直接证明顶层 FI 约化到低 sector。

## 不使用传播子置换的精确推导

令 $T=\int d^dL/(i\pi^{d/2})(L^2-m^2+i0)^{-1}$，$T_2=\partial_{m^2}T$。动量 IBP

$$
0=\int\partial_{L_1}\cdot\frac{L_1}{D_1D_2D_3}
=(d-3)I_{111}-2m^2(I_{211}+I_{112})+TT_2.
$$

推导使用

$$
2L_1\cdot(L_1+L_2)=D_1+D_3-D_2+2m^2,
$$

以及缩掉 $D_1$ 后出现的无质量 tadpole 为零。另一方面，只有一个非零尺度，因而

$$
m^2\partial_{m^2}I_{111}=(d-3)I_{111},\qquad
\partial_{m^2}I_{111}=I_{211}+I_{112}.
$$

这里没有用 $I_{211}=I_{112}$。两式代入得到

$$
\boxed{I_{111}=\frac{TT_2}{d-3}
=\frac{d-2}{2(d-3)m^2}T^2.}
$$

因此本例不是 FI 恒等于零，而是顶层没有独立主积分。

## 为什么固定 t 的 J 不能照搬这个积分简化

用 Cheng-Wu 规范 $\alpha_1+\alpha_3=1$，令 $u=\alpha_1$、$x=\alpha_2\in[0,\infty)$。则

$$
U=u(1-u)+x,\quad F=m^2U,
$$

$$
I_{111}=-\Gamma(3-d)(m^2)^{d-3}
\int_0^1du\int_0^\infty dx\,[u(1-u)+x]^{-d/2}.
$$

$x$ 积分是一个初等幂积分，再积 $u$ 是 Beta 函数：

$$
I_{111}=-\frac{\Gamma(3-2h)}{h-1}(m^2)^{2h-3}
\mathrm B(2-h,2-h),\quad h=d/2,
$$

与 tadpole 约化完全相同。

但是 $t=\alpha_1/(\alpha_1+\alpha_2)=u/(u+x)$。固定 $t$ 后，$x=u(1-t)/t$，不再是可独立积分的 $[0,\infty)$ 变量。上面的初等积分简化不能在该切片内直接使用。事实上 $J_{111}(t)$ 是非零的 Appell 函数：

$$
J_{111}=-\frac{\Gamma(3-2h)}{2-h}(m^2)^{2h-3}
 t^{2h-3}[t(1-t)]^{-h}
F_1\left(1;h,3-2h;3-h;
1-\frac1{t(1-t)},-\frac{1-t}{t}\right).
$$

这具体说明了：完整积分可自由执行的单尺度 IBP/参数积分，在固定 $t$ 时会改变约束，因而不会变成 J 的点态零关系。

## 检验

`count.wl` 是独立 Blade maximal-cut 计数；`request.json`/`j_result.json` 是独立 J 输入输出。`check.py` 对 Gamma 结果和 tadpole 结果核对，并对两个 $t$ 点用原始 z 积分核对 Appell 表达式，见 `checks.json`。所有表达式按维数正规化解析延拓；完整真空参数积分可先在 $2<\Re d<3$ 的收敛区域计算。
