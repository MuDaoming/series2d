# 此四传播子家族的解析积分

## 对象与结论

沿用 [家族定义](../README.md) 的物理归一化，取

$$
D_1=\ell_1^2+i0,\quad D_2=\ell_2^2+i0,\quad
D_3=(\ell_2+Q)^2+i0,\quad D_4=(\ell_1+\ell_2+P)^2+i0,
$$

$$
I_\nu^{(d)}=\int\frac{d^d\ell_1}{i\pi^{d/2}}
\frac{d^d\ell_2}{i\pi^{d/2}}\prod_{i=1}^4D_i^{-\nu_i}.
$$

记 $s_P=P^2,s_Q=Q^2,s_R=(Q-P)^2,h=d/2,v=\sum_i\nu_i$。三个 branch 分别为 $\{D_1\},\{D_2,D_3\},\{D_4\}$。采用

$$
X_1=(1-z)t,\quad X_2=(1-z)(1-t),\quad X_3=z,
\qquad I_\nu=\int_0^1J_\nu(t)\,dt.
$$

这里六维基底是 $B=B_{z=0;1111}$、$J_{1011}$、$J_{1101}$、$J_{1111}$、$J_{1112}$、$J_{1121}$。前面三个已得到完整解析式；后面三个目前只完成了 $z$ 积分，仍有一个含 Gauss 超几何函数的积分。不能把后者称为完整闭式解。

以下幂与积分均采用由 Feynman $i0$ 给出的解析延拓。对于本例的正不变量，$(-s-i0)^\alpha=e^{-i\pi\alpha}s^\alpha$。

## 两个三传播子函数

令 $a=t(1-t)$，$L=a+(1-a)z$。三传播子无质量 sunrise 的 $U=(1-z)L$，$F=-sX_1X_2X_3$，因而

$$
J_s(t)=-\Gamma(3-d)(-s-i0)^{d-3}a^{d-3}
\int_0^1z^{d-3}(1-z)^{h-2}L^{3-3h}\,dz.
$$

作变换

$$
w=\frac zL,\quad z=\frac{aw}{1-(1-a)w},\quad
L=\frac a{1-(1-a)w},\quad
1-z=\frac{1-w}{1-(1-a)w},\quad
 dz=\frac{a\,dw}{[1-(1-a)w]^2}.
$$

所有 $1-(1-a)w$ 的幂抵消，积分成为 Euler Beta 积分。定义 $\mathrm B(x,y)=\Gamma(x)\Gamma(y)/\Gamma(x+y)$，则

$$
C_s=-\Gamma(3-d)(-s-i0)^{d-3}\mathrm B(d-2,h-1),
$$

$$
\boxed{J_{1011}=C_{s_R}[t(1-t)]^{h-2},\qquad
J_{1101}=C_{s_P}[t(1-t)]^{h-2}.}
$$

因此两者的比值为与 $t$ 无关的 $(s_P/s_R)^{d-3}$（本例正不变量）。这是固定 $d$、固定运动学下的解析常数比值，通常不属于约化所用的运动学和 $d$ 的有理函数系数域，不能据此直接把它们计为一个有理约化主积分。

再积掉 $t$ 得到熟知的 Gamma 函数形式，此处也直接由两个 Beta 积分推出：

$$
I_s=-\Gamma(3-d)(-s-i0)^{d-3}
\frac{\Gamma(h-1)^3}{\Gamma(3h-3)}.
$$

## 边界函数

$B$ 是 $J_{1111}$ 的 $z=0$ 系数函数，不是将整个 $J$ 的积分简单设为 $z=0$。其归一化为

$$
B(t)=\Gamma(4-d)X_2U^{4-3h}
\int_0^1[-UX_2s_Qy(1-y)-i0]^{d-4}\,dy
\quad (X_2=1-t,\ U=t(1-t)).
$$

于是

$$
C_B=\Gamma(4-d)(-s_Q-i0)^{d-4}\mathrm B(d-3,d-3),
\qquad \boxed{B(t)=C_Bt^{-h}(1-t)^{h-3}.}
$$

三个函数的对数导数为

$$
\frac{J_s'}{J_s}=(h-2)\frac{1-2t}{t(1-t)},\qquad
\frac{B'}B=\frac{h-3t}{t(t-1)}.
$$

它们均与保存的有限域 DE 完全一致。在 $d=26/7$，两个 $J_s$ 的幂为 $-1/7$，$B$ 的两端幂分别为 $-13/7,-8/7$。

边界函数还有一个有理函数原函数关系：

$$
g_B(t)=\frac{t(1-t)(2-h-t)}{(1-h)(2-h)},\qquad
\frac d{dt}[g_B(t)B(t)]=B(t).
$$

按独立端点指数的 Beta 函数解析延拓，

$$
\int_0^1B(t)\,dt=C_B\mathrm B(1-h,h-2)=0,
$$

因为分母为 $\Gamma(-1)$，而一般 $h$ 下分子有限。但此处没有普通积分的共同收敛带：$t=0$ 要求 $\Re h<1$，$t=1$ 要求 $\Re h>2$。因此这个零是指定解析正则化后的结果，不能用未经正则化的端点代入来证明。这也只解释了此边界函数的正则化周期；尚不能代替其余顶层关系的端点分析。

## 三个四传播子函数：解析完成一重积分

对所有 $\nu_i>0$，令 $y$ 为 branch 2 内 $D_3$ 的参数，$D_2$ 的参数为 $1-y$。定义

$$
q(y)=s_Qy(1-y),\qquad
H(y)=s_P-(s_P+s_Q-s_R)y+s_Qy^2.
$$

FBI 多项式及其归一化为

$$
F=-X_1X_2X_3H(y)-UX_2q(y)
=-(1-z)^2(1-t)[Lq(y)+tzH(y)],
$$

$$
\mathcal I_\nu^{(d)}(X)=\frac{(-1)^v\Gamma(v-d)}{\prod_i\Gamma(\nu_i)}
\int_0^1y^{\nu_3-1}(1-y)^{\nu_2-1}(F-i0)^{d-v}\,dy,
$$

$$
J_\nu(t)=\int_0^1(1-z)X_1^{\nu_1-1}X_2^{\nu_2+\nu_3-1}
X_3^{\nu_4-1}U^{v-3h}\mathcal I_\nu^{(d)}(X)\,dz.
$$

用同一个 $w=z/L$ 变换，$z$ 积分归结为

$$
\int_0^1w^{\nu_4-1}(1-w)^{h-\nu_4-1}
\left(1+\frac{tH(y)}{q(y)}w\right)^{d-v}dw
=\mathrm B(\nu_4,h-\nu_4)
{}_2F_1\left(v-d,\nu_4;h;-\frac{tH(y)}{q(y)}\right).
$$

最终

$$
\boxed{\begin{aligned}
J_\nu(t)={}&\frac{(-1)^v\Gamma(v-d)\Gamma(h-\nu_4)}
{\Gamma(\nu_1)\Gamma(\nu_2)\Gamma(\nu_3)\Gamma(h)}
t^{\nu_1+\nu_4-h-1}(1-t)^{h-\nu_1-1}\\
&\times\int_0^1 y^{\nu_3-1}(1-y)^{\nu_2-1}
[-q(y)-i0]^{d-v}
{}_2F_1\left(v-d,\nu_4;h;-\frac{tH(y)}{q(y)}\right)dy.
\end{aligned}}
$$

分别代入 $1111,1112,1121$ 即为三个顶层对象。此式减少了一重积分，但剩余 $y$ 积分尚未求出完整闭式。不将未经端点分析的 $t$ 幂级数当作完整解。

## 可复现检验

运行 `python check_formulas.py`（依赖 mpmath、sympy），输出 [checks.json](checks.json)。不调用 Blade 或远程作业。

- 精确符号检查 $U,F$ 分解、$w$ 变换、$g_B$ 原函数关系。
- 将三个对数导数与已有 DE 在质数 $1000003$ 下作多项式恒等检验。
- 50 位精度下，两个维度、两个 $t$ 点检验低 sector 的原始 $z$ 积分和边界的原始 $y$ 积分。
- 对三个顶层对象，在 $d=11/2$、$t=0.37,y=0.29$ 的收敛区域比较原始 $z$ 积分与超几何表达式；这是固定 $y$ 的内层积分检验，不是对完整顶层 $J$ 的数值检验。一般 $d$ 的公式通过解析延拓使用。
- 检查低 sector 完整 FI 的 Gamma 函数恒等式，并记录边界 Beta 正则化的趋零样本。
