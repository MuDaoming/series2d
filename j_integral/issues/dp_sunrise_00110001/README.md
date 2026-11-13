# DP 三传播子、三 branch：等质量 on-shell sunrise

## 从 FI 定义独立家族

原 DP 中取 $D_3,D_4,D_8$，令 $L_1=\ell_1-p_2,L_2=\ell_2,p=p_4$，重编号为

$$
D_1=L_1^2-m^2+i0,\quad D_2=L_2^2-m^2+i0,\quad
D_3=(L_1+L_2+p)^2-m^2+i0,\qquad p^2=m^2.
$$

DP 的数值为 $m^2=p^2=1$。定义

$$
I_{\nu_1\nu_2\nu_3}^{(d)}=\int\frac{d^dL_1}{i\pi^{d/2}}
\frac{d^dL_2}{i\pi^{d/2}}D_1^{-\nu_1}D_2^{-\nu_2}D_3^{-\nu_3}.
$$

每个 branch 恰好一条传播子，因此 FBI 内部没有连续积分参数。记 $h=d/2,v=\nu_1+\nu_2+\nu_3$。在 $X_1+X_2+X_3=1$ 上，

$$
U=X_1X_2+X_1X_3+X_2X_3,
$$

$$
F=m^2U-p^2X_1X_2X_3
=m^2(X_1+X_2)(X_1+X_3)(X_2+X_3).
$$

FBI 在此仅为

$$
\mathcal I_\nu(X)=\frac{(-1)^v\Gamma(v-d)}{\prod_i\Gamma(\nu_i)}(F-i0)^{d-v}.
$$

取 $X_1=(1-z)t,X_2=(1-z)(1-t),X_3=z$，则

$$
J_\nu(t)=\frac{(-1)^v\Gamma(v-2h)}{\prod_i\Gamma(\nu_i)}
\int_0^1(1-z)X_1^{\nu_1-1}X_2^{\nu_2-1}z^{\nu_3-1}
U^{v-3h}(F-i0)^{2h-v}\,dz,
\qquad I_\nu=\int_0^1J_\nu(t)\,dt.
$$

## 完全解析完成 J 的积分

记 $a=t(1-t)$，则

$$
U=(1-z)[a+(1-a)z],
$$

$$
F=m^2a(1-z)
\left(1+\frac{t}{1-t}z\right)
\left(1+\frac{1-t}{t}z\right).
$$

所以被积函数恰好是三个线性因子的幂乘以 Euler 权重。这对应 Lauricella 超几何函数 $F_D^{(3)}$，不是数值积分的近似。其定义为

$$
F_D^{(3)}(A;B_0,B_1,B_2;C;x_0,x_1,x_2)
=\sum_{k_0,k_1,k_2\ge0}
\frac{(A)_{k_0+k_1+k_2}\prod_{j=0}^2(B_j)_{k_j}x_j^{k_j}}
{(C)_{k_0+k_1+k_2}\prod_{j=0}^2k_j!},
$$

先在 $|x_j|<1$ 定义，再按原始积分作解析延拓，$(A)_k=\Gamma(A+k)/\Gamma(A)$。令

$$
x_0=1-\frac1a,\quad x_1=-\frac{t}{1-t},\quad
x_2=-\frac{1-t}{t}.
$$

对任意正整数指标，有完整解析结果

$$
\boxed{\begin{aligned}
J_\nu(t)={}&\frac{(-1)^v\Gamma(v-2h)\Gamma(\nu_1+\nu_2-h)}
{\Gamma(\nu_1)\Gamma(\nu_2)\Gamma(v-h)}
(m^2)^{2h-v}t^{\nu_1-1}(1-t)^{\nu_2-1}a^{-h}\\
&\times F_D^{(3)}\left(\nu_3;
3h-v,v-2h,v-2h;v-h;x_0,x_1,x_2\right).
\end{aligned}}
$$

特别地，两个普通 J 主积分为

$$
\boxed{J_{111}=-\frac{\Gamma(3-2h)\Gamma(2-h)}{\Gamma(3-h)}
(m^2)^{2h-3}a^{-h}
F_D^{(3)}(1;3h-3,3-2h,3-2h;3-h;x_0,x_1,x_2),}
$$

$$
\boxed{J_{112}=\frac{\Gamma(4-2h)\Gamma(2-h)}{\Gamma(4-h)}
(m^2)^{2h-4}a^{-h}
F_D^{(3)}(2;3h-4,4-2h,4-2h;4-h;x_0,x_1,x_2).}
$$

这两个式子没有剩余积分，但仍是特殊函数表达式，尚未进一步化成更低阶超几何函数或 Gamma 函数。能写出解析式不等于两者存在有理函数系数的约化关系。

## 边界函数与 FI 的低 sector

第三个主积分是边界函数

$$
\boxed{B_{z=0;111}(t)=-\Gamma(3-d)(m^2)^{d-3}[t(1-t)]^{-h}.}
$$

其最后一重积分不为零。按 Beta 函数解析延拓，定义单圈 tadpole

$$
T(m^2)=\int\frac{d^dL}{i\pi^{d/2}}\frac1{L^2-m^2+i0}
=-\Gamma(1-h)(m^2)^{h-1},
$$

则

$$
\boxed{\int_0^1B(t)dt
=-\Gamma(3-2h)(m^2)^{2h-3}
\frac{\Gamma(1-h)^2}{\Gamma(2-2h)}
=\frac{d-2}{m^2}T(m^2)^2.}
$$

此等式可先在 $\Re h<1$ 的端点收敛区域验证，再作解析延拓。它说明这里的边界函数确实携带非零 FI 低 sector，不能删去。

原 Blade 清单显示的是三传播子 sector 有一个 FI 主积分 $I_{111}$；这不表示整个独立家族只有一个 FI 主积分。缩掉任一传播子得到两个 tadpole 的乘积。将等质量、传播子置换对称性纳入后，低 sector 还有 $T^2$。

现已从精确的 $z$ 全导数关系、端点项以及传播子置换对称性推出 $I_{112}$ 到 $I_{111},T^2$ 的解析约化，详见 [J 与 FI 的关系](analytic/J_and_FI_relations.md)。该推导不依赖特殊函数求值。

## 独立数值 IBP 和解析检验

输入见 [family_request.json](input/family_request.json)。运行

```sh
j_integral/src/IBP/build/j-ibp \
  --input j_integral/issues/dp_sunrise_00110001/input/family_request.json \
  --output j_integral/issues/dp_sunrise_00110001/analytic/numeric_ibp.json
python j_integral/issues/dp_sunrise_00110001/analytic/check.py
```

在 $p=1000003,d=26/7,t=17$ 的独立数值点得到 $B_{z=0;111},J_{111},J_{112}$，DE 闭合。输出状态明确是此点的候选主积分，单点闭合不等于任意指标完备性的证明。

解析检验用 $d=2.6$，此时两个原始 $z$ 积分均收敛。在 $t=0.37,0.5$ 比较原始积分和独立超几何级数计算。级数用 $z\mapsto1-z$ 的 Euler 变换，将三个参数变成 $x_j/(x_j-1)\in(0,1)$，再按总阶数求和，避免在收敛域外直接求原始级数。结果见 [checks.json](analytic/checks.json)。项目维度 $d=26/7$ 的两个积分也收敛；另在 $t=0.37$ 直接核验两式，用 $1-z=u^{14}$ 消除较强的可积端点奇性。
