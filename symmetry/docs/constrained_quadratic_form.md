# FBI Branch 约束下二次型的对称规范表示

## 1. 问题

FBI 的参数分母写成

$$
F(y;X)=\frac12 y^T R(X)y.
$$

参数按照 branch 分组：

$$
y=
\begin{pmatrix}
y_1\\
\vdots\\
y_B
\end{pmatrix},
\qquad
y_b\in K^{n_b},
$$

并满足每个 branch 独立的仿射约束

$$
\mathbf 1_b^T y_b=1,
\qquad b=1,\ldots,B,
$$

其中

$$
\mathbf 1_b=(1,\ldots,1)^T\in K^{n_b}.
$$

因此 FBI 真正依赖的是 $F$ 在约束面

$$
\mathcal A=
\left\{
y\ \middle|\ 
\mathbf 1_b^T y_b=1,\ b=1,\ldots,B
\right\}
$$

上的限制，而不是矩阵 $R$ 在整个参数空间中的具体表示。

两个不同矩阵 $R$ 和 $R'$ 可能满足

$$
\frac12y^TRy
=
\frac12y^TR'y
\qquad
\text{对所有 }y\in\mathcal A,
$$

尽管

$$
R\ne R'.
$$

例如，向 $F$ 加入

$$
\sum_{b=1}^B
\left(\mathbf 1_b^Ty_b-1\right)H_b(y)
$$

不会改变其在 $\mathcal A$ 上的值，却通常会改变表面的二次型矩阵。

所以，直接比较

$$
P^TRP=R'
$$

是过强的。需要构造一个只依赖于 $F|_{\mathcal A}$、同时不选择具体消元参数的表示。

## 2. Branch 重心与零和子空间

### 2.1 Branch 重心

对第 $b$ 个 branch，定义重心向量

$$
c_b=\frac1{n_b}\mathbf 1_b.
$$

显然

$$
\mathbf 1_b^T c_b=1.
$$

将所有 branch 重心组合为

$$
c=
\begin{pmatrix}
c_1\\
\vdots\\
c_B
\end{pmatrix}.
$$

$c$ 是约束面 $\mathcal A$ 上一个不依赖传播子顺序的特殊点。

### 2.2 零和子空间

定义第 $b$ 个 branch 的零和子空间

$$
V_b=
\left\{
u_b\in K^{n_b}\ \middle|\ 
\mathbf 1_b^Tu_b=0
\right\}.
$$

对应投影矩阵为

$$
H_b=
I_{n_b}
-\frac1{n_b}\mathbf 1_b\mathbf 1_b^T.
$$

它满足

$$
H_b^2=H_b,
\qquad
H_b^T=H_b,
\qquad
H_b\mathbf 1_b=0.
$$

将各 branch 投影组合成分块对角矩阵

$$
H=\operatorname{diag}(H_1,\ldots,H_B).
$$

对应的总零和子空间为

$$
V=
\bigoplus_{b=1}^B V_b
=
\operatorname{Im}(H).
$$

## 3. 约束面参数化

任意 $y\in\mathcal A$ 可以唯一写成

$$
y=c+u,
\qquad u\in V.
$$

证明如下。

如果 $y\in\mathcal A$，定义

$$
u=y-c.
$$

对每个 branch，

$$
\mathbf 1_b^Tu_b
=
\mathbf 1_b^Ty_b-\mathbf 1_b^Tc_b
=
1-1=0.
$$

所以 $u\in V$。

反过来，如果 $u\in V$，则

$$
\mathbf 1_b^T(c_b+u_b)=1+0=1,
$$

所以 $c+u\in\mathcal A$。

由于 $u\in\operatorname{Im}(H)$，还有

$$
Hu=u.
$$

## 4. 二次型在约束面上的表示

把

$$
y=c+u
$$

代入

$$
F(y)=\frac12y^TRy.
$$

得到

$$
\begin{aligned}
F(c+u)
&=
\frac12(c+u)^TR(c+u)\\
&=
\frac12c^TRc
+u^TRc
+\frac12u^TRu.
\end{aligned}
$$

因为 $Hu=u$ 且 $H^T=H$，

$$
u^TRc=u^THRc,
$$

并且

$$
u^TRu=u^THRHu.
$$

定义

$$
\boxed{
Q=HRH
}
$$

$$
\boxed{
L=HRc
}
$$

$$
\boxed{
F_0=\frac12c^TRc
}
$$

则约束面上的二次型可以写成

$$
\boxed{
F(c+u)
=
F_0+u^TL+\frac12u^TQu,
\qquad u\in V.
}
$$

因此，三元组

$$
\boxed{
\mathcal R(R)=(Q,L,F_0)
}
$$

给出了 $F$ 在 branch 约束面上的表示。

## 5. 等价判据

设

$$
F(y)=\frac12y^TRy,
\qquad
F'(y)=\frac12y^TR'y.
$$

令

$$
\Delta R=R-R'.
$$

两者在约束面上相等，当且仅当

$$
\frac12y^T\Delta R\,y=0
\qquad
\text{对所有 }y\in\mathcal A.
$$

代入 $y=c+u$：

$$
\frac12(c+u)^T\Delta R(c+u)
=
\Delta F_0+u^T\Delta L+\frac12u^T\Delta Q\,u,
$$

其中

$$
\Delta Q=H\Delta R H,
$$

$$
\Delta L=H\Delta R c,
$$

$$
\Delta F_0=\frac12c^T\Delta R c.
$$

如果

$$
\Delta Q=0,\qquad
\Delta L=0,\qquad
\Delta F_0=0,
$$

那么显然两个二次型在 $\mathcal A$ 上相等。

反过来，如果上式对所有 $u\in V$ 为零，则这是 $u$ 上的恒零二次多项式，所以其二次、一次和常数部分分别为零。因此

$$
\boxed{
F|_{\mathcal A}=F'|_{\mathcal A}
\iff
\begin{cases}
HRH=HR'H,\\
HRc=HR'c,\\
c^TRc=c^TR'c.
\end{cases}
}
$$

等价地，

$$
\boxed{
F|_{\mathcal A}=F'|_{\mathcal A}
\iff
\mathcal R(R)=\mathcal R(R').
}
$$

## 6. 为什么不需要选择被消去的参数

常见做法是在每个 branch 中消去一个变量，例如

$$
y_{b,n_b}
=
1-\sum_{i=1}^{n_b-1}y_{b,i}.
$$

这种做法本身正确，但需要指定“最后一个”参数。传播子置换可能把被消去参数映射到未被消去参数，使比较时出现额外仿射变量变换。

三元组 $(Q,L,F_0)$ 不选择任何具体参数作为 dependent variable。

- $c_b$ 对 branch 中所有参数完全对称；
- $H_b$ 只区分重心方向和零和方向；
- 整个约束面由 $c+V$ 表示；
- 任意 branch 内参数置换都保持该分解。

因此它适合用于参数置换 symmetry 的 canonicalization。

## 7. 传播子置换下的协变性

设 $P$ 是保持 branch 划分的传播子置换矩阵。对每个 branch，$P$ 只重排该 branch 内的参数。

由于

$$
P^T\mathbf 1_b=\mathbf 1_b,
$$

所以

$$
P^Tc=c.
$$

同时

$$
P^THP=H.
$$

令

$$
R'=P^TRP.
$$

则

$$
\begin{aligned}
Q'
&=HR'H\\
&=HP^TRPH\\
&=P^THRHP\\
&=P^TQP,
\end{aligned}
$$

以及

$$
\begin{aligned}
L'
&=HR'c\\
&=HP^TRPc\\
&=P^THRc\\
&=P^TL.
\end{aligned}
$$

常数部分满足

$$
F_0'
=
\frac12c^TR'c
=
\frac12c^TP^TRPc
=
\frac12c^TRc
=
F_0.
$$

所以

$$
\boxed{
(Q,L,F_0)
\xrightarrow{P}
(P^TQP,\;P^TL,\;F_0).
}
$$

这是标准的矩阵、向量和标量置换作用，适合继续使用 Pak 型规范化。

## 8. Branch 置换

如果还允许 branch 置换 $\tau$，则传播子置换矩阵 $P$ 会把整个 branch block 映射到另一个 branch block。

设 source branch $b$ 有 $n_b$ 个活跃参数。只有在

$$
n_b=n_{\tau(b)}
$$

时才能建立 branch 参数双射。

在重新排列后的 canonical branch 顺序中重新构造

$$
c_\tau,\qquad H_\tau.
$$

或者等价地使用整体置换矩阵 $P$：

$$
c_\tau=P^Tc,
\qquad
H_\tau=P^THP.
$$

同时 branch 参数 $X$ 按 $\tau$ 变换。在当前项目的二维坐标中，这对应

$$
(X,Y)\mapsto g_\tau(X,Y).
$$

因此两个 sector 的正确比较对象是：

$$
\mathcal R_s(X,Y)
=
(Q_s(X,Y),L_s(X,Y),F_{0,s}(X,Y)),
$$

以及

$$
\mathcal R_t(g_\tau(X,Y)).
$$

## 9. 用于 Pak 规范化

对固定 branch 变换，传播子置换 $P$ 对三元组的作用为

$$
Q\mapsto P^TQP,
$$

$$
L\mapsto P^TL,
$$

$$
F_0\mapsto F_0.
$$

可以构造 canonical 序列，例如

$$
\mathcal S(Q,L,F_0)
=
\left(
F_0,\,
L_1,Q_{11},\,
L_2,Q_{21},Q_{22},\,
\ldots,\,
L_n,Q_{n1},\ldots,Q_{nn}
\right).
$$

Pak 算法逐个确定传播子顺序：

1. 尝试一个尚未使用且 branch 相容的传播子；
2. 把对应的 $L_i$、$Q_{ij}$ 和 $Q_{ii}$ 加入当前前缀；
3. 比较完整前缀；
4. 只保留字典序最小的并列候选；
5. 重复直到所有传播子位置确定。

由于 $(Q,L,F_0)$ 在传播子置换下按普通矩阵和向量方式变换，原先针对矩阵 $R$ 的 Pak 前缀逻辑可以直接推广。

差别只是 canonical token 从

$$
R_{ij}
$$

改为

$$
F_0,\quad L_i,\quad Q_{ij}.
$$

## 10. 最终 symmetry 判据

设 source sector 为 $s$，target sector 为 $t$，传播子映射为 $P_\sigma$，branch 置换诱导二维坐标变换 $g_\tau$。

分别计算

$$
\mathcal R_s(X,Y)
=
(Q_s,L_s,F_{0,s})
$$

和

$$
\mathcal R_t(g_\tau(X,Y))
=
(Q_t^\tau,L_t^\tau,F_{0,t}^\tau).
$$

则正确的 FBI sector symmetry 判据为

$$
\boxed{
P_\sigma^TQ_sP_\sigma=Q_t^\tau,
}
$$

$$
\boxed{
P_\sigma^TL_s=L_t^\tau,
}
$$

$$
\boxed{
F_{0,s}=F_{0,t}^\tau.
}
$$

这等价于

$$
\frac12y^TR_sy
=
\frac12y'^TR_ty'
$$

在所有 branch 参数约束成立时相等。

它不要求

$$
P_\sigma^TR_sP_\sigma=R_t^\tau.
$$

因此可以识别那些 $R$ 表面不同、但对应同一个约束后 FBI 分母的 symmetry。

## 11. 系数域限制

重心和投影矩阵需要使用

$$
\frac1{n_b}.
$$

所以系数域 $K$ 中必须满足

$$
n_b\ne0.
$$

在特征零的有理函数域中总是成立。

在有限域 $\mathbb F_p$ 中，需要要求

$$
p\nmid n_b
$$

对所有 branch 成立。当前项目使用的大素数远大于 branch 大小，因此通常没有问题。

## 12. 小例子

考虑单个 branch，有两个参数：

$$
y_1+y_2=1.
$$

此时

$$
c=
\frac12
\begin{pmatrix}
1\\1
\end{pmatrix},
\qquad
H=
\frac12
\begin{pmatrix}
1&-1\\
-1&1
\end{pmatrix}.
$$

设

$$
R=
\begin{pmatrix}
a&b\\
b&d
\end{pmatrix}.
$$

则

$$
F_0=\frac18(a+2b+d),
$$

$$
L=
\frac14
\begin{pmatrix}
a-d\\
d-a
\end{pmatrix},
$$

而 $Q=HRH$ 只保留沿 $y_1-y_2$ 方向的二次信息。

交换 $y_1\leftrightarrow y_2$ 后，

$$
a\leftrightarrow d,
$$

$Q$ 做同步行列置换，$L$ 的两个分量交换，$F_0$ 不变。

整个表示不需要提前决定消去 $y_1$ 还是 $y_2$。

## 13. 实现结论

当前 symmetry 模块中直接比较 $R$ 的做法应替换为：

```text
R_s
  ↓
根据 sector 的活跃 branch 构造 c 和 H
  ↓
计算 Q = H R_s H
       L = H R_s c
       F0 = 1/2 c^T R_s c
  ↓
对 (Q,L,F0) 做 Pak canonicalization
  ↓
按 canonical key 建立 sector orbit
  ↓
用 (Q,L,F0) 的三个精确等式验证 mapping
```

该方法：

- 不选择被消去的 $y$；
- 与每个 branch 的约束完全等价；
- 对传播子置换协变；
- 可以继续使用 Pak 前缀剪枝；
- 不需要枚举全部 $\sigma$。
