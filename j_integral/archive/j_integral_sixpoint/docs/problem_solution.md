# 两圈六点例子的完整 J 连接

## 1. 原始积分和输入来源

采用 [2412.21054，第 IV 节图 5](https://arxiv.org/html/2412.21054v1) 的两圈六点例子。原作者程序在本机的 `Documents/DCT/gitlab/dct/example/e6_FBI.txt` 保存了传播子和数值点；原文件已复制到 `inputs/e6_FBI_original.txt`，路径和 SHA256 记在 `inputs/provenance.json`。不从图形猜测传播子，也不将浮点近似直接用于有限域。

定义

$$
I_{\nu}^{d}=\int\frac{d^d\ell_1}{i\pi^{d/2}}\frac{d^d\ell_2}{i\pi^{d/2}}\prod_{a=1}^7D_a^{-\nu_a},\qquad \nu_a\geq0.
$$

$d$ 是时空维度，$\ell_1,\ell_2$ 是圈动量，$p_1,\ldots,p_5$ 是输入的五个独立外动量，$p_6=p_1+p_2-p_3-p_4-p_5$。传播子的 $+i0$ 处方在下列简写中省略：

$$
\begin{aligned}
D_1&=\ell_1^2,&D_2&=(\ell_1+p_3)^2,\\
D_3&=\ell_2^2,&D_4&=(\ell_2+p_5)^2,\\
D_5&=(\ell_2+p_1+p_2-p_3-p_4)^2,\\
D_6&=(\ell_1+\ell_2+p_1)^2,&D_7&=(\ell_1+\ell_2+p_1-p_4)^2.
\end{aligned}
$$

全部内部质量为零。外腿满足 $p_1^2=p_2^2=1$，$p_3^2=\cdots=p_6^2=0$。从论文的精确 Mandelstam 值恢复原程序的标量积：

$$
\begin{gathered}
p_1\cdot p_2=3/2,\quad p_1\cdot p_3=-1/3,\quad p_1\cdot p_4=-3/8,\quad p_1\cdot p_5=-2/5,\\
p_2\cdot p_3=-5/11,\quad p_2\cdot p_4=-3/7,\quad p_2\cdot p_5=-6/13,\\
p_3\cdot p_4=1/30,\quad p_3\cdot p_5=1/34,\quad p_4\cdot p_5=-10242329/2042040.
\end{gathered}
$$

最后一项由 $p_6^2=0$ 求出；其余对角元由在壳条件给出。`family.json` 的 `gram` 就是 $p_i\cdot p_j$，不是 $2p_i\cdot p_j$。

## 2. FBI、J 和边界对象

三个 branch 为 $\mathcal B_1=\{1,2\}$、$\mathcal B_2=\{3,4,5\}$、$\mathcal B_3=\{6,7\}$。只在三个 branch 都有正指数时定义本次普通 J。令 $v=\sum_a\nu_a$、$\sigma_b=\sum_{a\in\mathcal B_b}\nu_a$，并将费曼参数写成 $x_a=X_by_a$；每个 branch 内 $\sum y_a=1$，$\sum_bX_b=1$。

沿用 [共同数学定义](../../j_integral_bootstrap/docs/problem_solution.md) 中包含 Gamma 归一化的 FBI $\mathcal I_{\nu}^{d}(X)$。$U$ 是第一 Symanzik 多项式，$F$ 是替换后的第二 Symanzik 多项式，FBI 的被积核为 $F^{d-v}\prod_a y_a^{\nu_a-1}$。使用

$$
X_1=(1-z)t,\quad X_2=(1-z)(1-t),\quad X_3=z,\qquad 0\leq t,z\leq1,
$$

并定义

$$
G_\nu(t,z)=(1-z)\prod_bX_b^{\sigma_b-1}U^{v-3d/2}\mathcal I_\nu^d(X),\qquad
J_\nu(t)=\int_0^1G_\nu(t,z)\,dz.
$$

于是原始三 branch 积分为 $I_\nu^d=\int_0^1J_\nu(t)\,dt$。本次只建立第一重积分后的连接，不将它称为原始 FI 全约化。

程序中 `Bz0{nu}` 是同一 FBI 主分量对应 $G_\nu$ 在 $z=0$ 的正则化有限项；`Bz1{nu}` 在当前端点处方下为零。边界独立性在物理解析关系及精确 corner 恒等式取商后判断，不能任意附加自由边界常数。

## 3. 实际测试集合和结果

固定 $p=1000003$、$d=26/7$，以 $\tau=t-2/5$ 保存有理系数。这里 $t$ 仍是符号变量，并未固定在某一个数值点。指数集合为所有三个 branch 活跃的 sector，且

$$
\operatorname{dot}(\nu)=\sum_{a:\nu_a>0}(\nu_a-1)\leq2.
$$

共有 63 个 sector、1111 个普通 J 目标。FBI 原始 corner 基有 43 个，去掉一条精确重复后保留 42 个；加入两个端点候选共 84 个，因此数值后端处理 1195 个对象。最后得到 79 个独立主对象：78 个普通 J 和 1 个边界函数。具体指数逐项列在 `inputs/basis.json`，不能将这 79 与论文在固定 branch 参数下的 FBI 数目 43 混为一谈。

全部 79 个主对象的导数都在该基内，完整连接

$$
\frac{dM(t)}{dt}=A(t)M(t),\qquad A\in\operatorname{Mat}_{79}(\mathbb F_p(t)),
$$

有 1164 个非零条目。第一次重构使用 $t=17,\ldots,112$ 共 96 点及内部保留点；少数条目先乘已有分母，再对剩余有理函数重构。没有使用现成的 J DE。

冻结连接后，在此前未用于重构的 $t=151,\ldots,158$ 八个完整数值点验证全部矩阵分量，包含隐含零：49928 次比较全部通过。结果见 `results/connection.json` 和 `results/verification_report.json`。这证明了本次固定维度、运动学、质数下的完整有理 $t$ 连接通过独立检验；尚不是一般 $d$ 或一般运动学的符号答案。

另用 $p=1000033$、$t=23$，将普通 J 测试集合提高到 dot≤3，并固定相同的 79 个主对象作约化和导数闭合检查，报告位于 `results/dot3_prime2.json`。这是额外的基底覆盖与闭合测试，不是将第一质数下的系数直接拿去比较另一质数。

## 4. 原包路线单独测试

`official_route/run_upstream.py` 不修改上游源码，调用原 `learn_anchor`、`prepare_setting`、`solve_de`。数值接口只返回当前 $t$ 点的 J-IBP 系数。原包的一次完整有理 anchor 是它明确要求的学习输入；它用来学习分母及次数，后续待求耦合系数不作为已知答案传入。

该路线显式允许共同归一化 $\widetilde M=[t(t-1)]^3M$，因此对角连接加 $3(1/t+1/(t-1))$；输出时按这一已知公式恢复物理 J 的连接。原包路线必须另外完成、比较，不能把已经完成的直接重构当成原包也已成功。

原包路线现已完成全部 47 个依赖块。其物理连接与独立重构在全部 6241 个矩阵条目上精确相同，见 `upstream_exact_comparison_v2.json`。随后又生成此前完全未使用的 $t=211,\ldots,218$ 八点，49928 个分量全部通过，见 `upstream_fresh_verification_report.json`。原包输出使用其 `SolvedDE.context` 中的真实主基顺序；适配层已加入非平凡 sector 重排回归，避免集合正确但行列标错。

79 维连接的分子、分母最高次数均为 50。唯一边界主函数是 $Bz0_{\{0,1,1,0,1,0,1\}}$；top sector 的五个普通 J 主函数为

$$
\begin{gathered}
J_{\{1,1,1,1,1,1,1\}},\quad J_{\{1,1,1,1,1,1,2\}},\quad J_{\{1,1,1,1,1,2,1\}},\\
J_{\{1,1,1,1,2,1,1\}},\quad J_{\{1,1,1,2,1,1,1\}}.
\end{gathered}
$$

第二质数的 dot≤3 测试覆盖 2852 个普通 J 目标，全部可在同一 79 维基中表示，导数未闭合数为零。
