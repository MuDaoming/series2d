# 从 J 的约化和 DE 得到原始 FI：等质量 vac 的 IBP 验证

## 1. 目标与结果

本测试直接处理最后一重 $t$ 积分，不再引入辅助变形参数。

在三个质量平方均为 1、$d=7/13$、$p=1000003$ 的 2loopvac 上，对全部

$$
\nu_i\ge1,\qquad\sum_{i=1}^{3}(\nu_i-1)\le3
$$

的 20 个 FI 进行约化。只读取已经算出的 J 约化和 DE，搜索 $t$ 全导数证书，结果约到 $I_{111}$ 与 $I_{110}$。

20 个目标的有理函数 IBP 证书全部通过精确有限域检查；与 Blade 从原始圈动量传播子独立计算的结果模 p 逐项比较，20/20 一致，无差异。这是完整家族的验证，没有把下 sector 在 topsector cut 上设为零。

## 2. 对象、归一化和边界

以 $I_{\boldsymbol\nu}=\operatorname{FI}_{\boldsymbol\nu}^{d}$ 表示原始 FI。采用前文相同的 $G$ 和 $J$ 定义，令

$$
h=d/2,\quad a=t(1-t),\quad L=a+(1-a)z,\quad\kappa=-\Gamma(3-d).
$$

等质量 vac corner 为

$$
G_{111}(t,z)=\kappa(1-z)^{1-h}L^{-h},\quad
J(t)=\int_0^1G_{111}(t,z)dz,\quad
B(t)=G_{111}(t,0)=\kappa a^{-h}.
$$

上一阶段已获得并独立验证

$$
\boldsymbol{J}(t)=\begin{pmatrix}J(t)\\B(t)\end{pmatrix},\qquad
\boldsymbol{J}'=A(t)\boldsymbol{J},
$$

$$
A(t)=\begin{pmatrix}
\dfrac{(1-2t)(2-h)}{1-a}&-\dfrac{1-2t}{1-a}\\
0&-\dfrac{h(1-2t)}a
\end{pmatrix}.
$$

定义两种积分 period

$$
P_J=\int_0^1J(t)dt=I_{111},\qquad P_B=\int_0^1B(t)dt.
$$

后者不是额外的未知 FI 主积分，其归一化由 Beta/Gamma 恒等式直接确定：

$$
P_B=\kappa\frac{\Gamma(1-h)^2}{\Gamma(2-2h)}
=(d-2)I_{110},\qquad I_{110}=\Gamma(1-h)^2.
$$

剩余 $t$ 端点必须单独处理。对 $\operatorname{Re}h<0$ 的公共收敛域，

$$
B(0)=B(1)=0,
$$

但

$$
\boxed{J(0)=J(1)=\kappa\operatorname B(1-h,2-h)=\tfrac12P_B\ne0.}
$$

再按维度解析延拓。有限域上的 $d=7/13$ 不用来判断这些极限的大小。
因此不能把所有 $t$ 全导数边界项一律丢弃；也不能因为上一重 z 积分的某个端点采用解析零边界，就认为 J 的 t 端点也为零。

## 3. 有理 ansatz 的 IBP

从 J 约化读取每个目标的有理行

$$
J_{\boldsymbol\nu}=\boldsymbol r_{\boldsymbol\nu}(t)\boldsymbol J(t).
$$

对每个目标，同时求有限域常数行 $\boldsymbol c_{\boldsymbol\nu}$ 和有理行 $\boldsymbol g_{\boldsymbol\nu}(t)$，使

$$
\boxed{\boldsymbol r_{\boldsymbol\nu}
=\boldsymbol c_{\boldsymbol\nu}+\boldsymbol g_{\boldsymbol\nu}'
+\boldsymbol g_{\boldsymbol\nu}A.}
$$

本轮取 $\tau=t-2/5$，从已知 DE 选择 $D(\tau)=1-a(t)$（允许有限域整体常数归一化），设每个 g 分量的分母为 $D^2$、分子次数不超过 12。通分后比较多项式系数，解有限域线性系统。所有参数由输入给出，C++ 引擎中没有 vac 的质量或运动学。

解出后直接在有理函数域验证上述恒等式，不以有限阶系数通过替代积分证书。随后计算

$$
I_{\boldsymbol\nu}=c_JP_J+c_BP_B
+\big[g_JJ+g_BB\big]_0^1
=c_JP_J+\left(c_B+\frac{g_J(1)-g_J(0)}2\right)P_B.
$$

本轮证书在两个 t 端点都正则，$D(0)=D(1)=1$，所以不需要奇异证书的额外极限算法。
最终用 $P_B=(d-2)I_{110}$ 转回原始 FI 基。结果文件分别保存全导数证书、边界贡献、积分后的 period 系数和物理 FI 系数。

需要强调：两个 period 由上一阶段的 J 基给出；本轮验证的是约化算法和边界映射。没有声称当前程序已经为任意家族自动发现最小 FI 基或自动求出端点映射。

## 4. Blade 独立对照

`tools/check_vac_blade.wl` 从原始传播子

$$
\ell_1^2-1,\quad\ell_2^2-1,\quad(\ell_1+\ell_2)^2-1
$$

调用本地 Blade，完整家族不加 cut，目标为相同的 20 个指标，优选 $I_{111},I_{110}$。Blade 返回两个 masters。

Blade 的内部 master 标签先依据优选积分的实际约化转换回物理积分；未消去的内部对象或剩余表达式不能计为成功。Blade 的精确有理系数再模 $p=1000003$，与 C++ 的结果比较。

例如在固定 $d=7/13$ 下，以下有理数形式由独立 Blade 输出提供，其模 p 与新方法相同：

$$
I_{211}=-\frac{32}{39}I_{111},
$$

$$
I_{311}=\frac{1552}{1521}I_{111}+\frac{361}{2028}I_{110},
\qquad
I_{221}=\frac{608}{1521}I_{111}-\frac{361}{2028}I_{110}.
$$

这里没有将一次模 p 比较冒充对任意 d 的有理函数重构。

## 5. 文件与复现

- `native/fi_ibp.cpp`：从 J DE 和目标有理行建立最后一重积分的 IBP，精确检验证书，并加入显式端点映射。
- `examples/vac_fi_ibp.json`：证书 ansatz、端点坐标、端点值到 period 的映射、period 到原始 FI 的归一化。
- `native/compare_blade.cpp`：比较目标覆盖、prime、两个物理 FI 系数及未约化对象状态。
- `tools/check_vac_blade.wl`：独立调用 Blade；输出位于 `results/vac_fi_ibp/`，缓存位于同一目录，未修改 Blade 安装。
- `tests/data/vac_blade_dot3_modp.json`：实际 Blade 运行的 20 个模 p 结果，保留作独立回归基准。
- `tests/test_fi_ibp.sh`：从 J 展开重新跑到 FI 约化，再比较保存的 Blade 结果。

```bash
bash j_integral/tests/test_fi_ibp.sh
```

需要重新生成独立 Blade 数据时运行：

```bash
wolframscript -file j_integral/tools/check_vac_blade.wl
```

原始本轮 J 搜索采用 160 阶、训练至 100 阶、关系次数 6；回归采用 80 阶、训练至 50 阶，也得到相同 20 个 FI 约化。FI IBP 本身不依赖展开截断阶数，只依赖已经取得的 J 有理约化、DE 和端点映射。
