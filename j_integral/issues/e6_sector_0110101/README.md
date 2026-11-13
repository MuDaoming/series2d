# e6 四传播子家族：J 与 FI 主积分数目不同

## 1. 本问题研究什么

只取e6的原sector `0110101`，也就是原来的D2、D3、D5、D7。它有三个branch，传播子数为 **(1,2,1)**。在原e6的已有记录中，该sector有2个FI主积分，而当前J基中有3个普通J主积分。

本目录先把家族、归一化、输入和已有证据独立固定下来。要研究的核心是：这一个sector多出的J如何在最后一重t积分中处理，边界对象是否独立，以及原有J约化是否还遗漏物理关系。不预先把“最后积分一定消掉它”作为结论。

## 2. 独立的四传播子定义

原e6使用圈动量k1、k2。作平移并定义两个外动量：

$$
\ell_1=k_1+p_3,\qquad \ell_2=k_2,
$$

$$
P=p_1-p_3-p_4,\qquad Q=p_1+p_2-p_3-p_4=P+p_2.
$$

原来的四条传播子现在重新编号为

$$
D_1=\ell_1^2+i0,\qquad
D_2=\ell_2^2+i0,
$$

$$
D_3=(\ell_2+Q)^2+i0,\qquad
D_4=(\ell_1+\ell_2+P)^2+i0.
$$

全部内部质量为零。编号对应关系是

| 本目录编号 | 原e6编号 | branch（输入编号） |
| --- | --- | --- |
| D1 | D2 | 0 |
| D2 | D3 | 1 |
| D3 | D5 | 1 |
| D4 | D7 | 2 |

数学上三个branch分别为 $\mathcal B_1=\{1\}$、$\mathcal B_2=\{2,3\}$、$\mathcal B_3=\{4\}$。输入JSON使用从0开始的branch编号，不改变上述分组。

原始FI定义为

$$
I_{\nu_1\nu_2\nu_3\nu_4}^{d}
=\int\frac{d^d\ell_1}{i\pi^{d/2}}
      \frac{d^d\ell_2}{i\pi^{d/2}}
      \frac1{D_1^{\nu_1}D_2^{\nu_2}D_3^{\nu_3}D_4^{\nu_4}}.
$$

本目录考虑非负整数指标；0表示删除传播子。topsector为四个指标均为正。没有额外辅助传播子或不可约分子。

## 3. 运动学及本次数值设置

家族只有三个独立外部不变量：

$$
s_P=P^2,\qquad s_Q=Q^2,\qquad s_R=(Q-P)^2,
\qquad P\cdot Q=\frac{s_P+s_Q-s_R}{2}.
$$

从原e6的完整Gram矩阵精确得到

$$
s_P=\frac{149}{60},\qquad
s_Q=\frac{38113}{4620},\qquad
s_R=1,\qquad
P\cdot Q=\frac{22483}{4620}.
$$

因此不需要继续携带原来的五个外动量。新的family输入使用两个独立外动量P、Q及其2×2 Gram矩阵。

当前复现设置为

$$
d=\frac{26}{7},\qquad p=1000003,\qquad t=17,
\qquad \tau=t-\frac25.
$$

$p$是有限域质数；维数和运动学以精确有理数输入。t是J的剩余变量，17只是本次数值IBP检查点；保存的DE与关系证书仍保留完整t有理依赖，系数按tau升幂排列。

## 4. FBI 与 J 的明确定义

固定使用

$$
X_1=(1-z)t,\qquad X_2=(1-z)(1-t),\qquad X_3=z.
$$

三个branch的指标之和为

$$
\sigma_1=\nu_1,\qquad \sigma_2=\nu_2+\nu_3,\qquad \sigma_3=\nu_4,
\qquad v=\nu_1+\nu_2+\nu_3+\nu_4.
$$

在topsector，branch1与3各只有一条传播子，内部参数固定为1；branch2令D3的参数为y、D2的参数为 $1-y$。圈动量配方给出

$$
A=\begin{pmatrix}X_1+X_3&X_3\\X_3&X_2+X_3\end{pmatrix},
\qquad U=\det A=X_1X_2+X_1X_3+X_2X_3,
$$

$$
B=\begin{pmatrix}X_3P\\X_2yQ+X_3P\end{pmatrix},
\qquad C=X_2ys_Q+X_3s_P,
\qquad F=B^{\mathsf T}\operatorname{adj}(A)B-UC.
$$

等价地，

$$
F=-X_1X_2X_3\left[s_P-(s_P+s_Q-s_R)y+s_Qy^2\right]
-UX_2s_Qy(1-y).
$$

保留Gamma及符号归一化的FBI为

$$
\mathcal I_{\boldsymbol\nu}^{d}(X)
=\frac{(-1)^v\Gamma(v-d)}{\prod_{i=1}^4\Gamma(\nu_i)}
\int_0^1dy\,
\frac{(1-y)^{\nu_2-1}y^{\nu_3-1}}{[F(y;X)-i0]^{v-d}}.
$$

上述单个y积分公式针对四个指标均为正；删除D2或D3时应按剩余活跃传播子重定义参数测度，不能直接代入Gamma(0)。

$$
G_{\boldsymbol\nu}(t,z)
=(1-z)X_1^{\nu_1-1}X_2^{\nu_2+\nu_3-1}X_3^{\nu_4-1}
U^{v-3d/2}\mathcal I_{\boldsymbol\nu}^{d}(X),
$$

$$
J_{\boldsymbol\nu}(t)=\int_0^1G_{\boldsymbol\nu}(t,z)dz,
\qquad I_{\boldsymbol\nu}^{d}=\int_0^1J_{\boldsymbol\nu}(t)dt.
$$

端点处理配置暂沿用原e6的 `fbi_analytic` / `regulated_zero`。这记录当前程序所用条件，不代替本问题对物理端点的进一步推导。最后的t端点不能因为z端点的配置而直接置零。

## 5. 已有记录与待确认范围

原e6中该topsector的FI主对象为

$$
I_{1111},\qquad I_{2111}.
$$

继承到本家族的低sector还有 $I_{1011}$、$I_{1101}$，合计4个FI候选主对象。这是从原e6的Blade清单筛选、重编号得到的记录，尚未在独立四传播子输入上重新运行Blade确认最小性。

现有J的完整微分闭包为6个对象：

$$
\mathbf J=(Bz0_{1111},J_{1011},J_{1101},J_{1111},J_{1112},J_{1121})^{\mathsf T}.
$$

其中同一topsector的三个普通J为 $J_{1111},J_{1112},J_{1121}$；两个低sector普通J与一个边界函数必须保留。`Bz0`表示规范化corner FBI的z=0正则化端点，不是J在t=0的值。

`reference/j_de.json`保留6维物理DE，未施加AB测试中的共同规范变换。它从已经验证的e6完整连接抽出，并作4指标重编号，没有引入新的系数拟合。新的四传播子输入在t=17直接计算的36个DE系数须与之比较。

已有一个精确t-IBP恒等式（模 $p=1000003$）：

$$
\frac{d}{dt}(g_1J_{1111}+g_2J_{1112}+g_3J_{1121})
=J_{1111}+734168J_{1112}+56454J_{1121}
+\sum_{a\in\mathrm{lower}}R_a(t)J_a(t).
$$

lower包括两个低sector J及Bz0。全部g和R记录在 `reference/t_ibp_certificate.json`，这是原e6证书的精确重编号；其t端点项仍保留，没有用它宣布FI约化已经完成。

## 6. 接下来要研究的问题

1. 用这个独立四传播子输入确认FI和J的主积分计数，检查不同dot、维数及非奇异运动学点。
2. 从上述显式FBI定义推导端点函数，核对当前物理边界条件及是否遗漏J函数关系。
3. 完整处理t-IBP恒等式中的低sector积分和端点项，得到真正的FI约化并与Blade比较。
4. 区分“固定t时的J函数独立性”和“最后积分后的FI独立性”，不要把全导数误当作固定t时的零关系。

本次只建立研究对象并核对输入等价性，不预先宣称所有差异都已解释。

## 7. 文件与运行

- `input/family_request.json`：独立四传播子家族；保留sP/sQ/sR符号及本次精确赋值；可直接交给通用J-IBP程序。
- `input/index_map.json`：原e6到本家族的指标、动量和不变量映射。
- `reference/j_de.json`、`fi_masters.json`：继承的DE和FI候选主基。
- `reference/t_ibp_certificate.json`：重编号后的已有t-IBP证书。
- `reference/provenance.json`：来源文件指纹及计数范围。
- `check_definition.py`：检查外动量变换、Gram等价性及本次36个数值DE系数。

从项目根目录运行：

```sh
j_integral/src/IBP/build/j-ibp \
  --input j_integral/issues/e6_sector_0110101/input/family_request.json \
  --output /tmp/e6_fourprop.json
python3 j_integral/issues/e6_sector_0110101/check_definition.py /tmp/e6_fourprop.json
```

## 解析积分尝试

见 [解析公式与推导](analytic/README.md)：两个低 sector 的 $J$ 和边界函数已完全积分，三个顶层 $J$ 暂得到完成 $z$ 积分的超几何表示；检验脚本及结果同时保存。
