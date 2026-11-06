# Analyticity Bootstrap 原包阅读与 DP 接入诊断

本文针对用户下载的 `extern/analytic_bootstrap-master`（0.3.0），记录实际阅读的模块、数学对象及可复现测试。原包副本未修改。论文是 [Analyticity Bootstrap of Feynman Integrals, arXiv:2610.01154v1](https://arxiv.org/html/2610.01154v1)。生产前端的独立实验另见 `official_route/`；本文重点是原包模构造和拟合端。

## 论文中的问题与我们的输入

固定 $p=1000003$、$d=26/7$ 和运动学，将一个 sector 的待求 J 主对象记为 $T$，已经求完的低层对象记为 $S$：

$$
T'=A_{11}T+C S,\qquad S'=A_{22}S.
$$

已知 $A_{11}$、$A_{22}$，求 $C$。它们均是 $t$ 的有理矩阵；实现用 $\tau=t-2/5$。这里 $S,T$ 都是积分掉 $z$ 后的 J 及必要边界对象，不是 FBI。实际数值 IBP 提供普通点的 $C(t_a)$，不是 Taylor 系数。

对 $C_{ij}$ 先指定分母 $D_{ij}$ 和分子次数 cap；将满足解析约束的有理矩阵写成 $C=\sum_a\kappa_a W_a$，最终由数值 IBP 固定 $\kappa_a$。因此“需要多少点”同时依赖包络和解析约束，不能仅凭单变量或主积分数量推断。

论文的解析区域有限性是前提，适当规范变换可以恢复它。共同标量变换 $\widetilde T=hT,\widetilde S=hS$ 给出

$$
\widetilde A_{11}=A_{11}+\frac{h'}h\mathbf1,\qquad
\widetilde A_{22}=A_{22}+\frac{h'}h\mathbf1,\qquad
\widetilde C=C.
$$

故改变 $h$ 后仍可拟合同一份 $C(t_a)$。这种变换和增加分母次数、增加导数阶数是三个不同动作，应分别测试。

## 原包模构造的真实步骤

入口 `analytic_bootstrap/pl.py:build_pl_ansatz` 先把分母包络按不可约因子 $q$ 分解。每个因子取公共极点上界 $k$，在 $\mathbb F_p[\tau]/(q^k)$ 中计算，最后合并不同因子的主部及多项式部分。

1. `backends/pl/factor.py:pure_sub_space` 计算论文的 $N_{A_{22},Q}$，$Q=q^k$。`native_generators` 或 Laurent 路径生成
   $$r_i^{[0]}=e_i,\qquad r_i^{[m+1]}=\partial_\tau r_i^{[m]}+r_i^{[m]}A_{22}.$$
   随后 Hasse pole-state 消元及 `project_global_numerator_module` 得到局部分子行。它实际包含导数递推、辅助极点状态消元、规范化和投影，并非一次黑箱“构造空间”。
2. `extension.py:quotient_module` 构造论文的低层商 $V_{A_{22}}$。当前实现保留分子矩阵的 RREF chart 和自由列；乘 $q$ 的作用按需读取，避免创建稠密投影和移位矩阵。
3. `topmodule.py:top_pole_module` 从 $A_{11}$ 的主部构造 $P_{A_{11}}$：提取 $Q A_{11}\bmod Q$ 的行，做截断环上的 Popov 分解，保留生成元 valuation、原行坐标和 syzygy。
4. `extension.py:extension_hom` 计算 $\operatorname{Hom}_{\mathbb F_p[\tau]/(Q)}(P_{A_{11}},V_{A_{22}})$。对 valuation 为 $v_i$ 的源生成元，需要商上乘 $q^{k-v_i}$ 的核；相同指数共用核，但各生成元仍有独立参数。原包有懒计算移位代表和分块 Hom map，没有把所有对象稠密展开。
5. `factor.py:assemble_factor_allowed` 将上述兼容条件与各矩阵条目的允许极点/支撑相交；纯低层因子还可提前按支撑组处理。`factor_basis_apart` 或融合坐标输出得到 $W_a$ 的局部主部系数。
6. `cap_structure.py:impose_caps_structural` 用无穷远 moments 施加全局分子 cap。主部的多个因子会在无穷远相互抵消，因此不能逐因子各自截断代替这一步。结果用隐式核保存，不必显式构造所有独立参数的稠密基。

`FactorBlock.metadata` 已原生记录上述步骤的时长；接入层应保存这些信息，而不是只输出“空间构造耗时”。

## 原包拟合与我们添加的部分

`sampling_batch.py:native_joint_fit` 直接联立次数约束 $Kz=0$ 和各点 $E(t_a)z=C(t_a)$，其中 $z$ 是局部主部及正则参数。它不先反复构造 $E(t_a)\ker K$。按单目标私有列优先、共享 Hom 列最后排序，并让点值方程先于次数约束进入原生稀疏消元。

`tests/native_fit.py:joint_fit` 参考这个已有实现，增加已知有理条目的恒等式约束。已知条目分母整除输入分母，且分子满足 cap；给足 cap+1 个普通点后，它们固定的是有理恒等式。原包的联立消元思想不是我们新设计的算法。

原生消元返回相容见证后，两种实现均检查原始方程残差、局部参数与独立参数往返、独立 ansatz 求值。秩不足只表示相容但尚未唯一，不能当作完成。主脚本额外留内部检查点；最终 $t=151,\ldots,158$ 的全局盲点与学习完全隔离。

原包还有独立 `ProductionFitter` 的原生后端；是否实际使用原生库必须检查，不能因 API 名称包含 production 就假定已走最快路径。该路径的安装与实测由 `official_route/` 单独记录。

## 原包的 learning 到底学什么

`learning.py:learn_profiles` 已实现 source-order、target-prefix、Laurent precision/window、cap layout 和 denominator transport 的学习及冻结。`optimized_learning` 开启 quotient 学习、相邻复用、分子缓存和图/布局缓存。DE 模式没有 reduction 模式的 target-prefix 学习，不能把两种模式的成本和功能混为一谈。

`upstream/envelope.py:learn` 从一个已知有理 anchor 的真实非零条目读出分母因子幂、分子相对分母的 growth 和零掩码；`transport` 按对角 witness 的结构签名把这些模板传到其他 setting。因子签名不可区分时默认拒绝，或显式允许以组内最大幂给上界。它不从单个数值 $t$ 点猜出完整未知有理函数的分母。

因此目前 DP 自写的候选分母尝试不是原包 source-order learning 的替代：它解决的是原包 anchor 所需的真实包络尚未拿到的问题。但在包络确定后继续手写这些 profile 学习没有必要，应优先使用原包接口。

## 导数阶数不是已证明的完备上界

`pl.py:_order_records` 默认在低层矩阵于该因子有极点时取 $R=k$，否则 $R=0$；允许显式 `orders` 或 `extra_order`。返回 `analytic_completeness=False`。

原包 `tests/test_public_pl.py:resonant_spec` 本身就给出反例：$k=1$，默认 $R=1$ 得零维；显式 $R=9$ 才找到一维。`structural_orders.py` 的学习是在给定 reference order 内找更小可复用阶数；代码要求 $s_*\le R_{\rm prod}\le R_{\rm ref}$。所以它不会自动发现任意超过默认 reference 的迟发共振。不能把“学习所得阶数稳定”说成理论完备。

## 实际 DP 失败块的隔离实验

可复现入口为 `diagnosis/profile_block.py`，输出 `results/dp_diagnosis/profile_block.json`。它复制主连接检查点，只在独立 snapshot 上运行，没有修改主驱动或上游源码，也未获取新 IBP 点。测试对象为

$$
\begin{aligned}
T=(&J_{10111011},J_{10111012},J_{10111021},\\
   &J_{10112011},J_{10121011},J_{10111013})^T,
\end{aligned}
$$

下层秩 66，共 396 个耦合条目；其中 120 个已有有理结果。采用原主驱动候选包络和 PL-structural 窗口，38 个训练点，3 个内部留出点。固定分母、cap 和规范因子，仅改变 $R$：

| extra_order | 解析约束后参数数 | 构造秒数 | 数值 $\kappa$ 方程 |
| --- | ---: | ---: | --- |
| 0 | 7175 | 5.31 | 不相容 |
| 4 | 7175 | 9.42 | 不相容 |
| 12 | 7175 | 25.40 | 不相容 |

这一实验排除了“只需在这个候选包络上多算至额外 12 阶就能修复”的说法，不能据此证明任意阶完备。默认时长按论文对象分解：

| 对象/操作 | 秒 |
| --- | ---: |
| $N_{A_{22},Q}$，全部因子之和 | 2.795 |
| $V_{A_{22}}$ 的 quotient chart | 0.713 |
| $P_{A_{11}}$ | 0.116 |
| $\operatorname{Hom}(P,V)$ | 0.077 |
| 条目极点/支撑相交 | 0.847 |
| 有理主部坐标输出 | 0.603 |
| 全局分子次数限制 | 0.008 |
| 生成点值、已知条目和次数的 $\kappa$ 方程 | 8.213 |
| 原生消元这些 $\kappa$ 方程 | 1.352 |

提高 extra_order 到 12 后，仅 $N$ 步骤已耗 22.05 秒，而参数数不变。因此这个块不应继续盲加导数阶数。当前更大的可避免成本是错误包络/解析规范导致重复构造和重复生成方程，Hom 本身不是瓶颈。

## 排除低层输入与 cut 对角不一致

`diagnosis/check_inputs.py` 将所有已知非零函数对 $t=17,\ldots,57$ 的可用学习行逐项求值；未读取最终盲点。结果在 `results/dp_diagnosis/input_validation.json`：

| 输入 | 非零有理条目 | 实际比较次数 | 不一致 |
| --- | ---: | ---: | ---: |
| snapshot 的已完成连接 | 7805 | 340823 | 0 |
| cut 同 sector 输入 | 1998 | 78262 | 0 |
| exceptional 对角输入 | 289 | 21002 | 0 |

没有一次求值落在分母零点。这排除了这些已知函数与当前学习数据不匹配；还不是对整个未知连接的最终验证。

## 后续规范与已知函数隔离

`diagnosis/capture_block.py` 保存了冻结的 A、B、包络、已知函数及训练点到 `results/dp_diagnosis/capsule.json`，重现实验不再读取全局巨型数据文件。仅加入原有 120 个已知有理条目时，相容、增加秩 473；实际未知条目的点值才使系统不相容。

`diagnosis/normalize_capsule.py` 保持同一包络，进一步隔离规范：

| 相对于初始端点三次规范的改变 | 参数数 | 结果 |
| --- | ---: | --- |
| 端点再乘三次，即共六次 | 7706 | 不相容 |
| 端点再乘九次，即共十二次 | 8186 | 不相容 |
| 其余候选因子各乘三次 | 46463 | 相容但 38 点仅秩 12320 |

最后一项严重欠定；它不能证明失败根因是非端点规范。放松足够多约束后，有限样本也可能被错误分母包络容纳。这一点必须和找到真实有理见证区分。端点幂升高也不能脱离有限导数阶数讨论理论完备性。

`diagnosis/check_new_coefficients.py` 进一步检查独立极点学习路线新增的真实 C。当前可用 189 个函数的真实分母均整除冻结候选分母，分子均满足原 cap；没有发现候选外因子或更高阶极点。这个结论仅针对已经重建的条目，不能外推到尚未重建的条目。

原包 `test_public_pl.py` 的四项回归已原样运行并通过，包括迟发共振反例；记录见 `results/dp_diagnosis/upstream_regression.log`。

补充检查也包含已完成连接的隐含零：共 6,772,797 次零分量比较，无遗漏非零项，见 `known_zeros.json`。合并补充路线后，已知真 C 增至 276 条，以其精确分母和 cap 构造的 4207 维 ansatz 仍容纳全部这些有理恒等式。尚未知条目的数值不相容不能从这些已知函数上直接归因。

冻结胶囊可用 `diagnosis/profile_capsule.py` 重放 R 对照，不依赖后续变化的主检查点。原包 `test_structural_extension.py` 两项回归也已通过，覆盖被自动消去的 annihilator syzygy 与必须保留的 mixed syzygy，包含线性及二次不可约因子和不同有限域。
