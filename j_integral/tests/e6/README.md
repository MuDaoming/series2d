# E6 的原包 AB 测试

FI定义为 $I_{\nu}^d=\int d^d k_1/(i\pi^{d/2})\,d^d k_2/(i\pi^{d/2})\prod_iD_i^{-\nu_i}$，其中 $D_i=q_i^2-m_i^2+i0$。传播子按下表编号，branch编号使用输入的0、1、2。

| i | $q_i$ | $m_i^2$ | branch |
| --- | --- | --- | --- |
| 1 | $k1$ | $0$ | 0 |
| 2 | $k1+p3$ | $0$ | 0 |
| 3 | $k2$ | $0$ | 1 |
| 4 | $k2+p5$ | $0$ | 1 |
| 5 | $k2+p1+p2-p3-p4$ | $0$ | 1 |
| 6 | $k1+k2+p1$ | $0$ | 2 |
| 7 | $k1+k2+p1-p4$ | $0$ | 2 |

独立外动量依次为 $p1$, $p2$, $p3$, $p4$, $p5$。全部精确标量积 $p_a\cdot p_b$ 为：

| |$p1$|$p2$|$p3$|$p4$|$p5$|
|---|---|---|---|---|---|
|$p1$|$1$|$3/2$|$-1/3$|$-3/8$|$-2/5$|
|$p2$|$3/2$|$1$|$-5/11$|$-3/7$|$-6/13$|
|$p3$|$-1/3$|$-5/11$|$0$|$1/30$|$1/34$|
|$p4$|$-3/8$|$-3/7$|$1/30$|$0$|$-10242329/2042040$|
|$p5$|$-2/5$|$-6/13$|$1/34$|$-10242329/2042040$|$0$|

固定 $d=26/7$、$p=1000003$。J使用 $X_1=(1-z)t$、$X_2=(1-z)(1-t)$、$X_3=z$，端点定义为输入中的fbi_analytic/regulated_zero，与已有物理DE一致。

完整主基有 **79** 个对象，分成 **47** 个原包依赖块；物理DE有 **1164** 个非零项。主基逐项在 `input/basis.json`，不是只做topsector cut。

原始家族来自 `archive/j_integral_sixpoint/inputs/family.json`；学习所需完整DE来自 `archive/j_integral_sixpoint/reference/connection_mod1000003.json`，SHA在 `input/provenance.json`。本次不重新获取anchor，也不把anchor的非对角系数提供给数值拟合。

运行 `bash j_integral/tests/e6/run.sh RUN_NAME`。变量平移、共同规范、原包命令、接口权限和记录位置见 [共同说明](../README.md)。

e6为原两圈六点例子，额外外动量满足 $p_6=p_1+p_2-p_3-p_4-p_5$；这里的Gram输入已经消去该动量。原输入来源保留在 `archive/j_integral_sixpoint/inputs/provenance.json`。

首次原包CLI运行 `original_cli_001` 已完成：学习5.00秒、求解18.07秒、真实数值IBP10次（6个不同训练点，含块检查共8点，再加2个最终检查点），无失败。79×79全部6241项与已有有理DE一致，撤回共同规范后的物理DE也逐条目完全相同。可追踪结果在 `reference/`，完整运行记录在 `runs/original_cli_001/`。

FI与J主积分数目的进一步比较及一个明确的t-IBP关系见 [analysis/README.md](analysis/README.md)。
