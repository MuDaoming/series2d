# DP 原始 FI 的 maximal-cut 约化与 Blade 对照

## 1. 范围和当前结果

本测试只处理 DP 的原始 FI maximal cut：8 条物理传播子全部被 cut，下 sector 在这个商空间中为零。没有运行普通完整家族约化。

输入为 `examples/dp_topcut.json`，$d=26/7$、$p=1000003$，与 J 层测试使用相同质量和运动学。

必须区分三个计数：

| 对象 | 当前数目 |
|---|---:|
| 旧一般变域 F(delta) 的 topsector 对象 | 43 |
| 修正后 J(t)/边界闭合基 | 8（7 J + 1 边界） |
| 原始 FI maximal-cut master 数 | 15 |

最后一行由本轮 Blade maximal-cut 分析得到，使用 IntegralOrdering=2、PreferredMI=Dot，15 个 master 均无 rank，最高 dot 为 2。

## 2. J 层数据的独立检查

J 层分别在 $(p,t_0)=(1000003,2/5)$、$(1000033,3/7)$ 下展开 200、220 阶，均找到同一组 8 个对象并获得闭合 DE。随后在每个 prime 的 8 个新 t 点直接通过已知 FBI 的 z 积分证书，检查 82 条约化、8 条 DE 和 45 个零对象，全部通过。

这些检查不使用训练点 Taylor 系数。详细定义、计数修正原因及复现命令见 `dp_test.md`。

## 3. 从 J 到 FI 的 IBP

对全部总 dot≤2 的 45 个 FI，使用修正后的 J 约化和 DE，搜索

$$
\boldsymbol r_\nu
=\sum_a c_{\nu a}\boldsymbol r_{\eta_a}
+\boldsymbol g_\nu'+\boldsymbol g_\nu A
+\boldsymbol h_\nu^{\mathrm{boundary}}.
$$

$c_{\nu a}$ 为有限域常数，不依赖 t。边界模的有理余项保留在输出中，只有在声明的原始 FI cut 投影后才被视为下 sector；不能只把一个不带权的边界 period 设为零，却遗漏带有理权重的边界项。

`fi_ibp.cpp` 验证边界函数确实构成 DE 的不变子模，投影后检查每条有理证书，再回到完整 J 基检查余项全部落在指定边界模中。t 端点也按本次 cut 的 pinched-branch 约定处理。这不是完整、未 cut 家族的边界求解。

本轮证书公分母取 J DE 和目标约化分母的最小公倍式，次数为 39，证书分子次数上限 64。

### 不使用 Blade 基的独立选择

将全部 45 个目标作为候选 FI masters，先消去全导数方向，再选择剩余独立目标。程序自动保留 15 个 FI，45/45 个有理证书均通过；本机最后一重 IBP 阶段实测约 2.56 秒。

配置：`examples/dp_fi_ibp_cut.json`；输出：`results/dp_fi_ibp/fi_reduction_auto.json`。

### 转到 Blade 的同一组基

仅读取 Blade 选出的 15 个 master 的指标，重新在这组基上求解；不读取 Blade 的约化系数。结果仍是 45/45 个证书通过，实测约 1.66 秒。

配置：`examples/dp_fi_ibp_blade_basis.json`；输出：`results/dp_fi_ibp/fi_reduction_blade_basis.json`。

这两个计时都不包含先前 J 展开、DE 搜索或 Blade 计算，不能作为完整流程的总耗时。

## 4. Blade 范围和对照状态

`tools/check_dp_blade_cut.wl` 使用原始圈动量传播子，额外补齐 3 个 ISP，但目标和最终基都没有 rank。第一批只计算 corner 和 8 个单 dot，共 9 个目标；maximal-cut 分析选出的 dot2 masters 也传入作为优选基。

完成 master 分析和 cut IBP 生成后，复用其 8100 条方程（3281 个变量），用 Blade 相同的 FiniteFlow 稀疏求解器直接数值求解，绕过 BLReduce 后续重复构造符号拓扑的流程。方程、master 选择及数值求解都不读取 J 路线的约化系数。

最终独立求解覆盖了全部 45 个 dot≤2 目标，而不仅是最初 9 个小批目标。求解器得到 15 个独立变量、30 个依赖目标，15 个基全部无 rank。先将精确数值有理系数模 p，再逐目标逐系数比较：**45/45 一致，零差异**。其中 15 个是基的恒等行，30 个为非 master 目标的约化。

独立脚本为 `tools/check_dp_blade_generated_ibp.wl`；Blade 方程生成约 274 秒，复用已生成方程后的数值求解和导出约 2.72 秒。原先 BLReduce 的冗余小批后处理已停止，没有启动普通全家族约化。

实际比较报告为 `results/dp_fi_ibp/comparison_dot2.json`，记录 `full_target_coverage=true`、`matched=45`、`mismatches=[]`。独立基准保存在 `tests/data/dp_blade_dot2_cut_modp.json`，并由端到端回归重新核对。

## 5. 复现

```bash
# J 数据的重跑见 dp_test.md。
j_integral/build/fi_ibp j_integral/examples/dp_fi_ibp_cut.json \
  /tmp/dp_fi_auto.json

# 只运行 Blade 的 cut 小批对照，不是普通全家族约化。
wolframscript -file j_integral/tools/check_dp_blade_cut.wl

# 在上述方程生成后即可复用它们求解全部 45 个目标；不需等待冗余小批流程。
wolframscript -file j_integral/tools/check_dp_blade_generated_ibp.wl

j_integral/build/fi_ibp j_integral/examples/dp_fi_ibp_blade_basis.json \
  /tmp/dp_fi_blade_basis.json
j_integral/build/compare_blade /tmp/dp_fi_blade_basis.json \
  j_integral/results/dp_fi_ibp/blade_modp_dot2_generated_ibp.json \
  /tmp/dp_cut_comparison.json

# 使用已保存的独立 Blade 基准从头回归，无需重跑 Blade。
bash j_integral/tests/test_dp_fi_ibp.sh
```
