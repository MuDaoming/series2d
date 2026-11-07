# 原 AB 学习及生产接口核对

本目录隔离测试上游 `analytic_bootstrap` 0.3.0，不修改原包，也不修改主运行检查点。

## 官方生产链的数据契约

`analytic_bootstrap.production.run` 先准备 Blade 系统，通过 `upstream.prepared_anchor.acquire_selected` 在一个质数、一个其它运动学参数赋值上重构**完整的单变量有理 DE**。这里的 anchor 不是一个数值 $t$ 点：变量 $t$ 仍保持符号。

随后 `upstream.pipeline.learn_anchor` 从这个有理结果切出 $A_{11}$、$A_{22}$、$A_{12}$，对实际 $A_{12}$ 调用 `_bounds`，读取其约分后的分母与分子次数，并用 `learn` 学习分母向其它质数或运动学赋值的传输模板。`prepare_setting` 在后续 setting 获取对角块、传输分母模板、构造 AB 候选，并由 `ProductionFitter` 向数值 backend 请求系数值、确定参数及检查保留点。

因此，直接替换数值 backend 是可行的；但仅替换它，不能免去首份有理 anchor 的获取。当前项目从零构造第一份完整有理 $J$ DE，不能把尚未完成的答案输入 `learn_anchor`。官方的低层 `learn_profiles` 则可以对给定 $A_{11},A_{22},Q$ 与分子次数上限立即运行，并不要求已知未知 $A_{12}$ 的系数。

## 源导数阶的实际边界

`pl._order_records` 给每个存在 $A_{22}$ 极点的因子默认 $R=k$；不存在该极点则 $R=0$。`extra_order` 显式将前者改为 $R=k+\mathrm{extra\_order}$。

`learning.learn_profiles` 调用 `structural_orders.learn_order_profile`，后者在此有限参考 $R$ 内比较源分子模的前缀，选取能恢复该有限参考模的较小阶。它不会自行搜索 $R>k$，也不证明无限导数闭包。上游返回 `analytic_completeness=False` 与此一致。

## 当前 DP 失败块实测

`probe_orders.py` 复制当前已完成行到独立检查点，截获第一个待处理块及其输入，再调用**原包** `learn_profiles`、`build_pl_ansatz`。联合已知恒等式与数值约束仍沿用项目 `native_fit.joint_fit`；加 `--official-fit` 则额外直接调用原包 `ProductionFitter`。

每项约束对应

$$
A_{12}(t_s)=\sum_a\kappa_a V_a(t_s),
$$

而不是 $J(t_s)$ 的数值。数据只来自 $t=17,\ldots,57$ 的学习点；不读取最终盲验证目录。

首次测试块为以 `J{1,0,1,1,1,0,1,1}` 开头的 6 行，下层秩 66，396 个活跃条目，其中 120 个条目已有有理结果可固定。初始候选包络为 47474 个标量分子系数。

| 原包源阶选项 | AB 参数数 | 构造秒数 | 数值方程相容 |
| --- | ---: | ---: | --- |
| $R=k$ | 7175 | 5.82 | 否 |
| $R=k+1$ | 7175 | 6.81 | 否 |

官方 `learn_profiles` 耗时 5.77 秒，其中源模前缀测量 5.70 秒，分母传输模板学习 0.074 秒。63 个因子中 54 个具有正参考阶，9 个可选用更低生产阶。以上说明 $k\to k+1$ 对这个固定包络没有扩大最终 AB 候选；不能把该失败直接归因于缺少一阶源导数。

这些结果只是失败原因诊断，不是完整 DP DE 或约化已经完成的声明。详细机器结果在 `../results/dp_official/order*.json`。

补充：$R=k+3$ 仍为 7175 维，构造 8.03 秒，项目联合恒等式拟合仍不相容。独立使用原 `ProductionFitter` 的原生增量后端，在 0.413 秒内也抛出 `InconsistentSampleError`，与项目联合拟合的结论一致。这次直接拟合只输入数值点，不输入固定有理条目，故两者检验对象不同但排除结论相同。

上游增量拟合库默认没有随 PL 核一起编译。若缺少 `upstream_evaluate.so` 和 `upstream_fit.so`，原 `ProductionFitter(backend='native')` 会回退；因此实际复用此组件须先运行原包的 `analytic_bootstrap.upstream.build_native`，并设置 `ANALYTIC_BOOTSTRAP_UPSTREAM_KERNEL_DIR`。本次库与编译指纹保存在 `results/dp_official/upstream_native/`。首次缺库试验因进入慢回退而主动终止，没有将它当作算法所需耗时。

复现原生测试：

```sh
ANALYTIC_BOOTSTRAP_UPSTREAM_KERNEL_DIR="$PWD/j_integral_bootstrap/results/dp_official/upstream_native" \
ANALYTIC_BOOTSTRAP_FIT_THREADS=1 \
python j_integral_bootstrap/official_route/probe_orders.py --extra 3 --official-fit
```

实际传入 `learn_profiles` 的生产 `order_profile` 也重新用于 `build_pl_ansatz`；它仍得到 7175 个参数并拒绝同一数值数据。profile 的作用是减少源模构造工作，不会凭空修复一个不包含真实 $A_{12}$ 的有限候选。

## 廉价否决错误候选的接口

`preflight.preflight(ansatz, samples, max_points=32, kernel_directory=...)` 先用原包 `ProductionFitter` 增量加入少量真实数值约束。返回三种状态：

- `incompatible`：原包精确有限域消元证明此前缀已经矛盾，可直接跳过完整 fixed-entry 联合拟合。
- `compatible`：已检查所有传入的样本且获得唯一解；仍须进入完整恒等式约束和独立保留点验证。
- `undetermined`：此前缀相容但未定，或库不可用、执行错误；继续原来的完整流程，不能据此否决候选。

这里只有 `InconsistentSampleError` 被解释为数学矛盾。残差证书失败、数据格式错误、缺库等都不会误报候选不相容。使用 `backend='stream'` 显式禁止原 `backend='native'` 的慢回退。环境变量只在构造原生 fitter 时临时设置，退出时恢复；并行不同配置应使用不同进程。线程数沿用原包的 `ANALYTIC_BOOTSTRAP_FIT_THREADS`，并行任务可设置为 1。

此筛选先加入 $Cz=0$ 与少数 $E(t_s)z=y_s$，没有加入已经重构条目的有理恒等式。若此较弱系统已矛盾，再加入恒等式也必矛盾；因此先筛后做完整联合约束是安全的。相容前缀则不能证明包络正确，不能省去后续步骤。

`test_preflight.py` 对归档的 6×66 失败块及一个已完成 DP sector 使用真实学习点回归；不调用最终盲验证数据。

实测：失败 6×66 候选的前 3 点仍相容，rank 为 1177/7175；第 32 点首次暴露矛盾，原生筛选总耗时约 0.63 秒，因此默认预算采用 32 点，而非误认为 3 点即可否决。成功 DP 三行块的前 3 点返回未定（128/137），33 点全部检查后唯一相容，约 0.075 秒。此成功例的三个主对象为 `J{0,0,1,0,1,1,1,1}`、`J{0,0,1,0,1,1,1,2}`、`J{0,0,1,0,1,1,2,1}`。这些时长均不包含既有 AB ansatz 的构造。

主运行可这样接入：

```python
screen = preflight(ansatz, train_samples, kernel_directory=native_dir)
if screen['status'] == 'incompatible':
    # 记录本包络失败，直接进入下一个候选。
    continue
# compatible 和 undetermined 都继续加入 fixed-entry 的完整恒等式。
fit = joint_fit(ansatz, train_samples, fixed, identity_points)
```

生产构建现已纳入 `tests/build_kernels.py`：在 PL 核之外，调用未修改的原包 `upstream.build_native`，将增量 fitter 库及编译指纹写入稳定的 `build/upstream/`。两套核统一使用 `/usr/local` 的 FLINT。生产调用传入 `kernel_directory=ROOT/'build/upstream'`，不依赖上述历史结果目录。回归也已切换到稳定目录；另显式传入不存在的库目录，确认小于 512 维的成功 ansatz 仍只返回未定，不会误拒或进入慢回退。

## 完整 J 后端适配

`upstream_j.py` 分开实现 `AnchorBackend`、`DiagonalBackend` 和 `NumericBackend` 三种原包能力。`run_upstream.py` 直接调用原 `learn_anchor → prepare_setting → solve_de`；不经过项目自写的候选包络搜索或联合拟合。原包源码未改。

`run_case.py` 和 `cases/` 提供输入驱动的两阶段运行；`cold_anchor.py` 用数值 J-IBP 获取原包所要求的第一份完整有理 anchor。DP 的数值请求仅积分 485 个主基父指数并求 496 个导数，不顺带处理完整 dot≤2 测试集合。

两圈真空已实际跑通冷启动和原包生产两个阶段。2412.21054 的六点例子已完成原包全部47层，其79维物理连接与独立重构精确一致，并通过另取的8个全新完整数值点。见 `../../j_integral_sixpoint/docs/problem_solution.md`。

天河调度、环境前提和目前未提交的状态见 `tianhe/README.md`；不要把本地机械预检通过误读成远端正在跑。
