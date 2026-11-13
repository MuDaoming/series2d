# DP原始FI的Blade原包约化规模

与e6使用相同口径：固定 $d=26/7$ 及已有DP运动学，不建立运动学微分方程、不做cut，约化8个物理传播子所有非空sector中无分子、总dot不超过2的4095个目标。总dot定义为 $\sum_i\max(\nu_i-1,0)$；3个辅助传播子的目标指标均为0。

`input/fi_family.json` 来自 `j_integral/archive/j_integral_bootstrap/reference/fi_master_counts/dp_input.json`，列出完整的11个分母、运动学和维数。`input/all_targets.json` 明确列出4095个目标。此前只有DP topsector cut的Blade约化缓存，不能用它代替完整家族的规模。

使用公开入口 `BLFamilyDefine`、`BLMaximalCutMasters`、`BLReduce`，不自己选择种子。除dot优先主积分排序和单线程资源设置外，使用原包默认流程，与e6保持一致。`Cut` 显式全0并在启动脚本中断言。原包完成后自动执行只读的 `collect.wl`，逐个系统统计：初始非零方程数、不同FI数、筛选后方程数、筛选方程涉及的不同FI数。不同FI按完整指标去重；各子任务计数之和不冒充跨任务去重数。

## 后台运行

本机systemd用户服务 `dp-fi-native-dot2.service`，单线程，内存上限12 GiB、禁止该服务使用swap。整个服务进程组统一管理，避免中止Mathematica后残留子进程继续写缓存。启动日期2026-11-13；这次不使用天河/Slurm。

```bash
systemctl --user show dp-fi-native-dot2.service -p ActiveState -p SubState -p Result -p MemoryPeak
cat j_integral/tests/dp/fi_vs_j_dot2/fi/blade_native/run_status.json
```

日志和结果位于 `fi/blade_native/`：

- `run.log`、`cache/`：原包日志与方程缓存。
- `masters.json`：原包识别的主积分。
- `completion.json`、`reductions.wl`：原包完整约化完成后的信息及结果。
- `collect.log`、`system_sizes.json`：自动统计过程及规模结果。
- `run_status.json`：worker阶段；若服务被OOM终止，此文件可能停留在旧阶段，应以systemd的 `Result` 为准。

worker不会轮询；前一阶段结束后自动执行后一阶段。此次完整公开入口也可能做多素数有理数重构，因此总耗时不能叫作缓存就绪后的一次有限域求值时间。若触及内存上限，保留失败证据和缓存，再评估转移计算，不擅自修改播种参数缩小系统。
