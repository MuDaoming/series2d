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

## J的对应规模测试

`input/j_common_prime.json` 保留既有DP家族、branch选择和496个J/边界主对象。将目标替换为全部2952个保留三个branch的dot≤2积分，并设置 `derivatives=[]`，避免将求DE的开销混入纯约化。单点为 $t=17$、$d=26/7$、$p=9223372036854775783$。调用当前通用C++接口 `j_integral/src/IBP/build/j-ibp`，其中包含与e6相同的矩阵规模诊断。

`j/run.py` 记录真实进程耗时及接口返回的 `diagnostics`；完整输出为 `j/result.json`、诊断日志为 `j/run.log`、汇总为 `j/statistics.json`。本机服务名 `dp-j-ibp-size.service`，内存上限4 GiB、不使用swap；与FI任务独立运行。

比较范围必须区分：FI运行4095个目标，J运行其中保留三个branch的2952个共同目标。其余1143个FI目标缺少branch；DP有质量，不能像e6那样直接称这些目标全部为零。最终结果需同时报告完整FI系统范围与这一共同目标范围。J的z关系只是在FBI解析递推及Hermite约化之后提交的补充方程，不冒充全部原始IBP方程；还需另外报告FBI矩阵操作数和最大形状。

J单点已完成：106.0326秒，496个J/边界主对象。FBI部分147个sector、115个FBI主对象、141次小矩阵消元，总343行，最大7×7；矩阵形状为45个1×1、48个2×2、27个3×3、8个4×4、1个5×5、12个7×7。z方向提交799条非零关系，215条独立，5695个辅助积分坐标、34530个非零系数。最终目标/边界表示共3182个向量、4954个形式坐标、秩496。完整诊断见 `j/statistics.json`，精简统计见 `j/summary.json`；完整单点结果压缩保存在 `j/result.json.gz`。这些矩阵行列属于不同算法阶段，不相互相加冒充一个Laporta矩阵。

2026-11-13 17:11:44（北京时间）按约定复查：FI服务active/running，主PID60211，累计CPU455.46秒，当前内存约0.74 GiB；已完成maximal-cut主积分识别（547个候选，不冒充完整约化后的最终主积分数），开始原包rank=2、dot=2方程生成。J服务已成功退出，日志确认 `masters=496 unclosed=0`。状态快照保存在 `status_checks/20261113T091144Z.json`，由只读的 `check_status.py` 生成。
