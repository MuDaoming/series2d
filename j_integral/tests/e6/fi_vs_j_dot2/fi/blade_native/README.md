# Blade原包完整约化的规模测试

与手动生成方程的 `../numeric.wl` 区分，此目录只通过公开入口 `BLFamilyDefine`、`BLMaximalCutMasters` 和 `BLReduce` 驱动Blade。使用同一个 `../../input/fi_family.json` 家族和 `../../input/all_targets.json` 的1695个目标积分，维数、运动学完全数值化。

保留dot优先的主积分排序（`IntegralOrdering=2`、`PreferredMI="Dot"`）以及单线程资源设置。播种、`BlackBoxRank`、`FilterLevel`、`MaxIncrement`、任务拆分等使用Blade自身默认值；不定义 `GenSeeds`，不调用手动的 `GenIds`。固定数值参数时Blade自动走常规IBP求解分支，不搜索解析参数的块三角形式。

运行：

```bash
wolframscript -file j_integral/tests/e6/fi_vs_j_dot2/fi/blade_native/run.wl > j_integral/tests/e6/fi_vs_j_dot2/fi/blade_native/run.log 2>&1
wolframscript -file j_integral/tests/e6/fi_vs_j_dot2/fi/blade_native/collect.wl
```

`run.log` 及 `cache/` 保留原包日志与方程；未修改Blade原包。`collect.wl` 只读取缓存，对各个系统统计初始非零方程数、其中不同FI数、原包筛选后的方程数以及筛选方程实际涉及的不同FI数，保存到 `system_sizes.json`。不同FI按完整指标去重，不预先用IBP/对称关系合并。跨子任务的相同积分不可把各子任务计数相加后声称全局去重数。

`results/data` 是Blade筛选后保存的单个有限域数值约化结果，可用于本次规模测试和模素数正确性验证。公开入口之后可能继续进行多素数的有理数重构；那部分时间和素数数目应与一次有限域约化区分。

## 完成结果

原包完整约化完成，1695个目标全部返回。一个整体系统，无自动加大种子的第二轮：初始1294356条非零方程、497083个不同FI；原包筛选后64818条方程、72134个不同FI。计数由 `collect.wl` 读取原包文件得到，见 `system_sizes.json`，不是手动生成的替代方程。原包另有优先主积分的符号标签，不与FI指标混为一谈。

生成436.005秒；首次加载、学习、筛选及保存数值约化462.611秒；后续重新加载筛选方程并进行有理数重构117.299秒。完整公开入口 `BLReduce` 约1020.675秒，连同主积分识别的整个脚本1151.405秒。以上包含前处理和多素数重构，不能当成缓存就绪后一次有限域求值的耗时。原包日志的重复空参数取样计时也不能作为独立新点IBP的可靠计时，本次未据此报告单点速度。

成功任务的各阶段原包日志以 `.gz` 保存，完整有理约化为 `reductions.wl.gz`。先前一次错误地设置了maximal cut的启动已废弃；其残留进程在重试时造成了提前读取未完成文件，相关无效输出已隔离，不用于上述统计。正确任务重新生成完整方程，完整性检查、数值阶段和有理数阶段均完成。

`verification.json` 记录全部1695个目标约化系数的模素数核对通过；显式合并同一个等质量无质量sunrise别名后，与此前已验证的数值结果逐项一致。`verify.wl` 使用返回列表末尾的62个优先主积分结果，将原包 `BLInternalEX` 标签还原为FI，再比较系数。
