# 原包生产机制补充核对与实际接入

本次继续阅读 `production.py`、`reduction_production.py`、`upstream/pipeline.py`、`target_prefix.py`、`learning_fused.py`、`learning_context.py`、`support_profile.py`、`backends/pl/graph_learning.py` 和实际 native 入口，并结合生产合并、early-support、DE Hom/support 及统一学习文档核对。下面区分已经实测并接入的改善与当前不应优先投入的机制。

## 已落实：此前没有构建完整优化 runtime

原项目 `build/kernels` 只有六个基础核。原包 `laurent_production._fast` 在缺少 persistent worker 或共因子库时，进入 `legacy_laurent`；后者还检查 fused 可执行文件。此前 fused 文件也缺失，因此运行的是兼容基础实现。这不是缺一个 Python 配置，而是构建缺少生产算法实际使用的可执行文件与库。

现在由 `build_production_runtime.py` 将**未修改的原包构建脚本**复制到临时目录，原始 C/Python 源码只读引用；补齐基础、fused、高次、坐标、source/target 学习、图和秩相关文件，共 38 个产物。完整构建在临时树中结束后逐文件原子发布到 `build/production`。FLINT 版本从安装头文件读取，不改变已有 FLINT 安装。

复现入口：

```sh
python j_integral_bootstrap/tests/build_kernels.py --production
```

原生默认 persistent 分支并非在本 DP 输入上可用：默认提前支撑分支报 `early support dependence`；关闭提前支撑后，同一 worker 又以 SIGFPE 退出。证据分别保留在 `results/dp_official/runtime_default_failure.log`、`runtime_noearly_failure.log`。没有把错误当作候选不相容，也没有修改上游 C 算法。

简单、已验证的可用分支是关闭 persistent 和提前支撑，使用原包非 persistent fused 构造。`reference_space.py` 在完整产物存在时默认选择该分支：

```text
ANALYTIC_BOOTSTRAP_KERNEL_DIR = .../build/production
ANALYTIC_BOOTSTRAP_PERSISTENT_PURE = 0
ANALYTIC_BOOTSTRAP_EARLY_SUPPORT = 0
```

显式用户环境优先；生产文件不完整时退回既有基础目录。没有覆盖正在执行的基础核。

`compare_runtime.py` 对归档失败 6×66 块和真实成功三行块比较了**每个因子的完整支撑内空间**。使用原包 `FactorSpaceReference` 做精确空间包含检验，不仅比较参数数目。

| 实际对象 | 基础参数数 | fused 参数数 | 完整因子空间 | fused 构造秒数 | 真实数值结论 |
| --- | ---: | ---: | --- | ---: | --- |
| 失败 6×66 块 | 7175 | 7175 | 63 个全部相同 | 2.588 | 第 32 点仍证明不相容 |
| 已完成三行 sector | 137 | 137 | 16 个全部相同 | 0.126 | 33 点唯一相容 |

此前基础路径同一大块构造约 5.8 秒；这个比较是实际替换实现的收益，不是 AB 数学候选发生变化。详细结果为 `runtime_comparison.json`。

## 源阶、目标前缀并不减少最终未知量

原包有两种不同的有限阶学习。

源阶学习比较完整锚定源模，寻找在当前有限参考阶 $R$ 内得到同一个模的较短前缀。目标前缀则比较源模与指定 reduction 支撑相交后的空间。后者可能更早稳定，`learn_fused` 将两种测量放在同一学习会话，复用源分子模、相邻阶结果与图计算。

重要的程序约束是：目标前缀只适用于 reduction，要求 top 主部为零。`require_reduction_prefix` 明确拒绝 DE。不能把 DE 每一行当独立 reduction 使用此前缀，因为 DE 的非零 $A_{11}$ 让不同 top 行通过共享 Hom 参数相联。

这两种学习都要求保持最终候选空间相同，因此其直接收益是少做源模递推，**不是**凭空减少最终参数数目或所需 $t$ 点。它们不能修复不包含真实答案的分母/次数上限。

为验证它在当前项目上的作用，`probe_reduction_prefix.py` 使用实际已完成目标 `J{1,0,1,0,1,1,0,3}`、其 15 维闭合下层 DE 及已重建系数的精确分母、分子次数。原包 `learn_profiles` 实际调用了 `learn_fused`：

- 普通构造 0.079 秒；学习 0.089 秒；带 profile 的生产构造 0.064 秒。
- 24 份源分子模缓存全部消费，峰值约 22 KB；说明融合学习确实运行，不是只调用空接口。
- 12 个因子的目标前缀与默认阶相同，完整空间逐因子相等，最终均为 38 个参数。

这个实际局部 anchor 测试没有节省点数，单次构造也不能收回额外学习成本。对更高 dot/更稀疏支撑可再测试，但目前不能声称此机制能解决 DP 剩余未知系数。结果在 `reduction_prefix.json`。

## 已有部分有理条目能替代哪一部分 anchor

完整官方 `learn_anchor` 的 anchor 是固定其它参数与质数后、保留单变量的完整有理结果，不是一个数值 $t$。但其低层组件不要求整个问题完成：只要一个目标行或一个闭合 sector 的结果已重建，就能用该局部结果生成精确分母、次数上限、零掩码及其 profile；上述实际 reduction 测试已验证这个用法。

对于尚未完成的行，只有部分条目已知时，不能把其它位置填零交给 `learn_anchor`，也不能把已知条目的极点阶直接当作其它条目的确定阶。原包分母传输模板记录的是某个**已知目标**相对于已知连接因子签名的幂次，再运输到其它质数或运动学 setting；它没有“已知一个条目就推出同一行所有未知条目分母”的定理。

因此当前可立即利用 41→49 点新增数据的重点仍是：先重建更多条目，缩小每个实际条目的 cap 与分母，再把这些恒等式放进联合拟合。这个操作会降低剩余未知维数与需求点数。按行/分量保存已验证的幂次是实质信息；用所有已知条目的全局分母乘积覆盖所有位置，反而会膨胀问题。原包自动模板可以复用在已经完成的局部目标，不能替代未知包络的验证。

## support 与 graph 计划的优先级

`support_profile` 保存允许/禁止支撑消元的索引与主元计划；数值矩阵必须在当前输入重算并验证，失败回到普通求解。`GraphLearning` 记录源模闭包的 pivot decisions；它甚至只在存在足够重的 closure 调用时才实际挑选重放计划。两者都不改变解析候选，也不减少拟合参数或点数。

当前主要是在**同一个质数的一份首解**上不断修改未知包络。包络、支撑和 lower DE 随块变化，计划不能无条件跨候选复用。在当前耗时量级上，为单次 2.6 秒构造先做一次专门学习不优先；更值得在同一局部 anchor 完成后、多个质数/运动学点重复运行时接入。

此外 GraphLearning 要求 persistent graph worker；本输入的普通 persistent worker 尚有可复现故障，故不能为了套用图计划而重新打开未经验证的分支。完整产物已经可构建，今后可以限定小案例单独修复，但不阻塞当前 DP 完整结果。
