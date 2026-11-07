# 原包一次完整 solve 与逐 sector 调用的点数审计

本审计只检查本地已保存的 production checks、数值输出文件清单及 routing 事件，并运行一个轻量二 sector 数学模型。没有请求新的 DP IBP，没有读取 fresh R，也未修改 `official_route/run_upstream.py`。

## 原包实际行为

`PointContext.solve_de` 本身已经按依赖顺序迭代全部 sector。成功处理各 sector 后，在函数末尾调用 `_check_complete_de`。当 `shared_de_samples=True`、`global_holdouts=2`、`holdouts=2` 时，即使指定内部参数 `_sector_limit=k` 而仅处理前缀，这次调用也会进行两个新的前缀检查。

完成全部 sector 后 `self.de` 被缓存；此后的 `solve_de()` 直接返回缓存，不会重新取点。问题因此是第一次逐 sector 的 149 次前缀推进，而不是完成后数次读取结果。

旧驱动用 `k=1,...,149` 逐次调用内部 `_sector_limit`。前 148 次各产生两个额外前缀检查点，共 296 个；最后一次执行真正完整 DE 的两个全局检查。一次普通 `solve_de()` 只需最后那次全局检查。原包行为本身与其“内部恢复边界”语义一致，是驱动把恢复接口当作进度回调造成了额外工作。

## 不用 mock 替代求解的二 sector 回归

`tests/test_upstream_single_solve.py` 实际运行原包 `learn_anchor(learn_profiles=True)`、`prepare_setting`、PL 构造、拟合及检查。模型为

$$
S'=0,\qquad T'=\frac{2x}{x^2+1}S+\frac{2x}{x^2+1}T,
$$

固定 $p=1000003$。anchor 返回该模型完整有理 DE；数值 backend 只在请求的 x 上返回点值，每个 x 的完整评价缓存一次。仅对 `_check_complete_de` 计数，不替换它的实现。

| 调用方式 | 完整/前缀检查函数调用 | 前缀数值点 | 不同数值点总数 |
| --- | ---: | ---: | ---: |
| 一次 `solve_de()` | 1 | 0 | 6 |
| 每个 sector 调一次 `_sector_limit` | 2 | 2 | 8 |

两个结果的有理矩阵均与模型严格一致，sector 的训练和检查记录完全相同。一次求解的最终点为 145229、114498；循环的第一对前缀点恰为这两个，而循环的最终点变成 492800、12240。额外前缀检查污染了全局检查所用的 seen 集，导致最终随机点改变。重复读取已经完成的 `solve_de()` 没有任何新数值调用。

三个 unittest 全部通过。模型完整记录在 `results/dp_diagnosis/upstream_single_solve_model.json`，测试日志为 `upstream_single_solve.log`。

## DP 本地真实 checks 的独立审计

输入为 `results/dp_original_warm/production/connection.json` 的 149 条检查记录，不是由“每块最多多少点”估计。对所有 sector 的点集合做并集：

| 集合 | 不同点数 |
| --- | ---: |
| 所有 `training_points` 的并集 | 16 |
| 训练点与 sector 检查点的并集 | 18 |
| 与所有 sector 点均独立的最终全局检查 | 2 |
| 旧运行真正参与主生产任务的点总数 | 20 |

sector 检查点自身的并集也有 18 个，其中 16 个在其他 sector 被用作训练；不能把所有 sector 的点数相加。20 是 anchor/profile 已给定后该 production 的不同主取样点数，不包含取得 anchor 的成本，也不等于 backend 原始 evaluate 调用次数。

再逐个解析本地保存的 422 个成功数值输出路径：

$$422=20+106+292+4.$$

- 20 个完整输出，对应上述全部主点；20 个主点都有完整 family 的输出。
- 106 个子族输出，仍落在这 20 个主点上，是路由层在同一参数值上的额外子族计算。
- 292 个额外完整输出、4 个额外子族输出，对应 296 个额外不同参数点。

因此成功输出对应的不同参数值为 $20+296=316$；316 不是拟合需要的主点数。文件清单的成功输出数恰为 connection 中记录的 `native_new_points=422`。

routing 事件日志只记录子族及 fallback 情形，并不逐次记录所有直接完整求值；因此工具明确把事件覆盖标为 partial，不声称仅凭该日志逐个证明全部 296 个额外点的任务标签。296 与旧驱动 $148\times2$ 的实际调用语义及文件清单的额外不同点数完全一致。

## 修正驱动后不能复用旧最终点来冒充一次新运行

本地 routing 事件记录的第一对前缀参数值为 $t=945629,233564$。它们与 18 个 block 取样点不重合；没有中间前缀污染时，最终检查的随机流将首先尝试这两个参数值。旧运行的最终检查点不同。

所以“正常只需 20 个点”是点数结论，不能据此直接拿旧运行的 20 个完整输出声称修正后的运行无需任何新数值工作。新的最终全局检查必须按修正后原包的实际请求完成。尤其第一对前缀点原来只有子族输出，不能冒充新的完整 DE 检查数据。

所有 audit JSON 中的 `*_tau` 是原包求值坐标，`*_t` 是绝对参数 $t=\tau+2/5$；这里的 945629、233564 是绝对 t，不能与 tau 字段混用。

## 复现

```bash
python -m unittest discover -s j_integral_bootstrap/tests -p test_upstream_single_solve.py -v
python j_integral_bootstrap/diagnosis/audit_upstream_sampling.py \
  j_integral_bootstrap/results/dp_original_warm/production/connection.json \
  --numeric-file-list j_integral_bootstrap/results/dp_original_warm/production/numeric_file_paths.txt \
  --events j_integral_bootstrap/results/dp_original_warm/production/restricted_events.jsonl \
  --output j_integral_bootstrap/results/dp_diagnosis/upstream_sampling_audit.json
```

报告保存了 source connection 的 SHA256、实际点集合和输出分类；不改变已有结果。
