# 按 sector 缩小原生数值 IBP 接口

实现为 `restricted_numeric.SectorRestrictedNumericBackend`，继承现有完整 `NumericBackend`，不改原 AB 包、不改 `upstream_j.py` 与 `run_upstream.py`。

## 输入与映射

原包仍请求全部 master 列和当前 sector 的导数行。后端取这些目标的传播子支持并集 $S$，保留原独立基底中支持包含于 $S$ 的对象，按原传播子顺序压缩指数。保留对象的所有 J parent 都写入 targets，边界标签 `Bz0`、`Bz1` 保持端点身份。family 的运动学、维数、端点和三条 branch 约定全部沿用输入。

运行同一个 `global_point`，并要求保留子基的**所有导数**闭合；它同时检查该基底独立、所需 J/边界对象可由该基底表示。成功后按标签反向映射，零填回原包请求的全部列，添加原有规范因子的对角 shift，并保持请求 digest 与原列顺序。

还可以传入 `closure_by_task=closure_tasks(bundle)`。此函数仅提取原 bundle 的 sector master 名称和 `de_templates[task]['lower_indices']`，不会取用任何未知耦合 $C$ 的系数。这时保留基底是此闭包指定的对象，传播子支持也包含闭包中边界代表的支持。这样允许边界对称代表位于 top 的物理支持之外。压缩后的基底仍由原生程序逐点检验，不能因已有结构声明而跳过数值闭合验证。

## 明确的回退与错误语义

以下情况记录到 `numeric/restricted_events.jsonl`，然后调用原完整后端：

- 目标或结构闭包已经使用所有传播子；
- 不能同时保留 branch 0、1、2；
- 输入含不支持的传播子 cut、额外指数布局；
- 原生程序明确报告保留基底 dependent 或不张成目标/导数；
- 域、取样点与基底身份均正确，但保留导数不闭合。

任意其它原生失败、异常退出、JSON/域/标签错误会继续抛出，不能伪装成数学不相容或悄悄接受零。全族可能把低 sector 积分约化到另一个支持上的对称代表，因此“按支持保留的基底”并不总是完备；本实现对此安全回退，不擅自增加新基或改变原包的 basis。

## 实际验证

测试完全不读取 `dp_fresh_verification`。真实对照来自 DP 原学习完整点及 e6 原完整点，取 $t=17,18$。

- `test_restricted_numeric.py`：DP/e6 共 30 个物理支持，每个 2 点，共 51750 个 full 列分量，全同；测试还反转请求列序，并使用 gauge=3。
- `test_restricted_closure.py`：DP 支持 `[0,1,2,3,6,7]` 的完整 53 行 ×496 列 ×2 点，共 52576 项，全同。FBI 秩 16，实际新计算单点约 0.49/0.39 秒。
- 同一测试使用 e6 原包真实 bundle 的 block12 与 block42 结构闭包，共 948 项全同。block42 保留 33 维、6 个传播子，实际单点约 0.21/0.16 秒；block12 保留 6 维、4 个传播子，约 0.02 秒。
- 其它可用 e6 五传播子子族保留 4—14 维，实际新点约 0.007—0.06 秒。
- 更大的不完备支持确实触发 full 回退，其回填结果也逐项通过，未把失败当作成功加速。

机器记录：`results/dp_official/restricted_backend/{regression,closure_regression}.json`。新计算与缓存复用的耗时必须区别；后续重跑低 sector 已缓存时不能把缓存读取耗时称为原生求解时间。

## 调用

```python
from restricted_numeric import SectorRestrictedNumericBackend, closure_tasks
numeric = SectorRestrictedNumericBackend(
    family_path, basis_path, parent_targets_path, binary_path, work_path,
    point_directories=full_point_directories, gauge=3,
    closure_by_task=closure_tasks(bundle),  # 可选，仅名称和lower_indices
)
```

也可直接使用新入口 `official_route/run_restricted.py`，参数与原 `run_upstream.py` 相同；加 `--closure-from-bundle` 时，它在原学习/加载 bundle 完成后提取上述结构。未加该选项则只按请求目标支持限制。原学习流程、对角输入、`prepare_setting`、holdout 与完整验证逻辑不变。
