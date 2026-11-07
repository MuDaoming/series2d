# 用原包 AB 的第二组接口测试 J 的 DE

两个例子：`dp/`、`e6/`。使用已有完整J DE作为学习anchor，新 `../src/IBP/build/j-ibp` 作为唯一数值IBP程序。均固定 $p=1000003$、$d=26/7$ 及输入中的精确运动学，只重构单变量依赖。

## 调用流程

`dp/run.sh` 与 `e6/run.sh` 中的数学求解命令直接来自原包：

```sh
python -m analytic_bootstrap.upstream learn \
  --config RUN/backend.json --definition input/definition.json \
  --setting input/setting.json --output RUN/bundle.json --timings RUN/learn_timings.json

python -m analytic_bootstrap.upstream solve-de \
  --config RUN/backend.json --bundle RUN/bundle.json \
  --setting input/setting.json --output RUN/solution.json --timings RUN/solve_timings.json
```

没有自写sector循环、部分恢复接口调用、AB候选构造、拟合器或取点调度。由原包的一次 `solve_de()` 完成其内部依赖遍历及最终检查。两个命令的先后顺序只来自学习结果供求解使用这一依赖。

`backend_worker.py` 只实现原包 `CommandBackend` 的JSON-lines请求与响应：

- anchor能力读取 `anchor.json.gz`，返回学习设置的完整有理矩阵；
- cut能力读取 `diagonals.json`，只返回所请求块的有理对角矩阵；
- numeric能力只读取FI/J定义、主基和归一化信息，调用新的 `src/IBP`，返回原包要求的行列顺序。它不读取anchor或非对角参考系数。

相同设置、相同数值 $t$ 的完整IBP结果在所有会话间复用。缓存从空目录开始；没有旧点seed，没有子族尝试或失败回退，也没有接口自行选择新数值点。

## 输入

每例的 `input/` 包含：

| 文件 | 内容 |
| --- | --- |
| `family_request.json` | 完整FI家族、Gram矩阵、质量、J的分支与端点、固定数值设置，以及新IBP的请求模板 |
| `basis.json` | 物理J/B主基及其顺序 |
| `definition.json` | 原包所需主对象分块、块之间依赖、变量名及归一化 |
| `setting.json` | 模质数及其它参数赋值；本测试其它参数已写入家族，故assignments为空 |
| `anchor.json.gz` | 给学习接口的完整有理DE，按原包规范顺序排列 |
| `diagonals.json` | 从同一设置已有DE抽出的对角块，生产过程不能由它读到非对角答案 |
| `normalization.json` | 物理参数与原包变量的转换、共同规范因子 |
| `provenance.json` | 原始FI/DE输入来源、指纹、秩及非零项数 |

`prepare_inputs.py` 只把已有数据整理为上述输入，不计算新的IBP或DE。块由已有DE依赖图的强连通分量确定，不把只有同一传播子sector的对象错误拆开；原包自身执行合法的依赖排序。

### 变量及共同规范

归档的有理函数以 $\tau=t-2/5$ 表示。原包的变量为 `tau`，它询问 $\tau=x$ 时，接口将新数值IBP的物理参数设为

$$
t=x+2/5\pmod p.
$$

这只是坐标平移，没有Taylor展开。沿用已有正确性测试的共同基变换

$$
\widetilde{M}(t)=[t(t-1)]^3 M(t),\qquad
\widetilde A=A+3\left(\frac1t+\frac1{t-1}\right)\mathbf1.
$$

anchor、对角块和数值响应都使用同一规范；物理J数值IBP仍由新程序直接计算。此变换不输入任何未知非对角耦合。输出需按 `normalization.json` 撤回共同因子后解释为物理J的DE。

本次学习与求解使用同一维数、运动学和质数，因此目标是验证正确接入及实际成本，不声称已经完成跨运动学或跨质数结构迁移。

## 执行

```sh
make -C j_integral/src/IBP
bash j_integral/tests/e6/run.sh my_run
bash j_integral/tests/dp/run.sh my_run
```

已有运行目录不覆盖。环境由 `runtime.sh` 指向现有原版 `analytic_bootstrap` 源码、python-flint及已构建原生库。没有改动原包源码；禁用的是已知有问题的可选persistent实现，使用其现有nonpersistent原生实现。所有内核线程数为1，单设置逐点取样使用原包默认的逐点接口。

## 记录及验证

原包自带的记录：

- `bundle.json`：学习结果；
- `learn_timings.json`、`solve_timings.json`：阶段与总耗时；
- `solution.json`：完整结果、每块实际训练点/保留点及最终独立检查点。

接口补充的记录：

- `run.json`：主机、起止时间、Python与包位置、关键环境、全部输入/接口/执行文件SHA256；
- `numeric/events.jsonl`：原包的请求及缓存命中，不新增求值；
- `numeric/receipts/`：每次真正调用数值IBP的退出码、耗时、CPU时间、输出SHA及闭合检查；
- `numeric/points/`、`numeric/logs/`：真实新IBP输入、结果和日志。

`verify_result.py` 在求解结束后作有理函数恒等比较，包括隐含零，并核对真实IBP输出数等于原包实际需要的不同数值点数、没有前缀检查任务。成功才写 `verification.json` 和 `scientific_complete=true`。数值工作和多项式恒等比较分别记录，不将实际IBP次数和后端响应次数混淆。

原始运行保存在各例的 `runs/`（Git忽略）；最终摘要和已验证DE会另外保存为可追踪参考结果。

`export_physical.py` 仅撤销已声明的共同规范，把物理J的DE保存为 `physical_de.json.gz`；不进行新拟合或新取样。

## 完成结果

| 例子 | 主对象 | 原包学习 | 原包求解 | 真实数值IBP | 失败 | 有理结果比较 |
| --- | ---: | ---: | ---: | ---: | ---: | --- |
| DP | 496 | 652.68秒 | 1311.97秒 | 20 | 0 | 全部246016项一致 |
| e6 | 79 | 5.00秒 | 18.07秒 | 10 | 0 | 全部6241项一致 |

DP包含最终核对总共约34分钟，e6约24秒。原包学习和求解的时间不包含此前取得anchor的成本；两例均复用已有anchor，数值取样从空缓存开始。撤回共同规范后，两例物理DE亦与各自参考完全一致。摘要为 `results_summary.json`，可追踪验证报告和物理DE在各例 `reference/`。
