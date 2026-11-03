# 部分分支积分 J(t) 的约化试验

从已知 FBI 关系直接生成完整积分掉一个分支变量后的 $J(t)$ 展开，搜索约化和 DE；不预先输入 $J$ 的 DE。

旧 δ 展开与当前 t 展开的逐步推导及区别见 [delta_vs_t_expansion](docs/delta_vs_t_expansion.md)。

数学定义见 [problem_solution](docs/problem_solution.md)，程序结构见 [code_structure](docs/code_structure.md)，实际结果见 [testing](docs/testing.md)。

生产程序为 C++17，依赖 FLINT、GMP、nlohmann/json 头文件。当前机器已安装。

```bash
make -C j_integral/native all
make -C j_integral/native test

j_integral/build/j_series j_integral/examples/db_topcut.json \
  /tmp/db_series.json --max-dot 2 --boundaries --order 200 --verify-all
j_integral/build/j_search /tmp/db_series.json /tmp/db_reduction.json \
  --degree 6 --train 120 --de-degree 16
```

`--max-dot 2` 是全部总 dot 不超过 2 的 topsector 指标，DB 共 28 个 J；`--boundaries` 再加入 56 个端点对象。输入中的运动学、维度、展开点和 cut 均未写入引擎。

补齐固定点积分的核方向关系后，DB topsector cut 的 J 层为 3 个 J 加 1 个正则端点函数，共 4 维。此前的五维表示有一个冗余。旧 Python 程序只保留作低阶回归。

最后一重积分的初步结果：等质量 vac 的全部 dot≤3 共 20 个 FI，利用 J 的 DE 做 IBP 后，与原始传播子 Blade 约化模 p 全部一致。推导和边界项见 [vac_fi_ibp](docs/vac_fi_ibp.md)。回归命令：

```bash
bash j_integral/tests/test_fi_ibp.sh
```

DP 更新后的计数和验证见 [dp_test](docs/dp_test.md)：J/边界为 8 个，原始 FI 的 cut 为 15 个，两者不同。全部 45 个 dot≤2 的 FI 与 Blade 生成的 maximal-cut IBP 独立解模 p 全部一致，包含 30 个非 master 目标；详见 [dp_fi_ibp](docs/dp_fi_ibp.md)。

```bash
bash j_integral/tests/test_dp_fi_ibp.sh
```
