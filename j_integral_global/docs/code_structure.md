# 程序结构

| 文件 | 作用 |
|---|---|
| `native/global_point.cpp` | 一个有限域 t 点的完整 FBI、所有 J、物理端点、选基、约化与 DE 点值；支持固定基、冻结候选及检查点恢复 |
| `native/triangular.hpp` | 一元有理函数、逐不可约因子 Hermite、完整三角耦合传递、Laurent-log 端点及物理边界商空间 |
| `native/reconstruct.cpp` | 从已完成且闭合的多点数据重构约化与 DE 系数，保留 8 个验证点；不做物理积分求值 |
| `native/audit_fbi.cpp` | 用一次双数导数检查完整多 MFBI 连接的相容性 |
| `native/probe.cpp` | sector、MFBI、共同分母、三角稀疏性和端点极点的前期成本探针 |
| `tests/run_reconstruction.py` | 外层动态调度；一个 C++ worker 对应一个 t 点；冻结输入/二进制、点锁、退出与闭合检查、恢复和自适应增加样本 |
| `tests/verify_reconstruction.py` | 八个新点从头验证完整多 MFBI 的约化和 DE；含故意错误系数的拒绝检查 |
| `tests/test_*.cpp`、`test_dispatch.py` | 数学与执行行为回归 |
| `results/summary.json` | 已完成结果、未完成项及输入指纹 |

共享的 `j_integral/native/prepare.hpp` 只增加可选的点值、一次双数和 FBI 恒等对象功能。原输入不启用新选项时保持原默认逻辑，不覆盖原 cut 测试结果。拓扑、质量和 Gram 矩阵均来自输入 JSON。

`examples/dp_physical.json` 明确启用本例已证明的物理边界条件；`examples/dp.json` 保留一般端点，用于对照。诊断阶段的旧候选及升 dot 输入不作为生产默认。
