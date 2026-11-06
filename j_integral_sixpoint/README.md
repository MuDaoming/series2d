# 2412.21054 两圈六点例子的第一重积分

独立于 DP 的输入及结果目录。共用现有通用 C++ J-IBP、单变量重构及原 AB 包接口，不另写此拓扑专用约化程序。

数学定义、传播子和待验证范围见 [问题与结果](docs/problem_solution.md)。输入为 `inputs/family.json`、`inputs/targets_dot2.json`、`inputs/basis.json`。

完整有限域连接在 `results/connection.json`，最终八个独立点检查在 `results/verification_report.json`。另运行原包 `learn_anchor → prepare_setting → solve_de`，结果及日志在 `results/upstream/`，两条路线必须分别报告。
