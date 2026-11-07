# DP 原包路线的天河准备

这里不修改 `extern/analytic_bootstrap-master`。实际调用链为：

1. `run_case --stage anchor` → `cold_anchor` → 通用 `collect_points` → 原 C++ J-IBP。只对固定主基的 485 个不同父指数作积分，计算 496 个导数；不再顺带计算全部 2952 个 dot≤2 目标。由数值点直接取得一次完整的有理 anchor，不使用 AB 候选猜测。
2. `run_case --stage production` → 原 `learn_anchor` → 原 `prepare_setting` → 原 `PointContext.solve_de`。保留原包的 profile 学习、由下至上的 sector 注册、数值增量求解和保留点检查。只有 IBP、anchor、对角块三种输入能力由本项目适配。
3. 物理连接按原包 `SolvedDE.context['masters']` 的真实行列顺序导出，并撤销显式共同归一化。逐条目与 anchor 作有理函数恒等比较，比较通过才写阶段成功凭证。

原包 `learn_anchor` 的契约要求首份完整单变量有理 DE；它不是从一个数值 $t$ 点直接学出所有分母。`cold_anchor` 替代原包中由 Blade 负责的这一次标准重构，不改其后续学习与生产算法。对单个固定参数点而言，这一次 anchor 的成本必须计入，不能说 AB 已免费省去了它。

## 实际调度审查

每个 $t$ 点是独立数值单元。当前需求批次的点由真实 `ThreadPoolExecutor` 动态派发，短任务结束后立即领取下一点；只有出现未重构系数时才增加下一批点，避免盲目计算到最大点数。`test_collect_dispatch.py` 实际执行该收集器，验证第三点在较慢第一点结束前启动、输出唯一、再次运行不重算。

拟议单节点 anchor 作业申请 48 个单核槽位、131072 MiB；每点按 2048 MiB 预算，48 点并发约 98304 MiB，剩余给协调和 IO。并发受实际 CPU affinity 检查，不能在只有一个可用核的 launcher 下启动 48 个 worker。所有软件线程数为 1，纯 C++/Python，不启动 Mathematica，因此不套用 Mathematica 的启动错峰。

anchor 完整通过后释放其分配；再提交原包生产作业，申请一个核、65536 MiB。原包单 setting 的逐层依赖是串行的，不保留 48 核资源等待这个阶段，也不拆成多个长期低利用节点。当前没有两个同时需要占用的计算节点。每个阶段记录 job、主机、PID、affinity、输入/脚本哈希、起止、程序退出码、验证退出码及科学完成状态。

`submission_manifest.json` 已通过技能机械预检；实际调度审查的证据是上述调用链、源文件和动态派发回归，不以机械 pass 代替审查。

## 上线前必须完成的环境核实

用户重新连接后，SSH 已恢复。已部署到隔离目录 `/fs2/home/mdm/Projects/series2d/j_ab_production_20261107`，25 个实际源码文件哈希一致，私人 python-flint 0.9.0 加载及三项适配器回归通过。anchor 作业 `12288434` 已提交，原包生产作业 `12288510` 以 `afterok:12288434` 排队；观测时分别为 Priority 和 Dependency，尚未进入数值计算。

恢复连接后先用既有 `pipeline/tools/remote_tianhe.py` 核实节点、当前分配/队列、配额及私人运行环境。远端系统 Python 为3.6.8，已选择 `/fs2/software/python/3.10_manbaforge_4.14.0-0/bin/python3`。私人 FLINT2.9 位于 `/fs2/home/mdm/.local`；nlohmann 头文件和 manylinux cp310-abi3 的 python-flint 已隔离部署。GCC8.5 的 C++17 filesystem 在构建时显式链接 `libstdc++fs`。不得把本地 CPython 扩展或 Linux 动态二进制直接当作可用的天河版本。

提交前在远端原生构建/核实 C++ worker、重构器及 AB 原生内核，确认输入、完整源码和运行库已部署且哈希一致。`tests/build_kernels.py --production` 提供本机原包运行库构建入口；现在支持 `--prefix` / `FLINT_PREFIX` 指定实际安装；SBATCH 明确要求 `AB_FLINT_PREFIX`，不把本机 `/usr/local` 套到天河。`build_backend.py` 和原子发布的原包运行库构建已在本机实际执行验证。检查日志父目录存在且可写。

再重新执行 skill 的 `workflow.py prepare`。两阶段均通过 `workflow.py submit` 提交；第二阶段的 `--dependency=afterok:12288434` 已写入资源清单并核对调度器。第一阶段退出0要求完整 anchor 科学标志；第二阶段还会再次检查完整输入，因此依赖等待不占计算节点，也不需人工等候后再提交。

## 已有实际验证

相同通用接口已在两圈真空例子跑通两个阶段；收集器、stage launcher、原包学习与生产、阶段凭证及 anchor 恒等比较均实际执行。

论文六点例子完整运行原包 47 层，79×79 物理连接与独立直接重构的全部 6241 项作为有理函数相同，另在全新的八个数值点比较 49928 个分量全部通过。相应记录在 `j_integral_sixpoint/results/`。这验证了接入方法，不等于天河环境已经可用或 DP 原包作业已经完成。

## 2026-11-07 恢复记录（取代上面的初次排队状态）

初次 anchor 作业12288434在7秒内因计算节点没有 `/usr/bin/g++` 而失败，未进入数值正文；其未启动的依赖作业12288510经身份核对后取消。故障已通过技能 incident 机制记录，未改公共软件。

新的 warm 作业 **12289927** 使用 `warm.sbatch` / `warm_manifest.json`，已在cn81132运行。共享GCC11.3对filesystem及FLINT的编译、链接、执行探针已经通过；当前继续构建私人运行库。复用独立新点验证通过的496维anchor，SHA256为 `dafcbfcb4ba7a3af89daf32823f83281d709e673d2fd5f9035291d60994a657b`，不重复冷重构。新工作目录为 `results/dp_original_warm`。

生产入口选 `run_restricted.py --closure-from-bundle`，只替换数值IBP接口；原包学习与自下而上求解不变。子族接口的105274个真实分量已与全族结果逐项核对；不满足子族闭合时记录原因并回退全族计算。原版流程全部结束并与anchor逐项有理恒等比较通过，才会认定生产完成。

## 已完成与当前后续任务

恢复作业 **12289940** 已 `COMPLETED 0:0`，149块全部完成，与独立完整连接的246016矩阵元逐项有理恒等检查全部通过；整作业3:50:01，峰值RSS3116192KiB。私人cc别名问题已按技能记录并通过实际生产验证。旧12289927只完成C++构建、未做科学计算，不代表成功生产；12289940复用了那些有效二进制。

当前补候选J约化的作业是 **12296043**，`remaining_r.sbatch`：R-only全族1022父指数，400个数值点，48动态单核worker。其科学验收要求所有400点有有效回执及完整1252对象输出；取样完成之后还需有理重构及独立点检查，不能以队列状态宣布完整约化成功。
