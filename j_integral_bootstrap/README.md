# J(t) 的 Analyticity Bootstrap 测试

在原来的 $t,z$ 坐标下，复用已有数值 J-IBP 作为取样器，测试
arXiv:2610.01154 的解析性约束。上游 0.3.0 源码由用户提供 ZIP，独立保存在
`extern/analytic_bootstrap-master/`，不修改原始副本。

- [对象、算法与端点归一化](docs/problem_solution.md)
- [测试结果、比较口径与运行命令](docs/results.md)
- [上游身份和接口](docs/upstream.md)

已完成给定参考包络及对角块的验证：vacuum、DB/DP top cut，及六传播子
44 维全局 J 系统。44 维例子使用实际 C++ IBP 取样，重建全部未知耦合块，
并用所得 DE 完成 284 个非主目标的约化；结果逐项与独立旧结果精确比较。

这尚不是从头完成 496 维 DP 全局系统。`experiments/period_probe.cpp`
是被否决的任意 MFBI 周期基底实验，会扩大对象集合，不能作为物理 J 约化结果。
