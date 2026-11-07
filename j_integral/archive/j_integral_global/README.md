# DP 第一重 J 的全局约化

继续使用原 t,z 参数化，处理未 cut FBI 系统的低 sector 耦合。本目录输入、程序、诊断和结果独立保存。

当前已在有限域点处理全部 147 个普通三 branch sector、2952 个 dot≤2 的 J。纳入 FBI 恒等对象与物理端点关系后，候选基为 **481 个 J＋15 个边界函数，共 496 个**，单点 DE 闭合。**完整 DP 的 t 有理重构尚未完成，496 暂不称为最终主积分数。**

定义和推导见 [problem_solution](docs/problem_solution.md)，验证状态及区别于旧诊断的说明见 [testing](docs/testing.md)，生产并行的实际审查见 [workload_review](docs/workload_review.md)。

```bash
make -C j_integral_global/native -j2 all test
j_integral_global/build/global_point j_integral_global/examples/dp_physical.json \
  j_integral_global/results/production_probe/point17.json --point 17 \
  --targets-file j_integral_global/examples/dp_targets_dot2.json \
  --basis-file j_integral_global/examples/dp496_basis.json --values 1 --compact 1
```

`dp_physical.json` 显式启用本例已有解析证明的物理端点条件。一般 FBI 解使用 `dp.json`，其 596 维闭合表示包含额外端点自由度，不能当成物理最小基。

完整重构入口为 `tests/run_reconstruction.py`：每个 C++ worker 只负责一个 t 点，外层动态调度、冻结输入和执行文件、保存每点退出码与闭合验收，再调用 C++ 有理重构。完整 DP 生产取样拟在 Tianhe 执行；当前 SSH 入口不可达，未提交作业。不要把“可运行入口”误认为生产结果已完成。

较小的六传播子未 cut 家族已完成 44 维物理 J 系统的 t 重构和独立新点检验，记录见测试文档。
