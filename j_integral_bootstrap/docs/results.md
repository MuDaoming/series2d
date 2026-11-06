# 第一轮测试记录

## 1. 范围与比较口径

本轮先验证上游算法与 J 系统的兼容性。分母包络、分子次数上界、已知对角块
来自旧结果；因此下面不是 496 维 DP 全局系统的从头重建。

基准包括：二圈 vacuum、DB 与 DP 的 topsector cut，以及从 DP 提取的
六传播子小型全局例子（44 个 J 主函数，保留所有低层耦合）。
vacuum 的两个主对象及 DB/DP cut 中的边界对象沿用旧测试定义，不能把
旧 cut 的主函数数目当作完整物理全局主函数数目。

只比较相同包络下的线性分子系数参数。论文在其他例子中的收益不保证
在本例达到同样倍数；本轮没有加入额外无穷远区域约束。

## 2. 已知连接下的约化

| 系统 | 非主目标数 | 原参数 | Bootstrap 参数 | 训练点 | 留出点 |
| --- | ---: | ---: | ---: | ---: | ---: |
| vacuum | 14 | 91 | 61 | 4 | 8 |
| DB top cut | 37 | 260 | 170 | 3 | 8 |
| DP top cut | 82 | 3651 | 986 | 5 | 8 |
| 六传播子全局 44 维 | 284 | 23642 | 6683 | 6 | 8 |

四例都逐项通过精确有理函数恒等式核对。前三例的初始取样来自参考函数求值。
DP cut 在第二质数 $p=1000033$ 下复验，仍为 3651→986 个参数、5 个训练点，
精确有理函数比较通过。
44 维例子还使用现有 C++ 原生 J-IBP 的实际取样，所有 12496 个矩阵分量
（含结构零）通过精确比较。

44 维约化使用本轮逐块重新构造的完整连接，输出
`results/minimal44_reduction_native.json`，包含主函数、约化系数与连接。
其 PL 构造约 0.62 秒，拟合约 4.50 秒；取样文件复用同一轮 DE 验证的实际 C++ 数据。
这些时间不包含包络、对角块获取与新 C++ 取样的成本。

在相同已知分母条件下，普通逐系数多项式插值的最大分子系数数目分别为
7、4、15、37；所以应与这些数目比较，而不能把旧的无分母先验有理重构
所用 25、9、89、121 个训练点全部算作解析性约束本身的收益。

## 3. 不输入未知耦合块的逐层 DE 验证

44 维系统的非零依赖图有 19 个强连通块。每层只输入对角块和已完成低层连接，
构造、拟合本层未知耦合块；下一层的低层连接来自刚刚算出的结果。

原归一化下，在 $J_{\{1,0,1,1,1,0\}}$ 层发现真实系数不属于候选空间。
已用独立原生取样排除“参考结果与 C++ 不一致”，并用局部 jet、PL 两种路线复现。
逐因子诊断定位在 $t=1$，残数解释见 [算法文档](problem_solution.md)。

统一以 $h=t-1$ 变换已知的对角与低层连接后：

- 全部 19 个块通过，包括实际 C++ 数值取样。
- 未知耦合块的原参数 3509，约束后 1217。
- 各层最多 10 个训练点，各有 8 个留出点；不同层共用取样。
- 最终 44×44 连接与参考逐项精确相等。
- 结果保存于 `results/minimal44_sector_de_native_t1.json`。

两项 44 维测试共使用本轮生成的 18 个原生取样点，C++ 文件内记录的计算时间合计
约 11.80 秒，不含子进程启动；重复运行会复用这些文件，不能将缓存运行称作新取样。
原生取样身份文件包含程序二进制、family 输入 SHA256、质数及完整基底顺序。

把所有未知分母因子乘进规范变换的对照虽然也通过，但参数变为 3509→3509，
没有留下简化效果。该对照不作为采用方案。

## 4. 可复现命令

在项目根目录执行：

```bash
python3 j_integral_bootstrap/tests/build_kernels.py
python3 j_integral_bootstrap/tests/sector_de.py \
  j_integral_global/results/reconstruction/minimal44/reduction.json \
  j_integral_bootstrap/results/minimal44_sector_de_native_t1.json \
  --gauge 1 --gauge-policy t1 \
  --config j_integral_global/examples/minimal44_jonly.json
python3 j_integral_bootstrap/tests/reference_space.py \
  j_integral_global/results/reconstruction/minimal44/reduction.json \
  j_integral_bootstrap/results/minimal44_reduction_native.json \
  --config j_integral_global/examples/minimal44_jonly.json \
  --connection j_integral_bootstrap/results/minimal44_sector_de_native_t1.json
```

`python-flint==0.9.0` 使用本目录 `.deps/` 的隔离安装。外部源码和原生内核
目录不纳入主 Git；上游提交与 ZIP 指纹见 [来源记录](upstream.md)。

## 5. 当前尚未完成的环节

完整 DP 的 496 个主对象在已有两个普通点的 DE 非零图中分成 149 个强连通块，
最大块 9 维；J 主函数分布在 124 个 sector。这个统计是数值结构观察，
不是依赖图完备性证明。

maximal-cut 标量 FBI 的第一点测试用时约 1.52 秒：124 个 J sector 中，
99 个的 2024 个同 sector DE 分量全部与全局数值结果相同。
其余 25 个不能直接使用这条路径：8 个没有保留的 cut MFBI，15 个在去掉 cut
边界后原 J 集合不再独立，2 个原 J 集合不张成 cut 导数。程序明确报告未完成，
不把这些情形强行拼成全局对角输入。
在 $p=1000003,t=23$ 和第二质数 $p=1000033,t=19$ 两个设置复验，
同样的 99 个 sector、每点 2024 个标量分量全部通过，25 个未完成分类不变。

另做了按活跃传播子数优先的基底选择对照，在三个独立设置（含第二质数）
仍选出同一组 496 个主对象，因此这个差异不能归咎于原基底的排序。

已检查按规范 MFBI 的一重积分周期组织基底能否处理上述边界扩张。
任意多项式／有理周期试验已经否决：它没有保持原 J 的闭合子空间。
44 维例子被扩大成 60 或 69 个对象，DP 的多项式试验出现 607 个对象且
12 个导数不闭合。不能把这些数目称作新的 J 主积分统计，也不能据此增加生产基底。
源码保存在 `experiments/period_probe.cpp`，汇总为 `results/period_probe_summary.json`；
它是失败方案记录，不属于默认编译或生产入口。
完整 DP 仍需确定实际对角块、获取包络及次数上界、验证端点归一化、逐层求解耦合。
不能把上述 44 维验证称为 496 维全局 DE 已完成，也没有产生新的全局原始 FI 约化。
