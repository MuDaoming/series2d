# DP topsector：修正后的计数与验证

## 1. 当前结论

补齐系数积分消元核产生的关系后，当前 DP topsector cut 得到 **8 个 J/边界主对象：7 个 J 加 1 个正则端点函数**。

$$
\boxed{\{J_{11111122},J_{11111113},J_{11112111},J_{11111211},
J_{11111121},J_{11111112},J_{11111111},B^{(z)}_{0,11111111}\}}.
$$

这不是原始 FI 的 master 数。旧一般变域 F(delta) 家族的 43 个对象（15 FI、24 BFI、4 BBFI）只是另一家族的比较基准。

## 2. 输入和目标集合

输入 `examples/dp_topcut.json` 的传播子、质量、Gram 矩阵来自 `mmaexample/dp.nb` 的 DP 定义，三个 branch 大小为 $(3,3,2)$，$d=26/7$，cut 为 11111111。

主测试为全部 $\nu_i\ge1$、$\sum_i(\nu_i-1)\le2$，共 45 个 J；每个目标再生成两个 z 端点对象，共 90 个边界对象。零对象单列，不计入主函数。

z=0 的正则端点显式保留；z=1 使用已声明的解析零边界约定。J 定义含完整 Jacobian、分支权重、U 幂和 FBI 归一化。

## 3. 展开、搜索和两组取样

| 参数 | 主测试 | 第二取样 |
|---|---:|---:|
| prime | 1000003 | 1000033 |
| t0 | 2/5 | 3/7 |
| 展开到 | 200 | 220 |
| 训练到 | 120 | 140 |
| 关系留出 | 121–200 | 141–220 |
| 多项式关系次数上限 | 8 | 8 |
| J/边界主对象数 | 8 | 8 |
| DE | 全部闭合 | 全部闭合 |

两组得到相同基。主测试训练关系拒绝数为 0，DE 验证到 199 阶；第二取样验证到 219 阶。每个展开系数的有理积分证书都执行精确检查，并在输出中保留额外的 period completion 证书。

修正后的 DE 最高关系次数为 **37**，八行的搜索次数为 37、37、23、23、23、23、21、5（与基的排列对应时以实际输出为准）。这些是多项式关系次数，不是微分方程导数阶数。

## 4. 不使用训练级数的独立点检查

`native/verify_j_points.cpp` 重新从原始输入构造已知 FBI 关系，直接在新的 t 值计算完整 z 系数积分及其精确总导数证书；它不读取训练点的 Taylor 系数。

两个 prime 各检验以下 8 点：

$$
1/7,2/7,3/7,4/7,5/7,6/7,2/9,4/9.
$$

每点通过 82 条有理约化检查、8 条 DE 检查、45 个零对象检查；每点所需 52 个直接积分证书均通过。两组报告为：

- `results/dp_topcut/independent_point_checks.json`
- `results/dp_topcut/independent_second_point_checks.json`

这些是独立点的精确有限域证书，不等于对一般运动学或一般 t 的完整符号证明。原始 FI 的独立 Blade cut 对照单独记录，不将上述检查冒充 Blade 验证。

## 5. 为什么此前的 9 个被修正为 8 个

旧程序保留了奇异消元映射的余项，却未利用其核方向产生的低次精确关系。本例无穷远 n=2 的核方向给出 g=z^2；最高次项抵消后，剩余低次项的完整积分与端点值之间还有一条关系。

输出 `period_completion.certificates` 保存生成函数、分母幂次和原始关系。程序先精确验证该关系，再统一化简每个系数的固定点坐标。加入该关系后，原来的 9 个候选中有一个冗余；先前 9 维、最高次数 46 的记录只适用于冗余表示，不能再引用为最终计数。

修正并未改变 FI/J 的数学定义，也不是通过删掉某个积分坐标强行降维；每次消去都有精确总导数和端点证书。

## 6. 新旧 DE 的比较

| 指标 | 旧 F(delta) | 修正后 J/边界 |
|---|---:|---:|
| 维数 | 43 | 8 |
| 逐行通分后的最高关系次数 | 79 | 37 |

旧值来自 `search/runs/dp/bc1/de/Ax.wl` 在相应有限域约分后逐行通分统计；旧搜索上限 m=100、展开阶数 6000 均不是实际 DE 次数。

## 7. 复现

```bash
j_integral/build/j_series j_integral/examples/dp_topcut.json \
  /tmp/dp.json --max-dot 2 --boundaries --order 200 --verify-all
j_integral/build/j_search /tmp/dp.json /tmp/dp_search.json \
  --degree 8 --train 120 --de-degree 60
make -C j_integral/native ../build/verify_j_points
j_integral/build/verify_j_points j_integral/examples/dp_topcut.json \
  /tmp/dp_search.json /tmp/dp_points.json
```

## 8. 原始 FI 的 Blade cut 对照

`tools/check_dp_blade_cut.wl` 只处理原始 FI 的 maximal cut，使用 IntegralOrdering=2 和 PreferredMI=Dot；物理传播子之外补齐 3 个 ISP，但最终基要求没有 rank。
最初小批仅含 corner 和 8 个单 dot；方程生成后，复用全部 8100 条 maximal-cut IBP、用 Blade 相同的 FiniteFlow 求解器独立求出全部 45 个 dot≤2 目标。最终原始 FI 基为 15 个 dot-only masters，全部 45 个目标与 J 路线的 FI 约化模 p 一致（其中 30 个为非 master 目标）。未启动完整未 cut 的全家族约化。详见 `dp_fi_ibp.md`。
