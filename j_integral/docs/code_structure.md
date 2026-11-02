# 程序结构与数据契约

## 1. 生产程序：C++

| 文件 | 职责 |
|---|---|
| `native/prepare.hpp` | JSON 家族输入；有限域 U、R、sector 分类、已知 FBI 递推；重定义和连接；相容性检查；固定多项式递推输入 |
| `native/series.cpp` | 有限宽度多项式 LRR；循环缓存；FLINT 多项式运算；固定极点分层积分和证书；正则端点 BFI 展开 |
| `native/search.cpp` | 多坐标有限域关系搜索；训练与留出检查；有理函数约化；候选基和 DE 搜索 |
| `native/tests.cpp` | 分治进制转换、矩阵方向、端点、奇异消元余项、vac 精确 DE 恒等式测试 |
| `native/check_series.cpp` | 与保存的慢原型系数作精确有限域比较 |
| `native/Makefile` | 构建；生产代码不依赖 Python |
| `examples/` | 拓扑、运动学、维度、prime、展开点、目标和端点约定 |
| `results/` | 可再生运行输出，默认不加入 Git |

C++ 构造的 FBI 基列满足 $f'=A f$，目标为行向量 $r f$。
总导数证书必须用 $g'+gA$，不能转置。拓扑及运动学只来自输入。
FBI 前处理可计算其已知连接，但从不读取或预先输入 $J$ 的约化/DE。

## 2. 输入与输出

`j_series INPUT OUTPUT --max-dot 2 --boundaries --order 200 --verify-all`：

- `two_loop_family` 输入直接给传播子、三个 branch、Gram 矩阵及质量；指标非负。
- 可输入已知的 `fbi_connection`，用于独立连接测试。
- `--max-dot D` 生成 topsector 中 $\sum_i(\nu_i-1)\le D$ 的全部指标；不设置时使用输入的 `targets`。
- `--boundaries` 为每个目标同时生成两端的 BFI 型对象；正则端点取值，`regulated_zero` 端点按输入约定为零。
- `--prime P --t0 A/B` 改变取样，无需修改引擎。
- 输出每个对象、每一阶在固定点积分/端点常数上的坐标，以及输入目标指标、计时和证书检查数。
- 默认对每个目标前五阶和最后一阶检查完整证书；`--verify-all` 检查全部阶数。

`j_search SERIES OUTPUT --degree 6 --train 120 --de-degree 16`：

- 搜索多项式关系并在更高阶逐项验证，再计算有理函数约化。
- 高 dot 对象优先消去，优先保留简单 $J$ 和 corner 边界。
- 零对象单独输出；master 列表不包含零边界。
- DE 从候选基展开的导数搜索得到，并保留验证阶数。
- 训练核的某些基向量失败时，求整个训练核与留出条件的交；结果报告失败基向量数，不能将使用过的留出阶数继续称为独立验证。
- `--j-only` 排除显式边界对象，用于比较 J-only 的闭合性。

## 3. 缓存和主要复杂度

生产循环不使用符号有理函数：保存 $F_n=H^n\Phi_n$、$V_n=H^{n+1}v_n$ 的多项式分子，仅保留递推需要的有限历史窗口。
系数积分用 $q$ 进制分解，高阶极点消元只更新邻层；消元映射按层缓存，跨目标复用。
当前仍有随展开阶数增长的多项式和极点层，尚不是线性标量成本。

## 4. 适用范围和失败处理

当前前端支持三个 branch、非负指标的二圈家族及可选 cut。后端要求固定点连接有限点 Fuchs 型、无穷远无多项式部分；不满足时明确报错。一般奇异端点的非零渐近系数尚未实现。
端点模式和所有有限域输入必须显式给出。证书失败、非普通点、除零和不能输出文件均须失败，不能静默生成结果。

## 5. 参考原型和 Git

`src/*.py`、`tools/run.py`、`tests/test_core.py` 是此前的慢版正确性原型；不作为生产入口。保留它们便于独立低阶回归。其依赖单列在 `requirements.txt`。
使用根仓库，不另建嵌套仓库；本模块之外的用户 staged/unstaged 修改保持不动。

## 6. 最后一重 FI 积分的初步实现

`native/fi_ibp.cpp` 读取 J 的有理约化和闭合 DE，设有理证书 ansatz 并比较多项式系数，精确检验 $r=c+g'+gA$，然后使用输入的端点映射和 period 归一化输出 FI 关系。端点映射必须显式给出；本轮只验证了等质量 vac，不自动推导一般家族端点。
`native/compare_blade.cpp` 与独立原始传播子 Blade 数据逐目标比较。完整测试定义及适用边界见 `vac_fi_ibp.md`。

## 7. 独立点检查与 FI 基选择

`native/verify_j_points.cpp` 不读取训练级数，直接在新 t 点构造 FBI z 连接并完成系数积分，以精确积分证书检查 J 的有理约化、DE 和零对象。当前验证入口限于 scalar FBI 连接；不将其描述为一般矩阵连接的完整符号证明。

`fi_ibp.cpp` 支持 `master_integrands="all_targets"` 和 `prefer_derivatives=true`，先消去全导数方向再选择剩余 FI 基；也可通过 `master_integrands_from` 读取外部给定的 dot-only master 指标。外部约化系数不是求解输入。

`compare_blade --allow-subset` 明确标识只核对 Blade 已计算的小批目标，不把小批结果声称为全部 FI 目标已被 Blade 检查。
