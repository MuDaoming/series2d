# E6：固定运动学的全局 dot<=2 数值约化，比较 FI 与 J

本任务不建立 FI 的运动学 DE，不把任何物理不变量改成变量。完整家族及固定标量积沿用 [e6 输入](../input/family_request.json)，$d=26/7$。目标为七个物理传播子所有非空 sector、无 rank、$\sum_i\max(\nu_i-1,0)\le2$ 的1695个积分。六个辅助传播子的指标全部为0。

其中584个缺少某 branch 的目标含无标度子圈而为0；其余1111个保留三个 branch，是 J 接口可处理的共同集合。FI 实测全部1695个，没有只做topsector cut。1695个目标明确列在 `input/all_targets.json`，1111个共同目标列在 `input/common_targets.json`。

## 一次 IBP 的实际规模

FI 与 J 的规模测试都使用 $p=9223372036854775783$（FiniteFlow默认第0个质数）。

| 方法 | 实际使用的代数结构 |
| --- | --- |
| FI 数值 IBP | 生成267120条非零方程，129722个积分变量，4326889个非零系数；择行后120817条，裁剪未用方程后保留54419条 |
| 原坐标 J79 | FBI闭式公式涉及71次小矩阵消元，共227行，最大系数矩阵7×7；z方向流式提交94条非零关系，其中54条独立，涉及501个形式坐标 |
| 解析分离坐标 J61 | 同样71次FBI小矩阵消元；z方向69条非零关系，54条独立，涉及465个形式坐标 |

FBI小矩阵按形状分为24个1×1、19个2×2、7个3×3、1个4×4、20个7×7。它们是约束 Hessian/KKT 公式中的矩阵求逆或消元，不能叫作71次完整数值IBP。

**J没有生成与FI相同类型的Laporta大方程组，不能仅用某一个数字假装二者是同型矩阵。**其FBI阶段使用已知解析递推，在 $\mathbb F_p(z)$ 上处理上述小矩阵；z阶段使用Hermite消去极点及稀疏关系的流式消元，在 $\mathbb F_p$ 上约化其形式坐标。局部闭式递推和Hermite操作没有逐条展开为一个巨大的IBP方程表，故94/69只指实际提交到全局z/边界关系消元器的非零行，不是所有操作的总次数。

最后的表示式选择还有1195个向量（1111个目标加84个辅助端点向量），原坐标涉及368个形式坐标、秩79；分离坐标涉及312个、秩61。这些目标表示式不重复计作IBP关系。形式坐标的列数会随有限域因子分解及取样点变化，不能当作主积分数；所有采样点的范围也保存在comparison.json，z关系行数和独立秩在本次全部采样中保持不变。

FI的267120是送入求解流程的生成方程数，54419是FiniteFlow裁剪后保留的独立方程数；原始变量数129722不冒充裁剪后的列数，后者本次未单独记录。

## 数值 IBP 次数：直接重构和原包 AB

J多点重构/AB采用原有已知DE和学习数据相同的 $p=1000003$，变量为 $\tau=t-2/5$；物理维数、质量、运动学不变。这里与上节63位同质数的规模测量分开记录，不隐瞒两个有限域的区别。

| 方法 | 原坐标 J79 | 分离坐标 J61 |
| --- | ---: | ---: |
| 不利用分离公式的通用直接有理重构，实际IBP | 148 =144训练+4检验 | 14 =10训练+4检验 |
| 原包AB的约化阶段，新增IBP | 7 =5训练+2检验 | 6 =4训练+2检验 |
| 原包默认完整reduce调用，含DE阶段 | 17 =10 DE+7约化 | 10 =4 DE+6约化 |
| 本次获取完整约化学习锚点的成本 | 148（复用上述直接重构结果） | 1（一个数值点+已证明的分离公式） |
| 包含本次约化锚点成本 | 165 | 11 |

最后一行**不包含历史DE锚点的获取成本**；原坐标使用项目已有DE。分离坐标DE由解析幂函数直接构造，无需先做数值DE重构。不能将17或7称为从裸家族冷启动取得全部约化的总成本。

对这个特殊e6，已证明分离公式后，J61在任意t的完整约化只需**1次数值IBP**，无需再做通用有理重构或AB。因此14和10是刻意测试通用重构流程的实际成本，并非使用解析结构后的必要成本。

直接重构使用同一批完整1111行矩阵数据，逐项做未知分子/分母的平衡有理重构、以新增点检验，按批次增加采样数。报告实际调用数，不把事后知道次数后的最小插值点数冒充实耗。原矩阵最复杂元 $\deg N+\deg D=137$（69+68），分离矩阵为4。

AB调用原包 `analytic_bootstrap.upstream learn`、`reduce`，没有改动原包的学习、拟合、选点或sector调度。`backend_worker.py` 只是JSON-lines接口：读取声明的anchor/cut，或调用公共J-IBP程序。每个t的完整数值结果在原包各会话间复用；缓存命中、prepare请求和每块检查均不重复算成IBP。两组最终运行均没有坏点或失败数值IBP。

在需要DE数值检查的点，接口额外请求DE，first-jet会使用第二个FBI模型，小矩阵消元变为142次；**约化阶段只请求约化，不重复计算DE**，其规模就是上表71次。z关系数量不变。第一次接口诊断曾重复计算DE，现场单独保存在本地 `j_ab/diagnostic_all_derivatives/`，不混入最终运行的17次计数。

## 主积分计数与正确性

- FI数值IBP使用已有Blade的62个代表；其中一对完全相同的massless sunrise可以显式合并，物理去重后61。
- 原J系统79维（含一个边界函数），分离后61维。
- 在同一个63位质数下，**全部1111个FI约化与分离J按解析因子换算的结果逐项相等**；其余584行也全部为0。
- 原J的10502个非零有理矩阵元、分离J的8396个有理矩阵元，AB结果分别与直接重构/解析构造的anchor逐项一致，并通过原包独立数值holdout。
- 两套直接重构各保留4个未参与该次拟合的数值点，完整矩阵通过检查。
- 新增的求解规模诊断不改变数值算法；在此前的e6 dot<=2输入上核对全部约化和DE条目未变化。

结构和数目以当前目标集合、固定输入与实际日志为证据；不据此宣称任意高dot和任意家族的完备性定理。

## 原包接口的一处已记录问题

最初尝试第二层 `reduce --existing-de`，原包在residue witness读取未初始化的diagonal块时抛出 `KeyError: block0`，尚未进行数值IBP。没有修改原包来绕过；随后按原包默认 `reduce` 自动完成DE阶段再约化，所以将DE和约化次数分别列出。失败日志留在本地 `j_ab/existing_de_failure.log`。

FI求解结束时，部分日志元数据的Mathematica `AssociateTo` 使用了不正确的参数个数。矩阵和裁剪计数已先保存，完整结果通过验证；元数据从这些保存结果恢复，没有重复求解。入口脚本已修正这个记录错误，原日志压缩留存。

## 文件及验证

[comparison.json](comparison.json) 是本次汇总，包含精确计数、维度、质数、实际源文件指纹和验证范围。`reference/` 保存压缩的完整FI结果、J代表点、AB/直接重构结果、实际采样事件和独立检验点。约化后的巨大临时缓存、重复点、可再生bundle和编译二进制不提交git，仍保留在本地。

仅用标准Python重新检查已提交证据：

```sh
python j_integral/tests/e6/fi_vs_j_dot2/verify_saved.py
```

重新计算请使用新的运行目录/checkout，避免覆盖已完成现场。依赖和构建：

```sh
make -C j_integral/src/IBP
make -C j_integral/tests/e6/fi_vs_j_dot2
source j_integral/tests/runtime.sh
python j_integral/tests/e6/fi_vs_j_dot2/prepare_inputs.py
wolframscript -file j_integral/tests/e6/fi_vs_j_dot2/fi/numeric.wl
python j_integral/tests/e6/fi_vs_j_dot2/sample_direct.py
python j_integral/tests/e6/fi_vs_j_dot2/prepare_ab.py
```

AB的两个数学命令（分离版相同，改用 `j_ab_separated` 配置）：

```sh
python -m analytic_bootstrap.upstream learn \
 --config j_integral/tests/e6/fi_vs_j_dot2/j_ab/backend.json \
 --definition j_integral/tests/e6/fi_vs_j_dot2/j_ab/definition.json \
 --setting j_integral/tests/e6/input/setting.json \
 --output j_integral/tests/e6/fi_vs_j_dot2/j_ab/bundle.json
python -m analytic_bootstrap.upstream reduce \
 --config j_integral/tests/e6/fi_vs_j_dot2/j_ab/backend.json \
 --bundle j_integral/tests/e6/fi_vs_j_dot2/j_ab/bundle.json \
 --setting j_integral/tests/e6/input/setting.json --task dot2 \
 --output j_integral/tests/e6/fi_vs_j_dot2/j_ab/result.json
```

`sample_direct.py --separated` 测分离版通用直接重构；`prepare_ab_separated.py` 从一个数值点和分离公式构造该版anchor。`verify_results.py` 从完整本地现场重新汇总和压缩证据；已提交的 `verify_saved.py` 不依赖未提交缓存。
