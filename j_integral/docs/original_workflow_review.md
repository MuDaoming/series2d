# 原有方法的数学层次与实现核对

## 1. 本文范围与研究目标

本文记录在设计新方法前，对原项目数学说明和核心实现的重新核对。
它不是新的 J 算法设计，也不宣称已经得到通用的 J 展开程序。

用户明确的最终目标只有一个：原始二圈费曼积分的约化。
FBI 表示、高阶展开、辅助参数、cut、BL 和 DE 结构约束都服务于这一目标。

新研究对象已经确定为

$$
X_1=(1-z)t,\qquad X_2=(1-z)(1-t),\qquad X_3=z,
$$

$$
J_\nu(t)=\int_0^1 dz\,(1-z)
 H_\nu((1-z)t,(1-z)(1-t),z),\qquad
I_\nu=\int_0^1dt\,J_\nu(t).
$$

其中 H 包含原始分支被积函数的全部因子；(1-z) 为变量变换的 Jacobian。
不能把 J 与旧程序的 BFI 或局部 FBI 混为一谈。

## 2. 必须区分的数学层次和参数

| 层次 | 对象 | 已知信息或任务 |
|---|---|---|
| FBI | 固定分支参数的积分 I_nu^Delta(X,Y) | 有现成的 IBP、维度迁移和求导公式 |
| MFBI | FBI 层的主积分 | Case 分类与对称性确定；其常点初值作为展开输入 |
| 外层对象 | FI(delta)、BFI(delta)、BBFI(delta) | 先生成展开，随后搜索关系和 masters |
| 原始积分 | FI(delta=1) | 旧路线最终需要返回的对象 |
| 新对象 | 完整积掉 z 后的 J_nu(t) | 尚需设计其直接展开和关系搜索方法 |

Delta 是 FBI 的维度参数。`SeriesSolver` 缓存键中的 `delta` 实际对应这一
维度参数，不是控制外层积分域的辅助 delta。`config.d` 是费曼积分维度；
两圈时 FBI 基准维度 Delta=L*d/2=d。

FBI masters 的数目、搜索矩阵未知多项式系数的数目、外层积分 masters
的数目，是三种不同的计数。

## 3. FBI 家族分析：S、Case 与主积分

数学来源：`docs/background/B001_feynman_integral_and_fbi.md`、
`docs/problem_solution/PS002_fbi_reduction.md`、`expand/docs/problem_solution.md`。

实现：`expand/src/family.tpp`、`expand/src/sector.tpp`。

S 的分支关联块给出传播子的 branch，右下块是二次型矩阵 R。
`Family::findSectors` 枚举覆盖全部 branch 的 sector；不覆盖全部 branch 的
sector 在当前这套展开框架中不作为有效 sector。

`Sector::rowReduce` 同时计算 S 的行化简和行操作矩阵。
非奇异时行操作矩阵就是 S 的逆；奇异时还需寻找零空间。
求得 C_b、z_alpha，并令 C=sum_b C_b。

当前 C++ 使用 Case 0/1/2/3；根目录的部分说明使用 Case 1/2/3/4。
应按条件而不是只按编号对照：

| 条件 | C++ Case | 处理 |
|---|---:|---|
| S 非奇异、C 非零 | 0 | 含一个 corner MFBI；IBP 降 dot，维度迁移对齐 Delta |
| S 非奇异、C 为零 | 1 | 约化到子 sector |
| S 奇异、C 非零 | 2 | 通过降指标约化 |
| S 奇异、C 为零 | 3 | 通过非零 z 分量转移指标，进入更低 sector |

`Sector` 将所需有理系数保存成多项式分子和公分母，供后续 LRR 使用。
这一步不是完成外层费曼积分约化。

## 4. expand：约化与求导结合的直接逐阶展开

数学来源：`expand/docs/problem_solution.md` 第 3—6 节、PS003 第 4 节。
核心实现：`expand/src/series_solver.tpp`。

FBI 的求导公式是

$$
\partial_X I_\nu^\Delta
=-\frac12\sum_{i,j}(\partial_X R_{ij})
 f_{ij}(\nu) I_{\nu+e_i+e_j}^{\Delta+1},
$$

其中非对角项 f_ij=nu_i*nu_j，对角项 f_ii=nu_i*(nu_i+1)。
Y 方向类似。求和限于当前 sector 的活跃传播子。

直接逐阶计算的依赖为：

1. 给定 MFBI 的零阶初值。
2. 主积分的 d+1 阶由求导公式使用非主 FBI 的 d 阶。
3. 非主 FBI 的 d 阶通过已有 FBI 关系递归计算。
4. 所需 FBI key 的已经完成阶数被缓存，不重复求解。

实际调用链：

```text
getFBISeries(nu, Delta, needDeg)
  -> 规范化 nu，查询缓存及已经完成阶数
  -> master key: solveMasterAtDeg
       -> p>0: solveMasterCoeffX
       -> p=0: solveMasterCoeffY
       -> 请求 shifted FBI 的低一阶
  -> nonmaster key: reduceFBIAtDeg
       -> Case 0/1/2/3 的已知关系
       -> 递归请求源 FBI 的所需阶数
       -> solveLRRAtDeg
```

所有关系统一为

$$
D(X,Y)g(X,Y)=\sum_i N_i(X,Y)f_i(X,Y).
$$

常点 D_00 非零时，

$$
g_{pq}=D_{00}^{-1}\left(
\sum_i\sum_{a,b}(N_i)_{ab}(f_i)_{p-a,q-b}
-\sum_{(a,b)\ne(0,0)}D_{ab}g_{p-a,q-b}
\right).
$$

这里遍历的是固定多项式的单项式，不是两个不断增长的级数之间的卷积。
`polySeriesCoeff`、`sumPolySeriesCoeff`、`solveLRRAtDeg` 直接实现这一点。

本路线不需要先插值得到完整的 FBI 约化系数或 MFBI 连接矩阵，
更不需要先得到外层 FI 的约化或 DE。

## 5. 重定义及复杂度

重定义为 I_tilde=U^pow_U I。对于中间维度迁移，`Redefinition`
按相对基准维度的整数偏移进行奇偶对齐，使不同项之间的 U 幂差是整数。
`applyRatioFactors` 通分、消去负整数幂，将关系保持为多项式系数 LRR。

求导时还要包含重定义产生的 dlog(U) 项；`setRedefinition` 预计算
dR/dX*U^L、dR/dY*U^L 及各 master 的 pow_U。

两圈时最终被积函数为

$$
G_\nu=P_\nu\widetilde I_\nu,\qquad P_\nu=\mathcal J W_\nu,
$$

不再需要额外将 U 的非整数幂级数与 FBI 级数相乘。
W 的分支幂由各 branch 内指标之和减一确定；代码见
`IntegrandExpander::buildFIPolynomial`。

令 N 为展开阶数，W_key 为一个 FBI key 的分母和源项分子单项式总数。
固定 family、目标集合及递归依赖后，系数工作的模型是

$$
O\!\left(N^2\sum_{key}W_{key}\right),
$$

另加固定结构的 master 求导项及实现开销。
这是已有算法的复杂度，不只是系数数量给出的希望。

`expand/experiments/expand2d_profile/analysis.md` 还有实际性能证据：
主要成本在 FBI 递推，投影很小；小 key 集合中展开 500 到 1000 阶的耗时
比接近 4。递推多项式重建、map 和规范化造成额外开销，不改变这一主项。

## 6. 从二维系数到 delta 级数

实现：`integrand_expander.tpp`、`delta_projector.tpp`、`fi_pipeline.tpp`。

原有 delta 只进入积分上下限。平移到常点后，被积函数本身不显含 delta。
若 G=sum c_pq X^p Y^q，FI 的单项式贡献为

$$
c_{pq}
\frac{(1-a)^{p+1}-(-a)^{p+1}}{p+1}
\frac{(1-b)^{q+1}-(-b)^{q+1}}{q+1}
\delta^{p+q+2}.
$$

BFI 固定一个移动边界、积分另一个变量，贡献到 p+q+1 阶；
BBFI 在移动角点取值，贡献到 p+q 阶。九种 head/边界对象共享同一 G_nu。

当前统一输出 delta^0 到 delta^deg，保留 FI 前两个零和 BFI 前一个零。
不能把被截断的输出直接当作已经精确求出了 delta=1 的完整积分。

`expand/experiments/derivative_with_delta/derivative_with_delta.md`
通过 Leibniz 法则解释了 FI 导数为何产生四个 BFI；
这是旧变域对象的边界结构，不是新 J 的定义。

## 7. reconstruct：显式重构的另一条路线

来源：`reconstruct/docs/problem_and_workflow.md`、`code_structure.md`。
已核对的核心实现包括 `fbi_reducer`、`fbi_interpolater`、`de_builder`、
`de_interpolater`、`diffeq` 和 `integratedseries`。

1. `FBIReducer` 在有限域数值点使用相同 FBI 关系，得到到 MFBI 的系数向量。
2. `FBIInterpolater` 通过 Firefly 重构 X,Y 有理函数约化系数。
3. `DEBuilder` 用 FBI 求导公式和数值约化构造 MFBI 连接，并处理 U 重定义。
4. `DEInterpolater` 插值得到显式 A_X、A_Y。
5. `DiffSystem` 对显式连接在常点逐阶求二维级数。

这一路线与 expand 可以交叉验证，但不构成 expand 的前置依赖。
旧工具和旧文档的某些输出约定、调试代码与当前统一投影不同，不能混读。

## 8. search：展开之后才进行的关系搜索

来源：`search/docs/problem_solution.md`、`code_structure.md`、
`ideal_absorbing_poly_search.md` 及对应核心实现。

给定多个 FBI 初值投影 b，搜索

$$
\sum_\alpha\sum_{k=0}^m c_{\alpha k}\delta^k
I_\alpha^{(b)}(\delta)=0.
$$

矩阵行是 (b,n)，列是 (alpha,k)，矩阵元为 a_alpha,b,n-k。
这里的未知量是多项式系数；零空间向量才对应积分关系。

`poly_relation_searcher` 的入口负责 train/check 分割。
`RelationSearcher` 核心本身只是组矩阵和消元；
`relation_certifier.hpp` 检查留出行是否扩张训练行空间。

之后 `integral_solver` 对关系的多项式系数在指定 delta 值求值，
再次做积分层消元，得到 masters 和约化。它不对截断积分级数直接求和。

还存在两种 absorb 实现：

- `absorbing_relation_searcher.tpp`：不同 degree 仍会构造矩阵，
  但删去已知 leading relation 可覆盖的变量。
- `incremental_absorbing_relation_searcher.tpp`：按列 oracle 生成移位系数列，
  增量维护独立列及原变量表达式，并加入关系正规化和留出检查。

旧关系的单项式倍数只在全部移位项仍属于 ansatz 时才可用于 absorb。
Absorb 不等于完整的最小关系模计算。

`search_de` 对已经生成的级数直接求导，使用
d_n=(n+1)c_(n+1)，再搜索导数到给定 master 基的关系。
因此外层 DE 是展开和关系搜索的结果；它不是展开生成器的输入。
当前代码还允许导数目标和右侧基来自不同 sector，以生成矩形交叉块。

## 9. cut、对称性与 pipeline

来源：`pipeline/docs/bl_reduction_cut_logic.md`、`plan.md`，
以及 driver、zero probe、cut search/solve、g_series 等入口。

FBI 初值投影用于分离独立贡献。cut 不应与传播子 sector 简单等同：
并非每个 sector 都有独立 MFBI，且对称 orbit 只保留代表。

`expand` 的 cut 模式要求初值恰有一个非零分量。递归中只保留支持集
包含所选 cut（或相应对称 orbit）的源项。normal 模式不执行此截断。

为了使 master 搜索完整，pipeline 还处理：

- 无自身 MFBI cut 的 sector：把对象送到最近的有效低层 cut。
- 某些 dot 移位对象在自身 cut 上为零：不能据此全局删除，需送低层搜索。
- 低层剥离需要高层 masters 在当前 cut 上的展开，不能只准备低层自身目标。

symmetry 的不同材料也须区分：

- C++ 当前实现对 R 做 branch 置换与 Pak 型规范化，再精确核验矩阵关系。
- `constrained_quadratic_form.md` 提出更完整的约束面表示
  (HRH, HRc, c^T R c/2)；不能把该方案当作已在上述 C++ 中实现。
- Wolfram 的 `SectorSymmetry.wl` 对 U,F 联合规范化；pipeline 的生成器
  进一步筛选不跨 branch、适用于现有展开接口的映射。

## 10. bl_sector_reduce：已知 masters 后的分层重构

来源：该模块的三份数学/结构说明，及 `contribution`、`reducer`、
`approximant_basis`、`sector_tree` 的实现。

当前层实际约化的是扣除已知高层贡献后的 C_s(A)，不是直接对完整 A
反复做同一个拟合。目标为

$$
P_0 C_s(A)=\sum_j P_j M_j,\qquad P_0\ne0.
$$

degree bound 为 m 时使用

$$
K=(r+1)(m+1)+K_{safety}+K_{cert}.
$$

这里 safety/cert 阶数都参与 approximant 计算，不是另设的留出区间。
当前算法从单位多项式基出发，逐阶计算 discrepancy、消元并将 pivot 行
乘以 delta。它不是完整的 BL 分治 shifted minimal approximant basis。
这一区别已在 `bl_algorithm_notes.md` 中明确说明。

虽然类名是 SectorTree，当前 `ancestorsOf` 实现按处理顺序寻找所有
支持集包含当前层的高层（含 orbit 判断），不是只沿一个人为父节点链。

## 11. 后续 DE 结构约束与导数模方法

来源：`pipeline/docs/de_constrained_reduction_method.md`，以及各工具的
README 和 Wolfram/C++ 核心实现。

这些方法的输入已经包含主积分 DE；其作用是压缩后续约化系数的搜索空间。

设 M 的解析局部 jet 为 M^(q)=B_q c，约化极点部分为
sum_(r=1)^d c_r M/x^r，则负幂消失要求

$$
\sum_{r=s}^{d}c_r B_{r-s}=0,\qquad s=1,\ldots,d.
$$

不同极点层必须联合求约束，不能只独立地使用留数或每层的一阶条件。

- 一次因子工具构造完整有限阶 Taylor 方程、投影 jet、检查稳定性，
  再求 Toeplitz 约束的核。
- 不可约因子工具使用 F_p[x]/(d)，构造局部坐标 x(t) 并回到基础域计数；
  proper numerator 到局部各极点层的变换也需要计入。
- 无穷远需计入目标积分自己的允许发散度。
- 低层约化还会继承高层系数的分母；只有对角 DE 一般不足以确定继承结构。
- `series_basis_solver` 只求已给定结构基的坐标，并在剩余级数上验证。
- `de_derivative_basis` 检验多项式乘主积分导数的 ansatz；
  `de_derivative_module_basis` 从 B_(n+1)=B_n'+B_n A 出发做多项式模约化。
- Frobenius 和局部规范化工具用于分析已知 DE，不能反过来作为未知 J
  的直接展开已经解决的证据。

## 12. search_like_blade 的设计与实际状态

来源：其问题文档、结构文档、plan、absorb 说明，以及 `slb_tool.cpp`
和 `src/` 中相应接口。

目标设计是 G1|G2 分块、增量搜索、leading absorb、模板学习、精简重拟合、
收集和按依赖求解。输入约束仍是移位级数系数，不是 Blade 的点值数据库。

当前可运行 `slb_tool` 是全局 dense 原型；`src/pipeline.cpp` 的正式接口
尚为未实现，不能把目录中出现接口或文档设计当成实际计算已接通。
`search/experiments/search_like_blade` 另有已实现的逐行增量搜索实验，
与顶层该原型不是同一个程序。

## 13. 验证依据与本次阅读结论

现有验证包括：C/z 与转换检查、零阶 cache 检查、二维展开交叉比较、
FI/BFI/BBFI 投影检查、额外阶数关系检查，以及具体 DB/DP 算例记录。
dp_planar 的历史故障表明，大有理数转换到有限域时不能先缩为机器整数；
当前 converter 通过完整整数字符串进入 FLINT。

本次只阅读并记录上述方法，没有运行新的物理算例，没有修改原模块代码。

对于新 J，已知的是研究对象及最终目标。尚未从旧代码证明：

1. 如何在没有 J 约化或 DE 的前提下，通用地生成其完整 t 展开；
2. 完整 z 积分后的每阶系数能否用固定成本处理；
3. 对一般 family，端点选择、全部展开分支及初值数据如何确定。

旧程序的常点 LRR 依赖非零常数项，不能不经推导就用于退化端点。
单个等质量 vac 的显式公式只可作为例子或验证依据，不能替代通用算法。
后续设计必须以这些明确边界为起点，避免先假设 J 的 DE 已知，
也避免把被积函数的一维展开等同于完整 J 的展开。
