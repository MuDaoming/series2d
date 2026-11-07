# 通用两圈 J 的数值 IBP

在一个给定的 $t=a$、维数 $d$、运动学和质数 $p$ 下，计算 $J$ 积分及其一次 $t$ 导数到指定主基的约化系数。系数属于 $\mathbb F_p$。程序从 FI 家族构造 FBI 关系，再做完整的 $z$ 方向 IBP，保留低 sector 与边界项。

**输入不需要已有的 J 约化、J 的 DE 或 anchor。** 程序不做 AB 学习、$t$ 有理重构、Taylor 展开或最后一重 $t$ 积分。

## 1. 对象和归一化

一个家族有两个圈动量 $\ell_1,\ell_2$、$E$ 个独立外动量 $p_a$、$N$ 条传播子：

$$
q_i=c_{i1}\ell_1+c_{i2}\ell_2+\sum_{a=1}^E e_{ia}p_a,
\qquad D_i=q_i^2-m_i^2+i0.
$$

`loop` 给出 $c_{i1},c_{i2}$，`shift` 给出 $e_{ia}$，`mass2` 是 $m_i^2$，`gram[a][b]` 是 $p_a\cdot p_b$。点积采用输入中的符号约定，不在程序中额外改变运动学符号。

原始 FI 为

$$
I_{\boldsymbol\nu}^{d}
=\int\frac{d^d\ell_1}{i\pi^{d/2}}
      \frac{d^d\ell_2}{i\pi^{d/2}}
      \prod_{i:\nu_i>0}D_i^{-\nu_i}.
$$

本版本要求 $\nu_i\ge0$，不处理不可约分子或负指标。相同圈动量二次结构的传播子组成一个 branch。输入的 `j.branches[i]` 为0、1、2，对应数学上的 $\mathcal B_1,\mathcal B_2,\mathcal B_3$。普通 J 的每个目标都必须保留三个 branch。

令 $v=\sum_i\nu_i$、$\sigma_b=\sum_{i\in\mathcal B_b}\nu_i$，在每个非空 branch 内令 $x_i=X_by_i$，约束为

$$
\sum_{b=1}^3X_b=1,\qquad
\sum_{i\in\mathcal B_b,\,\nu_i>0}y_i=1.
$$

以 $[dy]_{\boldsymbol\nu}$ 表示每条活跃传播子的 $dy_i$ 以及上述三个 branch 内的 Dirac 约束。配方

$$
\sum_i X_{b(i)}y_i(q_i^2-m_i^2)
=\ell^{\mathsf T}A\ell+2\ell\cdot B+C
$$

后定义

$$
U=\det A,\qquad F=B^{\mathsf T}\operatorname{adj}(A)B-UC.
$$

FBI 为

$$
\mathcal I_{\boldsymbol\nu}^{d}(X)
=\frac{(-1)^v\Gamma(v-d)}{\prod_{i:\nu_i>0}\Gamma(\nu_i)}
\int[dy]_{\boldsymbol\nu}
\frac{\prod_{i:\nu_i>0}y_i^{\nu_i-1}}{(F-i0)^{v-d}}.
$$

当前固定使用

$$
X_1=(1-z)t,\qquad X_2=(1-z)(1-t),\qquad X_3=z.
$$

包含 Jacobian、branch 权重和 U 幂的完整对象是

$$
G_{\boldsymbol\nu}(t,z)
=(1-z)\prod_{b=1}^3X_b^{\sigma_b-1}
U^{v-3d/2}\mathcal I_{\boldsymbol\nu}^{d}(X),
$$

$$
J_{\boldsymbol\nu}(t)=\int_0^1G_{\boldsymbol\nu}(t,z)\,dz,
\qquad I_{\boldsymbol\nu}^{d}=\int_0^1J_{\boldsymbol\nu}(t)\,dt.
$$

不同指标的 Gamma 因子及符号没有任意删除。`J{1,2,0,...}` 的指标顺序始终与输入传播子顺序一致。

## 2. 输入：仅有 FI 家族够不够

**FI 家族决定物理对象，但单次数值 IBP 还需要 J 的分支选择、数值点和待约化对象。** 这些都通过同一个 JSON 输入，不写入程序。

| 字段 | 内容 | 是否需要 |
| --- | --- | --- |
| `schema_version` | 当前为1 | 是 |
| `family.loops` | 两个圈动量名称 | 是 |
| `family.external_momenta` | 独立外动量名称；真空为 `[]` | 是 |
| `family.gram` | 外动量 Gram 矩阵；真空为 `[]` | 是 |
| `family.propagators` | 每条传播子的 `loop`、`shift`、`mass2` | 是 |
| `j.coordinates` | 当前为 `"tz"` | 可省略，默认tz |
| `j.branches` | 每条传播子的branch编号0/1/2 | 建议显式给定；标准单位圈动量路由可推断 |
| `j.boundary` | 下述端点定义 | 是 |
| `j.cut_sector` | 可选的0/1传播子cut掩码 | 否；默认全局保留低sector |
| `point.prime` | 奇质数 $p$ | 是 |
| `point.dimension` | 本次 $d$ 的精确值或已赋值表达式 | 是 |
| `point.t` | 本次 $t=a$，应为非奇异常点 | 是 |
| `point.parameters` | 家族中符号参数的精确赋值 | 使用符号参数时需要 |
| `request.targets` | 待约化 J 的非负整数指标数组 | 与自动生成二选一 |
| `request.generate_targets` | `max_dot` 和 `scope=all/top` | 与显式目标二选一 |
| `request.basis` | 固定的 J/B 主对象名称数组 | 否；省略则在当前候选中选基 |
| `request.derivatives` | `"basis"`、名称数组或 `[]` | 默认对全部主基求导；`[]`表示只约化 |
| `request.include_boundary_reductions` | 是否输出边界候选的约化 | 默认true |
| `options.completion` | 补充IBP关系的有限阶数 | 默认16 |
| `options.germ_order` | 端点Laurent-log局部解的截断阶数 | 默认16 |

`targets` 和 `basis` 可以直接给数组，也可以给相对于请求文件的 JSON 数组文件名。`derivatives` 除 `"basis"` 外也接受数组文件名。家族、目标、基底分别改变只需换输入，不改源码。

精确数值用整数或字符串，例如 `"26/7"`。支持参数表达式的加、减、乘、除、整数幂和括号；先按 `point.parameters` 代入，再在有限域中计算。不接受浮点数。分母在模 $p$ 下不能为零。

省略 `j.branches` 时，仅对路由 $\pm\ell_1$、$\pm\ell_2$ 和 $\pm\ell_1\pm\ell_2$ 自动分为0、1、2。第三组必须具有相同的圈动量二次型，混用 $\ell_1+\ell_2$ 和 $\ell_1-\ell_2$ 会被拒绝。其它路由必须明确分组；核仍检查每组二次型相同。

`generate_targets.max_dot=k` 表示每个活跃sector内 $\sum_{i:\nu_i>0}(\nu_i-1)\le k$。`scope=all` 枚举所有覆盖三个branch的sector；`top`只取全部传播子活跃的sector。省略某传播子即指标为0，不是负指标。

### 端点与边界主对象

`j.boundary` 的示例为：

```json
{"z0": "fbi_analytic", "z1": "regulated_zero"}
```

- `z0="regularized"`：保留一般耦合FBI局部解，端点取原函数乘FBI的 $z^0(\log z)^0$ 系数。
- `z0="fbi_analytic"`：额外声明物理corner FBI在 $z=0$ 解析，使规范化 $f_i$ 从 $z^{n_{3,i}-1}$ 起始；还合并精确相同的物理端点。它是需要物理依据的输入条件，不是程序从有限域数值猜出的性质，也不在所有家族中默认启用。
- `z1="regulated_zero"`：使用维数正规化/解析延拓下的端点零条件。本版本只实现这一 $z=1$ 处方。

`Bz0{...}` 是对应规范化corner MFBI的正则化端点值，`Bz1{...}` 在当前处方中为零。B 对象来自 FBI，不是把同指标 J 在 $t=0$ 取值。它们与 J 一起约化；不能为降低主积分数而随意删除。

## 3. 输出

```json
{
  "schema_version": 1,
  "status": "complete",
  "point": {"prime": 1000003, "dimension": "7/13", "t": "17", "parameters": {}},
  "basis": ["J{...}", "Bz0{...}"],
  "basis_status": "candidate_at_this_point",
  "reductions": {"J{...}": {"J{...}": 123, "Bz0{...}": 456}},
  "derivatives": {"J{...}": {"J{...}": 789}},
  "unclosed_derivatives": [],
  "diagnostics": {}
}
```

输出中的系数是整数 $0,\ldots,p-1$。未列出的系数为0，所有行列都使用完整对象名称。例如

$$
J_\nu(a)=\sum_{M\in\mathrm{basis}}c_{\nu M}(a)M(a),
\qquad
\partial_tM(a)=\sum_{N\in\mathrm{basis}}A_{MN}(a)N(a).
$$

这里是积分恒等式系数的有限域专门化，**不计算积分函数本身的数值**。

`status=complete` 只表示本次请求的约化成功且请求的导数闭合；不表示已在所有维数、运动学及任意dot下证明最小主积分集合。未给固定基底时，基底是当前数值点、当前目标和有限completion阶数下的候选，应通过增加目标/阶数及其它常点检查稳定性。给固定基底时，程序检查其独立性及对本次目标的张成性；不够完备便拒绝或返回 `unclosed`。

退出码为0（请求完成）、2（有未闭合导数）、1（输入、奇异点、固定基底或其它计算错误）。错误写到stderr，失败时不把半成品当作成功输出。输出文件先写临时文件，再原子替换。

`diagnostics` 记录FBI秩、基底大小、补充关系数、物理边界关系数、标量恒等式检查数和各阶段耗时。不同数值点之间不隐式重算、不自动取其它点、不读历史结果。

## 4. 编译与运行

需要C++17、FLINT2.x、GMP和nlohmann/json头文件：

```sh
make -C j_integral/src/IBP
j_integral/src/IBP/build/j-ibp \
  --input j_integral/src/IBP/examples/vacuum.json \
  --output /tmp/vacuum_ibp.json
```

用 `FLINT_PREFIX=/path/to/flint` 指定安装目录；必要时通过 `CPPFLAGS` 加入JSON头文件目录。`--output -` 将唯一的结果JSON写到stdout，诊断进度始终在stderr。

同一程序可分别运行 `examples/vacuum.json`、`double_box.json`、`dp.json`。DP示例使用归档中已验证的496维固定主基作为**基底输入**，不输入任何已有DE系数。测试数据不写入算法。

```sh
make -C j_integral/src/IBP test
```

## 5. 算法与代码来源

将规范化corner MFBI写成列向量 $\mathbf f$，先由通用FBI公式得到

$$
G_\nu=r_\nu\mathbf f,\quad
\partial_z\mathbf f=A_z\mathbf f,\quad
\partial_t\mathbf f=A_t\mathbf f.
$$

在给定 $t=a$ 后，系数属于 $\mathbb F_p(z)$。对有理行向量 $g$ 使用

$$
\int_0^1(\partial_zg+gA_z)\mathbf f\,dz=[g\mathbf f]_0^1.
$$

按sector的三角结构逐分量做有理函数Hermite约化，消去一个分量时把非对角项完整传到低sector，保留端点项及共振关系。缓存相同标量积分，避免同一对象重复降极点。导数使用FBI求导公式和一次双数计算；不会把已取值的系数直接当常数求导。

主要成本是FBI系统构造及这些关于z的有理函数约化，取决于FBI秩、极点结构和目标数；不宣称它只因有一个剩余参数就具有与目标数无关的常数成本。

`kernel/` 从归档的已验证通用C++实现迁入，保留算法，重整本地include关系；新的入口负责输入验证、参数代入、目标/固定基底父对象补齐及输出约定：

| 新文件 | 原来源 |
| --- | --- |
| `kernel/polynomial.cpp` | `archive/j_integral/native/series.cpp` |
| `kernel/fbi.hpp` | `archive/j_integral/native/prepare.hpp` |
| `kernel/ibp.hpp` | `archive/j_integral_methods/native/direct_j.cpp` |
| `kernel/z_integrator.hpp` | `archive/j_integral_global/native/triangular.hpp` |
| `kernel/point.cpp` | `archive/j_integral_global/native/global_point.cpp` |

新程序仅编译本目录代码，不include、执行或读取archive中的程序和结果。archived DP参考结果只由测试显式读取用于比较。

## 6. 本次整理的验证

2026-11-07：普通接口测试通过；另用新入口在 $t=23$ 计算DP完整496维数值DE，与归档中已经独立重构并验证的有理连接逐项比较，包含零项的246016个矩阵元全部一致。连同真空、DB及接口检查共6项测试通过。参考校验只发生在测试中，生产程序不读取已有DE。

复现包括DP的检查：

```sh
python3 j_integral/src/IBP/tests/test_interface.py --reference
```

这项测试显式依赖archive保存的参考文件；普通编译和生产计算没有这个依赖。真空测试还检查固定基底下dot2升到dot3、另一个常点及只请求部分导数。新增的接口并未宣称自动证明任意家族的全局最小主基。
