# 两圈六点 e6：J(t) 数值 IBP 与原包 AB 测试

本包复现一个固定运动学设置下的完整 $J(t)$ 微分方程。已有的完整 DE 用作学习 anchor；随后直接调用原版 AB 的 `learn` 和 `solve-de` 命令，由新写的通用 J-IBP 程序提供实际数值约化。没有自写 AB 求解流程或逐 sector 调度。

原始例子来自 arXiv:2412.21054 的两圈六点家族。家族的具体传播子、质量、Gram矩阵与归一化见 `j_integral/tests/e6/README.md` 和 `input/family_request.json`。

## 1. 包内有什么

```text
README.md
setup.sh / build.sh / run_e6.sh
requirements.txt
vendor/analytic_bootstrap/           原版 AB 0.3.0 的所需源码和原生构建脚本
j_integral/src/IBP/                  通用 J 数值 IBP 的 C++ 源码
j_integral/tests/
  backend_worker.py                 原包JSON协议与J-IBP之间的薄接口
  prepare_run.py                   生成配置、记录输入指纹
  verify_result.py                 有理函数恒等比较、实际IBP次数统计
  export_physical.py                撤销共同规范，导出物理J的DE
  runtime.sh
  e6/
    README.md
    input/                         家族、主基、anchor、对角块与原包definition
    reference/                     已验证物理DE、计时、次数与次数统计
    run.sh                         直接执行原包learn和solve-de
```

不打包 FLINT、GMP、MPFR、python-flint安装、Blade、Mathematica、虚拟环境或本机编译的二进制。AB源码没有改动，原始文件指纹保存在 `UPSTREAM_SOURCE_SHA256.json`。不包含DP的AB测试、历史archive或旧数值取样缓存。

## 2. 外部依赖与编译

已验证环境为Linux、Python3.10、GCC11、原生FLINT2.9、python-flint0.9.0。需要已有：

- Python3.10及以上；使用 `setup.sh` 时还需Python的venv/ensurepip组件（Debian/Ubuntu通常为 `python3-venv`）；
- C/C++编译器（C++17、OpenMP）、make、pkg-config；
- FLINT开发头文件和库、GMP、MPFR；
- nlohmann/json C++头文件；
- Python包 `python-flint==0.9.0`。

**此例运行不需要 Blade 或 Mathematica**，因为学习anchor和对角块已经提供。

已有Python环境及依赖时：

```sh
export FLINT_PREFIX=/your/flint/prefix   # 例如 /usr/local 或 /usr
bash build.sh > build.log 2>&1
bash run_e6.sh first_run
```

编译脚本默认尝试 `/usr/local`，不存在时尝试 `/usr`。在自定义Python环境中可用 `PYTHON=/path/to/python bash build.sh`，运行前激活同一环境。`cc` 必须可调用，因为原包的构建脚本使用这个命令名。

如果希望创建独立Python虚拟环境（只安装Python依赖，不下载或打包系统库）：

```sh
FLINT_PREFIX=/your/flint/prefix bash setup.sh > build.log 2>&1
bash run_e6.sh first_run
```

## 3. 数学对象和输入

固定 $d=26/7$、$p=1000003$，其余运动学全部在家族输入中精确给出。三个branch分别含2、3、2条传播子。普通J定义为

$$
X_1=(1-z)t,\quad X_2=(1-z)(1-t),\quad X_3=z,
$$

$$
J_\nu(t)=\int_0^1 dz\,(1-z)
\prod_{b=1}^{3}X_b^{\sigma_b-1}
U^{|\nu|-3d/2}\mathcal I_\nu^d(X),
\qquad I_\nu^d=\int_0^1J_\nu(t)dt.
$$

$\sigma_b$ 是第b个branch中指标之和，$\mathcal I$ 是保留标准Gamma及符号归一化的FBI。完整定义见 `j_integral/src/IBP/README.md`。端点按输入中的 `fbi_analytic` / `regulated_zero` 处方处理。

主对象共79个：78个普通J和1个边界函数；这是全部已测试三branch sector的耦合系统，不只是topsector cut。完整主基见 `input/basis.json`。

| 输入文件 | 内容 |
| --- | --- |
| `family_request.json` | 完整FI/J定义及本次数值IBP请求模板 |
| `definition.json` | 原包变量、主对象分块及依赖 |
| `setting.json` | 学习/生产的质数；其它参数已写入家族，assignments为空 |
| `anchor.json.gz` | 学习所需的完整有理DE |
| `diagonals.json` | 同一设置下的有理对角块；不向数值拟合提供非对角答案 |
| `normalization.json` | 坐标及共同规范转换 |
| `provenance.json` | 原始输入来源和指纹；其中历史路径仅作来源记录，不是运行依赖 |

有理函数以 $\tau=t-2/5$ 保存。AB询问 $\tau=x$ 时，接口调用物理 $t=x+2/5\pmod p$ 的数值IBP。共同规范为

$$
\widetilde M=[t(t-1)]^3M,\qquad
\widetilde A=A+3\left(\frac1t+\frac1{t-1}\right)\mathbf1.
$$

学习、对角块和数值响应使用同一规范；最后导出时撤销该规范。没有Taylor展开，也不计算积分函数本身的数值。

## 4. 如何确认跑对

新运行保存在 `j_integral/tests/e6/runs/first_run/`，已有目录不会覆盖。取样缓存从空目录开始，不使用参考DE来返回数值响应。

成功后查看：

- `verification.json`：`status="complete"`、`all_rational_entries_equal=true`；
- `learn_timings.json`、`solve_timings.json`：原包自带的时间记录；
- `solution.json`：原包有理DE及训练/检查点记录；
- `numeric/receipts/`：每次真正执行数值IBP的记录；
- `physical_de.json.gz`：撤销共同规范后的物理J微分方程。

本例已有运行记录：学习约5秒，求解约18秒，真实数值IBP **10次**（6个拟合点、加块检查共8点、最终检查再加2点），无失败。79×79共6241项有理恒等一致。耗时依机器而变；学习和求解不包含此前获得anchor的成本。

本包只做这个固定设置、单个质数的重算验证；不等同于一般维数/运动学或多质数有理数系数重构。

## 5. 原程序调用与结果格式

`e6/run.sh` 实际执行：

```sh
python -m analytic_bootstrap.upstream learn ...
python -m analytic_bootstrap.upstream solve-de ...
```

后端接口只做输入输出转换及相同数值点的结果缓存。原包自己选择数值点、构造候选、拟合、按依赖顺序求解并验证。

物理DE文件的 `de[row].rhs[column]` 给出对应矩阵元；`numerator_tau`、`denominator_tau` 按 $\tau$ 升幂排列，系数属于 $\mathbb F_p$，未列出的项为0。最复杂元的分子、分母次数均为50。参考结果位于 `reference/`。

文件完整性可在解压目录运行 `sha256sum -c MANIFEST.sha256` 检查。仅打包源码和输入数据，全部原生程序在合作者机器上重新编译。
