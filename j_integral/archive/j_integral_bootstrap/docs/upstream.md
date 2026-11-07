# Analytic Bootstrap 源码与接入记录

## 来源

用户提供的 `/home/mudaoming/Downloads/analytic_bootstrap-master.zip` 已解压到
`j_integral_bootstrap/extern/analytic_bootstrap-master/`，保留原始目录与源码。
该目录不纳入主仓库；我们自己的接口、文档与测试放在其外面。

- README 版本：0.3.0。
- ZIP 注释中的上游提交：`babd5ce45b0a27cc4c8e211fa09abd9ffc74ace9`。
- ZIP SHA256：`0ef2d6a02b3e62775d1f2135b879a58b3a8f4c24dd64e6dadba51d18e7f570f8`。
- 论文：https://arxiv.org/abs/2610.01154 。
- 上游地址：https://github.com/lixfrank/analytic_bootstrap 。本次使用本地 ZIP。

## 已核对的接口

`README.md`、`docs/CONTRACT.md`、`analytic_bootstrap/spec.py`、
`modes.py`、`ansatz.py` 与 `backends/polynomial_module.py` 表明：

1. `run` 是面向费曼积分 family 的完整生产入口，默认数值后端使用 Blade。
2. `BootstrapSpec.from_inputs`、`build_de_ansatz`、`build_reduction_ansatz`
   可以直接接收单变量有理连接、分母包络、分子次数上界，无需先走 Blade family 入口。
3. `RationalAnsatz.evaluate` 和 `assemble` 分别用于数值拟合与生成有理函数结果。
4. 未知 DE 耦合块的构造不要求输入该耦合块。已知低 sector 连接 $B(t)$
   与当前 sector 对角连接 $A(t)$ 后，构造 $C(t)$ 的候选空间：

$$
\frac{d}{dt}
\begin{pmatrix}S\\T\end{pmatrix}
=
\begin{pmatrix}B&0\\C&A\end{pmatrix}
\begin{pmatrix}S\\T\end{pmatrix}.
$$

这里 $S$ 是已知低 sector 主函数列，$T$ 是当前 sector 主函数列；
主函数顺序必须为 sub、top。用少量数值约化数据确定 $C$ 的剩余参数，
解出本层后才能把它加入下一层的已知连接。

## 对当前项目的接入方向

保留 $X_1=(1-z)t$、$X_2=(1-z)(1-t)$、$X_3=z$，
$J_\nu(t)=\int_0^1G_\nu(t,z)\,dz$。
新方法作用在 $t$ 的有理约化系数上。现有有限域一重积分约化器提供
取样数据，包含物理边界函数；不把完整未知的全局 DE 当作起始输入。

先用已有完整连接的小例子核对约束空间、数值拟合与独立点验证，
再测试逐 sector 构造未知耦合块。分母包络、次数上界的获得方式与成本
必须单独记录；不能把读取完整参考结果取得包络的验证实验称为从头求解。

上游区分有限代数恒等式、有限窗口稳定与解析完备性证书。
测试结果应沿用这个区分，窗口稳定本身不能替代约化结果验证。

后续已完成 vacuum、DB/DP cut 及 44 维全局 J 系统的受控测试。
范围、结果与仍未完成的完整 DP 环节见 [测试记录](results.md)。

## 本机接口测试

`python-flint==0.9.0` 安装到本目录的 `.deps/`，未修改系统 Python。
以下命令通过上游 `test_sampling.py` 的 5 个测试，覆盖增量拟合、独立点、
错误样本拒绝、列顺序以及未知 DE 块的共享参数拟合：

```bash
PYTHONPATH="$PWD/j_integral_bootstrap/.deps:$PWD/j_integral_bootstrap/extern/analytic_bootstrap-master" python3 -m unittest discover -s j_integral_bootstrap/extern/analytic_bootstrap-master/tests -p test_sampling.py -v
```

四路线比较 `test_complete_reduction_routes.py` 尚未完成：首次运行约 40 秒
未返回结果，主进程与子进程均几乎不消耗 CPU，已停止，需另查运行依赖或进程等待。
上述 5 项是上游接口测试，不是本项目 J 积分验证。
