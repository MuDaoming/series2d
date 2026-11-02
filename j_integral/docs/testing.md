# 当前验证记录

## 1. 权威计数与验证范围

以下为补齐固定点积分核方向关系后的结果。此前 DB 的 5 个、DP 的 9 个 J/边界对象属于冗余表示，不再作为当前 master 计数。

| 算例 | J/边界主对象 | 原始 FI 基 | 独立 FI 对照 |
|---|---:|---:|---|
| 等质量 2loopvac 完整家族 | 2：1 J + 1 边界 | 2：sunset + tadpole 乘积 | dot≤3 的 20 个 FI，Blade 模 p 20/20 一致 |
| DB topsector cut 111111 | 4：3 J + 1 边界 | 0 | dot≤2 的 28 个 FI，Blade 全部约零，28/28 一致 |
| DB 低 sector cut 101111 | 2：1 J + 1 边界 | 1 | dot≤2 的 21 个 FI，Blade 21/21 一致 |
| DP topsector cut 11111111 | 8：7 J + 1 边界 | 15 | dot≤2 的 45 个 FI，Blade 生成的 cut IBP 独立解 45/45 一致 |

这些数属于不同积分对象和不同 cut，不能互相替代。特别是 J 主函数数目不是原始 FI 主积分数的上界：最后一重积分中的有理权重还需要做 IBP。

## 2. 被验证的对象

定义

$$
X_1=(1-z)t,\quad X_2=(1-z)(1-t),\quad X_3=z,
\qquad J_\nu(t)=\int_0^1G_\nu(t,z)dz,
\qquad I_\nu=\int_0^1J_\nu(t)dt.
$$

G 包含 Jacobian、分支权重、U 的幂和 FBI 的完整指标归一化。正则 z=0 端点 $B^{(z)}_{0,\nu}(t)=G_\nu(t,0)$ 显式保留；奇异 z=1 端点使用输入声明的解析零边界约定。该约定不能扩展成所有 t 端点也自动为零。

原始 FI cut 的最后一步 IBP 另外保存完整 J 基上的边界模余项；只在明确的 pinched-branch cut 商空间内将其视为下 sector。完整、未 cut 家族则需提供非零端点及下 sector 的映射，vac 已独立推导并验证。

## 3. vac 的独立解析检查

令 $h=d/2$、$a=t(1-t)$、$L=a+(1-a)z$、$\kappa=-\Gamma(3-d)$。等质量、质量平方 1 的 corner 为

$$
G_{111}=\kappa(1-z)^{1-h}L^{-h},\quad B=\kappa a^{-h}.
$$

由精确总导数

$$
\partial_z[(1-z)G]=-(2-h)G+(1-a)\partial_aG
$$

得到闭合 DE

$$
J'_{111}=\frac{1-2t}{1-a}[(2-h)J_{111}-B],\qquad
B'=-h\frac{1-2t}{a}B.
$$

dot≤1 的 80 阶和 dot≤2 的 160 阶搜索都给出此基；独立有理证书通过。完整 FI 测试进一步覆盖 dot≤3 的 20 个目标，非零 t 端点按 Beta 积分映射处理。推导见 `vac_fi_ibp.md`。

## 4. 消元核的回归

旧程序在消元矩阵奇异时只保留余项，未利用核方向生成额外 period 关系。修正后，对核向量生成 g，精确检查

$$
\int( g'+gA_z)f_0\,dz=[gf_0]_0^1,
$$

把所得 period/端点坐标关系统一应用到所有系数。

`native/tests.cpp` 有一个专门的无穷远 n=2 核回归：未完成时，已知全导数仍返回冗余 period；加入关系后只剩正确端点值。其他检查包括 q 进分治转换、矩阵方向、非零正则端点、奇异有限极点余项、vac 精确 DE 证书。

当前检查核方向的次数上限显式记入输出，不把固定上限解释为所有家族的完备性证明。

## 5. DP 的 J 与 FI 验证

J 层目标：所有 $\nu_i\ge1$、$\sum_i(\nu_i-1)\le2$，45 个 J 加 90 个端点对象。

- prime=1000003、t0=2/5，200 阶，训练到 120 阶。
- prime=1000033、t0=3/7，220 阶，训练到 140 阶。
- 两组均得 8 个对象，DE 闭合，最高关系次数为 37。
- 两个 prime 分别在 8 个新 t 点，不读取训练级数，直接由 FBI z 积分证书检查每点 82 条约化、8 条 DE、45 个零对象，全部通过。

FI 层：利用 J 约化和 DE、正确的 sector cut，自动选择 15 个 FI masters，45 个目标的完整有理证书全部通过。再仅根据 Blade 的 dot-only master 指标换到同一组基，独立比较全部 45 个目标系数，模 p 45/45 一致。其中 15 个为基的恒等行，30 个为非 master 目标的约化。

Blade 对照来自其已生成的 8100 条 maximal-cut IBP，使用 Blade 相同的 FiniteFlow 稀疏求解器求解；方程和求解中不读取 J 路线的约化系数。不是完整未 cut 的 DP 家族约化。详见 `dp_test.md` 和 `dp_fi_ibp.md`。

## 6. DB 非零低 sector

选择原 DB 的 101111，即保留原传播子 {1,3,4,5,6}；局部五指标与原六指标必须按此顺序映射。该输入只依赖 q=p1+p2，q²=1。

完整 J 基是 $\{J_{11111},B^{(z)}_{0,11111}\}$。在原始 FI 的该 sector cut 上，全部 dot≤2 的 21 个 FI 约到 corner，原六指标为 $I_{101111}$；Blade 采用 IntegralOrdering=2、PreferredMI=Dot，独立得到同一无 rank 基及全部系数。

另两个已检查的零基情况：原 DB topsector 111111 的普通 Blade cut IBP 全部为零，并由 J 路线复现；111101 的 maximal-cut master 分析为空，故未将其作为非零正确性测试。

## 7. 性能口径

本机一次记录：DP 200 阶、全部 45 个 J、全系数证书检查的展开约 20 秒；从已知 J 约化/DE 到 45 个 FI 的最后一重 IBP 约 2.56 秒，指定 Blade 基时约 1.66 秒。最后两个数字不包括 J 展开和 Blade 方程生成。

当前仍处理随阶数增长的多项式，不是 O(M) 标量算法；没有宣称百万阶一秒内完成。

## 8. 一键回归

生产入口为 C++；Python 原型只作历史低阶参考。

```bash
make -C j_integral/native all test
bash j_integral/tests/test_fi_ibp.sh
bash j_integral/tests/test_db_fi_ibp.sh
bash j_integral/tests/test_dp_fi_ibp.sh
```

三个脚本从展开重新跑到 FI 约化，并与实际 Blade 运行留下的独立模 p 基准比较。DP 回归还检查 8 个不参与训练的新点。

固定运动学与有限域下的这些验证，不等于一般运动学参数的有理重构，也不等于已完成未 cut 的 DP 全 sector 约化。
