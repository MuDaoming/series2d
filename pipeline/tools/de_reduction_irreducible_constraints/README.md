# `de_reduction_irreducible_constraints`

该工具计算不可约多项式因子 $d(x)$ 处，微分方程允许的联合约化极点
分子空间。

加载：

```wl
Get["pipeline/tools/de_reduction_irreducible_constraints/de_reduction_irreducible_constraints.wl"];
```

主要入口：

```wl
result =
  DEIrreduciblePoleConstraints`
    ComputeIrreduciblePoleConstraintStructures[
      A,
      x,
      d,
      maxDI,
      "Modulus" -> p
    ];
```

输入：

- `A`：方程 $dM/dx=A(x)M$ 的矩阵；
- `x`：方程变量；
- `d`：$F_p[x]$ 中的首一不可约因子；
- `maxDI`：需要计算到的约化极点阶数；
- `"Modulus"`：有限域素数。

`result["Structures"][dI]` 同时给出：

- `AllowedDimensionK`：在 $K_d=F_p[x]/(d)$ 上的维数；
- `AllowedDimensionFp`：限制标量到 $F_p$ 后的维数；
- `ProperNumeratorConstraintMatrixFp`：标准 Apart proper numerator
  满足的基础域线性约束；
- `AllowedProperNumeratorBasisFp`：允许空间的一组基础域基；
- proper numerator 到局部极点层的三角变换及其可逆性检查。

检验一个实际约化：

```wl
ok =
  DEIrreduciblePoleConstraints`
    ProperNumeratorTupleContainedQ[
      tuple,
      result["Structures"][dI]
    ];
```

`tuple` 的维数必须为 `{dI, Length[A]}`，极点层顺序是

```wl
{numeratorAtPowerDI, ..., numeratorAtPower1}
```

每个元素是在变量 `x` 中、次数小于 `Degree[d,x]` 的 proper numerator
多项式。

对于 `dDE == dI == 1`，可以启用直接留数分支：

```wl
"FuchsianResidueOnly" -> True
```

此时 `"EquationOrderLimit"` 用作正整数共振的有限扫描上限。若
`dDE != 1` 或 `dI != 1`，工具会自动使用完整局部 Taylor jet 递推。

理论、计数口径和测试结果见：

`pipeline/runs/db/de_reduction_pole_constraints/docs/high_degree_irreducible_factor_constraints.md`
