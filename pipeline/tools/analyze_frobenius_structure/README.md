# Local Frobenius structure

`analyze_frobenius_structure.wl` wraps DERun's
`FindAsymptoticBehavior` and `CalcTaylor2Num`.

For a differential equation

$$
\frac{dM}{dx}=A(x)M,
$$

it determines the true local power-log structures at a specified point.
If DERun uses a candidate label such as $\mu=-1$ but the first nonzero
coefficient occurs at series order $n=1$, the reported exponent is corrected
to

$$
\mu_{\mathrm{true}}=-1+1=0.
$$

## Usage

```wl
Get["pipeline/tools/analyze_frobenius_structure/analyze_frobenius_structure.wl"];

result = LocalFrobeniusStructure`AnalyzeLocalFrobeniusStructure[
  A,
  delta,
  delta0,
  "Masters" -> masterLabels
];
```

The most useful fields are:

- `result["ActualExponentCounts"]`: number of independent free parameters
  for each true exponent;
- `result["RecurrenceFreeParameters"]`: the candidate exponent and series
  order at which each free variable is created by the coefficient recurrence;
- `result["FreeParameters"]`: true exponent, log power, and first appearance
  of a convenient global coordinate basis after full propagation;
- `result["PerIntegralStructures"]`: actual structures present in every
  integral;
- `result["ParameterProfiles"]`: propagation of each independent parameter
  through all integrals;
- `result["CompleteQ"]`: whether the retained coefficient data span the full
  solution space.

At a simple pole, compare against the residue spectrum with:

```wl
comparison =
  LocalFrobeniusStructure`CompareSimplePoleWithResidue[
    A, delta, delta0, result
  ];
```
