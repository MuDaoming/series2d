# DE derivative module basis

`de_derivative_module_basis.wl` computes the fixed-order module

$$
\mathcal L_r=\sum_{n=0}^r\mathbb F_p[\delta]^{1\times N}B_n,
\qquad M^{(n)}=B_nM.
$$

It clears one common denominator and applies polynomial row operations until
the nonzero rows are in weak Popov form.  The returned `RationalBasis` is a
uniform DE-only basis; it does not use any reduced target integral.

Important output keys are `ModuleRank`, `RowDegrees`, `RationalBasis`, and
`CommonDenominator`.

`CoordinatesInDerivativeModule[result, targetRow, delta]` performs exact
polynomial module division over the same finite field.  It reports membership,
the polynomial coordinates, their degrees, and the actual scalar coefficient
count.
