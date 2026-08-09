# DE-structured series reducer

This directory contains the finite-field linear solver used after a cut-specific
DE has supplied a basis of allowed rational coefficient rows.

`series_basis_solver` solves

$$
I(\delta)=\sum_s a_s F_s(\delta)
$$

from a prefix of the supplied series and verifies the result on every remaining
coefficient.  It does not discover denominators or DE structures; those are
prepared by the run-specific Mathematica front end.

Build with:

```bash
make -C pipeline/tools/de_structured_series_reducer
```

Run with:

```bash
series_basis_solver INPUT SOLUTION TRAIN_EQUATIONS
```

The input contains `P`, `ORDER`, one `TARGET` row, and a list of pre-evaluated
`COLUMN` series.  The output records every solved coordinate and the full-series
verification result.
