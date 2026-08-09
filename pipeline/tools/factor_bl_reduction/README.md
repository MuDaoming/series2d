# factor_bl_reduction

Factor one explicitly selected successful BL sector reduction over `GF(p)`.

```bash
factor_bl_reduction \
  <reduction_output> \
  <object_label> \
  <sector> \
  <output_path>
```

Example:

```bash
pipeline/tools/factor_bl_reduction/factor_bl_reduction \
  pipeline/runs/dp/bl_reduce_FI31111111_closed_valid31/output \
  'FI{3,1,1,1,1,1,1,1}' \
  '{1,1,1,1,1,1,1,1}' \
  /tmp/FI311_top_factor_structure
```

The object and sector must match exactly one block. The selected block must be
successful, nonzero, non-master, and contain a denominator and reduction terms.

The output contains:

- monic irreducible factors over the prime field from the BL output header;
- the numerator/denominator factor powers for every master coefficient after
  cancelling common factors;
- the direct polynomial-ansatz unknown count;
- the shared-factor structured unknown count;
- the corresponding minimum number of series coefficients and expansion
  degree by equation counting.

`rank_check_required=1` means that the reported minimum is an equation-counting
minimum. A rank check is still required to prove that those coefficients are
independent for the supplied series.

The wrapper uses an explicit `WOLFRAM_KERNEL` when set, then `math -script`,
and finally `wolframscript`. The same wrapper is used locally and on Tianhe.
