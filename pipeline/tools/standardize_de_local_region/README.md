# Local DE standardization

`StandardizeDEAtPoint[A, delta, a]` uses the local coordinate
`eta = delta-a`.  With `a -> Infinity` it uses `eta=1/delta`, including the
Jacobian factor in the transformed connection.

It then calls DERun's rational Fuchs reduction and eigenvalue normalization.
The result is local: different points generally have different gauge matrices.
It does not claim that one rational gauge standardizes every point at once.
