feat: generalize search reductions to integral tags

Extend search integral labels from bare nu vectors to structured FI/BFI/BBFI
tags with explicit boundary data. Parse FI{...}, BFI[XU]{...}, and
BBFI[XU,YD]{...}, while keeping bare {...} as legacy FI input.

Use the full integral tag for series matching, relation variables, coefficient
expansion, reduction variables, formatting, and master selection. Order
integrals by head first, with FI more complex than BFI and BFI more complex
than BBFI, then by nu complexity; boundary tags only provide deterministic
tie-breaking.

Rename the FI-only second-stage solver and reduction classes to integral
solver/reduction names, and update tool documentation and output sections for
general integral reductions.

Validate that FI-only vac, dp, and dp_planar reductions reproduce the previous
master sets and reductions. Run mixed FI/BFI/BBFI searches for dp and
dp_planar, obtaining 43 and 32 master integrals respectively.