# debug

This directory stores reusable debug templates and archived debug cases.

Scope:
- Not a strict regression test suite.
- Used for investigation, diagnosis, and recording conclusions.

Workflow (minimal):
1. Use `templates/gen_S_and_Cz.wl` to prepare S inputs and C/z references.
2. Run expand side (cachedeg0 and targetdeg20).
3. Run reconstruct side (series2d executable).
4. Compare outputs and record findings in a case README.

Layout:
- `templates/`: reusable templates (currently placeholders).
- `cases/<case_name>/`: one folder per debug incident.
  - `input/`: fixed case inputs.
  - `output/expand/`: outputs from expand side.
  - `output/reconstruct/`: outputs from reconstruct side.
  - `output/ref_wl/`: references generated from Mathematica.
  - `README.md`: incident summary, commands, conclusion.
