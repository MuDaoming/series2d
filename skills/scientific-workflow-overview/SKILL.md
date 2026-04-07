---
name: scientific-workflow-overview
description: Overview of the three-phase workflow for implementing scientific/professional problem solutions as code. Use when starting a new scientific computation project, planning the overall workflow, or deciding which phase to enter next. Triggers on phrases like "start a new scientific project", "plan the implementation workflow", "what phase should we do next".
---

# Scientific Problem → Code: Workflow Overview

Three-phase workflow for turning a mathematically-defined professional problem into working code.

## Phases

| Phase | Role | Input | Output | Skill |
|:---:|:---|:---|:---|:---|
| a | Domain PhD (can fill in derivations) | Expert draft | `problem_solution.md` | `scientific-problem-definition` |
| b | Architect (math-literate) | Phase a doc | `code_structure.md` | `scientific-code-architecture` |
| c | Programmer (math-literate) | Phase b doc (+ phase a for lookup) | Source code | `scientific-implementation` |

## Flow

```
Expert draft ──→ [Phase a] ──→ problem_solution.md
                                      │
                               Expert reviews math correctness  ← REQUIRED
                                      │
                                      ▼
                              [Phase b] ──→ code_structure.md
                                      │
                                      ▼
                              [Phase c] ──→ Source code
```

## Key Constraints

1. **Phase a output must be self-contained for a non-expert reader**: every formula used must be explicitly written out; no "easily verified" or "it follows that" jumps.
2. **Phase b must cross-reference phase a**: every function must cite the specific section/formula in `problem_solution.md` it implements.
3. **Expert math review after phase a is mandatory**: LLM may introduce errors during normalization.
4. **Implementation-specific derivations stay in the same document** (as later sections of `problem_solution.md`), not in separate files.

## Document Organization

```
project/
├── docs/
│   ├── background.md                   # Domain background (from expert, optional)
│   └── problem_solution.md                # Standalone problem+solution docs (optional)
├── subproject/
│   └── docs/
│       ├── problem_solution.md   # Phase a output
│       └── code_structure.md         # Phase b output
```

## Cross-Reference Convention

Phase b references phase a using:

```markdown
This function implements `problem_solution.md` §3.2 (the recurrence formula for ...).
```

Phase c implementer reads only the cited section — no need to understand the full derivation chain.
