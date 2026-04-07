---
name: scientific-implementation
description: Phase c of the scientific workflow. Implement source code by following a code_structure.md architecture document. Role is a math-literate programmer who can read formulas and translate them to code but does NOT design architecture or derive new formulas. Use when code_structure.md is ready and source code needs to be written, or when existing code needs modification following an updated architecture doc. Triggers on "implement the code", "write the source code", "code this module", "implement the architecture".
---

# Phase c: Code Implementation

Role: **Math-literate programmer** — can read formulas and translate them to code; does not design architecture or derive new formulas.

## Input

- `code_structure.md` (phase b output) — **primary reference**, follow it step by step
- `problem_solution.md` (phase a output) — **lookup reference**, consult only when `code_structure.md` cites a specific section

## Reading Strategy

1. Read `code_structure.md` §3 (architecture overview) to understand the class hierarchy and data flow.
2. For each class to implement, read its specification in `code_structure.md` §4.
3. When the spec says "implements `problem_solution.md` §X.Y", open that section and read **only the cited formula**. Do not read surrounding derivation.
4. Translate the formula to code following the symbol table in `code_structure.md`.

## Implementation Order

Bottom-up: implement classes with no dependencies first, then classes that depend on them. The specific order is determined by the class hierarchy in `code_structure.md` §3.

## Formula Translation Rules

### Match symbols to code names using the symbol table

`code_structure.md` provides a symbol table mapping math symbols to code names. Use it consistently. For example, if the table says $A_{ij}$ maps to `A[i][j]`, use exactly that.

### One formula, one function

If `code_structure.md` says a function implements formula §X.Y, that function should contain **exactly** the computation described by that formula. Do not combine multiple formulas into one function, and do not split one formula across functions, unless `code_structure.md` explicitly says to.

### Respect preconditions

If `code_structure.md` states a precondition, either:
- Add a runtime check (`assert` or `throw`) at the function entry, or
- Document in a comment that the caller guarantees this

### Loop structure mirrors formula structure

When a formula contains a summation $\sum_{i,j}$ over indices, the code should have corresponding nested loops over those indices. Keep the loop structure as close to the formula structure as possible.

### Index offsets: be explicit

Mathematical indices often start from 1; code indices from 0. When the formula uses 1-based indexing, comment the offset:

```cpp
// i in formula is 1-based; idx is 0-based (idx = i - 1)
for (int idx = 0; idx < n; ++idx) { ... }
```

## Coding Conventions

Follow whatever conventions are specified in `code_structure.md` §2. For any non-trivial computation, cite the formula section:

```cpp
/// Implements problem_solution.md §3.2: [name of formula]
void solve(...) { ... }
```

## Compile-Test Cycle

After implementing each class or coherent group:
1. Write a minimal compilation test (instantiate the class, call basic methods)
2. If `code_structure.md` provides expected input/output examples, verify against them
3. Fix compilation errors before moving to the next class

## Anti-Patterns

| ❌ Don't | ✅ Do |
|:---|:---|
| Redesign the architecture while implementing | Follow `code_structure.md` as given; flag concerns to the architect |
| Read all of `problem_solution.md` to "understand the context" | Read only the cited §X.Y when looking up a formula |
| Add optimization not described in the architecture | Implement exactly what's specified; optimize later |
| Skip precondition checks for "obvious" invariants | Add checks or comments for every stated precondition |
| Use different variable names than the symbol table | Match the symbol table exactly |
