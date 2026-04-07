---
name: scientific-problem-definition
description: Phase a of the scientific workflow. Normalize an expert's draft into a rigorous problem_solution.md document. Role is a domain PhD who can read, organize, formalize, and fill in mathematical content following the expert's overall logic. Use when the expert provides a rough draft of a problem definition and solution, or when an existing problem_solution.md needs to be extended with implementation-specific derivations. Triggers on "write problem document", "normalize the math doc", "formalize the solution", "extend the problem document".
---

# Phase a: Problem Definition & Solution Document

Role: **Domain PhD** — can read, organize, and formalize mathematical content. Can fill in intermediate derivation steps, fix errors, and write solutions for simple sub-problems, as long as the overall logic and key derivations follow the expert's draft. Can discuss alternative approaches when the expert requests it. Does NOT initiate new research directions or override the expert's key design decisions.

## Input

Expert provides one or more of:
- Rough draft describing the problem and solution (possibly incomplete, symbols inconsistent)
- Background knowledge references
- Verbal explanation of the approach

## Output

A single `problem_solution.md` with the structure below.

## Document Structure

Use numbered sections. Typical outline:

```
§1  Research objective (what is being computed, final output form)
§2  Mathematical structures (definitions of all objects; may reference background docs)
§3  Solution Part I (first major component of the solution)
§4  Solution Part II (second major component, if any)
§...  Solution Part III, IV, ... (as needed)
§K  Implementation-specific extensions (adaptations required by the chosen implementation strategy)
§N-1 Complete computation flow (pseudocode-level summary)
§N  References
```

Notes:
- §2 may reference definitions from `background.md` documents. Either include the referenced definitions inline or provide explicit cross-references (e.g., "see `background.md` §2.3 for the definition of ...").
- The solution sections (§3, §4, ...) should each cover a logically coherent part of the solution. Their titles should be descriptive.
- §K (implementation extensions) covers any new definitions, reformulations, or derivations that arise from the choice of a specific computational strategy. Structure:
  ```
  §K.1 Motivation (why this adaptation is needed)
  §K.2 New definitions
  §K.3 How each previous relation changes
  §K.4 Simplifications for the specific problem instance
  ```

## Normalization Rules

### Definitions
- Every mathematical object used must have an explicit definition with its domain, type, and notation.
- Use a single notation per concept throughout. If the expert's draft uses two different symbols for the same concept, pick one and state the relationship.
- Notation must be rigorous: subscripts, superscripts, and parameters should not be omitted for brevity.

### Formulas
- Every formula that the solution depends on must be **explicitly written out** in full.
- Forbidden phrases: "it is easy to show", "one can verify", "by a similar argument", "it follows that" (without showing what follows).
- If the expert writes "WLOG assume [condition]", either add the general case or state explicitly that this is a restriction and specify the transformation that makes it equivalent to the general case.

### Computational Flow
- The final "complete computation flow" section should be a step-by-step procedure readable as pseudocode.
- Each step references back to the relevant section (e.g., "Step 2: apply the reduction formula from §3.2").

## After Completion

**Mandatory**: ask the expert to review all formulas and logical chains for mathematical correctness. The expert does NOT need to review formatting, style, or section organization — only whether the math is right.

## Anti-Patterns

| ❌ Don't | ✅ Do |
|:---|:---|
| Invent key formulas or change the expert's overall approach | Fill in intermediate steps and solve simple sub-problems within the expert's framework; ask the expert for key formulas |
| Write "by standard results" | Write the result explicitly |
| Use different symbols for the same object in different sections | Normalize symbols globally |
| Omit subscripts/superscripts for brevity | Write full notation with all indices and parameters |
