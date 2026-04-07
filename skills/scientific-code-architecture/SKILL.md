---
name: scientific-code-architecture
description: Phase b of the scientific workflow. Design code architecture from a problem_solution.md document. Role is a math-literate software architect who can READ formulas but NOT derive new ones. Use when a problem_solution.md is ready and code architecture needs to be designed, or when an existing code_structure.md needs to be updated after the problem document changed. Triggers on "design the code architecture", "write code_structure.md", "update the architecture", "map formulas to code".
---

# Phase b: Code Architecture Design

Role: **Math-literate architect** — can read and understand given formulas; cannot derive new ones. Produces a `code_structure.md` that a junior programmer can follow.

## Input

- `problem_solution.md` (phase a output, expert-reviewed)

## Output

- `code_structure.md` — complete code architecture document

## Document Structure

```
§1  Overview (project goal, core method, one-paragraph summary)
§2  Code conventions (file layout, naming, template design)
§3  Architecture (class hierarchy, data flow, call graph)
§4  Class-by-class specification (the bulk of the document)
§5  Verification strategy (how to test correctness)
§6  References
```

## Architecture Principles

- **Single Responsibility (SOLID)**: each class owns one mathematical object or one algorithm. Do not mix data representation with solving logic.
- **Minimal Interface (SOLID)**: expose only what downstream callers need. Internal helpers stay private.
- **Match the Math**: class hierarchy should mirror the mathematical structure in `problem_solution.md` §2. Do not introduce abstractions without a counterpart in the math.
- **KISS**: prefer straightforward loops over clever metaprogramming. Scientific code is read by people who think in formulas, not design patterns.

## Architecture Process

### Step 1: Identify computational objects → classes/types

Read `problem_solution.md` §2 (mathematical structures). Each major mathematical object typically maps to a class or type. Create a mapping table:

```markdown
| Math object (from §2) | Class/Type name | Notes |
|:---|:---|:---|
| [object A] | `ClassA` | ... |
| [object B] | `ClassB<T>` | ... |
```

### Step 2: Identify algorithms → solver/pipeline classes

Read `problem_solution.md` §3–§K (solution parts + implementation extensions). Each algorithm or major computational step maps to a class or function:

```markdown
| Algorithm (from §X) | Class/Function name | Notes |
|:---|:---|:---|
| [algorithm 1] | `Solver` | Core computation |
| [algorithm 2] | `Transformer` | Pre/post-processing |
| [end-to-end flow] | `runPipeline()` | Top-level entry |
```

### Step 3: Identify critical complexity

Scan the formulas in `problem_solution.md`. Mark functions as **critical** if:
- The formula has multiple cases or conditional branches
- Multiple terms are summed with different index structures
- There are sign conventions or normalization factors that are easy to get wrong
- The formula involves a recurrence where the current value depends on previously computed values

For critical functions: provide full pseudocode in `code_structure.md`.
For non-critical functions: provide function signature + one-line description.

## Class Specification Format

For each key class, write:

1. **Header declaration** — write the actual `.hpp` content including all member variables and member function signatures. This is the authoritative interface that the implementer codes against.
2. **Annotations** — for each member and method, add comments explaining purpose, formula references, and pre/post-conditions.

Template:

````markdown
### X.Y ClassName

**File**: `include/class_name.hpp`, `src/class_name.cpp`

**Purpose**: [one sentence]

```cpp
// include/class_name.hpp
template <typename T>
class ClassName {
public:
    // --- Construction ---
    ClassName(ParamType param);   // Initialize from [what]

    // --- Core interface ---
    /// Implements problem_solution.md §X.Y ([formula name])
    /// Pre: [precondition]
    /// Post: [postcondition]
    ReturnType compute(ArgType arg);

    // --- Accessors ---
    const MemberType& getData() const;

private:
    MemberType data_;   // [what it stores, ref to problem doc if needed]
    int count_;         // [meaning]
};
```

**Critical method pseudocode** (for methods marked critical in Step 3):
- [Step-by-step pseudocode referencing §X.Y]
````

For non-key classes (simple wrappers, IO utilities), a brief signature table suffices — no need for full header declarations.

## Formula-to-Code Mapping Rules

### Every function that implements a formula MUST cite the source

```markdown
/// Implements problem_solution.md §3.2: [name of the formula]
void solve(OutputType& result, const InputType& input, ...);
```

### Map math symbols to code names consistently

Establish a symbol table in the document:

```markdown
| Math symbol | Code name | Type |
|:---|:---|:---|
| $x$ | `x` | `double` |
| $A_{ij}$ | `A[i][j]` | `Matrix` |
| ... | ... | ... |
```

### Interface contracts (preconditions/postconditions)

Write explicitly for every public method:

```markdown
- **Pre**: [what must be true about inputs]
- **Post**: [what is guaranteed about outputs]
```

## Granularity Guide

| Category | Header declaration | Pseudocode | Notes |
|:---|:---|:---|:---|
| Core algorithm class | Full `.hpp` with all members and methods | Yes, for critical methods | The bulk of the document |
| Mathematical data class | Full `.hpp` | No (semantics in comments suffice) | Key storage decisions documented |
| Pipeline / top-level glue | Function signatures only | Execution order + wiring diagram | "construct → configure → solve → output" |
| Infrastructure (IO, conversion) | Brief signature table | No | Minimal detail |

**Rule of thumb**: if the implementer might get it wrong without seeing the exact interface, write the full header. If it's obvious, a brief signature suffices.

## Verification Strategy

The `code_structure.md` should include a §5 describing how to verify correctness. For each core class or algorithm, specify:

- **Expected input/output examples**: at least one concrete small-scale example with known answer that the implementer can use as a sanity check.
- **Sub-module independent verification**: how each module can be tested in isolation before integration.
- **Regression baselines**: reference data or known analytical results that the full pipeline should reproduce. Specify the comparison criterion (exact match, relative error tolerance, etc.).
- **Incremental verification order**: which modules to verify first, and how to confirm each layer before building the next (mirrors the bottom-up implementation order).

## Anti-Patterns

| ❌ Don't | ✅ Do |
|:---|:---|
| Write complete language-specific implementation | Write pseudocode + function signatures |
| Describe a formula without citing `problem_solution.md` | Always cite §X.Y |
| Design complexity beyond what the formulas require | Match architecture to formula structure |
| Omit preconditions for functions that can fail | State preconditions explicitly |
| Assume the reader understands the math | Cite the exact formula; the reader looks it up |
