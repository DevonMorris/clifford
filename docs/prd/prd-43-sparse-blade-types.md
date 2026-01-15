# PRD-43: Sparse Blade Types for Subspace Constraints

**Status**: Draft
**Goal**: Enable code generation for types that use a subset of blades within a grade

## Problem Statement

The current codegen assumes that each type spans **all** blades of its specified grades. However, Conformal Geometric Algebra (and potentially other algebras) has geometric objects that occupy only a **subset** of blades within a grade:

| Grade | Full Type | Sparse Type | Constraint |
|-------|-----------|-------------|------------|
| 3 | Circle (10 blades) | Line (6 blades) | g = 0 (passes through infinity) |
| 4 | Sphere (5 blades) | Plane (4 blades) | u = 0 (passes through infinity) |
| 2 | Dipole (10 blades) | FlatPoint (4 blades) | v = 0, m = 0 |

This is fundamentally different from the algebraic constraints we currently handle (Study condition, Plücker condition), which express relationships between non-zero coefficients. These are **subspace constraints** where certain coefficients are identically zero.

## Current Limitation

Currently, the TOML type specification:
```toml
[types.Circle]
grades = [3]
fields = ["gx", "gy", "gz", "gw", "vx", "vy", "vz", "mx", "my", "mz"]
```

The codegen assumes `len(fields) == sum(binomial(n, k) for k in grades)`. If we specify:
```toml
[types.Line]
grades = [3]
fields = ["vx", "vy", "vz", "mx", "my", "mz"]  # Only 6 fields
```

The codegen has no way to know **which** 6 of the 10 grade-3 blades these fields correspond to.

## Proposed Solution

### Option A: Explicit Blade Mapping (Recommended)

Add an optional `blades` field that explicitly maps fields to blade indices:

```toml
[types.Line]
grades = [3]
fields = ["vx", "vy", "vz", "mx", "my", "mz"]
blades = ["e415", "e425", "e435", "e235", "e315", "e125"]  # Explicit blade names
```

The codegen would:
1. Parse blade names to compute blade indices
2. Generate a sparse struct with only the specified fields
3. Generate conversions that properly map to/from full multivectors

### Option B: Parent Type Reference

Define sparse types as restrictions of full types:

```toml
[types.Line]
parent = "Circle"
zero_fields = ["gx", "gy", "gz", "gw"]
```

The codegen would:
1. Look up the parent type's blade mapping
2. Remove the zero fields
3. Generate a sparse struct

### Option C: Blade Index Specification

Use numeric blade indices directly:

```toml
[types.Line]
grades = [3]
fields = ["vx", "vy", "vz", "mx", "my", "mz"]
blade_indices = [17, 18, 20, 22, 21, 19]  # Internal blade indices
```

## Recommendation

**Option A (Explicit Blade Mapping)** is recommended because:
1. Self-documenting - blade names make the geometry clear
2. No dependency on parent type ordering
3. Works for any sparse subset, not just "zero some fields"
4. Easier to verify against reference materials

## Implementation Details

### TOML Schema Extension

```toml
[types.TypeName]
grades = [3]                    # Required: grades present
fields = ["f1", "f2", ...]      # Required: field names
blades = ["e123", "e124", ...]  # Optional: explicit blade mapping
                                # If omitted, assumes all blades of grades in order
description = "..."             # Optional
```

### Blade Name Parsing

Blade names follow the pattern `e{i}{j}...` where indices are 1-based:
- `e1` = basis vector 1
- `e12` = e1 ∧ e2
- `e415` = e4 ∧ e1 ∧ e5

The parser must:
1. Extract indices from the name
2. Compute the blade's internal index (bitmap representation)
3. Verify the blade belongs to the specified grades

### Conversion Generation

For sparse types, conversions must handle the mapping:

```rust
impl From<Line<T>> for Multivector<T, Conformal3> {
    fn from(line: Line<T>) -> Self {
        let mut mv = Multivector::zero();
        mv.set(Blade::from_indices(&[4, 1, 5]), line.vx());  // e415
        mv.set(Blade::from_indices(&[4, 2, 5]), line.vy());  // e425
        // ... etc
        mv
    }
}

impl TryFrom<Multivector<T, Conformal3>> for Line<T> {
    type Error = &'static str;

    fn try_from(mv: Multivector<T, Conformal3>) -> Result<Self, Self::Error> {
        // Check that non-Line blades are zero
        if mv.get(Blade::from_indices(&[4, 2, 3])).abs() > T::EPSILON {
            return Err("Not a line: gx component non-zero");
        }
        // ... etc
        Ok(Line::new_unchecked(...))
    }
}
```

### Product Generation

Products involving sparse types have two options:

1. **Promote to full type**: Convert sparse to full, compute product, convert back
2. **Direct sparse product**: Generate optimized code that only computes relevant blades

Option 1 is simpler and recommended for initial implementation. Option 2 is an optimization for later.

## CGA Type Definitions

With this feature, conformal3.toml would define:

```toml
# Full grade-3 type
[types.Circle]
grades = [3]
fields = ["gx", "gy", "gz", "gw", "vx", "vy", "vz", "mx", "my", "mz"]

# Sparse grade-3 type (Line is a Circle through infinity)
[types.Line]
grades = [3]
fields = ["vx", "vy", "vz", "mx", "my", "mz"]
blades = ["e415", "e425", "e435", "e235", "e315", "e125"]
description = "Line (circle through infinity, g=0)"

# Full grade-4 type
[types.Sphere]
grades = [4]
fields = ["u", "x", "y", "z", "w"]

# Sparse grade-4 type (Plane is a Sphere through infinity)
[types.Plane]
grades = [4]
fields = ["x", "y", "z", "w"]
blades = ["e4235", "e4315", "e4125", "e3215"]
description = "Plane (sphere through infinity, u=0)"

# Sparse grade-2 type
[types.FlatPoint]
grades = [2]
fields = ["px", "py", "pz", "pw"]
blades = ["e15", "e25", "e35", "e45"]
description = "Flat point (v=0, m=0)"
```

## Blade Ordering in CGA

For reference, the blade ordering in Cl(4,1,0) with basis (e1, e2, e3, e4, e5):

### Grade 2 (10 blades)
| Index | Blade | Circle Field | Line Uses? |
|-------|-------|--------------|------------|
| 3 | e12 | mx | Yes |
| 5 | e13 | my | Yes |
| 6 | e23 | mz | Yes |
| 9 | e14 | ? | No |
| 10 | e24 | ? | No |
| 12 | e34 | ? | No |
| 17 | e15 | vx | Yes |
| 18 | e25 | vy | Yes |
| 20 | e35 | vz | Yes |
| 24 | e45 | ? | No |

(Note: Exact mapping depends on basis convention - must verify against conformalgeometricalgebra.org)

## Implementation Plan

1. **Update spec parser** to handle optional `blades` field
2. **Add blade name parser** to convert "e123" to blade index
3. **Update TypeGenerator** to handle sparse types
4. **Update ConversionsGenerator** for sparse ↔ multivector
5. **Update TraitsGenerator** for products involving sparse types
6. **Add tests** for sparse type generation
7. **Update conformal3.toml** with Line, Plane, FlatPoint types

## Verification

1. `cargo build` succeeds with sparse types
2. Conversion from Line to Multivector sets only the correct blades
3. Conversion from Multivector to Line fails if non-Line blades are non-zero
4. Products of sparse types give correct results

## Future Work

- Optimized product generation for sparse × sparse
- Automatic detection of sparse results (e.g., Plane ∧ Plane might give a Line)
- Subtype relationships (Line is-a Circle conceptually)

## References

- [Conformal GA Wiki - Primitives](https://conformalgeometricalgebra.org/wiki/index.php?title=Primitives)
- [Conformal GA Wiki - Operations](https://conformalgeometricalgebra.org/wiki/index.php?title=Operations)
