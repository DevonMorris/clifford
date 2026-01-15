# PRD-37: Implement Minkowski Plane via Code Generation

**Status**: Draft
**Goal**: Implement the Minkowski plane as Cl(1,1,0) using the existing code generation infrastructure

## Problem Statement

The library lacks a 2D spacetime algebra. The Minkowski plane Cl(1,1,0) is a fundamental algebra for understanding 2D special relativity and hyperbolic geometry. It provides a simple setting to explore indefinite metrics before moving to higher-dimensional spacetime algebras.

## Mathematical Foundation

### Minkowski Plane as Cl(1,1,0)

In Cl(1,1,0) with basis vectors e1, e2 where e1² = +1 and e2² = -1:

| Grade | Blades | Description |
|-------|--------|-------------|
| 0 | 1 | Scalar |
| 1 | e1, e2 | Vectors (spacelike and timelike) |
| 2 | e12 | Bivector (pseudoscalar) |

Total: 4 basis blades (2² = 4).

### Metric Properties

The indefinite metric gives three types of vectors:
- **Spacelike**: v² > 0 (e.g., e1)
- **Timelike**: v² < 0 (e.g., e2)
- **Null/Lightlike**: v² = 0 (e.g., e1 + e2)

### Norm and Conjugate

For the Minkowski plane, the norm uses **Clifford conjugate**:
```
For a = s + xe1 + ye2 + be12:
conjugate(a) = s - xe1 - ye2 - be12  (negates grades 1 and 2)

a * conjugate(a) = s² - x² + y² - b²
```

Note: This norm is **indefinite** - it can be positive, negative, or zero for non-zero elements.

### Relation to Split-Complex Numbers

The even subalgebra (grades 0 and 2) is isomorphic to the hyperbolic numbers:
- Grade 0: scalar
- Grade 2: e12 where (e12)² = e1·e2·e1·e2 = -e1²·e2² = -1·(-1) = +1

So the bivector squares to +1, giving split-complex structure.

## Proposed Solution

### TOML Specification

Create `algebras/minkowski2.toml`:

```toml
# Minkowski Plane
# Cl(1,1,0) - 2D spacetime algebra

[algebra]
name = "minkowski2"
module_path = "minkowski/dim2"
description = "Minkowski plane Cl(1,1,0)"

[signature]
positive = ["e1"]
negative = ["e2"]
zero = []

[norm]
# Minkowski plane uses Clifford conjugate for norm
# Note: This gives an indefinite norm (can be positive, negative, or zero)
primary_involution = "clifford_conjugate"

[types.Scalar]
grades = [0]
description = "Scalar (grade-0 element)"

[types.Vector]
grades = [1]
description = "Vector in Minkowski plane (can be spacelike, timelike, or null)"
fields = ["x", "t"]

[types.Bivector]
grades = [2]
description = "Bivector (pseudoscalar, squares to +1)"
fields = ["xt"]

[types.Eventor]
grades = [0, 2]
description = "Even-graded element (scalar + bivector), isomorphic to hyperbolic numbers"
fields = ["s", "xt"]

[types.Multivector]
grades = [0, 1, 2]
description = "Full multivector in the Minkowski plane"
fields = ["s", "x", "t", "xt"]
```

### Type Definitions

The codegen will generate:

1. **Scalar<T>**: Grade-0 element
2. **Vector<T>**: Grade-1 element with x (spacelike) and t (timelike) components
3. **Bivector<T>**: Grade-2 element (pseudoscalar)
4. **Eventor<T>**: Even subalgebra (grades 0, 2), isomorphic to hyperbolic numbers
5. **Multivector<T>**: Full algebra element (all grades)

### Generated Operations

The codegen will automatically generate:
- Geometric product with correct signs for indefinite metric
- Wedge and antiwedge products
- Involutions: reverse, clifford_conjugate, grade_involution
- IndefiniteNormed trait (since p > 0 and q > 0)
- Approx traits for floating-point comparison

### Wrapper Types

Since Cl(1,1,0) is indefinite (p=1, q=1, r=0), the codegen should generate:
- `Unit<T>` aliases (implements Normed)
- `Proper<T>`, `Spacelike<T>`, `Null<T>` aliases (implements IndefiniteNormed)

## Implementation Plan

### Step 1: Create TOML Specification

Create `algebras/minkowski2.toml` with the specification above.

### Step 2: Add Cl1_1_0 Signature

Add to `src/signature/euclidean.rs`:
```rust
pub struct Cl1_1_0;

impl Signature for Cl1_1_0 {
    type NumBlades = typenum::U4; // 2^2 = 4
    const P: usize = 1;
    const Q: usize = 1;
    const R: usize = 0;

    fn metric(i: usize) -> i8 {
        match i {
            0 => 1,  // e1² = +1
            1 => -1, // e2² = -1
            _ => panic!("invalid basis index")
        }
    }
}

pub type Minkowski2 = Cl1_1_0;
```

### Step 3: Add to Build System

Update `build.rs` to include minkowski2:

```rust
(
    "minkowski2",
    "algebras/minkowski2.toml",
    "src/specialized/minkowski/dim2/generated",
),
```

### Step 4: Create Module Structure

Create `src/specialized/minkowski/mod.rs` and `src/specialized/minkowski/dim2/mod.rs`.

### Step 5: Build and Verify

Run `cargo build` to generate the Minkowski plane types and verify compilation.

## Files to Modify

1. `algebras/minkowski2.toml` (new)
2. `build.rs` - Add minkowski2 to ALGEBRAS
3. `src/signature/euclidean.rs` - Add Cl1_1_0 signature
4. `src/signature/mod.rs` - Export Cl1_1_0 and Minkowski2
5. `src/specialized/mod.rs` - Add minkowski module
6. `src/specialized/minkowski/mod.rs` (new)
7. `src/specialized/minkowski/dim2/mod.rs` (new)

## Verification

1. `cargo build` - Codegen and compilation succeed
2. `cargo nextest run` - All tests pass
3. `cargo clippy` - No warnings
4. Verify generated types have correct structure:
   - `Vector<T>` with fields x, t
   - Geometric product follows indefinite metric rules
   - IndefiniteNormed trait is implemented

## Success Criteria

1. Minkowski plane types are generated via codegen
2. Geometric product correctly handles mixed signature (e1²=+1, e2²=-1)
3. Norm correctly computes indefinite values
4. Proper, Spacelike, Null wrapper aliases are generated
5. All standard traits are implemented

## Future Enhancements

After initial implementation, consider adding:

1. **Lorentz boost helpers**: `from_rapidity()`, `boost()`
2. **Light cone utilities**: `is_spacelike()`, `is_timelike()`, `is_null()`
3. **Hyperbolic angle functions**: `rapidity()`, `proper_time()`
4. **Connection to hyperbolic numbers**: Explicit isomorphism for even subalgebra

## Notes

- The indefinite metric means norm² can be negative - this is expected behavior
- Null vectors (v² = 0) are important for light cone geometry
- The even subalgebra is isomorphic to hyperbolic numbers Cl(1,0,0)
- This is the simplest indefinite Clifford algebra and serves as a stepping stone to Cl(1,3,0) (Minkowski spacetime)
