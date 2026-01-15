# PRD-38: Implement Dual Quaternions via Code Generation

**Status**: Draft
**Goal**: Implement dual quaternions as Cl(0,2,1) using the existing code generation infrastructure

## Problem Statement

Dual quaternions are essential for representing rigid body transformations (rotation + translation) in a unified algebraic framework. They're widely used in robotics, computer graphics, and skeletal animation. The library currently lacks this important algebra.

## Mathematical Foundation

### Dual Quaternions as Cl(0,2,1)

In Cl(0,2,1) with basis vectors e1, e2 (both square to -1) and e3 (squares to 0):

| Grade | Blades | Count | Description |
|-------|--------|-------|-------------|
| 0 | 1 | 1 | Scalar |
| 1 | e1, e2, e3 | 3 | Vectors |
| 2 | e12, e13, e23 | 3 | Bivectors |
| 3 | e123 | 1 | Trivector (pseudoscalar) |

Total: 8 basis blades (2³ = 8).

### Basis Properties

- `e1² = -1` (quaternion i direction)
- `e2² = -1` (quaternion j direction)
- `e3² = 0` (dual/nilpotent direction)
- `e12 = e1·e2` (quaternion k direction, squares to -1)
- `e13, e23` (dual quaternion components)
- `e123` (pseudoscalar, squares to 0)

### Structure as Dual Quaternions

A dual quaternion can be written as `q = q_r + ε·q_d` where:
- `q_r` is the "real" quaternion part (grades 0, 1[e1,e2], 2[e12])
- `q_d` is the "dual" quaternion part (grades 1[e3], 2[e13,e23], 3[e123])
- `ε = e3` is the dual unit with `ε² = 0`

### Norm and Conjugate

For dual quaternions, the norm uses **Clifford conjugate**:
```
conjugate negates grades 1, 2, and 3
q * conjugate(q) = |q_r|² + ε·(q_r·conjugate(q_d) + q_d·conjugate(q_r))
```

Since this is a degenerate algebra (r=1), types will implement `DegenerateNormed`.

### Applications

- **Rigid body transformations**: Rotation and translation in one operation
- **Screw motion**: Natural representation of screw axes
- **Skeletal animation**: Efficient bone transformations with linear blending
- **Robotics**: Forward/inverse kinematics

## Proposed Solution

### TOML Specification

Create `algebras/dualquat.toml`:

```toml
# Dual Quaternions
# Cl(0,2,1) - Quaternions extended with dual numbers

[algebra]
name = "dualquat"
module_path = "dualquat"
description = "Dual quaternions Cl(0,2,1)"

[signature]
positive = []
negative = ["e1", "e2"]
zero = ["e3"]

[norm]
# Dual quaternions use Clifford conjugate
# Degenerate norm due to null basis vector
primary_involution = "clifford_conjugate"

[types.Scalar]
grades = [0]
description = "Scalar (grade-0 element)"

[types.Vector]
grades = [1]
description = "Vector with quaternion (i,j) and dual (eps) components"
fields = ["i", "j", "eps"]

[types.Bivector]
grades = [2]
description = "Bivector (k and dual components)"
fields = ["k", "ei", "ej"]

[types.Trivector]
grades = [3]
description = "Trivector (pseudoscalar)"
fields = ["ek"]

[types.Quaternion]
grades = [0, 2]
description = "Pure quaternion part (scalar + k bivector)"
fields = ["w", "k"]

[types.DualQuaternion]
grades = [0, 1, 2, 3]
description = "Full dual quaternion for rigid transformations"
fields = ["w", "i", "j", "eps", "k", "ei", "ej", "ek"]
```

### Type Definitions

The codegen will generate:

1. **Scalar<T>**: Grade-0 element
2. **Vector<T>**: Grade-1 (i, j from quaternion; eps from dual)
3. **Bivector<T>**: Grade-2 (k from quaternion; ei, ej from dual)
4. **Trivector<T>**: Grade-3 (pseudoscalar, ek)
5. **Quaternion<T>**: Even subalgebra of the quaternion part (grades 0, 2[k only])
6. **DualQuaternion<T>**: Full 8-component dual quaternion

### Generated Operations

The codegen will automatically generate:
- Geometric product with correct signs for mixed signature
- Wedge and antiwedge products
- Involutions: reverse, clifford_conjugate, grade_involution
- DegenerateNormed trait (since r > 0)
- Bulk and Unitized wrapper aliases
- Approx traits for floating-point comparison

### Wrapper Types

Since Cl(0,2,1) is degenerate (r=1), the codegen will generate:
- `Bulk<T>` aliases (bulk normalization)
- `Unitized<T>` aliases (weight normalization)

## Implementation Plan

### Step 1: Create TOML Specification

Create `algebras/dualquat.toml` with the specification above.

### Step 2: Add Cl0_2_1 Signature

Add to `src/signature/euclidean.rs`:
```rust
pub struct Cl0_2_1;

impl Signature for Cl0_2_1 {
    type NumBlades = typenum::U8; // 2^3 = 8
    const P: usize = 0;
    const Q: usize = 2;
    const R: usize = 1;

    fn metric(i: usize) -> i8 {
        match i {
            0 => -1, // e1² = -1
            1 => -1, // e2² = -1
            2 => 0,  // e3² = 0
            _ => panic!("invalid basis index")
        }
    }
}

pub type DualQuaternion3 = Cl0_2_1;
```

### Step 3: Add to Build System

Update `build.rs` to include dualquat:

```rust
(
    "dualquat",
    "algebras/dualquat.toml",
    "src/specialized/dualquat/generated",
),
```

### Step 4: Create Module Structure

Create `src/specialized/dualquat/mod.rs`.

### Step 5: Build and Verify

Run `cargo build` to generate the dual quaternion types and verify compilation.

## Files to Modify

1. `algebras/dualquat.toml` (new)
2. `build.rs` - Add dualquat to ALGEBRAS
3. `src/signature/euclidean.rs` - Add Cl0_2_1 signature
4. `src/signature/mod.rs` - Export Cl0_2_1 and DualQuaternion3
5. `src/specialized/mod.rs` - Add dualquat module
6. `src/specialized/dualquat/mod.rs` (new)

## Verification

1. `cargo build` - Codegen and compilation succeed
2. `cargo nextest run` - All tests pass
3. `cargo clippy` - No warnings
4. Verify generated types have correct structure:
   - `DualQuaternion<T>` with 8 fields
   - Geometric product follows dual quaternion multiplication rules
   - DegenerateNormed trait is implemented

## Success Criteria

1. Dual quaternion types are generated via codegen
2. Geometric product correctly handles mixed signature (e1²=e2²=-1, e3²=0)
3. Bulk and Unitized wrapper aliases are generated
4. DegenerateNormed trait is implemented
5. All standard traits are implemented

## Future Enhancements

After initial implementation, consider adding:

1. **Rigid transformation helpers**: `from_rotation_translation()`, `to_rotation_translation()`
2. **Screw motion**: `from_screw_axis()`, `screw_interpolate()`
3. **DLB (Dual quaternion Linear Blending)**: For skeletal animation
4. **Conversion to/from**: Matrices, separate rotation+translation

## Notes

- The degenerate direction (e3² = 0) encodes translation information
- Pure rotations have zero dual part
- Pure translations have zero quaternion part (except scalar)
- Unit dual quaternions satisfy: `q * conjugate(q) = 1`
- The algebra is 8-dimensional, matching the 8 DOF of SE(3) with scale
