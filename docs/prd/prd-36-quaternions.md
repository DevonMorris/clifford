# PRD-36: Implement Quaternions via Code Generation

**Status**: Complete
**Goal**: Implement quaternions as Cl(0,2,0) using the existing code generation infrastructure

## Problem Statement

Quaternions are fundamental in computer graphics, robotics, and physics for representing 3D rotations. While the library has 3D Euclidean rotors (even subalgebra of Cl(3,0,0)), it lacks a standalone quaternion implementation.

Quaternions can be elegantly represented as the full Clifford algebra Cl(0,2,0) - two basis vectors that both square to -1.

## Mathematical Foundation

### Quaternion as Cl(0,2,0)

In Cl(0,2,0) with basis vectors e1, e2 where e1² = e2² = -1:

| Grade | Blades | Quaternion Component |
|-------|--------|---------------------|
| 0 | 1 | w (real/scalar) |
| 1 | e1, e2 | x (i), y (j) |
| 2 | e12 | z (k) |

The quaternion basis maps as:
- `1 → 1` (scalar)
- `i → e1`
- `j → e2`
- `k → e12 = e1·e2`

### Verification of Quaternion Properties

With e1² = e2² = -1:

```
i² = e1² = -1 ✓
j² = e2² = -1 ✓
k² = (e12)² = e1·e2·e1·e2 = -e1·e1·e2·e2 = -(-1)(-1) = -1 ✓

ij = e1·e2 = e12 = k ✓
ji = e2·e1 = -e12 = -k ✓
jk = e2·(e12) = e2·e1·e2 = -e1·e2·e2 = e1 = i ✓
kj = (e12)·e2 = e1·e2·e2 = -e1 = -i ✓
ki = (e12)·e1 = e1·e2·e1 = -e1·e1·e2 = e2 = j ✓
ik = e1·(e12) = e1·e1·e2 = -e2 = -j ✓

ijk = e1·e2·e12 = e1·e2·e1·e2 = -1 ✓
```

### Norm and Conjugate

For quaternion q = w + xi + yj + zk:

The **Clifford conjugate** negates grades 1 and 2:
```
conjugate(q) = w - xi - yj - zk
```

The norm squared is:
```
q · conjugate(q) = w² + x² + y² + z²
```

This matches the standard quaternion norm.

## Proposed Solution

### TOML Specification

Create `algebras/quaternion.toml`:

```toml
# Quaternions
# Cl(0,2,0) - 4D hypercomplex numbers for rotations

[algebra]
name = "quaternion"
module_path = "quaternion"
description = "Quaternions Cl(0,2,0)"

[signature]
positive = []
negative = ["e1", "e2"]
zero = []

[norm]
# Quaternions use Clifford conjugate: q * conjugate(q) = |q|²
primary_involution = "clifford_conjugate"

# ============================================================
# Types
# ============================================================

[types.Scalar]
grades = [0]
description = "Scalar (grade-0 element)"

[types.Imaginary]
grades = [1]
description = "Pure imaginary quaternion (i and j components)"
fields = ["i", "j"]

[types.Bivector]
grades = [2]
description = "Bivector component (k = ij)"
fields = ["k"]

[types.Quaternion]
grades = [0, 1, 2]
description = "Quaternion: w + xi + yj + zk"
fields = ["w", "i", "j", "k"]
```

### Type Definitions

The codegen will generate:

1. **Scalar<T>**: Grade-0 element (real part)
2. **Imaginary<T>**: Grade-1 element (i, j components)
3. **Bivector<T>**: Grade-2 element (k component)
4. **Quaternion<T>**: Full quaternion (all grades)

### Generated Operations

The codegen will automatically generate:

- **Geometric product**: Full quaternion multiplication
- **Addition/Subtraction**: Component-wise
- **Scalar multiplication**: Scale all components
- **Involutions**: reverse, clifford_conjugate, grade_involution
- **Normed trait**: norm(), norm_squared(), normalize()
- **Approx traits**: AbsDiffEq, RelativeEq, UlpsEq

### Versor Detection

Since Quaternion has grades [0, 1, 2], it spans all grades in the algebra. The codegen should recognize that unit quaternions can be used for rotations via the sandwich product, though the standard quaternion rotation formula uses a different convention.

## Implementation Plan

### Step 1: Create TOML Specification

Create `algebras/quaternion.toml` with the specification above.

### Step 2: Add to Build System

Update `build.rs` to include quaternion:

```rust
const ALGEBRAS: &[(&str, &str, &str)] = &[
    // ... existing algebras ...
    (
        "quaternion",
        "algebras/quaternion.toml",
        "src/specialized/quaternion/generated",
    ),
];
```

### Step 3: Create Module Structure

Create `src/specialized/quaternion/mod.rs`:

```rust
//! Quaternions Cl(0,2,0)
//!
//! Quaternions are 4-dimensional hypercomplex numbers commonly used
//! for representing 3D rotations.

mod generated;

pub use generated::*;
```

### Step 4: Update Specialized Module

Add to `src/specialized/mod.rs`:

```rust
pub mod quaternion;
```

### Step 5: Build and Verify

Run `cargo build` to generate the quaternion types and verify compilation.

## Files to Modify

1. `algebras/quaternion.toml` (new)
2. `build.rs` - Add quaternion to ALGEBRAS
3. `src/specialized/mod.rs` - Add quaternion module
4. `src/specialized/quaternion/mod.rs` (new)

## Verification

1. `cargo build` - Codegen and compilation succeed
2. `cargo nextest run` - All tests pass
3. `cargo clippy` - No warnings
4. Verify generated types have correct structure:
   - `Quaternion<T>` with fields w, i, j, k
   - Geometric product follows quaternion multiplication rules
   - Norm uses Clifford conjugate

## Success Criteria

1. Quaternion types are generated via codegen
2. Quaternion multiplication matches standard formula: (w1 + v1)(w2 + v2) = w1w2 - v1·v2 + w1v2 + w2v1 + v1×v2
3. Norm is computed correctly: |q|² = w² + x² + y² + z²
4. Unit quaternions can be constructed via normalize()
5. All standard traits (Debug, Clone, PartialEq, etc.) are implemented

## Future Enhancements

After initial implementation, consider adding:

1. **Rotation helpers**: `from_axis_angle()`, `to_axis_angle()`
2. **SLERP**: Spherical linear interpolation
3. **Euler angle conversions**: To/from roll-pitch-yaw
4. **Matrix conversions**: To/from rotation matrices
5. **Interop with 3D Euclidean rotors**: Conversion functions

These would be hand-written extensions in the quaternion module, not codegen.

## Notes

- Quaternions in Cl(0,2,0) use ALL grades (0, 1, 2), unlike 3D rotors which are even-graded
- The quaternion i, j, k map to e1, e2, e12 respectively
- This is a complete 4D algebra (2^2 = 4 basis blades)
