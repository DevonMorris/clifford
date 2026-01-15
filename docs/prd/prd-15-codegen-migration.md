# PRD-15: Codegen Migration for Euclidean and Projective Algebras

**Status**: Complete
**Goal**: Replace hand-rolled implementations with generated code while preserving the current public API

## Sub-PRDs

This PRD is broken into four sub-PRDs, one per algebra:

| Sub-PRD | Algebra | Types | Status | Notes |
|---------|---------|-------|--------|-------|
| [PRD-15.1](prd-15.1-euclidean3-migration.md) | Euclidean 3D | Vector, Bivector, Trivector, Rotor, Even | Complete | Most complex, validates pattern |
| [PRD-15.2](prd-15.2-euclidean2-migration.md) | Euclidean 2D | Vector, Bivector, Rotor | Complete | Simpler, confirms pattern |
| [PRD-15.3](prd-15.3-projective3-migration.md) | Projective 3D | Point, Line, Plane, Motor, Flector | Complete | Study + Plücker constraints |
| [PRD-15.4](prd-15.4-projective2-migration.md) | Projective 2D | Point, Line, Motor | Complete | Confirms PGA pattern |

**Recommended order**: 15.1 → 15.2 → 15.3 → 15.4

Each sub-PRD should be implemented as a separate PR with full test coverage before proceeding to the next.

### Note on Existing TOML Files

The `algebras/` directory contains TOML files created during PRD-14 as **test fixtures for the codegen tool**. These files may differ from the specifications in these PRDs because:

1. **Module paths**: Test files use `generated::*` paths; PRDs specify actual `specialized/*` paths
2. **Constraint format**: PRDs use simplified notation for clarity; implementation uses actual parser format

During implementation, the existing TOMLs will be **updated or replaced** to match these PRD specifications. The PRDs describe the desired end state, not the current test fixtures.

## Background

The clifford-codegen framework (PRD-14) provides automatic generation of geometric algebra types, products, and traits from TOML specifications. We now want to migrate the existing hand-rolled implementations in `specialized/euclidean` and `specialized/projective` to use generated code. This provides:

1. **Correctness guarantees** - Algebraic formulas are derived symbolically, not hand-coded
2. **Consistency** - All algebras follow the same patterns
3. **Maintainability** - Changes to generation logic automatically propagate
4. **Extensibility** - Adding new algebras becomes declarative

### Scope

| Module | Status | Migration |
|--------|--------|-----------|
| `euclidean::dim2` | Hand-rolled | **This PRD** |
| `euclidean::dim3` | Hand-rolled | **This PRD** |
| `projective::dim2` | Hand-rolled | **This PRD** |
| `projective::dim3` | Hand-rolled | **This PRD** |
| `conformal::dim3` | Hand-rolled | Future PRD |

### Key Constraint: API Preservation

The existing public API must be preserved. Users should not need to change their code after migration. This means:

- Same type names at same module paths
- Same method signatures
- Same trait implementations
- Same conversion patterns

## Current State Analysis

### What Codegen Generates

The codegen framework produces:

1. **Type definitions** - Structs with private fields
2. **Accessors** - `x()`, `y()`, `z()`, etc.
3. **Basic constructors** - `new()`, `zero()`, `unit_*()`
4. **Norm operations** - `norm()`, `norm_squared()`, `normalized()`, `scale()`
5. **Products** - Geometric, outer, inner, scalar, sandwich as free functions
6. **Standard traits** - `Clone`, `Copy`, `Debug`, `PartialEq`, `Default`
7. **Approx traits** - `AbsDiffEq`, `RelativeEq`, `UlpsEq`
8. **Conversions** - `From`/`Into` for Multivector

### What Hand-Rolled Code Has (Beyond Codegen)

#### Euclidean Types

| Type | Additional Methods | Notes |
|------|-------------------|-------|
| `Vector` | `dot()`, `wedge()`, `cross()`, `perp()`, `geometric()` | Cross product is 3D only |
| `Bivector` | `reverse()`, `dual()` | Dual is 3D only |
| `Trivector` | `reverse()` | 3D only |
| `Rotor` | `from_angle()`, `from_angle_plane()`, `from_vectors()`, `rotate()`, `compose()`, `inverse()`, `lerp()`, `slerp()`, `angle()` | Core rotation API |
| `Even` | `to_rotor()` | Alias of Rotor |

#### Projective Types

| Type | Additional Methods | Notes |
|------|-------------------|-------|
| `Point` | `origin()`, `ideal()`, `from_homogeneous()`, `is_ideal()`, `is_finite()`, `to_cartesian()`, `attitude()`, `bulk_norm()`, `weight_norm()`, `geometric_norm()`, `join()`, `distance()`, `midpoint()` | Homogeneous coordinates |
| `Line` | `from_implicit()`, `x_axis()`, `y_axis()`, `z_axis()`, `from_point_and_direction()`, `from_plucker()`, `direction()`, `moment()`, `normal()`, `attitude()`, `distance_from_origin()`, `meet()`, `distance_to_point()`, `project()`, `reflect()`, `angle()`, `unitized()`, `is_parallel()`, `intersects()`, `plucker_inner()`, `satisfies_plucker_condition()`, `plucker_residual()` | Plücker constraints |
| `Plane` | `from_normal_and_distance()`, `xy()`, `xz()`, `yz()`, `normal()`, `distance()`, `attitude()`, `unitized()` | 3D only |
| `Motor` | `identity()`, `from_rotation()`, `from_translation()`, `from_rotation_around()`, `from_axis_angle()`, `from_line()`, `compose()`, `inverse()`, `lerp()`, `rotation_angle()`, `translation()`, `transform_point()`, `transform_line()`, `transform_plane()`, `is_unitized()`, `satisfies_study_condition()`, `study_residual()` | Study condition |
| `Flector` | `from_plane()`, `from_plane_through_origin()`, `reflect_xy()`, `reflect_xz()`, `reflect_yz()`, `point_part()`, `plane_part()`, `is_pure_reflection()`, `transform_point()`, `transform_line()`, `transform_plane()` | 3D only |

#### Additional Modules (Not Generated)

- `arbitrary.rs` - Proptest `Arbitrary` implementations
- `nalgebra.rs` - nalgebra type conversions
- `rerun.rs` - Rerun visualization support
- `conversions.rs` - Multivector conversions with tolerance validation
- `errors.rs` - Domain-specific error types

## Architecture

### Module Structure

```
src/specialized/
├── euclidean/
│   ├── mod.rs                    # Re-exports
│   ├── dim2/
│   │   ├── mod.rs                # Re-exports all public types
│   │   ├── generated/            # NEW: Generated code (do not edit)
│   │   │   ├── mod.rs
│   │   │   ├── types.rs          # Generated type definitions
│   │   │   ├── products.rs       # Generated product functions
│   │   │   ├── traits.rs         # Generated trait impls
│   │   │   └── arbitrary.rs      # Generated proptest support
│   │   ├── extensions.rs         # NEW: Domain-specific methods
│   │   ├── ops.rs                # Operator overloads (hand-written)
│   │   ├── conversions.rs        # Multivector conversions (hand-written)
│   │   ├── nalgebra.rs           # nalgebra support (hand-written)
│   │   └── rerun.rs              # Visualization (hand-written)
│   └── dim3/
│       └── (same structure)
└── projective/
    ├── mod.rs
    ├── errors.rs                 # Shared error types
    ├── dim2/
    │   └── (same structure)
    └── dim3/
        └── (same structure)
```

### Code Organization Strategy

**Generated code** (`generated/`):
- Type struct definitions with private fields
- Public accessors
- Basic constructors (`new()`, `zero()`)
- Norm operations
- Product functions (not operator overloads)
- Basic trait implementations
- Proptest `Arbitrary` implementations

**Extension code** (`extensions.rs`):
- Domain-specific constructors (`from_angle()`, `from_vectors()`, etc.)
- Domain-specific operations (`rotate()`, `transform_point()`, etc.)
- Geometric relationships (`join()`, `meet()`, `distance()`, etc.)
- Interpolation (`lerp()`, `slerp()`)
- Constraint checking methods

**Hand-written modules**:
- `ops.rs` - Operator overloads using generated products
- `conversions.rs` - Multivector conversions with validation
- `nalgebra.rs` - nalgebra conversions
- `rerun.rs` - Visualization

### Type Re-export Pattern

```rust
// src/specialized/euclidean/dim3/mod.rs

// Re-export generated types
pub use generated::types::{Vector, Bivector, Trivector, Rotor, Even};

// Extensions are implemented directly on the generated types
mod extensions;  // impl blocks extending generated types

// Operator overloads
mod ops;

// Other modules
mod conversions;
#[cfg(any(test, feature = "proptest-support"))]
pub mod arbitrary;
#[cfg(any(feature = "nalgebra-0_32", feature = "nalgebra-0_33", feature = "nalgebra-0_34"))]
mod nalgebra;
#[cfg(feature = "rerun")]
mod rerun;

mod generated;
```

## Implementation Phases

### Phase 1: Create Algebra Specifications

Create/update TOML specifications for all target algebras.

#### 1.1 Euclidean 2D (`algebras/euclidean2.toml`)

```toml
[algebra]
name = "euclidean2"
module_path = "euclidean::dim2"
description = "2D Euclidean Geometric Algebra"

[signature]
positive = ["e1", "e2"]

[blades]
e1 = "x"
e2 = "y"
e12 = "xy"

[types.Vector]
grades = [1]
description = "2D vector"
fields = ["x", "y"]

[types.Bivector]
grades = [2]
description = "2D bivector (pseudoscalar)"
fields = ["xy"]

[types.Rotor]
grades = [0, 2]
description = "2D rotation element"
fields = ["s", "xy"]

[types.Even]
grades = [0, 2]
description = "Even subalgebra"
fields = ["s", "xy"]
alias_of = "Rotor"

[products.geometric]
Vector_Vector = "Rotor"
Rotor_Rotor = "Rotor"

[products.outer]
Vector_Vector = "Bivector"

[products.scalar]
Vector_Vector = "T"

[options]
generate_serde = true
generate_arbitrary = true
generate_tests = true
```

#### 1.2 Euclidean 3D (`algebras/euclidean3.toml`)

Already exists - verify it matches current implementation.

#### 1.3 Projective 2D (`algebras/projective2.toml`)

```toml
[algebra]
name = "projective2"
module_path = "projective::dim2"
description = "2D Projective Geometric Algebra (Point-based)"

[signature]
positive = ["e1", "e2"]
zero = ["e0"]

[blades]
e1 = "e1"
e2 = "e2"
e0 = "e0"
e12 = "e12"
e20 = "e20"
e01 = "e01"
e012 = "e012"

[types.Point]
grades = [1]
description = "2D point in homogeneous coordinates"
fields = ["e1", "e2", "e0"]

[types.Line]
grades = [2]
description = "2D line"
fields = ["e12", "e20", "e01"]

[types.Motor]
grades = [0, 2]
description = "2D rigid transformation"
fields = ["s", "e12", "e20", "e01"]

[types.Motor.constraints.study]
expression = "s*e012 = 0"
solve_for = "e012"
enforce = "project"

[products.geometric]
Point_Point = "Motor"
Motor_Motor = "Motor"

[products.regressive]
Point_Point = "Line"

[products.outer]
Line_Line = "Point"

[products.sandwich]
Motor_Point = "Point"
Motor_Line = "Line"

[options]
generate_serde = true
generate_arbitrary = true
generate_tests = true
```

#### 1.4 Projective 3D (`algebras/projective3.toml`)

```toml
[algebra]
name = "projective3"
module_path = "projective::dim3"
description = "3D Projective Geometric Algebra (Point-based)"

[signature]
positive = ["e1", "e2", "e3"]
zero = ["e0"]

[blades]
e1 = "e1"
e2 = "e2"
e3 = "e3"
e0 = "e0"
e23 = "e23"
e31 = "e31"
e12 = "e12"
e01 = "e01"
e02 = "e02"
e03 = "e03"
e123 = "e123"
e032 = "e032"
e013 = "e013"
e021 = "e021"
e0123 = "e0123"

[types.Point]
grades = [1]
description = "3D point in homogeneous coordinates"
fields = ["e1", "e2", "e3", "e0"]

[types.Line]
grades = [2]
description = "3D Plücker line"
fields = ["e01", "e02", "e03", "e23", "e31", "e12"]

[types.Line.constraints.plucker]
expression = "e01*e23 + e02*e31 + e03*e12 = 0"
# No solve_for - this is a geometric constraint, not constructive

[types.Plane]
grades = [3]
description = "3D plane"
fields = ["e032", "e013", "e021", "e123"]

[types.Motor]
grades = [0, 2, 4]
description = "3D rigid transformation"
fields = ["s", "e23", "e31", "e12", "e01", "e02", "e03", "e0123"]

[types.Motor.constraints.study]
expression = "s*e0123 - e23*e01 - e31*e02 - e12*e03 = 0"
solve_for = "e0123"
enforce = "project"

[types.Flector]
grades = [1, 3]
description = "3D reflection/glide"
fields = ["e1", "e2", "e3", "e0", "e032", "e013", "e021", "e123"]

[products.geometric]
Motor_Motor = "Motor"

[products.regressive]
Point_Point = "Line"
Point_Line = "Plane"
Point_Point_Point = "Plane"

[products.outer]
Plane_Plane = "Line"
Plane_Line = "Point"
Line_Line = "Point"

[products.sandwich]
Motor_Point = "Point"
Motor_Line = "Line"
Motor_Plane = "Plane"
Flector_Point = "Point"
Flector_Line = "Line"
Flector_Plane = "Plane"

[options]
generate_serde = true
generate_arbitrary = true
generate_tests = true
```

### Phase 2: Generate Code

For each algebra:

```bash
cargo run --package clifford-codegen -- generate algebras/euclidean2.toml \
    -o src/specialized/euclidean/dim2/generated/ --force

cargo run --package clifford-codegen -- generate algebras/euclidean3.toml \
    -o src/specialized/euclidean/dim3/generated/ --force

cargo run --package clifford-codegen -- generate algebras/projective2.toml \
    -o src/specialized/projective/dim2/generated/ --force

cargo run --package clifford-codegen -- generate algebras/projective3.toml \
    -o src/specialized/projective/dim3/generated/ --force
```

### Phase 3: Create Extension Modules

Create `extensions.rs` files that add domain-specific methods to generated types.

#### 3.1 Euclidean Extensions Pattern

```rust
// src/specialized/euclidean/dim3/extensions.rs

use super::generated::types::{Vector, Bivector, Trivector, Rotor};
use super::generated::products;
use crate::scalar::Float;

// ============================================================================
// Vector extensions
// ============================================================================

impl<T: Float> Vector<T> {
    /// Dot product (inner product returning scalar).
    #[inline]
    pub fn dot(&self, other: &Self) -> T {
        products::scalar_vector_vector(self, other)
    }

    /// Wedge product (outer product).
    #[inline]
    pub fn wedge(&self, other: &Self) -> Bivector<T> {
        products::outer_vector_vector(self, other)
    }

    /// Geometric product.
    #[inline]
    pub fn geometric(&self, other: &Self) -> Rotor<T> {
        products::geometric_vector_vector(self, other)
    }

    /// Cross product (3D only, via Hodge dual of wedge).
    #[inline]
    pub fn cross(&self, other: &Self) -> Self {
        self.wedge(other).dual()
    }
}

// ============================================================================
// Rotor extensions
// ============================================================================

impl<T: Float> Rotor<T> {
    /// Identity rotor (no rotation).
    #[inline]
    pub fn identity() -> Self {
        Self::new(T::one(), T::zero(), T::zero(), T::zero())
    }

    /// Create from rotation angle and plane.
    pub fn from_angle_plane(angle: T, plane: Bivector<T>) -> Self {
        let half = angle / T::TWO;
        let (sin_half, cos_half) = (half.sin(), half.cos());
        let plane_norm = plane.normalized();
        Self::new(
            cos_half,
            -sin_half * plane_norm.xy(),
            -sin_half * plane_norm.xz(),
            -sin_half * plane_norm.yz(),
        )
    }

    /// Create rotation from vector a to vector b.
    pub fn from_vectors(a: Vector<T>, b: Vector<T>) -> Self {
        // 1 + ba gives rotor rotating a to b (when normalized)
        let one = T::one();
        let ba = b.geometric(&a);
        Self::new(
            one + ba.s(),
            ba.xy(),
            ba.xz(),
            ba.yz(),
        ).normalized()
    }

    /// Apply rotation to vector: R v R̃
    #[inline]
    pub fn rotate(&self, v: Vector<T>) -> Vector<T> {
        products::sandwich_rotor_vector(self, &v)
    }

    /// Compose rotations: self then other = other * self
    #[inline]
    pub fn compose(&self, other: &Self) -> Self {
        products::geometric_rotor_rotor(other, self)
    }

    /// Inverse rotation (reverse for unit rotors).
    #[inline]
    pub fn inverse(&self) -> Self {
        let norm_sq = self.norm_squared();
        let rev = self.reverse();
        Self::new(
            rev.s() / norm_sq,
            rev.xy() / norm_sq,
            rev.xz() / norm_sq,
            rev.yz() / norm_sq,
        )
    }

    /// Spherical linear interpolation.
    pub fn slerp(&self, other: &Self, t: T) -> Self {
        // Implementation uses logarithm approach
        // ...
    }
}
```

#### 3.2 Projective Extensions Pattern

```rust
// src/specialized/projective/dim3/extensions.rs

use super::generated::types::{Point, Line, Plane, Motor, Flector};
use super::generated::products;
use crate::scalar::Float;

// ============================================================================
// Point extensions
// ============================================================================

impl<T: Float> Point<T> {
    /// Point at the origin.
    #[inline]
    pub fn origin() -> Self {
        Self::new(T::zero(), T::zero(), T::zero(), T::one())
    }

    /// Ideal point (point at infinity) in the given direction.
    #[inline]
    pub fn ideal(dx: T, dy: T, dz: T) -> Self {
        Self::new(dx, dy, dz, T::zero())
    }

    /// Returns true if this is an ideal point (w ≈ 0).
    #[inline]
    pub fn is_ideal(&self) -> bool {
        self.e0().abs() < T::epsilon()
    }

    /// Cartesian x-coordinate (requires finite point).
    #[inline]
    pub fn x(&self) -> T {
        self.e1() / self.e0()
    }

    /// Cartesian y-coordinate (requires finite point).
    #[inline]
    pub fn y(&self) -> T {
        self.e2() / self.e0()
    }

    /// Cartesian z-coordinate (requires finite point).
    #[inline]
    pub fn z(&self) -> T {
        self.e3() / self.e0()
    }

    /// Join of two points: the line through them.
    pub fn join(&self, other: &Point<T>) -> Line<T> {
        products::regressive_point_point(self, other)
    }

    /// Euclidean distance to another point.
    pub fn distance(&self, other: &Point<T>) -> T {
        // Implementation...
    }
}

// ============================================================================
// Motor extensions
// ============================================================================

impl<T: Float> Motor<T> {
    /// Identity motor (no transformation).
    #[inline]
    pub fn identity() -> Self {
        Self::new_unchecked(
            T::one(),
            T::zero(), T::zero(), T::zero(),
            T::zero(), T::zero(), T::zero(),
            T::zero(),
        )
    }

    /// Pure translation motor.
    pub fn from_translation(dx: T, dy: T, dz: T) -> Self {
        let half = T::one() / T::TWO;
        Self::new_unchecked(
            T::one(),
            T::zero(), T::zero(), T::zero(),
            dx * half, dy * half, dz * half,
            T::zero(),
        )
    }

    /// Pure rotation around an axis through origin.
    pub fn from_axis_angle(axis: &crate::specialized::euclidean::dim3::Vector<T>, angle: T) -> Self {
        let half = angle / T::TWO;
        let (sin_half, cos_half) = (half.sin(), half.cos());
        let axis_norm = axis.normalized();
        // Rotation bivector is dual of axis
        Self::new_unchecked(
            cos_half,
            -sin_half * axis_norm.x(),
            -sin_half * axis_norm.y(),
            -sin_half * axis_norm.z(),
            T::zero(), T::zero(), T::zero(),
            T::zero(),
        )
    }

    /// Transform a point: M P M̃
    #[inline]
    pub fn transform_point(&self, p: &Point<T>) -> Point<T> {
        products::sandwich_motor_point(self, p)
    }

    /// Transform a line: M L M̃
    #[inline]
    pub fn transform_line(&self, l: &Line<T>) -> Line<T> {
        products::sandwich_motor_line(self, l)
    }

    /// Transform a plane: M Π M̃
    #[inline]
    pub fn transform_plane(&self, p: &Plane<T>) -> Plane<T> {
        products::sandwich_motor_plane(self, p)
    }

    /// Compose motors: self then other.
    #[inline]
    pub fn compose(&self, other: &Motor<T>) -> Motor<T> {
        products::geometric_motor_motor(other, self)
    }

    /// Returns true if motor satisfies Study condition.
    pub fn satisfies_study_condition(&self, epsilon: T) -> bool {
        let residual = self.study_residual();
        residual.abs() < epsilon
    }

    /// Study condition residual: s*e0123 - (e23*e01 + e31*e02 + e12*e03)
    pub fn study_residual(&self) -> T {
        self.s() * self.e0123()
            - self.e23() * self.e01()
            - self.e31() * self.e02()
            - self.e12() * self.e03()
    }
}
```

### Phase 4: Update Operator Overloads

Update `ops.rs` to use generated product functions.

```rust
// src/specialized/euclidean/dim3/ops.rs

use super::generated::types::{Vector, Bivector, Rotor};
use super::generated::products;
use crate::scalar::Float;
use std::ops::{Add, Sub, Neg, Mul, BitXor};

// ============================================================================
// Vector operators
// ============================================================================

impl<T: Float> Neg for Vector<T> {
    type Output = Self;
    fn neg(self) -> Self {
        Self::new(-self.x(), -self.y(), -self.z())
    }
}

impl<T: Float> Add for Vector<T> {
    type Output = Self;
    fn add(self, rhs: Self) -> Self {
        Self::new(self.x() + rhs.x(), self.y() + rhs.y(), self.z() + rhs.z())
    }
}

impl<T: Float> Mul<T> for Vector<T> {
    type Output = Self;
    fn mul(self, rhs: T) -> Self {
        Self::new(self.x() * rhs, self.y() * rhs, self.z() * rhs)
    }
}

/// Wedge product via ^ operator
impl<T: Float> BitXor<Vector<T>> for Vector<T> {
    type Output = Bivector<T>;
    fn bitxor(self, rhs: Vector<T>) -> Bivector<T> {
        products::outer_vector_vector(&self, &rhs)
    }
}

// Scalar * Vector (commutative)
impl Mul<Vector<f64>> for f64 {
    type Output = Vector<f64>;
    fn mul(self, rhs: Vector<f64>) -> Vector<f64> {
        rhs * self
    }
}

impl Mul<Vector<f32>> for f32 {
    type Output = Vector<f32>;
    fn mul(self, rhs: Vector<f32>) -> Vector<f32> {
        rhs * self
    }
}

// ============================================================================
// Rotor operators
// ============================================================================

impl<T: Float> Mul<Rotor<T>> for Rotor<T> {
    type Output = Self;
    fn mul(self, rhs: Self) -> Self {
        products::geometric_rotor_rotor(&self, &rhs)
    }
}
```

### Phase 5: Update Conversions

Update `conversions.rs` to work with generated types.

```rust
// src/specialized/euclidean/dim3/conversions.rs

use super::generated::types::{Vector, Bivector, Trivector, Rotor};
use crate::blade::Blade;
use crate::multivector::Multivector;
use crate::scalar::Float;
use crate::signature::Euclidean3;

/// Tolerance for conversion validation.
pub const CONVERSION_TOLERANCE: f64 = 1e-10;

// ============================================================================
// Vector conversions
// ============================================================================

impl<T: Float> From<Vector<T>> for Multivector<T, Euclidean3> {
    fn from(v: Vector<T>) -> Self {
        let mut mv = Multivector::zero();
        mv.set(Blade::basis_vector(0), v.x()); // e1
        mv.set(Blade::basis_vector(1), v.y()); // e2
        mv.set(Blade::basis_vector(2), v.z()); // e3
        mv
    }
}

impl<T: Float> TryFrom<Multivector<T, Euclidean3>> for Vector<T> {
    type Error = ConversionError;

    fn try_from(mv: Multivector<T, Euclidean3>) -> Result<Self, Self::Error> {
        // Verify only grade-1 components are non-zero
        let tolerance = T::from_f64(CONVERSION_TOLERANCE);

        // Check scalar is zero
        if mv.scalar_part().abs() > tolerance {
            return Err(ConversionError::NonZeroGrade(0));
        }
        // Check grade-2 is zero
        // ... etc

        Ok(Self::new(
            mv.get(Blade::basis_vector(0)),
            mv.get(Blade::basis_vector(1)),
            mv.get(Blade::basis_vector(2)),
        ))
    }
}

/// Unchecked conversion for AD compatibility (branch-free).
impl<T: Float> Vector<T> {
    pub fn from_multivector_unchecked(mv: &Multivector<T, Euclidean3>) -> Self {
        Self::new(
            mv.get(Blade::basis_vector(0)),
            mv.get(Blade::basis_vector(1)),
            mv.get(Blade::basis_vector(2)),
        )
    }
}
```

### Phase 6: Migrate Tests

Update tests to verify generated code matches expected behavior.

```rust
// tests/euclidean_codegen_consistency.rs

use clifford::specialized::euclidean::dim3::{Vector, Bivector, Rotor};
use clifford::multivector::Multivector;
use clifford::signature::Euclidean3;
use proptest::prelude::*;

proptest! {
    /// Generated geometric product matches generic Multivector
    #[test]
    fn vector_geometric_matches_multivector(
        ax in -100.0f64..100.0, ay in -100.0f64..100.0, az in -100.0f64..100.0,
        bx in -100.0f64..100.0, by in -100.0f64..100.0, bz in -100.0f64..100.0,
    ) {
        let a = Vector::new(ax, ay, az);
        let b = Vector::new(bx, by, bz);

        // Specialized product
        let spec = a.geometric(&b);

        // Generic product
        let mv_a = Multivector::<f64, Euclidean3>::from(a);
        let mv_b = Multivector::<f64, Euclidean3>::from(b);
        let gen = &mv_a * &mv_b;
        let gen_rotor = Rotor::from_multivector_unchecked(&gen);

        prop_assert!(abs_diff_eq!(spec, gen_rotor, epsilon = 1e-10));
    }

    /// Generated sandwich product matches explicit R * v * R̃
    #[test]
    fn rotor_sandwich_matches_explicit(
        r in any::<UnitRotor<f64>>(),
        v in any::<Vector<f64>>(),
    ) {
        // Using sandwich product function
        let sandwich = r.rotate(v);

        // Explicit: R * v * R̃
        let mv_r = Multivector::<f64, Euclidean3>::from(*r);
        let mv_v = Multivector::<f64, Euclidean3>::from(v);
        let explicit = &(&mv_r * &mv_v) * &mv_r.reverse();
        let explicit_vec = Vector::from_multivector_unchecked(&explicit);

        prop_assert!(abs_diff_eq!(sandwich, explicit_vec, epsilon = 1e-10));
    }
}
```

## Verification Checklist

### Per-Algebra Verification

For each algebra (euclidean2, euclidean3, projective2, projective3):

- [ ] TOML specification created/updated
- [ ] Code generated successfully
- [ ] Extensions module created with all domain methods
- [ ] Operators updated to use generated products
- [ ] Conversions updated and tested
- [ ] All existing tests pass
- [ ] New consistency tests pass (generated vs generic)
- [ ] Benchmarks show no regression
- [ ] Documentation builds without warnings

### Global Verification

- [ ] `cargo fmt` passes
- [ ] `cargo clippy --all-features` passes
- [ ] `cargo doc --all-features --no-deps` passes
- [ ] `cargo test --all-features` passes
- [ ] `cargo deny check` passes
- [ ] All public API preserved (no breaking changes)
- [ ] nalgebra conversions work correctly
- [ ] Proptest arbitrary implementations work
- [ ] Rerun visualization works (if feature enabled)

## Migration Order

Recommended order to minimize risk:

1. **Euclidean 3D** - Most complex, best for finding issues
2. **Euclidean 2D** - Simpler, validates pattern
3. **Projective 3D** - Complex with constraints
4. **Projective 2D** - Validates PGA pattern

Each migration should be a separate PR with full test coverage.

## Risks and Mitigations

| Risk | Impact | Mitigation |
|------|--------|------------|
| Generated code differs from hand-rolled | Breaking changes | Extensive consistency tests against Multivector |
| Performance regression | User impact | Benchmark before/after, profile hot paths |
| Missing methods in extensions | Compile errors | Audit all public methods before migration |
| Constraint handling differs | Incorrect geometry | Property tests for Study/Plücker conditions |
| nalgebra conversion breaks | User impact | Test all nalgebra conversions explicitly |

## Future Work

After this PRD is complete:

1. **PRD-16**: Migrate `conformal::dim3` using same pattern
2. **PRD-17**: Add codegen support for additional products (left contraction, etc.)
3. **PRD-18**: Generate operator overloads automatically
4. **PRD-19**: Generate arbitrary implementations automatically

## Summary

| Phase | Deliverables | Effort |
|-------|--------------|--------|
| 1 | TOML specifications | Low |
| 2 | Generated code | Low (automated) |
| 3 | Extension modules | Medium |
| 4 | Operator overloads | Low |
| 5 | Conversions | Medium |
| 6 | Tests | Medium |

Total: ~4 algebras × ~3 person-days = ~12 person-days estimated.
