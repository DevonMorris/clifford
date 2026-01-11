# PRD-14: Geometric Algebra Code Generator

**Status**: Draft
**Goal**: Eliminate manual implementation of specialized GA types by generating optimized, strongly-typed code from algebraic specifications

## Sub-PRDs

This PRD is implemented through the following sub-PRDs:

| PRD | Title | Status |
|-----|-------|--------|
| [PRD-14.1](prd-14.1-blade-algebra.md) | Blade Algebra Engine | Complete |
| [PRD-14.2](prd-14.2-specification.md) | Specification Format & Parser | Complete |
| [PRD-14.3](prd-14.3-type-generation.md) | Type & Struct Generation | Complete |
| [PRD-14.4](prd-14.4-product-generation.md) | Product Code Generation | Complete |
| [PRD-14.5](prd-14.5-constrained-types.md) | Constrained Types (Unit, Normalized) | Superseded by 14.11 |
| [PRD-14.6](prd-14.6-traits-tests.md) | Trait Impls & Test Generation | Complete |
| [PRD-14.7](prd-14.7-cli-integration.md) | CLI Tool & Build Integration | Complete |
| [PRD-14.8](prd-14.8-constraint-engine.md) | Geometric Constraint Engine | Complete |
| [PRD-14.9](prd-14.9-entity-discovery.md) | Entity Discovery | Complete |
| [PRD-14.10](prd-14.10-product-inference.md) | Product Output Inference | Complete |
| [PRD-14.11](prd-14.11-constraint-simplification.md) | Constraint Simplification | Complete |
| [PRD-14.12](prd-14.12-symbolica-verification.md) | Symbolic Constraint Verification | Draft |
| [PRD-14.13](prd-14.13-constraint-constructors.md) | Constraint-Enforcing Constructors | Draft |
| [PRD-14.14](prd-14.14-flexible-constraints.md) | Flexible Constraint System | Draft |

## Problem Statement

### Current Pain Points

1. **Tedious Manual Implementation**: Each algebra (Euclidean, PGA, CGA) requires manually implementing the same patterns:
   - Type definitions for each grade combination
   - Geometric, inner, outer, regressive products
   - Sandwich products for transformations
   - Unary operations (reverse, dual, conjugate)
   - Norm and normalization
   - Operator overloads (`Add`, `Sub`, `Mul`, `Neg`, etc.)
   - Trait implementations (`Debug`, `Clone`, `PartialEq`, `AbsDiffEq`, etc.)
   - Conversions to/from `Multivector`
   - Arbitrary implementations for property testing
   - nalgebra interop

2. **Error-Prone**: Manual algebraic derivations are prone to sign errors, index mistakes, and inconsistencies. Even with SymPy derivations, transcribing formulas to Rust introduces errors.

3. **Multivector Limitations**:
   - **No semantic typing**: A `Multivector` doesn't distinguish between a rotor, motor, or random element
   - **Performance penalty**: Carries zeros for unused blades (e.g., a 3D rotor uses 4 of 8 coefficients)
   - **Index opacity**: Users must know bitmask-to-blade mappings

4. **Combinatorial Explosion**: For an n-dimensional algebra:
   - 2^n possible grade combinations (types)
   - (2^n)^2 possible binary products between types
   - Each product has unique coefficient formulas

5. **Constraint Boilerplate**: Types with invariants (unit rotors, normalized lines) require:
   - Wrapper types (`UnitRotor<T>`)
   - Validated constructors
   - `Deref`/`AsRef` implementations
   - Arbitrary implementations with constraint satisfaction
   - Re-implementation of operations that preserve constraints

### The Insight

Once you define:
1. **Signature**: Which basis vectors square to +1, -1, or 0
2. **Type map**: Which grade combinations get named types
3. **Constraints**: Which types have invariants (unit, normalized, null)

...the rest is **entirely formulaic**:
- All products are determined by the signature
- Coefficient formulas follow from the geometric product definition
- Return types are determined by grade selection rules
- Constraint-preserving operations are identifiable algebraically

## Solution: clifford-codegen

A code generator that:
1. Takes an **algebra specification** (signature + type definitions + constraints)
2. Generates **complete, optimized Rust modules**
3. Produces **constrained wrapper types** with proper APIs
4. Produces code that is **verified against `Multivector`** via property tests
5. Allows **semantic layering** on top of generated algebraic operations

### Key Principles

1. **Correctness by Construction**: Generated code matches `Multivector` by design
2. **Zero Runtime Overhead**: All type information erased at compile time
3. **Complete Algebra**: Every valid product is generated
4. **Type-Safe Constraints**: Unit/normalized types enforced at compile time
5. **Customizable Semantics**: Users add meaning via extension traits

## Constrained Types System

### The Problem with Unconstrained Types

Many geometric algebra operations only make sense for elements satisfying constraints:

```rust
// Rotation only works with UNIT rotors
fn rotate(&self, v: Vector<T>) -> Vector<T>  // self must satisfy: ||R|| = 1

// Normalization requires NON-ZERO elements
fn normalized(&self) -> Self  // self must satisfy: ||x|| ≠ 0

// CGA points must be NULL vectors
fn from_euclidean(p: Vec3<T>) -> Point<T>  // result satisfies: p · p = 0
```

Without type-level enforcement, users can:
- Pass unnormalized rotors to `rotate()`, getting wrong results
- Call `normalized()` on zero vectors, getting NaN
- Construct invalid CGA points

### Constraint Categories

| Constraint | Meaning | Examples |
|------------|---------|----------|
| `unit` | `‖x‖ = 1` | Rotor, Motor, UnitVector |
| `nonzero` | `‖x‖ ≠ 0` | NonZeroVector, NonZeroBivector |
| `normalized` | Has canonical form (context-dependent) | Point (w=1), Line (moment² + direction² = 1) |
| `null` | `x · x = 0` | CGA Point, CGA PointPair |
| `flat` | Contains e∞ (infinite radius) | CGA Plane, CGA Line |
| `ideal` | Lives in ideal (null) subspace | PGA Translator (no Euclidean bivector) |

### Generated Type Hierarchy

For a type with a constraint, the generator produces both:

```rust
// Base type: no constraint, all operations available
pub struct Rotor<T: Float> {
    s: T,
    xy: T,
    xz: T,
    yz: T,
}

// Constrained wrapper: enforces ||R|| = 1
#[derive(Clone, Copy, Debug)]
#[repr(transparent)]
pub struct UnitRotor<T: Float>(Rotor<T>);
```

### Constrained Type API Pattern

```rust
impl<T: Float> UnitRotor<T> {
    // ========================================
    // Constructors (all guarantee constraint)
    // ========================================

    /// Creates a unit rotor from angle and plane.
    /// Always produces a valid unit rotor.
    pub fn from_angle_plane(angle: T, plane: Bivector<T>) -> Self {
        let half = angle / T::TWO;
        Self(Rotor::new(
            half.cos(),
            plane.xy() * half.sin(),
            plane.xz() * half.sin(),
            plane.yz() * half.sin(),
        ))
    }

    /// Creates the identity rotation.
    pub fn identity() -> Self {
        Self(Rotor::new(T::one(), T::zero(), T::zero(), T::zero()))
    }

    /// Attempts to normalize an arbitrary rotor.
    /// Returns `None` if the rotor is zero.
    pub fn try_from_rotor(r: Rotor<T>) -> Option<Self> {
        let norm = r.norm();
        if norm < T::epsilon() {
            None
        } else {
            Some(Self(Rotor::new(
                r.s() / norm,
                r.xy() / norm,
                r.xz() / norm,
                r.yz() / norm,
            )))
        }
    }

    /// Creates from components without validation.
    ///
    /// # Safety (Logical)
    /// Caller must ensure the components form a unit rotor.
    /// Use for internal operations known to preserve the constraint.
    pub fn new_unchecked(s: T, xy: T, xz: T, yz: T) -> Self {
        Self(Rotor::new(s, xy, xz, yz))
    }

    // ========================================
    // Constraint-Preserving Operations
    // ========================================

    /// Composes two unit rotors. Result is always unit.
    pub fn compose(self, other: Self) -> Self {
        // UnitRotor * UnitRotor = UnitRotor (algebraic fact)
        Self::new_unchecked(/* generated formula */)
    }

    /// Returns the inverse (which equals the reverse for unit rotors).
    pub fn inverse(self) -> Self {
        // For unit rotor: R⁻¹ = R̃
        Self(self.0.reverse())
    }

    /// Interpolates between rotors. Result is always unit.
    pub fn slerp(self, other: Self, t: T) -> Self {
        // Slerp preserves unit constraint
        Self::new_unchecked(/* generated formula */)
    }

    // ========================================
    // Operations returning base type
    // ========================================

    /// Scales the rotor, breaking the unit constraint.
    pub fn scale(self, s: T) -> Rotor<T> {
        self.0.scale(s)  // Returns Rotor, not UnitRotor
    }

    // ========================================
    // Access to inner value
    // ========================================

    /// Returns the underlying rotor.
    pub fn into_inner(self) -> Rotor<T> {
        self.0
    }

    /// Returns a reference to the underlying rotor.
    pub fn as_rotor(&self) -> &Rotor<T> {
        &self.0
    }
}

// Deref for convenient access to Rotor methods
impl<T: Float> Deref for UnitRotor<T> {
    type Target = Rotor<T>;
    fn deref(&self) -> &Self::Target {
        &self.0
    }
}

// From for ergonomic conversion
impl<T: Float> From<UnitRotor<T>> for Rotor<T> {
    fn from(u: UnitRotor<T>) -> Self {
        u.0
    }
}
```

### Constraint Preservation Rules

The generator knows which operations preserve constraints:

| Operation | Preserves Unit? | Preserves NonZero? |
|-----------|-----------------|---------------------|
| `a * b` (both unit) | ✓ | ✓ |
| `a * b` (one unit) | ✗ | ✓ |
| `a + b` | ✗ | ✗ |
| `a - b` | ✗ | ✗ |
| `reverse(a)` | ✓ | ✓ |
| `dual(a)` | ✗ | ✓ |
| `sandwich(unit, x)` | depends on x | ✓ |
| `slerp(unit, unit, t)` | ✓ | ✓ |
| `inverse(unit)` | ✓ | ✓ |

### Arbitrary Implementations for Constrained Types

```rust
// generated/euclidean3/arbitrary.rs

impl<T: Float + Debug> Arbitrary for UnitRotor<T>
where
    Rotor<T>: Arbitrary,
{
    type Parameters = ();
    type Strategy = BoxedStrategy<Self>;

    fn arbitrary_with(_: Self::Parameters) -> Self::Strategy {
        // Generate angle and unit bivector, construct via from_angle_plane
        (-std::f64::consts::PI..std::f64::consts::PI)
            .prop_flat_map(|angle| {
                any::<UnitBivector<T>>().prop_map(move |plane| {
                    UnitRotor::from_angle_plane(T::from_f64(angle), plane.into_inner())
                })
            })
            .boxed()
    }
}

impl<T: Float + Debug> Arbitrary for NonZeroVector<T>
where
    Vector<T>: Arbitrary,
{
    type Parameters = ();
    type Strategy = BoxedStrategy<Self>;

    fn arbitrary_with(_: Self::Parameters) -> Self::Strategy {
        any::<Vector<T>>()
            .prop_filter_map("zero vector", |v| {
                if v.norm_squared() > T::epsilon() {
                    Some(NonZeroVector(v))
                } else {
                    None
                }
            })
            .boxed()
    }
}
```

### Constraint Specification in TOML

```toml
[types.Rotor]
grades = [0, 2]
description = "Even subalgebra element"

[types.Rotor.constraints.unit]
# Defines UnitRotor<T> wrapper
norm = "euclidean"  # ||R||² = s² + b·b
constructors = [
    { name = "from_angle_plane", params = ["angle: T", "plane: Bivector<T>"] },
    { name = "from_vectors", params = ["a: Vector<T>", "b: Vector<T>"] },
    { name = "identity", params = [] },
]
# Operations that preserve unit constraint when inputs are unit
preserving_ops = ["compose", "inverse", "reverse", "slerp"]

[types.Vector]
grades = [1]
description = "Grade-1 vector"

[types.Vector.constraints.unit]
norm = "euclidean"
constructors = [
    { name = "unit_x", params = [] },
    { name = "unit_y", params = [] },
    { name = "unit_z", params = [] },
]
preserving_ops = []  # No binary ops preserve unit for vectors

[types.Vector.constraints.nonzero]
# Defines NonZeroVector<T> wrapper
condition = "norm_squared > epsilon"
constructors = []  # Only try_from
preserving_ops = ["scale_nonzero"]  # scale by nonzero scalar
```

### CGA Null Constraint Example

```toml
[types.Point]
grades = [1]
description = "Conformal point (null vector representation)"

[types.Point.constraints.null]
# x · x = 0 (null vector)
condition = "inner_self == 0"
constructors = [
    { name = "from_euclidean", params = ["x: T", "y: T", "z: T"] },
    { name = "origin", params = [] },
]
# Conformal transformations preserve null property
preserving_ops = ["sandwich_by_rotor", "sandwich_by_translator", "sandwich_by_motor"]
```

## Specification Format

### Complete Algebra Definition (TOML)

```toml
[algebra]
name = "euclidean3"
module_path = "euclidean::dim3"
description = "3D Euclidean Geometric Algebra"

[signature]
positive = ["e1", "e2", "e3"]  # Square to +1
negative = []                   # Square to -1
zero = []                       # Square to 0

# Blade naming conventions
[blades]
e1 = "x"
e2 = "y"
e3 = "z"
e12 = "xy"
e13 = "xz"
e23 = "yz"
e123 = "xyz"

# Type definitions
[types.Scalar]
grades = [0]
description = "Grade-0 scalar"

[types.Vector]
grades = [1]
description = "Grade-1 vector: e₁, e₂, e₃"
fields = ["x", "y", "z"]

[types.Vector.constraints.unit]
norm = "euclidean"
constructors = ["unit_x", "unit_y", "unit_z"]

[types.Vector.constraints.nonzero]
condition = "norm_squared > epsilon"

[types.Bivector]
grades = [2]
description = "Grade-2 bivector: e₁₂, e₁₃, e₂₃"
fields = ["xy", "xz", "yz"]

[types.Bivector.constraints.unit]
norm = "euclidean"
constructors = ["unit_xy", "unit_xz", "unit_yz"]

[types.Trivector]
grades = [3]
description = "Grade-3 pseudoscalar: e₁₂₃"
fields = ["xyz"]

[types.Even]
grades = [0, 2]
description = "Even subalgebra (scalar + bivector)"
fields = ["s", "xy", "xz", "yz"]

[types.Rotor]
grades = [0, 2]
description = "Rotation element"
fields = ["s", "xy", "xz", "yz"]
alias_of = "Even"  # Same storage, different semantic

[types.Rotor.constraints.unit]
norm = "euclidean"
constructors = [
    "identity",
    "from_angle_plane(angle: T, plane: Bivector<T>)",
    "from_angle_axis(angle: T, axis: Vector<T>)",
    "from_vectors(a: Vector<T>, b: Vector<T>)",
]
preserving_ops = ["compose", "inverse", "slerp"]

[types.Odd]
grades = [1, 3]
description = "Odd subalgebra (vector + trivector)"
fields = ["x", "y", "z", "xyz"]

[types.Full]
grades = [0, 1, 2, 3]
description = "Full multivector (all grades)"
fields = ["s", "x", "y", "z", "xy", "xz", "yz", "xyz"]

# Semantic operations (hand-written extensions)
[semantic]
# These methods will be generated as stubs for manual implementation
rotor_methods = ["rotate(v: Vector<T>) -> Vector<T>"]
```

### PGA3D Example

```toml
[algebra]
name = "pga3"
module_path = "projective::dim3"
description = "3D Projective Geometric Algebra"

[signature]
positive = ["e1", "e2", "e3"]
negative = []
zero = ["e0"]

[blades]
# 1 component
e0 = "e0"
e1 = "e1"
e2 = "e2"
e3 = "e3"
# 2 components
e01 = "e01"
e02 = "e02"
e03 = "e03"
e12 = "e12"
e31 = "e31"
e23 = "e23"
# 3 components
e021 = "e021"
e013 = "e013"
e032 = "e032"
e123 = "e123"
# 4 component
e0123 = "e0123"

[types.Plane]
grades = [1]
description = "Plane: d·e₀ + a·e₁ + b·e₂ + c·e₃"
fields = ["d", "a", "b", "c"]

[types.Plane.constraints.normalized]
# Normal vector has unit length: a² + b² + c² = 1
condition = "a*a + b*b + c*c == 1"
constructors = [
    "from_normal_distance(normal: euclidean::dim3::Vector<T>, d: T)",
    "from_point_normal(point: Point<T>, normal: euclidean::dim3::Vector<T>)",
]

[types.Line]
grades = [2]
description = "Plücker line (6 components)"
fields = ["e01", "e02", "e03", "e12", "e31", "e23"]

[types.Line.constraints.normalized]
# Plücker condition satisfied and direction normalized
condition = "direction_norm == 1"
constructors = [
    "from_points(a: Point<T>, b: Point<T>)",
    "from_planes(a: Plane<T>, b: Plane<T>)",
]

[types.Point]
grades = [3]
description = "Homogeneous point: x·e₀₂₃ + y·e₀₃₁ + z·e₀₁₂ + w·e₁₂₃"
fields = ["x", "y", "z", "w"]

[types.Point.constraints.normalized]
# w = 1 (Euclidean point, not at infinity)
condition = "w == 1"
constructors = [
    "from_xyz(x: T, y: T, z: T)",
    "origin()",
]

[types.Point.constraints.nonzero]
# Any non-zero point (includes ideal points)
condition = "norm_squared > epsilon"

[types.IdealPoint]
grades = [3]
description = "Point at infinity (direction): e₀₂₃, e₀₃₁, e₀₁₂"
fields = ["x", "y", "z"]
# w = 0 implicitly

[types.Motor]
grades = [0, 2]
description = "Rigid transformation (8 components)"
fields = ["s", "e01", "e02", "e03", "e12", "e31", "e23", "e0123"]

[types.Motor.constraints.unit]
# Motor norm = 1
norm = "motor"  # s² + e12² + e31² + e23² = 1, with ideal part constraint
constructors = [
    "identity()",
    "from_rotor(r: UnitRotor<T>)",
    "from_translator(t: Translator<T>)",
    "from_rotor_translator(r: UnitRotor<T>, t: Translator<T>)",
    "from_line_angle(line: Line<T>, angle: T)",
    "from_line_distance(line: Line<T>, distance: T)",
]
preserving_ops = ["compose", "inverse", "slerp"]

[types.Rotor]
grades = [0, 2]
description = "Pure rotation (4 Euclidean components only)"
fields = ["s", "e12", "e31", "e23"]
# e01, e02, e03, e0123 implicitly zero

[types.Rotor.constraints.unit]
norm = "euclidean"
constructors = [
    "identity()",
    "from_angle_plane(angle: T, plane: euclidean::dim3::Bivector<T>)",
]
preserving_ops = ["compose", "inverse"]

[types.Translator]
grades = [0, 2]
description = "Pure translation"
fields = ["s", "e01", "e02", "e03"]
# s = 1, e12 = e31 = e23 = e0123 = 0 implicitly

[types.Translator.constraints.unit]
# s = 1 always
constructors = [
    "identity()",
    "from_vector(v: euclidean::dim3::Vector<T>)",
    "from_xyz(x: T, y: T, z: T)",
]
preserving_ops = ["compose", "inverse"]
```

## Generated Output

### Module Structure

```
generated/
└── euclidean3/
    ├── mod.rs              # Module root with re-exports
    ├── types.rs            # Base type definitions (Rotor, Vector, etc.)
    ├── constrained.rs      # Constrained wrappers (UnitRotor, NonZeroVector, etc.)
    ├── products.rs         # All binary products
    ├── unary.rs            # Reverse, dual, conjugate, norms
    ├── ops.rs              # Operator overloads for base types
    ├── ops_constrained.rs  # Operator overloads for constrained types
    ├── conversions.rs      # To/from Multivector
    ├── approx.rs           # AbsDiffEq, RelativeEq, UlpsEq
    ├── arbitrary.rs        # proptest Arbitrary impls (base + constrained)
    ├── nalgebra.rs         # nalgebra interop (optional)
    └── tests.rs            # Consistency tests against Multivector
```

### Generated Base Type

```rust
// generated/euclidean3/types.rs

/// 3D rotor: scalar + bivector (even subalgebra).
///
/// Represents a rotation element in geometric algebra.
/// For guaranteed unit norm, use [`UnitRotor`].
///
/// # Components
///
/// | Field | Blade | Description |
/// |-------|-------|-------------|
/// | `s`   | 1     | Scalar part |
/// | `xy`  | e₁₂   | XY-plane component |
/// | `xz`  | e₁₃   | XZ-plane component |
/// | `yz`  | e₂₃   | YZ-plane component |
#[derive(Clone, Copy, Debug, PartialEq)]
#[repr(C)]
pub struct Rotor<T: Float> {
    s: T,
    xy: T,
    xz: T,
    yz: T,
}

impl<T: Float> Rotor<T> {
    /// Creates a new rotor from components.
    #[inline]
    pub fn new(s: T, xy: T, xz: T, yz: T) -> Self {
        Self { s, xy, xz, yz }
    }

    // Accessors
    #[inline] pub fn s(&self) -> T { self.s }
    #[inline] pub fn xy(&self) -> T { self.xy }
    #[inline] pub fn xz(&self) -> T { self.xz }
    #[inline] pub fn yz(&self) -> T { self.yz }

    /// Creates the identity rotor (scalar = 1, bivector = 0).
    #[inline]
    pub fn identity() -> Self {
        Self::new(T::one(), T::zero(), T::zero(), T::zero())
    }

    /// Creates the zero rotor.
    #[inline]
    pub fn zero() -> Self {
        Self::new(T::zero(), T::zero(), T::zero(), T::zero())
    }

    /// Returns the squared norm: s² + xy² + xz² + yz².
    #[inline]
    pub fn norm_squared(&self) -> T {
        self.s * self.s + self.xy * self.xy + self.xz * self.xz + self.yz * self.yz
    }

    /// Returns the norm.
    #[inline]
    pub fn norm(&self) -> T {
        self.norm_squared().sqrt()
    }

    /// Returns the reverse: s - xy·e₁₂ - xz·e₁₃ - yz·e₂₃.
    #[inline]
    pub fn reverse(&self) -> Self {
        Self::new(self.s, -self.xy, -self.xz, -self.yz)
    }

    /// Attempts to normalize this rotor.
    /// Returns `None` if the norm is too small.
    #[inline]
    pub fn try_normalize(&self) -> Option<Self> {
        let n = self.norm();
        if n < T::epsilon() {
            None
        } else {
            Some(Self::new(self.s / n, self.xy / n, self.xz / n, self.yz / n))
        }
    }

    /// Normalizes this rotor, panicking if zero.
    #[inline]
    pub fn normalize(&self) -> Self {
        self.try_normalize().expect("cannot normalize zero rotor")
    }
}
```

### Generated Constrained Type

```rust
// generated/euclidean3/constrained.rs

/// Unit rotor with guaranteed norm = 1.
///
/// This wrapper type ensures the rotor always has unit norm,
/// which is required for proper rotation operations.
///
/// # Construction
///
/// Use one of the safe constructors that guarantee unit norm:
///
/// ```
/// use clifford::specialized::euclidean::dim3::{UnitRotor, Bivector};
///
/// // From angle and plane (always unit)
/// let r = UnitRotor::from_angle_plane(0.5, Bivector::unit_xy());
///
/// // From two vectors (always unit)
/// let r = UnitRotor::from_vectors(Vector::unit_x(), Vector::unit_y());
///
/// // Try to normalize an existing rotor
/// let r = UnitRotor::try_from_rotor(some_rotor)?;
/// ```
///
/// # Operations
///
/// Unit rotors can be composed and inverted, always producing unit rotors:
///
/// ```
/// let r3: UnitRotor<f64> = r1.compose(r2);  // Still unit
/// let r_inv: UnitRotor<f64> = r1.inverse(); // Still unit
/// ```
#[derive(Clone, Copy, Debug, PartialEq)]
#[repr(transparent)]
pub struct UnitRotor<T: Float>(Rotor<T>);

impl<T: Float> UnitRotor<T> {
    // ========================================
    // Constructors (all guarantee unit norm)
    // ========================================

    /// Creates the identity rotation.
    #[inline]
    pub fn identity() -> Self {
        Self(Rotor::identity())
    }

    /// Creates a unit rotor from angle and rotation plane.
    ///
    /// The plane bivector does not need to be normalized; only its
    /// direction matters.
    ///
    /// # Arguments
    ///
    /// * `angle` - Rotation angle in radians
    /// * `plane` - Bivector defining the rotation plane
    #[inline]
    pub fn from_angle_plane(angle: T, plane: Bivector<T>) -> Self {
        let half = angle / T::TWO;
        let (sin_half, cos_half) = (half.sin(), half.cos());
        let plane_norm = plane.norm();

        if plane_norm < T::epsilon() {
            return Self::identity();
        }

        let scale = sin_half / plane_norm;
        Self(Rotor::new(
            cos_half,
            plane.xy() * scale,
            plane.xz() * scale,
            plane.yz() * scale,
        ))
    }

    /// Creates a unit rotor that rotates vector `a` to vector `b`.
    #[inline]
    pub fn from_vectors(a: Vector<T>, b: Vector<T>) -> Self {
        // ... generated implementation ...
    }

    /// Attempts to create a unit rotor by normalizing the given rotor.
    ///
    /// Returns `None` if the rotor's norm is too small.
    #[inline]
    pub fn try_from_rotor(r: Rotor<T>) -> Option<Self> {
        r.try_normalize().map(Self)
    }

    /// Creates a unit rotor without checking the constraint.
    ///
    /// # Safety (Logical)
    ///
    /// The caller must ensure `s² + xy² + xz² + yz² = 1`.
    /// This is useful for internal operations known to preserve unit norm.
    #[inline]
    pub fn new_unchecked(s: T, xy: T, xz: T, yz: T) -> Self {
        Self(Rotor::new(s, xy, xz, yz))
    }

    // ========================================
    // Constraint-Preserving Operations
    // ========================================

    /// Composes two rotations: `self` followed by `other`.
    ///
    /// The result is always a unit rotor.
    #[inline]
    pub fn compose(self, other: Self) -> Self {
        // Generated: unit * unit = unit
        Self::new_unchecked(
            // ... generated formula ...
        )
    }

    /// Returns the inverse rotation.
    ///
    /// For unit rotors, the inverse equals the reverse.
    #[inline]
    pub fn inverse(self) -> Self {
        Self(self.0.reverse())
    }

    /// Spherical linear interpolation between two rotors.
    ///
    /// The result is always a unit rotor.
    #[inline]
    pub fn slerp(self, other: Self, t: T) -> Self {
        // ... generated implementation ...
    }

    // ========================================
    // Access
    // ========================================

    /// Returns the underlying rotor.
    #[inline]
    pub fn into_inner(self) -> Rotor<T> {
        self.0
    }

    /// Returns a reference to the underlying rotor.
    #[inline]
    pub fn as_rotor(&self) -> &Rotor<T> {
        &self.0
    }
}

// Deref to Rotor for accessing fields and non-mutating methods
impl<T: Float> std::ops::Deref for UnitRotor<T> {
    type Target = Rotor<T>;

    #[inline]
    fn deref(&self) -> &Self::Target {
        &self.0
    }
}

impl<T: Float> AsRef<Rotor<T>> for UnitRotor<T> {
    #[inline]
    fn as_ref(&self) -> &Rotor<T> {
        &self.0
    }
}

impl<T: Float> From<UnitRotor<T>> for Rotor<T> {
    #[inline]
    fn from(u: UnitRotor<T>) -> Self {
        u.0
    }
}

// Unit rotor multiplication produces unit rotor
impl<T: Float> std::ops::Mul for UnitRotor<T> {
    type Output = Self;

    #[inline]
    fn mul(self, rhs: Self) -> Self {
        self.compose(rhs)
    }
}
```

### Generated Products with Constraint Awareness

```rust
// generated/euclidean3/products.rs

/// Geometric product: Rotor × Rotor → Rotor
#[inline]
pub fn geometric_rotor_rotor<T: Float>(a: Rotor<T>, b: Rotor<T>) -> Rotor<T> {
    // Generated formula
    Rotor::new(
        a.s()*b.s() - a.xy()*b.xy() - a.xz()*b.xz() - a.yz()*b.yz(),
        a.s()*b.xy() + a.xy()*b.s() + a.xz()*b.yz() - a.yz()*b.xz(),
        a.s()*b.xz() + a.xz()*b.s() - a.xy()*b.yz() + a.yz()*b.xy(),
        a.s()*b.yz() + a.yz()*b.s() + a.xy()*b.xz() - a.xz()*b.xy(),
    )
}

/// Geometric product: UnitRotor × UnitRotor → UnitRotor
///
/// This specialized version knows the output is unit because
/// the product of two unit rotors is always unit.
#[inline]
pub fn geometric_unit_rotor_unit_rotor<T: Float>(a: UnitRotor<T>, b: UnitRotor<T>) -> UnitRotor<T> {
    UnitRotor::new_unchecked(
        a.s()*b.s() - a.xy()*b.xy() - a.xz()*b.xz() - a.yz()*b.yz(),
        a.s()*b.xy() + a.xy()*b.s() + a.xz()*b.yz() - a.yz()*b.xz(),
        a.s()*b.xz() + a.xz()*b.s() - a.xy()*b.yz() + a.yz()*b.xy(),
        a.s()*b.yz() + a.yz()*b.s() + a.xy()*b.xz() - a.xz()*b.xy(),
    )
}

/// Sandwich product: UnitRotor × Vector × UnitRotor† → Vector
///
/// Applies the rotation to a vector.
#[inline]
pub fn sandwich_unit_rotor_vector<T: Float>(r: UnitRotor<T>, v: Vector<T>) -> Vector<T> {
    // Optimized sandwich formula (not r * v * r.reverse())
    // Generated from symbolic expansion
}
```

## Semantic Layer

Generated code provides the **algebraic foundation**. Semantic meaning is added via extension traits or wrapper methods in separate (non-generated) files:

```rust
// src/specialized/euclidean/dim3/semantic.rs (hand-written, not generated)

use crate::generated::euclidean3::{UnitRotor, Vector, sandwich_unit_rotor_vector};

impl<T: Float> UnitRotor<T> {
    /// Rotates a vector by this rotation.
    ///
    /// Computes the sandwich product R̃ v R, which rotates v
    /// in the plane and by the angle encoded in this rotor.
    ///
    /// # Example
    ///
    /// ```
    /// use clifford::specialized::euclidean::dim3::{UnitRotor, Vector, Bivector};
    /// use std::f64::consts::FRAC_PI_2;
    ///
    /// let r = UnitRotor::from_angle_plane(FRAC_PI_2, Bivector::unit_xy());
    /// let v = Vector::unit_x();
    /// let rotated = r.rotate(v);
    /// // rotated ≈ Vector::unit_y()
    /// ```
    #[inline]
    pub fn rotate(&self, v: Vector<T>) -> Vector<T> {
        sandwich_unit_rotor_vector(*self, v)
    }
}
```

## Architecture

```
clifford-codegen/
├── Cargo.toml
├── src/
│   ├── main.rs                 # CLI entry point
│   ├── lib.rs                  # Library interface
│   │
│   ├── spec/                   # Specification parsing (PRD-14.2)
│   │   ├── mod.rs
│   │   ├── algebra.rs          # Top-level algebra definition
│   │   ├── signature.rs        # Signature parsing
│   │   ├── types.rs            # Type definitions
│   │   └── constraints.rs      # Constraint definitions
│   │
│   ├── algebra/                # Algebraic computations (PRD-14.1)
│   │   ├── mod.rs
│   │   ├── blade.rs            # Blade representation & ordering
│   │   ├── basis_product.rs    # Basis blade multiplication
│   │   ├── sign.rs             # Sign computation
│   │   ├── grade.rs            # Grade utilities
│   │   └── metric.rs           # Metric signature handling
│   │
│   ├── codegen/                # Code generation
│   │   ├── mod.rs
│   │   ├── types.rs            # Type generation (PRD-14.3)
│   │   ├── constrained.rs      # Constrained types (PRD-14.5)
│   │   ├── products.rs         # Product generation (PRD-14.4)
│   │   ├── unary.rs            # Unary ops generation (PRD-14.4)
│   │   ├── traits.rs           # Trait impl generation (PRD-14.6)
│   │   ├── arbitrary.rs        # Arbitrary generation (PRD-14.6)
│   │   └── tests.rs            # Test generation (PRD-14.6)
│   │
│   └── emit/                   # Rust code emission
│       ├── mod.rs
│       ├── rust.rs             # Rust AST building
│       └── format.rs           # rustfmt integration
│
├── algebras/                   # Bundled algebra specifications
│   ├── euclidean2.toml
│   ├── euclidean3.toml
│   ├── pga2.toml
│   ├── pga3.toml
│   ├── cga2.toml
│   └── cga3.toml
│
└── tests/
    ├── algebra_tests.rs        # Blade algebra tests
    ├── codegen_tests.rs        # Code generation tests
    └── integration/            # Full integration tests
        ├── euclidean3.rs
        └── pga3.rs
```

## Verification Strategy

### 1. Blade Algebra Verification (PRD-14.1)

```rust
proptest! {
    // Verify sign computation is symmetric/antisymmetric correctly
    #[test]
    fn basis_product_sign_anticommutes_for_vectors(i in 0usize..3, j in 0usize..3) {
        prop_assume!(i != j);
        let blade_i = 1 << i;
        let blade_j = 1 << j;
        let (sign_ij, _) = basis_product(blade_i, blade_j, euclidean_metric);
        let (sign_ji, _) = basis_product(blade_j, blade_i, euclidean_metric);
        prop_assert_eq!(sign_ij, -sign_ji);
    }
}
```

### 2. Product Consistency (PRD-14.4)

```rust
proptest! {
    #[test]
    fn geometric_product_matches_multivector(
        a in any::<Rotor<f64>>(),
        b in any::<Rotor<f64>>(),
    ) {
        let gen_result = geometric_rotor_rotor(a, b);
        let mv_a = Multivector::<f64, Euclidean3>::from(a);
        let mv_b = Multivector::<f64, Euclidean3>::from(b);
        let mv_result = mv_a * mv_b;
        let mv_as_rotor = Rotor::try_from(mv_result).unwrap();
        prop_assert!(abs_diff_eq!(gen_result, mv_as_rotor, epsilon = ABS_DIFF_EQ_EPS));
    }
}
```

### 3. Constraint Preservation (PRD-14.5)

```rust
proptest! {
    #[test]
    fn unit_rotor_composition_preserves_unit(
        a in any::<UnitRotor<f64>>(),
        b in any::<UnitRotor<f64>>(),
    ) {
        let result = a.compose(b);
        prop_assert!(abs_diff_eq!(result.norm(), 1.0, epsilon = ABS_DIFF_EQ_EPS));
    }

    #[test]
    fn unit_rotor_inverse_preserves_unit(r in any::<UnitRotor<f64>>()) {
        let inv = r.inverse();
        prop_assert!(abs_diff_eq!(inv.norm(), 1.0, epsilon = ABS_DIFF_EQ_EPS));
    }

    #[test]
    fn unit_rotor_slerp_preserves_unit(
        a in any::<UnitRotor<f64>>(),
        b in any::<UnitRotor<f64>>(),
        t in 0.0f64..1.0,
    ) {
        let result = a.slerp(b, t);
        prop_assert!(abs_diff_eq!(result.norm(), 1.0, epsilon = ABS_DIFF_EQ_EPS));
    }
}
```

## Success Criteria

1. **Correctness**: All generated products match `Multivector` within floating-point tolerance
2. **Completeness**: Every valid product combination is generated
3. **Type Safety**: Constrained types enforce invariants at compile time
4. **Performance**: Generated code performs identically to hand-written (verified via benchmarks)
5. **Ergonomics**: Semantic layer provides intuitive API
6. **Maintainability**: Adding new algebra takes < 1 hour

## Open Questions

1. **Proc macro vs CLI**: Should we pursue a proc macro for compile-time generation?
2. **Grade selection**: How to handle products that produce new grade combinations not in the type map?
3. **Constraint composition**: Can we automatically infer which operations preserve which constraints?
4. **SIMD generation**: Should the generator have SIMD-aware output modes?
5. **Documentation generation**: How much mathematical documentation should be auto-generated?

## References

- [ganja.js](https://github.com/enkimute/ganja.js) - JavaScript GA with code generation
- [Klein](https://github.com/jeremyong/klein) - C++ PGA library with SIMD, uses templates
- [geometric_algebra](https://github.com/weshoke/geometric_algebra) - Rust GA with const generics
- [galern](https://github.com/Checkmate50/galern) - Lean 4 GA formalization
