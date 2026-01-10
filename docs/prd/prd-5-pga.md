# PRD-5: Projective Geometric Algebra (PGA)

**Status**: In Progress (Phases 1-3 Complete)
**Goal**: Point-based PGA for representing points, lines, planes, and rigid transformations

## Background

Projective Geometric Algebra (PGA) extends Euclidean GA by adding a degenerate (null) basis vector `e₀` that squares to zero (`e₀² = 0`). This enables elegant representation of:

- **Points** as grade-1 vectors (homogeneous coordinates)
- **Lines** as grade-2 bivectors
- **Planes** as grade-3 trivectors (in 3D)
- **Rigid body transforms** (rotation + translation) as motors

### Point-Based vs Plane-Based PGA

There are two equivalent formulations of PGA, related by the Hodge dual:

| Formulation | Points | Lines | Planes | Convention |
|-------------|--------|-------|--------|------------|
| **Point-based** | Grade 1 | Grade 2 | Grade 3 | Vectors are points |
| Plane-based | Grade 3 | Grade 2 | Grade 1 | Vectors are planes |

**This PRD uses the point-based formulation** for consistency with:
- Euclidean GA where vectors represent directions/positions
- nalgebra's conventions where `Point3` and `Vector3` are grade-1 objects

### Signature

| Algebra | Signature | Dimension | Blades |
|---------|-----------|-----------|--------|
| 2D PGA | `Cl(2,0,1)` | 3 | 8 |
| 3D PGA | `Cl(3,0,1)` | 4 | 16 |

Basis vectors:
- 2D: `e₁, e₂, e₀` where `e₁² = e₂² = 1`, `e₀² = 0`
- 3D: `e₁, e₂, e₃, e₀` where `e₁² = e₂² = e₃² = 1`, `e₀² = 0`

## Implementation Strategy

**Incremental approach with verification at each step:**

1. **Phase 1: Generic Multivector with PGA Signatures**
   - Add `Projective2` and `Projective3` signatures
   - Property test mathematical laws (associativity, distributivity, etc.)
   - Test against nalgebra for consistency where applicable
   - Commit when all tests pass

2. **Phase 2: Specialized 2D Types** (`specialized/projective/dim2/`)
   - Point, Line, Motor types
   - Test against Multivector for consistency
   - Test operations against nalgebra equivalents
   - Commit when all tests pass

3. **Phase 3: Specialized 3D Types** (`specialized/projective/dim3/`)
   - Point, Line, Plane, Motor types
   - Test against Multivector for consistency
   - Test operations against nalgebra equivalents
   - Commit when all tests pass

## Phase 1: PGA Signatures

### Deliverables

#### 1. Signature Types (`src/signature/projective.rs`)

```rust
/// 2D Projective GA: Cl(2,0,1)
///
/// Point-based formulation where:
/// - Grade 1 (vectors): Points in homogeneous coordinates
/// - Grade 2 (bivectors): Lines
/// - Grade 3 (pseudoscalar): Oriented area
///
/// Basis: e₁, e₂, e₀ where e₀² = 0 (null/degenerate)
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct Projective2;

impl Signature for Projective2 {
    type NumBlades = typenum::U8; // 2³ = 8

    const P: usize = 2;
    const Q: usize = 0;
    const R: usize = 1;

    fn metric(i: usize) -> i8 {
        match i {
            0 | 1 => 1,  // e₁² = e₂² = 1
            2 => 0,      // e₀² = 0
            _ => panic!("basis index {i} out of range for Projective2"),
        }
    }
}

/// 3D Projective GA: Cl(3,0,1)
///
/// Point-based formulation where:
/// - Grade 1 (vectors): Points in homogeneous coordinates
/// - Grade 2 (bivectors): Lines
/// - Grade 3 (trivectors): Planes
/// - Grade 4 (pseudoscalar): Oriented volume
///
/// Basis: e₁, e₂, e₃, e₀ where e₀² = 0 (null/degenerate)
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct Projective3;

impl Signature for Projective3 {
    type NumBlades = typenum::U16; // 2⁴ = 16

    const P: usize = 3;
    const Q: usize = 0;
    const R: usize = 1;

    fn metric(i: usize) -> i8 {
        match i {
            0 | 1 | 2 => 1,  // e₁² = e₂² = e₃² = 1
            3 => 0,          // e₀² = 0
            _ => panic!("basis index {i} out of range for Projective3"),
        }
    }
}
```

### Phase 1 Testing

#### Mathematical Properties (proptest)

```rust
use crate::signature::{Projective2, Projective3};

proptest! {
    // ========================================================================
    // Algebraic properties for PGA2
    // ========================================================================

    #[test]
    fn pga2_geometric_product_associative(
        a in any::<Multivector<f64, Projective2>>(),
        b in any::<Multivector<f64, Projective2>>(),
        c in any::<Multivector<f64, Projective2>>(),
    ) {
        let lhs = &(&a * &b) * &c;
        let rhs = &a * &(&b * &c);
        prop_assert!(abs_diff_eq!(lhs, rhs, epsilon = ABS_DIFF_EQ_EPS));
    }

    #[test]
    fn pga2_geometric_product_distributive(
        a in any::<Multivector<f64, Projective2>>(),
        b in any::<Multivector<f64, Projective2>>(),
        c in any::<Multivector<f64, Projective2>>(),
    ) {
        let lhs = &a * &(&b + &c);
        let rhs = &(&a * &b) + &(&a * &c);
        prop_assert!(abs_diff_eq!(lhs, rhs, epsilon = ABS_DIFF_EQ_EPS));
    }

    #[test]
    fn pga2_null_vector_squares_to_zero(s in -100.0f64..100.0) {
        // e₀ component only
        let null_vec: Multivector<f64, Projective2> = {
            let mut mv = Multivector::zero();
            mv.set(Blade::basis_vector(2), s); // e₀
            mv
        };
        let sq = &null_vec * &null_vec;
        prop_assert!(sq.is_zero(ABS_DIFF_EQ_EPS));
    }

    #[test]
    fn pga2_reverse_involutory(a in any::<Multivector<f64, Projective2>>()) {
        prop_assert!(abs_diff_eq!(a.reverse().reverse(), a, epsilon = ABS_DIFF_EQ_EPS));
    }

    #[test]
    fn pga2_outer_associative(
        a in any::<Multivector<f64, Projective2>>(),
        b in any::<Multivector<f64, Projective2>>(),
        c in any::<Multivector<f64, Projective2>>(),
    ) {
        let lhs = a.outer(&b).outer(&c);
        let rhs = a.outer(&b.outer(&c));
        prop_assert!(abs_diff_eq!(lhs, rhs, epsilon = ABS_DIFF_EQ_EPS));
    }

    // ========================================================================
    // Algebraic properties for PGA3
    // ========================================================================

    #[test]
    fn pga3_geometric_product_associative(
        a in any::<Multivector<f64, Projective3>>(),
        b in any::<Multivector<f64, Projective3>>(),
        c in any::<Multivector<f64, Projective3>>(),
    ) {
        let lhs = &(&a * &b) * &c;
        let rhs = &a * &(&b * &c);
        prop_assert!(abs_diff_eq!(lhs, rhs, epsilon = ABS_DIFF_EQ_EPS));
    }

    #[test]
    fn pga3_null_vector_squares_to_zero(s in -100.0f64..100.0) {
        // e₀ component only
        let null_vec: Multivector<f64, Projective3> = {
            let mut mv = Multivector::zero();
            mv.set(Blade::basis_vector(3), s); // e₀
            mv
        };
        let sq = &null_vec * &null_vec;
        prop_assert!(sq.is_zero(ABS_DIFF_EQ_EPS));
    }

    #[test]
    fn pga3_reverse_antimorphism(
        a in any::<Multivector<f64, Projective3>>(),
        b in any::<Multivector<f64, Projective3>>(),
    ) {
        let lhs = (&a * &b).reverse();
        let rhs = &b.reverse() * &a.reverse();
        prop_assert!(abs_diff_eq!(lhs, rhs, epsilon = ABS_DIFF_EQ_EPS));
    }
}
```

#### nalgebra Consistency Tests

Test that PGA operations are consistent with nalgebra where applicable:

```rust
#[cfg(any(feature = "nalgebra-0_33", feature = "nalgebra-0_34"))]
mod nalgebra_consistency {
    use super::*;
    use nalgebra as na;

    proptest! {
        /// Euclidean subspace dot product matches nalgebra
        #[test]
        fn pga3_euclidean_dot_matches_nalgebra(
            ax in -100.0f64..100.0, ay in -100.0f64..100.0, az in -100.0f64..100.0,
            bx in -100.0f64..100.0, by in -100.0f64..100.0, bz in -100.0f64..100.0,
        ) {
            // Pure Euclidean vectors (no e₀ component)
            let pga_a: Multivector<f64, Projective3> =
                Multivector::vector(&[ax, ay, az, 0.0]);
            let pga_b: Multivector<f64, Projective3> =
                Multivector::vector(&[bx, by, bz, 0.0]);

            let na_a = na::Vector3::new(ax, ay, az);
            let na_b = na::Vector3::new(bx, by, bz);

            let pga_dot = pga_a.inner(&pga_b).scalar_part();
            let na_dot = na_a.dot(&na_b);

            prop_assert!(abs_diff_eq!(pga_dot, na_dot, epsilon = ABS_DIFF_EQ_EPS));
        }

        /// Euclidean subspace geometric product gives correct scalar + bivector
        #[test]
        fn pga3_euclidean_vectors_product(
            ax in -10.0f64..10.0, ay in -10.0f64..10.0, az in -10.0f64..10.0,
            bx in -10.0f64..10.0, by in -10.0f64..10.0, bz in -10.0f64..10.0,
        ) {
            let pga_a: Multivector<f64, Projective3> =
                Multivector::vector(&[ax, ay, az, 0.0]);
            let pga_b: Multivector<f64, Projective3> =
                Multivector::vector(&[bx, by, bz, 0.0]);

            let na_a = na::Vector3::new(ax, ay, az);
            let na_b = na::Vector3::new(bx, by, bz);

            // Geometric product of vectors: ab = a·b + a∧b
            let product = &pga_a * &pga_b;

            // Scalar part should be dot product
            let expected_scalar = na_a.dot(&na_b);
            prop_assert!(abs_diff_eq!(
                product.scalar_part(),
                expected_scalar,
                epsilon = ABS_DIFF_EQ_EPS
            ));

            // Should have no vector or trivector parts
            prop_assert!(product.grade_select(1).is_zero(ABS_DIFF_EQ_EPS));
            prop_assert!(product.grade_select(3).is_zero(ABS_DIFF_EQ_EPS));
        }
    }
}
```

### Files to Create (Phase 1)

- `src/signature/projective.rs` - Signature definitions
- `src/algebra/arbitrary.rs` - Add Arbitrary impls for PGA multivectors
- Update `src/signature/mod.rs` - Re-export new signatures

### Phase 1 Verification

- [x] `cargo fmt` passes
- [x] `cargo clippy --all-features` passes
- [x] `cargo doc --all-features --no-deps` passes
- [x] `cargo test --all-features` - all property tests pass
- [x] `cargo deny check` passes

---

## Phase 2: Specialized 2D Types

### Module Structure

```
src/specialized/projective/
├── mod.rs              # Module root
└── dim2/
    ├── mod.rs          # Public API
    ├── types.rs        # Point, Line, Motor
    ├── ops.rs          # Operations
    ├── conversions.rs  # To/from Multivector
    ├── arbitrary.rs    # Proptest support
    └── nalgebra.rs     # nalgebra conversions (feature-gated)
```

### Type Definitions (`specialized/projective/dim2/types.rs`)

```rust
use crate::scalar::Float;

/// A point in 2D PGA (grade 1 vector).
///
/// In point-based PGA, a point is represented as a homogeneous coordinate:
/// `P = x·e₁ + y·e₂ + w·e₀`
///
/// where `(x/w, y/w)` are the Cartesian coordinates when `w ≠ 0`.
/// Points with `w = 0` represent ideal points (points at infinity).
///
/// # Example
///
/// ```
/// use clifford::specialized::projective::dim2::Point;
///
/// // Point at (3, 4)
/// let p = Point::new(3.0, 4.0);
/// assert_eq!(p.x(), 3.0);
/// assert_eq!(p.y(), 4.0);
///
/// // Point at infinity in direction (1, 0)
/// let ideal = Point::ideal(1.0, 0.0);
/// assert!(ideal.is_ideal(1e-10));
/// ```
#[derive(Clone, Copy, Debug, PartialEq)]
#[repr(C)]
pub struct Point<T: Float> {
    /// Coefficient of e₁ (x-coordinate × w)
    pub e1: T,
    /// Coefficient of e₂ (y-coordinate × w)
    pub e2: T,
    /// Coefficient of e₀ (homogeneous weight)
    pub e0: T,
}

impl<T: Float> Point<T> {
    /// Creates a finite point at Cartesian coordinates (x, y).
    #[inline]
    pub fn new(x: T, y: T) -> Self {
        Self { e1: x, e2: y, e0: T::one() }
    }

    /// Creates a point from homogeneous coordinates.
    #[inline]
    pub fn from_homogeneous(e1: T, e2: T, e0: T) -> Self {
        Self { e1, e2, e0 }
    }

    /// Creates an ideal point (point at infinity) in the given direction.
    #[inline]
    pub fn ideal(dx: T, dy: T) -> Self {
        Self { e1: dx, e2: dy, e0: T::zero() }
    }

    /// Origin point (0, 0).
    #[inline]
    pub fn origin() -> Self {
        Self::new(T::zero(), T::zero())
    }

    /// Returns the x-coordinate (requires w ≠ 0).
    #[inline]
    pub fn x(&self) -> T {
        self.e1 / self.e0
    }

    /// Returns the y-coordinate (requires w ≠ 0).
    #[inline]
    pub fn y(&self) -> T {
        self.e2 / self.e0
    }

    /// Returns true if this is an ideal point (point at infinity).
    #[inline]
    pub fn is_ideal(&self, epsilon: T) -> bool {
        self.e0.abs() < epsilon
    }

    /// Normalizes the homogeneous coordinates so w = 1 (if finite).
    pub fn normalize(&self) -> Option<Self> {
        if self.e0.abs() < T::epsilon() {
            None
        } else {
            Some(Self {
                e1: self.e1 / self.e0,
                e2: self.e2 / self.e0,
                e0: T::one(),
            })
        }
    }
}

/// A line in 2D PGA (grade 2 bivector).
///
/// In point-based PGA, a line is a bivector `L = a·e₁₂ + b·e₂₀ + c·e₀₁`
/// representing the line `ax + by + c = 0` in implicit form.
///
/// # Example
///
/// ```
/// use clifford::specialized::projective::dim2::{Point, Line};
///
/// let p1 = Point::new(0.0, 0.0);
/// let p2 = Point::new(1.0, 1.0);
///
/// // Line through two points
/// let line = p1.join(&p2);
/// ```
#[derive(Clone, Copy, Debug, PartialEq)]
#[repr(C)]
pub struct Line<T: Float> {
    /// Coefficient of e₁₂ (related to line's distance from origin)
    pub e12: T,
    /// Coefficient of e₂₀ (x-component of normal, scaled)
    pub e20: T,
    /// Coefficient of e₀₁ (y-component of normal, scaled)
    pub e01: T,
}

impl<T: Float> Line<T> {
    /// Creates a line from bivector coefficients.
    #[inline]
    pub fn new(e12: T, e20: T, e01: T) -> Self {
        Self { e12, e20, e01 }
    }

    /// Creates a line from implicit equation ax + by + c = 0.
    #[inline]
    pub fn from_implicit(a: T, b: T, c: T) -> Self {
        Self { e12: c, e20: a, e01: b }
    }

    /// Creates the x-axis (y = 0).
    #[inline]
    pub fn x_axis() -> Self {
        Self::from_implicit(T::zero(), T::one(), T::zero())
    }

    /// Creates the y-axis (x = 0).
    #[inline]
    pub fn y_axis() -> Self {
        Self::from_implicit(T::one(), T::zero(), T::zero())
    }

    /// Returns the normal direction (a, b) of the line ax + by + c = 0.
    #[inline]
    pub fn normal(&self) -> (T, T) {
        (self.e20, self.e01)
    }
}

/// A motor (rigid transformation) in 2D PGA.
///
/// Motors represent rigid body transformations (rotation + translation).
/// In 2D PGA, a motor is an even-grade element: `M = s + d·e₁₂ + tx·e₂₀ + ty·e₀₁`
///
/// # Example
///
/// ```
/// use clifford::specialized::projective::dim2::{Point, Motor};
/// use std::f64::consts::FRAC_PI_2;
///
/// // 90° rotation around origin
/// let rotation = Motor::from_rotation(FRAC_PI_2);
///
/// // Translation by (1, 2)
/// let translation = Motor::from_translation(1.0, 2.0);
///
/// // Compose: first rotate, then translate
/// let combined = translation.compose(&rotation);
/// ```
#[derive(Clone, Copy, Debug, PartialEq)]
#[repr(C)]
pub struct Motor<T: Float> {
    /// Scalar part (cos(θ/2) for rotation)
    pub s: T,
    /// Coefficient of e₁₂ (sin(θ/2) for rotation)
    pub e12: T,
    /// Coefficient of e₂₀ (translation x component / 2)
    pub e20: T,
    /// Coefficient of e₀₁ (translation y component / 2)
    pub e01: T,
}

impl<T: Float> Motor<T> {
    /// Creates the identity motor (no transformation).
    #[inline]
    pub fn identity() -> Self {
        Self {
            s: T::one(),
            e12: T::zero(),
            e20: T::zero(),
            e01: T::zero(),
        }
    }

    /// Creates a pure rotation motor around the origin.
    #[inline]
    pub fn from_rotation(angle: T) -> Self {
        let half = angle / T::TWO;
        Self {
            s: half.cos(),
            e12: half.sin(),
            e20: T::zero(),
            e01: T::zero(),
        }
    }

    /// Creates a pure translation motor.
    #[inline]
    pub fn from_translation(dx: T, dy: T) -> Self {
        Self {
            s: T::one(),
            e12: T::zero(),
            e20: dx / T::TWO,
            e01: dy / T::TWO,
        }
    }

    /// Transforms a point: M P M̃
    pub fn transform_point(&self, p: &Point<T>) -> Point<T>;

    /// Transforms a line: M L M̃
    pub fn transform_line(&self, l: &Line<T>) -> Line<T>;

    /// Composes two motors: M₂ ∘ M₁ = M₂ M₁
    pub fn compose(&self, other: &Self) -> Self;

    /// Returns the inverse motor: M⁻¹
    pub fn inverse(&self) -> Self;

    /// Linear interpolation between motors.
    pub fn lerp(&self, other: &Self, t: T) -> Self;
}
```

### Operations (`specialized/projective/dim2/ops.rs`)

```rust
impl<T: Float> Point<T> {
    /// Join of two points: the line through them.
    ///
    /// P₁ ∨ P₂ (regressive product in point-based PGA)
    pub fn join(&self, other: &Point<T>) -> Line<T>;

    /// Distance to another point (Euclidean distance).
    pub fn distance(&self, other: &Point<T>) -> T;
}

impl<T: Float> Line<T> {
    /// Meet of two lines: their intersection point.
    ///
    /// L₁ ∧ L₂ (outer product in point-based PGA)
    pub fn meet(&self, other: &Line<T>) -> Point<T>;

    /// Distance from a point to this line.
    pub fn distance_to_point(&self, p: &Point<T>) -> T;

    /// Project a point onto this line.
    pub fn project(&self, p: &Point<T>) -> Point<T>;
}
```

### Phase 2 Testing

```rust
proptest! {
    // ========================================================================
    // Consistency with Multivector
    // ========================================================================

    #[test]
    fn point_roundtrip(
        x in -100.0f64..100.0,
        y in -100.0f64..100.0
    ) {
        let point = Point::new(x, y);
        let mv = Multivector::<f64, Projective2>::from(point);
        let back = Point::try_from(mv).unwrap();
        prop_assert!(abs_diff_eq!(point.e1, back.e1, epsilon = ABS_DIFF_EQ_EPS));
        prop_assert!(abs_diff_eq!(point.e2, back.e2, epsilon = ABS_DIFF_EQ_EPS));
        prop_assert!(abs_diff_eq!(point.e0, back.e0, epsilon = ABS_DIFF_EQ_EPS));
    }

    #[test]
    fn join_consistency(
        x1 in -10.0f64..10.0, y1 in -10.0f64..10.0,
        x2 in -10.0f64..10.0, y2 in -10.0f64..10.0,
    ) {
        let p1 = Point::new(x1, y1);
        let p2 = Point::new(x2, y2);

        // Specialized join
        let spec_line = p1.join(&p2);

        // Generic regressive product
        let mv1 = Multivector::<f64, Projective2>::from(p1);
        let mv2 = Multivector::<f64, Projective2>::from(p2);
        let gen_line = mv1.regressive(&mv2);
        let gen_as_line = Line::try_from(gen_line).unwrap();

        // Should match (up to scale)
        // Normalize both for comparison
        let spec_norm = (spec_line.e12.powi(2) + spec_line.e20.powi(2) + spec_line.e01.powi(2)).sqrt();
        let gen_norm = (gen_as_line.e12.powi(2) + gen_as_line.e20.powi(2) + gen_as_line.e01.powi(2)).sqrt();

        if spec_norm > ABS_DIFF_EQ_EPS && gen_norm > ABS_DIFF_EQ_EPS {
            let spec_normalized = Line::new(
                spec_line.e12 / spec_norm,
                spec_line.e20 / spec_norm,
                spec_line.e01 / spec_norm,
            );
            let gen_normalized = Line::new(
                gen_as_line.e12 / gen_norm,
                gen_as_line.e20 / gen_norm,
                gen_as_line.e01 / gen_norm,
            );
            // Lines equal up to sign
            let same = abs_diff_eq!(spec_normalized, gen_normalized, epsilon = ABS_DIFF_EQ_EPS);
            let opposite = abs_diff_eq!(
                spec_normalized,
                Line::new(-gen_normalized.e12, -gen_normalized.e20, -gen_normalized.e01),
                epsilon = ABS_DIFF_EQ_EPS
            );
            prop_assert!(same || opposite);
        }
    }

    #[test]
    fn motor_transform_consistency(
        angle in -std::f64::consts::PI..std::f64::consts::PI,
        tx in -10.0f64..10.0, ty in -10.0f64..10.0,
        px in -10.0f64..10.0, py in -10.0f64..10.0,
    ) {
        let motor = Motor::from_translation(tx, ty).compose(&Motor::from_rotation(angle));
        let point = Point::new(px, py);

        // Specialized transform
        let spec_result = motor.transform_point(&point);

        // Generic sandwich product
        let mv_motor = Multivector::<f64, Projective2>::from(motor);
        let mv_point = Multivector::<f64, Projective2>::from(point);
        let gen_result = mv_motor.sandwich(&mv_point);
        let gen_as_point = Point::try_from(gen_result).unwrap();

        prop_assert!(abs_diff_eq!(
            spec_result.normalize().unwrap(),
            gen_as_point.normalize().unwrap(),
            epsilon = ABS_DIFF_EQ_EPS
        ));
    }

    // ========================================================================
    // Motor properties
    // ========================================================================

    #[test]
    fn motor_preserves_distance(
        motor in any::<UnitMotor2<f64>>(),
        p1 in any::<Point<f64>>(),
        p2 in any::<Point<f64>>(),
    ) {
        let d_before = p1.distance(&p2);
        let t1 = motor.transform_point(&p1);
        let t2 = motor.transform_point(&p2);
        let d_after = t1.distance(&t2);
        prop_assert!(abs_diff_eq!(d_before, d_after, epsilon = ABS_DIFF_EQ_EPS));
    }

    #[test]
    fn motor_composition_associative(
        m1 in any::<Motor<f64>>(),
        m2 in any::<Motor<f64>>(),
        m3 in any::<Motor<f64>>(),
    ) {
        let lhs = m1.compose(&m2).compose(&m3);
        let rhs = m1.compose(&m2.compose(&m3));
        prop_assert!(abs_diff_eq!(lhs, rhs, epsilon = ABS_DIFF_EQ_EPS));
    }

    #[test]
    fn motor_inverse(
        motor in any::<UnitMotor2<f64>>(),
        point in any::<Point<f64>>(),
    ) {
        let transformed = motor.transform_point(&point);
        let back = motor.inverse().transform_point(&transformed);
        prop_assert!(abs_diff_eq!(
            point.normalize().unwrap(),
            back.normalize().unwrap(),
            epsilon = ABS_DIFF_EQ_EPS
        ));
    }

    // ========================================================================
    // Join/Meet duality
    // ========================================================================

    #[test]
    fn join_meet_duality(
        x1 in -10.0f64..10.0, y1 in -10.0f64..10.0,
        x2 in -10.0f64..10.0, y2 in -10.0f64..10.0,
        x3 in -10.0f64..10.0, y3 in -10.0f64..10.0,
    ) {
        let p1 = Point::new(x1, y1);
        let p2 = Point::new(x2, y2);
        let p3 = Point::new(x3, y3);

        // Line through p1 and p2
        let line12 = p1.join(&p2);
        // Line through p1 and p3
        let line13 = p1.join(&p3);

        // Meet of two lines through p1 should give p1
        let meet = line12.meet(&line13);

        if let (Some(p1_norm), Some(meet_norm)) = (p1.normalize(), meet.normalize()) {
            prop_assert!(abs_diff_eq!(p1_norm, meet_norm, epsilon = ABS_DIFF_EQ_EPS));
        }
    }
}
```

### nalgebra Consistency Tests (Phase 2)

```rust
#[cfg(any(feature = "nalgebra-0_33", feature = "nalgebra-0_34"))]
mod nalgebra_tests {
    proptest! {
        #[test]
        fn point_to_nalgebra_roundtrip(x in -100.0f64..100.0, y in -100.0f64..100.0) {
            let pga_point = Point::new(x, y);
            let na_point: na::Point2<f64> = pga_point.try_into().unwrap();
            let back: Point<f64> = na_point.into();

            prop_assert!(abs_diff_eq!(pga_point.x(), back.x(), epsilon = ABS_DIFF_EQ_EPS));
            prop_assert!(abs_diff_eq!(pga_point.y(), back.y(), epsilon = ABS_DIFF_EQ_EPS));
        }

        #[test]
        fn motor_matches_isometry2(
            angle in -std::f64::consts::PI..std::f64::consts::PI,
            tx in -10.0f64..10.0, ty in -10.0f64..10.0,
            px in -10.0f64..10.0, py in -10.0f64..10.0,
        ) {
            // PGA motor
            let motor = Motor::from_translation(tx, ty).compose(&Motor::from_rotation(angle));
            let pga_point = Point::new(px, py);
            let pga_result = motor.transform_point(&pga_point);

            // nalgebra Isometry2
            let iso = na::Isometry2::new(na::Vector2::new(tx, ty), angle);
            let na_point = na::Point2::new(px, py);
            let na_result = iso * na_point;

            prop_assert!(abs_diff_eq!(pga_result.x(), na_result.x, epsilon = ABS_DIFF_EQ_EPS));
            prop_assert!(abs_diff_eq!(pga_result.y(), na_result.y, epsilon = ABS_DIFF_EQ_EPS));
        }
    }
}
```

### Files to Create (Phase 2)

- `src/specialized/projective/mod.rs`
- `src/specialized/projective/dim2/mod.rs`
- `src/specialized/projective/dim2/types.rs`
- `src/specialized/projective/dim2/ops.rs`
- `src/specialized/projective/dim2/conversions.rs`
- `src/specialized/projective/dim2/arbitrary.rs`
- `src/specialized/projective/dim2/nalgebra.rs`
- Update `src/specialized/mod.rs`

### Phase 2 Verification

- [x] `cargo fmt` passes
- [x] `cargo clippy --all-features` passes
- [x] `cargo doc --all-features --no-deps` passes
- [x] `cargo test --all-features` - all property tests pass
- [x] `cargo deny check` passes
- [x] Operations match generic Multivector
- [x] Operations match nalgebra equivalents

---

## Phase 3: Specialized 3D Types

### Module Structure

```
src/specialized/projective/dim3/
├── mod.rs          # Public API
├── types.rs        # Point, Line, Plane, Motor
├── ops.rs          # Operations
├── conversions.rs  # To/from Multivector
├── arbitrary.rs    # Proptest support
└── nalgebra.rs     # nalgebra conversions (feature-gated)
```

### Type Definitions (`specialized/projective/dim3/types.rs`)

```rust
/// A point in 3D PGA (grade 1 vector).
///
/// In point-based PGA, a point is represented as a homogeneous coordinate:
/// `P = x·e₁ + y·e₂ + z·e₃ + w·e₀`
///
/// where `(x/w, y/w, z/w)` are the Cartesian coordinates when `w ≠ 0`.
/// Points with `w = 0` represent ideal points (points at infinity).
#[derive(Clone, Copy, Debug, PartialEq)]
#[repr(C)]
pub struct Point<T: Float> {
    /// Coefficient of e₁
    pub e1: T,
    /// Coefficient of e₂
    pub e2: T,
    /// Coefficient of e₃
    pub e3: T,
    /// Coefficient of e₀ (homogeneous weight)
    pub e0: T,
}

impl<T: Float> Point<T> {
    /// Creates a finite point at Cartesian coordinates (x, y, z).
    #[inline]
    pub fn new(x: T, y: T, z: T) -> Self {
        Self { e1: x, e2: y, e3: z, e0: T::one() }
    }

    /// Creates a point from homogeneous coordinates.
    #[inline]
    pub fn from_homogeneous(e1: T, e2: T, e3: T, e0: T) -> Self {
        Self { e1, e2, e3, e0 }
    }

    /// Creates an ideal point (point at infinity) in the given direction.
    #[inline]
    pub fn ideal(dx: T, dy: T, dz: T) -> Self {
        Self { e1: dx, e2: dy, e3: dz, e0: T::zero() }
    }

    /// Origin point (0, 0, 0).
    #[inline]
    pub fn origin() -> Self {
        Self::new(T::zero(), T::zero(), T::zero())
    }

    /// Returns true if this is an ideal point.
    #[inline]
    pub fn is_ideal(&self, epsilon: T) -> bool {
        self.e0.abs() < epsilon
    }
}

/// A line in 3D PGA (grade 2 bivector).
///
/// A line in 3D PGA has 6 components representing a Plücker line.
#[derive(Clone, Copy, Debug, PartialEq)]
#[repr(C)]
pub struct Line<T: Float> {
    /// Coefficient of e₁₂
    pub e12: T,
    /// Coefficient of e₁₃
    pub e13: T,
    /// Coefficient of e₂₃
    pub e23: T,
    /// Coefficient of e₀₁
    pub e01: T,
    /// Coefficient of e₀₂
    pub e02: T,
    /// Coefficient of e₀₃
    pub e03: T,
}

impl<T: Float> Line<T> {
    /// Creates a line through two points.
    pub fn through_points(p1: &Point<T>, p2: &Point<T>) -> Self;

    /// Returns the direction vector of the line.
    pub fn direction(&self) -> (T, T, T);

    /// Returns a point on the line.
    pub fn point_on_line(&self) -> Option<Point<T>>;
}

/// A plane in 3D PGA (grade 3 trivector).
///
/// In point-based PGA, a plane is represented as a trivector.
#[derive(Clone, Copy, Debug, PartialEq)]
#[repr(C)]
pub struct Plane<T: Float> {
    /// Coefficient of e₁₂₃
    pub e123: T,
    /// Coefficient of e₀₁₂
    pub e012: T,
    /// Coefficient of e₀₁₃
    pub e013: T,
    /// Coefficient of e₀₂₃
    pub e023: T,
}

impl<T: Float> Plane<T> {
    /// Creates a plane from implicit equation ax + by + cz + d = 0.
    pub fn from_implicit(a: T, b: T, c: T, d: T) -> Self;

    /// Creates the xy-plane (z = 0).
    #[inline]
    pub fn xy() -> Self {
        Self::from_implicit(T::zero(), T::zero(), T::one(), T::zero())
    }

    /// Creates the xz-plane (y = 0).
    #[inline]
    pub fn xz() -> Self {
        Self::from_implicit(T::zero(), T::one(), T::zero(), T::zero())
    }

    /// Creates the yz-plane (x = 0).
    #[inline]
    pub fn yz() -> Self {
        Self::from_implicit(T::one(), T::zero(), T::zero(), T::zero())
    }

    /// Returns the normal direction (a, b, c).
    pub fn normal(&self) -> (T, T, T);

    /// Signed distance from origin.
    pub fn distance_from_origin(&self) -> T;
}

/// A motor (rigid transformation) in 3D PGA.
///
/// Motors represent rigid body transformations (rotation + translation).
/// A motor is an even-grade element with 8 components.
#[derive(Clone, Copy, Debug, PartialEq)]
#[repr(C)]
pub struct Motor<T: Float> {
    /// Scalar part
    pub s: T,
    /// Bivector parts (6 components)
    pub e12: T,
    pub e13: T,
    pub e23: T,
    pub e01: T,
    pub e02: T,
    pub e03: T,
    /// Pseudoscalar part (e₀₁₂₃)
    pub e0123: T,
}

impl<T: Float> Motor<T> {
    /// Creates the identity motor.
    pub fn identity() -> Self;

    /// Creates a pure rotation around an axis through the origin.
    pub fn from_axis_angle(axis: (T, T, T), angle: T) -> Self;

    /// Creates a pure translation.
    pub fn from_translation(dx: T, dy: T, dz: T) -> Self;

    /// Creates a rotation around an arbitrary line.
    pub fn from_line_angle(line: &Line<T>, angle: T) -> Self;

    /// Transforms a point.
    pub fn transform_point(&self, p: &Point<T>) -> Point<T>;

    /// Transforms a line.
    pub fn transform_line(&self, l: &Line<T>) -> Line<T>;

    /// Transforms a plane.
    pub fn transform_plane(&self, p: &Plane<T>) -> Plane<T>;

    /// Composes two motors.
    pub fn compose(&self, other: &Self) -> Self;

    /// Returns the inverse motor.
    pub fn inverse(&self) -> Self;

    /// Screw linear interpolation.
    pub fn slerp(&self, other: &Self, t: T) -> Self;
}
```

### Operations (`specialized/projective/dim3/ops.rs`)

```rust
impl<T: Float> Point<T> {
    /// Join of two points: line through them.
    pub fn join(&self, other: &Point<T>) -> Line<T>;

    /// Join of three points: plane through them.
    pub fn join_plane(&self, p2: &Point<T>, p3: &Point<T>) -> Plane<T>;

    /// Distance to another point.
    pub fn distance(&self, other: &Point<T>) -> T;
}

impl<T: Float> Line<T> {
    /// Meet of two lines: point of intersection (if they intersect).
    pub fn meet(&self, other: &Line<T>) -> Option<Point<T>>;

    /// Join of line and point: plane through them.
    pub fn join(&self, p: &Point<T>) -> Plane<T>;

    /// Distance from a point to the line.
    pub fn distance_to_point(&self, p: &Point<T>) -> T;

    /// Closest point on line to given point.
    pub fn closest_point(&self, p: &Point<T>) -> Point<T>;
}

impl<T: Float> Plane<T> {
    /// Meet of two planes: line of intersection.
    pub fn meet(&self, other: &Plane<T>) -> Line<T>;

    /// Meet of three planes: point of intersection.
    pub fn meet_point(&self, p2: &Plane<T>, p3: &Plane<T>) -> Point<T>;

    /// Meet of plane and line: point of intersection.
    pub fn meet_line(&self, l: &Line<T>) -> Point<T>;

    /// Signed distance from a point to the plane.
    pub fn distance_to_point(&self, p: &Point<T>) -> T;

    /// Project a point onto the plane.
    pub fn project(&self, p: &Point<T>) -> Point<T>;
}
```

### Phase 3 Testing

Similar to Phase 2, with additional tests for:

```rust
proptest! {
    // Plane operations
    #[test]
    fn three_points_define_plane(
        p1 in any::<Point<f64>>(),
        p2 in any::<Point<f64>>(),
        p3 in any::<Point<f64>>(),
    ) {
        let plane = p1.join_plane(&p2, &p3);
        // All three points should lie on the plane
        prop_assert!(plane.distance_to_point(&p1).abs() < ABS_DIFF_EQ_EPS);
        prop_assert!(plane.distance_to_point(&p2).abs() < ABS_DIFF_EQ_EPS);
        prop_assert!(plane.distance_to_point(&p3).abs() < ABS_DIFF_EQ_EPS);
    }

    #[test]
    fn plane_meet_gives_line(
        plane1 in any::<Plane<f64>>(),
        plane2 in any::<Plane<f64>>(),
    ) {
        let line = plane1.meet(&plane2);
        // Line should lie on both planes
        if let Some(point) = line.point_on_line() {
            prop_assert!(plane1.distance_to_point(&point).abs() < ABS_DIFF_EQ_EPS);
            prop_assert!(plane2.distance_to_point(&point).abs() < ABS_DIFF_EQ_EPS);
        }
    }

    // Motor with Isometry3
    #[test]
    fn motor_matches_isometry3(
        axis_x in -1.0f64..1.0, axis_y in -1.0f64..1.0, axis_z in -1.0f64..1.0,
        angle in -std::f64::consts::PI..std::f64::consts::PI,
        tx in -10.0f64..10.0, ty in -10.0f64..10.0, tz in -10.0f64..10.0,
        px in -10.0f64..10.0, py in -10.0f64..10.0, pz in -10.0f64..10.0,
    ) {
        // Normalize axis
        let len = (axis_x*axis_x + axis_y*axis_y + axis_z*axis_z).sqrt();
        if len < ABS_DIFF_EQ_EPS {
            return Ok(());
        }
        let axis = (axis_x/len, axis_y/len, axis_z/len);

        // PGA motor
        let rotation = Motor::from_axis_angle(axis, angle);
        let translation = Motor::from_translation(tx, ty, tz);
        let motor = translation.compose(&rotation);
        let pga_point = Point::new(px, py, pz);
        let pga_result = motor.transform_point(&pga_point);

        // nalgebra Isometry3
        let na_axis = na::Unit::new_normalize(na::Vector3::new(axis.0, axis.1, axis.2));
        let rotation = na::UnitQuaternion::from_axis_angle(&na_axis, angle);
        let translation = na::Translation3::new(tx, ty, tz);
        let iso = na::Isometry3::from_parts(translation, rotation);
        let na_point = na::Point3::new(px, py, pz);
        let na_result = iso * na_point;

        prop_assert!(abs_diff_eq!(pga_result.x(), na_result.x, epsilon = ABS_DIFF_EQ_EPS));
        prop_assert!(abs_diff_eq!(pga_result.y(), na_result.y, epsilon = ABS_DIFF_EQ_EPS));
        prop_assert!(abs_diff_eq!(pga_result.z(), na_result.z, epsilon = ABS_DIFF_EQ_EPS));
    }
}
```

### Files to Create (Phase 3)

- `src/specialized/projective/dim3/mod.rs`
- `src/specialized/projective/dim3/types.rs`
- `src/specialized/projective/dim3/ops.rs`
- `src/specialized/projective/dim3/conversions.rs`
- `src/specialized/projective/dim3/arbitrary.rs`
- `src/specialized/projective/dim3/nalgebra.rs`

### Phase 3 Verification

- [x] `cargo fmt` passes
- [x] `cargo clippy --all-features` passes
- [x] `cargo doc --all-features --no-deps` passes
- [x] `cargo test --all-features` - all property tests pass
- [x] `cargo deny check` passes
- [x] Operations match generic Multivector
- [x] Operations match nalgebra equivalents

---

## nalgebra Interoperability

### Conversions (feature-gated)

```rust
// Point conversions
impl<T: Float + na::Scalar> From<na::Point2<T>> for Point<T>;
impl<T: Float + na::Scalar> TryFrom<Point<T>> for na::Point2<T>;

impl<T: Float + na::Scalar> From<na::Point3<T>> for Point<T>;
impl<T: Float + na::Scalar> TryFrom<Point<T>> for na::Point3<T>;

// Motor conversions
impl<T: Float + na::RealField> From<na::Isometry2<T>> for Motor<T>;
impl<T: Float + na::RealField> From<Motor<T>> for na::Isometry2<T>;

impl<T: Float + na::RealField> From<na::Isometry3<T>> for Motor<T>;
impl<T: Float + na::RealField> From<Motor<T>> for na::Isometry3<T>;
```

### Error Types

```rust
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum NalgebraConversionError {
    /// Point has w ≈ 0 (point at infinity)
    PointAtInfinity,
    /// Line is degenerate (zero)
    DegenerateLine,
    /// Motor is not normalized
    UnnormalizedMotor,
}
```

---

## Benchmarks

### Files to Create

- `benches/projective.rs` - PGA operation benchmarks
- `benches/reports/projective/dim2/` - 2D benchmark SVGs
- `benches/reports/projective/dim3/` - 3D benchmark SVGs

### Key Benchmarks

```rust
// Motor operations
fn bench_motor3_transform_point(c: &mut Criterion);
fn bench_motor3_compose(c: &mut Criterion);

// Meet/Join operations
fn bench_point_join_line(c: &mut Criterion);
fn bench_plane_meet_line(c: &mut Criterion);

// Comparison with nalgebra
fn bench_motor_vs_isometry(c: &mut Criterion);
```

### Expected Performance Targets

| Operation | Target | Notes |
|-----------|--------|-------|
| Motor transform point | < 15 ns | Sandwich product optimization |
| Motor compose | < 20 ns | Even subalgebra product |
| Point join (line) | < 5 ns | Regressive product |
| Plane meet (line) | < 5 ns | Outer product |
| Motor vs Isometry3 | Comparable | Similar computational cost |

---

## Summary

| Phase | Scope | Key Deliverables | Status |
|-------|-------|------------------|--------|
| 1 | Generic PGA | Signatures, property tests | ✅ Complete |
| 2 | 2D Specialized | Point, Line, Motor, nalgebra tests | ✅ Complete |
| 3 | 3D Specialized | Point, Line, Plane, Motor, nalgebra tests | ✅ Complete |
| 4 | Geometric Constraints | Study condition, Plücker condition, invariants | Pending |

Each phase must pass all verification checks before proceeding to the next. Commits happen at the end of each phase after tests pass.

---

## Phase 4: Enforce Geometric Constraints as Invariants

**Status**: Pending

**Reference**: https://rigidgeometricalgebra.org/wiki/index.php?title=Geometric_constraint

### Background

PGA elements must satisfy geometric constraints to represent valid geometric objects. These constraints should be enforced as type-level invariants where possible, and verified via property tests.

### Key Constraints

| Element | Constraint | Meaning |
|---------|------------|---------|
| **Motor** | `s·e0123 - e23·e01 - e31·e02 - e12·e03 = 0` | Study condition: proper rigid transformation |
| **Unit Motor** | `s² + e23² + e31² + e12² = 1` | Normalized rotation part |
| **Line** | `e23·e01 + e31·e02 + e12·e03 = 0` | Plücker condition: valid 3D line |
| **Unit Line** | `e23² + e31² + e12² = 1` | Normalized direction |
| **Flector** | (depends on form) | Valid reflection/glide |

### Implementation Approach

1. **Wrapper Types for Constrained Elements**
   - `UnitMotor<T>` - Motor satisfying study condition + unit norm
   - `UnitLine<T>` - Line satisfying Plücker condition + unit direction
   - Constructors verify constraints, unsafe constructors for trusted input

2. **Property Tests for Invariants**
   - Verify constraints hold after construction
   - Verify constraints preserved by operations (composition, inverse)
   - Test constraint satisfaction in `Arbitrary` implementations

3. **Constraint Verification Methods**
   ```rust
   impl<T: Float> Motor<T> {
       /// Returns true if this motor satisfies the study condition.
       pub fn satisfies_study_condition(&self, epsilon: T) -> bool;

       /// Returns true if the rotation part is normalized.
       pub fn is_normalized(&self, epsilon: T) -> bool;
   }

   impl<T: Float> Line<T> {
       /// Returns true if this line satisfies the Plücker condition.
       pub fn satisfies_plucker_condition(&self, epsilon: T) -> bool;
   }
   ```

### Deliverables

- [ ] Add `satisfies_study_condition()` to Motor
- [ ] Add `satisfies_plucker_condition()` to Line
- [ ] Add property tests verifying constraints preserved by operations
- [ ] Document constraints in rustdoc with references to wiki
- [ ] Consider `UnitMotor` and `UnitLine` wrapper types (if ergonomic)

### Phase 4 Verification

- [ ] `cargo fmt` passes
- [ ] `cargo clippy --all-features` passes
- [ ] `cargo doc --all-features --no-deps` passes
- [ ] `cargo test --all-features` - constraint property tests pass
- [ ] `cargo deny check` passes
