# PRD-6.5: Translator and Rotor Versors

**Status**: Pending
**Parent**: PRD-6 (Conformal Geometric Algebra)
**Depends On**: PRD-6.4 (Flat Points)
**Goal**: Implement translator and rotor versors for rigid body transformations

## Reference

- https://conformalgeometricalgebra.org/wiki/index.php?title=Translation
- https://conformalgeometricalgebra.org/wiki/index.php?title=Rotation

## Background

### Versors in CGA

Versors are the geometric product of vectors that represent transformations. In CGA:
- **Translator**: Translation by a displacement vector
- **Rotor**: Rotation around a line through the origin
- **Motor**: Combined rotation and translation (translator * rotor)

All versors transform geometric objects via the **sandwich product**: `T ⋈ X ⋈ T̃`

### Translator

The translator for displacement τ = (τₓ, τᵧ, τᵤ) is:

```
T = 1 - (τ/2)·e∞
  = 1 - (τₓ/2)·e₁∞ - (τᵧ/2)·e₂∞ - (τᵤ/2)·e₃∞
```

where e∞ = e₋ + e₊ is the point at infinity.

### Rotor

The rotor for rotation by angle θ around an axis through the origin:

```
R = cos(θ/2) + sin(θ/2)·B
```

where B is the bivector representing the plane of rotation.

### Motor (Translator × Rotor)

A motor combines rotation and translation:
```
M = T · R
```

Composition follows: `M₁ · M₂` applies M₂ first, then M₁.

## Deliverables

### 1. Translator Type (`dim3/translator.rs`)

```rust
use crate::scalar::Float;
use super::{Point, Sphere, Plane, Circle, Line};

/// A translator in 3D CGA.
///
/// Translates geometric objects by a displacement vector τ = (τₓ, τᵧ, τᵤ).
///
/// # Formula
///
/// ```text
/// T = 1 - (τ/2)·e∞ = 1 - (τₓ/2)·e₁∞ - (τᵧ/2)·e₂∞ - (τᵤ/2)·e₃∞
/// ```
///
/// # Example
///
/// ```
/// use clifford::specialized::conformal::dim3::{Point, Translator};
///
/// // Translate by (1, 2, 3)
/// let t = Translator::new(1.0, 2.0, 3.0);
///
/// // Transform a point
/// let p = Point::new(0.0, 0.0, 0.0);
/// let q = t.transform_point(&p);
///
/// assert!((q.x() - 1.0).abs() < 1e-10);
/// assert!((q.y() - 2.0).abs() < 1e-10);
/// assert!((q.z() - 3.0).abs() < 1e-10);
/// ```
///
/// Reference: https://conformalgeometricalgebra.org/wiki/index.php?title=Translation
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Translator<T: Float> {
    /// Scalar part (always 1 for unit translator)
    s: T,
    /// Coefficient of e₁∞ (half of x displacement)
    e1i: T,
    /// Coefficient of e₂∞ (half of y displacement)
    e2i: T,
    /// Coefficient of e₃∞ (half of z displacement)
    e3i: T,
}

impl<T: Float> Translator<T> {
    /// Creates a translator for displacement (dx, dy, dz).
    #[inline]
    pub fn new(dx: T, dy: T, dz: T) -> Self {
        let half = T::one() / T::TWO;
        Self {
            s: T::one(),
            e1i: -dx * half,
            e2i: -dy * half,
            e3i: -dz * half,
        }
    }

    /// Identity translator (no displacement).
    #[inline]
    pub fn identity() -> Self {
        Self {
            s: T::one(),
            e1i: T::zero(),
            e2i: T::zero(),
            e3i: T::zero(),
        }
    }

    /// Returns the displacement vector (dx, dy, dz).
    #[inline]
    pub fn displacement(&self) -> (T, T, T) {
        let two = T::TWO;
        (-self.e1i * two, -self.e2i * two, -self.e3i * two)
    }

    /// Returns the reverse (inverse) translator.
    ///
    /// For a translator, T̃ = T⁻¹ = 1 + (τ/2)·e∞
    #[inline]
    pub fn reverse(&self) -> Self {
        Self {
            s: self.s,
            e1i: -self.e1i,
            e2i: -self.e2i,
            e3i: -self.e3i,
        }
    }

    /// Composes two translators: self applied first, then other.
    ///
    /// Translation is commutative: T₁ · T₂ = T₂ · T₁
    pub fn compose(&self, other: &Self) -> Self {
        // Translations add
        let (dx1, dy1, dz1) = self.displacement();
        let (dx2, dy2, dz2) = other.displacement();
        Self::new(dx1 + dx2, dy1 + dy2, dz1 + dz2)
    }

    /// Transforms a round point: T ⋈ P ⋈ T̃
    pub fn transform_point(&self, p: &Point<T>) -> Point<T> {
        // Derived from SymPy
        let (dx, dy, dz) = self.displacement();
        Point::new(p.x() + dx, p.y() + dy, p.z() + dz)
    }

    /// Transforms a sphere.
    pub fn transform_sphere(&self, s: &Sphere<T>) -> Sphere<T> {
        // Sphere center translates, radius unchanged
        todo!("Derive from SymPy: derive_translator_transform_sphere()")
    }

    /// Transforms a plane.
    pub fn transform_plane(&self, p: &Plane<T>) -> Plane<T> {
        todo!("Derive from SymPy: derive_translator_transform_plane()")
    }

    /// Transforms a circle.
    pub fn transform_circle(&self, c: &Circle<T>) -> Circle<T> {
        todo!("Derive from SymPy: derive_translator_transform_circle()")
    }

    /// Transforms a line.
    pub fn transform_line(&self, l: &Line<T>) -> Line<T> {
        todo!("Derive from SymPy: derive_translator_transform_line()")
    }

    /// Component accessors.
    #[inline]
    pub fn s(&self) -> T { self.s }
    #[inline]
    pub fn e1i(&self) -> T { self.e1i }
    #[inline]
    pub fn e2i(&self) -> T { self.e2i }
    #[inline]
    pub fn e3i(&self) -> T { self.e3i }
}
```

### 2. Rotor Type (`dim3/rotor.rs`)

```rust
use crate::scalar::Float;
use super::{Point, Sphere, Plane, Circle, Line};

/// A rotor in 3D CGA (rotation around an axis through the origin).
///
/// # Formula
///
/// ```text
/// R = cos(θ/2) + sin(θ/2)·B
/// ```
///
/// where B is the bivector representing the plane of rotation.
///
/// # Example
///
/// ```
/// use clifford::specialized::conformal::dim3::{Point, Rotor};
/// use std::f64::consts::FRAC_PI_2;
///
/// // 90° rotation around z-axis
/// let r = Rotor::from_axis_angle((0.0, 0.0, 1.0), FRAC_PI_2);
///
/// // Rotate point (1, 0, 0) -> (0, 1, 0)
/// let p = Point::new(1.0, 0.0, 0.0);
/// let q = r.transform_point(&p);
///
/// assert!(q.x().abs() < 1e-10);
/// assert!((q.y() - 1.0).abs() < 1e-10);
/// ```
///
/// Reference: https://conformalgeometricalgebra.org/wiki/index.php?title=Rotation
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Rotor<T: Float> {
    /// Scalar part: cos(θ/2)
    s: T,
    /// Coefficient of e₁₂: sin(θ/2) · axis_z
    e12: T,
    /// Coefficient of e₁₃: -sin(θ/2) · axis_y
    e13: T,
    /// Coefficient of e₂₃: sin(θ/2) · axis_x
    e23: T,
}

impl<T: Float> Rotor<T> {
    /// Creates a rotor from axis (unit vector) and angle.
    pub fn from_axis_angle(axis: (T, T, T), angle: T) -> Self {
        let half = angle / T::TWO;
        let c = half.cos();
        let s = half.sin();
        Self {
            s: c,
            e12: s * axis.2,   // z-component
            e13: -s * axis.1,  // -y-component
            e23: s * axis.0,   // x-component
        }
    }

    /// Creates a rotor from Euler angles (ZYX convention).
    pub fn from_euler_angles(roll: T, pitch: T, yaw: T) -> Self {
        let rz = Self::from_axis_angle((T::zero(), T::zero(), T::one()), yaw);
        let ry = Self::from_axis_angle((T::zero(), T::one(), T::zero()), pitch);
        let rx = Self::from_axis_angle((T::one(), T::zero(), T::zero()), roll);
        rz.compose(&ry).compose(&rx)
    }

    /// Identity rotor (no rotation).
    #[inline]
    pub fn identity() -> Self {
        Self {
            s: T::one(),
            e12: T::zero(),
            e13: T::zero(),
            e23: T::zero(),
        }
    }

    /// Rotation around x-axis by angle.
    pub fn rotate_x(angle: T) -> Self {
        Self::from_axis_angle((T::one(), T::zero(), T::zero()), angle)
    }

    /// Rotation around y-axis by angle.
    pub fn rotate_y(angle: T) -> Self {
        Self::from_axis_angle((T::zero(), T::one(), T::zero()), angle)
    }

    /// Rotation around z-axis by angle.
    pub fn rotate_z(angle: T) -> Self {
        Self::from_axis_angle((T::zero(), T::zero(), T::one()), angle)
    }

    /// Returns the rotation angle.
    pub fn angle(&self) -> T {
        T::TWO * self.s.acos()
    }

    /// Returns the rotation axis (unit vector).
    pub fn axis(&self) -> Option<(T, T, T)> {
        let sin_half = (T::one() - self.s * self.s).sqrt();
        if sin_half.abs() < T::epsilon() {
            None // Identity rotation, axis undefined
        } else {
            Some((
                self.e23 / sin_half,
                -self.e13 / sin_half,
                self.e12 / sin_half,
            ))
        }
    }

    /// Returns the reverse (inverse) rotor.
    #[inline]
    pub fn reverse(&self) -> Self {
        Self {
            s: self.s,
            e12: -self.e12,
            e13: -self.e13,
            e23: -self.e23,
        }
    }

    /// Composes two rotors: self applied first, then other.
    pub fn compose(&self, other: &Self) -> Self {
        // Quaternion-like composition
        todo!("Derive from SymPy: derive_rotor_composition()")
    }

    /// Transforms a round point: R ⋈ P ⋈ R̃
    pub fn transform_point(&self, p: &Point<T>) -> Point<T> {
        todo!("Derive from SymPy: derive_rotor_transform_point()")
    }

    /// Transforms a sphere.
    pub fn transform_sphere(&self, s: &Sphere<T>) -> Sphere<T> {
        todo!("Derive from SymPy")
    }

    /// Transforms a plane.
    pub fn transform_plane(&self, p: &Plane<T>) -> Plane<T> {
        todo!("Derive from SymPy")
    }

    /// Transforms a circle.
    pub fn transform_circle(&self, c: &Circle<T>) -> Circle<T> {
        todo!("Derive from SymPy")
    }

    /// Transforms a line.
    pub fn transform_line(&self, l: &Line<T>) -> Line<T> {
        todo!("Derive from SymPy")
    }

    /// Spherical linear interpolation between rotors.
    pub fn slerp(&self, other: &Self, t: T) -> Self {
        todo!("Implement slerp")
    }

    /// Component accessors.
    #[inline]
    pub fn s(&self) -> T { self.s }
    #[inline]
    pub fn e12(&self) -> T { self.e12 }
    #[inline]
    pub fn e13(&self) -> T { self.e13 }
    #[inline]
    pub fn e23(&self) -> T { self.e23 }
}
```

### 3. Motor Type (`dim3/motor.rs`)

```rust
use crate::scalar::Float;
use super::{Point, Sphere, Plane, Circle, Line, Translator, Rotor};

/// A motor in 3D CGA (combined rotation and translation).
///
/// M = T · R where T is a translator and R is a rotor.
///
/// # Example
///
/// ```
/// use clifford::specialized::conformal::dim3::{Point, Motor, Translator, Rotor};
/// use std::f64::consts::FRAC_PI_2;
///
/// // Rotate 90° around z, then translate by (1, 0, 0)
/// let r = Rotor::rotate_z(FRAC_PI_2);
/// let t = Translator::new(1.0, 0.0, 0.0);
/// let m = Motor::from_translator_rotor(&t, &r);
///
/// // (1, 0, 0) -> rotate to (0, 1, 0) -> translate to (1, 1, 0)
/// let p = Point::new(1.0, 0.0, 0.0);
/// let q = m.transform_point(&p);
/// ```
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Motor<T: Float> {
    /// Scalar part
    s: T,
    /// Grade-2 bivector parts (6 components)
    e12: T,
    e13: T,
    e23: T,
    e1i: T,  // e₁∞ - translation in x
    e2i: T,  // e₂∞ - translation in y
    e3i: T,  // e₃∞ - translation in z
    /// Pseudoscalar part (e₁₂₃∞)
    e123i: T,
}

impl<T: Float> Motor<T> {
    /// Creates a motor from translator and rotor.
    ///
    /// M = T · R (translator applied after rotor)
    pub fn from_translator_rotor(t: &Translator<T>, r: &Rotor<T>) -> Self {
        todo!("Derive from SymPy: derive_motor_from_tr()")
    }

    /// Creates a motor from rotor and translator.
    ///
    /// M = R · T (rotor applied after translator)
    pub fn from_rotor_translator(r: &Rotor<T>, t: &Translator<T>) -> Self {
        todo!("Derive from SymPy")
    }

    /// Identity motor (no transformation).
    pub fn identity() -> Self {
        Self {
            s: T::one(),
            e12: T::zero(),
            e13: T::zero(),
            e23: T::zero(),
            e1i: T::zero(),
            e2i: T::zero(),
            e3i: T::zero(),
            e123i: T::zero(),
        }
    }

    /// Creates a pure translation motor.
    pub fn from_translation(dx: T, dy: T, dz: T) -> Self {
        Self::from_translator_rotor(&Translator::new(dx, dy, dz), &Rotor::identity())
    }

    /// Creates a pure rotation motor.
    pub fn from_rotation(axis: (T, T, T), angle: T) -> Self {
        Self::from_translator_rotor(&Translator::identity(), &Rotor::from_axis_angle(axis, angle))
    }

    /// Decomposes into translator and rotor components.
    pub fn decompose(&self) -> (Translator<T>, Rotor<T>) {
        todo!("Derive from SymPy")
    }

    /// Returns the reverse (inverse) motor.
    pub fn reverse(&self) -> Self {
        todo!("Derive from SymPy")
    }

    /// Composes two motors: self applied first, then other.
    pub fn compose(&self, other: &Self) -> Self {
        todo!("Derive from SymPy: derive_motor_composition()")
    }

    /// Transforms a point.
    pub fn transform_point(&self, p: &Point<T>) -> Point<T> {
        todo!("Derive from SymPy")
    }

    /// Transforms a sphere.
    pub fn transform_sphere(&self, s: &Sphere<T>) -> Sphere<T> {
        todo!("Derive from SymPy")
    }

    /// Transforms a plane.
    pub fn transform_plane(&self, p: &Plane<T>) -> Plane<T> {
        todo!("Derive from SymPy")
    }

    /// Transforms a circle.
    pub fn transform_circle(&self, c: &Circle<T>) -> Circle<T> {
        todo!("Derive from SymPy")
    }

    /// Transforms a line.
    pub fn transform_line(&self, l: &Line<T>) -> Line<T> {
        todo!("Derive from SymPy")
    }
}
```

## SymPy Derivations

Add to `derivations/src/clifford_derivations/cga.py`:

```python
@with_timeout(300)
def derive_translator_transform_point():
    """Derive T ⋈ P ⋈ T̃ for translator."""
    dx, dy, dz = symbols('dx dy dz')
    px, py, pz = symbols('px py pz')

    # Translator: T = 1 - (τ/2)·e∞
    # Point: P = px·e₁ + py·e₂ + pz·e₃ + e₀ + ½(px² + py² + pz²)·e∞

    # Result should be: P' at (px + dx, py + dy, pz + dz)

    print("Translator transform point:")
    print("  Result: (px + dx, py + dy, pz + dz)")
    print("  Translation simply adds to coordinates")


@with_timeout(300)
def derive_rotor_composition():
    """Derive R₁ · R₂ rotor composition."""
    # Similar to quaternion multiplication
    s1, e12_1, e13_1, e23_1 = symbols('s1 e12_1 e13_1 e23_1')
    s2, e12_2, e13_2, e23_2 = symbols('s2 e12_2 e13_2 e23_2')

    # TODO: Compute full composition
    print("TODO: Implement rotor composition derivation")


@with_timeout(300)
def derive_rotor_transform_point():
    """Derive R ⋈ P ⋈ R̃ for rotor transforming a point."""
    # This is the 3D rotation formula
    print("TODO: Implement rotor transform point derivation")


@with_timeout(600)
def derive_motor_composition():
    """Derive M₁ · M₂ motor composition."""
    # Full 8-component motor product
    print("TODO: Implement motor composition derivation")
```

## Property Tests

```rust
proptest! {
    // ================================================================
    // Translator tests
    // ================================================================

    #[test]
    fn translator_transform_matches_nalgebra(
        dx in -10.0f64..10.0, dy in -10.0f64..10.0, dz in -10.0f64..10.0,
        px in -10.0f64..10.0, py in -10.0f64..10.0, pz in -10.0f64..10.0,
    ) {
        let t = Translator::new(dx, dy, dz);
        let p = Point::new(px, py, pz);
        let q = t.transform_point(&p);

        let na_t = na::Translation3::new(dx, dy, dz);
        let na_p = na::Point3::new(px, py, pz);
        let na_q = na_t * na_p;

        prop_assert!(abs_diff_eq!(q.x(), na_q.x, epsilon = ABS_DIFF_EQ_EPS));
        prop_assert!(abs_diff_eq!(q.y(), na_q.y, epsilon = ABS_DIFF_EQ_EPS));
        prop_assert!(abs_diff_eq!(q.z(), na_q.z, epsilon = ABS_DIFF_EQ_EPS));
    }

    #[test]
    fn translator_composition_is_additive(
        dx1 in -10.0f64..10.0, dy1 in -10.0f64..10.0, dz1 in -10.0f64..10.0,
        dx2 in -10.0f64..10.0, dy2 in -10.0f64..10.0, dz2 in -10.0f64..10.0,
    ) {
        let t1 = Translator::new(dx1, dy1, dz1);
        let t2 = Translator::new(dx2, dy2, dz2);
        let t12 = t1.compose(&t2);

        let (dx, dy, dz) = t12.displacement();
        prop_assert!(abs_diff_eq!(dx, dx1 + dx2, epsilon = ABS_DIFF_EQ_EPS));
        prop_assert!(abs_diff_eq!(dy, dy1 + dy2, epsilon = ABS_DIFF_EQ_EPS));
        prop_assert!(abs_diff_eq!(dz, dz1 + dz2, epsilon = ABS_DIFF_EQ_EPS));
    }

    #[test]
    fn translator_inverse(
        dx in -10.0f64..10.0, dy in -10.0f64..10.0, dz in -10.0f64..10.0,
        px in -10.0f64..10.0, py in -10.0f64..10.0, pz in -10.0f64..10.0,
    ) {
        let t = Translator::new(dx, dy, dz);
        let p = Point::new(px, py, pz);

        let q = t.transform_point(&p);
        let back = t.reverse().transform_point(&q);

        prop_assert!(abs_diff_eq!(back.x(), p.x(), epsilon = ABS_DIFF_EQ_EPS));
        prop_assert!(abs_diff_eq!(back.y(), p.y(), epsilon = ABS_DIFF_EQ_EPS));
        prop_assert!(abs_diff_eq!(back.z(), p.z(), epsilon = ABS_DIFF_EQ_EPS));
    }

    // ================================================================
    // Rotor tests
    // ================================================================

    #[test]
    fn rotor_preserves_distance(
        axis_x in -1.0f64..1.0, axis_y in -1.0f64..1.0, axis_z in -1.0f64..1.0,
        angle in -std::f64::consts::PI..std::f64::consts::PI,
        p1 in any::<Point<f64>>(),
        p2 in any::<Point<f64>>(),
    ) {
        let len = (axis_x*axis_x + axis_y*axis_y + axis_z*axis_z).sqrt();
        if len < ABS_DIFF_EQ_EPS {
            return Ok(());
        }
        let axis = (axis_x/len, axis_y/len, axis_z/len);

        let r = Rotor::from_axis_angle(axis, angle);
        let q1 = r.transform_point(&p1);
        let q2 = r.transform_point(&p2);

        let d_before = p1.distance(&p2);
        let d_after = q1.distance(&q2);

        prop_assert!(abs_diff_eq!(d_before, d_after, epsilon = ABS_DIFF_EQ_EPS));
    }

    #[test]
    fn rotor_matches_nalgebra_quaternion(
        axis_x in -1.0f64..1.0, axis_y in -1.0f64..1.0, axis_z in -1.0f64..1.0,
        angle in -std::f64::consts::PI..std::f64::consts::PI,
        px in -10.0f64..10.0, py in -10.0f64..10.0, pz in -10.0f64..10.0,
    ) {
        let len = (axis_x*axis_x + axis_y*axis_y + axis_z*axis_z).sqrt();
        if len < ABS_DIFF_EQ_EPS {
            return Ok(());
        }
        let axis = (axis_x/len, axis_y/len, axis_z/len);

        let r = Rotor::from_axis_angle(axis, angle);
        let p = Point::new(px, py, pz);
        let q = r.transform_point(&p);

        let na_axis = na::Unit::new_normalize(na::Vector3::new(axis.0, axis.1, axis.2));
        let na_r = na::UnitQuaternion::from_axis_angle(&na_axis, angle);
        let na_p = na::Point3::new(px, py, pz);
        let na_q = na_r * na_p;

        prop_assert!(abs_diff_eq!(q.x(), na_q.x, epsilon = ABS_DIFF_EQ_EPS));
        prop_assert!(abs_diff_eq!(q.y(), na_q.y, epsilon = ABS_DIFF_EQ_EPS));
        prop_assert!(abs_diff_eq!(q.z(), na_q.z, epsilon = ABS_DIFF_EQ_EPS));
    }

    // ================================================================
    // Motor tests
    // ================================================================

    #[test]
    fn motor_matches_isometry3(
        axis_x in -1.0f64..1.0, axis_y in -1.0f64..1.0, axis_z in -1.0f64..1.0,
        angle in -std::f64::consts::PI..std::f64::consts::PI,
        dx in -10.0f64..10.0, dy in -10.0f64..10.0, dz in -10.0f64..10.0,
        px in -10.0f64..10.0, py in -10.0f64..10.0, pz in -10.0f64..10.0,
    ) {
        let len = (axis_x*axis_x + axis_y*axis_y + axis_z*axis_z).sqrt();
        if len < ABS_DIFF_EQ_EPS {
            return Ok(());
        }
        let axis = (axis_x/len, axis_y/len, axis_z/len);

        let r = Rotor::from_axis_angle(axis, angle);
        let t = Translator::new(dx, dy, dz);
        let m = Motor::from_translator_rotor(&t, &r);
        let p = Point::new(px, py, pz);
        let q = m.transform_point(&p);

        let na_axis = na::Unit::new_normalize(na::Vector3::new(axis.0, axis.1, axis.2));
        let na_r = na::UnitQuaternion::from_axis_angle(&na_axis, angle);
        let na_t = na::Translation3::new(dx, dy, dz);
        let na_iso = na::Isometry3::from_parts(na_t, na_r);
        let na_p = na::Point3::new(px, py, pz);
        let na_q = na_iso * na_p;

        prop_assert!(abs_diff_eq!(q.x(), na_q.x, epsilon = ABS_DIFF_EQ_EPS));
        prop_assert!(abs_diff_eq!(q.y(), na_q.y, epsilon = ABS_DIFF_EQ_EPS));
        prop_assert!(abs_diff_eq!(q.z(), na_q.z, epsilon = ABS_DIFF_EQ_EPS));
    }
}
```

## Files to Create/Modify

### New Files
- `src/specialized/conformal/dim3/translator.rs`
- `src/specialized/conformal/dim3/rotor.rs`
- `src/specialized/conformal/dim3/motor.rs`

### Modified Files
- `src/specialized/conformal/dim3/mod.rs` - Export new types
- `derivations/src/clifford_derivations/cga.py` - Add derivations

## Verification Checklist

- [ ] `cargo fmt` passes
- [ ] `cargo clippy --all-features` passes
- [ ] `cargo doc --all-features --no-deps` passes
- [ ] `cargo test --all-features` passes
- [ ] `cargo deny check` passes
- [ ] Translator matches nalgebra Translation3
- [ ] Rotor matches nalgebra UnitQuaternion
- [ ] Motor matches nalgebra Isometry3
- [ ] Distance preserved under rotations

## Dependencies

- PRD-6.4 (Flat Points) - must be complete

## Next Steps

After this PRD is complete, proceed to PRD-6.6 (Dilator and Inversor Versors).
