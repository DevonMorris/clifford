# PRD-6.6: Dilator and Inversor Versors

**Status**: Pending
**Parent**: PRD-6 (Conformal Geometric Algebra)
**Depends On**: PRD-6.5 (Translator and Rotor Versors)
**Goal**: Implement conformal-specific transformations: dilation and inversion

## Reference

- https://conformalgeometricalgebra.org/wiki/index.php?title=Dilation
- https://conformalgeometricalgebra.org/wiki/index.php?title=Inversion

## Background

### Dilator (Uniform Scaling)

The dilator performs **uniform scaling** by a factor σ about a center point. Unlike Euclidean scaling, this is a conformal transformation that maps circles to circles and spheres to spheres.

```
D = (1-σ)/2 · att(S) + (1+σ)/2 · 𝟙
```

where `att(S)` is the attitude of a sphere centered at the scaling center.

For scaling centered at the origin:
```
D = cosh(λ/2) + sinh(λ/2)·e₀∞
```
where σ = e^λ (scale factor as exponential).

### Inversor (Sphere Inversion)

The inversor performs **inversion in a sphere**, mapping points P to P' such that:
- P' lies on the ray from sphere center through P
- |CP| · |CP'| = r² (where C is center, r is radius)

This is a fundamental conformal transformation that:
- Maps circles/lines to circles/lines
- Maps spheres/planes to spheres/planes
- Is its own inverse (applying twice returns the original)

### Transversor (Special Conformal)

The transversor is a combination of inversion, translation, and inversion. It's the conformal analog of translation at infinity.

## Deliverables

### 1. Dilator Type (`dim3/dilator.rs`)

```rust
use crate::scalar::Float;
use super::{Point, Sphere, Plane, Circle, Line};

/// A dilator in 3D CGA (uniform scaling).
///
/// Performs uniform scaling by factor σ about a center point.
///
/// # Formula
///
/// For scaling centered at origin:
/// ```text
/// D = cosh(λ/2) + sinh(λ/2)·e₀∞
/// ```
/// where σ = e^λ is the scale factor.
///
/// # Example
///
/// ```
/// use clifford::specialized::conformal::dim3::{Point, Dilator};
///
/// // Scale by factor 2 centered at origin
/// let d = Dilator::from_scale(2.0);
///
/// // Point (1, 0, 0) -> (2, 0, 0)
/// let p = Point::new(1.0, 0.0, 0.0);
/// let q = d.transform_point(&p);
///
/// assert!((q.x() - 2.0).abs() < 1e-10);
/// ```
///
/// Reference: https://conformalgeometricalgebra.org/wiki/index.php?title=Dilation
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Dilator<T: Float> {
    /// cosh(λ/2) where σ = e^λ
    cosh_half: T,
    /// sinh(λ/2) where σ = e^λ
    sinh_half: T,
    /// Center of dilation (None = origin)
    center: Option<Point<T>>,
}

impl<T: Float> Dilator<T> {
    /// Creates a dilator for scaling by factor σ centered at origin.
    pub fn from_scale(factor: T) -> Self {
        // σ = e^λ => λ = ln(σ)
        let lambda = factor.ln();
        let half = lambda / T::TWO;
        Self {
            cosh_half: half.cosh(),
            sinh_half: half.sinh(),
            center: None,
        }
    }

    /// Creates a dilator for scaling by factor σ centered at a point.
    pub fn from_center_scale(center: &Point<T>, factor: T) -> Self {
        let lambda = factor.ln();
        let half = lambda / T::TWO;
        Self {
            cosh_half: half.cosh(),
            sinh_half: half.sinh(),
            center: Some(*center),
        }
    }

    /// Identity dilator (scale factor 1).
    pub fn identity() -> Self {
        Self {
            cosh_half: T::one(),
            sinh_half: T::zero(),
            center: None,
        }
    }

    /// Returns the scale factor σ.
    pub fn scale_factor(&self) -> T {
        // σ = e^λ = e^(2·atanh(sinh/cosh))
        // For small angles: σ ≈ (1 + sinh_half) / (1 - sinh_half)
        let lambda = T::TWO * (self.sinh_half / self.cosh_half).atanh();
        lambda.exp()
    }

    /// Returns the reverse (inverse) dilator.
    ///
    /// The inverse scales by 1/σ.
    pub fn reverse(&self) -> Self {
        Self {
            cosh_half: self.cosh_half,
            sinh_half: -self.sinh_half,
            center: self.center,
        }
    }

    /// Composes two dilators (multiplication of scale factors).
    pub fn compose(&self, other: &Self) -> Self {
        // For same-center dilators, factors multiply
        // For different centers, this is more complex
        todo!("Derive from SymPy: derive_dilator_composition()")
    }

    /// Transforms a point.
    pub fn transform_point(&self, p: &Point<T>) -> Point<T> {
        todo!("Derive from SymPy: derive_dilator_transform_point()")
    }

    /// Transforms a sphere (scales radius and translates center).
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

### 2. Inversor Type (`dim3/inversor.rs`)

```rust
use crate::scalar::Float;
use super::{Point, Sphere, Plane, Circle, Line};

/// An inversor in 3D CGA (inversion in a sphere).
///
/// Inverts geometry through a sphere: maps point P to P' where
/// |CP| · |CP'| = r² (C = center, r = radius).
///
/// # Properties
///
/// - Self-inverse: applying twice returns the original
/// - Maps circles/lines to circles/lines
/// - Maps spheres/planes to spheres/planes
/// - The center point maps to infinity
///
/// # Example
///
/// ```
/// use clifford::specialized::conformal::dim3::{Point, Sphere, Inversor};
///
/// // Inversion in unit sphere at origin
/// let inv = Inversor::from_sphere(&Sphere::from_center_radius(0.0, 0.0, 0.0, 1.0));
///
/// // Point at distance 2 maps to distance 1/2
/// let p = Point::new(2.0, 0.0, 0.0);
/// let q = inv.transform_point(&p);
///
/// assert!((q.x() - 0.5).abs() < 1e-10);
/// ```
///
/// Reference: https://conformalgeometricalgebra.org/wiki/index.php?title=Inversion
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Inversor<T: Float> {
    /// The sphere of inversion
    sphere: Sphere<T>,
}

impl<T: Float> Inversor<T> {
    /// Creates an inversor for the given sphere.
    pub fn from_sphere(sphere: &Sphere<T>) -> Self {
        Self { sphere: *sphere }
    }

    /// Creates an inversor for a unit sphere at origin.
    pub fn unit_sphere() -> Self {
        Self::from_sphere(&Sphere::from_center_radius(
            T::zero(), T::zero(), T::zero(), T::one()
        ))
    }

    /// Creates an inversor for a sphere with given center and radius.
    pub fn from_center_radius(cx: T, cy: T, cz: T, r: T) -> Self {
        Self::from_sphere(&Sphere::from_center_radius(cx, cy, cz, r))
    }

    /// Returns the sphere of inversion.
    pub fn sphere(&self) -> &Sphere<T> {
        &self.sphere
    }

    /// Inversion is self-inverse: reverse returns itself.
    pub fn reverse(&self) -> Self {
        *self
    }

    /// Transforms a point (inversion formula).
    ///
    /// For sphere with center C and radius r:
    /// P' = C + r²(P - C) / |P - C|²
    pub fn transform_point(&self, p: &Point<T>) -> Point<T> {
        todo!("Derive from SymPy: derive_inversor_transform_point()")
    }

    /// Transforms a sphere.
    ///
    /// Spheres map to spheres (or planes if passing through center).
    pub fn transform_sphere(&self, s: &Sphere<T>) -> Sphere<T> {
        todo!("Derive from SymPy")
    }

    /// Transforms a plane.
    ///
    /// Planes map to spheres through the inversion center
    /// (or planes if passing through center).
    pub fn transform_plane(&self, plane: &Plane<T>) -> Sphere<T> {
        todo!("Derive from SymPy")
    }

    /// Transforms a circle.
    pub fn transform_circle(&self, c: &Circle<T>) -> Circle<T> {
        todo!("Derive from SymPy")
    }

    /// Transforms a line.
    ///
    /// Lines map to circles through the inversion center
    /// (or lines if passing through center).
    pub fn transform_line(&self, l: &Line<T>) -> Circle<T> {
        todo!("Derive from SymPy")
    }
}
```

### 3. Transversor Type (`dim3/transversor.rs`)

```rust
use crate::scalar::Float;
use super::{Point, Translator, Inversor};

/// A transversor in 3D CGA (special conformal transformation).
///
/// The transversor is a combination of:
/// 1. Inversion in unit sphere at origin
/// 2. Translation by vector τ
/// 3. Inversion in unit sphere at origin
///
/// This is the "translation at infinity" - it maps the point at
/// infinity in one direction to infinity in another direction.
///
/// # Example
///
/// ```
/// use clifford::specialized::conformal::dim3::{Point, Transversor};
///
/// let trans = Transversor::new(1.0, 0.0, 0.0);
/// let p = Point::new(1.0, 0.0, 0.0);
/// let q = trans.transform_point(&p);
/// ```
///
/// Reference: https://conformalgeometricalgebra.org/wiki/index.php?title=Transversion
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Transversor<T: Float> {
    /// Translation vector τ
    tau_x: T,
    tau_y: T,
    tau_z: T,
}

impl<T: Float> Transversor<T> {
    /// Creates a transversor with translation vector τ.
    pub fn new(tau_x: T, tau_y: T, tau_z: T) -> Self {
        Self { tau_x, tau_y, tau_z }
    }

    /// Identity transversor (τ = 0).
    pub fn identity() -> Self {
        Self::new(T::zero(), T::zero(), T::zero())
    }

    /// Returns the translation vector.
    pub fn tau(&self) -> (T, T, T) {
        (self.tau_x, self.tau_y, self.tau_z)
    }

    /// Returns the reverse (inverse) transversor.
    pub fn reverse(&self) -> Self {
        Self::new(-self.tau_x, -self.tau_y, -self.tau_z)
    }

    /// Transforms a point.
    ///
    /// K * P = I * T * I * P where I = inversion, T = translation
    pub fn transform_point(&self, p: &Point<T>) -> Point<T> {
        // P' = (P + τ|P|²) / (1 + 2τ·P + |τ|²|P|²)
        todo!("Derive from SymPy: derive_transversor_transform()")
    }
}
```

## SymPy Derivations

Add to `derivations/src/clifford_derivations/cga.py`:

```python
@with_timeout(300)
def derive_dilator_transform_point():
    """Derive D ⋈ P ⋈ D̃ for dilator centered at origin.

    For scale factor σ = e^λ, the dilator is:
    D = cosh(λ/2) + sinh(λ/2)·e₀∞

    Result: Point at (x, y, z) maps to (σx, σy, σz)
    """
    sigma = symbols('sigma', positive=True)
    x, y, z = symbols('x y z')

    # For origin-centered dilation: P' = σ·P
    print("Dilator transform (origin-centered):")
    print(f"  (x, y, z) -> (σx, σy, σz)")


@with_timeout(300)
def derive_inversor_transform_point():
    """Derive sphere inversion of a point.

    For sphere with center C = (cx, cy, cz) and radius r:
    P' = C + r²(P - C) / |P - C|²
    """
    cx, cy, cz, r = symbols('cx cy cz r')
    px, py, pz = symbols('px py pz')

    # Vector from center to point
    dx = px - cx
    dy = py - cy
    dz = pz - cz
    dist_sq = dx**2 + dy**2 + dz**2

    # Inverted position
    factor = r**2 / dist_sq
    qx = cx + dx * factor
    qy = cy + dy * factor
    qz = cz + dz * factor

    print("Inversor transform:")
    print(f"  qx = {rust_code(qx)}")
    print(f"  qy = {rust_code(qy)}")
    print(f"  qz = {rust_code(qz)}")


@with_timeout(300)
def derive_transversor_transform():
    """Derive special conformal transformation.

    K(P) = (P + τ|P|²) / (1 + 2τ·P + |τ|²|P|²)
    """
    tx, ty, tz = symbols('tx ty tz')
    px, py, pz = symbols('px py pz')

    p_sq = px**2 + py**2 + pz**2
    tau_sq = tx**2 + ty**2 + tz**2
    tau_dot_p = tx*px + ty*py + tz*pz

    denom = 1 + 2*tau_dot_p + tau_sq * p_sq

    qx = (px + tx * p_sq) / denom
    qy = (py + ty * p_sq) / denom
    qz = (pz + tz * p_sq) / denom

    print("Transversor transform:")
    print(f"  qx = {rust_code(qx)}")
    print(f"  qy = {rust_code(qy)}")
    print(f"  qz = {rust_code(qz)}")
```

## Property Tests

```rust
proptest! {
    // ================================================================
    // Dilator tests
    // ================================================================

    #[test]
    fn dilator_scales_distance_from_origin(
        factor in 0.1f64..10.0,
        px in -10.0f64..10.0, py in -10.0f64..10.0, pz in -10.0f64..10.0,
    ) {
        let d = Dilator::from_scale(factor);
        let p = Point::new(px, py, pz);
        let q = d.transform_point(&p);

        // Distance from origin should scale by factor
        let dist_before = (px*px + py*py + pz*pz).sqrt();
        let dist_after = (q.x()*q.x() + q.y()*q.y() + q.z()*q.z()).sqrt();

        prop_assert!(abs_diff_eq!(dist_after, factor * dist_before, epsilon = ABS_DIFF_EQ_EPS));
    }

    #[test]
    fn dilator_inverse(
        factor in 0.1f64..10.0,
        px in -10.0f64..10.0, py in -10.0f64..10.0, pz in -10.0f64..10.0,
    ) {
        let d = Dilator::from_scale(factor);
        let p = Point::new(px, py, pz);

        let q = d.transform_point(&p);
        let back = d.reverse().transform_point(&q);

        prop_assert!(abs_diff_eq!(back.x(), p.x(), epsilon = ABS_DIFF_EQ_EPS));
        prop_assert!(abs_diff_eq!(back.y(), p.y(), epsilon = ABS_DIFF_EQ_EPS));
        prop_assert!(abs_diff_eq!(back.z(), p.z(), epsilon = ABS_DIFF_EQ_EPS));
    }

    #[test]
    fn dilator_preserves_angles(
        factor in 0.1f64..10.0,
        p1 in any::<Point<f64>>(),
        p2 in any::<Point<f64>>(),
        p3 in any::<Point<f64>>(),
    ) {
        // Conformal transformations preserve angles
        let d = Dilator::from_scale(factor);

        let q1 = d.transform_point(&p1);
        let q2 = d.transform_point(&p2);
        let q3 = d.transform_point(&p3);

        // Angle at p2 between p1-p2-p3
        let angle_before = angle_at_vertex(&p1, &p2, &p3);
        let angle_after = angle_at_vertex(&q1, &q2, &q3);

        prop_assert!(abs_diff_eq!(angle_before, angle_after, epsilon = ABS_DIFF_EQ_EPS));
    }

    // ================================================================
    // Inversor tests
    // ================================================================

    #[test]
    fn inversor_is_self_inverse(
        cx in -5.0f64..5.0, cy in -5.0f64..5.0, cz in -5.0f64..5.0,
        r in 0.1f64..5.0,
        px in -10.0f64..10.0, py in -10.0f64..10.0, pz in -10.0f64..10.0,
    ) {
        // Skip if point is at center (maps to infinity)
        let dist = ((px-cx).powi(2) + (py-cy).powi(2) + (pz-cz).powi(2)).sqrt();
        if dist < 0.1 {
            return Ok(());
        }

        let inv = Inversor::from_center_radius(cx, cy, cz, r);
        let p = Point::new(px, py, pz);

        let q = inv.transform_point(&p);
        let back = inv.transform_point(&q);

        prop_assert!(abs_diff_eq!(back.x(), p.x(), epsilon = ABS_DIFF_EQ_EPS));
        prop_assert!(abs_diff_eq!(back.y(), p.y(), epsilon = ABS_DIFF_EQ_EPS));
        prop_assert!(abs_diff_eq!(back.z(), p.z(), epsilon = ABS_DIFF_EQ_EPS));
    }

    #[test]
    fn inversor_product_of_distances(
        r in 0.1f64..5.0,
        px in 0.1f64..10.0, py in -10.0f64..10.0, pz in -10.0f64..10.0,
    ) {
        // For inversion in sphere at origin with radius r:
        // |P| · |P'| = r²

        let inv = Inversor::from_center_radius(0.0, 0.0, 0.0, r);
        let p = Point::new(px, py, pz);
        let q = inv.transform_point(&p);

        let dist_p = (px*px + py*py + pz*pz).sqrt();
        let dist_q = (q.x()*q.x() + q.y()*q.y() + q.z()*q.z()).sqrt();

        prop_assert!(abs_diff_eq!(dist_p * dist_q, r * r, epsilon = ABS_DIFF_EQ_EPS));
    }

    #[test]
    fn inversor_preserves_angles(
        r in 0.1f64..5.0,
        p1 in any::<Point<f64>>(),
        p2 in any::<Point<f64>>(),
        p3 in any::<Point<f64>>(),
    ) {
        // Conformal transformations preserve angles
        let inv = Inversor::from_center_radius(0.0, 0.0, 0.0, r);

        // Skip if any point too close to origin
        if p1.distance(&Point::origin()) < 0.1
            || p2.distance(&Point::origin()) < 0.1
            || p3.distance(&Point::origin()) < 0.1 {
            return Ok(());
        }

        let q1 = inv.transform_point(&p1);
        let q2 = inv.transform_point(&p2);
        let q3 = inv.transform_point(&p3);

        let angle_before = angle_at_vertex(&p1, &p2, &p3);
        let angle_after = angle_at_vertex(&q1, &q2, &q3);

        prop_assert!(abs_diff_eq!(angle_before, angle_after, epsilon = ABS_DIFF_EQ_EPS));
    }

    // ================================================================
    // Transversor tests
    // ================================================================

    #[test]
    fn transversor_inverse(
        tx in -1.0f64..1.0, ty in -1.0f64..1.0, tz in -1.0f64..1.0,
        px in -5.0f64..5.0, py in -5.0f64..5.0, pz in -5.0f64..5.0,
    ) {
        let k = Transversor::new(tx, ty, tz);
        let p = Point::new(px, py, pz);

        let q = k.transform_point(&p);
        let back = k.reverse().transform_point(&q);

        prop_assert!(abs_diff_eq!(back.x(), p.x(), epsilon = ABS_DIFF_EQ_EPS));
        prop_assert!(abs_diff_eq!(back.y(), p.y(), epsilon = ABS_DIFF_EQ_EPS));
        prop_assert!(abs_diff_eq!(back.z(), p.z(), epsilon = ABS_DIFF_EQ_EPS));
    }
}
```

## Files to Create/Modify

### New Files
- `src/specialized/conformal/dim3/dilator.rs`
- `src/specialized/conformal/dim3/inversor.rs`
- `src/specialized/conformal/dim3/transversor.rs`

### Modified Files
- `src/specialized/conformal/dim3/mod.rs` - Export new types
- `derivations/src/clifford_derivations/cga.py` - Add derivations

## Verification Checklist

- [ ] `cargo fmt` passes
- [ ] `cargo clippy --all-features` passes
- [ ] `cargo doc --all-features --no-deps` passes
- [ ] `cargo test --all-features` passes
- [ ] `cargo deny check` passes
- [ ] Dilator scales distances correctly
- [ ] Inversor is self-inverse
- [ ] |P| · |P'| = r² for inversion
- [ ] Angles preserved by all conformal transforms

## Dependencies

- PRD-6.5 (Translator and Rotor Versors) - must be complete

## Next Steps

After this PRD is complete, proceed to PRD-6.7 (Operations and Integration).
