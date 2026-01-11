//! Type definitions for 3D Conformal Geometric Algebra.

use approx::{AbsDiffEq, RelativeEq, UlpsEq};

use crate::scalar::Float;

/// A round point in 3D CGA (null vector).
///
/// Represents a Euclidean point `(x, y, z)` embedded in the conformal model:
///
/// ```text
/// P = x·e₁ + y·e₂ + z·e₃ + e₀ + ½(x² + y² + z²)·e∞
/// ```
///
/// The null constraint `P · P = 0` is maintained by construction.
///
/// # Internal Representation
///
/// Internally, we store coefficients in the orthogonal basis `(e₊, e₋)` rather
/// than the null basis `(e₀, e∞)`:
/// - `ep`: coefficient of `e₊`
/// - `em`: coefficient of `e₋`
///
/// The relationship is:
/// - `e₀ = (e₋ - e₊) / 2`
/// - `e∞ = e₋ + e₊`
///
/// # Example
///
/// ```
/// use clifford::specialized::conformal::dim3::Point;
///
/// // Create a point at (1, 2, 3)
/// let p = Point::<f64>::new(1.0, 2.0, 3.0);
///
/// // Verify null constraint
/// assert!(p.is_null(1e-10));
///
/// // Distance between points
/// let q = Point::new(4.0, 6.0, 3.0);
/// let d = p.distance(&q);
/// assert!((d - 5.0).abs() < 1e-10);  // 3-4-5 triangle
/// ```
///
/// # References
///
/// - <https://conformalgeometricalgebra.org/wiki/index.php?title=Round_point>
#[derive(Clone, Copy, Debug, PartialEq)]
#[repr(C)]
pub struct Point<T: Float> {
    /// Coefficient of `e₁` (x-coordinate).
    e1: T,
    /// Coefficient of `e₂` (y-coordinate).
    e2: T,
    /// Coefficient of `e₃` (z-coordinate).
    e3: T,
    /// Coefficient of `e₊` (positive conformal basis).
    ep: T,
    /// Coefficient of `e₋` (negative conformal basis).
    em: T,
}

impl<T: Float> Point<T> {
    /// Creates a point at Euclidean coordinates `(x, y, z)`.
    ///
    /// The conformal embedding ensures `P · P = 0`.
    ///
    /// # Formula
    ///
    /// ```text
    /// P = x·e₁ + y·e₂ + z·e₃ + e₀ + ½(x² + y² + z²)·e∞
    /// ```
    ///
    /// In orthogonal basis:
    /// ```text
    /// ep = (r² - 1) / 2
    /// em = (r² + 1) / 2
    /// ```
    /// where `r² = x² + y² + z²`.
    ///
    /// # Example
    ///
    /// ```
    /// use clifford::specialized::conformal::dim3::Point;
    ///
    /// let p = Point::<f64>::new(1.0, 2.0, 3.0);
    /// assert!((p.x() - 1.0).abs() < 1e-10);
    /// ```
    #[inline]
    pub fn new(x: T, y: T, z: T) -> Self {
        // P = x·e₁ + y·e₂ + z·e₃ + e₀ + ½(x² + y² + z²)·e∞
        // where e₀ = (e₋ - e₊)/2, e∞ = e₋ + e₊
        //
        // ep coefficient: -1/2 + (x² + y² + z²)/2 = (r² - 1)/2
        // em coefficient: 1/2 + (x² + y² + z²)/2 = (r² + 1)/2
        let r_sq = x * x + y * y + z * z;
        let half = T::one() / T::TWO;
        Self {
            e1: x,
            e2: y,
            e3: z,
            ep: (r_sq - T::one()) * half,
            em: (r_sq + T::one()) * half,
        }
    }

    /// Creates a point from raw conformal components (unchecked).
    ///
    /// Use this when you know the components satisfy the null constraint,
    /// or for internal operations where the constraint will be restored.
    ///
    /// # Safety
    ///
    /// This does not verify `P · P = 0`. Use [`Point::new`] for safe construction.
    ///
    /// # Panics
    ///
    /// The accessor methods [`x()`](Self::x), [`y()`](Self::y), [`z()`](Self::z)
    /// will panic or return infinity if the point has zero weight (i.e., `em - ep = 0`).
    /// Points created via [`Point::new`] always have unit weight.
    #[inline]
    pub fn from_conformal_unchecked(e1: T, e2: T, e3: T, ep: T, em: T) -> Self {
        Self { e1, e2, e3, ep, em }
    }

    /// The origin point `(0, 0, 0)`.
    ///
    /// # Example
    ///
    /// ```
    /// use clifford::specialized::conformal::dim3::Point;
    ///
    /// let o = Point::<f64>::origin();
    /// assert!((o.x()).abs() < 1e-10);
    /// assert!((o.y()).abs() < 1e-10);
    /// assert!((o.z()).abs() < 1e-10);
    /// ```
    #[inline]
    pub fn origin() -> Self {
        Self::new(T::zero(), T::zero(), T::zero())
    }

    /// Extracts the Euclidean x-coordinate.
    ///
    /// For a properly normalized point (weight = 1), this returns the x-coordinate.
    /// For non-unit weight points, the coordinate is divided by the weight.
    #[inline]
    pub fn x(&self) -> T {
        // Weight is the coefficient of e₀, which equals (em - ep)
        // We divide by the weight to get Cartesian coordinates
        let weight = self.em - self.ep;
        self.e1 / weight
    }

    /// Extracts the Euclidean y-coordinate.
    #[inline]
    pub fn y(&self) -> T {
        let weight = self.em - self.ep;
        self.e2 / weight
    }

    /// Extracts the Euclidean z-coordinate.
    #[inline]
    pub fn z(&self) -> T {
        let weight = self.em - self.ep;
        self.e3 / weight
    }

    /// Returns the weight (coefficient of `e₀`).
    ///
    /// For points created via [`Point::new`], this is always 1.
    /// Non-unit weights occur when points are scaled or combined.
    #[inline]
    pub fn weight(&self) -> T {
        self.em - self.ep
    }

    /// Checks if this is a valid null vector (`P · P = 0`).
    ///
    /// All properly constructed conformal points are null.
    ///
    /// # Example
    ///
    /// ```
    /// use clifford::specialized::conformal::dim3::Point;
    ///
    /// let p = Point::new(1.0, 2.0, 3.0);
    /// assert!(p.is_null(1e-10));
    /// ```
    #[inline]
    pub fn is_null(&self, epsilon: T) -> bool {
        // P · P = e1² + e2² + e3² + ep² - em²
        //       = r² + ep² - em²
        //
        // For numerical stability, we use the identity:
        //   ep² - em² = (ep - em)(ep + em)
        //
        // For a valid conformal point with weight 1:
        //   ep - em = -1, ep + em = r², so ep² - em² = -r²
        //   Thus P · P = r² - r² = 0
        let r_sq = self.e1 * self.e1 + self.e2 * self.e2 + self.e3 * self.e3;
        let ep_sq_minus_em_sq = (self.ep - self.em) * (self.ep + self.em);
        let norm_sq = r_sq + ep_sq_minus_em_sq;

        // Use relative tolerance scaled by magnitude for numerical stability
        // with large coordinates
        let scale = T::one() + r_sq;
        norm_sq.abs() < epsilon * scale
    }

    /// Euclidean distance to another point.
    ///
    /// Uses the inner product formula: `d² = -2(P₁ · P₂)` for unit-weight points.
    ///
    /// # Formula
    ///
    /// ```text
    /// P₁ · P₂ = e1₁·e1₂ + e2₁·e2₂ + e3₁·e3₂ + ep₁·ep₂ - em₁·em₂
    /// d² = -2(P₁ · P₂) / (w₁ · w₂)
    /// ```
    ///
    /// # Example
    ///
    /// ```
    /// use clifford::specialized::conformal::dim3::Point;
    ///
    /// let p1 = Point::<f64>::new(0.0, 0.0, 0.0);
    /// let p2 = Point::new(3.0, 4.0, 0.0);
    /// let d = p1.distance(&p2);
    /// assert!((d - 5.0).abs() < 1e-10);  // 3-4-5 triangle
    /// ```
    pub fn distance(&self, other: &Point<T>) -> T {
        // For numerical stability, compute the Euclidean distance directly
        // rather than using the CGA inner product formula (which can have
        // precision issues with large coordinate values).
        //
        // For conformal points, we need to normalize by weights first
        let w1 = self.weight();
        let w2 = other.weight();

        let x1 = self.e1 / w1;
        let y1 = self.e2 / w1;
        let z1 = self.e3 / w1;
        let x2 = other.e1 / w2;
        let y2 = other.e2 / w2;
        let z2 = other.e3 / w2;

        let dx = x1 - x2;
        let dy = y1 - y2;
        let dz = z1 - z2;

        (dx * dx + dy * dy + dz * dz).sqrt()
    }

    /// Inner product with another point.
    ///
    /// For conformal points: `P₁ · P₂ = -½|p₁ - p₂|²` when both have unit weight.
    #[inline]
    pub fn inner(&self, other: &Point<T>) -> T {
        self.e1 * other.e1 + self.e2 * other.e2 + self.e3 * other.e3 + self.ep * other.ep
            - self.em * other.em
    }

    /// Accessor for `e₁` component.
    #[inline]
    pub fn e1(&self) -> T {
        self.e1
    }

    /// Accessor for `e₂` component.
    #[inline]
    pub fn e2(&self) -> T {
        self.e2
    }

    /// Accessor for `e₃` component.
    #[inline]
    pub fn e3(&self) -> T {
        self.e3
    }

    /// Accessor for `e₊` component.
    #[inline]
    pub fn ep(&self) -> T {
        self.ep
    }

    /// Accessor for `e₋` component.
    #[inline]
    pub fn em(&self) -> T {
        self.em
    }
}

impl<T: Float> Default for Point<T> {
    fn default() -> Self {
        Self::origin()
    }
}

// ============================================================================
// Sphere type
// ============================================================================

/// A sphere in 3D CGA.
///
/// Represents a sphere with center `(cx, cy, cz)` and radius `r`.
///
/// # Representation
///
/// In the dual representation: `S* = C - ½r²·e∞` where `C` is the center point.
///
/// # Real vs Imaginary
///
/// - `r² > 0`: Real sphere with positive radius
/// - `r² = 0`: Point sphere (degenerate to a point)
/// - `r² < 0`: Imaginary sphere (no real surface points)
///
/// # Example
///
/// ```
/// use clifford::specialized::conformal::dim3::{Point, Sphere};
///
/// // Unit sphere at origin
/// let sphere = Sphere::<f64>::from_center_radius(0.0, 0.0, 0.0, 1.0);
///
/// // Point on sphere surface
/// let p = Point::new(1.0, 0.0, 0.0);
/// assert!(sphere.contains(&p, 1e-10));
///
/// // Point inside sphere
/// let q = Point::new(0.5, 0.0, 0.0);
/// assert!(sphere.signed_distance(&q) < 0.0);
/// ```
///
/// # References
///
/// - <https://conformalgeometricalgebra.org/wiki/index.php?title=Sphere>
#[derive(Clone, Copy, Debug, PartialEq)]
#[repr(C)]
pub struct Sphere<T: Float> {
    /// Center x-coordinate.
    cx: T,
    /// Center y-coordinate.
    cy: T,
    /// Center z-coordinate.
    cz: T,
    /// Squared radius (can be negative for imaginary spheres).
    r_sq: T,
}

impl<T: Float> Sphere<T> {
    /// Creates a sphere from center coordinates and radius.
    ///
    /// # Example
    ///
    /// ```
    /// use clifford::specialized::conformal::dim3::Sphere;
    ///
    /// let sphere = Sphere::<f64>::from_center_radius(1.0, 2.0, 3.0, 5.0);
    /// assert!((sphere.radius() - 5.0).abs() < 1e-10);
    /// ```
    #[inline]
    pub fn from_center_radius(cx: T, cy: T, cz: T, r: T) -> Self {
        Self {
            cx,
            cy,
            cz,
            r_sq: r * r,
        }
    }

    /// Creates a sphere from center coordinates and squared radius.
    ///
    /// This allows creating imaginary spheres with `r_sq < 0`.
    #[inline]
    pub fn from_center_radius_squared(cx: T, cy: T, cz: T, r_sq: T) -> Self {
        Self { cx, cy, cz, r_sq }
    }

    /// Creates a unit sphere centered at the origin.
    #[inline]
    pub fn unit() -> Self {
        Self::from_center_radius(T::zero(), T::zero(), T::zero(), T::one())
    }

    /// Returns the center as Euclidean coordinates.
    #[inline]
    pub fn center(&self) -> (T, T, T) {
        (self.cx, self.cy, self.cz)
    }

    /// Returns the center x-coordinate.
    #[inline]
    pub fn cx(&self) -> T {
        self.cx
    }

    /// Returns the center y-coordinate.
    #[inline]
    pub fn cy(&self) -> T {
        self.cy
    }

    /// Returns the center z-coordinate.
    #[inline]
    pub fn cz(&self) -> T {
        self.cz
    }

    /// Returns the squared radius.
    ///
    /// - Positive: real sphere
    /// - Zero: point sphere
    /// - Negative: imaginary sphere
    #[inline]
    pub fn radius_squared(&self) -> T {
        self.r_sq
    }

    /// Returns the radius.
    ///
    /// # Panics
    ///
    /// Panics if the sphere is imaginary (`r² < 0`).
    #[inline]
    pub fn radius(&self) -> T {
        assert!(
            self.r_sq >= T::zero(),
            "Cannot take radius of imaginary sphere"
        );
        self.r_sq.sqrt()
    }

    /// Returns true if this is a real sphere (`r² > 0`).
    #[inline]
    pub fn is_real(&self, epsilon: T) -> bool {
        self.r_sq > epsilon
    }

    /// Returns true if this is a point sphere (`r² ≈ 0`).
    #[inline]
    pub fn is_point(&self, epsilon: T) -> bool {
        self.r_sq.abs() < epsilon
    }

    /// Returns true if this is an imaginary sphere (`r² < 0`).
    #[inline]
    pub fn is_imaginary(&self, epsilon: T) -> bool {
        self.r_sq < -epsilon
    }

    /// Returns true if the point lies on the sphere surface.
    ///
    /// # Example
    ///
    /// ```
    /// use clifford::specialized::conformal::dim3::{Point, Sphere};
    ///
    /// let sphere = Sphere::from_center_radius(0.0, 0.0, 0.0, 1.0);
    /// let on_surface = Point::new(1.0, 0.0, 0.0);
    /// assert!(sphere.contains(&on_surface, 1e-10));
    /// ```
    #[inline]
    pub fn contains(&self, p: &Point<T>, epsilon: T) -> bool {
        self.signed_distance(p).abs() < epsilon
    }

    /// Signed distance from point to sphere surface.
    ///
    /// - Negative: inside sphere
    /// - Zero: on surface
    /// - Positive: outside sphere
    ///
    /// # Example
    ///
    /// ```
    /// use clifford::specialized::conformal::dim3::{Point, Sphere};
    ///
    /// let sphere = Sphere::from_center_radius(0.0, 0.0, 0.0, 1.0);
    ///
    /// // Inside
    /// let p_in = Point::new(0.5, 0.0, 0.0);
    /// assert!(sphere.signed_distance(&p_in) < 0.0);
    ///
    /// // Outside
    /// let p_out = Point::new(2.0, 0.0, 0.0);
    /// assert!(sphere.signed_distance(&p_out) > 0.0);
    /// ```
    ///
    /// # Panics
    ///
    /// Panics if the sphere is imaginary (`r² < 0`).
    pub fn signed_distance(&self, p: &Point<T>) -> T {
        // Distance from point to center minus radius
        let px = p.x();
        let py = p.y();
        let pz = p.z();

        let dx = px - self.cx;
        let dy = py - self.cy;
        let dz = pz - self.cz;

        let dist_to_center = (dx * dx + dy * dy + dz * dz).sqrt();
        // Use radius() which panics with a clear message for imaginary spheres
        dist_to_center - self.radius()
    }

    /// Returns the center as a conformal [`Point`].
    #[inline]
    pub fn center_point(&self) -> Point<T> {
        Point::new(self.cx, self.cy, self.cz)
    }
}

impl<T: Float> Default for Sphere<T> {
    fn default() -> Self {
        Self::unit()
    }
}

// ============================================================================
// approx trait implementations
// ============================================================================

impl<T: Float> AbsDiffEq for Point<T> {
    type Epsilon = T;

    fn default_epsilon() -> Self::Epsilon {
        T::default_epsilon()
    }

    fn abs_diff_eq(&self, other: &Self, epsilon: Self::Epsilon) -> bool {
        T::abs_diff_eq(&self.e1, &other.e1, epsilon)
            && T::abs_diff_eq(&self.e2, &other.e2, epsilon)
            && T::abs_diff_eq(&self.e3, &other.e3, epsilon)
            && T::abs_diff_eq(&self.ep, &other.ep, epsilon)
            && T::abs_diff_eq(&self.em, &other.em, epsilon)
    }
}

impl<T: Float> RelativeEq for Point<T> {
    fn default_max_relative() -> Self::Epsilon {
        T::default_max_relative()
    }

    fn relative_eq(
        &self,
        other: &Self,
        epsilon: Self::Epsilon,
        max_relative: Self::Epsilon,
    ) -> bool {
        T::relative_eq(&self.e1, &other.e1, epsilon, max_relative)
            && T::relative_eq(&self.e2, &other.e2, epsilon, max_relative)
            && T::relative_eq(&self.e3, &other.e3, epsilon, max_relative)
            && T::relative_eq(&self.ep, &other.ep, epsilon, max_relative)
            && T::relative_eq(&self.em, &other.em, epsilon, max_relative)
    }
}

impl<T: Float> UlpsEq for Point<T> {
    fn default_max_ulps() -> u32 {
        T::default_max_ulps()
    }

    fn ulps_eq(&self, other: &Self, epsilon: Self::Epsilon, max_ulps: u32) -> bool {
        T::ulps_eq(&self.e1, &other.e1, epsilon, max_ulps)
            && T::ulps_eq(&self.e2, &other.e2, epsilon, max_ulps)
            && T::ulps_eq(&self.e3, &other.e3, epsilon, max_ulps)
            && T::ulps_eq(&self.ep, &other.ep, epsilon, max_ulps)
            && T::ulps_eq(&self.em, &other.em, epsilon, max_ulps)
    }
}

impl<T: Float> AbsDiffEq for Sphere<T> {
    type Epsilon = T;

    fn default_epsilon() -> Self::Epsilon {
        T::default_epsilon()
    }

    fn abs_diff_eq(&self, other: &Self, epsilon: Self::Epsilon) -> bool {
        T::abs_diff_eq(&self.cx, &other.cx, epsilon)
            && T::abs_diff_eq(&self.cy, &other.cy, epsilon)
            && T::abs_diff_eq(&self.cz, &other.cz, epsilon)
            && T::abs_diff_eq(&self.r_sq, &other.r_sq, epsilon)
    }
}

impl<T: Float> RelativeEq for Sphere<T> {
    fn default_max_relative() -> Self::Epsilon {
        T::default_max_relative()
    }

    fn relative_eq(
        &self,
        other: &Self,
        epsilon: Self::Epsilon,
        max_relative: Self::Epsilon,
    ) -> bool {
        T::relative_eq(&self.cx, &other.cx, epsilon, max_relative)
            && T::relative_eq(&self.cy, &other.cy, epsilon, max_relative)
            && T::relative_eq(&self.cz, &other.cz, epsilon, max_relative)
            && T::relative_eq(&self.r_sq, &other.r_sq, epsilon, max_relative)
    }
}

impl<T: Float> UlpsEq for Sphere<T> {
    fn default_max_ulps() -> u32 {
        T::default_max_ulps()
    }

    fn ulps_eq(&self, other: &Self, epsilon: Self::Epsilon, max_ulps: u32) -> bool {
        T::ulps_eq(&self.cx, &other.cx, epsilon, max_ulps)
            && T::ulps_eq(&self.cy, &other.cy, epsilon, max_ulps)
            && T::ulps_eq(&self.cz, &other.cz, epsilon, max_ulps)
            && T::ulps_eq(&self.r_sq, &other.r_sq, epsilon, max_ulps)
    }
}

// ============================================================================
// Tests
// ============================================================================

#[cfg(test)]
mod tests {
    use super::*;
    use crate::test_utils::ABS_DIFF_EQ_EPS;
    use approx::abs_diff_eq;
    use proptest::prelude::*;

    // ========================================================================
    // Point null constraint
    // ========================================================================

    proptest! {
        #[test]
        fn point_is_null(
            x in -100.0f64..100.0,
            y in -100.0f64..100.0,
            z in -100.0f64..100.0,
        ) {
            let p = Point::new(x, y, z);
            prop_assert!(p.is_null(ABS_DIFF_EQ_EPS));
        }

        #[test]
        fn point_roundtrip_coordinates(
            x in -100.0f64..100.0,
            y in -100.0f64..100.0,
            z in -100.0f64..100.0,
        ) {
            let p = Point::new(x, y, z);
            prop_assert!(abs_diff_eq!(p.x(), x, epsilon = ABS_DIFF_EQ_EPS));
            prop_assert!(abs_diff_eq!(p.y(), y, epsilon = ABS_DIFF_EQ_EPS));
            prop_assert!(abs_diff_eq!(p.z(), z, epsilon = ABS_DIFF_EQ_EPS));
        }

        #[test]
        fn point_unit_weight(
            x in -100.0f64..100.0,
            y in -100.0f64..100.0,
            z in -100.0f64..100.0,
        ) {
            let p = Point::new(x, y, z);
            prop_assert!(abs_diff_eq!(p.weight(), 1.0, epsilon = ABS_DIFF_EQ_EPS));
        }
    }

    // ========================================================================
    // Distance
    // ========================================================================

    proptest! {
        #[test]
        fn distance_is_euclidean(
            x1 in -10.0f64..10.0, y1 in -10.0f64..10.0, z1 in -10.0f64..10.0,
            x2 in -10.0f64..10.0, y2 in -10.0f64..10.0, z2 in -10.0f64..10.0,
        ) {
            let p1 = Point::new(x1, y1, z1);
            let p2 = Point::new(x2, y2, z2);

            let cga_dist = p1.distance(&p2);
            let euclidean_dist = ((x2-x1).powi(2) + (y2-y1).powi(2) + (z2-z1).powi(2)).sqrt();

            prop_assert!(abs_diff_eq!(cga_dist, euclidean_dist, epsilon = ABS_DIFF_EQ_EPS));
        }

        #[test]
        fn distance_is_symmetric(
            x1 in -10.0f64..10.0, y1 in -10.0f64..10.0, z1 in -10.0f64..10.0,
            x2 in -10.0f64..10.0, y2 in -10.0f64..10.0, z2 in -10.0f64..10.0,
        ) {
            let p1 = Point::new(x1, y1, z1);
            let p2 = Point::new(x2, y2, z2);

            let d12 = p1.distance(&p2);
            let d21 = p2.distance(&p1);

            prop_assert!(abs_diff_eq!(d12, d21, epsilon = ABS_DIFF_EQ_EPS));
        }

        #[test]
        fn distance_to_self_is_zero(
            x in -100.0f64..100.0,
            y in -100.0f64..100.0,
            z in -100.0f64..100.0,
        ) {
            let p = Point::new(x, y, z);
            let d = p.distance(&p);
            prop_assert!(abs_diff_eq!(d, 0.0, epsilon = ABS_DIFF_EQ_EPS));
        }
    }

    // ========================================================================
    // Sphere containment
    // ========================================================================

    proptest! {
        #[test]
        fn sphere_contains_surface_points(
            cx in -10.0f64..10.0, cy in -10.0f64..10.0, cz in -10.0f64..10.0,
            r in 0.1f64..10.0,
            theta in 0.0f64..std::f64::consts::TAU,
            phi in 0.0f64..std::f64::consts::PI,
        ) {
            let sphere = Sphere::from_center_radius(cx, cy, cz, r);
            let p = Point::new(
                cx + r * phi.sin() * theta.cos(),
                cy + r * phi.sin() * theta.sin(),
                cz + r * phi.cos(),
            );
            prop_assert!(sphere.contains(&p, ABS_DIFF_EQ_EPS));
        }

        #[test]
        fn sphere_center_inside(
            cx in -10.0f64..10.0, cy in -10.0f64..10.0, cz in -10.0f64..10.0,
            r in 0.1f64..10.0,
        ) {
            let sphere = Sphere::from_center_radius(cx, cy, cz, r);
            let center = Point::new(cx, cy, cz);
            let signed_dist = sphere.signed_distance(&center);
            prop_assert!(signed_dist < 0.0);
            prop_assert!(abs_diff_eq!(signed_dist, -r, epsilon = ABS_DIFF_EQ_EPS));
        }

        #[test]
        fn sphere_outside_positive_distance(
            cx in -5.0f64..5.0, cy in -5.0f64..5.0, cz in -5.0f64..5.0,
            r in 0.1f64..2.0,
            dx in 3.0f64..10.0,
        ) {
            let sphere = Sphere::from_center_radius(cx, cy, cz, r);
            // Point definitely outside (dx > r guaranteed by ranges)
            let p = Point::new(cx + dx, cy, cz);
            let signed_dist = sphere.signed_distance(&p);
            prop_assert!(signed_dist > 0.0);
        }
    }

    // ========================================================================
    // Specific examples
    // ========================================================================

    #[test]
    fn origin_is_at_origin() {
        let o = Point::<f64>::origin();
        assert!(abs_diff_eq!(o.x(), 0.0, epsilon = ABS_DIFF_EQ_EPS));
        assert!(abs_diff_eq!(o.y(), 0.0, epsilon = ABS_DIFF_EQ_EPS));
        assert!(abs_diff_eq!(o.z(), 0.0, epsilon = ABS_DIFF_EQ_EPS));
        assert!(o.is_null(ABS_DIFF_EQ_EPS));
    }

    #[test]
    fn distance_3_4_5_triangle() {
        let p1 = Point::new(0.0, 0.0, 0.0);
        let p2 = Point::new(3.0, 4.0, 0.0);
        let d = p1.distance(&p2);
        assert!(abs_diff_eq!(d, 5.0, epsilon = ABS_DIFF_EQ_EPS));
    }

    #[test]
    fn unit_sphere_contains_unit_points() {
        let sphere = Sphere::<f64>::unit();

        // Points on unit sphere
        assert!(sphere.contains(&Point::new(1.0, 0.0, 0.0), ABS_DIFF_EQ_EPS));
        assert!(sphere.contains(&Point::new(0.0, 1.0, 0.0), ABS_DIFF_EQ_EPS));
        assert!(sphere.contains(&Point::new(0.0, 0.0, 1.0), ABS_DIFF_EQ_EPS));
        assert!(sphere.contains(&Point::new(-1.0, 0.0, 0.0), ABS_DIFF_EQ_EPS));
    }

    #[test]
    fn sphere_types() {
        let real = Sphere::from_center_radius(0.0f64, 0.0, 0.0, 1.0);
        assert!(real.is_real(ABS_DIFF_EQ_EPS));
        assert!(!real.is_point(ABS_DIFF_EQ_EPS));
        assert!(!real.is_imaginary(ABS_DIFF_EQ_EPS));

        let point = Sphere::from_center_radius(0.0f64, 0.0, 0.0, 0.0);
        assert!(!point.is_real(ABS_DIFF_EQ_EPS));
        assert!(point.is_point(ABS_DIFF_EQ_EPS));
        assert!(!point.is_imaginary(ABS_DIFF_EQ_EPS));

        let imaginary = Sphere::from_center_radius_squared(0.0f64, 0.0, 0.0, -1.0);
        assert!(!imaginary.is_real(ABS_DIFF_EQ_EPS));
        assert!(!imaginary.is_point(ABS_DIFF_EQ_EPS));
        assert!(imaginary.is_imaginary(ABS_DIFF_EQ_EPS));
    }
}
