//! Type definitions for 3D Conformal Geometric Algebra.

use approx::{AbsDiffEq, RelativeEq, UlpsEq};

use crate::scalar::Float;
use crate::specialized::euclidean::dim3::Vector;

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

    /// Returns the weight for coordinate normalization.
    ///
    /// This returns `em - ep`, which equals `2·e₀` where `e₀ = (e₋ - e₊)/2`.
    /// For points created via [`Point::new`], this is always 1.
    /// Non-unit weights occur when points are scaled or combined.
    ///
    /// The coordinate accessors [`x()`](Self::x), [`y()`](Self::y), [`z()`](Self::z)
    /// divide by this weight to extract Euclidean coordinates.
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

    /// Returns the center as a Euclidean vector.
    #[inline]
    pub fn center(&self) -> Vector<T> {
        Vector::new(self.cx, self.cy, self.cz)
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
    ///
    /// # Panics
    ///
    /// Panics if the sphere is imaginary (`r² < 0`).
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
// Plane type
// ============================================================================

/// A plane in 3D CGA (flat sphere - grade-4 quadvector).
///
/// Represents the plane `ax + by + cz + d = 0` in implicit form.
/// A plane is a sphere that contains the point at infinity.
///
/// # Example
///
/// ```
/// use clifford::specialized::conformal::dim3::{Point, Plane};
///
/// // xy-plane (z = 0)
/// let plane = Plane::<f64>::from_normal_distance(0.0, 0.0, 1.0, 0.0);
///
/// // Point on plane
/// let p = Point::new(1.0, 2.0, 0.0);
/// assert!(plane.contains(&p, 1e-10));
/// ```
///
/// # References
///
/// - <https://conformalgeometricalgebra.org/wiki/index.php?title=Plane>
#[derive(Clone, Copy, Debug, PartialEq)]
#[repr(C)]
pub struct Plane<T: Float> {
    // We store the implicit form (a, b, c, d) for simplicity.
    // The CGA representation can be derived when needed.
    /// Normal x-component.
    a: T,
    /// Normal y-component.
    b: T,
    /// Normal z-component.
    c: T,
    /// Distance term.
    d: T,
}

impl<T: Float> Plane<T> {
    /// Creates a plane from normalized normal `(a, b, c)` and signed distance `d`.
    ///
    /// The plane equation is `ax + by + cz + d = 0`.
    ///
    /// # Note
    ///
    /// The normal `(a, b, c)` should be a unit vector for proper distance calculations.
    #[inline]
    pub fn from_normal_distance(a: T, b: T, c: T, d: T) -> Self {
        Self { a, b, c, d }
    }

    /// Creates a plane from three non-collinear points.
    ///
    /// The normal direction follows the right-hand rule for `p1 -> p2 -> p3`.
    /// Returns `None` if the points are collinear (no unique plane exists).
    ///
    /// # Example
    ///
    /// ```
    /// use clifford::specialized::conformal::dim3::{Point, Plane};
    ///
    /// let p1 = Point::new(0.0, 0.0, 0.0);
    /// let p2 = Point::new(1.0, 0.0, 0.0);
    /// let p3 = Point::new(0.0, 1.0, 0.0);
    /// let plane = Plane::from_three_points(&p1, &p2, &p3).unwrap();
    /// assert!(plane.contains(&p1, 1e-10));
    /// ```
    pub fn from_three_points(p1: &Point<T>, p2: &Point<T>, p3: &Point<T>) -> Option<Self> {
        // Compute normal via cross product
        let v1 = (p2.x() - p1.x(), p2.y() - p1.y(), p2.z() - p1.z());
        let v2 = (p3.x() - p1.x(), p3.y() - p1.y(), p3.z() - p1.z());

        // Cross product: v1 × v2
        let nx = v1.1 * v2.2 - v1.2 * v2.1;
        let ny = v1.2 * v2.0 - v1.0 * v2.2;
        let nz = v1.0 * v2.1 - v1.1 * v2.0;

        // Normalize
        let len = (nx * nx + ny * ny + nz * nz).sqrt();

        // Collinear points: normal has zero length
        if len < T::epsilon() {
            return None;
        }

        let (a, b, c) = (nx / len, ny / len, nz / len);

        // d = -(a*x + b*y + c*z) for point p1 on plane
        let d = -(a * p1.x() + b * p1.y() + c * p1.z());

        Some(Self { a, b, c, d })
    }

    /// Creates a plane from a point and normal direction.
    ///
    /// The normal is automatically normalized.
    ///
    /// # Panics
    ///
    /// Panics if the normal vector has zero length.
    pub fn from_point_normal(point: &Point<T>, normal: &Vector<T>) -> Self {
        let len = normal.norm();
        assert!(
            len > T::epsilon(),
            "Cannot create plane from zero-length normal"
        );
        let (a, b, c) = (normal.x() / len, normal.y() / len, normal.z() / len);

        let d = -(a * point.x() + b * point.y() + c * point.z());
        Self { a, b, c, d }
    }

    /// The xy-plane (z = 0).
    #[inline]
    pub fn xy() -> Self {
        Self::from_normal_distance(T::zero(), T::zero(), T::one(), T::zero())
    }

    /// The xz-plane (y = 0).
    #[inline]
    pub fn xz() -> Self {
        Self::from_normal_distance(T::zero(), T::one(), T::zero(), T::zero())
    }

    /// The yz-plane (x = 0).
    #[inline]
    pub fn yz() -> Self {
        Self::from_normal_distance(T::one(), T::zero(), T::zero(), T::zero())
    }

    /// Returns the unit normal vector.
    #[inline]
    pub fn normal(&self) -> Vector<T> {
        Vector::new(self.a, self.b, self.c)
    }

    /// Returns the signed distance from origin to the plane.
    #[inline]
    pub fn distance_from_origin(&self) -> T {
        self.d
    }

    /// Signed distance from a point to the plane.
    ///
    /// - Positive: point is on the side the normal points to
    /// - Zero: point is on the plane
    /// - Negative: point is on the opposite side
    #[inline]
    pub fn signed_distance(&self, p: &Point<T>) -> T {
        self.a * p.x() + self.b * p.y() + self.c * p.z() + self.d
    }

    /// Returns true if the point lies on the plane.
    #[inline]
    pub fn contains(&self, p: &Point<T>, epsilon: T) -> bool {
        self.signed_distance(p).abs() < epsilon
    }

    /// Projects a point onto the plane.
    pub fn project(&self, p: &Point<T>) -> Point<T> {
        let dist = self.signed_distance(p);
        Point::new(
            p.x() - dist * self.a,
            p.y() - dist * self.b,
            p.z() - dist * self.c,
        )
    }

    /// Accessor for normal x-component.
    #[inline]
    pub fn a(&self) -> T {
        self.a
    }

    /// Accessor for normal y-component.
    #[inline]
    pub fn b(&self) -> T {
        self.b
    }

    /// Accessor for normal z-component.
    #[inline]
    pub fn c(&self) -> T {
        self.c
    }

    /// Accessor for distance term.
    #[inline]
    pub fn d(&self) -> T {
        self.d
    }
}

impl<T: Float> Default for Plane<T> {
    fn default() -> Self {
        Self::xy()
    }
}

// ============================================================================
// Circle type
// ============================================================================

/// A circle in 3D CGA (grade-3 trivector).
///
/// Represents a circle as the intersection of two spheres, or equivalently
/// as the set of points equidistant from a center in a given plane.
///
/// # Example
///
/// ```
/// use clifford::specialized::conformal::dim3::{Point, Circle};
///
/// // Circle from three points
/// let p1 = Point::<f64>::new(1.0, 0.0, 0.0);
/// let p2 = Point::new(0.0, 1.0, 0.0);
/// let p3 = Point::new(-1.0, 0.0, 0.0);
/// let circle = Circle::from_three_points(&p1, &p2, &p3).unwrap();
///
/// // Verify radius
/// assert!((circle.radius() - 1.0).abs() < 1e-10);
/// ```
///
/// # References
///
/// - <https://conformalgeometricalgebra.org/wiki/index.php?title=Circle>
#[derive(Clone, Copy, Debug, PartialEq)]
#[repr(C)]
pub struct Circle<T: Float> {
    // We store center, radius, and normal for simplicity.
    /// Center x-coordinate.
    cx: T,
    /// Center y-coordinate.
    cy: T,
    /// Center z-coordinate.
    cz: T,
    /// Radius.
    r: T,
    /// Normal x-component (normalized).
    nx: T,
    /// Normal y-component.
    ny: T,
    /// Normal z-component.
    nz: T,
}

impl<T: Float> Circle<T> {
    /// Creates a circle from center, radius, and normal direction.
    ///
    /// The normal defines the plane containing the circle and is automatically normalized.
    ///
    /// # Panics
    ///
    /// Panics if the normal vector has zero length.
    pub fn from_center_radius_normal(center: &Vector<T>, radius: T, normal: &Vector<T>) -> Self {
        let len = normal.norm();
        assert!(
            len > T::epsilon(),
            "Cannot create circle from zero-length normal"
        );
        let (nx, ny, nz) = (normal.x() / len, normal.y() / len, normal.z() / len);

        Self {
            cx: center.x(),
            cy: center.y(),
            cz: center.z(),
            r: radius,
            nx,
            ny,
            nz,
        }
    }

    /// Creates a circle from three points.
    ///
    /// The three points define a unique circle. Returns `None` if the points
    /// are collinear (no unique circle exists).
    ///
    /// # Algorithm
    ///
    /// Uses perpendicular bisector intersection to find the circumcenter.
    /// See `derivations/src/clifford_derivations/cga.py::derive_circumcenter()`.
    ///
    /// # Example
    ///
    /// ```
    /// use clifford::specialized::conformal::dim3::{Point, Circle};
    ///
    /// let p1 = Point::new(1.0, 0.0, 0.0);
    /// let p2 = Point::new(0.0, 1.0, 0.0);
    /// let p3 = Point::new(-1.0, 0.0, 0.0);
    ///
    /// let circle = Circle::from_three_points(&p1, &p2, &p3).unwrap();
    /// assert!(circle.contains(&p1, 1e-10));
    /// ```
    pub fn from_three_points(p1: &Point<T>, p2: &Point<T>, p3: &Point<T>) -> Option<Self> {
        // Compute plane normal via cross product of edge vectors
        let v1 = (p2.x() - p1.x(), p2.y() - p1.y(), p2.z() - p1.z());
        let v2 = (p3.x() - p1.x(), p3.y() - p1.y(), p3.z() - p1.z());

        let nx = v1.1 * v2.2 - v1.2 * v2.1;
        let ny = v1.2 * v2.0 - v1.0 * v2.2;
        let nz = v1.0 * v2.1 - v1.1 * v2.0;

        let len = (nx * nx + ny * ny + nz * nz).sqrt();

        // Collinear points: normal has zero length
        if len < T::epsilon() {
            return None;
        }

        let (nx, ny, nz) = (nx / len, ny / len, nz / len);

        // Find circumcenter using perpendicular bisector intersection
        // Algorithm from derivations/src/clifford_derivations/cga.py
        let ax = p1.x();
        let ay = p1.y();
        let az = p1.z();
        let bx = p2.x();
        let by = p2.y();
        let bz = p2.z();
        let cx = p3.x();
        let cy = p3.y();
        let cz = p3.z();

        // Midpoints
        let m1 = ((ax + bx) / T::TWO, (ay + by) / T::TWO, (az + bz) / T::TWO);

        // Edge directions
        let d1 = (bx - ax, by - ay, bz - az);

        // Perpendicular direction in the plane: p1 = d1 × n
        let perp1 = (
            d1.1 * nz - d1.2 * ny,
            d1.2 * nx - d1.0 * nz,
            d1.0 * ny - d1.1 * nx,
        );

        // Delta midpoint: dm = m2 - m1 = (C - B) / 2
        let dm = ((cx - bx) / T::TWO, (cy - by) / T::TWO, (cz - bz) / T::TWO);

        // Edge direction d2
        let d2 = (cx - ax, cy - ay, cz - az);

        // Perpendicular direction: p2 = d2 × n
        let perp2 = (
            d2.1 * nz - d2.2 * ny,
            d2.2 * nx - d2.0 * nz,
            d2.0 * ny - d2.1 * nx,
        );

        // Numerator: ((m2 - m1) × p2) · n
        let dm_cross_p2 = (
            dm.1 * perp2.2 - dm.2 * perp2.1,
            dm.2 * perp2.0 - dm.0 * perp2.2,
            dm.0 * perp2.1 - dm.1 * perp2.0,
        );
        let numer = dm_cross_p2.0 * nx + dm_cross_p2.1 * ny + dm_cross_p2.2 * nz;

        // Denominator: (p1 × p2) · n
        let p1_cross_p2 = (
            perp1.1 * perp2.2 - perp1.2 * perp2.1,
            perp1.2 * perp2.0 - perp1.0 * perp2.2,
            perp1.0 * perp2.1 - perp1.1 * perp2.0,
        );
        let denom = p1_cross_p2.0 * nx + p1_cross_p2.1 * ny + p1_cross_p2.2 * nz;

        // Degenerate case: perpendicular bisectors are parallel or nearly parallel
        // Use a relative threshold to avoid numerical instability
        let denom_threshold = T::epsilon() * (T::one() + numer.abs());
        if denom.abs() < denom_threshold {
            return None;
        }

        let t = numer / denom;

        // Circumcenter
        let center_x = m1.0 + t * perp1.0;
        let center_y = m1.1 + t * perp1.1;
        let center_z = m1.2 + t * perp1.2;

        // Radius
        let dx = center_x - ax;
        let dy = center_y - ay;
        let dz = center_z - az;
        let r = (dx * dx + dy * dy + dz * dz).sqrt();

        Some(Self {
            cx: center_x,
            cy: center_y,
            cz: center_z,
            r,
            nx,
            ny,
            nz,
        })
    }

    /// Returns the center as a Euclidean vector.
    #[inline]
    pub fn center(&self) -> Vector<T> {
        Vector::new(self.cx, self.cy, self.cz)
    }

    /// Returns the center as a conformal [`Point`].
    #[inline]
    pub fn center_point(&self) -> Point<T> {
        Point::new(self.cx, self.cy, self.cz)
    }

    /// Returns the radius.
    #[inline]
    pub fn radius(&self) -> T {
        self.r
    }

    /// Returns the normal direction of the plane containing the circle.
    #[inline]
    pub fn normal(&self) -> Vector<T> {
        Vector::new(self.nx, self.ny, self.nz)
    }

    /// Returns true if this is a real circle (positive radius).
    #[inline]
    pub fn is_real(&self, epsilon: T) -> bool {
        self.r > epsilon
    }

    /// Returns true if the point lies on the circle.
    pub fn contains(&self, p: &Point<T>, epsilon: T) -> bool {
        // Point must be in the plane and at distance r from center
        let normal = self.normal();
        let plane = Plane::from_point_normal(&self.center_point(), &normal);
        if !plane.contains(p, epsilon) {
            return false;
        }

        let dx = p.x() - self.cx;
        let dy = p.y() - self.cy;
        let dz = p.z() - self.cz;
        let dist = (dx * dx + dy * dy + dz * dz).sqrt();

        (dist - self.r).abs() < epsilon
    }

    /// Returns the plane containing this circle.
    pub fn carrier_plane(&self) -> Plane<T> {
        let normal = self.normal();
        Plane::from_point_normal(&self.center_point(), &normal)
    }

    /// Accessor for center x-coordinate.
    #[inline]
    pub fn cx(&self) -> T {
        self.cx
    }

    /// Accessor for center y-coordinate.
    #[inline]
    pub fn cy(&self) -> T {
        self.cy
    }

    /// Accessor for center z-coordinate.
    #[inline]
    pub fn cz(&self) -> T {
        self.cz
    }
}

impl<T: Float> Default for Circle<T> {
    /// Returns the unit circle in the xy-plane centered at the origin.
    fn default() -> Self {
        Self::from_center_radius_normal(
            &Vector::new(T::zero(), T::zero(), T::zero()),
            T::one(),
            &Vector::new(T::zero(), T::zero(), T::one()),
        )
    }
}

// ============================================================================
// Line type
// ============================================================================

/// A line in 3D CGA (flat circle - grade-3 trivector through infinity).
///
/// Represents a line as a point and direction, or as the intersection of two planes.
///
/// # Example
///
/// ```
/// use clifford::specialized::conformal::dim3::{Point, Line};
///
/// // Line through two points
/// let p1 = Point::<f64>::new(0.0, 0.0, 0.0);
/// let p2 = Point::new(1.0, 1.0, 1.0);
/// let line = Line::from_two_points(&p1, &p2);
///
/// // Both points lie on line
/// assert!(line.contains(&p1, 1e-10));
/// assert!(line.contains(&p2, 1e-10));
/// ```
///
/// # References
///
/// - <https://conformalgeometricalgebra.org/wiki/index.php?title=Line>
#[derive(Clone, Copy, Debug, PartialEq)]
#[repr(C)]
pub struct Line<T: Float> {
    // We store a point on the line and the direction.
    /// Point on line x-coordinate.
    px: T,
    /// Point on line y-coordinate.
    py: T,
    /// Point on line z-coordinate.
    pz: T,
    /// Direction x-component (normalized).
    dx: T,
    /// Direction y-component.
    dy: T,
    /// Direction z-component.
    dz: T,
}

impl<T: Float> Line<T> {
    /// Creates a line from a point and direction.
    ///
    /// The direction is automatically normalized.
    ///
    /// # Panics
    ///
    /// Panics if the direction vector has zero length.
    pub fn from_point_direction(point: &Point<T>, direction: &Vector<T>) -> Self {
        let len = direction.norm();
        assert!(
            len > T::epsilon(),
            "Cannot create line from zero-length direction"
        );
        let (dx, dy, dz) = (
            direction.x() / len,
            direction.y() / len,
            direction.z() / len,
        );

        Self {
            px: point.x(),
            py: point.y(),
            pz: point.z(),
            dx,
            dy,
            dz,
        }
    }

    /// Creates a line from two points.
    pub fn from_two_points(p1: &Point<T>, p2: &Point<T>) -> Self {
        let direction = Vector::new(p2.x() - p1.x(), p2.y() - p1.y(), p2.z() - p1.z());
        Self::from_point_direction(p1, &direction)
    }

    /// The x-axis.
    #[inline]
    pub fn x_axis() -> Self {
        Self::from_point_direction(&Point::origin(), &Vector::unit_x())
    }

    /// The y-axis.
    #[inline]
    pub fn y_axis() -> Self {
        Self::from_point_direction(&Point::origin(), &Vector::unit_y())
    }

    /// The z-axis.
    #[inline]
    pub fn z_axis() -> Self {
        Self::from_point_direction(&Point::origin(), &Vector::unit_z())
    }

    /// Returns the unit direction vector.
    #[inline]
    pub fn direction(&self) -> Vector<T> {
        Vector::new(self.dx, self.dy, self.dz)
    }

    /// Returns a point on the line (the stored reference point).
    #[inline]
    pub fn point_on_line(&self) -> Point<T> {
        Point::new(self.px, self.py, self.pz)
    }

    /// Returns true if the point lies on the line.
    #[inline]
    pub fn contains(&self, p: &Point<T>, epsilon: T) -> bool {
        self.distance_to_point(p) < epsilon
    }

    /// Distance from a point to the line.
    pub fn distance_to_point(&self, p: &Point<T>) -> T {
        // Vector from line point to p
        let v = (p.x() - self.px, p.y() - self.py, p.z() - self.pz);

        // Cross product v × d
        let cross = (
            v.1 * self.dz - v.2 * self.dy,
            v.2 * self.dx - v.0 * self.dz,
            v.0 * self.dy - v.1 * self.dx,
        );

        // Distance is |v × d| / |d|, but d is normalized so |d| = 1
        (cross.0 * cross.0 + cross.1 * cross.1 + cross.2 * cross.2).sqrt()
    }

    /// Closest point on line to given point.
    pub fn closest_point(&self, p: &Point<T>) -> Point<T> {
        // Project p onto the line
        let v = (p.x() - self.px, p.y() - self.py, p.z() - self.pz);

        // t = (v · d) where d is already normalized
        let t = v.0 * self.dx + v.1 * self.dy + v.2 * self.dz;

        Point::new(
            self.px + t * self.dx,
            self.py + t * self.dy,
            self.pz + t * self.dz,
        )
    }

    /// Meet of line and plane: point of intersection.
    ///
    /// Returns `None` if line is parallel to plane.
    pub fn meet_plane(&self, plane: &Plane<T>) -> Option<Point<T>> {
        // Line: P(t) = p + t*d
        // Plane: a*x + b*y + c*z + d = 0
        // Substitute: a*(px + t*dx) + b*(py + t*dy) + c*(pz + t*dz) + d = 0
        // Solve for t: t = -(a*px + b*py + c*pz + d) / (a*dx + b*dy + c*dz)

        let denom = plane.a() * self.dx + plane.b() * self.dy + plane.c() * self.dz;

        if denom.abs() < T::epsilon() {
            // Line is parallel to plane
            return None;
        }

        let numer = -(plane.a() * self.px + plane.b() * self.py + plane.c() * self.pz + plane.d());
        let t = numer / denom;

        Some(Point::new(
            self.px + t * self.dx,
            self.py + t * self.dy,
            self.pz + t * self.dz,
        ))
    }

    /// Accessor for reference point x-coordinate.
    #[inline]
    pub fn px(&self) -> T {
        self.px
    }

    /// Accessor for reference point y-coordinate.
    #[inline]
    pub fn py(&self) -> T {
        self.py
    }

    /// Accessor for reference point z-coordinate.
    #[inline]
    pub fn pz(&self) -> T {
        self.pz
    }
}

impl<T: Float> Default for Line<T> {
    /// Returns the x-axis (line through origin along x direction).
    fn default() -> Self {
        Self::x_axis()
    }
}

// ============================================================================
// Dipole type
// ============================================================================

/// A dipole (point pair) in 3D CGA (grade-2 bivector).
///
/// Represents two points as a single algebraic object.
/// `D = P₁ ∧ P₂`
///
/// # Example
///
/// ```
/// use clifford::specialized::conformal::dim3::{Point, Dipole};
///
/// let p1 = Point::<f64>::new(0.0, 0.0, 0.0);
/// let p2 = Point::new(1.0, 0.0, 0.0);
/// let dipole = Dipole::from_two_points(&p1, &p2);
///
/// // Verify separation
/// assert!((dipole.separation() - 1.0).abs() < 1e-10);
/// ```
///
/// # References
///
/// - <https://conformalgeometricalgebra.org/wiki/index.php?title=Dipole>
#[derive(Clone, Copy, Debug, PartialEq)]
#[repr(C)]
pub struct Dipole<T: Float> {
    // We store the two points directly for simplicity.
    /// First point x-coordinate.
    x1: T,
    /// First point y-coordinate.
    y1: T,
    /// First point z-coordinate.
    z1: T,
    /// Second point x-coordinate.
    x2: T,
    /// Second point y-coordinate.
    y2: T,
    /// Second point z-coordinate.
    z2: T,
}

impl<T: Float> Dipole<T> {
    /// Creates a dipole from two points.
    ///
    /// `D = P₁ ∧ P₂`
    #[inline]
    pub fn from_two_points(p1: &Point<T>, p2: &Point<T>) -> Self {
        Self {
            x1: p1.x(),
            y1: p1.y(),
            z1: p1.z(),
            x2: p2.x(),
            y2: p2.y(),
            z2: p2.z(),
        }
    }

    /// Extracts the two points from the dipole.
    #[inline]
    pub fn extract_points(&self) -> (Point<T>, Point<T>) {
        (
            Point::new(self.x1, self.y1, self.z1),
            Point::new(self.x2, self.y2, self.z2),
        )
    }

    /// Distance between the two points.
    pub fn separation(&self) -> T {
        let dx = self.x2 - self.x1;
        let dy = self.y2 - self.y1;
        let dz = self.z2 - self.z1;
        (dx * dx + dy * dy + dz * dz).sqrt()
    }

    /// Midpoint of the two points.
    #[inline]
    pub fn midpoint(&self) -> Point<T> {
        Point::new(
            (self.x1 + self.x2) / T::TWO,
            (self.y1 + self.y2) / T::TWO,
            (self.z1 + self.z2) / T::TWO,
        )
    }

    /// Returns true if this dipole is real (two distinct points).
    #[inline]
    pub fn is_real(&self, epsilon: T) -> bool {
        self.separation() > epsilon
    }

    /// Accessor for first point x-coordinate.
    #[inline]
    pub fn x1(&self) -> T {
        self.x1
    }

    /// Accessor for first point y-coordinate.
    #[inline]
    pub fn y1(&self) -> T {
        self.y1
    }

    /// Accessor for first point z-coordinate.
    #[inline]
    pub fn z1(&self) -> T {
        self.z1
    }

    /// Accessor for second point x-coordinate.
    #[inline]
    pub fn x2(&self) -> T {
        self.x2
    }

    /// Accessor for second point y-coordinate.
    #[inline]
    pub fn y2(&self) -> T {
        self.y2
    }

    /// Accessor for second point z-coordinate.
    #[inline]
    pub fn z2(&self) -> T {
        self.z2
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

impl<T: Float> AbsDiffEq for Plane<T> {
    type Epsilon = T;

    fn default_epsilon() -> Self::Epsilon {
        T::default_epsilon()
    }

    fn abs_diff_eq(&self, other: &Self, epsilon: Self::Epsilon) -> bool {
        T::abs_diff_eq(&self.a, &other.a, epsilon)
            && T::abs_diff_eq(&self.b, &other.b, epsilon)
            && T::abs_diff_eq(&self.c, &other.c, epsilon)
            && T::abs_diff_eq(&self.d, &other.d, epsilon)
    }
}

impl<T: Float> RelativeEq for Plane<T> {
    fn default_max_relative() -> Self::Epsilon {
        T::default_max_relative()
    }

    fn relative_eq(
        &self,
        other: &Self,
        epsilon: Self::Epsilon,
        max_relative: Self::Epsilon,
    ) -> bool {
        T::relative_eq(&self.a, &other.a, epsilon, max_relative)
            && T::relative_eq(&self.b, &other.b, epsilon, max_relative)
            && T::relative_eq(&self.c, &other.c, epsilon, max_relative)
            && T::relative_eq(&self.d, &other.d, epsilon, max_relative)
    }
}

impl<T: Float> UlpsEq for Plane<T> {
    fn default_max_ulps() -> u32 {
        T::default_max_ulps()
    }

    fn ulps_eq(&self, other: &Self, epsilon: Self::Epsilon, max_ulps: u32) -> bool {
        T::ulps_eq(&self.a, &other.a, epsilon, max_ulps)
            && T::ulps_eq(&self.b, &other.b, epsilon, max_ulps)
            && T::ulps_eq(&self.c, &other.c, epsilon, max_ulps)
            && T::ulps_eq(&self.d, &other.d, epsilon, max_ulps)
    }
}

impl<T: Float> AbsDiffEq for Circle<T> {
    type Epsilon = T;

    fn default_epsilon() -> Self::Epsilon {
        T::default_epsilon()
    }

    fn abs_diff_eq(&self, other: &Self, epsilon: Self::Epsilon) -> bool {
        T::abs_diff_eq(&self.cx, &other.cx, epsilon)
            && T::abs_diff_eq(&self.cy, &other.cy, epsilon)
            && T::abs_diff_eq(&self.cz, &other.cz, epsilon)
            && T::abs_diff_eq(&self.r, &other.r, epsilon)
            && T::abs_diff_eq(&self.nx, &other.nx, epsilon)
            && T::abs_diff_eq(&self.ny, &other.ny, epsilon)
            && T::abs_diff_eq(&self.nz, &other.nz, epsilon)
    }
}

impl<T: Float> RelativeEq for Circle<T> {
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
            && T::relative_eq(&self.r, &other.r, epsilon, max_relative)
            && T::relative_eq(&self.nx, &other.nx, epsilon, max_relative)
            && T::relative_eq(&self.ny, &other.ny, epsilon, max_relative)
            && T::relative_eq(&self.nz, &other.nz, epsilon, max_relative)
    }
}

impl<T: Float> UlpsEq for Circle<T> {
    fn default_max_ulps() -> u32 {
        T::default_max_ulps()
    }

    fn ulps_eq(&self, other: &Self, epsilon: Self::Epsilon, max_ulps: u32) -> bool {
        T::ulps_eq(&self.cx, &other.cx, epsilon, max_ulps)
            && T::ulps_eq(&self.cy, &other.cy, epsilon, max_ulps)
            && T::ulps_eq(&self.cz, &other.cz, epsilon, max_ulps)
            && T::ulps_eq(&self.r, &other.r, epsilon, max_ulps)
            && T::ulps_eq(&self.nx, &other.nx, epsilon, max_ulps)
            && T::ulps_eq(&self.ny, &other.ny, epsilon, max_ulps)
            && T::ulps_eq(&self.nz, &other.nz, epsilon, max_ulps)
    }
}

impl<T: Float> AbsDiffEq for Line<T> {
    type Epsilon = T;

    fn default_epsilon() -> Self::Epsilon {
        T::default_epsilon()
    }

    fn abs_diff_eq(&self, other: &Self, epsilon: Self::Epsilon) -> bool {
        T::abs_diff_eq(&self.px, &other.px, epsilon)
            && T::abs_diff_eq(&self.py, &other.py, epsilon)
            && T::abs_diff_eq(&self.pz, &other.pz, epsilon)
            && T::abs_diff_eq(&self.dx, &other.dx, epsilon)
            && T::abs_diff_eq(&self.dy, &other.dy, epsilon)
            && T::abs_diff_eq(&self.dz, &other.dz, epsilon)
    }
}

impl<T: Float> RelativeEq for Line<T> {
    fn default_max_relative() -> Self::Epsilon {
        T::default_max_relative()
    }

    fn relative_eq(
        &self,
        other: &Self,
        epsilon: Self::Epsilon,
        max_relative: Self::Epsilon,
    ) -> bool {
        T::relative_eq(&self.px, &other.px, epsilon, max_relative)
            && T::relative_eq(&self.py, &other.py, epsilon, max_relative)
            && T::relative_eq(&self.pz, &other.pz, epsilon, max_relative)
            && T::relative_eq(&self.dx, &other.dx, epsilon, max_relative)
            && T::relative_eq(&self.dy, &other.dy, epsilon, max_relative)
            && T::relative_eq(&self.dz, &other.dz, epsilon, max_relative)
    }
}

impl<T: Float> UlpsEq for Line<T> {
    fn default_max_ulps() -> u32 {
        T::default_max_ulps()
    }

    fn ulps_eq(&self, other: &Self, epsilon: Self::Epsilon, max_ulps: u32) -> bool {
        T::ulps_eq(&self.px, &other.px, epsilon, max_ulps)
            && T::ulps_eq(&self.py, &other.py, epsilon, max_ulps)
            && T::ulps_eq(&self.pz, &other.pz, epsilon, max_ulps)
            && T::ulps_eq(&self.dx, &other.dx, epsilon, max_ulps)
            && T::ulps_eq(&self.dy, &other.dy, epsilon, max_ulps)
            && T::ulps_eq(&self.dz, &other.dz, epsilon, max_ulps)
    }
}

impl<T: Float> AbsDiffEq for Dipole<T> {
    type Epsilon = T;

    fn default_epsilon() -> Self::Epsilon {
        T::default_epsilon()
    }

    fn abs_diff_eq(&self, other: &Self, epsilon: Self::Epsilon) -> bool {
        T::abs_diff_eq(&self.x1, &other.x1, epsilon)
            && T::abs_diff_eq(&self.y1, &other.y1, epsilon)
            && T::abs_diff_eq(&self.z1, &other.z1, epsilon)
            && T::abs_diff_eq(&self.x2, &other.x2, epsilon)
            && T::abs_diff_eq(&self.y2, &other.y2, epsilon)
            && T::abs_diff_eq(&self.z2, &other.z2, epsilon)
    }
}

impl<T: Float> RelativeEq for Dipole<T> {
    fn default_max_relative() -> Self::Epsilon {
        T::default_max_relative()
    }

    fn relative_eq(
        &self,
        other: &Self,
        epsilon: Self::Epsilon,
        max_relative: Self::Epsilon,
    ) -> bool {
        T::relative_eq(&self.x1, &other.x1, epsilon, max_relative)
            && T::relative_eq(&self.y1, &other.y1, epsilon, max_relative)
            && T::relative_eq(&self.z1, &other.z1, epsilon, max_relative)
            && T::relative_eq(&self.x2, &other.x2, epsilon, max_relative)
            && T::relative_eq(&self.y2, &other.y2, epsilon, max_relative)
            && T::relative_eq(&self.z2, &other.z2, epsilon, max_relative)
    }
}

impl<T: Float> UlpsEq for Dipole<T> {
    fn default_max_ulps() -> u32 {
        T::default_max_ulps()
    }

    fn ulps_eq(&self, other: &Self, epsilon: Self::Epsilon, max_ulps: u32) -> bool {
        T::ulps_eq(&self.x1, &other.x1, epsilon, max_ulps)
            && T::ulps_eq(&self.y1, &other.y1, epsilon, max_ulps)
            && T::ulps_eq(&self.z1, &other.z1, epsilon, max_ulps)
            && T::ulps_eq(&self.x2, &other.x2, epsilon, max_ulps)
            && T::ulps_eq(&self.y2, &other.y2, epsilon, max_ulps)
            && T::ulps_eq(&self.z2, &other.z2, epsilon, max_ulps)
    }
}

// ============================================================================
// Flat Point type
// ============================================================================

/// A flat point in 3D CGA (grade-2 bivector).
///
/// Represents a point without radius/weight information, analogous to
/// a point in Projective Geometric Algebra (PGA).
///
/// # Representation
///
/// A flat point has the form:
/// ```text
/// p = x·e₁∞ + y·e₂∞ + z·e₃∞ + w·e₀∞
/// ```
///
/// where each component contains the factor e∞ (point at infinity).
///
/// # Relationship to Round Points
///
/// - `FlatPoint = RoundPoint ∧ e∞`
/// - `RoundPoint = FlatPoint ⌋ e₀` (left contraction)
///
/// Flat points are useful for:
/// - Line representation (outer product of flat points)
/// - Plane representation (outer product of flat points)
/// - When only position matters, not conformal weight
///
/// # Example
///
/// ```
/// use clifford::specialized::conformal::dim3::{Point, FlatPoint};
///
/// // Create a flat point at (1, 2, 3)
/// let fp = FlatPoint::<f64>::new(1.0, 2.0, 3.0);
///
/// // Coordinates round-trip
/// assert!((fp.x() - 1.0).abs() < 1e-10);
/// assert!((fp.y() - 2.0).abs() < 1e-10);
/// assert!((fp.z() - 3.0).abs() < 1e-10);
///
/// // Convert to round point
/// let rp: Point<f64> = fp.to_round();
/// assert!((rp.x() - 1.0).abs() < 1e-10);
/// ```
///
/// # References
///
/// - <https://conformalgeometricalgebra.org/wiki/index.php?title=Flat_point>
#[derive(Clone, Copy, Debug, PartialEq)]
#[repr(C)]
pub struct FlatPoint<T: Float> {
    /// Coefficient of e₁∞ (x-coordinate times e∞).
    e1i: T,
    /// Coefficient of e₂∞ (y-coordinate times e∞).
    e2i: T,
    /// Coefficient of e₃∞ (z-coordinate times e∞).
    e3i: T,
    /// Coefficient of e₀∞ (weight times e∞).
    e0i: T,
}

impl<T: Float> FlatPoint<T> {
    /// Creates a flat point at Euclidean coordinates (x, y, z).
    ///
    /// Creates a unit-weight flat point where `e0i = 1`.
    ///
    /// # Example
    ///
    /// ```
    /// use clifford::specialized::conformal::dim3::FlatPoint;
    ///
    /// let fp = FlatPoint::<f64>::new(1.0, 2.0, 3.0);
    /// assert!((fp.x() - 1.0).abs() < 1e-10);
    /// ```
    #[inline]
    pub fn new(x: T, y: T, z: T) -> Self {
        // FlatPoint = RoundPoint ∧ e∞
        // For unit weight: e1i = x, e2i = y, e3i = z, e0i = 1
        Self {
            e1i: x,
            e2i: y,
            e3i: z,
            e0i: T::one(),
        }
    }

    /// Creates a flat point from raw components (unchecked).
    ///
    /// # Panics
    ///
    /// The accessor methods [`x()`](Self::x), [`y()`](Self::y), [`z()`](Self::z)
    /// will panic or return infinity if the point has zero weight (i.e., `e0i = 0`).
    /// Points created via [`FlatPoint::new`] always have unit weight.
    #[inline]
    pub fn from_components_unchecked(e1i: T, e2i: T, e3i: T, e0i: T) -> Self {
        Self { e1i, e2i, e3i, e0i }
    }

    /// Creates a flat point from a round point.
    ///
    /// `FlatPoint = RoundPoint ∧ e∞`
    ///
    /// For a round point P = x·e₁ + y·e₂ + z·e₃ + e₀ + ½r²·e∞,
    /// the outer product with e∞ gives a flat point with:
    /// - e1i = x (from x·e₁ ∧ e∞)
    /// - e2i = y (from y·e₂ ∧ e∞)
    /// - e3i = z (from z·e₃ ∧ e∞)
    /// - e0i = weight (from e₀ ∧ e∞)
    #[inline]
    pub fn from_round(p: &Point<T>) -> Self {
        // The round point already stores x, y, z directly in e1, e2, e3
        // The weight is em - ep (which equals 1 for unit-weight points)
        Self {
            e1i: p.e1(),
            e2i: p.e2(),
            e3i: p.e3(),
            e0i: p.weight(),
        }
    }

    /// The origin as a flat point.
    #[inline]
    pub fn origin() -> Self {
        Self::new(T::zero(), T::zero(), T::zero())
    }

    /// Extracts the x-coordinate.
    ///
    /// # Panics
    ///
    /// Panics or returns infinity if the flat point has zero weight.
    #[inline]
    pub fn x(&self) -> T {
        self.e1i / self.e0i
    }

    /// Extracts the y-coordinate.
    ///
    /// # Panics
    ///
    /// Panics or returns infinity if the flat point has zero weight.
    #[inline]
    pub fn y(&self) -> T {
        self.e2i / self.e0i
    }

    /// Extracts the z-coordinate.
    ///
    /// # Panics
    ///
    /// Panics or returns infinity if the flat point has zero weight.
    #[inline]
    pub fn z(&self) -> T {
        self.e3i / self.e0i
    }

    /// Returns the weight (e₀∞ component).
    ///
    /// For flat points created via [`FlatPoint::new`], this is always 1.
    #[inline]
    pub fn weight(&self) -> T {
        self.e0i
    }

    /// Normalizes the flat point so weight = 1.
    ///
    /// Returns `None` for ideal points (at infinity) where weight ≈ 0.
    pub fn normalize(&self) -> Option<Self> {
        if self.e0i.abs() < T::epsilon() {
            None // Ideal point (at infinity)
        } else {
            Some(Self {
                e1i: self.e1i / self.e0i,
                e2i: self.e2i / self.e0i,
                e3i: self.e3i / self.e0i,
                e0i: T::one(),
            })
        }
    }

    /// Returns true if this is an ideal flat point (at infinity).
    #[inline]
    pub fn is_ideal(&self, epsilon: T) -> bool {
        self.e0i.abs() < epsilon
    }

    /// Converts to a round point.
    ///
    /// `RoundPoint = FlatPoint ⌋ e₀`
    ///
    /// # Panics
    ///
    /// Panics if the flat point has zero weight.
    pub fn to_round(&self) -> Point<T> {
        let (x, y, z) = (self.x(), self.y(), self.z());
        Point::new(x, y, z)
    }

    /// Accessor for `e₁∞` component.
    #[inline]
    pub fn e1i(&self) -> T {
        self.e1i
    }

    /// Accessor for `e₂∞` component.
    #[inline]
    pub fn e2i(&self) -> T {
        self.e2i
    }

    /// Accessor for `e₃∞` component.
    #[inline]
    pub fn e3i(&self) -> T {
        self.e3i
    }

    /// Accessor for `e₀∞` component.
    #[inline]
    pub fn e0i(&self) -> T {
        self.e0i
    }
}

impl<T: Float> Default for FlatPoint<T> {
    fn default() -> Self {
        Self::origin()
    }
}

impl<T: Float> From<Point<T>> for FlatPoint<T> {
    fn from(p: Point<T>) -> Self {
        FlatPoint::from_round(&p)
    }
}

impl<T: Float> From<FlatPoint<T>> for Point<T> {
    fn from(fp: FlatPoint<T>) -> Self {
        fp.to_round()
    }
}

impl<T: Float> AbsDiffEq for FlatPoint<T> {
    type Epsilon = T;

    fn default_epsilon() -> Self::Epsilon {
        T::default_epsilon()
    }

    fn abs_diff_eq(&self, other: &Self, epsilon: Self::Epsilon) -> bool {
        T::abs_diff_eq(&self.e1i, &other.e1i, epsilon)
            && T::abs_diff_eq(&self.e2i, &other.e2i, epsilon)
            && T::abs_diff_eq(&self.e3i, &other.e3i, epsilon)
            && T::abs_diff_eq(&self.e0i, &other.e0i, epsilon)
    }
}

impl<T: Float> RelativeEq for FlatPoint<T> {
    fn default_max_relative() -> Self::Epsilon {
        T::default_max_relative()
    }

    fn relative_eq(
        &self,
        other: &Self,
        epsilon: Self::Epsilon,
        max_relative: Self::Epsilon,
    ) -> bool {
        T::relative_eq(&self.e1i, &other.e1i, epsilon, max_relative)
            && T::relative_eq(&self.e2i, &other.e2i, epsilon, max_relative)
            && T::relative_eq(&self.e3i, &other.e3i, epsilon, max_relative)
            && T::relative_eq(&self.e0i, &other.e0i, epsilon, max_relative)
    }
}

impl<T: Float> UlpsEq for FlatPoint<T> {
    fn default_max_ulps() -> u32 {
        T::default_max_ulps()
    }

    fn ulps_eq(&self, other: &Self, epsilon: Self::Epsilon, max_ulps: u32) -> bool {
        T::ulps_eq(&self.e1i, &other.e1i, epsilon, max_ulps)
            && T::ulps_eq(&self.e2i, &other.e2i, epsilon, max_ulps)
            && T::ulps_eq(&self.e3i, &other.e3i, epsilon, max_ulps)
            && T::ulps_eq(&self.e0i, &other.e0i, epsilon, max_ulps)
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

    // ========================================================================
    // Plane tests
    // ========================================================================

    proptest! {
        #[test]
        fn plane_contains_origin_when_d_is_zero(
            a in -1.0f64..1.0, b in -1.0f64..1.0, c in -1.0f64..1.0,
        ) {
            let len = (a*a + b*b + c*c).sqrt();
            if len < ABS_DIFF_EQ_EPS {
                return Ok(());
            }
            let plane = Plane::from_normal_distance(a/len, b/len, c/len, 0.0);
            let origin = Point::origin();
            prop_assert!(plane.contains(&origin, ABS_DIFF_EQ_EPS));
        }

        #[test]
        fn plane_distance_matches_euclidean(
            a in -1.0f64..1.0, b in -1.0f64..1.0, c in -1.0f64..1.0,
            d in -10.0f64..10.0,
            px in -10.0f64..10.0, py in -10.0f64..10.0, pz in -10.0f64..10.0,
        ) {
            let len = (a*a + b*b + c*c).sqrt();
            if len < ABS_DIFF_EQ_EPS {
                return Ok(());
            }
            let (a, b, c) = (a/len, b/len, c/len);

            let plane = Plane::from_normal_distance(a, b, c, d);
            let p = Point::new(px, py, pz);

            let cga_dist = plane.signed_distance(&p);
            let euclidean_dist = a*px + b*py + c*pz + d;

            prop_assert!(abs_diff_eq!(cga_dist, euclidean_dist, epsilon = ABS_DIFF_EQ_EPS));
        }

        #[test]
        fn plane_projection_on_plane(
            a in -1.0f64..1.0, b in -1.0f64..1.0, c in -1.0f64..1.0,
            d in -10.0f64..10.0,
            px in -10.0f64..10.0, py in -10.0f64..10.0, pz in -10.0f64..10.0,
        ) {
            let len = (a*a + b*b + c*c).sqrt();
            if len < ABS_DIFF_EQ_EPS {
                return Ok(());
            }
            let (a, b, c) = (a/len, b/len, c/len);

            let plane = Plane::from_normal_distance(a, b, c, d);
            let p = Point::new(px, py, pz);
            let projected = plane.project(&p);

            prop_assert!(plane.contains(&projected, ABS_DIFF_EQ_EPS));
        }
    }

    #[test]
    fn plane_from_three_points() {
        let p1 = Point::new(0.0, 0.0, 0.0);
        let p2 = Point::new(1.0, 0.0, 0.0);
        let p3 = Point::new(0.0, 1.0, 0.0);
        let plane = Plane::from_three_points(&p1, &p2, &p3).unwrap();

        // All three points should be on the plane
        assert!(plane.contains(&p1, ABS_DIFF_EQ_EPS));
        assert!(plane.contains(&p2, ABS_DIFF_EQ_EPS));
        assert!(plane.contains(&p3, ABS_DIFF_EQ_EPS));

        // Normal should be (0, 0, 1) or (0, 0, -1)
        let normal = plane.normal();
        assert!(abs_diff_eq!(normal.x(), 0.0, epsilon = ABS_DIFF_EQ_EPS));
        assert!(abs_diff_eq!(normal.y(), 0.0, epsilon = ABS_DIFF_EQ_EPS));
        assert!(abs_diff_eq!(
            normal.z().abs(),
            1.0,
            epsilon = ABS_DIFF_EQ_EPS
        ));
    }

    #[test]
    fn plane_from_three_points_collinear_returns_none() {
        // Collinear points should return None
        let p1 = Point::<f64>::new(0.0, 0.0, 0.0);
        let p2 = Point::new(1.0, 0.0, 0.0);
        let p3 = Point::new(2.0, 0.0, 0.0);
        assert!(Plane::from_three_points(&p1, &p2, &p3).is_none());
    }

    // ========================================================================
    // Line tests
    // ========================================================================

    proptest! {
        #[test]
        fn line_contains_defining_points(
            x1 in -10.0f64..10.0, y1 in -10.0f64..10.0, z1 in -10.0f64..10.0,
            x2 in -10.0f64..10.0, y2 in -10.0f64..10.0, z2 in -10.0f64..10.0,
        ) {
            // Skip degenerate case
            let dx = x2 - x1;
            let dy = y2 - y1;
            let dz = z2 - z1;
            let len = (dx*dx + dy*dy + dz*dz).sqrt();
            if len < ABS_DIFF_EQ_EPS {
                return Ok(());
            }

            let p1 = Point::new(x1, y1, z1);
            let p2 = Point::new(x2, y2, z2);
            let line = Line::from_two_points(&p1, &p2);

            prop_assert!(line.contains(&p1, ABS_DIFF_EQ_EPS));
            prop_assert!(line.contains(&p2, ABS_DIFF_EQ_EPS));
        }

        #[test]
        fn line_closest_point_is_perpendicular(
            x1 in -10.0f64..10.0, y1 in -10.0f64..10.0, z1 in -10.0f64..10.0,
            dx in -1.0f64..1.0, dy in -1.0f64..1.0, dz in -1.0f64..1.0,
            px in -10.0f64..10.0, py in -10.0f64..10.0, pz in -10.0f64..10.0,
        ) {
            let len = (dx*dx + dy*dy + dz*dz).sqrt();
            if len < ABS_DIFF_EQ_EPS {
                return Ok(());
            }

            let p1 = Point::new(x1, y1, z1);
            let direction = Vector::new(dx, dy, dz);
            let line = Line::from_point_direction(&p1, &direction);
            let p = Point::new(px, py, pz);
            let closest = line.closest_point(&p);

            // Vector from closest to p should be perpendicular to line direction
            let to_p = Vector::new(p.x() - closest.x(), p.y() - closest.y(), p.z() - closest.z());
            let dir = line.direction();
            let dot = to_p.x() * dir.x() + to_p.y() * dir.y() + to_p.z() * dir.z();

            prop_assert!(abs_diff_eq!(dot, 0.0, epsilon = ABS_DIFF_EQ_EPS));
        }

        #[test]
        fn line_distance_is_nonnegative(
            x1 in -10.0f64..10.0, y1 in -10.0f64..10.0, z1 in -10.0f64..10.0,
            dx in -1.0f64..1.0, dy in -1.0f64..1.0, dz in -1.0f64..1.0,
            px in -10.0f64..10.0, py in -10.0f64..10.0, pz in -10.0f64..10.0,
        ) {
            let len = (dx*dx + dy*dy + dz*dz).sqrt();
            if len < ABS_DIFF_EQ_EPS {
                return Ok(());
            }

            let p1 = Point::new(x1, y1, z1);
            let direction = Vector::new(dx, dy, dz);
            let line = Line::from_point_direction(&p1, &direction);
            let p = Point::new(px, py, pz);

            prop_assert!(line.distance_to_point(&p) >= 0.0);
        }
    }

    #[test]
    fn line_meet_plane_intersection() {
        // Line along z-axis
        let line = Line::x_axis();
        // Plane at x = 1
        let plane = Plane::from_normal_distance(1.0, 0.0, 0.0, -1.0);

        let intersection = line.meet_plane(&plane);
        assert!(intersection.is_some());
        let p = intersection.unwrap();
        assert!(abs_diff_eq!(p.x(), 1.0, epsilon = ABS_DIFF_EQ_EPS));
        assert!(abs_diff_eq!(p.y(), 0.0, epsilon = ABS_DIFF_EQ_EPS));
        assert!(abs_diff_eq!(p.z(), 0.0, epsilon = ABS_DIFF_EQ_EPS));
    }

    #[test]
    fn line_parallel_to_plane_no_intersection() {
        // Line along x-axis
        let line = Line::x_axis();
        // xy-plane (z = 0), parallel to x-axis and containing it
        // Actually x-axis lies ON the xy-plane, let's use a parallel plane
        let plane = Plane::from_normal_distance(0.0, 0.0, 1.0, -1.0); // z = 1 plane

        let intersection = line.meet_plane(&plane);
        assert!(intersection.is_none());
    }

    // ========================================================================
    // Circle tests
    // ========================================================================

    proptest! {
        #[test]
        fn circle_from_center_radius_normal_contains_points(
            cx in -10.0f64..10.0, cy in -10.0f64..10.0, cz in -10.0f64..10.0,
            r in 0.1f64..10.0,
            nx in -1.0f64..1.0, ny in -1.0f64..1.0, nz in -1.0f64..1.0,
            theta in 0.0f64..std::f64::consts::TAU,
        ) {
            // Skip degenerate normals
            let norm_sq = nx * nx + ny * ny + nz * nz;
            prop_assume!(norm_sq > 0.01);
            let norm = norm_sq.sqrt();
            let nx = nx / norm;
            let ny = ny / norm;
            let nz = nz / norm;

            let circle = Circle::from_center_radius_normal(
                &Vector::new(cx, cy, cz),
                r,
                &Vector::new(nx, ny, nz),
            );

            // Generate a point on the circle using angle theta
            // Need two perpendicular vectors in the circle's plane
            let (u, v) = {
                // Find a vector not parallel to normal
                let candidate = if nx.abs() < 0.9 {
                    Vector::new(1.0, 0.0, 0.0)
                } else {
                    Vector::new(0.0, 1.0, 0.0)
                };
                // u = n × candidate (normalized)
                let ux = ny * candidate.z() - nz * candidate.y();
                let uy = nz * candidate.x() - nx * candidate.z();
                let uz = nx * candidate.y() - ny * candidate.x();
                let u_len = (ux * ux + uy * uy + uz * uz).sqrt();
                let u = Vector::new(ux / u_len, uy / u_len, uz / u_len);
                // v = n × u
                let vx = ny * u.z() - nz * u.y();
                let vy = nz * u.x() - nx * u.z();
                let vz = nx * u.y() - ny * u.x();
                (u, Vector::new(vx, vy, vz))
            };

            // Point on circle: center + r * (cos(theta) * u + sin(theta) * v)
            let px = cx + r * (theta.cos() * u.x() + theta.sin() * v.x());
            let py = cy + r * (theta.cos() * u.y() + theta.sin() * v.y());
            let pz = cz + r * (theta.cos() * u.z() + theta.sin() * v.z());
            let p = Point::new(px, py, pz);

            prop_assert!(circle.contains(&p, 1e-6));
        }

        #[test]
        fn circle_radius_matches_distance_from_center(
            cx in -10.0f64..10.0, cy in -10.0f64..10.0, cz in -10.0f64..10.0,
            r in 0.1f64..10.0,
        ) {
            let circle = Circle::from_center_radius_normal(
                &Vector::new(cx, cy, cz),
                r,
                &Vector::new(0.0, 0.0, 1.0),
            );

            prop_assert!(abs_diff_eq!(circle.radius(), r, epsilon = ABS_DIFF_EQ_EPS));

            let center = circle.center();
            prop_assert!(abs_diff_eq!(center.x(), cx, epsilon = ABS_DIFF_EQ_EPS));
            prop_assert!(abs_diff_eq!(center.y(), cy, epsilon = ABS_DIFF_EQ_EPS));
            prop_assert!(abs_diff_eq!(center.z(), cz, epsilon = ABS_DIFF_EQ_EPS));
        }

        #[test]
        fn circle_carrier_plane_contains_circle(
            cx in -10.0f64..10.0, cy in -10.0f64..10.0, cz in -10.0f64..10.0,
            r in 0.1f64..10.0,
            nx in -1.0f64..1.0, ny in -1.0f64..1.0, nz in -1.0f64..1.0,
        ) {
            // Skip degenerate normals
            let norm_sq = nx * nx + ny * ny + nz * nz;
            prop_assume!(norm_sq > 0.01);

            let circle = Circle::from_center_radius_normal(
                &Vector::new(cx, cy, cz),
                r,
                &Vector::new(nx, ny, nz),
            );

            let plane = circle.carrier_plane();
            let center_point = Point::new(cx, cy, cz);
            prop_assert!(plane.contains(&center_point, 1e-6));
        }
    }

    #[test]
    fn circle_from_three_points_xy() {
        // Unit circle in xy-plane
        let p1 = Point::<f64>::new(1.0, 0.0, 0.0);
        let p2 = Point::new(0.0, 1.0, 0.0);
        let p3 = Point::new(-1.0, 0.0, 0.0);
        let circle = Circle::from_three_points(&p1, &p2, &p3).unwrap();

        // Center should be at origin
        assert!(abs_diff_eq!(circle.cx(), 0.0, epsilon = ABS_DIFF_EQ_EPS));
        assert!(abs_diff_eq!(circle.cy(), 0.0, epsilon = ABS_DIFF_EQ_EPS));
        assert!(abs_diff_eq!(circle.cz(), 0.0, epsilon = ABS_DIFF_EQ_EPS));

        // Radius should be 1
        assert!(abs_diff_eq!(
            circle.radius(),
            1.0,
            epsilon = ABS_DIFF_EQ_EPS
        ));

        // Normal should be (0, 0, ±1)
        let normal = circle.normal();
        assert!(abs_diff_eq!(normal.x(), 0.0, epsilon = ABS_DIFF_EQ_EPS));
        assert!(abs_diff_eq!(normal.y(), 0.0, epsilon = ABS_DIFF_EQ_EPS));
        assert!(abs_diff_eq!(
            normal.z().abs(),
            1.0,
            epsilon = ABS_DIFF_EQ_EPS
        ));
    }

    #[test]
    fn circle_from_three_points_collinear_returns_none() {
        // Collinear points should return None
        let p1 = Point::<f64>::new(0.0, 0.0, 0.0);
        let p2 = Point::new(1.0, 0.0, 0.0);
        let p3 = Point::new(2.0, 0.0, 0.0);
        assert!(Circle::from_three_points(&p1, &p2, &p3).is_none());
    }

    #[test]
    fn circle_contains_defining_points() {
        let p1 = Point::<f64>::new(1.0, 0.0, 0.0);
        let p2 = Point::new(0.0, 1.0, 0.0);
        let p3 = Point::new(-1.0, 0.0, 0.0);
        let circle = Circle::from_three_points(&p1, &p2, &p3).unwrap();

        assert!(circle.contains(&p1, ABS_DIFF_EQ_EPS));
        assert!(circle.contains(&p2, ABS_DIFF_EQ_EPS));
        assert!(circle.contains(&p3, ABS_DIFF_EQ_EPS));
    }

    // ========================================================================
    // Dipole tests
    // ========================================================================

    proptest! {
        #[test]
        fn dipole_separation_matches_distance(
            x1 in -10.0f64..10.0, y1 in -10.0f64..10.0, z1 in -10.0f64..10.0,
            x2 in -10.0f64..10.0, y2 in -10.0f64..10.0, z2 in -10.0f64..10.0,
        ) {
            let p1 = Point::new(x1, y1, z1);
            let p2 = Point::new(x2, y2, z2);
            let dipole = Dipole::from_two_points(&p1, &p2);

            let separation = dipole.separation();
            let distance = p1.distance(&p2);

            prop_assert!(abs_diff_eq!(separation, distance, epsilon = ABS_DIFF_EQ_EPS));
        }

        #[test]
        fn dipole_midpoint_is_equidistant(
            x1 in -10.0f64..10.0, y1 in -10.0f64..10.0, z1 in -10.0f64..10.0,
            x2 in -10.0f64..10.0, y2 in -10.0f64..10.0, z2 in -10.0f64..10.0,
        ) {
            let p1 = Point::new(x1, y1, z1);
            let p2 = Point::new(x2, y2, z2);
            let dipole = Dipole::from_two_points(&p1, &p2);
            let mid = dipole.midpoint();

            let d1 = p1.distance(&mid);
            let d2 = p2.distance(&mid);

            prop_assert!(abs_diff_eq!(d1, d2, epsilon = ABS_DIFF_EQ_EPS));
        }

        #[test]
        fn dipole_roundtrip(
            x1 in -10.0f64..10.0, y1 in -10.0f64..10.0, z1 in -10.0f64..10.0,
            x2 in -10.0f64..10.0, y2 in -10.0f64..10.0, z2 in -10.0f64..10.0,
        ) {
            let p1 = Point::new(x1, y1, z1);
            let p2 = Point::new(x2, y2, z2);
            let dipole = Dipole::from_two_points(&p1, &p2);
            let (q1, q2) = dipole.extract_points();

            // Points should match exactly (no need to check swapped since we store directly)
            prop_assert!(abs_diff_eq!(q1.x(), p1.x(), epsilon = ABS_DIFF_EQ_EPS));
            prop_assert!(abs_diff_eq!(q1.y(), p1.y(), epsilon = ABS_DIFF_EQ_EPS));
            prop_assert!(abs_diff_eq!(q1.z(), p1.z(), epsilon = ABS_DIFF_EQ_EPS));
            prop_assert!(abs_diff_eq!(q2.x(), p2.x(), epsilon = ABS_DIFF_EQ_EPS));
            prop_assert!(abs_diff_eq!(q2.y(), p2.y(), epsilon = ABS_DIFF_EQ_EPS));
            prop_assert!(abs_diff_eq!(q2.z(), p2.z(), epsilon = ABS_DIFF_EQ_EPS));
        }
    }

    #[test]
    fn dipole_is_real_when_distinct() {
        let p1 = Point::new(0.0, 0.0, 0.0);
        let p2 = Point::new(1.0, 0.0, 0.0);
        let dipole = Dipole::from_two_points(&p1, &p2);
        assert!(dipole.is_real(ABS_DIFF_EQ_EPS));
    }

    #[test]
    fn dipole_is_degenerate_when_same_point() {
        let p = Point::new(1.0, 2.0, 3.0);
        let dipole = Dipole::from_two_points(&p, &p);
        assert!(!dipole.is_real(ABS_DIFF_EQ_EPS));
    }

    // ========================================================================
    // FlatPoint tests
    // ========================================================================

    proptest! {
        #[test]
        fn flat_point_roundtrip_coordinates(
            x in -100.0f64..100.0,
            y in -100.0f64..100.0,
            z in -100.0f64..100.0,
        ) {
            let fp = FlatPoint::new(x, y, z);
            prop_assert!(abs_diff_eq!(fp.x(), x, epsilon = ABS_DIFF_EQ_EPS));
            prop_assert!(abs_diff_eq!(fp.y(), y, epsilon = ABS_DIFF_EQ_EPS));
            prop_assert!(abs_diff_eq!(fp.z(), z, epsilon = ABS_DIFF_EQ_EPS));
        }

        #[test]
        fn flat_point_unit_weight(
            x in -100.0f64..100.0,
            y in -100.0f64..100.0,
            z in -100.0f64..100.0,
        ) {
            let fp = FlatPoint::new(x, y, z);
            prop_assert!(abs_diff_eq!(fp.weight(), 1.0, epsilon = ABS_DIFF_EQ_EPS));
        }

        #[test]
        fn flat_round_roundtrip(
            x in -100.0f64..100.0,
            y in -100.0f64..100.0,
            z in -100.0f64..100.0,
        ) {
            let rp = Point::new(x, y, z);
            let fp = FlatPoint::from_round(&rp);
            let rp2 = fp.to_round();

            prop_assert!(abs_diff_eq!(rp.x(), rp2.x(), epsilon = ABS_DIFF_EQ_EPS));
            prop_assert!(abs_diff_eq!(rp.y(), rp2.y(), epsilon = ABS_DIFF_EQ_EPS));
            prop_assert!(abs_diff_eq!(rp.z(), rp2.z(), epsilon = ABS_DIFF_EQ_EPS));
        }

        #[test]
        fn round_flat_roundtrip(
            x in -100.0f64..100.0,
            y in -100.0f64..100.0,
            z in -100.0f64..100.0,
        ) {
            let fp = FlatPoint::new(x, y, z);
            let rp = fp.to_round();
            let fp2 = FlatPoint::from_round(&rp);

            prop_assert!(abs_diff_eq!(fp.x(), fp2.x(), epsilon = ABS_DIFF_EQ_EPS));
            prop_assert!(abs_diff_eq!(fp.y(), fp2.y(), epsilon = ABS_DIFF_EQ_EPS));
            prop_assert!(abs_diff_eq!(fp.z(), fp2.z(), epsilon = ABS_DIFF_EQ_EPS));
        }

        #[test]
        fn flat_point_normalize_preserves_coordinates(
            x in -100.0f64..100.0,
            y in -100.0f64..100.0,
            z in -100.0f64..100.0,
            scale in 0.1f64..10.0,
        ) {
            // Create a scaled flat point
            let fp = FlatPoint::from_components_unchecked(x * scale, y * scale, z * scale, scale);
            let normalized = fp.normalize().expect("non-zero weight");

            prop_assert!(abs_diff_eq!(normalized.x(), x, epsilon = ABS_DIFF_EQ_EPS));
            prop_assert!(abs_diff_eq!(normalized.y(), y, epsilon = ABS_DIFF_EQ_EPS));
            prop_assert!(abs_diff_eq!(normalized.z(), z, epsilon = ABS_DIFF_EQ_EPS));
            prop_assert!(abs_diff_eq!(normalized.weight(), 1.0, epsilon = ABS_DIFF_EQ_EPS));
        }
    }

    #[test]
    fn flat_point_origin() {
        let fp = FlatPoint::<f64>::origin();
        assert!(abs_diff_eq!(fp.x(), 0.0, epsilon = ABS_DIFF_EQ_EPS));
        assert!(abs_diff_eq!(fp.y(), 0.0, epsilon = ABS_DIFF_EQ_EPS));
        assert!(abs_diff_eq!(fp.z(), 0.0, epsilon = ABS_DIFF_EQ_EPS));
        assert!(abs_diff_eq!(fp.weight(), 1.0, epsilon = ABS_DIFF_EQ_EPS));
    }

    #[test]
    fn flat_point_is_ideal_at_infinity() {
        // An ideal flat point has zero weight
        let ideal = FlatPoint::from_components_unchecked(1.0, 0.0, 0.0, 0.0f64);
        assert!(ideal.is_ideal(ABS_DIFF_EQ_EPS));

        // A normal flat point is not ideal
        let normal = FlatPoint::new(1.0, 2.0, 3.0);
        assert!(!normal.is_ideal(ABS_DIFF_EQ_EPS));
    }

    #[test]
    fn flat_point_from_round_preserves_coordinates() {
        let rp = Point::<f64>::new(1.0, 2.0, 3.0);
        let fp = FlatPoint::from_round(&rp);

        assert!(abs_diff_eq!(fp.x(), 1.0, epsilon = ABS_DIFF_EQ_EPS));
        assert!(abs_diff_eq!(fp.y(), 2.0, epsilon = ABS_DIFF_EQ_EPS));
        assert!(abs_diff_eq!(fp.z(), 3.0, epsilon = ABS_DIFF_EQ_EPS));
    }

    #[test]
    fn flat_point_from_trait() {
        let rp = Point::<f64>::new(1.0, 2.0, 3.0);
        let fp: FlatPoint<f64> = rp.into();

        assert!(abs_diff_eq!(fp.x(), 1.0, epsilon = ABS_DIFF_EQ_EPS));
        assert!(abs_diff_eq!(fp.y(), 2.0, epsilon = ABS_DIFF_EQ_EPS));
        assert!(abs_diff_eq!(fp.z(), 3.0, epsilon = ABS_DIFF_EQ_EPS));
    }

    #[test]
    fn flat_point_to_round_trait() {
        let fp = FlatPoint::<f64>::new(1.0, 2.0, 3.0);
        let rp: Point<f64> = fp.into();

        assert!(abs_diff_eq!(rp.x(), 1.0, epsilon = ABS_DIFF_EQ_EPS));
        assert!(abs_diff_eq!(rp.y(), 2.0, epsilon = ABS_DIFF_EQ_EPS));
        assert!(abs_diff_eq!(rp.z(), 3.0, epsilon = ABS_DIFF_EQ_EPS));
    }

    // ========================================================================
    // Translator tests
    // ========================================================================

    #[test]
    fn translator_identity_is_noop() {
        proptest!(|(p in any::<Point<f64>>())| {
            let t = Translator::<f64>::identity();
            let q = t.transform_point(&p);

            prop_assert!(abs_diff_eq!(p.x(), q.x(), epsilon = ABS_DIFF_EQ_EPS));
            prop_assert!(abs_diff_eq!(p.y(), q.y(), epsilon = ABS_DIFF_EQ_EPS));
            prop_assert!(abs_diff_eq!(p.z(), q.z(), epsilon = ABS_DIFF_EQ_EPS));
        });
    }

    #[test]
    fn translator_inverse_roundtrip() {
        proptest!(|(t in any::<Translator<f64>>(), p in any::<Point<f64>>())| {
            let q = t.transform_point(&p);
            let r = t.inverse().transform_point(&q);

            prop_assert!(abs_diff_eq!(p.x(), r.x(), epsilon = ABS_DIFF_EQ_EPS));
            prop_assert!(abs_diff_eq!(p.y(), r.y(), epsilon = ABS_DIFF_EQ_EPS));
            prop_assert!(abs_diff_eq!(p.z(), r.z(), epsilon = ABS_DIFF_EQ_EPS));
        });
    }

    #[test]
    fn translator_compose_matches_sequential() {
        proptest!(|(t1 in any::<Translator<f64>>(), t2 in any::<Translator<f64>>(), p in any::<Point<f64>>())| {
            // Apply t1 then t2 sequentially
            let q_seq = t2.transform_point(&t1.transform_point(&p));

            // Apply composed translator
            let composed = t1.compose(&t2);
            let q_composed = composed.transform_point(&p);

            prop_assert!(abs_diff_eq!(q_seq.x(), q_composed.x(), epsilon = ABS_DIFF_EQ_EPS));
            prop_assert!(abs_diff_eq!(q_seq.y(), q_composed.y(), epsilon = ABS_DIFF_EQ_EPS));
            prop_assert!(abs_diff_eq!(q_seq.z(), q_composed.z(), epsilon = ABS_DIFF_EQ_EPS));
        });
    }

    #[test]
    fn translator_preserves_null() {
        proptest!(|(t in any::<Translator<f64>>(), p in any::<Point<f64>>())| {
            let q = t.transform_point(&p);
            prop_assert!(q.is_null(ABS_DIFF_EQ_EPS));
        });
    }

    // ========================================================================
    // Rotor tests
    // ========================================================================

    #[test]
    fn rotor_identity_is_noop() {
        proptest!(|(p in any::<Point<f64>>())| {
            let r = Rotor::<f64>::identity();
            let q = r.transform_point(&p);

            prop_assert!(abs_diff_eq!(p.x(), q.x(), epsilon = ABS_DIFF_EQ_EPS));
            prop_assert!(abs_diff_eq!(p.y(), q.y(), epsilon = ABS_DIFF_EQ_EPS));
            prop_assert!(abs_diff_eq!(p.z(), q.z(), epsilon = ABS_DIFF_EQ_EPS));
        });
    }

    #[test]
    fn rotor_preserves_distance_from_origin() {
        proptest!(|(r in any::<Rotor<f64>>(), p in any::<Point<f64>>())| {
            let q = r.transform_point(&p);

            let dist_p = p.distance(&Point::origin());
            let dist_q = q.distance(&Point::origin());

            prop_assert!(abs_diff_eq!(dist_p, dist_q, epsilon = ABS_DIFF_EQ_EPS));
        });
    }

    #[test]
    fn rotor_preserves_null() {
        proptest!(|(r in any::<Rotor<f64>>(), p in any::<Point<f64>>())| {
            let q = r.transform_point(&p);
            prop_assert!(q.is_null(ABS_DIFF_EQ_EPS));
        });
    }

    #[test]
    fn rotor_compose_matches_sequential() {
        proptest!(|(r1 in any::<Rotor<f64>>(), r2 in any::<Rotor<f64>>(), p in any::<Point<f64>>())| {
            // Apply r1 then r2 sequentially
            let q_seq = r2.transform_point(&r1.transform_point(&p));

            // Apply composed rotor
            let composed = r1.compose(&r2);
            let q_composed = composed.transform_point(&p);

            prop_assert!(abs_diff_eq!(q_seq.x(), q_composed.x(), epsilon = ABS_DIFF_EQ_EPS));
            prop_assert!(abs_diff_eq!(q_seq.y(), q_composed.y(), epsilon = ABS_DIFF_EQ_EPS));
            prop_assert!(abs_diff_eq!(q_seq.z(), q_composed.z(), epsilon = ABS_DIFF_EQ_EPS));
        });
    }

    #[test]
    fn rotor_inverse_roundtrip() {
        proptest!(|(r in any::<Rotor<f64>>(), p in any::<Point<f64>>())| {
            let q = r.transform_point(&p);
            let p_back = r.inverse().transform_point(&q);

            prop_assert!(abs_diff_eq!(p.x(), p_back.x(), epsilon = ABS_DIFF_EQ_EPS));
            prop_assert!(abs_diff_eq!(p.y(), p_back.y(), epsilon = ABS_DIFF_EQ_EPS));
            prop_assert!(abs_diff_eq!(p.z(), p_back.z(), epsilon = ABS_DIFF_EQ_EPS));
        });
    }

    // ========================================================================
    // Dilator tests
    // ========================================================================

    #[test]
    fn dilator_identity_is_noop() {
        proptest!(|(p in any::<Point<f64>>())| {
            let d = Dilator::<f64>::identity();
            let q = d.transform_point(&p);

            prop_assert!(abs_diff_eq!(p.x(), q.x(), epsilon = ABS_DIFF_EQ_EPS));
            prop_assert!(abs_diff_eq!(p.y(), q.y(), epsilon = ABS_DIFF_EQ_EPS));
            prop_assert!(abs_diff_eq!(p.z(), q.z(), epsilon = ABS_DIFF_EQ_EPS));
        });
    }

    #[test]
    fn dilator_center_is_invariant() {
        proptest!(|(d in any::<Dilator<f64>>())| {
            let center = d.center();
            let q = d.transform_point(&center);

            prop_assert!(abs_diff_eq!(center.x(), q.x(), epsilon = ABS_DIFF_EQ_EPS));
            prop_assert!(abs_diff_eq!(center.y(), q.y(), epsilon = ABS_DIFF_EQ_EPS));
            prop_assert!(abs_diff_eq!(center.z(), q.z(), epsilon = ABS_DIFF_EQ_EPS));
        });
    }

    #[test]
    fn dilator_inverse_roundtrip() {
        proptest!(|(d in any::<Dilator<f64>>(), p in any::<Point<f64>>())| {
            let q = d.transform_point(&p);
            let r = d.inverse().transform_point(&q);

            prop_assert!(abs_diff_eq!(p.x(), r.x(), epsilon = ABS_DIFF_EQ_EPS));
            prop_assert!(abs_diff_eq!(p.y(), r.y(), epsilon = ABS_DIFF_EQ_EPS));
            prop_assert!(abs_diff_eq!(p.z(), r.z(), epsilon = ABS_DIFF_EQ_EPS));
        });
    }

    #[test]
    fn dilator_preserves_null() {
        proptest!(|(d in any::<Dilator<f64>>(), p in any::<Point<f64>>())| {
            let q = d.transform_point(&p);
            prop_assert!(q.is_null(ABS_DIFF_EQ_EPS));
        });
    }

    #[test]
    fn dilator_compose_matches_sequential() {
        proptest!(|(d1 in any::<Dilator<f64>>(), d2 in any::<Dilator<f64>>(), p in any::<Point<f64>>())| {
            // For same-center dilators: d1.compose(d2) should equal sequential application
            // Create d2 with same center as d1
            let d2_same_center = Dilator::from_center_scale(&d1.center(), d2.scale());

            let q_seq = d2_same_center.transform_point(&d1.transform_point(&p));
            let q_composed = d1.compose(&d2_same_center).transform_point(&p);

            prop_assert!(abs_diff_eq!(q_seq.x(), q_composed.x(), epsilon = ABS_DIFF_EQ_EPS));
            prop_assert!(abs_diff_eq!(q_seq.y(), q_composed.y(), epsilon = ABS_DIFF_EQ_EPS));
            prop_assert!(abs_diff_eq!(q_seq.z(), q_composed.z(), epsilon = ABS_DIFF_EQ_EPS));
        });
    }
}

// ============================================================================
// Translator type
// ============================================================================

/// A translator (translation versor) in 3D CGA.
///
/// Represents a translation by a displacement vector `τ = (τₓ, τᵧ, τz)`.
///
/// The translator versor is:
/// ```text
/// T = 1 - (1/2)·τ·e∞
/// ```
///
/// It transforms objects via the sandwich product `T ⊛ x ⊛ T̃`.
///
/// # Example
///
/// ```
/// use clifford::specialized::conformal::dim3::{Point, Translator};
///
/// // Create a translator that moves by (1, 2, 3)
/// let t = Translator::new(1.0_f64, 2.0, 3.0);
///
/// // Translate a point at the origin
/// let p = Point::<f64>::origin();
/// let q = t.transform_point(&p);
///
/// assert!((q.x() - 1.0).abs() < 1e-10);
/// assert!((q.y() - 2.0).abs() < 1e-10);
/// assert!((q.z() - 3.0).abs() < 1e-10);
/// ```
///
/// # References
///
/// - <https://conformalgeometricalgebra.org/wiki/index.php?title=Translation>
#[derive(Clone, Copy, Debug, PartialEq)]
#[repr(C)]
pub struct Translator<T: Float> {
    /// Translation in x direction.
    tx: T,
    /// Translation in y direction.
    ty: T,
    /// Translation in z direction.
    tz: T,
}

impl<T: Float> Translator<T> {
    /// Creates a translator with the given displacement vector.
    ///
    /// The translator will move points by `(tx, ty, tz)`.
    ///
    /// # Example
    ///
    /// ```
    /// use clifford::specialized::conformal::dim3::Translator;
    ///
    /// let t = Translator::new(1.0_f64, 2.0, 3.0);
    /// assert!((t.tx() - 1.0).abs() < 1e-10);
    /// ```
    #[inline]
    pub fn new(tx: T, ty: T, tz: T) -> Self {
        Self { tx, ty, tz }
    }

    /// Creates a translator from a displacement vector.
    ///
    /// # Example
    ///
    /// ```
    /// use clifford::specialized::conformal::dim3::Translator;
    /// use clifford::specialized::euclidean::dim3::Vector;
    ///
    /// let v = Vector::new(1.0_f64, 2.0, 3.0);
    /// let t = Translator::from_vector(&v);
    /// assert!((t.tx() - 1.0).abs() < 1e-10);
    /// ```
    #[inline]
    pub fn from_vector(v: &Vector<T>) -> Self {
        Self::new(v.x(), v.y(), v.z())
    }

    /// Returns the identity translator (no translation).
    ///
    /// # Example
    ///
    /// ```
    /// use clifford::specialized::conformal::dim3::{Point, Translator};
    ///
    /// let t = Translator::<f64>::identity();
    /// let p = Point::new(1.0, 2.0, 3.0);
    /// let q = t.transform_point(&p);
    ///
    /// assert!((q.x() - p.x()).abs() < 1e-10);
    /// ```
    #[inline]
    pub fn identity() -> Self {
        Self::new(T::zero(), T::zero(), T::zero())
    }

    /// Returns the x component of the translation.
    #[inline]
    pub fn tx(&self) -> T {
        self.tx
    }

    /// Returns the y component of the translation.
    #[inline]
    pub fn ty(&self) -> T {
        self.ty
    }

    /// Returns the z component of the translation.
    #[inline]
    pub fn tz(&self) -> T {
        self.tz
    }

    /// Returns the translation as a Euclidean vector.
    #[inline]
    pub fn as_vector(&self) -> Vector<T> {
        Vector::new(self.tx, self.ty, self.tz)
    }

    /// Returns the inverse translator (translates in the opposite direction).
    ///
    /// For translator `T`, the inverse `T⁻¹` satisfies `T * T⁻¹ = 1`.
    ///
    /// # Example
    ///
    /// ```
    /// use clifford::specialized::conformal::dim3::{Point, Translator};
    ///
    /// let t = Translator::new(1.0_f64, 2.0, 3.0);
    /// let t_inv = t.inverse();
    ///
    /// // Composing with inverse gives identity
    /// let composed = t.compose(&t_inv);
    /// assert!((composed.tx()).abs() < 1e-10);
    /// ```
    #[inline]
    pub fn inverse(&self) -> Self {
        Self::new(-self.tx, -self.ty, -self.tz)
    }

    /// Composes two translators: `self` applied first, then `other`.
    ///
    /// `a.compose(&b).transform(x) == b.transform(a.transform(x))`
    ///
    /// Since translations commute, the order doesn't matter for translators,
    /// but this convention is consistent with other versors.
    ///
    /// # Example
    ///
    /// ```
    /// use clifford::specialized::conformal::dim3::Translator;
    ///
    /// let t1 = Translator::new(1.0_f64, 0.0, 0.0);
    /// let t2 = Translator::new(0.0_f64, 2.0, 0.0);
    /// let composed = t1.compose(&t2);
    ///
    /// assert!((composed.tx() - 1.0).abs() < 1e-10);
    /// assert!((composed.ty() - 2.0).abs() < 1e-10);
    /// ```
    #[inline]
    pub fn compose(&self, other: &Self) -> Self {
        Self::new(self.tx + other.tx, self.ty + other.ty, self.tz + other.tz)
    }

    /// Transforms a point by this translator.
    ///
    /// Applies the sandwich product `T * P * T̃` where `T̃` is the reverse.
    ///
    /// # Example
    ///
    /// ```
    /// use clifford::specialized::conformal::dim3::{Point, Translator};
    ///
    /// let t = Translator::new(1.0_f64, 2.0, 3.0);
    /// let p = Point::new(0.0_f64, 0.0, 0.0);
    /// let q = t.transform_point(&p);
    ///
    /// assert!((q.x() - 1.0).abs() < 1e-10);
    /// assert!((q.y() - 2.0).abs() < 1e-10);
    /// assert!((q.z() - 3.0).abs() < 1e-10);
    /// ```
    #[inline]
    pub fn transform_point(&self, p: &Point<T>) -> Point<T> {
        // Translation simply adds the displacement to the Euclidean coordinates
        Point::new(p.x() + self.tx, p.y() + self.ty, p.z() + self.tz)
    }

    /// Transforms a sphere by this translator.
    ///
    /// Translates the sphere's center while preserving its radius.
    ///
    /// # Example
    ///
    /// ```
    /// use clifford::specialized::conformal::dim3::{Sphere, Translator};
    ///
    /// let t = Translator::new(1.0_f64, 2.0, 3.0);
    /// let s = Sphere::from_center_radius(0.0_f64, 0.0, 0.0, 5.0);
    /// let s2 = t.transform_sphere(&s);
    ///
    /// assert!((s2.cx() - 1.0).abs() < 1e-10);
    /// assert!((s2.radius() - 5.0).abs() < 1e-10);
    /// ```
    #[inline]
    pub fn transform_sphere(&self, s: &Sphere<T>) -> Sphere<T> {
        Sphere::from_center_radius(
            s.cx() + self.tx,
            s.cy() + self.ty,
            s.cz() + self.tz,
            s.radius(),
        )
    }

    /// Transforms a plane by this translator.
    ///
    /// The normal direction is preserved; only the distance changes.
    ///
    /// # Example
    ///
    /// ```
    /// use clifford::specialized::conformal::dim3::{Plane, Translator};
    ///
    /// // xy-plane at z=0
    /// let plane = Plane::<f64>::xy();
    /// let t = Translator::new(0.0_f64, 0.0, 5.0);
    /// let translated = t.transform_plane(&plane);
    ///
    /// // Now z=5 plane
    /// assert!((translated.d() + 5.0).abs() < 1e-10);
    /// ```
    #[inline]
    pub fn transform_plane(&self, plane: &Plane<T>) -> Plane<T> {
        // d' = d - (a*tx + b*ty + c*tz)
        let new_d = plane.d() - (plane.a() * self.tx + plane.b() * self.ty + plane.c() * self.tz);
        Plane::from_normal_distance(plane.a(), plane.b(), plane.c(), new_d)
    }
}

impl<T: Float> Default for Translator<T> {
    /// Returns the identity translator.
    fn default() -> Self {
        Self::identity()
    }
}

// ============================================================================
// Translator approx traits
// ============================================================================

impl<T: Float> AbsDiffEq for Translator<T> {
    type Epsilon = T;

    fn default_epsilon() -> Self::Epsilon {
        T::epsilon()
    }

    fn abs_diff_eq(&self, other: &Self, epsilon: Self::Epsilon) -> bool {
        T::abs_diff_eq(&self.tx, &other.tx, epsilon)
            && T::abs_diff_eq(&self.ty, &other.ty, epsilon)
            && T::abs_diff_eq(&self.tz, &other.tz, epsilon)
    }
}

impl<T: Float> RelativeEq for Translator<T> {
    fn default_max_relative() -> Self::Epsilon {
        T::epsilon()
    }

    fn relative_eq(
        &self,
        other: &Self,
        epsilon: Self::Epsilon,
        max_relative: Self::Epsilon,
    ) -> bool {
        T::relative_eq(&self.tx, &other.tx, epsilon, max_relative)
            && T::relative_eq(&self.ty, &other.ty, epsilon, max_relative)
            && T::relative_eq(&self.tz, &other.tz, epsilon, max_relative)
    }
}

impl<T: Float> UlpsEq for Translator<T> {
    fn default_max_ulps() -> u32 {
        4
    }

    fn ulps_eq(&self, other: &Self, epsilon: Self::Epsilon, max_ulps: u32) -> bool {
        T::ulps_eq(&self.tx, &other.tx, epsilon, max_ulps)
            && T::ulps_eq(&self.ty, &other.ty, epsilon, max_ulps)
            && T::ulps_eq(&self.tz, &other.tz, epsilon, max_ulps)
    }
}

// ============================================================================
// Rotor type
// ============================================================================

use crate::specialized::euclidean::dim3::{Bivector, Rotor as EuclideanRotor};

/// A rotor (rotation versor) in 3D CGA.
///
/// Represents a rotation around an axis through the origin. The rotation
/// is encoded as a unit quaternion in geometric algebra form:
///
/// ```text
/// R = cos(θ/2) + sin(θ/2)·B
/// ```
///
/// where `B` is the unit bivector representing the rotation plane.
///
/// # Example
///
/// ```
/// use clifford::specialized::conformal::dim3::{Point, Rotor};
/// use clifford::specialized::euclidean::dim3::Bivector;
/// use std::f64::consts::FRAC_PI_2;
///
/// // Create a 90° rotation around the z-axis
/// let r = Rotor::from_angle_plane(FRAC_PI_2, Bivector::unit_xy());
///
/// // Rotate a point on the x-axis
/// let p = Point::new(1.0, 0.0, 0.0);
/// let q = r.transform_point(&p);
///
/// assert!((q.x()).abs() < 1e-10);
/// assert!((q.y() - 1.0).abs() < 1e-10);
/// ```
///
/// # References
///
/// - <https://conformalgeometricalgebra.org/wiki/index.php?title=Rotation>
#[derive(Clone, Copy, Debug, PartialEq)]
#[repr(transparent)]
pub struct Rotor<T: Float>(EuclideanRotor<T>);

impl<T: Float> Rotor<T> {
    /// Creates a rotor from scalar and bivector parts.
    ///
    /// The rotor should be normalized for proper rotation behavior.
    ///
    /// # Example
    ///
    /// ```
    /// use clifford::specialized::conformal::dim3::Rotor;
    /// use clifford::specialized::euclidean::dim3::Bivector;
    ///
    /// // Identity rotor
    /// let r = Rotor::new(1.0_f64, Bivector::zero());
    /// assert!((r.s() - 1.0).abs() < 1e-10);
    /// ```
    #[inline]
    pub fn new(s: T, b: Bivector<T>) -> Self {
        // Use new_unchecked since we're passing verified components
        Self(EuclideanRotor::new_unchecked(s, b.xy(), b.xz(), b.yz()))
    }

    /// Creates a rotor from an angle and rotation plane.
    ///
    /// The plane bivector should be normalized. The rotation follows the
    /// right-hand rule: positive angles rotate from the first basis vector
    /// toward the second.
    ///
    /// # Example
    ///
    /// ```
    /// use clifford::specialized::conformal::dim3::Rotor;
    /// use clifford::specialized::euclidean::dim3::Bivector;
    /// use std::f64::consts::FRAC_PI_4;
    ///
    /// // 45° rotation in the xy-plane (around z-axis)
    /// let r = Rotor::from_angle_plane(FRAC_PI_4, Bivector::unit_xy());
    /// ```
    #[inline]
    pub fn from_angle_plane(angle: T, plane: Bivector<T>) -> Self {
        Self(EuclideanRotor::from_angle_plane(angle, plane))
    }

    /// Creates a rotor from an axis and angle.
    ///
    /// The axis vector should be normalized. The rotation follows the
    /// right-hand rule.
    ///
    /// # Example
    ///
    /// ```
    /// use clifford::specialized::conformal::dim3::Rotor;
    /// use clifford::specialized::euclidean::dim3::Vector;
    /// use std::f64::consts::FRAC_PI_2;
    ///
    /// // 90° rotation around the z-axis
    /// let r = Rotor::from_axis_angle(&Vector::unit_z(), FRAC_PI_2);
    /// ```
    #[inline]
    pub fn from_axis_angle(axis: &Vector<T>, angle: T) -> Self {
        // The rotation plane bivector is the dual of the axis vector
        // For axis (x, y, z), the dual bivector is (yz, -xz, xy) = (x*e23, y*e31, z*e12)
        // But we need the bivector such that rotating in it gives rotation around axis
        // The bivector is: axis_x * e23 + axis_y * e31 + axis_z * e12
        let plane = Bivector::new(axis.z(), -axis.y(), axis.x());
        Self(EuclideanRotor::from_angle_plane(angle, plane))
    }

    /// Returns the identity rotor (no rotation).
    ///
    /// # Example
    ///
    /// ```
    /// use clifford::specialized::conformal::dim3::{Point, Rotor};
    ///
    /// let r = Rotor::<f64>::identity();
    /// let p = Point::new(1.0, 2.0, 3.0);
    /// let q = r.transform_point(&p);
    ///
    /// assert!((q.x() - p.x()).abs() < 1e-10);
    /// ```
    #[inline]
    pub fn identity() -> Self {
        Self(EuclideanRotor::identity())
    }

    /// Returns the scalar part (grade 0).
    #[inline]
    pub fn s(&self) -> T {
        self.0.s()
    }

    /// Returns the bivector part (grade 2).
    #[inline]
    pub fn b(&self) -> Bivector<T> {
        self.0.b()
    }

    /// Returns the underlying Euclidean rotor.
    #[inline]
    pub fn as_euclidean(&self) -> &EuclideanRotor<T> {
        &self.0
    }

    /// Returns the reverse of this rotor.
    ///
    /// For a unit rotor, the reverse equals the inverse.
    #[inline]
    pub fn reverse(&self) -> Self {
        Self(self.0.reverse())
    }

    /// Returns the inverse of this rotor.
    ///
    /// For a unit rotor, this equals the reverse.
    /// For non-unit rotors, this properly scales by the norm squared.
    #[inline]
    pub fn inverse(&self) -> Self {
        Self(self.0.inverse())
    }

    /// Composes two rotors: `self` applied first, then `other`.
    ///
    /// `a.compose(&b).transform(x) == b.transform(a.transform(x))`
    ///
    /// # Example
    ///
    /// ```
    /// use clifford::specialized::conformal::dim3::Rotor;
    /// use clifford::specialized::euclidean::dim3::Bivector;
    /// use std::f64::consts::FRAC_PI_4;
    ///
    /// let r1 = Rotor::from_angle_plane(FRAC_PI_4, Bivector::unit_xy());
    /// let r2 = Rotor::from_angle_plane(FRAC_PI_4, Bivector::unit_xy());
    /// let composed = r1.compose(&r2);
    ///
    /// // Composed rotation is 90°
    /// ```
    #[inline]
    pub fn compose(&self, other: &Self) -> Self {
        // R_composed = R_other * R_self
        Self(other.0.compose(self.0))
    }

    /// Transforms a point by this rotor.
    ///
    /// Applies the sandwich product `R * P * R̃`.
    ///
    /// # Example
    ///
    /// ```
    /// use clifford::specialized::conformal::dim3::{Point, Rotor};
    /// use clifford::specialized::euclidean::dim3::Bivector;
    /// use std::f64::consts::FRAC_PI_2;
    ///
    /// let r = Rotor::from_angle_plane(FRAC_PI_2, Bivector::unit_xy());
    /// let p = Point::new(1.0, 0.0, 0.0);
    /// let q = r.transform_point(&p);
    ///
    /// assert!((q.x()).abs() < 1e-10);
    /// assert!((q.y() - 1.0).abs() < 1e-10);
    /// assert!((q.z()).abs() < 1e-10);
    /// ```
    #[inline]
    pub fn transform_point(&self, p: &Point<T>) -> Point<T> {
        // For rotation around origin, only the Euclidean coordinates change
        let v = Vector::new(p.x(), p.y(), p.z());
        let rotated = self.0.rotate(v);
        Point::new(rotated.x(), rotated.y(), rotated.z())
    }

    /// Transforms a sphere by this rotor.
    ///
    /// Rotates the sphere's center around the origin.
    ///
    /// # Example
    ///
    /// ```
    /// use clifford::specialized::conformal::dim3::{Sphere, Rotor};
    /// use clifford::specialized::euclidean::dim3::Bivector;
    /// use std::f64::consts::FRAC_PI_2;
    ///
    /// let r = Rotor::from_angle_plane(FRAC_PI_2, Bivector::unit_xy());
    /// let s = Sphere::from_center_radius(1.0, 0.0, 0.0, 5.0);
    /// let s2 = r.transform_sphere(&s);
    ///
    /// assert!((s2.cx()).abs() < 1e-10);
    /// assert!((s2.cy() - 1.0).abs() < 1e-10);
    /// ```
    #[inline]
    pub fn transform_sphere(&self, s: &Sphere<T>) -> Sphere<T> {
        let center = Vector::new(s.cx(), s.cy(), s.cz());
        let rotated = self.0.rotate(center);
        Sphere::from_center_radius(rotated.x(), rotated.y(), rotated.z(), s.radius())
    }

    /// Transforms a plane by this rotor.
    ///
    /// Rotates the plane's normal around the origin.
    #[inline]
    pub fn transform_plane(&self, plane: &Plane<T>) -> Plane<T> {
        // Rotate the normal vector
        let normal = Vector::new(plane.a(), plane.b(), plane.c());
        let rotated_normal = self.0.rotate(normal);

        // For a plane through the origin (d=0), just use the rotated normal
        // For other planes, we need to find a point on the plane, rotate it,
        // and compute the new distance
        if plane.d().abs() < T::epsilon() {
            Plane::from_normal_distance(
                rotated_normal.x(),
                rotated_normal.y(),
                rotated_normal.z(),
                T::zero(),
            )
        } else {
            // Find a point on the plane: p = -d * n (assuming unit normal)
            let point_on_plane = Vector::new(
                -plane.d() * plane.a(),
                -plane.d() * plane.b(),
                -plane.d() * plane.c(),
            );
            let rotated_point = self.0.rotate(point_on_plane);
            // New distance is -dot(rotated_point, rotated_normal)
            let new_d = -(rotated_point.x() * rotated_normal.x()
                + rotated_point.y() * rotated_normal.y()
                + rotated_point.z() * rotated_normal.z());
            Plane::from_normal_distance(
                rotated_normal.x(),
                rotated_normal.y(),
                rotated_normal.z(),
                new_d,
            )
        }
    }
}

impl<T: Float> Default for Rotor<T> {
    /// Returns the identity rotor.
    fn default() -> Self {
        Self::identity()
    }
}

impl<T: Float> From<EuclideanRotor<T>> for Rotor<T> {
    fn from(r: EuclideanRotor<T>) -> Self {
        Self(r)
    }
}

impl<T: Float> From<Rotor<T>> for EuclideanRotor<T> {
    fn from(r: Rotor<T>) -> Self {
        r.0
    }
}

// ============================================================================
// Rotor approx traits
// ============================================================================

impl<T: Float> AbsDiffEq for Rotor<T> {
    type Epsilon = T;

    fn default_epsilon() -> Self::Epsilon {
        T::epsilon()
    }

    fn abs_diff_eq(&self, other: &Self, epsilon: Self::Epsilon) -> bool {
        self.0.abs_diff_eq(&other.0, epsilon)
    }
}

impl<T: Float> RelativeEq for Rotor<T> {
    fn default_max_relative() -> Self::Epsilon {
        T::epsilon()
    }

    fn relative_eq(
        &self,
        other: &Self,
        epsilon: Self::Epsilon,
        max_relative: Self::Epsilon,
    ) -> bool {
        self.0.relative_eq(&other.0, epsilon, max_relative)
    }
}

impl<T: Float> UlpsEq for Rotor<T> {
    fn default_max_ulps() -> u32 {
        4
    }

    fn ulps_eq(&self, other: &Self, epsilon: Self::Epsilon, max_ulps: u32) -> bool {
        self.0.ulps_eq(&other.0, epsilon, max_ulps)
    }
}

// ============================================================================
// Dilator type
// ============================================================================

/// A dilator (uniform scaling versor) in 3D CGA.
///
/// Represents a uniform scaling by factor `σ` centered at a point.
///
/// The dilator versor for origin-centered scaling is:
/// ```text
/// D = cosh(λ/2) + sinh(λ/2)·e∞∧e₀
/// ```
/// where `λ = ln(σ)` and `σ` is the scale factor.
///
/// # Example
///
/// ```
/// use clifford::specialized::conformal::dim3::{Point, Dilator};
///
/// // Create a dilator that scales by factor 2 around the origin
/// let d = Dilator::from_scale(2.0_f64);
///
/// // Scale a point
/// let p = Point::new(1.0_f64, 2.0, 3.0);
/// let q = d.transform_point(&p);
///
/// assert!((q.x() - 2.0).abs() < 1e-10);
/// assert!((q.y() - 4.0).abs() < 1e-10);
/// assert!((q.z() - 6.0).abs() < 1e-10);
/// ```
///
/// # References
///
/// - <https://conformalgeometricalgebra.org/wiki/index.php?title=Dilation>
#[derive(Clone, Copy, Debug, PartialEq)]
#[repr(C)]
pub struct Dilator<T: Float> {
    /// Scale factor (σ).
    scale: T,
    /// Center x-coordinate.
    cx: T,
    /// Center y-coordinate.
    cy: T,
    /// Center z-coordinate.
    cz: T,
}

impl<T: Float> Dilator<T> {
    /// Creates a dilator with the given scale factor centered at the origin.
    ///
    /// # Panics
    ///
    /// Panics if `scale` is zero or negative.
    ///
    /// # Example
    ///
    /// ```
    /// use clifford::specialized::conformal::dim3::Dilator;
    ///
    /// let d = Dilator::from_scale(2.0_f64);
    /// assert!((d.scale() - 2.0).abs() < 1e-10);
    /// ```
    #[inline]
    pub fn from_scale(scale: T) -> Self {
        assert!(scale > T::zero(), "Scale factor must be positive");
        Self {
            scale,
            cx: T::zero(),
            cy: T::zero(),
            cz: T::zero(),
        }
    }

    /// Creates a dilator with the given scale factor centered at the specified point.
    ///
    /// # Panics
    ///
    /// Panics if `scale` is zero or negative.
    ///
    /// # Example
    ///
    /// ```
    /// use clifford::specialized::conformal::dim3::{Dilator, Point};
    ///
    /// let center = Point::new(1.0_f64, 2.0, 3.0);
    /// let d = Dilator::from_center_scale(&center, 2.0);
    ///
    /// // The center point is invariant under dilation
    /// let q = d.transform_point(&center);
    /// assert!((q.x() - center.x()).abs() < 1e-10);
    /// ```
    #[inline]
    pub fn from_center_scale(center: &Point<T>, scale: T) -> Self {
        assert!(scale > T::zero(), "Scale factor must be positive");
        Self {
            scale,
            cx: center.x(),
            cy: center.y(),
            cz: center.z(),
        }
    }

    /// Returns the identity dilator (scale factor 1).
    ///
    /// # Example
    ///
    /// ```
    /// use clifford::specialized::conformal::dim3::{Point, Dilator};
    ///
    /// let d = Dilator::<f64>::identity();
    /// let p = Point::new(1.0, 2.0, 3.0);
    /// let q = d.transform_point(&p);
    ///
    /// assert!((q.x() - p.x()).abs() < 1e-10);
    /// ```
    #[inline]
    pub fn identity() -> Self {
        Self {
            scale: T::one(),
            cx: T::zero(),
            cy: T::zero(),
            cz: T::zero(),
        }
    }

    /// Returns the scale factor.
    #[inline]
    pub fn scale(&self) -> T {
        self.scale
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

    /// Returns the center as a point.
    #[inline]
    pub fn center(&self) -> Point<T> {
        Point::new(self.cx, self.cy, self.cz)
    }

    /// Returns the inverse dilator (scales by 1/σ).
    ///
    /// # Example
    ///
    /// ```
    /// use clifford::specialized::conformal::dim3::{Point, Dilator};
    ///
    /// let d = Dilator::from_scale(2.0_f64);
    /// let d_inv = d.inverse();
    ///
    /// // Composing with inverse gives identity
    /// let composed = d.compose(&d_inv);
    /// assert!((composed.scale() - 1.0).abs() < 1e-10);
    /// ```
    #[inline]
    pub fn inverse(&self) -> Self {
        Self {
            scale: T::one() / self.scale,
            cx: self.cx,
            cy: self.cy,
            cz: self.cz,
        }
    }

    /// Composes two dilators with the same center.
    ///
    /// # Panics
    ///
    /// Panics in debug builds if the dilators have different centers.
    ///
    /// # Note
    ///
    /// Only valid if both dilators have the same center.
    /// For dilators with different centers, use sequential transformation.
    ///
    /// # Example
    ///
    /// ```
    /// use clifford::specialized::conformal::dim3::Dilator;
    ///
    /// let d1 = Dilator::from_scale(2.0_f64);
    /// let d2 = Dilator::from_scale(3.0_f64);
    /// let composed = d1.compose(&d2);
    ///
    /// assert!((composed.scale() - 6.0).abs() < 1e-10);
    /// ```
    #[inline]
    pub fn compose(&self, other: &Self) -> Self {
        debug_assert!(
            T::abs_diff_eq(&self.cx, &other.cx, T::default_epsilon())
                && T::abs_diff_eq(&self.cy, &other.cy, T::default_epsilon())
                && T::abs_diff_eq(&self.cz, &other.cz, T::default_epsilon()),
            "Dilator::compose requires both dilators to have the same center"
        );
        // For same-center dilators, scales multiply
        Self {
            scale: self.scale * other.scale,
            cx: self.cx,
            cy: self.cy,
            cz: self.cz,
        }
    }

    /// Transforms a point by this dilator.
    ///
    /// Scales the point relative to the center by the scale factor.
    ///
    /// # Example
    ///
    /// ```
    /// use clifford::specialized::conformal::dim3::{Point, Dilator};
    ///
    /// let d = Dilator::from_scale(2.0_f64);
    /// let p = Point::new(1.0_f64, 2.0, 3.0);
    /// let q = d.transform_point(&p);
    ///
    /// assert!((q.x() - 2.0).abs() < 1e-10);
    /// assert!((q.y() - 4.0).abs() < 1e-10);
    /// assert!((q.z() - 6.0).abs() < 1e-10);
    /// ```
    #[inline]
    pub fn transform_point(&self, p: &Point<T>) -> Point<T> {
        // Scale relative to center:
        // q = center + scale * (p - center)
        //   = center * (1 - scale) + scale * p
        let x = self.cx * (T::one() - self.scale) + self.scale * p.x();
        let y = self.cy * (T::one() - self.scale) + self.scale * p.y();
        let z = self.cz * (T::one() - self.scale) + self.scale * p.z();
        Point::new(x, y, z)
    }

    /// Transforms a sphere by this dilator.
    ///
    /// Scales both the center position and the radius.
    ///
    /// # Example
    ///
    /// ```
    /// use clifford::specialized::conformal::dim3::{Sphere, Dilator};
    ///
    /// let d = Dilator::from_scale(2.0_f64);
    /// let s = Sphere::from_center_radius(1.0_f64, 0.0, 0.0, 1.0);
    /// let s2 = d.transform_sphere(&s);
    ///
    /// assert!((s2.cx() - 2.0).abs() < 1e-10);
    /// assert!((s2.radius() - 2.0).abs() < 1e-10);
    /// ```
    #[inline]
    pub fn transform_sphere(&self, s: &Sphere<T>) -> Sphere<T> {
        // Scale center relative to dilator center
        let cx = self.cx * (T::one() - self.scale) + self.scale * s.cx();
        let cy = self.cy * (T::one() - self.scale) + self.scale * s.cy();
        let cz = self.cz * (T::one() - self.scale) + self.scale * s.cz();
        // Radius scales uniformly
        let r = self.scale * s.radius();
        Sphere::from_center_radius(cx, cy, cz, r)
    }

    /// Transforms a plane by this dilator.
    ///
    /// For planes through the dilation center, the plane is unchanged.
    /// For other planes, the distance from center scales.
    ///
    /// # Example
    ///
    /// ```
    /// use clifford::specialized::conformal::dim3::{Plane, Dilator};
    ///
    /// let d = Dilator::from_scale(2.0_f64);
    /// let plane = Plane::from_normal_distance(0.0_f64, 0.0, 1.0, -1.0);  // z=1 plane
    /// let scaled = d.transform_plane(&plane);
    ///
    /// // Distance doubles
    /// assert!((scaled.d() + 2.0).abs() < 1e-10);
    /// ```
    #[inline]
    pub fn transform_plane(&self, plane: &Plane<T>) -> Plane<T> {
        // The normal direction is unchanged
        // The distance from center scales
        // d_new = scale * d + (1 - scale) * (n · center)
        let center_dot_n = plane.a() * self.cx + plane.b() * self.cy + plane.c() * self.cz;
        let new_d = self.scale * plane.d() + (T::one() - self.scale) * center_dot_n;
        Plane::from_normal_distance(plane.a(), plane.b(), plane.c(), new_d)
    }
}

impl<T: Float> Default for Dilator<T> {
    /// Returns the identity dilator.
    fn default() -> Self {
        Self::identity()
    }
}

// ============================================================================
// Dilator approx traits
// ============================================================================

impl<T: Float> AbsDiffEq for Dilator<T> {
    type Epsilon = T;

    fn default_epsilon() -> Self::Epsilon {
        T::epsilon()
    }

    fn abs_diff_eq(&self, other: &Self, epsilon: Self::Epsilon) -> bool {
        T::abs_diff_eq(&self.scale, &other.scale, epsilon)
            && T::abs_diff_eq(&self.cx, &other.cx, epsilon)
            && T::abs_diff_eq(&self.cy, &other.cy, epsilon)
            && T::abs_diff_eq(&self.cz, &other.cz, epsilon)
    }
}

impl<T: Float> RelativeEq for Dilator<T> {
    fn default_max_relative() -> Self::Epsilon {
        T::epsilon()
    }

    fn relative_eq(
        &self,
        other: &Self,
        epsilon: Self::Epsilon,
        max_relative: Self::Epsilon,
    ) -> bool {
        T::relative_eq(&self.scale, &other.scale, epsilon, max_relative)
            && T::relative_eq(&self.cx, &other.cx, epsilon, max_relative)
            && T::relative_eq(&self.cy, &other.cy, epsilon, max_relative)
            && T::relative_eq(&self.cz, &other.cz, epsilon, max_relative)
    }
}

impl<T: Float> UlpsEq for Dilator<T> {
    fn default_max_ulps() -> u32 {
        4
    }

    fn ulps_eq(&self, other: &Self, epsilon: Self::Epsilon, max_ulps: u32) -> bool {
        T::ulps_eq(&self.scale, &other.scale, epsilon, max_ulps)
            && T::ulps_eq(&self.cx, &other.cx, epsilon, max_ulps)
            && T::ulps_eq(&self.cy, &other.cy, epsilon, max_ulps)
            && T::ulps_eq(&self.cz, &other.cz, epsilon, max_ulps)
    }
}
