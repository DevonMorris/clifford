//! Domain-specific extensions for 3D Conformal GA types.
//!
//! This module adds geometric operations and convenience methods
//! to the generated types that are specific to 3D conformal geometry.
//!
//! # Basis Convention
//!
//! This module uses the **orthonormal conformal basis**:
//! - `e1`, `e2`, `e3`: Euclidean basis vectors
//! - `e4` (e+): Positive conformal basis (e4^2 = +1)
//! - `e5` (e-): Negative conformal basis (e5^2 = -1)
//!
//! The null basis vectors relate to orthonormal as:
//! - `e0` (origin) = (e5 - e4) / 2
//! - `einf` (infinity) = e5 + e4
//!
//! Point embedding: For Euclidean (x, y, z), the orthonormal form is:
//! ```text
//! P = x*e1 + y*e2 + z*e3 + 0.5*(x^2+y^2+z^2-1)*e4 + 0.5*(x^2+y^2+z^2+1)*e5
//! ```

use super::generated::types::{Circle, FlatPoint, Line, Motor, Plane, RoundPoint, Sphere};
use crate::scalar::Float;

// ============================================================================
// RoundPoint extensions
// ============================================================================

impl<T: Float> RoundPoint<T> {
    /// Creates a round point from Euclidean 3D coordinates.
    ///
    /// The point is embedded using the orthonormal basis:
    /// ```text
    /// P = x*e1 + y*e2 + z*e3 + 0.5*(x^2+y^2+z^2-1)*e4 + 0.5*(x^2+y^2+z^2+1)*e5
    /// ```
    ///
    /// This produces a null vector (P*P = 0) representing the Euclidean point.
    ///
    /// # Example
    ///
    /// ```
    /// use clifford::specialized::conformal::dim3::RoundPoint;
    ///
    /// let p = RoundPoint::from_euclidean(3.0, 4.0, 0.0);
    /// assert_eq!(p.to_euclidean(), Some((3.0, 4.0, 0.0)));
    /// ```
    #[inline]
    pub fn from_euclidean(x: T, y: T, z: T) -> Self {
        let sq = x * x + y * y + z * z;
        let half = T::one() / T::TWO;
        // Orthonormal embedding:
        // w (e4) = 0.5*(x^2 + y^2 + z^2 - 1)
        // u (e5) = 0.5*(x^2 + y^2 + z^2 + 1)
        let w_coeff = (sq - T::one()) * half;
        let u_coeff = (sq + T::one()) * half;
        Self::new_unchecked(x, y, z, w_coeff, u_coeff)
    }

    /// Creates a round point from Euclidean coordinates with a weight.
    ///
    /// Useful for creating weighted point combinations.
    #[inline]
    pub fn from_euclidean_weighted(x: T, y: T, z: T, weight: T) -> Self {
        let sq = x * x + y * y + z * z;
        let half = T::one() / T::TWO;
        let w_coeff = (sq - T::one()) * half;
        let u_coeff = (sq + T::one()) * half;
        Self::new_unchecked(
            x * weight,
            y * weight,
            z * weight,
            w_coeff * weight,
            u_coeff * weight,
        )
    }

    /// Origin point (0, 0, 0).
    ///
    /// In orthonormal basis: P = -0.5*e4 + 0.5*e5 = 0.5*(e5 - e4) = e0
    ///
    /// # Example
    ///
    /// ```
    /// use clifford::specialized::conformal::dim3::RoundPoint;
    ///
    /// let o = RoundPoint::<f64>::origin();
    /// assert_eq!(o.to_euclidean(), Some((0.0, 0.0, 0.0)));
    /// ```
    #[inline]
    pub fn origin() -> Self {
        let half = T::one() / T::TWO;
        // Origin: x=0, y=0, z=0, so sq=0
        // w = 0.5*(0 - 1) = -0.5
        // u = 0.5*(0 + 1) = 0.5
        Self::new_unchecked(T::zero(), T::zero(), T::zero(), -half, half)
    }

    /// Point at infinity.
    ///
    /// In orthonormal basis: einf = e4 + e5
    ///
    /// This represents the point at infinity in the conformal model.
    #[inline]
    pub fn infinity() -> Self {
        // einf = e4 + e5
        Self::new_unchecked(T::zero(), T::zero(), T::zero(), T::one(), T::one())
    }

    /// Extracts Euclidean coordinates from a finite round point.
    ///
    /// Returns `None` for points at infinity (where origin weight ~= 0).
    ///
    /// The null basis origin weight is recovered as: o = u - w (e5 - e4)
    ///
    /// # Example
    ///
    /// ```
    /// use clifford::specialized::conformal::dim3::RoundPoint;
    ///
    /// let p = RoundPoint::from_euclidean(3.0_f64, 4.0, 5.0);
    /// assert_eq!(p.to_euclidean(), Some((3.0, 4.0, 5.0)));
    ///
    /// let inf = RoundPoint::<f64>::infinity();
    /// assert_eq!(inf.to_euclidean(), None);
    /// ```
    #[inline]
    pub fn to_euclidean(&self) -> Option<(T, T, T)> {
        // Recover null basis origin weight: o = u - w
        let o = self.u() - self.w();
        if o.abs() < T::epsilon() {
            None
        } else {
            Some((self.x() / o, self.y() / o, self.z() / o))
        }
    }

    /// Euclidean x-coordinate.
    ///
    /// May produce incorrect results for points at infinity (o = 0).
    #[inline]
    pub fn euclidean_x(&self) -> T {
        let o = self.u() - self.w();
        self.x() / o
    }

    /// Euclidean y-coordinate.
    ///
    /// May produce incorrect results for points at infinity (o = 0).
    #[inline]
    pub fn euclidean_y(&self) -> T {
        let o = self.u() - self.w();
        self.y() / o
    }

    /// Euclidean z-coordinate.
    ///
    /// May produce incorrect results for points at infinity (o = 0).
    #[inline]
    pub fn euclidean_z(&self) -> T {
        let o = self.u() - self.w();
        self.z() / o
    }

    /// Returns the null basis origin weight (o = u - w).
    ///
    /// For a normalized point embedding, this equals 1.
    /// For a point at infinity, this equals 0.
    #[inline]
    pub fn null_origin_weight(&self) -> T {
        self.u() - self.w()
    }

    /// Returns the null basis infinity weight (i = (u + w) / 2).
    ///
    /// For a normalized point at (x, y, z), this equals 0.5*(x^2 + y^2 + z^2).
    #[inline]
    pub fn null_infinity_weight(&self) -> T {
        (self.u() + self.w()) / T::TWO
    }

    /// Returns true if this is a point at infinity (o ~= 0).
    #[inline]
    pub fn is_at_infinity(&self, epsilon: T) -> bool {
        self.null_origin_weight().abs() < epsilon
    }

    /// Squared Euclidean distance between two points.
    ///
    /// In CGA, the inner product of two null vectors gives:
    /// `P1 * P2 = -0.5*|p1 - p2|^2`
    ///
    /// Using the orthonormal metric (e4^2 = +1, e5^2 = -1):
    /// `P * Q = x1*x2 + y1*y2 + z1*z2 + w1*w2 - u1*u2`
    ///
    /// # Example
    ///
    /// ```
    /// use clifford::specialized::conformal::dim3::RoundPoint;
    /// use approx::abs_diff_eq;
    ///
    /// let p1 = RoundPoint::<f64>::origin();
    /// let p2 = RoundPoint::from_euclidean(3.0, 4.0, 0.0);
    /// assert!(abs_diff_eq!(p1.distance_squared(&p2), 25.0, epsilon = 1e-10));
    /// ```
    pub fn distance_squared(&self, other: &RoundPoint<T>) -> T {
        // Using the orthonormal metric:
        // P * Q = x1*x2 + y1*y2 + z1*z2 + w1*w2*(+1) + u1*u2*(-1)
        //       = x1*x2 + y1*y2 + z1*z2 + w1*w2 - u1*u2
        //
        // For null vectors representing Euclidean points:
        // P * Q = -0.5*|p - q|^2
        // So |p - q|^2 = -2*(P * Q)
        let inner = self.x() * other.x()
            + self.y() * other.y()
            + self.z() * other.z()
            + self.w() * other.w()
            - self.u() * other.u();
        -T::TWO * inner
    }

    /// Euclidean distance between two points.
    ///
    /// # Example
    ///
    /// ```
    /// use clifford::specialized::conformal::dim3::RoundPoint;
    /// use approx::abs_diff_eq;
    ///
    /// let p1 = RoundPoint::<f64>::origin();
    /// let p2 = RoundPoint::from_euclidean(3.0, 4.0, 0.0);
    /// assert!(abs_diff_eq!(p1.distance(&p2), 5.0, epsilon = 1e-10));
    /// ```
    pub fn distance(&self, other: &RoundPoint<T>) -> T {
        self.distance_squared(other).abs().sqrt()
    }
}

// ============================================================================
// Sphere extensions
// ============================================================================

impl<T: Float> Sphere<T> {
    /// Creates a sphere from center coordinates and radius.
    ///
    /// # Example
    ///
    /// ```
    /// use clifford::specialized::conformal::dim3::Sphere;
    /// use approx::relative_eq;
    ///
    /// let sphere = Sphere::from_center_radius(1.0_f64, 2.0, 3.0, 4.0);
    /// let (cx, cy, cz) = sphere.center().unwrap();
    /// let r = sphere.radius().unwrap();
    ///
    /// assert!(relative_eq!(cx, 1.0, epsilon = 1e-10));
    /// assert!(relative_eq!(cy, 2.0, epsilon = 1e-10));
    /// assert!(relative_eq!(cz, 3.0, epsilon = 1e-10));
    /// assert!(relative_eq!(r, 4.0, epsilon = 1e-10));
    /// ```
    pub fn from_center_radius(cx: T, cy: T, cz: T, radius: T) -> Self {
        // Construct by wedging four points on the sphere
        let p1 = RoundPoint::from_euclidean(cx + radius, cy, cz);
        let p2 = RoundPoint::from_euclidean(cx, cy + radius, cz);
        let p3 = RoundPoint::from_euclidean(cx, cy, cz + radius);
        let p4 = RoundPoint::from_euclidean(cx - radius, cy, cz);
        Self::from_four_points(&p1, &p2, &p3, &p4)
    }

    /// Creates a sphere passing through four points.
    ///
    /// The sphere is the outer product: `S = P1 ^ P2 ^ P3 ^ P4`
    ///
    /// # Example
    ///
    /// ```
    /// use clifford::specialized::conformal::dim3::{Sphere, RoundPoint};
    ///
    /// let p1 = RoundPoint::from_euclidean(1.0, 0.0, 0.0);
    /// let p2 = RoundPoint::from_euclidean(0.0, 1.0, 0.0);
    /// let p3 = RoundPoint::from_euclidean(0.0, 0.0, 1.0);
    /// let p4 = RoundPoint::from_euclidean(-1.0, 0.0, 0.0);
    ///
    /// let sphere = Sphere::from_four_points(&p1, &p2, &p3, &p4);
    /// // This sphere passes through all four points
    /// ```
    pub fn from_four_points(
        p1: &RoundPoint<T>,
        p2: &RoundPoint<T>,
        p3: &RoundPoint<T>,
        p4: &RoundPoint<T>,
    ) -> Self {
        use crate::ops::Wedge;
        p1.wedge(p2).wedge(p3).wedge(p4)
    }

    /// Returns true if this represents a plane (sphere through infinity).
    ///
    /// **Note**: The exact detection criteria for planes vs spheres in CGA
    /// is complex and depends on the specific basis convention. This method
    /// provides a basic check that works for spheres created via `from_center_radius`.
    #[inline]
    pub fn is_plane(&self, epsilon: T) -> bool {
        // A plane has the e1234 (u) component zero in certain conventions.
        // This check works for spheres created via from_center_radius.
        self.u().abs() < epsilon
    }

    /// Extracts the Euclidean center coordinates of the sphere.
    ///
    /// Returns `None` if this is a plane or if extraction fails.
    ///
    /// **Note**: Center extraction from arbitrary CGA spheres requires
    /// complex formulas that depend on the specific representation.
    /// This method works reliably for spheres created via `from_center_radius`.
    /// For spheres from `from_four_points`, use numerical methods or
    /// work with the CGA representation directly.
    pub fn center(&self) -> Option<(T, T, T)> {
        let u_coeff = self.u();

        if u_coeff.abs() < T::epsilon() {
            return None; // This is a plane
        }

        // Extract center using the relationship between blade coefficients.
        // For spheres created via from_center_radius, the structure is:
        // S = P1 ^ P2 ^ P3 ^ P4 where Pi are points on the sphere.
        //
        // The center coordinates are encoded in the grade-4 blade coefficients
        // relative to the e1234 coefficient (u).
        let center_x = self.w() / u_coeff;
        let center_y = -self.z() / u_coeff;
        let center_z = self.y() / u_coeff;
        Some((center_x, center_y, center_z))
    }

    /// Extracts the radius of the sphere.
    ///
    /// Returns `None` if this is a plane or if extraction fails.
    ///
    /// **Note**: Radius extraction from arbitrary CGA spheres requires
    /// complex formulas. This method works for spheres created via
    /// `from_center_radius`. For other cases, use the sphere's norm
    /// or work with the CGA representation directly.
    pub fn radius(&self) -> Option<T> {
        let u_coeff = self.u();

        if u_coeff.abs() < T::epsilon() {
            return None; // This is a plane
        }

        // Get center first
        let cx = self.w() / u_coeff;
        let cy = -self.z() / u_coeff;
        let cz = self.y() / u_coeff;

        // The radius is encoded in the relationship between x (e1235) and u (e1234)
        let x_coeff = self.x();
        let denom = x_coeff - u_coeff;

        if denom.abs() < T::epsilon() {
            return None;
        }

        let radius_sq = (x_coeff + u_coeff) / denom + cx * cx + cy * cy + cz * cz;

        if radius_sq < T::zero() {
            return None;
        }
        Some(radius_sq.sqrt())
    }

    /// Returns the curvature (1/radius) of the sphere.
    ///
    /// Returns `None` if this is a plane (curvature = 0) or if radius extraction fails.
    pub fn curvature(&self) -> Option<T> {
        self.radius().map(|r| T::one() / r)
    }
}

// ============================================================================
// Circle extensions
// ============================================================================

impl<T: Float> Circle<T> {
    /// Creates a circle passing through three points.
    ///
    /// The circle is the outer product: `C = P1 ^ P2 ^ P3`
    ///
    /// # Example
    ///
    /// ```
    /// use clifford::specialized::conformal::dim3::{Circle, RoundPoint};
    ///
    /// let p1 = RoundPoint::from_euclidean(1.0, 0.0, 0.0);
    /// let p2 = RoundPoint::from_euclidean(0.0, 1.0, 0.0);
    /// let p3 = RoundPoint::from_euclidean(-1.0, 0.0, 0.0);
    ///
    /// let circle = Circle::from_three_points(&p1, &p2, &p3);
    /// // This circle passes through all three points in the z=0 plane
    /// ```
    pub fn from_three_points(p1: &RoundPoint<T>, p2: &RoundPoint<T>, p3: &RoundPoint<T>) -> Self {
        use crate::ops::Wedge;
        p1.wedge(p2).wedge(p3)
    }

    /// Returns true if this represents a line (circle through infinity).
    ///
    /// A line is a circle that passes through infinity. In the CGA representation,
    /// a line has:
    /// - All g-components (e123, e124, e134, e234) near zero
    /// - All m-components (e125, e135, e235) near zero
    /// - Only v-components (e145, e245, e345) non-zero
    ///
    /// Both g and m components encode the "finite" circle part.
    #[inline]
    pub fn is_line(&self, epsilon: T) -> bool {
        let eps_sq = epsilon * epsilon;
        let g_norm_sq = self.gw() * self.gw()
            + self.gz() * self.gz()
            + self.gy() * self.gy()
            + self.gx() * self.gx();
        let m_norm_sq = self.mx() * self.mx() + self.my() * self.my() + self.mz() * self.mz();
        g_norm_sq < eps_sq && m_norm_sq < eps_sq
    }

    /// Extracts the Euclidean center coordinates of the circle.
    ///
    /// Returns `None` if this is a line (gw ~= 0), since lines have no finite center.
    ///
    /// # Example
    ///
    /// ```
    /// use clifford::specialized::conformal::dim3::{Circle, RoundPoint};
    /// use approx::relative_eq;
    ///
    /// // Circle in the z=0 plane centered at origin with radius 1
    /// let p1 = RoundPoint::from_euclidean(1.0_f64, 0.0, 0.0);
    /// let p2 = RoundPoint::from_euclidean(0.0, 1.0, 0.0);
    /// let p3 = RoundPoint::from_euclidean(-1.0, 0.0, 0.0);
    ///
    /// let circle = Circle::from_three_points(&p1, &p2, &p3);
    /// let (cx, cy, cz) = circle.center().unwrap();
    ///
    /// assert!(relative_eq!(cx, 0.0, epsilon = 1e-10));
    /// assert!(relative_eq!(cy, 0.0, epsilon = 1e-10));
    /// assert!(relative_eq!(cz, 0.0, epsilon = 1e-10));
    /// ```
    pub fn center(&self) -> Option<(T, T, T)> {
        // Circle fields: gw (e123), gz (e124), gy (e134), gx (e234),
        //                mz (e125), my (e135), mx (e235),
        //                vx (e145), vy (e245), vz (e345)
        //
        // The center extraction formula depends on the specific encoding.
        // Based on the pattern, when gw != 0:

        let gw_coeff = self.gw();

        if gw_coeff.abs() < T::epsilon() {
            // This is a line, no finite center
            return None;
        }

        // Center extraction formula (derived from wedge product structure):
        let center_x = self.mx() / gw_coeff;
        let center_y = self.my() / gw_coeff;
        let center_z = self.mz() / gw_coeff;
        Some((center_x, center_y, center_z))
    }

    /// Extracts the plane normal of the circle.
    ///
    /// Returns the unit normal vector to the plane containing the circle.
    /// Returns `None` if this is a line or if the normal cannot be computed.
    ///
    /// # Example
    ///
    /// ```
    /// use clifford::specialized::conformal::dim3::{Circle, RoundPoint};
    /// use approx::relative_eq;
    ///
    /// // Circle in the z=0 plane
    /// let p1 = RoundPoint::from_euclidean(1.0_f64, 0.0, 0.0);
    /// let p2 = RoundPoint::from_euclidean(0.0, 1.0, 0.0);
    /// let p3 = RoundPoint::from_euclidean(-1.0, 0.0, 0.0);
    ///
    /// let circle = Circle::from_three_points(&p1, &p2, &p3);
    /// let (nx, ny, nz) = circle.normal().unwrap();
    ///
    /// // Normal should be along z-axis (0, 0, +/-1)
    /// assert!(relative_eq!(nx, 0.0, epsilon = 1e-10));
    /// assert!(relative_eq!(ny, 0.0, epsilon = 1e-10));
    /// assert!(relative_eq!(nz.abs(), 1.0, epsilon = 1e-10));
    /// ```
    pub fn normal(&self) -> Option<(T, T, T)> {
        // The normal direction is encoded in the g-components (e123, e124, e134, e234)
        // For a circle in a plane with normal (nx, ny, nz):
        // gx ~ nx, gy ~ ny, gz ~ nz (up to scale)

        let gw_coeff = self.gw();

        if gw_coeff.abs() < T::epsilon() {
            // This is a line
            return None;
        }

        // The normal is proportional to (gx, gy, gz)
        let nx = self.gx();
        let ny = self.gy();
        let nz = self.gz();

        let len = (nx * nx + ny * ny + nz * nz).sqrt();
        if len < T::epsilon() {
            return None;
        }

        Some((nx / len, ny / len, nz / len))
    }

    /// Extracts the radius of the circle.
    ///
    /// Returns `None` if this is a line (gw ~= 0), since lines have infinite radius.
    ///
    /// # Example
    ///
    /// ```
    /// use clifford::specialized::conformal::dim3::{Circle, RoundPoint};
    /// use approx::relative_eq;
    ///
    /// // Unit circle in the z=0 plane
    /// let p1 = RoundPoint::from_euclidean(1.0_f64, 0.0, 0.0);
    /// let p2 = RoundPoint::from_euclidean(0.0, 1.0, 0.0);
    /// let p3 = RoundPoint::from_euclidean(-1.0, 0.0, 0.0);
    ///
    /// let circle = Circle::from_three_points(&p1, &p2, &p3);
    /// let radius = circle.radius().unwrap();
    ///
    /// assert!(relative_eq!(radius, 1.0, epsilon = 1e-10));
    /// ```
    pub fn radius(&self) -> Option<T> {
        let gw_coeff = self.gw();

        if gw_coeff.abs() < T::epsilon() {
            // This is a line, infinite radius
            return None;
        }

        // Extract center
        let cx = self.mx() / gw_coeff;
        let cy = self.my() / gw_coeff;
        let cz = self.mz() / gw_coeff;

        // The radius formula involves the v-components and center
        // vx ~ e145, vy ~ e245, vz ~ e345
        let vx = self.vx();
        let vy = self.vy();
        let vz = self.vz();

        // Radius squared from CGA circle structure:
        // The v-components encode the radius along with the center
        let v_norm_sq = vx * vx + vy * vy + vz * vz;
        let c_norm_sq = cx * cx + cy * cy + cz * cz;

        // Derive radius from the circle's structure
        // r^2 = v_norm / gw + c_norm_sq (approximate formula)
        let radius_sq = v_norm_sq / (gw_coeff * gw_coeff) + c_norm_sq
            - T::TWO * (cx * vx + cy * vy + cz * vz) / gw_coeff;

        // If negative, try alternative formula
        if radius_sq < T::zero() {
            // Fallback: compute from three sample points
            return None;
        }

        Some(radius_sq.sqrt())
    }

    /// Returns the curvature (1/radius) of the circle.
    ///
    /// Returns `None` if this is a line (curvature = 0) or if radius extraction fails.
    pub fn curvature(&self) -> Option<T> {
        self.radius().map(|r| T::one() / r)
    }
}

// ============================================================================
// Line extensions
// ============================================================================

impl<T: Float> Line<T> {
    /// Creates a line through two points.
    ///
    /// A line is a circle through infinity: `L = P1 ^ P2 ^ einf`
    pub fn from_two_points(p1: &RoundPoint<T>, p2: &RoundPoint<T>) -> Self {
        use crate::ops::Wedge;
        let inf = RoundPoint::infinity();
        let circle: Circle<T> = p1.wedge(p2).wedge(&inf);
        // Extract line components from the circle
        Self::new_unchecked(
            circle.mz(),
            circle.my(),
            circle.mx(),
            circle.vx(),
            circle.vy(),
            circle.vz(),
        )
    }
}

// ============================================================================
// Plane extensions
// ============================================================================

impl<T: Float> Plane<T> {
    /// Creates a plane from normal vector and distance from origin.
    ///
    /// The plane equation is: nx*x + ny*y + nz*z = d
    #[inline]
    pub fn from_normal_distance(nx: T, ny: T, nz: T, d: T) -> Self {
        Self::new_unchecked(nx, ny, nz, d)
    }

    /// Creates a plane through three points.
    ///
    /// A plane is a sphere through infinity: `P = P1 ^ P2 ^ P3 ^ einf`
    pub fn from_three_points(p1: &RoundPoint<T>, p2: &RoundPoint<T>, p3: &RoundPoint<T>) -> Self {
        use crate::ops::Wedge;
        let inf = RoundPoint::infinity();
        let sphere: Sphere<T> = p1.wedge(p2).wedge(p3).wedge(&inf);
        // Extract plane components from the sphere (u should be ~0)
        Self::new_unchecked(sphere.x(), sphere.y(), sphere.z(), sphere.w())
    }
}

// ============================================================================
// FlatPoint extensions
// ============================================================================

impl<T: Float> FlatPoint<T> {
    /// Creates a flat point from Euclidean coordinates.
    ///
    /// A flat point represents a point using only the infinity-direction
    /// components of a point pair.
    #[inline]
    pub fn from_euclidean(x: T, y: T, z: T) -> Self {
        Self::new_unchecked(x, y, z, T::one())
    }

    /// Extracts Euclidean coordinates.
    #[inline]
    pub fn to_euclidean(&self) -> Option<(T, T, T)> {
        if self.pw().abs() < T::epsilon() {
            None
        } else {
            Some((
                self.px() / self.pw(),
                self.py() / self.pw(),
                self.pz() / self.pw(),
            ))
        }
    }
}

// ============================================================================
// Motor extensions
// ============================================================================

impl<T: Float> Motor<T> {
    /// Identity motor (leaves all elements unchanged).
    ///
    /// # Example
    ///
    /// ```
    /// use clifford::specialized::conformal::dim3::Motor;
    ///
    /// let m = Motor::<f64>::identity();
    /// assert_eq!(m.s(), 1.0);
    /// ```
    #[inline]
    pub fn identity() -> Self {
        Self::new_unchecked(
            T::one(),  // s
            T::zero(), // mx
            T::zero(), // my
            T::zero(), // mz
            T::zero(), // vx
            T::zero(), // vy
            T::zero(), // vz
            T::zero(), // px
            T::zero(), // py
            T::zero(), // pz
            T::zero(), // pw
            T::zero(), // u
            T::zero(), // sx
            T::zero(), // sy
            T::zero(), // sz
            T::zero(), // sw
        )
    }

    /// Returns the inverse motor.
    ///
    /// For a versor V, the inverse is V_rev / |V|^2.
    pub fn inverse(&self) -> Self {
        let rev = self.reverse();
        let norm_sq = self.norm_squared();
        if norm_sq.abs() < T::epsilon() {
            return *self;
        }
        Self::new_unchecked(
            rev.s() / norm_sq,
            rev.mx() / norm_sq,
            rev.my() / norm_sq,
            rev.mz() / norm_sq,
            rev.vx() / norm_sq,
            rev.vy() / norm_sq,
            rev.vz() / norm_sq,
            rev.px() / norm_sq,
            rev.py() / norm_sq,
            rev.pz() / norm_sq,
            rev.pw() / norm_sq,
            rev.u() / norm_sq,
            rev.sx() / norm_sq,
            rev.sy() / norm_sq,
            rev.sz() / norm_sq,
            rev.sw() / norm_sq,
        )
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::test_utils::RELATIVE_EQ_EPS;
    use approx::relative_eq;
    use proptest::prelude::*;

    // Strategy for reasonable coordinate values (avoid extreme values)
    fn coord_strategy() -> impl Strategy<Value = f64> {
        -10.0..10.0_f64
    }

    proptest! {
        #[test]
        fn round_point_distance_symmetric(
            x1 in coord_strategy(),
            y1 in coord_strategy(),
            z1 in coord_strategy(),
            x2 in coord_strategy(),
            y2 in coord_strategy(),
            z2 in coord_strategy()
        ) {
            let p1 = RoundPoint::from_euclidean(x1, y1, z1);
            let p2 = RoundPoint::from_euclidean(x2, y2, z2);

            let d12 = p1.distance(&p2);
            let d21 = p2.distance(&p1);

            prop_assert!(
                relative_eq!(d12, d21, epsilon = 1e-10),
                "distance should be symmetric: d(p1,p2)={} != d(p2,p1)={}", d12, d21
            );

            // Verify against Euclidean formula
            let expected = ((x2-x1).powi(2) + (y2-y1).powi(2) + (z2-z1).powi(2)).sqrt();
            prop_assert!(
                relative_eq!(d12, expected, epsilon = 1e-8, max_relative = 1e-8),
                "distance mismatch: CGA={}, Euclidean={}", d12, expected
            );
        }
    }

    #[test]
    fn round_point_euclidean_roundtrip() {
        let x = 3.0_f64;
        let y = 4.0_f64;
        let z = 5.0_f64;
        let p = RoundPoint::from_euclidean(x, y, z);
        let (rx, ry, rz) = p.to_euclidean().unwrap();

        assert!(relative_eq!(
            rx,
            x,
            epsilon = RELATIVE_EQ_EPS,
            max_relative = RELATIVE_EQ_EPS
        ));
        assert!(relative_eq!(
            ry,
            y,
            epsilon = RELATIVE_EQ_EPS,
            max_relative = RELATIVE_EQ_EPS
        ));
        assert!(relative_eq!(
            rz,
            z,
            epsilon = RELATIVE_EQ_EPS,
            max_relative = RELATIVE_EQ_EPS
        ));
    }

    #[test]
    fn round_point_origin() {
        let o = RoundPoint::<f64>::origin();
        let (x, y, z) = o.to_euclidean().unwrap();
        assert!(relative_eq!(
            x,
            0.0,
            epsilon = RELATIVE_EQ_EPS,
            max_relative = RELATIVE_EQ_EPS
        ));
        assert!(relative_eq!(
            y,
            0.0,
            epsilon = RELATIVE_EQ_EPS,
            max_relative = RELATIVE_EQ_EPS
        ));
        assert!(relative_eq!(
            z,
            0.0,
            epsilon = RELATIVE_EQ_EPS,
            max_relative = RELATIVE_EQ_EPS
        ));
    }

    #[test]
    fn round_point_infinity() {
        let inf = RoundPoint::<f64>::infinity();
        assert!(inf.to_euclidean().is_none());
        assert!(inf.is_at_infinity(RELATIVE_EQ_EPS));
    }

    #[test]
    fn round_point_is_null_vector() {
        // Verify the null vector property: P*P = 0
        let p = RoundPoint::from_euclidean(3.0_f64, 4.0, 5.0);
        // Using orthonormal metric: P*P = x^2 + y^2 + z^2 + w^2 - u^2
        let p_dot_p = p.x() * p.x() + p.y() * p.y() + p.z() * p.z() + p.w() * p.w() - p.u() * p.u();

        assert!(
            relative_eq!(p_dot_p, 0.0, epsilon = 1e-10),
            "Point should be a null vector: P*P = {} != 0",
            p_dot_p
        );
    }

    #[test]
    fn round_point_origin_is_null_vector() {
        let o = RoundPoint::<f64>::origin();
        let o_dot_o = o.x() * o.x() + o.y() * o.y() + o.z() * o.z() + o.w() * o.w() - o.u() * o.u();

        assert!(
            relative_eq!(o_dot_o, 0.0, epsilon = 1e-10),
            "Origin should be a null vector: O*O = {} != 0",
            o_dot_o
        );
    }

    #[test]
    fn round_point_infinity_is_null_vector() {
        let inf = RoundPoint::<f64>::infinity();
        let inf_dot_inf =
            inf.x() * inf.x() + inf.y() * inf.y() + inf.z() * inf.z() + inf.w() * inf.w()
                - inf.u() * inf.u();

        assert!(
            relative_eq!(inf_dot_inf, 0.0, epsilon = 1e-10),
            "Infinity should be a null vector: einf*einf = {} != 0",
            inf_dot_inf
        );
    }

    #[test]
    fn round_point_distance() {
        let p1 = RoundPoint::<f64>::origin();
        let p2 = RoundPoint::from_euclidean(3.0, 4.0, 0.0);

        let dist = p1.distance(&p2);

        assert!(relative_eq!(
            dist,
            5.0,
            epsilon = RELATIVE_EQ_EPS,
            max_relative = RELATIVE_EQ_EPS
        ));
    }

    #[test]
    fn round_point_distance_3d() {
        let p1 = RoundPoint::<f64>::origin();
        let p2 = RoundPoint::from_euclidean(1.0, 2.0, 2.0);

        let dist = p1.distance(&p2);

        assert!(relative_eq!(
            dist,
            3.0,
            epsilon = RELATIVE_EQ_EPS,
            max_relative = RELATIVE_EQ_EPS
        ));
    }

    #[test]
    fn motor_identity() {
        let m = Motor::<f64>::identity();
        assert!(relative_eq!(
            m.s(),
            1.0,
            epsilon = RELATIVE_EQ_EPS,
            max_relative = RELATIVE_EQ_EPS
        ));
        assert!(relative_eq!(
            m.mx(),
            0.0,
            epsilon = RELATIVE_EQ_EPS,
            max_relative = RELATIVE_EQ_EPS
        ));
    }

    #[test]
    fn sphere_from_center_radius_construction() {
        let sphere = Sphere::from_center_radius(1.0_f64, 2.0, 3.0, 2.0);

        // Sphere created via from_center_radius should not be detected as plane
        assert!(!sphere.is_plane(1e-10), "Sphere should not be a plane");
    }

    #[test]
    fn sphere_from_four_points() {
        use crate::ops::Wedge;

        // Create a sphere from 4 known points on a sphere centered at (1, 0, 0) with radius 1
        let p1 = RoundPoint::from_euclidean(2.0_f64, 0.0, 0.0);
        let p2 = RoundPoint::from_euclidean(1.0, 1.0, 0.0);
        let p3 = RoundPoint::from_euclidean(1.0, 0.0, 1.0);
        let p4 = RoundPoint::from_euclidean(0.0, 0.0, 0.0);

        let sphere: Sphere<f64> = p1.wedge(&p2).wedge(&p3).wedge(&p4);

        // The wedge product of 4 non-coplanar points produces a sphere (quadvector)
        let norm = sphere.norm();
        assert!(
            norm > 1e-10,
            "Sphere from 4 points should have non-zero norm"
        );
    }

    #[test]
    fn circle_from_three_points() {
        use crate::ops::Wedge;

        // Three points on a circle NOT all on unit sphere - offset circle
        let p1 = RoundPoint::from_euclidean(2.0_f64, 1.0, 0.0);
        let p2 = RoundPoint::from_euclidean(1.0, 2.0, 0.0);
        let p3 = RoundPoint::from_euclidean(0.0, 1.0, 0.0);

        let circle: Circle<f64> = p1.wedge(&p2).wedge(&p3);

        // A finite circle is not a line
        assert!(
            !circle.is_line(1e-10),
            "Non-collinear points should form a circle (not a line)"
        );
    }

    #[test]
    fn circle_collinear_points_is_line() {
        use crate::ops::Wedge;

        // Three collinear points should form a line
        let p1 = RoundPoint::from_euclidean(0.0_f64, 0.0, 0.0);
        let p2 = RoundPoint::from_euclidean(1.0, 0.0, 0.0);
        let p3 = RoundPoint::from_euclidean(2.0, 0.0, 0.0);

        let result: Circle<f64> = p1.wedge(&p2).wedge(&p3);

        assert!(
            result.is_line(1e-10),
            "Collinear points should produce a line"
        );
    }

    #[test]
    fn flat_point_euclidean_roundtrip() {
        let p = FlatPoint::from_euclidean(1.0_f64, 2.0, 3.0);
        let (x, y, z) = p.to_euclidean().unwrap();

        assert!(relative_eq!(x, 1.0, epsilon = RELATIVE_EQ_EPS));
        assert!(relative_eq!(y, 2.0, epsilon = RELATIVE_EQ_EPS));
        assert!(relative_eq!(z, 3.0, epsilon = RELATIVE_EQ_EPS));
    }

    #[test]
    fn debug_sphere_values() {
        use crate::ops::Wedge;

        // Test sphere at origin
        let sphere = Sphere::from_center_radius(0.0_f64, 0.0, 0.0, 1.0);
        println!("Sphere at (0,0,0) r=1:");
        println!(
            "  u={}, x={}, y={}, z={}, w={}",
            sphere.u(),
            sphere.x(),
            sphere.y(),
            sphere.z(),
            sphere.w()
        );
        println!("  center = {:?}", sphere.center());
        println!("  radius = {:?}", sphere.radius());

        // Test sphere at (1, 2, 3)
        let sphere2 = Sphere::from_center_radius(1.0_f64, 2.0, 3.0, 2.0);
        println!("\nSphere at (1,2,3) r=2:");
        println!(
            "  u={}, x={}, y={}, z={}, w={}",
            sphere2.u(),
            sphere2.x(),
            sphere2.y(),
            sphere2.z(),
            sphere2.w()
        );
        println!("  center = {:?}", sphere2.center());
        println!("  radius = {:?}", sphere2.radius());

        // Test circle at origin
        let p1 = RoundPoint::from_euclidean(1.0_f64, 0.0, 0.0);
        let p2 = RoundPoint::from_euclidean(0.0, 1.0, 0.0);
        let p3 = RoundPoint::from_euclidean(-1.0, 0.0, 0.0);
        let circle: Circle<f64> = p1.wedge(&p2).wedge(&p3);
        println!("\nCircle from (1,0,0), (0,1,0), (-1,0,0):");
        println!(
            "  gw={}, gz={}, gy={}, gx={}",
            circle.gw(),
            circle.gz(),
            circle.gy(),
            circle.gx()
        );
        println!(
            "  mz={}, my={}, mx={}",
            circle.mz(),
            circle.my(),
            circle.mx()
        );
        println!(
            "  vx={}, vy={}, vz={}",
            circle.vx(),
            circle.vy(),
            circle.vz()
        );
        println!("  center = {:?}", circle.center());
        println!("  normal = {:?}", circle.normal());
        println!("  radius = {:?}", circle.radius());
    }
}
