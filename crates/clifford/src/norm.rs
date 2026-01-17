//! Norm traits for geometric algebra types.
//!
//! This module provides a hierarchy of traits for computing norms and normalizing
//! elements across different geometric algebras. Different algebras have different
//! norm semantics:
//!
//! - **Euclidean**: Standard positive-definite norm
//! - **PGA**: Degenerate metric with bulk norm and weight norm
//! - **CGA**: Null vectors, normalization via e-infinity coefficient
//! - **Minkowski**: Indefinite metric (timelike/spacelike/lightlike)
//!
//! # Trait Hierarchy
//!
//! ```text
//!                    Normed
//!                   /      \
//!     DegenerateNormed    IndefiniteNormed
//!           |                    |
//!          PGA              Minkowski
//!
//!     ConformalNormed
//!           |
//!          CGA
//! ```
//!
//! # Example
//!
//! ```ignore
//! use clifford::norm::Normed;
//!
//! let v = Vector::new(3.0, 4.0, 0.0);
//! assert_eq!(v.norm(), 5.0);
//!
//! let unit = v.normalize();
//! assert!((unit.norm() - 1.0).abs() < 1e-10);
//! ```
//!
//! # References
//!
//! - [Rigid GA Wiki - Geometric norm](https://rigidgeometricalgebra.org/wiki/index.php?title=Geometric_norm)
//! - [Spacetime Algebra](https://en.wikipedia.org/wiki/Spacetime_algebra)
//! - [CGA Documentation](https://clifford.readthedocs.io/en/latest/tutorials/cga/index.html)

use crate::scalar::Float;
// Import num_traits::Float anonymously to bring trait methods into scope
// without conflicting with crate::scalar::Float
use num_traits::Float as _;
use num_traits::Zero as _;

// ============================================================================
// Base Normed Trait
// ============================================================================

/// Trait for types that have a well-defined norm.
///
/// In geometric algebra, the norm of an element `u` is computed from
/// `u * rev(u)` (geometric product with reverse). The scalar part gives
/// the squared norm, which may be negative for indefinite metrics.
///
/// This trait provides the foundation for normalization across all
/// algebra types, with specialized sub-traits for specific metrics.
///
/// # Mathematical Background
///
/// For a multivector `u`, the squared norm is defined as the scalar part
/// of `u * rev(u)`, where `rev(u)` is the reverse (reverses the order of
/// basis vectors in each blade).
///
/// For simple elements like vectors in Euclidean space:
/// - `v * rev(v) = v * v = |v|^2` (always positive)
///
/// For rotors (even-grade elements):
/// - `r * rev(r) = s^2 + B^2` where `s` is scalar, `B` is bivector part
///
/// # Example
///
/// ```ignore
/// use clifford::norm::Normed;
///
/// // Euclidean vector
/// let v = Vector::new(3.0, 4.0, 0.0);
/// assert_eq!(v.norm_squared(), 25.0);
/// assert_eq!(v.norm(), 5.0);
///
/// // Normalize to unit length
/// if let Some(unit) = v.try_normalize() {
///     assert!((unit.norm() - 1.0).abs() < 1e-10);
/// }
/// ```
pub trait Normed {
    /// The scalar type used for norm computations.
    type Scalar: Float;

    /// Returns the squared norm of this element.
    ///
    /// **Note**: For indefinite metrics (Minkowski, etc.), this can be negative.
    /// Use [`norm`](Self::norm) for the absolute magnitude.
    ///
    /// # Mathematical Definition
    ///
    /// `norm_squared(u) = scalar_part(u * rev(u))`
    fn norm_squared(&self) -> Self::Scalar;

    /// Returns the absolute norm of this element.
    ///
    /// Computed as `sqrt(|norm_squared()|)`. This is always non-negative,
    /// regardless of metric signature.
    ///
    /// # Example
    ///
    /// ```ignore
    /// let v = Vector::new(3.0, 4.0, 0.0);
    /// assert_eq!(v.norm(), 5.0);
    /// ```
    fn norm(&self) -> Self::Scalar {
        self.norm_squared().abs().sqrt()
    }

    /// Attempts to normalize this element to unit norm.
    ///
    /// Returns `None` if the norm is too small (would cause division by zero
    /// or numerical instability).
    ///
    /// # Example
    ///
    /// ```ignore
    /// let v = Vector::new(3.0, 4.0, 0.0);
    /// let unit = v.try_normalize().unwrap();
    /// assert!((unit.norm() - 1.0).abs() < 1e-10);
    ///
    /// let zero = Vector::new(0.0, 0.0, 0.0);
    /// assert!(zero.try_normalize().is_none());
    /// ```
    fn try_normalize(&self) -> Option<Self>
    where
        Self: Sized;

    /// Normalizes this element to unit norm.
    ///
    /// # Panics
    ///
    /// Panics if the norm is too small to normalize. Use [`try_normalize`](Self::try_normalize)
    /// for a non-panicking version.
    ///
    /// # Example
    ///
    /// ```ignore
    /// let v = Vector::new(3.0, 4.0, 0.0);
    /// let unit = v.normalize();
    /// assert!((unit.norm() - 1.0).abs() < 1e-10);
    /// ```
    fn normalize(&self) -> Self
    where
        Self: Sized,
    {
        self.try_normalize().expect("cannot normalize zero element")
    }

    /// Scales this element by a scalar factor.
    ///
    /// # Example
    ///
    /// ```ignore
    /// let v = Vector::new(1.0, 2.0, 3.0);
    /// let scaled = v.scale(2.0);
    /// assert_eq!(scaled.x(), 2.0);
    /// ```
    fn scale(&self, factor: Self::Scalar) -> Self
    where
        Self: Sized;

    /// Returns true if this element can be normalized (has non-zero norm).
    ///
    /// An element is normalizable if its norm is greater than machine epsilon.
    ///
    /// # Example
    ///
    /// ```ignore
    /// let v = Vector::new(1.0, 0.0, 0.0);
    /// assert!(v.is_normalizable());
    ///
    /// let zero = Vector::new(0.0, 0.0, 0.0);
    /// assert!(!zero.is_normalizable());
    /// ```
    fn is_normalizable(&self) -> bool
    where
        Self: Sized,
    {
        self.norm() > Self::Scalar::epsilon()
    }
}

// ============================================================================
// DegenerateNormed Trait (PGA)
// ============================================================================

/// Trait for types in algebras with degenerate metrics (e.g., PGA).
///
/// In Projective Geometric Algebra (PGA), the metric is degenerate because
/// one basis vector squares to zero (`e0^2 = 0`). This leads to two distinct
/// norms:
///
/// - **Bulk norm**: `||u||_bulk = sqrt(u . u)` - magnitude of non-degenerate part
/// - **Weight norm**: `||u||_weight = sqrt(u (.) u)` - magnitude via antidot product
///
/// The geometric norm combines these to give meaningful distances and angles.
///
/// # Unitization vs Normalization
///
/// In PGA, there are two distinct operations:
///
/// - **Normalize**: Divide by bulk norm, making `u . u = 1`
/// - **Unitize**: Divide by weight norm, making `u (.) u = 1`
///
/// For rigid body transformations (motors), bulk normalization is typically used.
/// For points and planes in homogeneous coordinates, weight normalization
/// (unitization) gives the standard form.
///
/// # Example
///
/// ```ignore
/// use clifford::norm::DegenerateNormed;
///
/// let motor = Motor::from_rotation_z(0.5);
///
/// // Bulk norm is the rotor part magnitude
/// assert!((motor.bulk_norm() - 1.0).abs() < 1e-10);
///
/// // Pure rotation has zero weight (no translation)
/// assert!((motor.weight_norm()).abs() < 1e-10);
/// ```
///
/// # Reference
///
/// [Rigid GA Wiki - Geometric norm](https://rigidgeometricalgebra.org/wiki/index.php?title=Geometric_norm)
pub trait DegenerateNormed: Normed {
    /// Returns the squared bulk norm (from dot product).
    ///
    /// The bulk norm measures the magnitude of the non-degenerate part of the element.
    ///
    /// `||u||_bulk^2 = u . u`
    fn bulk_norm_squared(&self) -> Self::Scalar;

    /// Returns the bulk norm.
    ///
    /// `||u||_bulk = sqrt(|u . u|)`
    fn bulk_norm(&self) -> Self::Scalar {
        self.bulk_norm_squared().abs().sqrt()
    }

    /// Returns the squared weight norm (from antidot product).
    ///
    /// The weight norm measures the magnitude via the antidot product,
    /// which involves the degenerate (projective) part.
    ///
    /// `||u||_weight^2 = u (.) u`
    fn weight_norm_squared(&self) -> Self::Scalar;

    /// Returns the weight norm.
    ///
    /// `||u||_weight = sqrt(|u (.) u|)`
    fn weight_norm(&self) -> Self::Scalar {
        self.weight_norm_squared().abs().sqrt()
    }

    /// Attempts to unitize this element by dividing by its weight norm.
    ///
    /// Unitization is distinct from normalization:
    /// - **Normalize**: divide by bulk norm (makes `u . u = 1`)
    /// - **Unitize**: divide by weight norm (makes `u (.) u = 1`)
    ///
    /// Returns `None` if the weight norm is too small.
    ///
    /// # Example
    ///
    /// ```ignore
    /// // Unitize a point to standard homogeneous form (w = 1)
    /// let point = Point::new(2.0, 4.0, 6.0, 2.0);  // x=2, y=4, z=6, w=2
    /// let unitized = point.try_unitize().unwrap();
    /// // Now w = 1, giving standard form (1, 2, 3)
    /// ```
    fn try_unitize(&self) -> Option<Self>
    where
        Self: Sized;

    /// Unitizes this element by dividing by its weight norm.
    ///
    /// # Panics
    ///
    /// Panics if the weight norm is too small. Use [`try_unitize`](Self::try_unitize)
    /// for a non-panicking version.
    fn unitize(&self) -> Self
    where
        Self: Sized,
    {
        self.try_unitize()
            .expect("cannot unitize element with zero weight")
    }
}

// ============================================================================
// IndefiniteNormed Trait (Minkowski)
// ============================================================================

/// Trait for types in algebras with indefinite metrics (e.g., Minkowski spacetime).
///
/// In Minkowski spacetime with signature `(+,-,-,-)`, the squared norm of a
/// vector can be positive, negative, or zero:
///
/// - **Timelike**: `v^2 > 0` - inside the light cone
/// - **Spacelike**: `v^2 < 0` - outside the light cone
/// - **Lightlike/Null**: `v^2 = 0` - on the light cone
///
/// This classification is fundamental to special relativity, where timelike
/// vectors represent possible 4-velocities and lightlike vectors represent
/// the paths of light rays.
///
/// # Example
///
/// ```ignore
/// use clifford::norm::{IndefiniteNormed, CausalCharacter};
///
/// // With (+---) signature: t^2 - x^2 - y^2 - z^2
/// let timelike = FourVector::new(2.0, 1.0, 0.0, 0.0);  // 4 - 1 = 3 > 0
/// let spacelike = FourVector::new(1.0, 2.0, 0.0, 0.0); // 1 - 4 = -3 < 0
/// let lightlike = FourVector::new(1.0, 1.0, 0.0, 0.0); // 1 - 1 = 0
///
/// assert!(timelike.is_timelike());
/// assert!(spacelike.is_spacelike());
/// assert!(lightlike.is_lightlike());
///
/// assert_eq!(timelike.causal_character(), CausalCharacter::Timelike);
/// ```
///
/// # Reference
///
/// [Spacetime Algebra](https://en.wikipedia.org/wiki/Spacetime_algebra)
pub trait IndefiniteNormed: Normed {
    /// Returns true if this element is timelike (`norm_squared > 0`).
    ///
    /// In Minkowski spacetime, timelike vectors are inside the light cone
    /// and represent possible 4-velocities of massive particles.
    fn is_timelike(&self) -> bool {
        self.norm_squared() > Self::Scalar::zero()
    }

    /// Returns true if this element is spacelike (`norm_squared < 0`).
    ///
    /// In Minkowski spacetime, spacelike vectors are outside the light cone
    /// and represent spatial separations.
    fn is_spacelike(&self) -> bool {
        self.norm_squared() < Self::Scalar::zero()
    }

    /// Returns true if this element is lightlike/null (`norm_squared ~= 0`).
    ///
    /// In Minkowski spacetime, lightlike vectors are on the light cone
    /// and represent the paths of massless particles (photons).
    fn is_lightlike(&self) -> bool {
        self.norm_squared().abs() < Self::Scalar::epsilon()
    }

    /// Returns the causal character of this element.
    ///
    /// This classifies the element as timelike, spacelike, or lightlike
    /// based on the sign of its squared norm.
    ///
    /// # Example
    ///
    /// ```ignore
    /// let v = FourVector::new(2.0, 1.0, 0.0, 0.0);
    /// match v.causal_character() {
    ///     CausalCharacter::Timelike => println!("Inside light cone"),
    ///     CausalCharacter::Spacelike => println!("Outside light cone"),
    ///     CausalCharacter::Lightlike => println!("On light cone"),
    /// }
    /// ```
    fn causal_character(&self) -> CausalCharacter {
        let ns = self.norm_squared();
        if ns.abs() < Self::Scalar::epsilon() {
            CausalCharacter::Lightlike
        } else if ns > Self::Scalar::zero() {
            CausalCharacter::Timelike
        } else {
            CausalCharacter::Spacelike
        }
    }
}

/// Classification of vectors in indefinite-metric spaces.
///
/// In Minkowski spacetime, vectors are classified by the sign of their
/// squared norm relative to the light cone.
///
/// # Physical Interpretation
///
/// - **Timelike**: Represents possible world-lines of massive particles.
///   The proper time along a timelike path is real.
/// - **Spacelike**: Represents spatial separations. No massive particle
///   can travel along a purely spacelike path.
/// - **Lightlike**: Represents the world-lines of massless particles
///   (photons). Also called "null" vectors.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum CausalCharacter {
    /// Inside the light cone (`v^2 > 0` with `+---` signature).
    ///
    /// Timelike vectors represent possible 4-velocities of massive particles.
    Timelike,

    /// Outside the light cone (`v^2 < 0` with `+---` signature).
    ///
    /// Spacelike vectors represent spatial separations between events.
    Spacelike,

    /// On the light cone (`v^2 = 0`).
    ///
    /// Lightlike (or null) vectors represent the paths of massless particles.
    Lightlike,
}

// ============================================================================
// ConformalNormed Trait (CGA)
// ============================================================================

/// Trait for types in Conformal Geometric Algebra (CGA).
///
/// In CGA, geometric entities are represented using null vectors (vectors with
/// zero norm). The algebra adds two extra dimensions: `e+` (or `e_o`, the origin)
/// and `e-` (or `e_inf`, infinity), which together form a Minkowski plane.
///
/// Points in CGA are null vectors, and normalization typically means setting
/// the `e_inf` coefficient to 1 rather than making the Euclidean norm equal to 1.
///
/// # Null Vectors
///
/// In CGA, a point `P` with Euclidean coordinates `(x, y, z)` is represented as:
///
/// ```text
/// P = e_o + x*e1 + y*e2 + z*e3 + (x^2 + y^2 + z^2)/2 * e_inf
/// ```
///
/// This vector is null: `P . P = 0`.
///
/// # Normalization
///
/// CGA normalization sets `e_inf` coefficient to 1, giving standard form.
/// This is different from Euclidean normalization (length = 1).
///
/// # Example
///
/// ```ignore
/// use clifford::norm::ConformalNormed;
///
/// let point = CGAPoint::from_euclidean(1.0, 2.0, 3.0);
///
/// // Points are null vectors
/// assert!(point.is_null());
///
/// // Normalize to standard form (e_inf = 1)
/// let normalized = point.try_normalize_cga().unwrap();
/// assert!((normalized.einf_coefficient() - 1.0).abs() < 1e-10);
/// ```
///
/// # Reference
///
/// [CGA Documentation](https://clifford.readthedocs.io/en/latest/tutorials/cga/index.html)
pub trait ConformalNormed: Normed {
    /// Returns the `e_inf` (infinity) coefficient of this element.
    ///
    /// In CGA, the `e_inf` basis vector represents the point at infinity.
    /// For normalized points, this coefficient equals 1.
    fn einf_coefficient(&self) -> Self::Scalar;

    /// Returns true if this element represents a null vector (has zero norm).
    ///
    /// In CGA, geometric entities like points, circles, and spheres are
    /// represented by null vectors or null vector combinations.
    fn is_null(&self) -> bool {
        self.norm_squared().abs() < Self::Scalar::epsilon()
    }

    /// Normalizes by setting the `e_inf` coefficient to 1.
    ///
    /// This is the standard CGA normalization for points and other
    /// geometric entities with an `e_inf` component.
    ///
    /// Returns `None` if the `e_inf` coefficient is too small.
    ///
    /// # Example
    ///
    /// ```ignore
    /// let point = CGAPoint::new(/* ... */);
    /// if let Some(normalized) = point.try_normalize_cga() {
    ///     assert!((normalized.einf_coefficient() - 1.0).abs() < 1e-10);
    /// }
    /// ```
    fn try_normalize_cga(&self) -> Option<Self>
    where
        Self: Sized;
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn causal_character_equality() {
        assert_eq!(CausalCharacter::Timelike, CausalCharacter::Timelike);
        assert_ne!(CausalCharacter::Timelike, CausalCharacter::Spacelike);
        assert_ne!(CausalCharacter::Spacelike, CausalCharacter::Lightlike);
    }

    #[test]
    fn causal_character_debug() {
        assert_eq!(format!("{:?}", CausalCharacter::Timelike), "Timelike");
        assert_eq!(format!("{:?}", CausalCharacter::Spacelike), "Spacelike");
        assert_eq!(format!("{:?}", CausalCharacter::Lightlike), "Lightlike");
    }

    #[test]
    fn causal_character_clone() {
        let c = CausalCharacter::Timelike;
        #[allow(clippy::clone_on_copy)] // Testing Clone trait explicitly
        let cloned = c.clone();
        assert_eq!(c, cloned);
    }

    #[test]
    fn causal_character_copy() {
        let c = CausalCharacter::Spacelike;
        let copied: CausalCharacter = c;
        assert_eq!(c, copied);
    }
}
