//! Geometry-specific wrapper types for normalized and constrained entities.
//!
//! This module provides wrapper types that guarantee specific normalization
//! or constraint properties at the type level. There are two kinds of wrappers:
//!
//! - **Normalization wrappers**: Scale elements to satisfy a norm constraint
//! - **Constraint wrappers**: Verify elements satisfy a property (no scaling)
//!
//! # Wrapper Types
//!
//! ## Euclidean Algebras
//!
//! | Wrapper | Type | Constraint | Use Case |
//! |---------|------|------------|----------|
//! | [`Unit<T>`] | Normalization | `norm() == 1` | Unit vectors, rotors |
//!
//! ## Projective GA (PGA)
//!
//! | Wrapper | Type | Constraint | Use Case |
//! |---------|------|------------|----------|
//! | [`Bulk<T>`] | Normalization | `bulk_norm() == 1` | Versors (motors, flectors) |
//! | [`Unitized<T>`] | Normalization | `weight_norm() == 1` | Finite points, planes |
//! | [`Ideal<T>`] | Constraint | `weight_norm() ≈ 0` | Points/planes at infinity |
//!
//! ## Minkowski (Indefinite) Algebras
//!
//! | Wrapper | Type | Constraint | Use Case |
//! |---------|------|------------|----------|
//! | [`Proper<T>`] | Normalization | timelike, `|norm²| == 1` | 4-velocities |
//! | [`Spacelike<T>`] | Normalization | spacelike, `|norm²| == 1` | Spatial directions |
//! | [`Null<T>`] | Constraint | `norm² ≈ 0` | Photon worldlines |
//!
//! # Normalization vs Constraint Wrappers
//!
//! **Normalization wrappers** scale the element to satisfy a norm:
//! - `Unit::try_new(v)` → divides by `norm()`, returns `None` if zero
//! - `Unitized::try_new(p)` → divides by `weight_norm()`, returns `None` if ideal
//!
//! **Constraint wrappers** verify a property without scaling:
//! - `Ideal::try_new(p)` → returns `Some` only if `weight_norm() ≈ 0`
//! - `Null::try_new(v)` → returns `Some` only if `norm² ≈ 0`
//!
//! # Example
//!
//! ```ignore
//! use clifford::wrappers::{Unit, Unitized, Ideal};
//! use clifford::specialized::euclidean::dim3::Vector;
//! use clifford::specialized::projective::dim3::Point;
//!
//! // Euclidean: normalize a vector
//! let v = Vector::new(3.0, 4.0, 0.0);
//! let unit = Unit::new_normalize(v);
//! assert!((unit.norm() - 1.0).abs() < 1e-10);
//!
//! // PGA: unitize a finite point
//! let p = Point::new(2.0, 4.0, 6.0, 2.0);
//! let unitized = Unitized::new_normalize(p);  // w = 1 form
//!
//! // PGA: constrain to ideal point (direction)
//! let dir = Point::new(1.0, 0.0, 0.0, 0.0);  // w = 0
//! let ideal = Ideal::try_new(dir).unwrap();
//! ```
//!
//! # Naming Rationale
//!
//! - **`Unit<T>`** - Standard mathematical term for norm = 1
//! - **`Bulk<T>`** - PGA terminology for normalized by non-degenerate part
//! - **`Unitized<T>`** - PGA terminology for standard homogeneous form (weight = 1)
//! - **`Ideal<T>`** - PGA terminology for elements at infinity (weight = 0)
//! - **`Proper<T>`** - Physics terminology for proper time/proper velocity (timelike)
//! - **`Spacelike<T>`** - Physics terminology for spatial separations
//! - **`Null<T>`** - Physics terminology for lightlike/null vectors
//!
//! # References
//!
//! - [Rigid GA Wiki - Geometric norm](https://rigidgeometricalgebra.org/wiki/index.php?title=Geometric_norm)
//! - [Rigid GA Wiki - Unitization](https://rigidgeometricalgebra.org/wiki/index.php?title=Unitization)

use core::fmt::Debug;
use core::hash::Hash;
use core::ops::Deref;

use crate::norm::{DegenerateNormed, IndefiniteNormed, Normed};
use crate::scalar::Float;
// Import num_traits for method access
use num_traits::Float as _;
use num_traits::One;

// ============================================================================
// Unit<T> - Euclidean Norm Wrapper
// ============================================================================

/// A wrapper that guarantees the inner value has unit Euclidean norm.
///
/// `Unit<T>` provides compile-time documentation that a value is normalized,
/// and enables type-safe APIs that require normalized inputs.
///
/// # Example
///
/// ```ignore
/// use clifford::wrappers::Unit;
/// use clifford::specialized::euclidean::dim3::Vector;
///
/// // Create a vector and normalize it
/// let v = Vector::new(3.0, 4.0, 0.0);
/// let unit = Unit::new_normalize(v);
///
/// // Unit vectors have norm = 1
/// assert!((unit.norm() - 1.0).abs() < 1e-10);
///
/// // Access methods via Deref
/// let x = unit.x();  // Works because Unit<T> derefs to T
/// ```
///
/// # When to Use
///
/// Use `Unit<T>` for:
/// - Euclidean vectors representing directions (not magnitudes)
/// - Rotation axes
/// - Normalized rotors/spinors
/// - Any quantity where magnitude should always be 1
#[derive(Clone, Copy)]
#[repr(transparent)]
pub struct Unit<T> {
    /// The wrapped value, guaranteed to satisfy the wrapper's constraint.
    inner: T,
}

impl<T> Unit<T> {
    /// Creates a `Unit<T>` without checking or enforcing normalization.
    ///
    /// # Safety
    ///
    /// The caller must ensure the inner value has unit norm. Using this
    /// with non-normalized values violates the type's invariant.
    ///
    /// Use this only when:
    /// - The value is known to be normalized (e.g., from a factory method)
    /// - Performance is critical and normalization has been verified elsewhere
    #[inline]
    pub fn new_unchecked(inner: T) -> Self {
        Self { inner }
    }

    /// Returns the inner value, consuming the wrapper.
    #[inline]
    pub fn into_inner(self) -> T {
        self.inner
    }

    /// Returns a reference to the inner value.
    #[inline]
    pub fn as_inner(&self) -> &T {
        &self.inner
    }
}

impl<T: Normed> Unit<T> {
    /// Creates a `Unit<T>` by normalizing the input.
    ///
    /// Returns `None` if the input has zero or near-zero norm.
    ///
    /// # Example
    ///
    /// ```ignore
    /// let v = Vector::new(3.0, 4.0, 0.0);
    /// let unit = Unit::try_new(v).unwrap();
    ///
    /// let zero = Vector::new(0.0, 0.0, 0.0);
    /// assert!(Unit::try_new(zero).is_none());
    /// ```
    #[inline]
    pub fn try_new(inner: T) -> Option<Self>
    where
        T: Sized,
    {
        inner
            .try_normalize()
            .map(|normalized| Self { inner: normalized })
    }

    /// Creates a `Unit<T>` by normalizing the input.
    ///
    /// # Panics
    ///
    /// Panics if the input has zero or near-zero norm.
    ///
    /// # Example
    ///
    /// ```ignore
    /// let v = Vector::new(3.0, 4.0, 0.0);
    /// let unit = Unit::new_normalize(v);
    /// assert!((unit.norm() - 1.0).abs() < 1e-10);
    /// ```
    #[inline]
    pub fn new_normalize(inner: T) -> Self
    where
        T: Sized,
    {
        Self::try_new(inner).expect("cannot normalize zero element")
    }
}

impl<T> Deref for Unit<T> {
    type Target = T;

    #[inline]
    fn deref(&self) -> &Self::Target {
        &self.inner
    }
}

impl<T> AsRef<T> for Unit<T> {
    #[inline]
    fn as_ref(&self) -> &T {
        &self.inner
    }
}

impl<T: Debug> Debug for Unit<T> {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        f.debug_tuple("Unit").field(&self.inner).finish()
    }
}

impl<T: PartialEq> PartialEq for Unit<T> {
    fn eq(&self, other: &Self) -> bool {
        self.inner == other.inner
    }
}

impl<T: Eq> Eq for Unit<T> {}

impl<T: Hash> Hash for Unit<T> {
    fn hash<H: core::hash::Hasher>(&self, state: &mut H) {
        self.inner.hash(state);
    }
}

#[cfg(feature = "serde")]
impl<T: serde::Serialize> serde::Serialize for Unit<T> {
    fn serialize<S>(&self, serializer: S) -> Result<S::Ok, S::Error>
    where
        S: serde::Serializer,
    {
        self.inner.serialize(serializer)
    }
}

#[cfg(feature = "serde")]
impl<'de, T: serde::Deserialize<'de>> serde::Deserialize<'de> for Unit<T> {
    fn deserialize<D>(deserializer: D) -> Result<Self, D::Error>
    where
        D: serde::Deserializer<'de>,
    {
        let inner = T::deserialize(deserializer)?;
        Ok(Self { inner })
    }
}

// ============================================================================
// Bulk<T> - PGA Versor Wrapper
// ============================================================================

/// A wrapper guaranteeing bulk norm = 1 (PGA rigid transformations).
///
/// In Projective Geometric Algebra (PGA), the bulk norm measures the
/// non-degenerate part of an element. For motors (versors), `bulk_norm() == 1`
/// guarantees a proper rigid body transformation (rotation + translation,
/// no scaling).
///
/// # Example
///
/// ```ignore
/// use clifford::wrappers::Bulk;
/// use clifford::specialized::projective::dim3::Motor;
///
/// // Create a motor and bulk-normalize it
/// let motor = Motor::from_translation(1.0, 2.0, 3.0);
/// let rigid = Bulk::new_normalize(motor);
///
/// // Bulk-normalized motors represent rigid transformations
/// assert!((rigid.bulk_norm() - 1.0).abs() < 1e-10);
/// ```
///
/// # When to Use
///
/// Use `Bulk<T>` for:
/// - PGA motors representing rigid body transformations
/// - PGA flectors representing rigid reflections
/// - Any PGA versor where scaling should be excluded
///
/// # Reference
///
/// [Rigid GA Wiki - Geometric norm](https://rigidgeometricalgebra.org/wiki/index.php?title=Geometric_norm)
#[derive(Clone, Copy)]
#[repr(transparent)]
pub struct Bulk<T> {
    /// The wrapped value, guaranteed to satisfy the wrapper's constraint.
    inner: T,
}

impl<T: DegenerateNormed> Bulk<T>
where
    T::Scalar: Float,
{
    /// Creates a `Bulk<T>` by normalizing by the bulk norm.
    ///
    /// Returns `None` if the bulk norm is zero or near-zero.
    ///
    /// # Example
    ///
    /// ```ignore
    /// let motor = Motor::from_rotation_z(0.5);
    /// let rigid = Bulk::try_new(motor).unwrap();
    /// ```
    #[inline]
    pub fn try_new(inner: T) -> Option<Self>
    where
        T: Sized,
    {
        let n = inner.bulk_norm();
        if n < T::Scalar::epsilon() {
            return None;
        }
        Some(Self {
            inner: inner.scale(<T::Scalar as One>::one() / n),
        })
    }

    /// Creates a `Bulk<T>` by normalizing by the bulk norm.
    ///
    /// # Panics
    ///
    /// Panics if the bulk norm is zero or near-zero.
    #[inline]
    pub fn new_normalize(inner: T) -> Self
    where
        T: Sized,
    {
        Self::try_new(inner).expect("cannot bulk-normalize zero element")
    }
}

impl<T> Bulk<T> {
    /// Creates a `Bulk<T>` without checking or enforcing normalization.
    ///
    /// # Safety
    ///
    /// The caller must ensure the inner value has unit bulk norm.
    #[inline]
    pub fn new_unchecked(inner: T) -> Self {
        Self { inner }
    }

    /// Returns the inner value, consuming the wrapper.
    #[inline]
    pub fn into_inner(self) -> T {
        self.inner
    }

    /// Returns a reference to the inner value.
    #[inline]
    pub fn as_inner(&self) -> &T {
        &self.inner
    }
}

impl<T> Deref for Bulk<T> {
    type Target = T;

    #[inline]
    fn deref(&self) -> &Self::Target {
        &self.inner
    }
}

impl<T> AsRef<T> for Bulk<T> {
    #[inline]
    fn as_ref(&self) -> &T {
        &self.inner
    }
}

impl<T: Debug> Debug for Bulk<T> {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        f.debug_tuple("Bulk").field(&self.inner).finish()
    }
}

impl<T: PartialEq> PartialEq for Bulk<T> {
    fn eq(&self, other: &Self) -> bool {
        self.inner == other.inner
    }
}

impl<T: Eq> Eq for Bulk<T> {}

impl<T: Hash> Hash for Bulk<T> {
    fn hash<H: core::hash::Hasher>(&self, state: &mut H) {
        self.inner.hash(state);
    }
}

#[cfg(feature = "serde")]
impl<T: serde::Serialize> serde::Serialize for Bulk<T> {
    fn serialize<S>(&self, serializer: S) -> Result<S::Ok, S::Error>
    where
        S: serde::Serializer,
    {
        self.inner.serialize(serializer)
    }
}

#[cfg(feature = "serde")]
impl<'de, T: serde::Deserialize<'de>> serde::Deserialize<'de> for Bulk<T> {
    fn deserialize<D>(deserializer: D) -> Result<Self, D::Error>
    where
        D: serde::Deserializer<'de>,
    {
        let inner = T::deserialize(deserializer)?;
        Ok(Self { inner })
    }
}

// ============================================================================
// Unitized<T> - PGA Homogeneous Wrapper (weight = 1)
// ============================================================================

/// A wrapper guaranteeing weight norm = 1 (standard homogeneous form).
///
/// In Projective Geometric Algebra (PGA), the weight norm measures the
/// projective/ideal (degenerate) part. For points, `weight_norm() == 1`
/// means the point is in standard homogeneous form with w = 1.
///
/// # Unitized vs Ideal
///
/// - **`Unitized<T>`**: Normalization wrapper - scales so weight_norm = 1 (finite elements)
/// - **`Ideal<T>`**: Constraint wrapper - verifies weight_norm ≈ 0 (elements at infinity)
///
/// # Example
///
/// ```ignore
/// use clifford::wrappers::Unitized;
/// use clifford::specialized::projective::dim3::Point;
///
/// // Create a homogeneous point and unitize to standard form
/// let p = Point::new(2.0, 4.0, 6.0, 2.0);  // Represents (1, 2, 3)
/// let std = Unitized::new_normalize(p);
///
/// // Standard form has weight norm = 1
/// assert!((std.weight_norm() - 1.0).abs() < 1e-10);
/// ```
///
/// # When to Use
///
/// Use `Unitized<T>` for:
/// - PGA points in standard homogeneous form (finite points)
/// - PGA planes in standard form
/// - Any element where the projective coordinate should be normalized to 1
///
/// # Unitization vs Normalization
///
/// In PGA terminology:
/// - **Normalize** = divide by bulk norm (use `Bulk<T>`)
/// - **Unitize** = divide by weight norm (use `Unitized<T>`)
///
/// # Reference
///
/// [Rigid GA Wiki - Unitization](https://rigidgeometricalgebra.org/wiki/index.php?title=Unitization)
#[derive(Clone, Copy)]
#[repr(transparent)]
pub struct Unitized<T> {
    /// The wrapped value, guaranteed to satisfy the wrapper's constraint.
    inner: T,
}

impl<T: DegenerateNormed> Unitized<T> {
    /// Creates a `Unitized<T>` by dividing by weight norm.
    ///
    /// Returns `None` if the weight norm is zero or near-zero
    /// (e.g., ideal/infinity points cannot be unitized).
    ///
    /// # Example
    ///
    /// ```ignore
    /// let p = Point::new(2.0, 4.0, 6.0, 2.0);
    /// let std = Unitized::try_new(p).unwrap();
    ///
    /// // Ideal points (at infinity) cannot be unitized
    /// let ideal = Point::new(1.0, 0.0, 0.0, 0.0);  // w = 0
    /// assert!(Unitized::try_new(ideal).is_none());
    /// ```
    #[inline]
    pub fn try_new(inner: T) -> Option<Self>
    where
        T: Sized,
    {
        inner.try_unitize().map(|u| Self { inner: u })
    }

    /// Creates a `Unitized<T>` by dividing by weight norm.
    ///
    /// # Panics
    ///
    /// Panics if the weight norm is zero or near-zero.
    #[inline]
    pub fn new_normalize(inner: T) -> Self
    where
        T: Sized,
    {
        Self::try_new(inner).expect("cannot unitize element with zero weight (ideal element)")
    }
}

impl<T> Unitized<T> {
    /// Creates a `Unitized<T>` without checking or enforcing normalization.
    ///
    /// # Safety
    ///
    /// The caller must ensure the inner value has unit weight norm.
    #[inline]
    pub fn new_unchecked(inner: T) -> Self {
        Self { inner }
    }

    /// Returns the inner value, consuming the wrapper.
    #[inline]
    pub fn into_inner(self) -> T {
        self.inner
    }

    /// Returns a reference to the inner value.
    #[inline]
    pub fn as_inner(&self) -> &T {
        &self.inner
    }
}

impl<T> Deref for Unitized<T> {
    type Target = T;

    #[inline]
    fn deref(&self) -> &Self::Target {
        &self.inner
    }
}

impl<T> AsRef<T> for Unitized<T> {
    #[inline]
    fn as_ref(&self) -> &T {
        &self.inner
    }
}

impl<T: Debug> Debug for Unitized<T> {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        f.debug_tuple("Unitized").field(&self.inner).finish()
    }
}

impl<T: PartialEq> PartialEq for Unitized<T> {
    fn eq(&self, other: &Self) -> bool {
        self.inner == other.inner
    }
}

impl<T: Eq> Eq for Unitized<T> {}

impl<T: Hash> Hash for Unitized<T> {
    fn hash<H: core::hash::Hasher>(&self, state: &mut H) {
        self.inner.hash(state);
    }
}

#[cfg(feature = "serde")]
impl<T: serde::Serialize> serde::Serialize for Unitized<T> {
    fn serialize<S>(&self, serializer: S) -> Result<S::Ok, S::Error>
    where
        S: serde::Serializer,
    {
        self.inner.serialize(serializer)
    }
}

#[cfg(feature = "serde")]
impl<'de, T: serde::Deserialize<'de>> serde::Deserialize<'de> for Unitized<T> {
    fn deserialize<D>(deserializer: D) -> Result<Self, D::Error>
    where
        D: serde::Deserializer<'de>,
    {
        let inner = T::deserialize(deserializer)?;
        Ok(Self { inner })
    }
}

// ============================================================================
// Ideal<T> - PGA Constraint Wrapper (weight = 0)
// ============================================================================

/// A constraint wrapper guaranteeing weight norm ≈ 0 (elements at infinity).
///
/// In Projective Geometric Algebra (PGA), elements with zero weight are called
/// "ideal" - they represent geometric entities at infinity. For example, an
/// ideal point represents a direction (point at infinity).
///
/// **Important**: This is a *constraint* wrapper, not a *normalization* wrapper.
/// It verifies that the element has near-zero weight but does NOT scale it.
/// You cannot normalize an ideal element by weight because weight = 0.
///
/// # Ideal vs Unitized
///
/// - **`Ideal<T>`**: Constraint wrapper - verifies weight_norm ≈ 0 (at infinity)
/// - **`Unitized<T>`**: Normalization wrapper - scales so weight_norm = 1 (finite)
///
/// # Example
///
/// ```ignore
/// use clifford::wrappers::Ideal;
/// use clifford::specialized::projective::dim3::Point;
///
/// // An ideal point (direction) has w = 0
/// let direction = Point::new(1.0, 0.0, 0.0, 0.0);  // Points in +x direction
/// let ideal = Ideal::try_new(direction).unwrap();
///
/// // Finite points cannot be wrapped as Ideal
/// let finite = Point::new(1.0, 2.0, 3.0, 1.0);  // w = 1
/// assert!(Ideal::try_new(finite).is_none());
/// ```
///
/// # When to Use
///
/// Use `Ideal<T>` for:
/// - PGA ideal points (directions, points at infinity)
/// - PGA ideal planes (plane at infinity)
/// - Type-safe APIs that require ideal/infinite elements
///
/// # Reference
///
/// [Rigid GA Wiki - Ideal elements](https://rigidgeometricalgebra.org/wiki/index.php?title=Projective_geometric_algebra)
#[derive(Clone, Copy)]
#[repr(transparent)]
pub struct Ideal<T> {
    /// The wrapped value, guaranteed to have weight ≈ 0.
    inner: T,
}

impl<T: DegenerateNormed> Ideal<T> {
    /// Creates an `Ideal<T>` if the element has near-zero weight.
    ///
    /// This is a constraint check, NOT a normalization. The element is
    /// wrapped as-is if its weight norm is below epsilon.
    ///
    /// Returns `None` if the weight norm is not near zero (finite element).
    ///
    /// # Example
    ///
    /// ```ignore
    /// // Ideal point (direction) - weight = 0
    /// let direction = Point::new(1.0, 0.0, 0.0, 0.0);
    /// let ideal = Ideal::try_new(direction).unwrap();
    ///
    /// // Finite point - weight ≠ 0, rejected
    /// let finite = Point::new(1.0, 2.0, 3.0, 1.0);
    /// assert!(Ideal::try_new(finite).is_none());
    /// ```
    #[inline]
    pub fn try_new(inner: T) -> Option<Self>
    where
        T: Sized,
    {
        if inner.weight_norm() < T::Scalar::epsilon() {
            Some(Self { inner })
        } else {
            None // Not ideal - has finite weight
        }
    }
}

impl<T> Ideal<T> {
    /// Creates an `Ideal<T>` without checking the weight constraint.
    ///
    /// # Safety
    ///
    /// The caller must ensure the inner value has near-zero weight norm.
    #[inline]
    pub fn new_unchecked(inner: T) -> Self {
        Self { inner }
    }

    /// Returns the inner value, consuming the wrapper.
    #[inline]
    pub fn into_inner(self) -> T {
        self.inner
    }

    /// Returns a reference to the inner value.
    #[inline]
    pub fn as_inner(&self) -> &T {
        &self.inner
    }
}

impl<T> Deref for Ideal<T> {
    type Target = T;

    #[inline]
    fn deref(&self) -> &Self::Target {
        &self.inner
    }
}

impl<T> AsRef<T> for Ideal<T> {
    #[inline]
    fn as_ref(&self) -> &T {
        &self.inner
    }
}

impl<T: Debug> Debug for Ideal<T> {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        f.debug_tuple("Ideal").field(&self.inner).finish()
    }
}

impl<T: PartialEq> PartialEq for Ideal<T> {
    fn eq(&self, other: &Self) -> bool {
        self.inner == other.inner
    }
}

impl<T: Eq> Eq for Ideal<T> {}

impl<T: Hash> Hash for Ideal<T> {
    fn hash<H: core::hash::Hasher>(&self, state: &mut H) {
        self.inner.hash(state);
    }
}

#[cfg(feature = "serde")]
impl<T: serde::Serialize> serde::Serialize for Ideal<T> {
    fn serialize<S>(&self, serializer: S) -> Result<S::Ok, S::Error>
    where
        S: serde::Serializer,
    {
        self.inner.serialize(serializer)
    }
}

#[cfg(feature = "serde")]
impl<'de, T: serde::Deserialize<'de>> serde::Deserialize<'de> for Ideal<T> {
    fn deserialize<D>(deserializer: D) -> Result<Self, D::Error>
    where
        D: serde::Deserializer<'de>,
    {
        let inner = T::deserialize(deserializer)?;
        Ok(Self { inner })
    }
}

// ============================================================================
// Proper<T> - Minkowski Timelike Wrapper
// ============================================================================

/// A wrapper guaranteeing the vector is timelike (proper).
///
/// In Minkowski spacetime with signature `(+,-,-,-)`, timelike vectors have
/// positive squared norm. They represent possible 4-velocities of massive
/// particles, proper time intervals, and other "proper" quantities.
///
/// # Example
///
/// ```ignore
/// use clifford::wrappers::Proper;
/// use clifford::specialized::minkowski::dim4::FourVector;
///
/// // Timelike vector (t² - x² - y² - z² > 0)
/// let v = FourVector::new(2.0, 1.0, 0.0, 0.0);  // 4 - 1 = 3 > 0
/// let proper = Proper::try_new(v).unwrap();
///
/// // Spacelike vectors are rejected
/// let spacelike = FourVector::new(1.0, 2.0, 0.0, 0.0);  // 1 - 4 < 0
/// assert!(Proper::try_new(spacelike).is_none());
/// ```
///
/// # When to Use
///
/// Use `Proper<T>` for:
/// - 4-velocities (must be timelike)
/// - Proper time intervals
/// - Any quantity that must be inside the light cone
///
/// # Reference
///
/// [Spacetime Algebra](https://en.wikipedia.org/wiki/Spacetime_algebra)
#[derive(Clone, Copy)]
#[repr(transparent)]
pub struct Proper<T> {
    /// The wrapped value, guaranteed to satisfy the wrapper's constraint.
    inner: T,
}

impl<T: IndefiniteNormed> Proper<T> {
    /// Creates a `Proper<T>` if the input is timelike.
    ///
    /// Returns `None` if the input is spacelike or lightlike.
    /// The input is also normalized.
    ///
    /// # Example
    ///
    /// ```ignore
    /// let timelike = FourVector::new(2.0, 1.0, 0.0, 0.0);
    /// let proper = Proper::try_new(timelike).unwrap();
    ///
    /// let spacelike = FourVector::new(1.0, 2.0, 0.0, 0.0);
    /// assert!(Proper::try_new(spacelike).is_none());
    /// ```
    #[inline]
    pub fn try_new(inner: T) -> Option<Self>
    where
        T: Sized,
    {
        if inner.is_timelike() {
            Some(Self {
                inner: inner.normalize(),
            })
        } else {
            None // Spacelike or lightlike
        }
    }
}

impl<T> Proper<T> {
    /// Creates a `Proper<T>` without checking the timelike property.
    ///
    /// # Safety
    ///
    /// The caller must ensure the inner value is timelike.
    #[inline]
    pub fn new_unchecked(inner: T) -> Self {
        Self { inner }
    }

    /// Returns the inner value, consuming the wrapper.
    #[inline]
    pub fn into_inner(self) -> T {
        self.inner
    }

    /// Returns a reference to the inner value.
    #[inline]
    pub fn as_inner(&self) -> &T {
        &self.inner
    }
}

impl<T> Deref for Proper<T> {
    type Target = T;

    #[inline]
    fn deref(&self) -> &Self::Target {
        &self.inner
    }
}

impl<T> AsRef<T> for Proper<T> {
    #[inline]
    fn as_ref(&self) -> &T {
        &self.inner
    }
}

impl<T: Debug> Debug for Proper<T> {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        f.debug_tuple("Proper").field(&self.inner).finish()
    }
}

impl<T: PartialEq> PartialEq for Proper<T> {
    fn eq(&self, other: &Self) -> bool {
        self.inner == other.inner
    }
}

impl<T: Eq> Eq for Proper<T> {}

impl<T: Hash> Hash for Proper<T> {
    fn hash<H: core::hash::Hasher>(&self, state: &mut H) {
        self.inner.hash(state);
    }
}

#[cfg(feature = "serde")]
impl<T: serde::Serialize> serde::Serialize for Proper<T> {
    fn serialize<S>(&self, serializer: S) -> Result<S::Ok, S::Error>
    where
        S: serde::Serializer,
    {
        self.inner.serialize(serializer)
    }
}

#[cfg(feature = "serde")]
impl<'de, T: serde::Deserialize<'de>> serde::Deserialize<'de> for Proper<T> {
    fn deserialize<D>(deserializer: D) -> Result<Self, D::Error>
    where
        D: serde::Deserializer<'de>,
    {
        let inner = T::deserialize(deserializer)?;
        Ok(Self { inner })
    }
}

// ============================================================================
// Spacelike<T> - Minkowski Spacelike Wrapper
// ============================================================================

/// A wrapper guaranteeing the vector is spacelike and normalized.
///
/// In Minkowski spacetime with signature `(+,-,-,-)`, spacelike vectors have
/// negative squared norm. They represent spatial separations between events
/// that are outside each other's light cones.
///
/// # Example
///
/// ```ignore
/// use clifford::wrappers::Spacelike;
/// use clifford::specialized::minkowski::dim4::FourVector;
///
/// // Spacelike vector (t² - x² - y² - z² < 0)
/// let v = FourVector::new(1.0, 2.0, 0.0, 0.0);  // 1 - 4 = -3 < 0
/// let spacelike = Spacelike::try_new(v).unwrap();
///
/// // Timelike vectors are rejected
/// let timelike = FourVector::new(2.0, 1.0, 0.0, 0.0);  // 4 - 1 > 0
/// assert!(Spacelike::try_new(timelike).is_none());
/// ```
///
/// # When to Use
///
/// Use `Spacelike<T>` for:
/// - Spatial separation vectors
/// - Space-like hypersurface normals
/// - Any quantity that must be outside the light cone
///
/// # Reference
///
/// [Spacetime Algebra](https://en.wikipedia.org/wiki/Spacetime_algebra)
#[derive(Clone, Copy)]
#[repr(transparent)]
pub struct Spacelike<T> {
    /// The wrapped value, guaranteed to be spacelike and normalized.
    inner: T,
}

impl<T: IndefiniteNormed> Spacelike<T> {
    /// Creates a `Spacelike<T>` if the input is spacelike.
    ///
    /// Returns `None` if the input is timelike or lightlike.
    /// The input is normalized by its magnitude (√|v·v|).
    ///
    /// # Example
    ///
    /// ```ignore
    /// let spacelike = FourVector::new(1.0, 2.0, 0.0, 0.0);  // 1 - 4 < 0
    /// let normalized = Spacelike::try_new(spacelike).unwrap();
    ///
    /// let timelike = FourVector::new(2.0, 1.0, 0.0, 0.0);
    /// assert!(Spacelike::try_new(timelike).is_none());
    /// ```
    #[inline]
    pub fn try_new(inner: T) -> Option<Self>
    where
        T: Sized,
    {
        if inner.is_spacelike() {
            // Normalize by magnitude: √|v·v|
            let mag = inner.norm(); // norm() returns √|norm_squared()|
            if mag > T::Scalar::epsilon() {
                Some(Self {
                    inner: inner.scale(<T::Scalar as One>::one() / mag),
                })
            } else {
                None
            }
        } else {
            None // Timelike or lightlike
        }
    }
}

impl<T> Spacelike<T> {
    /// Creates a `Spacelike<T>` without checking the spacelike property.
    ///
    /// # Safety
    ///
    /// The caller must ensure the inner value is spacelike and normalized.
    #[inline]
    pub fn new_unchecked(inner: T) -> Self {
        Self { inner }
    }

    /// Returns the inner value, consuming the wrapper.
    #[inline]
    pub fn into_inner(self) -> T {
        self.inner
    }

    /// Returns a reference to the inner value.
    #[inline]
    pub fn as_inner(&self) -> &T {
        &self.inner
    }
}

impl<T> Deref for Spacelike<T> {
    type Target = T;

    #[inline]
    fn deref(&self) -> &Self::Target {
        &self.inner
    }
}

impl<T> AsRef<T> for Spacelike<T> {
    #[inline]
    fn as_ref(&self) -> &T {
        &self.inner
    }
}

impl<T: Debug> Debug for Spacelike<T> {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        f.debug_tuple("Spacelike").field(&self.inner).finish()
    }
}

impl<T: PartialEq> PartialEq for Spacelike<T> {
    fn eq(&self, other: &Self) -> bool {
        self.inner == other.inner
    }
}

impl<T: Eq> Eq for Spacelike<T> {}

impl<T: Hash> Hash for Spacelike<T> {
    fn hash<H: core::hash::Hasher>(&self, state: &mut H) {
        self.inner.hash(state);
    }
}

#[cfg(feature = "serde")]
impl<T: serde::Serialize> serde::Serialize for Spacelike<T> {
    fn serialize<S>(&self, serializer: S) -> Result<S::Ok, S::Error>
    where
        S: serde::Serializer,
    {
        self.inner.serialize(serializer)
    }
}

#[cfg(feature = "serde")]
impl<'de, T: serde::Deserialize<'de>> serde::Deserialize<'de> for Spacelike<T> {
    fn deserialize<D>(deserializer: D) -> Result<Self, D::Error>
    where
        D: serde::Deserializer<'de>,
    {
        let inner = T::deserialize(deserializer)?;
        Ok(Self { inner })
    }
}

// ============================================================================
// Null<T> - Minkowski Lightlike Constraint Wrapper
// ============================================================================

/// A constraint wrapper guaranteeing the vector is lightlike/null.
///
/// In Minkowski spacetime, lightlike (or null) vectors have zero squared norm.
/// They represent the world-lines of massless particles (photons) and lie
/// exactly on the light cone.
///
/// **Important**: This is a *constraint* wrapper, not a *normalization* wrapper.
/// Null vectors cannot be normalized in the traditional sense because their
/// norm is zero. The wrapper verifies the constraint but does not scale.
///
/// # Example
///
/// ```ignore
/// use clifford::wrappers::Null;
/// use clifford::specialized::minkowski::dim4::FourVector;
///
/// // Lightlike vector (t² - x² - y² - z² = 0)
/// let photon = FourVector::new(1.0, 1.0, 0.0, 0.0);  // 1 - 1 = 0
/// let null = Null::try_new(photon).unwrap();
///
/// // Non-null vectors are rejected
/// let timelike = FourVector::new(2.0, 1.0, 0.0, 0.0);
/// assert!(Null::try_new(timelike).is_none());
/// ```
///
/// # When to Use
///
/// Use `Null<T>` for:
/// - Photon 4-momenta
/// - Light cone generators
/// - Any quantity that must lie on the light cone
///
/// # Reference
///
/// [Null vector](https://en.wikipedia.org/wiki/Null_vector)
#[derive(Clone, Copy)]
#[repr(transparent)]
pub struct Null<T> {
    /// The wrapped value, guaranteed to have norm ≈ 0.
    inner: T,
}

impl<T: IndefiniteNormed> Null<T> {
    /// Creates a `Null<T>` if the element is lightlike (norm² ≈ 0).
    ///
    /// This is a constraint check, NOT a normalization. The element is
    /// wrapped as-is if its squared norm is near zero.
    ///
    /// Returns `None` if the element is not lightlike.
    ///
    /// # Example
    ///
    /// ```ignore
    /// let photon = FourVector::new(1.0, 1.0, 0.0, 0.0);
    /// let null = Null::try_new(photon).unwrap();
    ///
    /// let timelike = FourVector::new(2.0, 1.0, 0.0, 0.0);
    /// assert!(Null::try_new(timelike).is_none());
    /// ```
    #[inline]
    pub fn try_new(inner: T) -> Option<Self>
    where
        T: Sized,
    {
        if inner.is_lightlike() {
            Some(Self { inner })
        } else {
            None // Not null
        }
    }
}

impl<T> Null<T> {
    /// Creates a `Null<T>` without checking the lightlike constraint.
    ///
    /// # Safety
    ///
    /// The caller must ensure the inner value is lightlike (norm² ≈ 0).
    #[inline]
    pub fn new_unchecked(inner: T) -> Self {
        Self { inner }
    }

    /// Returns the inner value, consuming the wrapper.
    #[inline]
    pub fn into_inner(self) -> T {
        self.inner
    }

    /// Returns a reference to the inner value.
    #[inline]
    pub fn as_inner(&self) -> &T {
        &self.inner
    }
}

impl<T> Deref for Null<T> {
    type Target = T;

    #[inline]
    fn deref(&self) -> &Self::Target {
        &self.inner
    }
}

impl<T> AsRef<T> for Null<T> {
    #[inline]
    fn as_ref(&self) -> &T {
        &self.inner
    }
}

impl<T: Debug> Debug for Null<T> {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        f.debug_tuple("Null").field(&self.inner).finish()
    }
}

impl<T: PartialEq> PartialEq for Null<T> {
    fn eq(&self, other: &Self) -> bool {
        self.inner == other.inner
    }
}

impl<T: Eq> Eq for Null<T> {}

impl<T: Hash> Hash for Null<T> {
    fn hash<H: core::hash::Hasher>(&self, state: &mut H) {
        self.inner.hash(state);
    }
}

#[cfg(feature = "serde")]
impl<T: serde::Serialize> serde::Serialize for Null<T> {
    fn serialize<S>(&self, serializer: S) -> Result<S::Ok, S::Error>
    where
        S: serde::Serializer,
    {
        self.inner.serialize(serializer)
    }
}

#[cfg(feature = "serde")]
impl<'de, T: serde::Deserialize<'de>> serde::Deserialize<'de> for Null<T> {
    fn deserialize<D>(deserializer: D) -> Result<Self, D::Error>
    where
        D: serde::Deserializer<'de>,
    {
        let inner = T::deserialize(deserializer)?;
        Ok(Self { inner })
    }
}

// ============================================================================
// Arbitrary Implementations
// ============================================================================

/// Proptest `Arbitrary` implementations for wrapper types.
///
/// These implementations generate valid wrapped values by:
/// 1. Generating a random inner value via `any::<T>()`
/// 2. Filtering/mapping to apply the wrapper's normalization
///
/// This requires the inner type `T` to implement both the appropriate
/// norm trait (`Normed`, `DegenerateNormed`, `IndefiniteNormed`) and `Arbitrary`.
#[cfg(any(test, feature = "proptest-support"))]
mod arbitrary_impl {
    use super::*;
    use proptest::prelude::*;
    use proptest::strategy::BoxedStrategy;

    impl<T> Arbitrary for Unit<T>
    where
        T: Normed + Arbitrary + Debug + 'static,
    {
        type Parameters = ();
        type Strategy = BoxedStrategy<Self>;

        fn arbitrary_with(_: Self::Parameters) -> Self::Strategy {
            any::<T>()
                .prop_filter_map("normalizable", |t| Unit::try_new(t))
                .boxed()
        }
    }

    impl<T> Arbitrary for Bulk<T>
    where
        T: DegenerateNormed + Arbitrary + Debug + 'static,
        T::Scalar: Float,
    {
        type Parameters = ();
        type Strategy = BoxedStrategy<Self>;

        fn arbitrary_with(_: Self::Parameters) -> Self::Strategy {
            any::<T>()
                .prop_filter_map("bulk-normalizable", |t| Bulk::try_new(t))
                .boxed()
        }
    }

    impl<T> Arbitrary for Unitized<T>
    where
        T: DegenerateNormed + Arbitrary + Debug + 'static,
    {
        type Parameters = ();
        type Strategy = BoxedStrategy<Self>;

        fn arbitrary_with(_: Self::Parameters) -> Self::Strategy {
            any::<T>()
                .prop_filter_map("unitizable (weight > 0)", |t| Unitized::try_new(t))
                .boxed()
        }
    }

    impl<T> Arbitrary for Ideal<T>
    where
        T: DegenerateNormed + Arbitrary + Debug + 'static,
    {
        type Parameters = ();
        type Strategy = BoxedStrategy<Self>;

        fn arbitrary_with(_: Self::Parameters) -> Self::Strategy {
            // Ideal elements have weight ≈ 0, which is rare in random data.
            // This filter will reject most samples, so Ideal<T> is expensive to generate.
            any::<T>()
                .prop_filter_map("ideal (weight ≈ 0)", |t| Ideal::try_new(t))
                .boxed()
        }
    }

    impl<T> Arbitrary for Proper<T>
    where
        T: IndefiniteNormed + Arbitrary + Debug + 'static,
    {
        type Parameters = ();
        type Strategy = BoxedStrategy<Self>;

        fn arbitrary_with(_: Self::Parameters) -> Self::Strategy {
            any::<T>()
                .prop_filter_map("timelike", |t| Proper::try_new(t))
                .boxed()
        }
    }

    impl<T> Arbitrary for Spacelike<T>
    where
        T: IndefiniteNormed + Arbitrary + Debug + 'static,
    {
        type Parameters = ();
        type Strategy = BoxedStrategy<Self>;

        fn arbitrary_with(_: Self::Parameters) -> Self::Strategy {
            any::<T>()
                .prop_filter_map("spacelike", |t| Spacelike::try_new(t))
                .boxed()
        }
    }

    impl<T> Arbitrary for Null<T>
    where
        T: IndefiniteNormed + Arbitrary + Debug + 'static,
    {
        type Parameters = ();
        type Strategy = BoxedStrategy<Self>;

        fn arbitrary_with(_: Self::Parameters) -> Self::Strategy {
            // Null elements have norm² ≈ 0, which is rare in random data.
            // This filter will reject most samples, so Null<T> is expensive to generate.
            any::<T>()
                .prop_filter_map("lightlike (norm² ≈ 0)", |t| Null::try_new(t))
                .boxed()
        }
    }
}

// ============================================================================
// Tests
// ============================================================================

#[cfg(test)]
mod tests {
    use super::*;

    // Basic structural tests for wrappers
    // Full integration tests require types implementing the norm traits

    #[test]
    fn unit_debug_format() {
        let unit = Unit { inner: 42i32 };
        assert_eq!(format!("{:?}", unit), "Unit(42)");
    }

    #[test]
    fn bulk_debug_format() {
        let bulk = Bulk { inner: 42i32 };
        assert_eq!(format!("{:?}", bulk), "Bulk(42)");
    }

    #[test]
    fn ideal_debug_format() {
        let ideal = Ideal { inner: 42i32 };
        assert_eq!(format!("{:?}", ideal), "Ideal(42)");
    }

    #[test]
    fn proper_debug_format() {
        let proper = Proper { inner: 42i32 };
        assert_eq!(format!("{:?}", proper), "Proper(42)");
    }

    #[test]
    fn unit_deref() {
        let unit = Unit { inner: 42i32 };
        assert_eq!(*unit, 42);
    }

    #[test]
    fn unit_as_ref() {
        let unit = Unit { inner: 42i32 };
        assert_eq!(unit.as_ref(), &42);
    }

    #[test]
    fn unit_into_inner() {
        let unit = Unit { inner: 42i32 };
        assert_eq!(unit.into_inner(), 42);
    }

    #[test]
    fn unit_equality() {
        let a = Unit { inner: 42i32 };
        let b = Unit { inner: 42i32 };
        let c = Unit { inner: 43i32 };
        assert_eq!(a, b);
        assert_ne!(a, c);
    }

    #[test]
    fn unit_clone() {
        let a = Unit { inner: 42i32 };
        let b = a.clone();
        assert_eq!(a, b);
    }

    #[test]
    fn unit_copy() {
        let a = Unit { inner: 42i32 };
        let b = a;
        assert_eq!(a, b);
    }

    #[test]
    fn unitized_debug_format() {
        let unitized = Unitized { inner: 42i32 };
        assert_eq!(format!("{:?}", unitized), "Unitized(42)");
    }

    #[test]
    fn unitized_deref() {
        let unitized = Unitized { inner: 42i32 };
        assert_eq!(*unitized, 42);
    }

    #[test]
    fn unitized_as_ref() {
        let unitized = Unitized { inner: 42i32 };
        assert_eq!(unitized.as_ref(), &42);
    }

    #[test]
    fn unitized_into_inner() {
        let unitized = Unitized { inner: 42i32 };
        assert_eq!(unitized.into_inner(), 42);
    }

    #[test]
    fn unitized_equality() {
        let a = Unitized { inner: 42i32 };
        let b = Unitized { inner: 42i32 };
        let c = Unitized { inner: 43i32 };
        assert_eq!(a, b);
        assert_ne!(a, c);
    }

    #[test]
    fn spacelike_debug_format() {
        let spacelike = Spacelike { inner: 42i32 };
        assert_eq!(format!("{:?}", spacelike), "Spacelike(42)");
    }

    #[test]
    fn spacelike_deref() {
        let spacelike = Spacelike { inner: 42i32 };
        assert_eq!(*spacelike, 42);
    }

    #[test]
    fn spacelike_into_inner() {
        let spacelike = Spacelike { inner: 42i32 };
        assert_eq!(spacelike.into_inner(), 42);
    }

    #[test]
    fn null_debug_format() {
        let null = Null { inner: 42i32 };
        assert_eq!(format!("{:?}", null), "Null(42)");
    }

    #[test]
    fn null_deref() {
        let null = Null { inner: 42i32 };
        assert_eq!(*null, 42);
    }

    #[test]
    fn null_into_inner() {
        let null = Null { inner: 42i32 };
        assert_eq!(null.into_inner(), 42);
    }

    #[test]
    fn null_equality() {
        let a = Null { inner: 42i32 };
        let b = Null { inner: 42i32 };
        let c = Null { inner: 43i32 };
        assert_eq!(a, b);
        assert_ne!(a, c);
    }
}
