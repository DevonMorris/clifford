//! Geometry-specific wrapper types for normalized entities.
//!
//! This module provides wrapper types that guarantee specific normalization
//! properties at the type level. Each geometry type has its own wrapper with
//! a clear name describing what it guarantees:
//!
//! | Wrapper | Constraint | Use Case |
//! |---------|------------|----------|
//! | [`Unit<T>`] | `norm() == 1` | Euclidean vectors, bivectors, rotors |
//! | [`Bulk<T>`] | `bulk_norm() == 1` | PGA versors (motors, flectors) |
//! | [`Ideal<T>`] | `weight_norm() == 1` | PGA homogeneous coords (points, planes) |
//! | [`Proper<T>`] | `is_timelike()` | Minkowski timelike vectors |
//!
//! # Example
//!
//! ```ignore
//! use clifford::wrappers::Unit;
//! use clifford::specialized::euclidean::dim3::Vector;
//!
//! let v = Vector::new(3.0, 4.0, 0.0);
//! let unit = Unit::new_normalize(v);
//! assert!((unit.norm() - 1.0).abs() < 1e-10);
//! ```
//!
//! # Naming Rationale
//!
//! - **`Unit<T>`** - Standard mathematical term for norm = 1
//! - **`Bulk<T>`** - PGA terminology for the non-degenerate part
//! - **`Ideal<T>`** - PGA terminology for the degenerate/projective part
//! - **`Proper<T>`** - Physics terminology for proper time/proper length

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
// Ideal<T> - PGA Homogeneous Wrapper
// ============================================================================

/// A wrapper guaranteeing weight norm = 1 (standard homogeneous form).
///
/// In Projective Geometric Algebra (PGA), the weight norm measures the
/// projective/ideal (degenerate) part. For points, `weight_norm() == 1`
/// means the point is in standard homogeneous form with w = 1.
///
/// # Example
///
/// ```ignore
/// use clifford::wrappers::Ideal;
/// use clifford::specialized::projective::dim3::Point;
///
/// // Create a homogeneous point and normalize to standard form
/// let p = Point::new(2.0, 4.0, 6.0, 2.0);  // Represents (1, 2, 3)
/// let std = Ideal::new_normalize(p);
///
/// // Standard form has weight norm = 1
/// assert!((std.weight_norm() - 1.0).abs() < 1e-10);
/// ```
///
/// # When to Use
///
/// Use `Ideal<T>` for:
/// - PGA points in standard homogeneous form
/// - PGA planes in standard form
/// - Any element where the projective coordinate should be normalized
///
/// # Unitization vs Normalization
///
/// In PGA terminology:
/// - **Normalize** = divide by bulk norm
/// - **Unitize** = divide by weight norm (what `Ideal` does)
#[derive(Clone, Copy)]
#[repr(transparent)]
pub struct Ideal<T> {
    /// The wrapped value, guaranteed to satisfy the wrapper's constraint.
    inner: T,
}

impl<T: DegenerateNormed> Ideal<T> {
    /// Creates an `Ideal<T>` by unitizing (dividing by weight norm).
    ///
    /// Returns `None` if the weight norm is zero or near-zero
    /// (e.g., ideal/infinity points).
    ///
    /// # Example
    ///
    /// ```ignore
    /// let p = Point::new(2.0, 4.0, 6.0, 2.0);
    /// let std = Ideal::try_new(p).unwrap();
    /// ```
    #[inline]
    pub fn try_new(inner: T) -> Option<Self>
    where
        T: Sized,
    {
        inner.try_unitize().map(|u| Self { inner: u })
    }

    /// Creates an `Ideal<T>` by unitizing (dividing by weight norm).
    ///
    /// # Panics
    ///
    /// Panics if the weight norm is zero or near-zero.
    #[inline]
    pub fn new_normalize(inner: T) -> Self
    where
        T: Sized,
    {
        Self::try_new(inner).expect("cannot weight-normalize ideal element")
    }
}

impl<T> Ideal<T> {
    /// Creates an `Ideal<T>` without checking or enforcing normalization.
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

    impl<T> Arbitrary for Ideal<T>
    where
        T: DegenerateNormed + Arbitrary + Debug + 'static,
    {
        type Parameters = ();
        type Strategy = BoxedStrategy<Self>;

        fn arbitrary_with(_: Self::Parameters) -> Self::Strategy {
            any::<T>()
                .prop_filter_map("weight-normalizable", |t| Ideal::try_new(t))
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
}
