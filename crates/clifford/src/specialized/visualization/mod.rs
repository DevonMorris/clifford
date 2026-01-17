//! Visualization helper types for geometric algebra primitives.
//!
//! This module provides wrapper types that disambiguate how geometric types
//! should be interpreted for visualization purposes.
//!
//! # Position vs Direction
//!
//! A [`Vector`](crate::specialized::euclidean::dim3::Vector) can represent either:
//! - A **position** (point in space)
//! - A **direction** (displacement/velocity)
//!
//! Use [`AsPosition`] and [`AsArrow`] wrappers to make intent explicit:
//!
//! ```ignore
//! use clifford::specialized::euclidean::dim3::Vector;
//! use clifford::specialized::visualization::{rerun, AsPosition, AsArrow};
//!
//! let v = Vector::new(1.0_f32, 2.0, 3.0);
//!
//! // Log as a point in space
//! rec.log("point", &rerun::Points3D::new([AsPosition(v)]))?;
//!
//! // Log as a direction arrow from origin
//! rec.log("arrow", &rerun::Arrows3D::from_vectors([AsArrow(v)]))?;
//! ```

/// Re-export of the rerun crate for convenience.
///
/// This allows users to access rerun types without adding a separate dependency.
pub use rerun_0_28 as rerun;

/// Wrapper to interpret a geometric type as a position (point in space).
///
/// Use this when logging vectors as points rather than directions.
///
/// # Example
///
/// ```ignore
/// use clifford::specialized::euclidean::dim3::Vector;
/// use clifford::specialized::visualization::AsPosition;
///
/// let point = Vector::new(1.0_f32, 2.0, 3.0);
/// rec.log("point", &rerun::Points3D::new([AsPosition(point)]))?;
/// ```
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct AsPosition<T>(pub T);

/// Wrapper to interpret a geometric type as a direction (arrow/displacement).
///
/// Use this when logging vectors as arrows rather than points.
///
/// # Example
///
/// ```ignore
/// use clifford::specialized::euclidean::dim3::Vector;
/// use clifford::specialized::visualization::AsArrow;
///
/// let direction = Vector::new(1.0_f32, 0.0, 0.0);
/// rec.log("arrow", &rerun::Arrows3D::from_vectors([AsArrow(direction)]))?;
/// ```
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct AsArrow<T>(pub T);

impl<T> AsPosition<T> {
    /// Unwraps the inner value.
    #[inline]
    pub fn into_inner(self) -> T {
        self.0
    }
}

impl<T> AsArrow<T> {
    /// Unwraps the inner value.
    #[inline]
    pub fn into_inner(self) -> T {
        self.0
    }
}

impl<T> AsRef<T> for AsPosition<T> {
    #[inline]
    fn as_ref(&self) -> &T {
        &self.0
    }
}

impl<T> AsRef<T> for AsArrow<T> {
    #[inline]
    fn as_ref(&self) -> &T {
        &self.0
    }
}

impl<T> core::ops::Deref for AsPosition<T> {
    type Target = T;

    #[inline]
    fn deref(&self) -> &Self::Target {
        &self.0
    }
}

impl<T> core::ops::Deref for AsArrow<T> {
    type Target = T;

    #[inline]
    fn deref(&self) -> &Self::Target {
        &self.0
    }
}
