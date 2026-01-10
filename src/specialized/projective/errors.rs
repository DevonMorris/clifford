//! Error types for PGA nalgebra conversions.
//!
//! This module provides shared error types used when converting between
//! clifford PGA types and nalgebra types.

use thiserror::Error;

/// Error type for point conversions to nalgebra types.
///
/// This error occurs when trying to convert a PGA point to a nalgebra point,
/// but the conversion cannot succeed because the PGA point is at infinity.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Error)]
pub enum PointConversionError {
    /// The point is at infinity (homogeneous weight ≈ 0).
    ///
    /// Ideal points (points at infinity) cannot be represented as
    /// nalgebra `Point2` or `Point3` types, which require finite coordinates.
    #[error("point is at infinity (w ≈ 0)")]
    PointAtInfinity,
}
