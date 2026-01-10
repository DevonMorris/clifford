//! Error types for PGA conversions and constraint violations.
//!
//! This module provides error types used when:
//! - Converting between clifford PGA types and nalgebra types
//! - Constructing PGA types that must satisfy geometric constraints

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

/// Error type for motor constraint violations.
///
/// Motors must satisfy the Study condition to represent valid rigid transformations.
/// See: <https://rigidgeometricalgebra.org/wiki/index.php?title=Geometric_constraint>
#[derive(Clone, Copy, Debug, PartialEq, Eq, Error)]
pub enum MotorConstraintError {
    /// Motor violates the Study condition.
    ///
    /// The Study condition requires: `s·e₀₁₂₃ - e₂₃·e₀₁ - e₃₁·e₀₂ - e₁₂·e₀₃ = 0`
    ///
    /// This ensures the motor represents a proper rigid transformation (rotation + translation).
    #[error("motor violates Study condition: s·e₀₁₂₃ ≠ v·m")]
    StudyConditionViolation,

    /// Motor has zero weight norm (degenerate).
    ///
    /// The weight norm `s² + e₂₃² + e₃₁² + e₁₂²` must be non-zero for a valid motor.
    #[error("motor has zero weight norm")]
    ZeroWeightNorm,
}

/// Error type for line constraint violations.
///
/// Lines in 3D PGA must satisfy the Plücker condition to represent valid geometric lines.
/// See: <https://rigidgeometricalgebra.org/wiki/index.php?title=Geometric_constraint>
#[derive(Clone, Copy, Debug, PartialEq, Eq, Error)]
pub enum LineConstraintError {
    /// Line violates the Plücker condition.
    ///
    /// The Plücker condition requires: `direction · moment = 0`
    /// i.e., `e₀₁·e₂₃ + e₀₂·e₃₁ + e₀₃·e₁₂ = 0`
    ///
    /// This ensures the line components describe a valid line in 3D space.
    #[error("line violates Plücker condition: direction·moment ≠ 0")]
    PluckerConditionViolation,

    /// Line has zero direction (degenerate).
    ///
    /// The direction `(e₀₁, e₀₂, e₀₃)` must be non-zero for a valid line.
    #[error("line has zero direction")]
    ZeroDirection,
}
