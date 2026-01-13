//! Algebraic product traits for Geometric Algebra.
//!
//! This module defines traits for the various products in Geometric Algebra.
//! These traits enable method-based APIs for products, improving discoverability
//! through IDE autocomplete.
//!
//! # Products Overview
//!
//! | Trait | Symbol | Description |
//! |-------|--------|-------------|
//! | [`GeometricProduct`] | `×` | Full geometric product |
//! | [`Wedge`] | `∧` | Exterior/wedge product (grade-raising) |
//! | [`Antiwedge`] | `∨` | Regressive product (antigrade-raising) |
//! | [`Inner`] | `·` | Symmetric inner product |
//! | [`LeftContract`] | `⌋` | Left contraction |
//! | [`RightContract`] | `⌊` | Right contraction |
//! | [`BulkContract`] | `a ∨ b★` | Bulk contraction |
//! | [`WeightContract`] | `a ∨ b☆` | Weight contraction |
//! | [`Sandwich`] | `v × x × ṽ` | Sandwich product |
//!
//! # Example
//!
//! ```ignore
//! use clifford::ops::{Wedge, Sandwich};
//! use clifford::specialized::projective::dim3::{Point, Motor};
//!
//! // Create a line through two points using wedge product
//! let line = p1.wedge(&p2);
//!
//! // Transform a point using sandwich product
//! let transformed = motor.sandwich(&point);
//! ```
//!
//! # Implementation
//!
//! These traits are implemented by the code generator for specialized types.
//! The implementations delegate to optimized free functions.

/// Geometric product.
///
/// The geometric product is the fundamental product in Geometric Algebra,
/// combining the inner and outer products: `ab = a·b + a∧b` for vectors.
///
/// # Example
///
/// ```ignore
/// use clifford::ops::GeometricProduct;
///
/// let result = a.geometric(&b);
/// ```
pub trait GeometricProduct<Rhs = Self> {
    /// The output type of the geometric product.
    type Output;

    /// Computes the geometric product `self * rhs`.
    fn geometric(&self, rhs: &Rhs) -> Self::Output;
}

/// Wedge (exterior) product.
///
/// The wedge product `a ∧ b` is grade-raising: `grade(a ∧ b) = grade(a) + grade(b)`.
/// It represents the "span" of two elements - two vectors wedged give the plane they span.
///
/// # Properties
///
/// - Antisymmetric: `a ∧ b = -(b ∧ a)`
/// - Associative: `(a ∧ b) ∧ c = a ∧ (b ∧ c)`
/// - `a ∧ a = 0` for vectors
///
/// # Example
///
/// ```ignore
/// use clifford::ops::Wedge;
///
/// // Line through two points
/// let line = p1.wedge(&p2);
///
/// // Plane through point and line
/// let plane = point.wedge(&line);
/// ```
pub trait Wedge<Rhs = Self> {
    /// The output type of the wedge product.
    type Output;

    /// Computes the wedge product `self ∧ rhs`.
    fn wedge(&self, rhs: &Rhs) -> Self::Output;
}

/// Antiwedge (regressive) product.
///
/// The antiwedge product `a ∨ b` is antigrade-raising. It's the dual of the wedge product
/// and computes the "meet" of geometric elements.
///
/// # Example
///
/// ```ignore
/// use clifford::ops::Antiwedge;
///
/// // Point where two lines meet
/// let point = line1.antiwedge(&line2);
///
/// // Line where two planes meet
/// let line = plane1.antiwedge(&plane2);
/// ```
pub trait Antiwedge<Rhs = Self> {
    /// The output type of the antiwedge product.
    type Output;

    /// Computes the antiwedge product `self ∨ rhs`.
    fn antiwedge(&self, rhs: &Rhs) -> Self::Output;
}

/// Symmetric inner product (Hestenes inner product).
///
/// The inner product reduces grade: `grade(a · b) = |grade(a) - grade(b)|`.
///
/// # Example
///
/// ```ignore
/// use clifford::ops::Inner;
///
/// let scalar = vector1.inner(&vector2);
/// ```
pub trait Inner<Rhs = Self> {
    /// The output type of the inner product.
    type Output;

    /// Computes the inner product `self · rhs`.
    fn inner(&self, rhs: &Rhs) -> Self::Output;
}

/// Left contraction.
///
/// The left contraction `a ⌋ b` projects `b` onto elements orthogonal to `a`.
/// Result is grade `grade(b) - grade(a)` when `grade(a) ≤ grade(b)`, else zero.
///
/// # Example
///
/// ```ignore
/// use clifford::ops::LeftContract;
///
/// // Contract a vector with a bivector
/// let result = vector.left_contract(&bivector);
/// ```
pub trait LeftContract<Rhs = Self> {
    /// The output type of the left contraction.
    type Output;

    /// Computes the left contraction `self ⌋ rhs`.
    fn left_contract(&self, rhs: &Rhs) -> Self::Output;
}

/// Right contraction.
///
/// The right contraction `a ⌊ b` is the dual of left contraction.
/// Result is grade `grade(a) - grade(b)` when `grade(b) ≤ grade(a)`, else zero.
///
/// # Example
///
/// ```ignore
/// use clifford::ops::RightContract;
///
/// let result = bivector.right_contract(&vector);
/// ```
pub trait RightContract<Rhs = Self> {
    /// The output type of the right contraction.
    type Output;

    /// Computes the right contraction `self ⌊ rhs`.
    fn right_contract(&self, rhs: &Rhs) -> Self::Output;
}

/// Bulk contraction (`a ∨ b★`).
///
/// The bulk contraction is the antiwedge of `a` with the bulk dual of `b`.
/// Used for projections in PGA.
///
/// # Example
///
/// ```ignore
/// use clifford::ops::BulkContract;
///
/// let projected = element.bulk_contract(&other);
/// ```
pub trait BulkContract<Rhs = Self> {
    /// The output type of the bulk contraction.
    type Output;

    /// Computes the bulk contraction `self ∨ rhs★`.
    fn bulk_contract(&self, rhs: &Rhs) -> Self::Output;
}

/// Weight contraction (`a ∨ b☆`).
///
/// The weight contraction is the antiwedge of `a` with the weight dual of `b`.
/// Used for projections in PGA.
///
/// # Example
///
/// ```ignore
/// use clifford::ops::WeightContract;
///
/// let projected = element.weight_contract(&other);
/// ```
pub trait WeightContract<Rhs = Self> {
    /// The output type of the weight contraction.
    type Output;

    /// Computes the weight contraction `self ∨ rhs☆`.
    fn weight_contract(&self, rhs: &Rhs) -> Self::Output;
}

/// Bulk expansion (`a ∧ b★`).
///
/// The bulk expansion is the wedge of `a` with the bulk dual of `b`.
///
/// # Example
///
/// ```ignore
/// use clifford::ops::BulkExpand;
///
/// let result = element.bulk_expand(&other);
/// ```
pub trait BulkExpand<Rhs = Self> {
    /// The output type of the bulk expansion.
    type Output;

    /// Computes the bulk expansion `self ∧ rhs★`.
    fn bulk_expand(&self, rhs: &Rhs) -> Self::Output;
}

/// Weight expansion (`a ∧ b☆`).
///
/// The weight expansion is the wedge of `a` with the weight dual of `b`.
///
/// # Example
///
/// ```ignore
/// use clifford::ops::WeightExpand;
///
/// let result = element.weight_expand(&other);
/// ```
pub trait WeightExpand<Rhs = Self> {
    /// The output type of the weight expansion.
    type Output;

    /// Computes the weight expansion `self ∧ rhs☆`.
    fn weight_expand(&self, rhs: &Rhs) -> Self::Output;
}

/// Scalar product.
///
/// The scalar product extracts only the grade-0 (scalar) part of the geometric product.
///
/// # Example
///
/// ```ignore
/// use clifford::ops::ScalarProduct;
///
/// let s: f64 = a.scalar_product(&b);
/// ```
pub trait ScalarProduct<Rhs = Self> {
    /// The scalar type.
    type Scalar;

    /// Computes the scalar product (grade-0 part of geometric product).
    fn scalar_product(&self, rhs: &Rhs) -> Self::Scalar;
}

/// Sandwich product (`v × x × ṽ`).
///
/// The sandwich product transforms `x` by the versor `v`. This is the fundamental
/// transformation operation in GA - rotors, motors, and flectors all transform
/// elements via sandwich products.
///
/// # Example
///
/// ```ignore
/// use clifford::ops::Sandwich;
///
/// // Transform a point by a motor
/// let transformed_point = motor.sandwich(&point);
///
/// // Rotate a vector by a rotor
/// let rotated = rotor.sandwich(&vector);
/// ```
pub trait Sandwich<Operand> {
    /// The output type (usually same as Operand).
    type Output;

    /// Computes the sandwich product `self × operand × rev(self)`.
    fn sandwich(&self, operand: &Operand) -> Self::Output;
}

/// Antisandwich product (`v ⟇ x ⟇ ṽ`).
///
/// The antisandwich product is the complement of the sandwich product,
/// used for transforming elements via the geometric antiproduct.
///
/// # Example
///
/// ```ignore
/// use clifford::ops::Antisandwich;
///
/// let transformed = flector.antisandwich(&element);
/// ```
pub trait Antisandwich<Operand> {
    /// The output type.
    type Output;

    /// Computes the antisandwich product.
    fn antisandwich(&self, operand: &Operand) -> Self::Output;
}

/// Antigeometric product.
///
/// The geometric antiproduct is the complement of the geometric product:
/// `a ⟇ b = complement(complement(a) × complement(b))`.
///
/// # Example
///
/// ```ignore
/// use clifford::ops::Antigeometric;
///
/// let result = a.antigeometric(&b);
/// ```
pub trait Antigeometric<Rhs = Self> {
    /// The output type of the antigeometric product.
    type Output;

    /// Computes the antigeometric product.
    fn antigeometric(&self, rhs: &Rhs) -> Self::Output;
}

// ============================================================================
// Unary Operations
// ============================================================================

/// Reverse operation.
///
/// The reverse `ã` reverses the order of basis vectors in each blade.
/// For a blade of grade k: `rev(e₁e₂...eₖ) = eₖ...e₂e₁ = (-1)^(k(k-1)/2) * blade`.
///
/// # Example
///
/// ```ignore
/// use clifford::ops::Reverse;
///
/// let reversed = bivector.reverse();
/// ```
pub trait Reverse {
    /// Computes the reverse.
    fn reverse(&self) -> Self;
}

/// Antireverse operation.
///
/// The antireverse `a̲` (tilde below) is the antiproduct complement of reverse.
///
/// # Example
///
/// ```ignore
/// use clifford::ops::Antireverse;
///
/// let antireversed = element.antireverse();
/// ```
pub trait Antireverse {
    /// Computes the antireverse.
    fn antireverse(&self) -> Self;
}

/// Right complement operation.
///
/// The right complement `ā` (bar above) satisfies: `a ∧ ā = 𝟙` (pseudoscalar).
///
/// # Example
///
/// ```ignore
/// use clifford::ops::RightComplement;
///
/// let complement = vector.right_complement();
/// ```
pub trait RightComplement {
    /// The output type of the right complement.
    type Output;

    /// Computes the right complement.
    fn right_complement(&self) -> Self::Output;
}

/// Left complement operation.
///
/// The left complement `a̱` (bar below) satisfies: `a̱ ∧ a = 𝟙` (pseudoscalar).
///
/// # Example
///
/// ```ignore
/// use clifford::ops::LeftComplement;
///
/// let complement = vector.left_complement();
/// ```
pub trait LeftComplement {
    /// The output type of the left complement.
    type Output;

    /// Computes the left complement.
    fn left_complement(&self) -> Self::Output;
}

/// Bulk dual operation.
///
/// The bulk dual `u★` is defined as: `u★ = ũ ⋙ 1` (geometric product with
/// pseudoscalar on right). For non-degenerate metrics, this equals the Hodge dual.
///
/// # Example
///
/// ```ignore
/// use clifford::ops::BulkDual;
///
/// let dual = vector.bulk_dual();
/// ```
pub trait BulkDual {
    /// The output type of the bulk dual.
    type Output;

    /// Computes the bulk dual.
    fn bulk_dual(&self) -> Self::Output;
}

/// Weight dual (antidual) operation.
///
/// The weight dual `u☆` is defined as: `u☆ = ũ ⋗ 1` (antigeometric product
/// with antiscalar on right). For degenerate metrics (PGA), this differs from bulk dual.
///
/// # Example
///
/// ```ignore
/// use clifford::ops::WeightDual;
///
/// let antidual = vector.weight_dual();
/// ```
pub trait WeightDual {
    /// The output type of the weight dual.
    type Output;

    /// Computes the weight dual (antidual).
    fn weight_dual(&self) -> Self::Output;
}
