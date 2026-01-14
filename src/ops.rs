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
//! | [`Wedge`] | `∧` | Exterior/wedge product (grade-raising) |
//! | [`Antiwedge`] | `∨` | Regressive product (antigrade-raising) |
//! | [`Join`] | `∨` | Alias for Wedge (PGA terminology) |
//! | [`Meet`] | `∧` | Alias for Antiwedge (PGA terminology) |
//! | [`Dot`] | `•` | Metric dot product (same-grade only) |
//! | [`Antidot`] | `⊚` | Metric antidot product |
//! | [`LeftContract`] | `⌋` | Left contraction |
//! | [`RightContract`] | `⌊` | Right contraction |
//! | [`BulkContract`] | `a ∨ b★` | Bulk contraction |
//! | [`WeightContract`] | `a ∨ b☆` | Weight contraction |
//! | [`Sandwich`] | `v × x × ṽ` | Sandwich product |
//! | [`Versor`] | `a ⟇ b` | Versor composition |
//!
//! # Type Safety
//!
//! All product traits in this module are **type-safe**: the output type is
//! determined by the input types and the product kind. This is achieved by:
//!
//! - **Grade-preserving products** (Sandwich, Antisandwich) return the same type
//! - **Grade-changing products** (Wedge, Antiwedge) have deterministic output grades
//! - **Scalar products** (Dot, Antidot) always return the scalar type
//!
//! The geometric product is **not** provided as a trait because it produces
//! mixed-grade outputs that cannot be represented by a single specialized type.
//! Use the grade-specific products instead.
//!
//! # Example
//!
//! ```ignore
//! use clifford::ops::{Wedge, Sandwich, Dot};
//! use clifford::specialized::projective::dim3::{Point, Motor};
//!
//! // Create a line through two points using wedge product
//! let line = p1.wedge(&p2);
//!
//! // Transform a point using sandwich product
//! let transformed = motor.sandwich(&point);
//!
//! // Compute dot product for angles/magnitudes
//! let cos_angle = v1.dot(&v2);
//! ```
//!
//! # Implementation
//!
//! These traits are implemented by the code generator for specialized types.
//! The implementations delegate to optimized free functions.

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

/// Dot product (metric inner product).
///
/// The dot product `a • b` is non-zero only when `grade(a) = grade(b)`.
/// It measures the "alignment" between same-grade elements using the metric.
///
/// Defined as: `a • b = aᵀ G b` where G is the metric exomorphism matrix.
///
/// # Properties
///
/// - Only same-grade elements produce non-zero results
/// - Symmetric: `a • b = b • a`
/// - Returns a scalar
///
/// # Example
///
/// ```ignore
/// use clifford::ops::Dot;
///
/// let cos_angle: f64 = vector1.dot(&vector2);
/// let bivector_magnitude_sq: f64 = bivector.dot(&bivector);
/// ```
///
/// # Reference
///
/// [RGA Dot Products](https://rigidgeometricalgebra.org/wiki/index.php?title=Dot_products)
pub trait Dot<Rhs = Self> {
    /// The scalar type.
    type Scalar;

    /// Computes the dot product `self • rhs`.
    fn dot(&self, rhs: &Rhs) -> Self::Scalar;
}

/// Antidot product (metric antiproduct inner product).
///
/// The antidot product `a ⊚ b` is the De Morgan dual of the dot product:
/// `a ⊚ b = ā • b̄` (complement dot complement).
///
/// Defined as: `a ⊚ b = aᵀ 𝔾 b` where 𝔾 is the metric antiexomorphism matrix.
///
/// # Properties
///
/// - Only same-antigrade elements produce non-zero results
/// - Symmetric
/// - Returns a scalar
///
/// # Example
///
/// ```ignore
/// use clifford::ops::Antidot;
///
/// let result: f64 = plane1.antidot(&plane2);
/// ```
///
/// # Reference
///
/// [RGA Dot Products](https://rigidgeometricalgebra.org/wiki/index.php?title=Dot_products)
pub trait Antidot<Rhs = Self> {
    /// The scalar type.
    type Scalar;

    /// Computes the antidot product `self ⊚ rhs`.
    fn antidot(&self, rhs: &Rhs) -> Self::Scalar;
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

// ============================================================================
// Versor Operations
// ============================================================================

/// Versor composition.
///
/// Versors are elements that represent transformations. This trait enables
/// composing versors where the result type may differ from the input types.
///
/// In Geometric Algebra, versors include:
///
/// - **Rotors**: Represent rotations (Euclidean and PGA)
/// - **Motors**: Represent rigid motions (rotation + translation) in PGA
/// - **Flectors**: Represent reflections and glide reflections in PGA
///
/// The `compose` method combines two versors into a single versor representing
/// the sequential application of both transformations: `(a.compose(b)).sandwich(x)`
/// is equivalent to `a.sandwich(b.sandwich(x))`.
///
/// # Mathematical Definition
///
/// For PGA versors, composition is the geometric product (which uses the
/// antigeometric structure internally). The output type depends on the inputs:
///
/// - Motor × Motor → Motor (even × even = even)
/// - Flector × Flector → Motor (odd × odd = even)
/// - Motor × Flector → Flector (even × odd = odd)
/// - Flector × Motor → Flector (odd × even = odd)
///
/// # Properties
///
/// - Associative: `a.compose(b).compose(c) = a.compose(b.compose(c))`
/// - Identity: The scalar `1` acts as identity
///
/// # Example
///
/// ```ignore
/// use clifford::ops::Versor;
///
/// // Compose two motors (returns Motor)
/// let combined_motion: Motor<f64> = motor1.compose(&motor2);
///
/// // Compose two flectors (returns Motor)
/// let motor: Motor<f64> = flector1.compose(&flector2);
///
/// // Compose motor with flector (returns Flector)
/// let flector: Flector<f64> = motor.compose(&flector);
/// ```
///
/// # Reference
///
/// [RGA Motors](https://rigidgeometricalgebra.org/wiki/index.php?title=Motor)
pub trait Versor<Rhs = Self> {
    /// The output type of the composition.
    type Output;

    /// Composes this versor with another, returning a versor representing
    /// the sequential application of both transformations.
    fn compose(&self, other: &Rhs) -> Self::Output;
}

// ============================================================================
// Join and Meet (aliases for Wedge and Antiwedge)
// ============================================================================

/// Join product (alias for wedge/exterior product).
///
/// The join `a ∨ b` is the wedge product `a ∧ b`. It "joins" two geometric
/// elements to create a higher-dimensional element spanning both.
///
/// In PGA terminology:
/// - Two points joined give the line through them
/// - A point and line joined give the plane containing both
/// - A point, line, and plane joined give the volume
///
/// # Relationship to Wedge
///
/// This trait is automatically implemented for any type implementing [`Wedge`].
/// It provides a more intuitive name for PGA operations.
///
/// # Example
///
/// ```ignore
/// use clifford::ops::Join;
///
/// // Line through two points
/// let line = p1.join(&p2);
///
/// // Plane through point and line
/// let plane = point.join(&line);
/// ```
///
/// # Reference
///
/// [RGA Wedge Product](https://rigidgeometricalgebra.org/wiki/index.php?title=Wedge_product)
pub trait Join<Rhs = Self> {
    /// The output type of the join.
    type Output;

    /// Computes the join `self ∨ rhs` (wedge product).
    fn join(&self, rhs: &Rhs) -> Self::Output;
}

/// Blanket implementation of Join for any type implementing Wedge.
impl<T, Rhs> Join<Rhs> for T
where
    T: Wedge<Rhs>,
{
    type Output = <T as Wedge<Rhs>>::Output;

    #[inline]
    fn join(&self, rhs: &Rhs) -> Self::Output {
        self.wedge(rhs)
    }
}

/// Meet product (alias for antiwedge/regressive product).
///
/// The meet `a ∧ b` is the antiwedge product `a ∨ b`. It finds the
/// "intersection" of two geometric elements.
///
/// In PGA terminology:
/// - Two planes meet at a line
/// - A plane and line meet at a point
/// - Two lines meet at a point (if they intersect)
///
/// # Relationship to Antiwedge
///
/// This trait is automatically implemented for any type implementing [`Antiwedge`].
/// It provides a more intuitive name for PGA operations.
///
/// # Example
///
/// ```ignore
/// use clifford::ops::Meet;
///
/// // Line where two planes intersect
/// let line = plane1.meet(&plane2);
///
/// // Point where plane and line intersect
/// let point = plane.meet(&line);
/// ```
///
/// # Reference
///
/// [RGA Antiwedge Product](https://rigidgeometricalgebra.org/wiki/index.php?title=Antiwedge_product)
pub trait Meet<Rhs = Self> {
    /// The output type of the meet.
    type Output;

    /// Computes the meet `self ∧ rhs` (antiwedge product).
    fn meet(&self, rhs: &Rhs) -> Self::Output;
}

/// Blanket implementation of Meet for any type implementing Antiwedge.
impl<T, Rhs> Meet<Rhs> for T
where
    T: Antiwedge<Rhs>,
{
    type Output = <T as Antiwedge<Rhs>>::Output;

    #[inline]
    fn meet(&self, rhs: &Rhs) -> Self::Output {
        self.antiwedge(rhs)
    }
}

/// Transform operation (alias for antisandwich product).
///
/// The transform `versor.transform(element)` is the antisandwich product
/// `versor ⊛ element ⊛ rev(versor)` where ⊛ is the geometric antiproduct.
///
/// In PGA, transformations (rotations, translations, reflections) are applied
/// via the antisandwich product. This trait provides a more intuitive name
/// for this operation.
///
/// # Implementation
///
/// Transform implementations are generated by codegen based on the algebra's signature:
/// - **Degenerate** (has zero elements like e0² = 0): delegates to [`Antisandwich`]
/// - **Non-degenerate** (no zero elements): delegates to [`Sandwich`]
///
/// This distinction exists because degenerate metrics cause sandwich terms to vanish,
/// requiring the antisandwich (which uses complements) to preserve information.
///
/// # Example
///
/// ```ignore
/// use clifford::ops::Transform;
///
/// // Apply motor (rotation + translation) to a point
/// let transformed_point = motor.transform(&point);
///
/// // Apply motor to a line
/// let transformed_line = motor.transform(&line);
///
/// // Apply motor to a plane
/// let transformed_plane = motor.transform(&plane);
/// ```
///
/// # Reference
///
/// [RGA Transformations](https://rigidgeometricalgebra.org/wiki/index.php?title=Transformations)
pub trait Transform<Operand> {
    /// The output type (same as operand for grade-preserving transforms).
    type Output;

    /// Transforms the operand using the appropriate sandwich product for this algebra.
    fn transform(&self, operand: &Operand) -> Self::Output;
}

// Transform implementations are generated by codegen based on signature:
// - Degenerate signature (has zero elements like e0² = 0): uses antisandwich
// - Non-degenerate signature (no zero elements): uses sandwich
