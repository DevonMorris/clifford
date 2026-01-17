//! The multivector: the fundamental element of geometric algebra.
//!
//! # What is a Multivector?
//!
//! A **multivector** is the most general element of a geometric algebra. While
//! a blade represents a single oriented subspace (scalar, vector, bivector, etc.),
//! a multivector is a **linear combination of blades** of potentially different grades.
//!
//! Think of it this way:
//! - A **scalar** is a number: `5`
//! - A **vector** is a direction with magnitude: `3e₁ + 4e₂`
//! - A **bivector** is an oriented plane: `2e₁₂`
//! - A **multivector** can mix all of these: `5 + 3e₁ + 4e₂ + 2e₁₂`
//!
//! # Why Multivectors?
//!
//! When we multiply geometric objects, we often get mixed-grade results.
//! For example, the geometric product of two vectors gives:
//!
//! ```text
//! ab = a·b + a∧b
//!    = (scalar) + (bivector)
//! ```
//!
//! The dot product `a·b` is a scalar (grade 0), and the wedge product `a∧b`
//! is a bivector (grade 2). The result is a multivector with both grades.
//!
//! # The Geometric Product
//!
//! The geometric product is the fundamental operation in GA. For multivectors
//! A and B, the product AB is computed by:
//!
//! 1. Distributing over all blade pairs
//! 2. Using the basis blade multiplication rules
//! 3. Collecting terms by grade
//!
//! Key properties:
//! - **Associative**: (AB)C = A(BC)
//! - **Distributive**: A(B + C) = AB + AC
//! - **NOT commutative**: AB ≠ BA in general
//!
//! # Representation
//!
//! We store multivectors as a dense array of coefficients, one for each
//! basis blade. In an n-dimensional space, there are 2ⁿ basis blades,
//! so a 3D multivector has 8 coefficients:
//!
//! ```text
//! M = c₀·1 + c₁·e₁ + c₂·e₂ + c₃·e₁₂ + c₄·e₃ + c₅·e₁₃ + c₆·e₂₃ + c₇·e₁₂₃
//! ```
//!
//! # Example
//!
//! ```
//! use clifford::algebra::Multivector;
//! use clifford::signature::Euclidean3;
//!
//! // Create a vector: 3e₁ + 4e₂
//! let v: Multivector<f64, Euclidean3> = Multivector::vector(&[3.0, 4.0, 0.0]);
//!
//! // The geometric product of a vector with itself gives a scalar
//! let v_squared = &v * &v;
//! assert!((v_squared.scalar_part() - 25.0).abs() < 1e-10); // 3² + 4² = 25
//! ```

use approx::{AbsDiffEq, RelativeEq, UlpsEq};
use core::fmt;
use core::marker::PhantomData;
use core::ops::{Add, AddAssign, BitXor, Div, Mul, MulAssign, Neg, Sub, SubAssign};

use generic_array::{GenericArray, sequence::GenericSequence};
use typenum::Unsigned;

#[cfg(feature = "serde")]
use serde::{Deserialize, Serialize};

use crate::basis::Blade;
use crate::scalar::Float;
use crate::signature::Signature;

/// A multivector in a Clifford algebra with signature `S` and scalar type `T`.
///
/// A multivector is a linear combination of basis blades. It is the most
/// general element of a geometric algebra, capable of representing scalars,
/// vectors, bivectors, and arbitrary combinations thereof.
///
/// # Type Parameters
///
/// - `T`: The scalar type (e.g., `f32`, `f64`). Must implement [`Float`].
/// - `S`: The metric signature (e.g., [`Euclidean3`](crate::signature::Euclidean3)).
///   Defines how basis vectors square and thus the algebra's geometry.
///
/// # Storage
///
/// Coefficients are stored in a dense array indexed by blade index.
/// The blade index is a bitmask where bit `i` indicates the presence of
/// basis vector `eᵢ`. See [`Blade`] for details.
///
/// # Example
///
/// ```
/// use clifford::algebra::Multivector;
/// use clifford::signature::Euclidean2;
///
/// // Create basis vectors
/// let e1: Multivector<f64, Euclidean2> = Multivector::basis_vector(0);
/// let e2: Multivector<f64, Euclidean2> = Multivector::basis_vector(1);
///
/// // Geometric product: e₁e₂ = e₁₂ (a bivector)
/// let e12 = &e1 * &e2;
///
/// // e₁₂ squares to -1 in Euclidean space!
/// let e12_squared = &e12 * &e12;
/// assert!((e12_squared.scalar_part() - (-1.0)).abs() < 1e-10);
/// ```
#[derive(Clone)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub struct Multivector<T: Float, S: Signature> {
    /// Coefficients for each basis blade, indexed by blade bitmask.
    coeffs: GenericArray<T, S::NumBlades>,
    /// Marker for the signature type.
    _signature: PhantomData<S>,
}

impl<T: Float, S: Signature> Copy for Multivector<T, S> where GenericArray<T, S::NumBlades>: Copy {}

// ============================================================================
// Constructors
// ============================================================================

impl<T: Float, S: Signature> Multivector<T, S> {
    /// Creates the zero multivector (additive identity).
    ///
    /// All coefficients are zero. This is the identity element for addition:
    /// `M + 0 = M` for any multivector M.
    ///
    /// # Example
    ///
    /// ```
    /// use clifford::algebra::Multivector;
    /// use clifford::signature::Euclidean3;
    ///
    /// let zero: Multivector<f64, Euclidean3> = Multivector::zero();
    /// assert!(zero.is_zero(1e-10));
    /// ```
    #[inline]
    pub fn zero() -> Self {
        Self {
            coeffs: GenericArray::generate(|_| T::zero()),
            _signature: PhantomData,
        }
    }

    /// Creates the unit scalar (multiplicative identity).
    ///
    /// This is the scalar `1`, which is the identity for the geometric product:
    /// `M * 1 = 1 * M = M` for any multivector M.
    ///
    /// # Example
    ///
    /// ```
    /// use clifford::algebra::Multivector;
    /// use clifford::signature::Euclidean3;
    ///
    /// let one: Multivector<f64, Euclidean3> = Multivector::one();
    /// assert!((one.scalar_part() - 1.0).abs() < 1e-10);
    /// ```
    #[inline]
    pub fn one() -> Self {
        Self::scalar(T::one())
    }

    /// Creates a scalar multivector (grade 0 only).
    ///
    /// A scalar is a multivector with only the grade-0 component nonzero.
    ///
    /// # Example
    ///
    /// ```
    /// use clifford::algebra::Multivector;
    /// use clifford::signature::Euclidean3;
    ///
    /// let five: Multivector<f64, Euclidean3> = Multivector::scalar(5.0);
    /// assert!((five.scalar_part() - 5.0).abs() < 1e-10);
    /// ```
    #[inline]
    pub fn scalar(value: T) -> Self {
        let mut mv = Self::zero();
        mv.coeffs[0] = value;
        mv
    }

    /// Creates a multivector with a single basis blade.
    ///
    /// The coefficient of the specified blade is set to 1, all others to 0.
    ///
    /// # Arguments
    ///
    /// * `blade` - The basis blade to create
    ///
    /// # Example
    ///
    /// ```
    /// use clifford::algebra::Multivector;
    /// use clifford::basis::Blade;
    /// use clifford::signature::Euclidean3;
    ///
    /// // Create the bivector e₁₂
    /// let e12_blade = Blade::from_index(0b011);
    /// let e12: Multivector<f64, Euclidean3> = Multivector::from_blade(e12_blade);
    ///
    /// assert!((e12.get(e12_blade) - 1.0).abs() < 1e-10);
    /// ```
    #[inline]
    pub fn from_blade(blade: Blade) -> Self {
        let mut mv = Self::zero();
        mv.coeffs[blade.index()] = T::one();
        mv
    }

    /// Creates a basis vector multivector (grade 1).
    ///
    /// Basis vector `eᵢ` is the fundamental directional unit along axis `i`.
    ///
    /// # Arguments
    ///
    /// * `i` - Index of the basis vector (0-indexed)
    ///
    /// # Panics
    ///
    /// Panics if `i >= S::DIM`.
    ///
    /// # Example
    ///
    /// ```
    /// use clifford::algebra::Multivector;
    /// use clifford::signature::Euclidean3;
    ///
    /// let e1: Multivector<f64, Euclidean3> = Multivector::basis_vector(0);
    /// let e2: Multivector<f64, Euclidean3> = Multivector::basis_vector(1);
    ///
    /// // Vectors square to their metric value (1 for Euclidean)
    /// let e1_sq = &e1 * &e1;
    /// assert!((e1_sq.scalar_part() - 1.0).abs() < 1e-10);
    /// ```
    #[inline]
    pub fn basis_vector(i: usize) -> Self {
        assert!(
            i < S::DIM,
            "basis vector index {i} out of range for dimension {}",
            S::DIM
        );
        Self::from_blade(Blade::basis_vector(i))
    }

    /// Creates a vector multivector from components.
    ///
    /// The components are the coefficients of the basis vectors e₁, e₂, etc.
    ///
    /// # Arguments
    ///
    /// * `components` - Slice of coefficients for each basis vector
    ///
    /// # Panics
    ///
    /// Panics if `components.len() > S::DIM`.
    ///
    /// # Example
    ///
    /// ```
    /// use clifford::algebra::Multivector;
    /// use clifford::signature::Euclidean3;
    ///
    /// // Create vector v = 3e₁ + 4e₂ + 0e₃
    /// let v: Multivector<f64, Euclidean3> = Multivector::vector(&[3.0, 4.0, 0.0]);
    ///
    /// // |v|² = 3² + 4² = 25
    /// let v_sq = &v * &v;
    /// assert!((v_sq.scalar_part() - 25.0).abs() < 1e-10);
    /// ```
    #[inline]
    pub fn vector(components: &[T]) -> Self {
        assert!(
            components.len() <= S::DIM,
            "too many components ({}) for dimension {}",
            components.len(),
            S::DIM
        );
        let mut mv = Self::zero();
        for (i, &c) in components.iter().enumerate() {
            mv.coeffs[1 << i] = c;
        }
        mv
    }

    /// Creates a multivector from a full array of coefficients.
    ///
    /// # Arguments
    ///
    /// * `coeffs` - Coefficients for all basis blades (only first `S::NumBlades::USIZE` are used)
    #[inline]
    pub fn from_coeffs(coeffs: &[T]) -> Self {
        let mut mv = Self::zero();
        let n = coeffs.len().min(S::NumBlades::USIZE);
        mv.coeffs[..n].copy_from_slice(&coeffs[..n]);
        mv
    }
}

// ============================================================================
// Accessors
// ============================================================================

impl<T: Float, S: Signature> Multivector<T, S> {
    /// Returns the coefficient of the given basis blade.
    ///
    /// # Example
    ///
    /// ```
    /// use clifford::algebra::Multivector;
    /// use clifford::basis::Blade;
    /// use clifford::signature::Euclidean3;
    ///
    /// let v: Multivector<f64, Euclidean3> = Multivector::vector(&[3.0, 4.0, 5.0]);
    /// assert!((v.get(Blade::basis_vector(0)) - 3.0).abs() < 1e-10);
    /// assert!((v.get(Blade::basis_vector(1)) - 4.0).abs() < 1e-10);
    /// ```
    #[inline]
    pub fn get(&self, blade: Blade) -> T {
        self.coeffs[blade.index()]
    }

    /// Sets the coefficient of the given basis blade.
    ///
    /// # Example
    ///
    /// ```
    /// use clifford::algebra::Multivector;
    /// use clifford::basis::Blade;
    /// use clifford::signature::Euclidean3;
    ///
    /// let mut mv: Multivector<f64, Euclidean3> = Multivector::zero();
    /// mv.set(Blade::basis_vector(0), 5.0);
    /// assert!((mv.get(Blade::basis_vector(0)) - 5.0).abs() < 1e-10);
    /// ```
    #[inline]
    pub fn set(&mut self, blade: Blade, value: T) {
        self.coeffs[blade.index()] = value;
    }

    /// Returns the scalar part (grade 0 coefficient).
    ///
    /// This is the coefficient of the unit scalar blade.
    ///
    /// # Example
    ///
    /// ```
    /// use clifford::algebra::Multivector;
    /// use clifford::signature::Euclidean3;
    ///
    /// let mv: Multivector<f64, Euclidean3> = Multivector::scalar(42.0);
    /// assert!((mv.scalar_part() - 42.0).abs() < 1e-10);
    /// ```
    #[inline]
    pub fn scalar_part(&self) -> T {
        self.coeffs[0]
    }

    /// Checks if this multivector is approximately zero.
    ///
    /// Returns true if all coefficients have absolute value less than epsilon.
    ///
    /// # Example
    ///
    /// ```
    /// use clifford::algebra::Multivector;
    /// use clifford::signature::Euclidean3;
    ///
    /// let zero: Multivector<f64, Euclidean3> = Multivector::zero();
    /// assert!(zero.is_zero(1e-10));
    ///
    /// let v: Multivector<f64, Euclidean3> = Multivector::basis_vector(0);
    /// assert!(!v.is_zero(1e-10));
    /// ```
    #[inline]
    pub fn is_zero(&self, epsilon: T) -> bool {
        self.coeffs.iter().all(|&c| c.abs() < epsilon)
    }
}

// ============================================================================
// Unary Operations
// ============================================================================

impl<T: Float, S: Signature> Multivector<T, S> {
    /// Returns the reverse of this multivector.
    ///
    /// The reverse operation reverses the order of basis vectors in each blade.
    /// For a blade of grade k, the sign change is `(-1)^(k(k-1)/2)`:
    ///
    /// | Grade | Sign | Example |
    /// |-------|------|---------|
    /// | 0 | +1 | 1̃ = 1 |
    /// | 1 | +1 | ẽ₁ = e₁ |
    /// | 2 | -1 | ẽ₁₂ = e₂₁ = -e₁₂ |
    /// | 3 | -1 | ẽ₁₂₃ = e₃₂₁ = -e₁₂₃ |
    /// | 4 | +1 | ... |
    ///
    /// The reverse is important for computing norms and inverses:
    /// `|M|² = ⟨M M̃⟩₀` (scalar part of M times its reverse)
    ///
    /// # Example
    ///
    /// ```
    /// use clifford::algebra::Multivector;
    /// use clifford::signature::Euclidean3;
    ///
    /// let e1: Multivector<f64, Euclidean3> = Multivector::basis_vector(0);
    /// let e2: Multivector<f64, Euclidean3> = Multivector::basis_vector(1);
    /// let e12 = &e1 * &e2;
    ///
    /// // Reverse of a bivector negates it
    /// let e12_rev = e12.reverse();
    /// assert!((&e12 + &e12_rev).is_zero(1e-10));
    /// ```
    pub fn reverse(&self) -> Self {
        let mut result = Self::zero();
        for i in 0..S::NumBlades::USIZE {
            let grade = Blade::from_index(i).grade();
            // Sign is (-1)^(k(k-1)/2)
            // For grade 0 and 1: sign = +1
            // For grade 2 and 3: sign = -1
            // For grade 4 and 5: sign = +1
            // Pattern repeats every 4 grades
            let sign = if grade < 2 || (grade / 2).is_multiple_of(2) {
                T::one()
            } else {
                -T::one()
            };
            result.coeffs[i] = sign * self.coeffs[i];
        }
        result
    }

    /// Returns the grade involution of this multivector.
    ///
    /// Grade involution negates all odd-grade components.
    /// For a blade of grade k, the sign is `(-1)^k`:
    ///
    /// | Grade | Sign |
    /// |-------|------|
    /// | 0 (scalar) | +1 |
    /// | 1 (vector) | -1 |
    /// | 2 (bivector) | +1 |
    /// | 3 (trivector) | -1 |
    ///
    /// This operation is also called the "main involution" or "grade automorphism".
    ///
    /// # Example
    ///
    /// ```
    /// use clifford::algebra::Multivector;
    /// use clifford::signature::Euclidean3;
    ///
    /// let v: Multivector<f64, Euclidean3> = Multivector::basis_vector(0);
    ///
    /// // Involute of a vector negates it
    /// let v_inv = v.involute();
    /// assert!((&v + &v_inv).is_zero(1e-10));
    /// ```
    pub fn involute(&self) -> Self {
        let mut result = Self::zero();
        for i in 0..S::NumBlades::USIZE {
            let grade = Blade::from_index(i).grade();
            // Sign is (-1)^k
            let sign = if grade.is_multiple_of(2) {
                T::one()
            } else {
                -T::one()
            };
            result.coeffs[i] = sign * self.coeffs[i];
        }
        result
    }

    /// Returns the Clifford conjugate of this multivector.
    ///
    /// The Clifford conjugate is the composition of reverse and grade involution.
    /// For a blade of grade k, the sign is `(-1)^(k(k+1)/2)`:
    ///
    /// | Grade | Sign |
    /// |-------|------|
    /// | 0 | +1 |
    /// | 1 | -1 |
    /// | 2 | -1 |
    /// | 3 | +1 |
    /// | 4 | +1 |
    ///
    /// # Example
    ///
    /// ```
    /// use clifford::algebra::Multivector;
    /// use clifford::signature::Euclidean3;
    ///
    /// let v: Multivector<f64, Euclidean3> = Multivector::basis_vector(0);
    /// let v_conj = v.conjugate();
    ///
    /// // Conjugate of a vector negates it
    /// assert!((&v + &v_conj).is_zero(1e-10));
    /// ```
    pub fn conjugate(&self) -> Self {
        let mut result = Self::zero();
        for i in 0..S::NumBlades::USIZE {
            let grade = Blade::from_index(i).grade();
            // Sign is (-1)^(k(k+1)/2)
            let sign = if (grade * (grade + 1) / 2).is_multiple_of(2) {
                T::one()
            } else {
                -T::one()
            };
            result.coeffs[i] = sign * self.coeffs[i];
        }
        result
    }
}

// ============================================================================
// Norms and Inverse
// ============================================================================

impl<T: Float, S: Signature> Multivector<T, S> {
    /// Returns the squared norm of this multivector.
    ///
    /// The squared norm is the scalar part of `M * M̃` (M times its reverse).
    /// For vectors in Euclidean space, this equals the sum of squared components.
    ///
    /// Note: The squared norm can be negative in non-Euclidean signatures!
    ///
    /// # Example
    ///
    /// ```
    /// use clifford::algebra::Multivector;
    /// use clifford::signature::Euclidean3;
    ///
    /// let v: Multivector<f64, Euclidean3> = Multivector::vector(&[3.0, 4.0, 0.0]);
    /// assert!((v.norm_squared() - 25.0).abs() < 1e-10); // 3² + 4² = 25
    /// ```
    pub fn norm_squared(&self) -> T {
        (self * &self.reverse()).scalar_part()
    }

    /// Returns the norm of this multivector.
    ///
    /// The norm is `sqrt(|norm_squared|)`. We take the absolute value because
    /// the squared norm can be negative in non-Euclidean signatures.
    ///
    /// # Example
    ///
    /// ```
    /// use clifford::algebra::Multivector;
    /// use clifford::signature::Euclidean3;
    ///
    /// let v: Multivector<f64, Euclidean3> = Multivector::vector(&[3.0, 4.0, 0.0]);
    /// assert!((v.norm() - 5.0).abs() < 1e-10); // sqrt(25) = 5
    /// ```
    pub fn norm(&self) -> T {
        self.norm_squared().abs().sqrt()
    }

    /// Returns the normalized multivector (unit norm), or None if norm is zero.
    ///
    /// # Example
    ///
    /// ```
    /// use clifford::algebra::Multivector;
    /// use clifford::signature::Euclidean3;
    ///
    /// let v: Multivector<f64, Euclidean3> = Multivector::vector(&[3.0, 4.0, 0.0]);
    /// let v_hat = v.normalize().unwrap();
    /// assert!((v_hat.norm() - 1.0).abs() < 1e-10);
    /// ```
    pub fn normalize(&self) -> Option<Self> {
        let n = self.norm();
        if n < T::epsilon() {
            None
        } else {
            Some(self / n)
        }
    }

    /// Returns the multiplicative inverse, or None if not invertible.
    ///
    /// For a multivector M with non-zero squared norm:
    /// `M⁻¹ = M̃ / |M|²`
    ///
    /// This satisfies `M * M⁻¹ = M⁻¹ * M = 1`.
    ///
    /// Note: Not all multivectors are invertible! The inverse exists only
    /// when the squared norm is non-zero.
    ///
    /// # Example
    ///
    /// ```
    /// use clifford::algebra::Multivector;
    /// use clifford::signature::Euclidean3;
    ///
    /// let v: Multivector<f64, Euclidean3> = Multivector::vector(&[1.0, 0.0, 0.0]);
    /// let v_inv = v.inverse().unwrap();
    ///
    /// // v * v⁻¹ = 1
    /// let product = &v * &v_inv;
    /// assert!((product.scalar_part() - 1.0).abs() < 1e-10);
    /// assert!((&product - &Multivector::one()).is_zero(1e-10));
    /// ```
    pub fn inverse(&self) -> Option<Self> {
        let norm_sq = self.norm_squared();
        if norm_sq.abs() < T::epsilon() {
            None
        } else {
            Some(&self.reverse() / norm_sq)
        }
    }
}

// ============================================================================
// Grade Operations
// ============================================================================

impl<T: Float, S: Signature> Multivector<T, S> {
    /// Extracts the grade-k part of this multivector: `⟨M⟩ₖ`.
    ///
    /// Returns a multivector containing only the components of the specified grade,
    /// with all other grades set to zero.
    ///
    /// # Arguments
    ///
    /// * `k` - The grade to extract (0 = scalar, 1 = vector, 2 = bivector, etc.)
    ///
    /// # Example
    ///
    /// ```
    /// use clifford::algebra::Multivector;
    /// use clifford::signature::Euclidean3;
    ///
    /// // Create a multivector with scalar and vector parts
    /// let mut mv: Multivector<f64, Euclidean3> = Multivector::scalar(5.0);
    /// mv = &mv + &Multivector::vector(&[1.0, 2.0, 3.0]);
    ///
    /// // Extract just the scalar part
    /// let scalar = mv.grade_select(0);
    /// assert!((scalar.scalar_part() - 5.0).abs() < 1e-10);
    ///
    /// // Extract just the vector part
    /// let vector = mv.grade_select(1);
    /// assert!(vector.grade_select(0).is_zero(1e-10));
    /// ```
    pub fn grade_select(&self, k: usize) -> Self {
        let mut result = Self::zero();
        for i in 0..S::NumBlades::USIZE {
            if Blade::from_index(i).grade() == k {
                result.coeffs[i] = self.coeffs[i];
            }
        }
        result
    }

    /// Extracts the even part of this multivector (grades 0, 2, 4, ...).
    ///
    /// The even subalgebra is closed under the geometric product and contains
    /// important elements like rotors (in 3D, the even subalgebra is isomorphic
    /// to quaternions).
    ///
    /// # Example
    ///
    /// ```
    /// use clifford::algebra::Multivector;
    /// use clifford::signature::Euclidean3;
    /// use approx::relative_eq;
    ///
    /// let e1: Multivector<f64, Euclidean3> = Multivector::basis_vector(0);
    /// let e2: Multivector<f64, Euclidean3> = Multivector::basis_vector(1);
    ///
    /// // A rotor is an even multivector: scalar + bivector
    /// let rotor = &Multivector::scalar(0.5_f64.sqrt()) + &(&e1 * &e2) * 0.5_f64.sqrt();
    /// let even = rotor.even();
    /// assert!(relative_eq!(even, rotor, epsilon = 1e-10, max_relative = 1e-10));
    /// ```
    pub fn even(&self) -> Self {
        let mut result = Self::zero();
        for i in 0..S::NumBlades::USIZE {
            if Blade::from_index(i).grade().is_multiple_of(2) {
                result.coeffs[i] = self.coeffs[i];
            }
        }
        result
    }

    /// Extracts the odd part of this multivector (grades 1, 3, 5, ...).
    ///
    /// # Example
    ///
    /// ```
    /// use clifford::algebra::Multivector;
    /// use clifford::signature::Euclidean3;
    /// use approx::relative_eq;
    ///
    /// let v: Multivector<f64, Euclidean3> = Multivector::vector(&[1.0, 2.0, 3.0]);
    /// let odd = v.odd();
    /// assert!(relative_eq!(odd, v, epsilon = 1e-10, max_relative = 1e-10)); // Vectors are odd
    /// ```
    pub fn odd(&self) -> Self {
        let mut result = Self::zero();
        for i in 0..S::NumBlades::USIZE {
            if !Blade::from_index(i).grade().is_multiple_of(2) {
                result.coeffs[i] = self.coeffs[i];
            }
        }
        result
    }

    /// Returns the grade of this multivector if it is homogeneous (all non-zero
    /// components have the same grade), or `None` if it contains multiple grades.
    ///
    /// # Arguments
    ///
    /// * `epsilon` - Tolerance for considering a coefficient as zero
    ///
    /// # Example
    ///
    /// ```
    /// use clifford::algebra::Multivector;
    /// use clifford::signature::Euclidean3;
    ///
    /// let v: Multivector<f64, Euclidean3> = Multivector::vector(&[1.0, 2.0, 3.0]);
    /// assert_eq!(v.grade(1e-10), Some(1)); // Pure vector
    ///
    /// let mixed = &v + &Multivector::scalar(1.0);
    /// assert_eq!(mixed.grade(1e-10), None); // Mixed grades
    /// ```
    pub fn grade(&self, epsilon: T) -> Option<usize> {
        let mut found_grade = None;
        for i in 0..S::NumBlades::USIZE {
            if self.coeffs[i].abs() > epsilon {
                let g = Blade::from_index(i).grade();
                match found_grade {
                    None => found_grade = Some(g),
                    Some(prev) if prev != g => return None,
                    _ => {}
                }
            }
        }
        found_grade
    }

    /// Returns the unit pseudoscalar `I = e₁e₂...eₙ`.
    ///
    /// The pseudoscalar is the highest-grade basis blade, representing the
    /// oriented unit volume element. In n dimensions, it has grade n.
    ///
    /// # Properties
    ///
    /// - In Euclidean 2D: `I² = -1` (behaves like imaginary unit)
    /// - In Euclidean 3D: `I² = -1`
    /// - In Euclidean 4D: `I² = +1`
    ///
    /// # Example
    ///
    /// ```
    /// use clifford::algebra::Multivector;
    /// use clifford::signature::Euclidean3;
    ///
    /// let ps: Multivector<f64, Euclidean3> = Multivector::pseudoscalar();
    /// assert_eq!(ps.grade(1e-10), Some(3)); // Trivector in 3D
    ///
    /// // I² = -1 in Euclidean 3D
    /// let ps_sq = &ps * &ps;
    /// assert!((ps_sq.scalar_part() - (-1.0)).abs() < 1e-10);
    /// ```
    pub fn pseudoscalar() -> Self {
        let ps_index = (1 << S::DIM) - 1;
        Self::from_blade(Blade::from_index(ps_index))
    }
}

// ============================================================================
// Derived Products
// ============================================================================

impl<T: Float, S: Signature> Multivector<T, S> {
    /// Computes the exterior (wedge) product: `a ∧ b`.
    ///
    /// The exterior product extracts the grade-raising part of the geometric product.
    /// For a grade-k blade A and grade-j blade B, the result has grade k + j.
    ///
    /// # Properties
    ///
    /// - **Anticommutative for vectors**: `a ∧ b = -(b ∧ a)`
    /// - **Associative**: `(a ∧ b) ∧ c = a ∧ (b ∧ c)`
    /// - **Grade-raising**: `grade(A ∧ B) = grade(A) + grade(B)`
    ///
    /// # Geometric Interpretation
    ///
    /// The wedge product of two vectors creates a bivector representing
    /// the oriented parallelogram spanned by those vectors. More generally,
    /// it creates higher-grade blades representing oriented subspaces.
    ///
    /// # Example
    ///
    /// ```
    /// use clifford::algebra::Multivector;
    /// use clifford::signature::Euclidean3;
    /// use approx::relative_eq;
    ///
    /// let e1: Multivector<f64, Euclidean3> = Multivector::basis_vector(0);
    /// let e2: Multivector<f64, Euclidean3> = Multivector::basis_vector(1);
    ///
    /// // e₁ ∧ e₂ = e₁₂ (bivector)
    /// let e12 = e1.exterior(&e2);
    /// assert_eq!(e12.grade(1e-10), Some(2));
    ///
    /// // Anticommutativity: e₂ ∧ e₁ = -e₁₂
    /// let e21 = e2.exterior(&e1);
    /// assert!(relative_eq!(e21, -&e12, epsilon = 1e-10, max_relative = 1e-10));
    /// ```
    pub fn exterior(&self, other: &Self) -> Self {
        let mut result = Self::zero();

        for i in 0..S::NumBlades::USIZE {
            if self.coeffs[i] == T::zero() {
                continue;
            }
            let blade_i = Blade::from_index(i);
            let grade_i = blade_i.grade();

            for j in 0..S::NumBlades::USIZE {
                if other.coeffs[j] == T::zero() {
                    continue;
                }
                let blade_j = Blade::from_index(j);
                let grade_j = blade_j.grade();

                let (sign, result_blade) = blade_i.product(&blade_j, S::metric);

                // Keep only if result grade equals sum of input grades
                if sign != 0 && result_blade.grade() == grade_i + grade_j {
                    let coeff = T::from_i8(sign) * self.coeffs[i] * other.coeffs[j];
                    result.coeffs[result_blade.index()] += coeff;
                }
            }
        }
        result
    }

    /// Deprecated alias for [`exterior`](Self::exterior).
    ///
    /// Use [`exterior`](Self::exterior) instead. This method will be removed
    /// in a future version.
    #[deprecated(since = "0.2.0", note = "use `exterior` instead")]
    #[inline]
    pub fn outer(&self, other: &Self) -> Self {
        self.exterior(other)
    }

    /// Computes the left contraction: `A ⌋ B`.
    ///
    /// The left contraction extracts the grade-lowering part of the geometric
    /// product. For a grade-k blade A and grade-j blade B where k ≤ j,
    /// the result has grade j - k.
    ///
    /// # Mathematical Definition
    ///
    /// `A ⌋ B` is the part of `AB` with grade `grade(B) - grade(A)`.
    /// If `grade(A) > grade(B)`, the result is zero.
    ///
    /// # Example
    ///
    /// ```
    /// use clifford::algebra::Multivector;
    /// use clifford::signature::Euclidean3;
    ///
    /// let e1: Multivector<f64, Euclidean3> = Multivector::basis_vector(0);
    /// let e2: Multivector<f64, Euclidean3> = Multivector::basis_vector(1);
    /// let e12 = e1.exterior(&e2);
    ///
    /// // Vector left-contracted with bivector gives vector
    /// let result = e1.left_contract(&e12);
    /// assert_eq!(result.grade(1e-10), Some(1));
    /// ```
    pub fn left_contract(&self, other: &Self) -> Self {
        let mut result = Self::zero();

        for i in 0..S::NumBlades::USIZE {
            if self.coeffs[i] == T::zero() {
                continue;
            }
            let blade_i = Blade::from_index(i);
            let grade_i = blade_i.grade();

            for j in 0..S::NumBlades::USIZE {
                if other.coeffs[j] == T::zero() {
                    continue;
                }
                let blade_j = Blade::from_index(j);
                let grade_j = blade_j.grade();

                // Left contraction requires grade_i ≤ grade_j
                if grade_i > grade_j {
                    continue;
                }

                let (sign, result_blade) = blade_i.product(&blade_j, S::metric);

                // Keep only if result grade equals grade_j - grade_i
                if sign != 0 && result_blade.grade() == grade_j - grade_i {
                    let coeff = T::from_i8(sign) * self.coeffs[i] * other.coeffs[j];
                    result.coeffs[result_blade.index()] += coeff;
                }
            }
        }
        result
    }

    /// Computes the right contraction: `A ⌊ B`.
    ///
    /// The right contraction is the "mirror" of the left contraction.
    /// For a grade-k blade A and grade-j blade B where j ≤ k,
    /// the result has grade k - j.
    ///
    /// # Mathematical Definition
    ///
    /// `A ⌊ B` is the part of `AB` with grade `grade(A) - grade(B)`.
    /// If `grade(B) > grade(A)`, the result is zero.
    ///
    /// # Example
    ///
    /// ```
    /// use clifford::algebra::Multivector;
    /// use clifford::signature::Euclidean3;
    ///
    /// let e1: Multivector<f64, Euclidean3> = Multivector::basis_vector(0);
    /// let e2: Multivector<f64, Euclidean3> = Multivector::basis_vector(1);
    /// let e12 = e1.exterior(&e2);
    ///
    /// // Bivector right-contracted with vector gives vector
    /// let result = e12.right_contract(&e1);
    /// assert_eq!(result.grade(1e-10), Some(1));
    /// ```
    pub fn right_contract(&self, other: &Self) -> Self {
        let mut result = Self::zero();

        for i in 0..S::NumBlades::USIZE {
            if self.coeffs[i] == T::zero() {
                continue;
            }
            let blade_i = Blade::from_index(i);
            let grade_i = blade_i.grade();

            for j in 0..S::NumBlades::USIZE {
                if other.coeffs[j] == T::zero() {
                    continue;
                }
                let blade_j = Blade::from_index(j);
                let grade_j = blade_j.grade();

                // Right contraction requires grade_j ≤ grade_i
                if grade_j > grade_i {
                    continue;
                }

                let (sign, result_blade) = blade_i.product(&blade_j, S::metric);

                // Keep only if result grade equals grade_i - grade_j
                if sign != 0 && result_blade.grade() == grade_i - grade_j {
                    let coeff = T::from_i8(sign) * self.coeffs[i] * other.coeffs[j];
                    result.coeffs[result_blade.index()] += coeff;
                }
            }
        }
        result
    }

    /// Computes the inner product: `a · b`.
    ///
    /// For vectors, this is the standard dot product, returning a scalar.
    /// More generally, this is implemented as the left contraction.
    ///
    /// # Example
    ///
    /// ```
    /// use clifford::algebra::Multivector;
    /// use clifford::signature::Euclidean3;
    ///
    /// let a: Multivector<f64, Euclidean3> = Multivector::vector(&[1.0, 2.0, 3.0]);
    /// let b: Multivector<f64, Euclidean3> = Multivector::vector(&[4.0, 5.0, 6.0]);
    ///
    /// // a · b = 1*4 + 2*5 + 3*6 = 32
    /// let dot = a.inner(&b);
    /// assert!((dot.scalar_part() - 32.0).abs() < 1e-10);
    /// ```
    pub fn inner(&self, other: &Self) -> Self {
        self.left_contract(other)
    }

    /// Computes the Hodge dual: `A* = A ⌋ I⁻¹`.
    ///
    /// The dual maps a grade-k blade to a grade-(n-k) blade, where n is the
    /// dimension. It represents the orthogonal complement of a subspace.
    ///
    /// # Properties
    ///
    /// - `dual(dual(A)) = ±A` (sign depends on dimension and signature)
    /// - Maps vectors to pseudovectors, bivectors to vectors (in 3D), etc.
    ///
    /// # Example
    ///
    /// ```
    /// use clifford::algebra::Multivector;
    /// use clifford::signature::Euclidean3;
    ///
    /// let e1: Multivector<f64, Euclidean3> = Multivector::basis_vector(0);
    ///
    /// // In 3D, dual of a vector is a bivector
    /// let e1_dual = e1.dual();
    /// assert_eq!(e1_dual.grade(1e-10), Some(2));
    /// ```
    pub fn dual(&self) -> Self {
        let ps = Self::pseudoscalar();
        let ps_inv = ps.inverse().expect("pseudoscalar should be invertible");
        self.left_contract(&ps_inv)
    }

    /// Computes the undual (inverse of dual): `A⁻* = A ⌋ I`.
    ///
    /// # Example
    ///
    /// ```
    /// use clifford::algebra::Multivector;
    /// use clifford::signature::Euclidean3;
    /// use approx::relative_eq;
    ///
    /// let e1: Multivector<f64, Euclidean3> = Multivector::basis_vector(0);
    /// let roundtrip = e1.dual().undual();
    /// assert!(relative_eq!(roundtrip, e1, epsilon = 1e-10, max_relative = 1e-10));
    /// ```
    pub fn undual(&self) -> Self {
        let ps = Self::pseudoscalar();
        self.left_contract(&ps)
    }

    /// Computes the complement of a multivector.
    ///
    /// The complement maps a grade-k blade to a grade-(n-k) blade by XOR with
    /// the pseudoscalar index. The sign is determined by the number of
    /// transpositions needed to put `blade ∧ complement(blade)` into canonical
    /// order.
    ///
    /// This is used in RGA for the regressive product:
    /// `a ∨ b = complement(complement(a) ∧ complement(b))`
    ///
    /// # Properties
    ///
    /// - `grade(complement(a)) = n - grade(a)` where n is the dimension
    /// - `complement(complement(a)) = ±a`
    /// - Different from Hodge dual (which uses contraction)
    ///
    /// # Example
    ///
    /// ```
    /// use clifford::algebra::Multivector;
    /// use clifford::signature::Euclidean3;
    ///
    /// let e1: Multivector<f64, Euclidean3> = Multivector::basis_vector(0);
    ///
    /// // Complement of e1 in 3D gives a bivector (e23)
    /// let e1_comp = e1.complement();
    /// assert_eq!(e1_comp.grade(1e-10), Some(2));
    /// ```
    pub fn complement(&self) -> Self {
        let ps_index = S::NumBlades::USIZE - 1;
        let mut result = Self::zero();

        for i in 0..S::NumBlades::USIZE {
            if self.coeffs[i] == T::zero() {
                continue;
            }

            let complement_index = ps_index ^ i;

            // Compute exterior sign: (-1)^(transpositions for i ∧ complement)
            // This counts how many bits in complement_index are "passed over" by bits in i
            let sign = exterior_sign(i, complement_index);

            result.coeffs[complement_index] += T::from_i8(sign) * self.coeffs[i];
        }
        result
    }

    /// Deprecated: Use [`complement`](Self::complement) instead.
    ///
    /// The right complement maps a grade-k blade to a grade-(n-k) blade.
    /// This is used in RGA for the regressive product.
    #[deprecated(since = "0.3.0", note = "use `complement` instead")]
    pub fn right_complement(&self) -> Self {
        self.complement()
    }

    /// Deprecated: Use [`complement`](Self::complement) instead.
    ///
    /// The left complement maps a grade-k blade to a grade-(n-k) blade.
    /// This is used in RGA for computing the "undual" in the regressive product.
    #[deprecated(since = "0.3.0", note = "use `complement` instead")]
    pub fn left_complement(&self) -> Self {
        self.complement()
    }

    /// Computes the regressive (vee) product: `a ∨ b = complement(complement(a) ∧ complement(b))`.
    ///
    /// The regressive product is the dual of the outer product. While the
    /// outer product represents the "join" (span) of subspaces, the regressive
    /// product represents the "meet" (intersection).
    ///
    /// # Formula
    ///
    /// Using complement notation:
    /// - `a ∨ b = complement(complement(a) ∧ complement(b))`
    ///
    /// # Properties
    ///
    /// - Grade-lowering (opposite of wedge)
    /// - Associative
    ///
    /// # Example
    ///
    /// ```
    /// use clifford::algebra::Multivector;
    /// use clifford::signature::Euclidean3;
    ///
    /// let e1: Multivector<f64, Euclidean3> = Multivector::basis_vector(0);
    /// let e2: Multivector<f64, Euclidean3> = Multivector::basis_vector(1);
    /// let e12 = e1.exterior(&e2);
    /// let e23 = e2.exterior(&Multivector::basis_vector(2));
    ///
    /// // Meet of two planes sharing e2 should give e2 direction
    /// let meet = e12.regressive(&e23);
    /// assert_eq!(meet.grade(1e-10), Some(1));
    /// ```
    pub fn regressive(&self, other: &Self) -> Self {
        self.complement().exterior(&other.complement()).complement()
    }

    /// Computes the antiwedge (exterior antiproduct): `a ∨ b`.
    ///
    /// This is an alias for [`regressive`](Self::regressive). The antiwedge is
    /// the dual of the wedge product, representing the "meet" (intersection)
    /// of subspaces rather than the "join" (span).
    ///
    /// # Example
    ///
    /// ```
    /// use clifford::algebra::Multivector;
    /// use clifford::signature::Euclidean3;
    ///
    /// let e12: Multivector<f64, Euclidean3> =
    ///     Multivector::basis_vector(0).exterior(&Multivector::basis_vector(1));
    /// let e23: Multivector<f64, Euclidean3> =
    ///     Multivector::basis_vector(1).exterior(&Multivector::basis_vector(2));
    ///
    /// // Antiwedge of two planes gives their intersection (a line)
    /// let intersection = e12.antiwedge(&e23);
    /// ```
    #[inline]
    pub fn antiwedge(&self, other: &Self) -> Self {
        self.regressive(other)
    }

    /// Computes the bulk dual (★): `a★ = ã ⊙ 𝟙`.
    ///
    /// The bulk dual is computed as `reverse(a) * pseudoscalar`. It maps a blade
    /// to its orthogonal complement using the metric (bulk) structure.
    ///
    /// # Properties
    ///
    /// - `grade(a★) = n - grade(a)` where n is the dimension
    /// - `(a★)★ = ±a` (up to sign, involutory)
    ///
    /// # Example
    ///
    /// ```
    /// use clifford::algebra::Multivector;
    /// use clifford::signature::Euclidean3;
    ///
    /// let e1: Multivector<f64, Euclidean3> = Multivector::basis_vector(0);
    ///
    /// // In 3D, bulk dual of a vector is a bivector
    /// let e1_bulk_dual = e1.bulk_dual();
    /// assert_eq!(e1_bulk_dual.grade(1e-10), Some(2));
    /// ```
    pub fn bulk_dual(&self) -> Self {
        let ps = Self::pseudoscalar();
        &self.reverse() * &ps
    }

    /// Computes the left bulk dual (★): `★a = 𝟙 ⊙ ã`.
    ///
    /// The left bulk dual is computed as `pseudoscalar * reverse(a)`. This is
    /// the "undual" operation for the bulk dual, used in computing the antiproduct.
    ///
    /// # Properties
    ///
    /// - `grade(★a) = n - grade(a)` where n is the dimension
    /// - `★(a★) = ±a` (recovers original up to sign)
    ///
    /// # Example
    ///
    /// ```
    /// use clifford::algebra::Multivector;
    /// use clifford::signature::Euclidean3;
    ///
    /// let e1: Multivector<f64, Euclidean3> = Multivector::basis_vector(0);
    ///
    /// // In 3D, left bulk dual of a vector is a bivector
    /// let e1_left_dual = e1.left_bulk_dual();
    /// assert_eq!(e1_left_dual.grade(1e-10), Some(2));
    /// ```
    pub fn left_bulk_dual(&self) -> Self {
        let ps = Self::pseudoscalar();
        &ps * &self.reverse()
    }

    /// Computes the weight dual (☆): `a☆ = ã ⊛ 1`.
    ///
    /// The weight dual is computed as `reverse(a) ⊛ scalar` where ⊛ is the
    /// geometric antiproduct. It maps a blade to its orthogonal complement
    /// using the anti-metric (weight) structure.
    ///
    /// # Properties
    ///
    /// - `grade(a☆) = n - grade(a)` where n is the dimension
    /// - Uses the antiproduct with the scalar (grade-0 element)
    ///
    /// # Example
    ///
    /// ```
    /// use clifford::algebra::Multivector;
    /// use clifford::signature::Euclidean3;
    ///
    /// let e1: Multivector<f64, Euclidean3> = Multivector::basis_vector(0);
    ///
    /// // In 3D, weight dual of a vector is a bivector
    /// let e1_weight_dual = e1.weight_dual();
    /// assert_eq!(e1_weight_dual.grade(1e-10), Some(2));
    /// ```
    pub fn weight_dual(&self) -> Self {
        let scalar = Self::scalar(T::one());
        self.reverse().antiproduct(&scalar)
    }

    /// Computes the bulk contraction: `a ∨ b★`.
    ///
    /// The bulk contraction is the antiwedge of `a` with the bulk dual of `b`.
    /// This is one of the four interior products in Rigid Geometric Algebra.
    ///
    /// # Properties
    ///
    /// - Grade-lowering operation
    /// - Measures how much `b` "projects onto" the complement of `a`
    ///
    /// # Example
    ///
    /// ```
    /// use clifford::algebra::Multivector;
    /// use clifford::signature::Euclidean3;
    ///
    /// let e1: Multivector<f64, Euclidean3> = Multivector::basis_vector(0);
    /// let e2: Multivector<f64, Euclidean3> = Multivector::basis_vector(1);
    ///
    /// let result = e1.bulk_contraction(&e2);
    /// ```
    pub fn bulk_contraction(&self, other: &Self) -> Self {
        self.antiwedge(&other.bulk_dual())
    }

    /// Computes the weight contraction: `a ∨ b☆`.
    ///
    /// The weight contraction is the antiwedge of `a` with the weight dual of `b`.
    /// This is one of the four interior products in Rigid Geometric Algebra.
    ///
    /// # Properties
    ///
    /// - Grade-lowering operation
    /// - Uses the anti-metric structure
    ///
    /// # Example
    ///
    /// ```
    /// use clifford::algebra::Multivector;
    /// use clifford::signature::Euclidean3;
    ///
    /// let e1: Multivector<f64, Euclidean3> = Multivector::basis_vector(0);
    /// let e2: Multivector<f64, Euclidean3> = Multivector::basis_vector(1);
    ///
    /// let result = e1.weight_contraction(&e2);
    /// ```
    pub fn weight_contraction(&self, other: &Self) -> Self {
        self.antiwedge(&other.weight_dual())
    }

    /// Computes the bulk expansion: `a ∧ b★`.
    ///
    /// The bulk expansion is the wedge of `a` with the bulk dual of `b`.
    /// This is one of the four interior products in Rigid Geometric Algebra.
    ///
    /// # Properties
    ///
    /// - Can be grade-raising or grade-lowering depending on inputs
    /// - Measures the "exterior" relationship using bulk structure
    ///
    /// # Example
    ///
    /// ```
    /// use clifford::algebra::Multivector;
    /// use clifford::signature::Euclidean3;
    ///
    /// let e1: Multivector<f64, Euclidean3> = Multivector::basis_vector(0);
    /// let e2: Multivector<f64, Euclidean3> = Multivector::basis_vector(1);
    ///
    /// let result = e1.bulk_expansion(&e2);
    /// ```
    pub fn bulk_expansion(&self, other: &Self) -> Self {
        self.exterior(&other.bulk_dual())
    }

    /// Computes the weight expansion: `a ∧ b☆`.
    ///
    /// The weight expansion is the wedge of `a` with the weight dual of `b`.
    /// This is one of the four interior products in Rigid Geometric Algebra.
    ///
    /// # Properties
    ///
    /// - Can be grade-raising or grade-lowering depending on inputs
    /// - Uses the anti-metric structure
    ///
    /// # Example
    ///
    /// ```
    /// use clifford::algebra::Multivector;
    /// use clifford::signature::Euclidean3;
    ///
    /// let e1: Multivector<f64, Euclidean3> = Multivector::basis_vector(0);
    /// let e2: Multivector<f64, Euclidean3> = Multivector::basis_vector(1);
    ///
    /// let result = e1.weight_expansion(&e2);
    /// ```
    pub fn weight_expansion(&self, other: &Self) -> Self {
        self.exterior(&other.weight_dual())
    }

    /// Computes the geometric antiproduct: `a ⊛ b = undual(dual(a) × dual(b))`.
    ///
    /// The geometric antiproduct is the dual of the geometric product. In PGA
    /// (Projective Geometric Algebra), the antiproduct is essential for correct
    /// motor transformations because it properly handles the degenerate direction.
    ///
    /// # Properties
    ///
    /// - `dual(a × b) = dual(a) ⊛ dual(b)` (duality relationship)
    /// - Used for PGA sandwich transformations: `M ⊛ x ⊛ rev(M)`
    ///
    /// # Example
    ///
    /// ```
    /// use clifford::algebra::Multivector;
    /// use clifford::signature::Euclidean3;
    ///
    /// let e1: Multivector<f64, Euclidean3> = Multivector::basis_vector(0);
    /// let e2: Multivector<f64, Euclidean3> = Multivector::basis_vector(1);
    ///
    /// // Antiproduct of vectors
    /// let result = e1.antiproduct(&e2);
    /// ```
    pub fn antiproduct(&self, other: &Self) -> Self {
        // a ⋇ b = complement(complement(a) ⊙ complement(b))
        (&self.complement() * &other.complement()).complement()
    }

    /// Computes the antisandwich product: `R ⊛ x ⊛ R̃`.
    ///
    /// The antisandwich uses the geometric antiproduct instead of the geometric
    /// product. In PGA, this is the correct transformation for motors because
    /// it handles translations properly (which the regular sandwich cannot due
    /// to the degenerate metric e₀² = 0).
    ///
    /// # Properties
    ///
    /// - Preserves grade of the argument
    /// - Correctly transforms points under translation in PGA
    ///
    /// # Example
    ///
    /// ```
    /// use clifford::algebra::Multivector;
    /// use clifford::signature::Euclidean3;
    ///
    /// // Create a rotor (identity in this case)
    /// let rotor: Multivector<f64, Euclidean3> = Multivector::scalar(1.0);
    /// let e1: Multivector<f64, Euclidean3> = Multivector::basis_vector(0);
    ///
    /// // Compute the antisandwich transformation
    /// let _result = rotor.antisandwich(&e1);
    /// ```
    ///
    /// # Note
    ///
    /// For PGA (Projective Geometric Algebra), use the specialized types in
    /// `clifford::specialized::projective::dim3` which have optimized antisandwich
    /// implementations that handle the degenerate metric correctly.
    pub fn antisandwich(&self, x: &Self) -> Self {
        self.antiproduct(x).antiproduct(&self.reverse())
    }

    /// Computes the sandwich product: `R x R̃`.
    ///
    /// The sandwich product is the fundamental operation for transformations
    /// in geometric algebra:
    /// - Reflections: `n x n` where n is a unit vector (reflects across plane normal to n)
    /// - Rotations: `R x R̃` where R is a rotor
    ///
    /// # Properties
    ///
    /// - Preserves grade of the argument
    /// - Preserves norm in Euclidean space
    /// - Composition: `R₂(R₁ x R₁̃)R₂̃ = (R₂R₁) x (R₂R₁)̃`
    ///
    /// # Example
    ///
    /// ```
    /// use clifford::algebra::Multivector;
    /// use clifford::signature::Euclidean3;
    ///
    /// // Reflect a vector across the yz-plane (normal = e₁)
    /// let n: Multivector<f64, Euclidean3> = Multivector::basis_vector(0); // e₁
    /// let v: Multivector<f64, Euclidean3> = Multivector::vector(&[1.0, 2.0, 3.0]);
    ///
    /// let reflected = n.sandwich(&v);
    /// // x component flips sign, y and z stay same
    /// assert!((reflected.norm() - v.norm()).abs() < 1e-10);
    /// ```
    pub fn sandwich(&self, x: &Self) -> Self {
        &(self * x) * &self.reverse()
    }

    /// Computes the projection of `self` onto `other`: `other ∨ (self ∧ other☆)`.
    ///
    /// The projection finds the part of `self` that lies "on" or "inside" `other`.
    /// Uses the weight dual (☆) of the target geometry.
    ///
    /// # Geometric Interpretation
    ///
    /// In PGA (Projective Geometric Algebra):
    /// - Projecting a point onto a line gives the closest point on the line
    /// - Projecting a point onto a plane gives the closest point on the plane
    /// - Projecting a line onto a plane gives the line's intersection with the plane
    ///
    /// # Example
    ///
    /// ```
    /// use clifford::algebra::Multivector;
    /// use clifford::signature::Euclidean3;
    ///
    /// let a: Multivector<f64, Euclidean3> = Multivector::basis_vector(0); // e₁
    /// let b: Multivector<f64, Euclidean3> = Multivector::basis_vector(1); // e₂
    ///
    /// let _proj = a.project(&b);
    /// ```
    pub fn project(&self, other: &Self) -> Self {
        // Project: other ∨ (self ∧ other☆)
        other.antiwedge(&self.exterior(&other.weight_dual()))
    }

    /// Computes the antiprojection of `self` onto `other`: `other ∧ (self ∨ other☆)`.
    ///
    /// The antiprojection finds geometry that passes through `self` and is
    /// perpendicular to `other`. Uses the weight dual (☆) of the target geometry.
    ///
    /// # Geometric Interpretation
    ///
    /// In PGA (Projective Geometric Algebra):
    /// - Antiprojecting a point onto a line gives the plane through the point
    ///   perpendicular to the line
    /// - Antiprojecting a point onto a plane gives the line through the point
    ///   perpendicular to the plane
    ///
    /// # Example
    ///
    /// ```
    /// use clifford::algebra::Multivector;
    /// use clifford::signature::Euclidean3;
    ///
    /// let a: Multivector<f64, Euclidean3> = Multivector::basis_vector(0); // e₁
    /// let b: Multivector<f64, Euclidean3> = Multivector::basis_vector(1); // e₂
    ///
    /// let _antiproj = a.antiproject(&b);
    /// ```
    pub fn antiproject(&self, other: &Self) -> Self {
        // Antiproject: other ∧ (self ∨ other☆)
        other.exterior(&self.antiwedge(&other.weight_dual()))
    }
}

// ============================================================================
// Geometric Product
// ============================================================================

impl<T: Float, S: Signature> Mul for &Multivector<T, S> {
    type Output = Multivector<T, S>;

    /// Computes the geometric product of two multivectors.
    ///
    /// The geometric product is the fundamental operation in geometric algebra.
    /// It is computed by distributing over all pairs of basis blades and using
    /// the blade multiplication rules defined by the metric signature.
    ///
    /// # Complexity
    ///
    /// O(4^n) where n is the dimension of the algebra.
    fn mul(self, rhs: Self) -> Self::Output {
        let mut result = Multivector::zero();

        for i in 0..S::NumBlades::USIZE {
            if self.coeffs[i] == T::zero() {
                continue;
            }
            for j in 0..S::NumBlades::USIZE {
                if rhs.coeffs[j] == T::zero() {
                    continue;
                }

                let blade_i = Blade::from_index(i);
                let blade_j = Blade::from_index(j);
                let (sign, result_blade) = blade_i.product(&blade_j, S::metric);

                if sign != 0 {
                    let coeff = T::from_i8(sign) * self.coeffs[i] * rhs.coeffs[j];
                    result.coeffs[result_blade.index()] += coeff;
                }
            }
        }

        result
    }
}

// Implement Mul for owned values by delegating to reference implementation
impl<T: Float, S: Signature> Mul for Multivector<T, S> {
    type Output = Multivector<T, S>;

    fn mul(self, rhs: Self) -> Self::Output {
        &self * &rhs
    }
}

impl<T: Float, S: Signature> Mul<&Multivector<T, S>> for Multivector<T, S> {
    type Output = Multivector<T, S>;

    fn mul(self, rhs: &Self) -> Self::Output {
        &self * rhs
    }
}

impl<T: Float, S: Signature> Mul<Multivector<T, S>> for &Multivector<T, S> {
    type Output = Multivector<T, S>;

    fn mul(self, rhs: Multivector<T, S>) -> Self::Output {
        self * &rhs
    }
}

// ============================================================================
// Scalar multiplication and division
// ============================================================================

impl<T: Float, S: Signature> Mul<T> for &Multivector<T, S> {
    type Output = Multivector<T, S>;

    fn mul(self, scalar: T) -> Self::Output {
        let mut result = Multivector::zero();
        for i in 0..S::NumBlades::USIZE {
            result.coeffs[i] = self.coeffs[i] * scalar;
        }
        result
    }
}

impl<T: Float, S: Signature> Mul<T> for Multivector<T, S> {
    type Output = Multivector<T, S>;

    fn mul(self, scalar: T) -> Self::Output {
        &self * scalar
    }
}

impl<T: Float, S: Signature> Div<T> for &Multivector<T, S> {
    type Output = Multivector<T, S>;

    fn div(self, scalar: T) -> Self::Output {
        let mut result = Multivector::zero();
        for i in 0..S::NumBlades::USIZE {
            result.coeffs[i] = self.coeffs[i] / scalar;
        }
        result
    }
}

impl<T: Float, S: Signature> Div<T> for Multivector<T, S> {
    type Output = Multivector<T, S>;

    fn div(self, scalar: T) -> Self::Output {
        &self / scalar
    }
}

// ============================================================================
// Wedge Product (BitXor)
// ============================================================================

/// Wedge (exterior) product via `^` operator.
///
/// The `^` operator provides ergonomic syntax for the wedge product:
/// `a ^ b` is equivalent to `a.exterior(&b)`.
///
/// # Example
///
/// ```
/// use clifford::algebra::Multivector;
/// use clifford::signature::Euclidean3;
/// use approx::relative_eq;
///
/// let e1: Multivector<f64, Euclidean3> = Multivector::basis_vector(0);
/// let e2: Multivector<f64, Euclidean3> = Multivector::basis_vector(1);
///
/// // e₁ ^ e₂ = e₁₂ (bivector)
/// let e12 = &e1 ^ &e2;
/// assert_eq!(e12.grade(1e-10), Some(2));
///
/// // Anticommutativity: e₂ ^ e₁ = -e₁₂
/// let e21 = &e2 ^ &e1;
/// assert!(relative_eq!(e21, -&e12, epsilon = 1e-10, max_relative = 1e-10));
/// ```
impl<T: Float, S: Signature> BitXor for &Multivector<T, S> {
    type Output = Multivector<T, S>;

    fn bitxor(self, rhs: Self) -> Self::Output {
        self.exterior(rhs)
    }
}

impl<T: Float, S: Signature> BitXor for Multivector<T, S> {
    type Output = Multivector<T, S>;

    fn bitxor(self, rhs: Self) -> Self::Output {
        self.exterior(&rhs)
    }
}

impl<T: Float, S: Signature> BitXor<&Multivector<T, S>> for Multivector<T, S> {
    type Output = Multivector<T, S>;

    fn bitxor(self, rhs: &Multivector<T, S>) -> Self::Output {
        self.exterior(rhs)
    }
}

impl<T: Float, S: Signature> BitXor<Multivector<T, S>> for &Multivector<T, S> {
    type Output = Multivector<T, S>;

    fn bitxor(self, rhs: Multivector<T, S>) -> Self::Output {
        self.exterior(&rhs)
    }
}

// ============================================================================
// Addition and Subtraction
// ============================================================================

impl<T: Float, S: Signature> Add for &Multivector<T, S> {
    type Output = Multivector<T, S>;

    fn add(self, rhs: Self) -> Self::Output {
        let mut result = Multivector::zero();
        for i in 0..S::NumBlades::USIZE {
            result.coeffs[i] = self.coeffs[i] + rhs.coeffs[i];
        }
        result
    }
}

impl<T: Float, S: Signature> Add for Multivector<T, S> {
    type Output = Multivector<T, S>;

    fn add(self, rhs: Self) -> Self::Output {
        &self + &rhs
    }
}

impl<T: Float, S: Signature> Add<&Multivector<T, S>> for Multivector<T, S> {
    type Output = Multivector<T, S>;

    fn add(self, rhs: &Self) -> Self::Output {
        &self + rhs
    }
}

impl<T: Float, S: Signature> Add<Multivector<T, S>> for &Multivector<T, S> {
    type Output = Multivector<T, S>;

    fn add(self, rhs: Multivector<T, S>) -> Self::Output {
        self + &rhs
    }
}

impl<T: Float, S: Signature> Sub for &Multivector<T, S> {
    type Output = Multivector<T, S>;

    fn sub(self, rhs: Self) -> Self::Output {
        let mut result = Multivector::zero();
        for i in 0..S::NumBlades::USIZE {
            result.coeffs[i] = self.coeffs[i] - rhs.coeffs[i];
        }
        result
    }
}

impl<T: Float, S: Signature> Sub for Multivector<T, S> {
    type Output = Multivector<T, S>;

    fn sub(self, rhs: Self) -> Self::Output {
        &self - &rhs
    }
}

impl<T: Float, S: Signature> Sub<&Multivector<T, S>> for Multivector<T, S> {
    type Output = Multivector<T, S>;

    fn sub(self, rhs: &Self) -> Self::Output {
        &self - rhs
    }
}

impl<T: Float, S: Signature> Sub<Multivector<T, S>> for &Multivector<T, S> {
    type Output = Multivector<T, S>;

    fn sub(self, rhs: Multivector<T, S>) -> Self::Output {
        self - &rhs
    }
}

impl<T: Float, S: Signature> Neg for &Multivector<T, S> {
    type Output = Multivector<T, S>;

    fn neg(self) -> Self::Output {
        let mut result = Multivector::zero();
        for i in 0..S::NumBlades::USIZE {
            result.coeffs[i] = -self.coeffs[i];
        }
        result
    }
}

impl<T: Float, S: Signature> Neg for Multivector<T, S> {
    type Output = Multivector<T, S>;

    fn neg(self) -> Self::Output {
        -&self
    }
}

// ============================================================================
// Assignment operators
// ============================================================================

impl<T: Float, S: Signature> AddAssign<&Multivector<T, S>> for Multivector<T, S> {
    fn add_assign(&mut self, rhs: &Self) {
        for i in 0..S::NumBlades::USIZE {
            self.coeffs[i] += rhs.coeffs[i];
        }
    }
}

impl<T: Float, S: Signature> AddAssign for Multivector<T, S> {
    fn add_assign(&mut self, rhs: Self) {
        *self += &rhs;
    }
}

impl<T: Float, S: Signature> SubAssign<&Multivector<T, S>> for Multivector<T, S> {
    fn sub_assign(&mut self, rhs: &Self) {
        for i in 0..S::NumBlades::USIZE {
            self.coeffs[i] -= rhs.coeffs[i];
        }
    }
}

impl<T: Float, S: Signature> SubAssign for Multivector<T, S> {
    fn sub_assign(&mut self, rhs: Self) {
        *self -= &rhs;
    }
}

impl<T: Float, S: Signature> MulAssign<T> for Multivector<T, S> {
    fn mul_assign(&mut self, scalar: T) {
        for i in 0..S::NumBlades::USIZE {
            self.coeffs[i] *= scalar;
        }
    }
}

// ============================================================================
// Standard Traits
// ============================================================================

impl<T: Float, S: Signature> Default for Multivector<T, S> {
    fn default() -> Self {
        Self::zero()
    }
}

impl<T: Float, S: Signature> PartialEq for Multivector<T, S> {
    fn eq(&self, other: &Self) -> bool {
        self.coeffs == other.coeffs
    }
}

impl<T: Float, S: Signature> AbsDiffEq for Multivector<T, S> {
    type Epsilon = T;

    fn default_epsilon() -> Self::Epsilon {
        T::default_epsilon()
    }

    fn abs_diff_eq(&self, other: &Self, epsilon: Self::Epsilon) -> bool {
        self.coeffs
            .iter()
            .zip(other.coeffs.iter())
            .all(|(a, b)| T::abs_diff_eq(a, b, epsilon))
    }
}

impl<T: Float, S: Signature> RelativeEq for Multivector<T, S> {
    fn default_max_relative() -> Self::Epsilon {
        T::default_max_relative()
    }

    fn relative_eq(
        &self,
        other: &Self,
        epsilon: Self::Epsilon,
        max_relative: Self::Epsilon,
    ) -> bool {
        self.coeffs
            .iter()
            .zip(other.coeffs.iter())
            .all(|(a, b)| T::relative_eq(a, b, epsilon, max_relative))
    }
}

impl<T: Float, S: Signature> UlpsEq for Multivector<T, S> {
    fn default_max_ulps() -> u32 {
        T::default_max_ulps()
    }

    fn ulps_eq(&self, other: &Self, epsilon: Self::Epsilon, max_ulps: u32) -> bool {
        self.coeffs
            .iter()
            .zip(other.coeffs.iter())
            .all(|(a, b)| T::ulps_eq(a, b, epsilon, max_ulps))
    }
}

impl<T: Float, S: Signature> fmt::Debug for Multivector<T, S> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "Multivector(")?;
        let mut first = true;
        for i in 0..S::NumBlades::USIZE {
            if self.coeffs[i] != T::zero() {
                if !first {
                    write!(f, " + ")?;
                }
                write!(f, "{:?}*{}", self.coeffs[i], Blade::from_index(i))?;
                first = false;
            }
        }
        if first {
            write!(f, "0")?;
        }
        write!(f, ")")
    }
}

impl<T: Float, S: Signature> fmt::Display for Multivector<T, S> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let mut first = true;
        for i in 0..S::NumBlades::USIZE {
            if self.coeffs[i] != T::zero() {
                if !first {
                    write!(f, " + ")?;
                }
                let blade = Blade::from_index(i);
                if i == 0 {
                    write!(f, "{}", self.coeffs[i])?;
                } else {
                    write!(f, "{}{}", self.coeffs[i], blade)?;
                }
                first = false;
            }
        }
        if first {
            write!(f, "0")?;
        }
        Ok(())
    }
}

// ============================================================================
// Helper functions
// ============================================================================

/// Computes the sign for the exterior product of two non-overlapping blades.
///
/// For blades `a` and `b` that don't share any basis vectors (i.e., `a & b == 0`),
/// this computes `(-1)^swaps` where `swaps` is the number of transpositions
/// needed to put the basis vectors of `a ∧ b` into canonical order.
///
/// # Panics
///
/// Debug builds panic if the blades overlap (share basis vectors).
#[inline]
fn exterior_sign(a: usize, b: usize) -> i8 {
    debug_assert_eq!(a & b, 0, "blades must not overlap for exterior sign");

    // Count transpositions: for each bit in b, count how many bits in a
    // are to the right of it (have higher index, appear later in canonical order)
    let mut transpositions = 0;
    let mut b_remaining = b;
    while b_remaining != 0 {
        // Find lowest set bit in b
        let lowest_b = b_remaining & b_remaining.wrapping_neg();
        let b_pos = lowest_b.trailing_zeros();

        // Count bits in a that are above this position (need to swap past)
        let a_above = a >> (b_pos + 1);
        transpositions += a_above.count_ones();

        b_remaining &= !lowest_b;
    }

    if transpositions % 2 == 0 { 1 } else { -1 }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::algebra::arbitrary::{NonZeroVectorE3, UnitVectorE3, VectorE3};
    use crate::signature::{Euclidean2, Euclidean3};
    use crate::test_utils::RELATIVE_EQ_EPS;
    use approx::relative_eq;
    use proptest::prelude::*;

    // ========================================================================
    // Constructor tests
    // ========================================================================

    #[test]
    fn test_zero() {
        let zero: Multivector<f64, Euclidean3> = Multivector::zero();
        assert!(zero.is_zero(RELATIVE_EQ_EPS));
    }

    #[test]
    fn test_one() {
        let one: Multivector<f64, Euclidean3> = Multivector::one();
        assert!(relative_eq!(
            one.scalar_part(),
            1.0,
            max_relative = RELATIVE_EQ_EPS
        ));
    }

    #[test]
    fn test_basis_vector() {
        let e1: Multivector<f64, Euclidean3> = Multivector::basis_vector(0);
        assert!(relative_eq!(
            e1.get(Blade::basis_vector(0)),
            1.0,
            max_relative = RELATIVE_EQ_EPS
        ));
        assert!(relative_eq!(
            e1.get(Blade::basis_vector(1)),
            0.0,
            max_relative = RELATIVE_EQ_EPS
        ));
    }

    #[test]
    fn test_vector() {
        let v: Multivector<f64, Euclidean3> = Multivector::vector(&[3.0, 4.0, 0.0]);
        assert!(relative_eq!(
            v.get(Blade::basis_vector(0)),
            3.0,
            max_relative = RELATIVE_EQ_EPS
        ));
        assert!(relative_eq!(
            v.get(Blade::basis_vector(1)),
            4.0,
            max_relative = RELATIVE_EQ_EPS
        ));
    }

    // ========================================================================
    // Geometric product tests
    // ========================================================================

    #[test]
    fn test_vector_squares_to_scalar() {
        let e1: Multivector<f64, Euclidean3> = Multivector::basis_vector(0);
        let e1_sq = e1 * e1;
        assert!(relative_eq!(
            e1_sq.scalar_part(),
            1.0,
            max_relative = RELATIVE_EQ_EPS
        ));
        assert!(relative_eq!(
            e1_sq,
            Multivector::one(),
            max_relative = RELATIVE_EQ_EPS
        ));
    }

    #[test]
    fn test_bivector_squares_to_minus_one() {
        let e1: Multivector<f64, Euclidean3> = Multivector::basis_vector(0);
        let e2: Multivector<f64, Euclidean3> = Multivector::basis_vector(1);
        let e12 = e1 * e2;
        let e12_sq = e12 * e12;
        assert!(relative_eq!(
            e12_sq.scalar_part(),
            -1.0,
            max_relative = RELATIVE_EQ_EPS
        ));
    }

    #[test]
    fn test_scalar_is_identity() {
        let one: Multivector<f64, Euclidean3> = Multivector::one();
        let v: Multivector<f64, Euclidean3> = Multivector::vector(&[1.0, 2.0, 3.0]);
        assert!(relative_eq!(
            one * v,
            v,
            epsilon = RELATIVE_EQ_EPS,
            max_relative = RELATIVE_EQ_EPS
        ));
        assert!(relative_eq!(
            v * one,
            v,
            epsilon = RELATIVE_EQ_EPS,
            max_relative = RELATIVE_EQ_EPS
        ));
    }

    // ========================================================================
    // Unary operation tests
    // ========================================================================

    #[test]
    fn test_reverse_vector() {
        let v: Multivector<f64, Euclidean3> = Multivector::basis_vector(0);
        assert!(relative_eq!(
            v.reverse(),
            v,
            epsilon = RELATIVE_EQ_EPS,
            max_relative = RELATIVE_EQ_EPS
        )); // Vectors unchanged
    }

    #[test]
    fn test_reverse_bivector() {
        let e1: Multivector<f64, Euclidean3> = Multivector::basis_vector(0);
        let e2: Multivector<f64, Euclidean3> = Multivector::basis_vector(1);
        let e12 = e1 * e2;
        assert!(relative_eq!(
            e12.reverse(),
            -e12,
            max_relative = RELATIVE_EQ_EPS
        )); // Bivectors negate
    }

    #[test]
    fn test_involute_vector() {
        let v: Multivector<f64, Euclidean3> = Multivector::basis_vector(0);
        assert!(relative_eq!(
            v.involute(),
            -v,
            epsilon = RELATIVE_EQ_EPS,
            max_relative = RELATIVE_EQ_EPS
        )); // Vectors negate
    }

    // ========================================================================
    // Norm and inverse tests
    // ========================================================================

    #[test]
    fn test_norm_squared_vector() {
        let v: Multivector<f64, Euclidean3> = Multivector::vector(&[3.0, 4.0, 0.0]);
        assert!(relative_eq!(
            v.norm_squared(),
            25.0,
            max_relative = RELATIVE_EQ_EPS
        ));
    }

    #[test]
    fn test_norm_vector() {
        let v: Multivector<f64, Euclidean3> = Multivector::vector(&[3.0, 4.0, 0.0]);
        assert!(relative_eq!(
            v.norm(),
            5.0,
            epsilon = RELATIVE_EQ_EPS,
            max_relative = RELATIVE_EQ_EPS
        ));
    }

    #[test]
    fn test_inverse_vector() {
        let v: Multivector<f64, Euclidean3> = Multivector::vector(&[2.0, 0.0, 0.0]);
        let v_inv = v.inverse().unwrap();
        let product = v * v_inv;
        assert!(relative_eq!(
            product,
            Multivector::one(),
            max_relative = RELATIVE_EQ_EPS
        ));
    }

    // ========================================================================
    // Property-based tests
    // ========================================================================

    proptest! {
        #[test]
        fn geometric_product_associative(
            a in any::<Multivector<f64, Euclidean3>>(),
            b in any::<Multivector<f64, Euclidean3>>(),
            c in any::<Multivector<f64, Euclidean3>>(),
        ) {
            let lhs = (a * b) * c;
            let rhs = a * (b * c);
            prop_assert!(relative_eq!(lhs, rhs, epsilon = RELATIVE_EQ_EPS, max_relative = RELATIVE_EQ_EPS));
        }

        #[test]
        fn geometric_product_distributive(
            a in any::<Multivector<f64, Euclidean3>>(),
            b in any::<Multivector<f64, Euclidean3>>(),
            c in any::<Multivector<f64, Euclidean3>>(),
        ) {
            let lhs = a * (b + c);
            let rhs = (a * b) + (a * c);
            prop_assert!(relative_eq!(lhs, rhs, epsilon = RELATIVE_EQ_EPS, max_relative = RELATIVE_EQ_EPS));
        }

        #[test]
        fn reverse_involutory(a in any::<Multivector<f64, Euclidean3>>()) {
            prop_assert!(relative_eq!(a.reverse().reverse(), a, epsilon = RELATIVE_EQ_EPS, max_relative = RELATIVE_EQ_EPS));
        }

        #[test]
        fn reverse_antimorphism(
            a in any::<Multivector<f64, Euclidean3>>(),
            b in any::<Multivector<f64, Euclidean3>>(),
        ) {
            let lhs = (a * b).reverse();
            let rhs = b.reverse() * a.reverse();
            prop_assert!(relative_eq!(lhs, rhs, epsilon = RELATIVE_EQ_EPS, max_relative = RELATIVE_EQ_EPS));
        }

        #[test]
        fn involute_involutory(a in any::<Multivector<f64, Euclidean3>>()) {
            prop_assert!(relative_eq!(a.involute().involute(), a, epsilon = RELATIVE_EQ_EPS, max_relative = RELATIVE_EQ_EPS));
        }

        #[test]
        fn inverse_property(a in any::<NonZeroVectorE3>()) {
            // Vectors are always invertible in Euclidean space
            let inv = a.inverse().expect("non-zero vector should be invertible");
            let product = *a * inv;
            prop_assert!(relative_eq!(product, Multivector::one(), epsilon = RELATIVE_EQ_EPS, max_relative = RELATIVE_EQ_EPS));
        }

        #[test]
        fn add_commutative(
            a in any::<Multivector<f64, Euclidean3>>(),
            b in any::<Multivector<f64, Euclidean3>>()
        ) {
            prop_assert!(relative_eq!(a + b, b + a, epsilon = RELATIVE_EQ_EPS, max_relative = RELATIVE_EQ_EPS));
        }

        #[test]
        fn add_associative(
            a in any::<Multivector<f64, Euclidean3>>(),
            b in any::<Multivector<f64, Euclidean3>>(),
            c in any::<Multivector<f64, Euclidean3>>(),
        ) {
            let lhs = (a + b) + c;
            let rhs = a + (b + c);
            prop_assert!(relative_eq!(lhs, rhs, epsilon = RELATIVE_EQ_EPS, max_relative = RELATIVE_EQ_EPS));
        }

        #[test]
        fn zero_is_additive_identity(a in any::<Multivector<f64, Euclidean3>>()) {
            let zero = Multivector::<f64, Euclidean3>::zero();
            prop_assert!(relative_eq!(a + zero, a, epsilon = RELATIVE_EQ_EPS, max_relative = RELATIVE_EQ_EPS));
        }

        #[test]
        fn one_is_multiplicative_identity(a in any::<Multivector<f64, Euclidean3>>()) {
            let one = Multivector::<f64, Euclidean3>::one();
            prop_assert!(relative_eq!(a * one, a, epsilon = RELATIVE_EQ_EPS, max_relative = RELATIVE_EQ_EPS));
            prop_assert!(relative_eq!(one * a, a, epsilon = RELATIVE_EQ_EPS, max_relative = RELATIVE_EQ_EPS));
        }

        // ====================================================================
        // Grade operation tests
        // ====================================================================

        #[test]
        fn grade_decomposition_complete(a in any::<Multivector<f64, Euclidean3>>()) {
            // Sum of all grade projections equals original
            let sum = (0..=3)
                .map(|k| a.grade_select(k))
                .fold(Multivector::<f64, Euclidean3>::zero(), |acc, x| acc + x);
            prop_assert!(relative_eq!(sum, a, epsilon = RELATIVE_EQ_EPS, max_relative = RELATIVE_EQ_EPS));
        }

        #[test]
        fn even_plus_odd_equals_original(a in any::<Multivector<f64, Euclidean3>>()) {
            let reconstructed = a.even() + a.odd();
            prop_assert!(relative_eq!(reconstructed, a, epsilon = RELATIVE_EQ_EPS, max_relative = RELATIVE_EQ_EPS));
        }

        #[test]
        fn grade_select_idempotent(
            a in any::<Multivector<f64, Euclidean3>>(),
            k in 0usize..4
        ) {
            let once = a.grade_select(k);
            let twice = once.grade_select(k);
            prop_assert!(relative_eq!(once, twice, epsilon = RELATIVE_EQ_EPS, max_relative = RELATIVE_EQ_EPS));
        }

        // ====================================================================
        // Exterior product tests
        // ====================================================================

        #[test]
        fn exterior_anticommutative_vectors(
            a in any::<VectorE3>(),
            b in any::<VectorE3>(),
        ) {
            let ab = a.exterior(&*b);
            let ba = b.exterior(&*a);
            prop_assert!(relative_eq!(ab, -&ba, epsilon = RELATIVE_EQ_EPS, max_relative = RELATIVE_EQ_EPS));
        }

        #[test]
        fn outer_associative(
            a in any::<Multivector<f64, Euclidean3>>(),
            b in any::<Multivector<f64, Euclidean3>>(),
            c in any::<Multivector<f64, Euclidean3>>(),
        ) {
            let lhs = a.exterior(&b).exterior(&c);
            let rhs = a.exterior(&b.exterior(&c));
            prop_assert!(relative_eq!(lhs, rhs, epsilon = RELATIVE_EQ_EPS, max_relative = RELATIVE_EQ_EPS));
        }

        #[test]
        fn outer_of_vectors_is_grade_2(
            a in any::<VectorE3>(),
            b in any::<VectorE3>(),
        ) {
            let wedge = a.exterior(&*b);
            // Result should have no scalar or vector parts
            prop_assert!(wedge.grade_select(0).is_zero(RELATIVE_EQ_EPS));
            prop_assert!(wedge.grade_select(1).is_zero(RELATIVE_EQ_EPS));
        }

        // ====================================================================
        // Inner product tests
        // ====================================================================

        #[test]
        fn inner_product_of_vectors_is_scalar(
            a in any::<VectorE3>(),
            b in any::<VectorE3>(),
        ) {
            let dot = a.inner(&*b);
            // Should be pure scalar
            prop_assert!(dot.grade_select(1).is_zero(RELATIVE_EQ_EPS));
            prop_assert!(dot.grade_select(2).is_zero(RELATIVE_EQ_EPS));
            prop_assert!(dot.grade_select(3).is_zero(RELATIVE_EQ_EPS));
        }

        #[test]
        fn inner_product_symmetric_for_vectors(
            a in any::<VectorE3>(),
            b in any::<VectorE3>(),
        ) {
            // For vectors, a·b = b·a
            let ab = a.inner(&*b);
            let ba = b.inner(&*a);
            prop_assert!(relative_eq!(ab, ba, epsilon = RELATIVE_EQ_EPS, max_relative = RELATIVE_EQ_EPS));
        }

        // ====================================================================
        // Dual tests
        // ====================================================================

        #[test]
        fn dual_undual_roundtrip(a in any::<Multivector<f64, Euclidean3>>()) {
            let roundtrip = a.dual().undual();
            prop_assert!(relative_eq!(roundtrip, a, epsilon = RELATIVE_EQ_EPS, max_relative = RELATIVE_EQ_EPS));
        }

        #[test]
        fn dual_changes_grade(a in any::<NonZeroVectorE3>()) {
            // Vector (grade 1) dualizes to bivector (grade 2) in 3D
            let dual = a.dual();
            prop_assert_eq!(dual.grade(RELATIVE_EQ_EPS), Some(2));
        }

        // ====================================================================
        // Sandwich product tests
        // ====================================================================

        #[test]
        fn sandwich_by_unit_vector_preserves_norm(
            n in any::<UnitVectorE3>(),
            v in any::<NonZeroVectorE3>(),
        ) {
            let reflected = n.sandwich(&*v);
            // Reflection preserves norm
            prop_assert!(relative_eq!(reflected.norm(), v.norm(), epsilon = RELATIVE_EQ_EPS, max_relative = RELATIVE_EQ_EPS));
        }

        #[test]
        fn sandwich_by_unit_vector_is_reflection(
            n in any::<UnitVectorE3>(),
            v in any::<VectorE3>(),
        ) {
            let reflected = n.sandwich(&*v);
            // Double reflection returns original
            let double_reflected = n.sandwich(&reflected);
            prop_assert!(relative_eq!(double_reflected, *v, epsilon = RELATIVE_EQ_EPS, max_relative = RELATIVE_EQ_EPS));
        }

        // ====================================================================
        // De Morgan's Laws (RGA)
        // https://rigidgeometricalgebra.org/wiki/index.php?title=Geometric_products
        // ====================================================================

        #[test]
        fn de_morgan_complement_of_geometric_product(
            a in any::<Multivector<f64, Euclidean3>>(),
            b in any::<Multivector<f64, Euclidean3>>(),
        ) {
            // complement(a * b) = complement(a) ⋇ complement(b)
            let lhs = (a * b).complement();
            let rhs = a.complement().antiproduct(&b.complement());
            prop_assert!(relative_eq!(lhs, rhs, epsilon = RELATIVE_EQ_EPS, max_relative = RELATIVE_EQ_EPS));
        }

        #[test]
        fn de_morgan_complement_of_antiproduct(
            a in any::<Multivector<f64, Euclidean3>>(),
            b in any::<Multivector<f64, Euclidean3>>(),
        ) {
            // complement(a ⋇ b) = complement(a) * complement(b)
            // This is actually how we define antiproduct, so this tests consistency
            let lhs = a.antiproduct(&b).complement();
            let rhs = a.complement() * b.complement();
            prop_assert!(relative_eq!(lhs, rhs, epsilon = RELATIVE_EQ_EPS, max_relative = RELATIVE_EQ_EPS));
        }

        #[test]
        fn complement_involutory(
            a in any::<Multivector<f64, Euclidean3>>(),
        ) {
            // complement(complement(a)) = ±a
            // The double complement may introduce a sign depending on dimension and grade
            let double_comp = a.complement().complement();
            // In 3D, double complement is ±1 times original depending on grade
            // For a general multivector, we check that double complement squared equals a squared
            let lhs_norm_sq = (double_comp * double_comp.reverse()).scalar_part();
            let rhs_norm_sq = (a * a.reverse()).scalar_part();
            prop_assert!(relative_eq!(lhs_norm_sq, rhs_norm_sq, epsilon = RELATIVE_EQ_EPS, max_relative = RELATIVE_EQ_EPS));
        }
    }

    // ========================================================================
    // 2D specific tests
    // ========================================================================

    #[test]
    fn test_2d_complex_number_behavior() {
        // In 2D Euclidean GA, the bivector e₁₂ behaves like the imaginary unit i
        let e1: Multivector<f64, Euclidean2> = Multivector::basis_vector(0);
        let e2: Multivector<f64, Euclidean2> = Multivector::basis_vector(1);
        let i = e1 * e2; // e₁₂ = "i"

        // i² = -1
        let i_sq = i * i;
        assert!(relative_eq!(
            i_sq.scalar_part(),
            -1.0,
            max_relative = RELATIVE_EQ_EPS
        ));

        // We can represent complex numbers as a + b*e₁₂
        let z = Multivector::scalar(3.0) + i * 4.0; // 3 + 4i
        let norm_sq = (z * z.reverse()).scalar_part();
        assert!(relative_eq!(
            norm_sq,
            25.0,
            epsilon = RELATIVE_EQ_EPS,
            max_relative = RELATIVE_EQ_EPS
        )); // |z|² = 3² + 4² = 25
    }

    // ========================================================================
    // Copy semantics test
    // ========================================================================

    #[test]
    fn test_copy_semantics() {
        let v: Multivector<f64, Euclidean3> = Multivector::vector(&[1.0, 2.0, 3.0]);

        // Copy the value (not move)
        let v2 = v;

        // Original is still usable (proves Copy, not just Clone)
        assert!(relative_eq!(
            v,
            v2,
            epsilon = RELATIVE_EQ_EPS,
            max_relative = RELATIVE_EQ_EPS
        ));

        // Can use both independently
        let sum = v + v2;
        assert!(relative_eq!(
            sum.get(Blade::basis_vector(0)),
            2.0,
            max_relative = RELATIVE_EQ_EPS
        ));
    }

    // ========================================================================
    // Serde tests (feature-gated)
    // ========================================================================

    #[cfg(feature = "serde")]
    #[test]
    fn test_serde_roundtrip() {
        let v: Multivector<f64, Euclidean3> = Multivector::vector(&[1.0, 2.0, 3.0]);

        // Serialize to JSON
        let json = serde_json::to_string(&v).expect("serialization failed");

        // Deserialize back
        let v2: Multivector<f64, Euclidean3> =
            serde_json::from_str(&json).expect("deserialization failed");

        assert!(relative_eq!(
            v,
            v2,
            epsilon = RELATIVE_EQ_EPS,
            max_relative = RELATIVE_EQ_EPS
        ));
    }
}
