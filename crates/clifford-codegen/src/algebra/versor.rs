//! Versor identification and classification.
//!
//! A versor in Geometric Algebra is an element that can be expressed as
//! a geometric product of invertible vectors. Versors are fundamental for
//! representing transformations via the sandwich product: `X' = V * X * rev(V)`.
//!
//! # Versor Types
//!
//! | Type | Grades | Parity | Transformation |
//! |------|--------|--------|----------------|
//! | Unit Vector | [1] | Odd | Reflection |
//! | Rotor (2D/3D) | [0, 2] | Even | Rotation |
//! | Motor (PGA) | [0, 2, 4] | Even | Rigid motion |
//! | Flector (PGA) | [1, 3] | Odd | Reflection + translation |
//!
//! # Properties
//!
//! 1. **Grade Parity**: All grades have the same parity (all even or all odd)
//! 2. **Versor Constraint**: `V * rev(V) = scalar` (or pseudoscalar for odd versors)
//! 3. **Closure**: Even * Even = Even, Even * Odd = Odd

/// Versor classification by grade parity.
///
/// Versors are classified as even or odd based on the parity of their grades.
/// Even versors (rotors, motors) preserve orientation; odd versors (reflectors,
/// flectors) reverse orientation.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum VersorParity {
    /// Even versor: grades 0, 2, 4, ... (rotors, motors).
    ///
    /// Even versors are products of an even number of vectors and preserve
    /// orientation under the sandwich product.
    Even,

    /// Odd versor: grades 1, 3, 5, ... (reflectors, flectors).
    ///
    /// Odd versors are products of an odd number of vectors and reverse
    /// orientation under the sandwich product.
    Odd,
}

impl VersorParity {
    /// Returns the parity value (0 for even, 1 for odd).
    #[inline]
    pub fn value(&self) -> usize {
        match self {
            VersorParity::Even => 0,
            VersorParity::Odd => 1,
        }
    }

    /// Returns the result of composing two versors.
    ///
    /// Even * Even = Even, Odd * Odd = Even, Even * Odd = Odd.
    #[inline]
    pub fn compose(self, other: VersorParity) -> VersorParity {
        if self.value() == other.value() {
            VersorParity::Even
        } else {
            VersorParity::Odd
        }
    }
}

/// Determines the versor parity of a type based on its grades.
///
/// Returns `Some(parity)` if all grades have the same parity, `None` otherwise.
///
/// # Examples
///
/// ```ignore
/// // Even versors
/// assert_eq!(versor_parity(&[0, 2]), Some(VersorParity::Even));    // Rotor
/// assert_eq!(versor_parity(&[0, 2, 4]), Some(VersorParity::Even)); // Motor
///
/// // Odd versors
/// assert_eq!(versor_parity(&[1]), Some(VersorParity::Odd));        // Vector
/// assert_eq!(versor_parity(&[1, 3]), Some(VersorParity::Odd));     // Flector
///
/// // Not a versor (mixed parity)
/// assert_eq!(versor_parity(&[0, 1, 2]), None);
/// ```
pub fn versor_parity(grades: &[usize]) -> Option<VersorParity> {
    if grades.is_empty() {
        return None;
    }

    let first_parity = grades[0] % 2;
    if grades.iter().all(|&g| g % 2 == first_parity) {
        Some(if first_parity == 0 {
            VersorParity::Even
        } else {
            VersorParity::Odd
        })
    } else {
        None
    }
}

/// Information about a verified versor type.
#[derive(Debug, Clone)]
pub struct VersorInfo {
    /// The parity of the versor (even or odd).
    pub parity: VersorParity,

    /// Whether this is a unit versor (`V * rev(V) = 1`).
    ///
    /// Unit versors have a simpler sandwich product formula since
    /// `V⁻¹ = rev(V)` for unit versors.
    pub is_unit: bool,

    /// Types that this versor can transform via sandwich product.
    ///
    /// A type can be transformed if the sandwich product `V * X * rev(V)`
    /// produces output of the same grades as the input.
    pub sandwich_targets: Vec<String>,
}

impl Default for VersorInfo {
    fn default() -> Self {
        Self {
            parity: VersorParity::Even,
            is_unit: false,
            sandwich_targets: Vec::new(),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn even_versor_identification() {
        // Scalar (grade 0)
        assert_eq!(versor_parity(&[0]), Some(VersorParity::Even));

        // Rotor (grades 0, 2)
        assert_eq!(versor_parity(&[0, 2]), Some(VersorParity::Even));

        // Motor (grades 0, 2, 4)
        assert_eq!(versor_parity(&[0, 2, 4]), Some(VersorParity::Even));

        // Even multivector (grades 0, 2, 4, 6)
        assert_eq!(versor_parity(&[0, 2, 4, 6]), Some(VersorParity::Even));
    }

    #[test]
    fn odd_versor_identification() {
        // Vector (grade 1)
        assert_eq!(versor_parity(&[1]), Some(VersorParity::Odd));

        // Flector (grades 1, 3)
        assert_eq!(versor_parity(&[1, 3]), Some(VersorParity::Odd));

        // Odd multivector (grades 1, 3, 5)
        assert_eq!(versor_parity(&[1, 3, 5]), Some(VersorParity::Odd));
    }

    #[test]
    fn mixed_parity_not_versor() {
        // Mixed grades are not versors
        assert_eq!(versor_parity(&[0, 1]), None);
        assert_eq!(versor_parity(&[0, 1, 2]), None);
        assert_eq!(versor_parity(&[1, 2]), None);
        assert_eq!(versor_parity(&[0, 2, 3]), None);
    }

    #[test]
    fn empty_grades_not_versor() {
        assert_eq!(versor_parity(&[]), None);
    }

    #[test]
    fn parity_composition() {
        // Even * Even = Even
        assert_eq!(
            VersorParity::Even.compose(VersorParity::Even),
            VersorParity::Even
        );

        // Odd * Odd = Even
        assert_eq!(
            VersorParity::Odd.compose(VersorParity::Odd),
            VersorParity::Even
        );

        // Even * Odd = Odd
        assert_eq!(
            VersorParity::Even.compose(VersorParity::Odd),
            VersorParity::Odd
        );

        // Odd * Even = Odd
        assert_eq!(
            VersorParity::Odd.compose(VersorParity::Even),
            VersorParity::Odd
        );
    }
}
