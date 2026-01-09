//! Blade indexing and grade utilities.
//!
//! This module provides functions for working with basis blade indices.
//! In this library, blades are indexed using their binary representation:
//! each bit indicates whether a particular basis vector is present.
//!
//! # Indexing Scheme
//!
//! For a basis blade, its index is a bitmask where bit `i` is set if
//! basis vector `e_i` is present. For example, in 3D:
//!
//! | Index | Binary | Blade |
//! |-------|--------|-------|
//! | 0 | `000` | `1` (scalar) |
//! | 1 | `001` | `e₁` |
//! | 2 | `010` | `e₂` |
//! | 3 | `011` | `e₁₂` |
//! | 4 | `100` | `e₃` |
//! | 5 | `101` | `e₁₃` |
//! | 6 | `110` | `e₂₃` |
//! | 7 | `111` | `e₁₂₃` |
//!
//! # Grade
//!
//! The **grade** of a blade is the number of basis vectors it contains,
//! which equals the number of set bits (popcount) in its index.

/// Returns the grade (number of basis vectors) of a blade given its index.
///
/// The grade equals the population count (number of 1-bits) in the binary
/// representation of the index.
///
/// # Arguments
///
/// * `index` - The blade index (bitmask of constituent basis vectors)
///
/// # Returns
///
/// The grade of the blade, from 0 (scalar) to N (pseudoscalar in N dimensions).
///
/// # Examples
///
/// ```
/// use clifford::basis::grade_of_blade;
///
/// assert_eq!(grade_of_blade(0b000), 0); // scalar
/// assert_eq!(grade_of_blade(0b001), 1); // e₁ (vector)
/// assert_eq!(grade_of_blade(0b011), 2); // e₁₂ (bivector)
/// assert_eq!(grade_of_blade(0b111), 3); // e₁₂₃ (trivector)
/// ```
#[inline]
pub const fn grade_of_blade(index: usize) -> usize {
    index.count_ones() as usize
}

/// Returns the number of blades of a given grade in a space of given dimension.
///
/// This computes the binomial coefficient `C(dim, grade)`, which counts
/// how many ways we can choose `grade` basis vectors from `dim` total.
///
/// # Arguments
///
/// * `dim` - The dimension of the vector space
/// * `grade` - The grade to count blades for
///
/// # Returns
///
/// The number of basis blades of the specified grade.
///
/// # Examples
///
/// ```
/// use clifford::basis::blades_of_grade;
///
/// // In 3D, there is 1 scalar, 3 vectors, 3 bivectors, 1 trivector
/// assert_eq!(blades_of_grade(3, 0), 1);
/// assert_eq!(blades_of_grade(3, 1), 3);
/// assert_eq!(blades_of_grade(3, 2), 3);
/// assert_eq!(blades_of_grade(3, 3), 1);
///
/// // In 4D, there are 6 bivectors (4 choose 2)
/// assert_eq!(blades_of_grade(4, 2), 6);
/// ```
#[inline]
pub const fn blades_of_grade(dim: usize, grade: usize) -> usize {
    if grade > dim {
        return 0;
    }
    binomial(dim, grade)
}

/// Computes the binomial coefficient C(n, k) = n! / (k! * (n-k)!).
///
/// Uses an iterative algorithm that avoids overflow for reasonable inputs.
///
/// # Arguments
///
/// * `n` - The total number of items
/// * `k` - The number of items to choose
///
/// # Returns
///
/// The number of ways to choose `k` items from `n` items.
#[inline]
pub const fn binomial(n: usize, k: usize) -> usize {
    if k > n {
        return 0;
    }
    if k == 0 || k == n {
        return 1;
    }

    // Use symmetry: C(n,k) = C(n, n-k)
    // Choose the smaller k to minimize iterations
    let k = if k > n - k { n - k } else { k };

    let mut result: usize = 1;
    let mut i: usize = 0;
    while i < k {
        result = result * (n - i) / (i + 1);
        i += 1;
    }
    result
}

/// Returns the index of the first blade of a given grade.
///
/// Blades are ordered by grade, then by index within grade.
/// This returns the starting index for a particular grade in that ordering.
///
/// # Arguments
///
/// * `dim` - The dimension of the vector space
/// * `grade` - The grade to find the start index for
///
/// # Returns
///
/// The cumulative count of all blades with grade less than `grade`.
///
/// # Examples
///
/// ```
/// use clifford::basis::grade_start_index;
///
/// // In 3D: grades 0,1,2,3 have 1,3,3,1 blades
/// assert_eq!(grade_start_index(3, 0), 0);  // scalars start at 0
/// assert_eq!(grade_start_index(3, 1), 1);  // vectors start at 1
/// assert_eq!(grade_start_index(3, 2), 4);  // bivectors start at 4
/// assert_eq!(grade_start_index(3, 3), 7);  // trivectors start at 7
/// ```
#[inline]
pub const fn grade_start_index(dim: usize, grade: usize) -> usize {
    let mut sum: usize = 0;
    let mut g: usize = 0;
    while g < grade {
        sum += blades_of_grade(dim, g);
        g += 1;
    }
    sum
}

#[cfg(test)]
mod tests {
    use super::*;
    use proptest::prelude::*;

    proptest! {
        /// Grade equals population count for any blade index.
        #[test]
        fn grade_is_popcount(index in 0usize..256) {
            prop_assert_eq!(grade_of_blade(index), index.count_ones() as usize);
        }

        /// Binomial coefficient symmetry: C(n,k) = C(n, n-k)
        #[test]
        fn binomial_symmetry(n in 0usize..20, k in 0usize..20) {
            if k <= n {
                prop_assert_eq!(binomial(n, k), binomial(n, n - k));
            }
        }

        /// Sum of binomial coefficients: sum of C(n,k) for k=0..n equals 2^n
        #[test]
        fn binomial_sum(n in 0usize..16) {
            let sum: usize = (0..=n).map(|k| binomial(n, k)).sum();
            prop_assert_eq!(sum, 1 << n);
        }
    }

    #[test]
    fn known_binomial_values() {
        // Pascal's triangle values
        assert_eq!(binomial(5, 0), 1);
        assert_eq!(binomial(5, 1), 5);
        assert_eq!(binomial(5, 2), 10);
        assert_eq!(binomial(5, 3), 10);
        assert_eq!(binomial(5, 4), 5);
        assert_eq!(binomial(5, 5), 1);
    }

    #[test]
    fn blades_of_grade_3d() {
        // 3D has 1 scalar, 3 vectors, 3 bivectors, 1 trivector
        assert_eq!(blades_of_grade(3, 0), 1);
        assert_eq!(blades_of_grade(3, 1), 3);
        assert_eq!(blades_of_grade(3, 2), 3);
        assert_eq!(blades_of_grade(3, 3), 1);
        assert_eq!(blades_of_grade(3, 4), 0); // no grade-4 blades in 3D
    }

    #[test]
    fn grade_start_indices() {
        // 3D cumulative: 0, 1, 4, 7, 8
        assert_eq!(grade_start_index(3, 0), 0);
        assert_eq!(grade_start_index(3, 1), 1);
        assert_eq!(grade_start_index(3, 2), 4);
        assert_eq!(grade_start_index(3, 3), 7);
    }
}
