//! Grade utility functions for blade algebra.
//!
//! In geometric algebra, the **grade** of a blade is the number of basis
//! vectors it contains. For example:
//! - Scalars have grade 0
//! - Vectors have grade 1
//! - Bivectors have grade 2
//! - etc.
//!
//! This module provides utilities for:
//! - Computing blade grades
//! - Enumerating blades by grade
//! - Computing binomial coefficients for grade counting

/// Returns the grade (number of basis vectors) of a blade.
///
/// The grade equals the number of 1-bits in the blade's bitmask index.
///
/// # Example
///
/// ```
/// use clifford_codegen::algebra::grade;
///
/// assert_eq!(grade(0b000), 0); // scalar
/// assert_eq!(grade(0b001), 1); // e1
/// assert_eq!(grade(0b011), 2); // e12
/// assert_eq!(grade(0b111), 3); // e123
/// ```
#[inline]
pub const fn grade(blade: usize) -> usize {
    blade.count_ones() as usize
}

/// Returns all blade indices of a given grade in an n-dimensional algebra.
///
/// Blades are returned in ascending index order (canonical ordering).
///
/// # Arguments
///
/// * `dim` - The dimension of the algebra (number of basis vectors)
/// * `g` - The target grade
///
/// # Returns
///
/// A vector of blade indices with the specified grade.
///
/// # Example
///
/// ```
/// use clifford_codegen::algebra::blades_of_grade;
///
/// // 3D algebra: grade-2 blades (bivectors)
/// let bivectors = blades_of_grade(3, 2);
/// assert_eq!(bivectors, vec![0b011, 0b101, 0b110]); // e12, e13, e23
///
/// // 3D algebra: grade-0 blade (scalar)
/// let scalars = blades_of_grade(3, 0);
/// assert_eq!(scalars, vec![0]);
///
/// // 3D algebra: grade-3 blade (pseudoscalar)
/// let pseudoscalar = blades_of_grade(3, 3);
/// assert_eq!(pseudoscalar, vec![0b111]); // e123
/// ```
pub fn blades_of_grade(dim: usize, g: usize) -> Vec<usize> {
    (0..(1 << dim)).filter(|&b| grade(b) == g).collect()
}

/// Returns all blade indices for the given grades in an n-dimensional algebra.
///
/// Blades are returned grouped by grade, with each grade's blades in
/// ascending index order.
///
/// # Arguments
///
/// * `dim` - The dimension of the algebra (number of basis vectors)
/// * `grades` - The target grades
///
/// # Example
///
/// ```
/// use clifford_codegen::algebra::blades_of_grades;
///
/// // 3D even subalgebra: grades 0 and 2
/// let even = blades_of_grades(3, &[0, 2]);
/// assert_eq!(even, vec![0, 0b011, 0b101, 0b110]); // 1, e12, e13, e23
/// ```
pub fn blades_of_grades(dim: usize, grades: &[usize]) -> Vec<usize> {
    grades
        .iter()
        .flat_map(|&g| blades_of_grade(dim, g))
        .collect()
}

/// Computes the binomial coefficient C(n, k) = n! / (k! * (n-k)!).
///
/// This gives the number of blades of grade k in an n-dimensional algebra.
///
/// # Example
///
/// ```
/// use clifford_codegen::algebra::binomial;
///
/// // Number of bivectors in 3D: C(3, 2) = 3
/// assert_eq!(binomial(3, 2), 3);
///
/// // Number of vectors in 4D: C(4, 1) = 4
/// assert_eq!(binomial(4, 1), 4);
///
/// // Number of scalars: always 1
/// assert_eq!(binomial(5, 0), 1);
///
/// // Edge cases
/// assert_eq!(binomial(5, 5), 1);
/// assert_eq!(binomial(5, 6), 0); // k > n
/// ```
pub const fn binomial(n: usize, k: usize) -> usize {
    if k > n {
        return 0;
    }
    if k == 0 || k == n {
        return 1;
    }

    // Use smaller k for efficiency: C(n, k) = C(n, n-k)
    let k = if k > n - k { n - k } else { k };

    let mut result = 1;
    let mut i = 0;
    while i < k {
        result = result * (n - i) / (i + 1);
        i += 1;
    }
    result
}

/// Grade selection for outer product (wedge).
///
/// The outer product `a ∧ b` has grade = grade(a) + grade(b), but only
/// if the result doesn't exceed the algebra dimension. Otherwise the
/// product is zero.
///
/// # Example
///
/// ```
/// use clifford_codegen::algebra::outer_grade;
///
/// // In 3D: bivector ∧ vector = trivector (grade 3)
/// assert_eq!(outer_grade(2, 1, 3), Some(3));
///
/// // In 3D: bivector ∧ bivector = 0 (grade 4 > dim 3)
/// assert_eq!(outer_grade(2, 2, 3), None);
/// ```
pub const fn outer_grade(grade_a: usize, grade_b: usize, dim: usize) -> Option<usize> {
    let result = grade_a + grade_b;
    if result <= dim { Some(result) } else { None }
}

/// Grade selection for left contraction.
///
/// The left contraction `a ⌋ b` has grade = grade(b) - grade(a), but only
/// if grade(a) ≤ grade(b). Otherwise the product is zero.
///
/// # Example
///
/// ```
/// use clifford_codegen::algebra::left_contraction_grade;
///
/// // vector ⌋ bivector = vector (grade 1)
/// assert_eq!(left_contraction_grade(1, 2), Some(1));
///
/// // bivector ⌋ vector = 0 (grade 2 > grade 1)
/// assert_eq!(left_contraction_grade(2, 1), None);
/// ```
pub const fn left_contraction_grade(grade_a: usize, grade_b: usize) -> Option<usize> {
    if grade_a <= grade_b {
        Some(grade_b - grade_a)
    } else {
        None
    }
}

/// Grade selection for inner product (symmetric).
///
/// The inner product `a · b` has grade = |grade(a) - grade(b)|.
///
/// # Example
///
/// ```
/// use clifford_codegen::algebra::inner_grade;
///
/// // vector · vector = scalar (grade 0)
/// assert_eq!(inner_grade(1, 1), 0);
///
/// // bivector · vector = vector (grade 1)
/// assert_eq!(inner_grade(2, 1), 1);
/// ```
pub const fn inner_grade(grade_a: usize, grade_b: usize) -> usize {
    grade_a.abs_diff(grade_b)
}

/// Returns all grades present in the geometric product of two blades.
///
/// The geometric product of a grade-a blade with a grade-b blade can
/// produce blades of grades from |a - b| to a + b, stepping by 2
/// (preserving parity).
///
/// # Example
///
/// ```
/// use clifford_codegen::algebra::geometric_grades;
///
/// // vector * vector = scalar + bivector (grades 0, 2)
/// assert_eq!(geometric_grades(1, 1, 4), vec![0, 2]);
///
/// // bivector * vector = vector + trivector (grades 1, 3)
/// assert_eq!(geometric_grades(2, 1, 4), vec![1, 3]);
/// ```
pub fn geometric_grades(grade_a: usize, grade_b: usize, dim: usize) -> Vec<usize> {
    let min_grade = grade_a.abs_diff(grade_b);
    let max_grade = (grade_a + grade_b).min(dim);

    (min_grade..=max_grade).step_by(2).collect()
}

/// Sign factor for the reverse of a k-blade.
///
/// The reverse operation flips the order of basis vectors in a blade,
/// introducing a sign of `(-1)^(k(k-1)/2)` where k is the grade.
///
/// The pattern is: `++--++--...` (repeating every 4 grades)
///
/// | Grade | k(k-1)/2 | Sign |
/// |-------|----------|------|
/// | 0     | 0        | +1   |
/// | 1     | 0        | +1   |
/// | 2     | 1        | -1   |
/// | 3     | 3        | -1   |
/// | 4     | 6        | +1   |
/// | 5     | 10       | +1   |
/// | 6     | 15       | -1   |
/// | 7     | 21       | -1   |
///
/// # Example
///
/// ```
/// use clifford_codegen::algebra::reverse_sign;
///
/// assert_eq!(reverse_sign(0), 1);  // scalar: +
/// assert_eq!(reverse_sign(1), 1);  // vector: +
/// assert_eq!(reverse_sign(2), -1); // bivector: -
/// assert_eq!(reverse_sign(3), -1); // trivector: -
/// assert_eq!(reverse_sign(4), 1);  // 4-vector: +
/// ```
#[inline]
pub const fn reverse_sign(grade: usize) -> i8 {
    // (-1)^(k(k-1)/2)
    let exponent = (grade * grade.saturating_sub(1)) / 2;
    if exponent.is_multiple_of(2) { 1 } else { -1 }
}

/// Sign factor for the antireverse of a k-blade in an n-dimensional algebra.
///
/// The antireverse is the dual of the reverse (or equivalently, the reverse
/// of the dual). For a k-blade in n dimensions, the sign is the reverse sign
/// of the dual grade (n - k).
///
/// # Arguments
///
/// * `grade` - The grade of the blade
/// * `dim` - The dimension of the algebra
///
/// # Example
///
/// ```
/// use clifford_codegen::algebra::antireverse_sign;
///
/// // In 3D algebra:
/// assert_eq!(antireverse_sign(0, 3), -1); // dual is grade 3, reverse sign is -1
/// assert_eq!(antireverse_sign(1, 3), -1); // dual is grade 2, reverse sign is -1
/// assert_eq!(antireverse_sign(2, 3), 1);  // dual is grade 1, reverse sign is +1
/// assert_eq!(antireverse_sign(3, 3), 1);  // dual is grade 0, reverse sign is +1
/// ```
#[inline]
pub const fn antireverse_sign(grade: usize, dim: usize) -> i8 {
    reverse_sign(dim - grade)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn grade_from_popcount() {
        assert_eq!(grade(0b0000), 0);
        assert_eq!(grade(0b0001), 1);
        assert_eq!(grade(0b0011), 2);
        assert_eq!(grade(0b0111), 3);
        assert_eq!(grade(0b1111), 4);
        assert_eq!(grade(0b10101), 3);
    }

    #[test]
    fn blades_of_grade_3d() {
        // Grade 0: scalar
        assert_eq!(blades_of_grade(3, 0), vec![0]);

        // Grade 1: vectors
        assert_eq!(blades_of_grade(3, 1), vec![1, 2, 4]);

        // Grade 2: bivectors
        assert_eq!(blades_of_grade(3, 2), vec![3, 5, 6]);

        // Grade 3: pseudoscalar
        assert_eq!(blades_of_grade(3, 3), vec![7]);

        // Grade 4: empty (exceeds dimension)
        assert_eq!(blades_of_grade(3, 4), Vec::<usize>::new());
    }

    #[test]
    fn blades_of_grades_even_subalgebra() {
        // Even subalgebra in 3D: grades 0 and 2
        let even = blades_of_grades(3, &[0, 2]);
        assert_eq!(even, vec![0, 3, 5, 6]);
    }

    #[test]
    fn binomial_coefficients() {
        // Pascal's triangle row 4: 1 4 6 4 1
        assert_eq!(binomial(4, 0), 1);
        assert_eq!(binomial(4, 1), 4);
        assert_eq!(binomial(4, 2), 6);
        assert_eq!(binomial(4, 3), 4);
        assert_eq!(binomial(4, 4), 1);

        // Edge cases
        assert_eq!(binomial(0, 0), 1);
        assert_eq!(binomial(5, 6), 0);
    }

    #[test]
    fn binomial_counts_blades() {
        // Number of blades of each grade equals binomial coefficient
        for dim in 0..=6 {
            for g in 0..=dim {
                let count = blades_of_grade(dim, g).len();
                assert_eq!(count, binomial(dim, g), "dim={}, grade={}", dim, g);
            }
        }
    }

    #[test]
    fn total_blades_is_power_of_two() {
        // Total blades = sum of C(n,k) for k=0..n = 2^n
        for dim in 0..=6 {
            let total: usize = (0..=dim).map(|g| binomial(dim, g)).sum();
            assert_eq!(total, 1 << dim);
        }
    }

    #[test]
    fn outer_grade_rules() {
        // Valid outer products
        assert_eq!(outer_grade(0, 1, 3), Some(1)); // scalar ∧ vector = vector
        assert_eq!(outer_grade(1, 1, 3), Some(2)); // vector ∧ vector = bivector
        assert_eq!(outer_grade(1, 2, 3), Some(3)); // vector ∧ bivector = trivector

        // Invalid (exceeds dimension)
        assert_eq!(outer_grade(2, 2, 3), None); // grade 4 > dim 3
        assert_eq!(outer_grade(3, 1, 3), None); // grade 4 > dim 3
    }

    #[test]
    fn left_contraction_grade_rules() {
        assert_eq!(left_contraction_grade(0, 2), Some(2)); // scalar ⌋ bivector = bivector
        assert_eq!(left_contraction_grade(1, 2), Some(1)); // vector ⌋ bivector = vector
        assert_eq!(left_contraction_grade(2, 2), Some(0)); // bivector ⌋ bivector = scalar
        assert_eq!(left_contraction_grade(3, 2), None); // trivector ⌋ bivector = 0
    }

    #[test]
    fn geometric_product_grades() {
        // vector * vector = scalar + bivector
        assert_eq!(geometric_grades(1, 1, 4), vec![0, 2]);

        // bivector * vector = vector + trivector
        assert_eq!(geometric_grades(2, 1, 4), vec![1, 3]);

        // bivector * bivector = scalar + bivector + 4-vector
        assert_eq!(geometric_grades(2, 2, 4), vec![0, 2, 4]);

        // Capped by dimension
        assert_eq!(geometric_grades(2, 2, 3), vec![0, 2]);
    }

    #[test]
    fn reverse_sign_pattern() {
        // Pattern: ++--++--...
        assert_eq!(reverse_sign(0), 1);
        assert_eq!(reverse_sign(1), 1);
        assert_eq!(reverse_sign(2), -1);
        assert_eq!(reverse_sign(3), -1);
        assert_eq!(reverse_sign(4), 1);
        assert_eq!(reverse_sign(5), 1);
        assert_eq!(reverse_sign(6), -1);
        assert_eq!(reverse_sign(7), -1);
        assert_eq!(reverse_sign(8), 1);
    }

    #[test]
    fn antireverse_is_dual_of_reverse() {
        // antireverse_sign(k, n) = reverse_sign(n - k)
        for dim in 1..=6 {
            for k in 0..=dim {
                let dual_grade = dim - k;
                assert_eq!(
                    antireverse_sign(k, dim),
                    reverse_sign(dual_grade),
                    "dim={}, k={}, dual_grade={}",
                    dim,
                    k,
                    dual_grade
                );
            }
        }
    }

    #[test]
    fn antireverse_3d() {
        // In 3D algebra
        assert_eq!(antireverse_sign(0, 3), -1); // dual grade 3 -> reverse_sign(3) = -1
        assert_eq!(antireverse_sign(1, 3), -1); // dual grade 2 -> reverse_sign(2) = -1
        assert_eq!(antireverse_sign(2, 3), 1); // dual grade 1 -> reverse_sign(1) = +1
        assert_eq!(antireverse_sign(3, 3), 1); // dual grade 0 -> reverse_sign(0) = +1
    }
}
