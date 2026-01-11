//! Property-based tests for blade algebra.
//!
//! These tests verify the algebraic properties of the blade algebra engine
//! using proptest to ensure correctness across a wide range of inputs.

#![cfg(test)]

use proptest::prelude::*;

use super::grade::grade;
use super::sign::basis_product;
use super::{Algebra, ProductTable};

/// Strategy for blade indices in a 4D algebra (0..16).
fn blade_4d() -> impl Strategy<Value = usize> {
    0usize..16
}

/// Strategy for blade indices in a 6D algebra (0..64).
fn blade_6d() -> impl Strategy<Value = usize> {
    0usize..64
}

/// Strategy for basis vector indices in 6D (0..6).
fn basis_6d() -> impl Strategy<Value = usize> {
    0usize..6
}

proptest! {
    /// The result blade of a product is always the XOR of the inputs.
    #[test]
    fn product_result_is_xor(a in blade_6d(), b in blade_6d()) {
        let (_, result) = basis_product(a, b, |_| 1);
        prop_assert_eq!(result, a ^ b);
    }

    /// The geometric product is associative: (a * b) * c = a * (b * c).
    #[test]
    fn product_is_associative(
        a in blade_4d(),
        b in blade_4d(),
        c in blade_4d(),
    ) {
        let metric = |_| 1i8; // Euclidean

        let (sign_ab, ab) = basis_product(a, b, metric);
        let (sign_abc1, abc1) = basis_product(ab, c, metric);
        let sign1 = sign_ab as i32 * sign_abc1 as i32;

        let (sign_bc, bc) = basis_product(b, c, metric);
        let (sign_abc2, abc2) = basis_product(a, bc, metric);
        let sign2 = sign_bc as i32 * sign_abc2 as i32;

        prop_assert_eq!(abc1, abc2, "result blades differ");
        prop_assert_eq!(sign1, sign2, "signs differ");
    }

    /// Scalar multiplication is identity: 1 * a = a * 1 = a.
    #[test]
    fn scalar_is_identity(a in blade_6d()) {
        let metric = |_| 1i8;

        let (sign1, result1) = basis_product(0, a, metric);
        prop_assert_eq!(sign1, 1);
        prop_assert_eq!(result1, a);

        let (sign2, result2) = basis_product(a, 0, metric);
        prop_assert_eq!(sign2, 1);
        prop_assert_eq!(result2, a);
    }

    /// Distinct basis vectors anticommute: eᵢeⱼ = -eⱼeᵢ for i ≠ j.
    #[test]
    fn basis_vectors_anticommute(i in basis_6d(), j in basis_6d()) {
        prop_assume!(i != j);
        let metric = |_| 1i8;

        let a = 1 << i;
        let b = 1 << j;

        let (sign_ab, result_ab) = basis_product(a, b, metric);
        let (sign_ba, result_ba) = basis_product(b, a, metric);

        prop_assert_eq!(result_ab, result_ba, "result blades should be the same");
        prop_assert_eq!(sign_ab, -sign_ba, "signs should be opposite");
    }

    /// The grade of the result is bounded: |grade(a) - grade(b)| ≤ grade(result) ≤ grade(a) + grade(b).
    #[test]
    fn product_grade_bounds(a in blade_6d(), b in blade_6d()) {
        let (sign, result) = basis_product(a, b, |_| 1);

        // Only check bounds for non-zero products
        if sign != 0 {
            let ga = grade(a);
            let gb = grade(b);
            let gr = grade(result);

            let min_grade = if ga >= gb { ga - gb } else { gb - ga };
            let max_grade = ga + gb;

            prop_assert!(gr >= min_grade, "grade {} < min {}", gr, min_grade);
            prop_assert!(gr <= max_grade, "grade {} > max {}", gr, max_grade);
        }
    }

    /// The result grade has the same parity as the sum of input grades.
    #[test]
    fn product_grade_parity(a in blade_6d(), b in blade_6d()) {
        let (_, result) = basis_product(a, b, |_| 1);

        let ga = grade(a);
        let gb = grade(b);
        let gr = grade(result);

        prop_assert_eq!((ga + gb) % 2, gr % 2, "parity mismatch");
    }

    /// A blade times itself squares to a scalar (grade 0).
    #[test]
    fn blade_squared_is_scalar(a in blade_6d()) {
        let (_, result) = basis_product(a, a, |_| 1);
        prop_assert_eq!(grade(result), 0, "blade squared should give scalar");
    }

    /// Degenerate basis vectors (metric 0) give zero products when contracted.
    #[test]
    fn degenerate_basis_gives_zero(i in basis_6d(), other in blade_6d()) {
        // Create a metric where basis i is degenerate
        let metric = |idx: usize| if idx == i { 0 } else { 1 };

        let basis_i = 1 << i;

        // e_i * e_i = 0
        let (sign, _) = basis_product(basis_i, basis_i, metric);
        prop_assert_eq!(sign, 0, "degenerate basis should square to 0");

        // Any blade containing e_i squared should be 0
        // (if the blade contains e_i and we multiply by itself)
        if (other >> i) & 1 == 1 {
            let (sign2, _) = basis_product(other, other, metric);
            prop_assert_eq!(sign2, 0, "blade with degenerate component should square to 0");
        }
    }

    /// ProductTable gives the same results as direct computation.
    #[test]
    fn product_table_matches_direct(a in blade_4d(), b in blade_4d()) {
        let algebra = Algebra::euclidean(4);
        let table = ProductTable::new(&algebra);

        let direct = algebra.basis_product(a, b);
        let lookup = table.geometric(a, b);

        prop_assert_eq!(direct, lookup, "table lookup differs from direct computation");
    }

    /// The reverse of a blade has sign (-1)^(k(k-1)/2) where k is the grade.
    ///
    /// For a blade A of grade k, the reverse Ã satisfies: A * Ã = ||A||² * scalar
    /// The sign pattern is: k=0,1: +, k=2,3: -, k=4,5: +, k=6,7: -, ...
    #[test]
    fn reverse_sign_pattern(blade in blade_6d()) {
        let g = grade(blade);
        let expected_positive = (g * (g.wrapping_sub(1)) / 2) % 2 == 0;

        // For Euclidean metric, blade * blade gives sign based on grade
        let (sign, result) = basis_product(blade, blade, |_| 1);

        // Result should be scalar
        prop_assert_eq!(result, 0, "blade * blade should give scalar");

        // Sign should follow the pattern
        if expected_positive {
            prop_assert!(sign == 1 || sign == 0, "expected positive sign, got {}", sign);
        } else {
            prop_assert!(sign == -1 || sign == 0, "expected negative sign, got {}", sign);
        }
    }
}

/// Tests for different metric signatures.
mod metric_tests {
    use super::*;

    proptest! {
        /// In Minkowski space, the time-like basis squares to -1.
        #[test]
        fn minkowski_signature(spatial in 0usize..3) {
            let algebra = Algebra::minkowski(3);
            let table = ProductTable::new(&algebra);

            // Spatial basis (indices 0, 1, 2) square to +1
            let spatial_blade = 1 << spatial;
            let (sign, result) = table.geometric(spatial_blade, spatial_blade);
            prop_assert_eq!(sign, 1, "spatial basis should square to +1");
            prop_assert_eq!(result, 0, "should give scalar");

            // Time-like basis (index 3) squares to -1
            let time_blade = 1 << 3;
            let (sign, result) = table.geometric(time_blade, time_blade);
            prop_assert_eq!(sign, -1, "time basis should square to -1");
            prop_assert_eq!(result, 0, "should give scalar");
        }

        /// In PGA, the degenerate basis squares to 0.
        #[test]
        fn pga_degenerate_signature(euclidean in 0usize..3) {
            let algebra = Algebra::pga(3);
            let table = ProductTable::new(&algebra);

            // Euclidean basis (indices 0, 1, 2) square to +1
            let euc_blade = 1 << euclidean;
            let (sign, result) = table.geometric(euc_blade, euc_blade);
            prop_assert_eq!(sign, 1, "euclidean basis should square to +1");
            prop_assert_eq!(result, 0, "should give scalar");

            // Degenerate basis (index 3) squares to 0
            let degen_blade = 1 << 3;
            let (sign, _) = table.geometric(degen_blade, degen_blade);
            prop_assert_eq!(sign, 0, "degenerate basis should square to 0");
        }

        /// In CGA, e+ squares to +1 and e- squares to -1.
        #[test]
        fn cga_null_signature(euclidean in 0usize..3) {
            let algebra = Algebra::cga(3);
            let table = ProductTable::new(&algebra);

            // Euclidean basis (indices 0, 1, 2) square to +1
            let euc_blade = 1 << euclidean;
            let (sign, result) = table.geometric(euc_blade, euc_blade);
            prop_assert_eq!(sign, 1, "euclidean basis should square to +1");
            prop_assert_eq!(result, 0);

            // e+ (index 3) squares to +1
            let e_plus = 1 << 3;
            let (sign, result) = table.geometric(e_plus, e_plus);
            prop_assert_eq!(sign, 1, "e+ should square to +1");
            prop_assert_eq!(result, 0);

            // e- (index 4) squares to -1
            let e_minus = 1 << 4;
            let (sign, result) = table.geometric(e_minus, e_minus);
            prop_assert_eq!(sign, -1, "e- should square to -1");
            prop_assert_eq!(result, 0);
        }
    }
}
