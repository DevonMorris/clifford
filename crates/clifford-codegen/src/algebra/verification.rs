//! Verification tests against the clifford library.
//!
//! These tests ensure that our blade algebra engine produces results
//! identical to the reference `clifford::Multivector` implementation.

#![cfg(test)]

use clifford::basis::Blade as CliffordBlade;
use clifford::prelude::Multivector;
use clifford::signature::{Euclidean2, Euclidean3};

use super::{Algebra, ProductTable};

/// Verifies that our product table matches clifford::Multivector for Euclidean 2D.
#[test]
fn matches_clifford_euclidean_2d() {
    let algebra = Algebra::euclidean(2);
    let table = ProductTable::new(&algebra);

    // 2D has 4 blades: 1, e1, e2, e12
    for a in 0..4 {
        for b in 0..4 {
            let (gen_sign, gen_result) = table.geometric(a, b);

            // Create Multivectors with single unit coefficients
            let mut mv_a = Multivector::<f64, Euclidean2>::zero();
            mv_a.set(CliffordBlade::from_index(a), 1.0);

            let mut mv_b = Multivector::<f64, Euclidean2>::zero();
            mv_b.set(CliffordBlade::from_index(b), 1.0);

            let mv_result = mv_a * mv_b;

            // Find the non-zero coefficient in the result
            let mut found_blade = None;
            let mut found_coeff = 0.0;
            for i in 0..4 {
                let coeff = mv_result.get(CliffordBlade::from_index(i));
                if coeff.abs() > 1e-10 {
                    assert!(
                        found_blade.is_none(),
                        "Multiple non-zero blades in result for {} * {}",
                        a,
                        b
                    );
                    found_blade = Some(i);
                    found_coeff = coeff;
                }
            }

            if gen_sign == 0 {
                assert!(
                    found_blade.is_none(),
                    "Expected zero product for {} * {}, got blade {} with coeff {}",
                    a,
                    b,
                    found_blade.unwrap_or(0),
                    found_coeff
                );
            } else {
                assert!(
                    found_blade.is_some(),
                    "Expected non-zero product for {} * {}, got zero",
                    a,
                    b
                );
                assert_eq!(
                    found_blade.unwrap(),
                    gen_result,
                    "Wrong result blade for {} * {}: expected {}, got {}",
                    a,
                    b,
                    gen_result,
                    found_blade.unwrap()
                );
                assert_eq!(
                    found_coeff as i8, gen_sign,
                    "Wrong sign for {} * {}: expected {}, got {}",
                    a, b, gen_sign, found_coeff
                );
            }
        }
    }
}

/// Verifies that our product table matches clifford::Multivector for Euclidean 3D.
#[test]
fn matches_clifford_euclidean_3d() {
    let algebra = Algebra::euclidean(3);
    let table = ProductTable::new(&algebra);

    // 3D has 8 blades
    for a in 0..8 {
        for b in 0..8 {
            let (gen_sign, gen_result) = table.geometric(a, b);

            let mut mv_a = Multivector::<f64, Euclidean3>::zero();
            mv_a.set(CliffordBlade::from_index(a), 1.0);

            let mut mv_b = Multivector::<f64, Euclidean3>::zero();
            mv_b.set(CliffordBlade::from_index(b), 1.0);

            let mv_result = mv_a * mv_b;

            // Find the non-zero coefficient
            let mut found_blade = None;
            let mut found_coeff = 0.0;
            for i in 0..8 {
                let coeff = mv_result.get(CliffordBlade::from_index(i));
                if coeff.abs() > 1e-10 {
                    assert!(
                        found_blade.is_none(),
                        "Multiple non-zero blades in result for {} * {}",
                        a,
                        b
                    );
                    found_blade = Some(i);
                    found_coeff = coeff;
                }
            }

            if gen_sign == 0 {
                assert!(
                    found_blade.is_none(),
                    "Expected zero product for {} * {}",
                    a,
                    b
                );
            } else {
                assert!(
                    found_blade.is_some(),
                    "Expected non-zero product for {} * {}",
                    a,
                    b
                );
                assert_eq!(
                    found_blade.unwrap(),
                    gen_result,
                    "Wrong result blade for {} * {}",
                    a,
                    b
                );
                assert_eq!(
                    found_coeff as i8, gen_sign,
                    "Wrong sign for {} * {}: expected {}, got {}",
                    a, b, gen_sign, found_coeff
                );
            }
        }
    }
}

/// Tests specific known products to verify signs are correct.
#[test]
fn known_products_euclidean_3d() {
    let algebra = Algebra::euclidean(3);
    let table = ProductTable::new(&algebra);

    // Define blade indices
    let e1 = 1;
    let e2 = 2;
    let e12 = 3;
    let e3 = 4;
    let e13 = 5;
    let e23 = 6;
    let e123 = 7;

    // Vectors square to +1
    assert_eq!(table.geometric(e1, e1), (1, 0));
    assert_eq!(table.geometric(e2, e2), (1, 0));
    assert_eq!(table.geometric(e3, e3), (1, 0));

    // Vectors anticommute
    assert_eq!(table.geometric(e1, e2), (1, e12));
    assert_eq!(table.geometric(e2, e1), (-1, e12));
    assert_eq!(table.geometric(e1, e3), (1, e13));
    assert_eq!(table.geometric(e3, e1), (-1, e13));
    assert_eq!(table.geometric(e2, e3), (1, e23));
    assert_eq!(table.geometric(e3, e2), (-1, e23));

    // Bivectors square to -1
    assert_eq!(table.geometric(e12, e12), (-1, 0));
    assert_eq!(table.geometric(e13, e13), (-1, 0));
    assert_eq!(table.geometric(e23, e23), (-1, 0));

    // Pseudoscalar squares to -1
    assert_eq!(table.geometric(e123, e123), (-1, 0));

    // Bivector * vector products
    // e12 * e3 = e123
    assert_eq!(table.geometric(e12, e3), (1, e123));

    // e12 * e1 = e1 e2 e1 = -e1 e1 e2 = -1 * e2 = -e2
    let (sign, result) = table.geometric(e12, e1);
    assert_eq!(result, e2);
    assert_eq!(sign, -1, "e12 * e1 should give -e2");
}

/// Verifies associativity using clifford::Multivector.
#[test]
fn associativity_matches_multivector() {
    let algebra = Algebra::euclidean(3);
    let table = ProductTable::new(&algebra);

    // Test several triple products
    let test_cases = [
        (1, 2, 4), // e1 * e2 * e3
        (3, 5, 6), // e12 * e13 * e23
        (1, 3, 7), // e1 * e12 * e123
        (7, 7, 7), // e123 * e123 * e123
        (2, 3, 4), // e2 * e12 * e3
        (0, 5, 3), // 1 * e13 * e12
        (6, 1, 2), // e23 * e1 * e2
        (4, 6, 7), // e3 * e23 * e123
    ];

    for (a, b, c) in test_cases {
        // Compute (a * b) * c using our engine
        let (sign_ab, ab) = table.geometric(a, b);
        let (sign_abc_left, result_left) = if sign_ab != 0 {
            let (s, r) = table.geometric(ab, c);
            (sign_ab as i32 * s as i32, r)
        } else {
            (0, 0)
        };

        // Compute a * (b * c) using our engine
        let (sign_bc, bc) = table.geometric(b, c);
        let (sign_abc_right, result_right) = if sign_bc != 0 {
            let (s, r) = table.geometric(a, bc);
            (sign_bc as i32 * s as i32, r)
        } else {
            (0, 0)
        };

        // Verify they match
        assert_eq!(
            result_left, result_right,
            "Associativity failed for ({}, {}, {}): left result {}, right result {}",
            a, b, c, result_left, result_right
        );
        assert_eq!(
            sign_abc_left, sign_abc_right,
            "Associativity sign failed for ({}, {}, {}): left sign {}, right sign {}",
            a, b, c, sign_abc_left, sign_abc_right
        );

        // Also verify against Multivector
        let mut mv_a = Multivector::<f64, Euclidean3>::zero();
        mv_a.set(CliffordBlade::from_index(a), 1.0);
        let mut mv_b = Multivector::<f64, Euclidean3>::zero();
        mv_b.set(CliffordBlade::from_index(b), 1.0);
        let mut mv_c = Multivector::<f64, Euclidean3>::zero();
        mv_c.set(CliffordBlade::from_index(c), 1.0);

        let mv_left = (mv_a * mv_b) * mv_c;
        let mv_right = mv_a * (mv_b * mv_c);

        // They should be equal
        for i in 0..8 {
            let blade = CliffordBlade::from_index(i);
            let diff = (mv_left.get(blade) - mv_right.get(blade)).abs();
            assert!(
                diff < 1e-10,
                "Multivector associativity failed for ({}, {}, {}) at blade {}: {} vs {}",
                a,
                b,
                c,
                i,
                mv_left.get(blade),
                mv_right.get(blade)
            );
        }
    }
}

/// Comprehensive test comparing all products against Multivector.
#[test]
fn all_products_match_multivector_3d() {
    let algebra = Algebra::euclidean(3);
    let table = ProductTable::new(&algebra);

    for a in 0..8 {
        for b in 0..8 {
            let (our_sign, our_result) = table.geometric(a, b);

            // Compute with Multivector
            let mut mv_a = Multivector::<f64, Euclidean3>::zero();
            mv_a.set(CliffordBlade::from_index(a), 1.0);

            let mut mv_b = Multivector::<f64, Euclidean3>::zero();
            mv_b.set(CliffordBlade::from_index(b), 1.0);

            let mv_result = mv_a * mv_b;

            // Extract result from Multivector
            let mut mv_sign = 0i8;
            let mut mv_blade = 0;
            for i in 0..8 {
                let coeff = mv_result.get(CliffordBlade::from_index(i));
                if coeff.abs() > 1e-10 {
                    mv_sign = coeff.signum() as i8;
                    mv_blade = i;
                }
            }

            // Compare
            assert_eq!(
                our_sign, mv_sign,
                "Sign mismatch for {} * {}: ours={}, mv={}",
                a, b, our_sign, mv_sign
            );
            if our_sign != 0 {
                assert_eq!(
                    our_result, mv_blade,
                    "Result mismatch for {} * {}: ours={}, mv={}",
                    a, b, our_result, mv_blade
                );
            }
        }
    }
}
