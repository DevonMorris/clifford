//! Sign computation for basis blade products.
//!
//! When multiplying basis blades, the sign comes from two sources:
//! 1. **Permutation sign**: The number of swaps needed to reorder into canonical form
//! 2. **Metric contribution**: The product of metric values for shared basis vectors
//!
//! # Example
//!
//! Computing `e₂ * e₁` in Euclidean space:
//! - Result blade: `e₂₁` which must be reordered to `e₁₂`
//! - One swap needed: `e₂₁ → -e₁₂`
//! - No shared basis vectors, so metric contributes `+1`
//! - Final sign: `-1`

use super::blade::MAX_DIM;

/// Computes the sign and result of multiplying basis blades.
///
/// Given two blades represented as bitmasks, computes:
/// - The result blade index (via XOR)
/// - The sign factor from permutation and metric
///
/// # Arguments
///
/// * `a` - First blade bitmask (bit `i` = 1 means `eᵢ` is present)
/// * `b` - Second blade bitmask
/// * `metric` - Function returning the metric value for each basis vector:
///   - `+1` for positive-square basis vectors (Euclidean)
///   - `-1` for negative-square basis vectors (Minkowski)
///   - `0` for null/degenerate basis vectors (PGA, CGA)
///
/// # Returns
///
/// A tuple `(sign, result)` where:
/// - `sign` is `-1`, `0`, or `+1`
/// - `result` is the blade index of the product
///
/// # Example
///
/// ```
/// use clifford_codegen::algebra::basis_product;
///
/// // Euclidean metric: all basis vectors square to +1
/// let euclidean = |_: usize| 1i8;
///
/// // e1 * e2 = e12 with sign +1
/// let (sign, result) = basis_product(0b01, 0b10, euclidean);
/// assert_eq!(sign, 1);
/// assert_eq!(result, 0b11);
///
/// // e2 * e1 = -e12 (anticommutative)
/// let (sign, result) = basis_product(0b10, 0b01, euclidean);
/// assert_eq!(sign, -1);
/// assert_eq!(result, 0b11);
///
/// // e1 * e1 = 1 (vectors square to scalar)
/// let (sign, result) = basis_product(0b01, 0b01, euclidean);
/// assert_eq!(sign, 1);
/// assert_eq!(result, 0);
/// ```
pub fn basis_product<F>(a: usize, b: usize, metric: F) -> (i8, usize)
where
    F: Fn(usize) -> i8,
{
    let result = a ^ b;

    // Count swaps needed to bring b's bits past a's bits
    // For each bit in b, count how many bits of a are to its left
    let mut swaps = 0u32;
    let mut temp_a = a;
    for i in 0..MAX_DIM {
        if (b >> i) & 1 == 1 {
            // Count how many bits of a are above position i
            swaps += (temp_a >> (i + 1)).count_ones();
        }
        // Remove bit i from temp_a if present (bits we've passed)
        temp_a &= !(1 << i);
    }

    // Compute metric contribution from shared bits (where both a and b have a 1)
    let shared = a & b;
    let mut metric_sign = 1i8;
    for i in 0..MAX_DIM {
        if (shared >> i) & 1 == 1 {
            metric_sign *= metric(i);
            if metric_sign == 0 {
                return (0, result);
            }
        }
    }

    // Final sign: (-1)^swaps * metric_sign
    let sign = if swaps.is_multiple_of(2) {
        metric_sign
    } else {
        -metric_sign
    };
    (sign, result)
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Euclidean metric for testing.
    fn euclidean(_i: usize) -> i8 {
        1
    }

    /// PGA metric: last basis (index 0 in our convention) is degenerate.
    fn pga(i: usize) -> i8 {
        if i == 0 { 0 } else { 1 }
    }

    #[test]
    fn scalar_times_anything() {
        // 1 * e1 = e1
        let (sign, result) = basis_product(0, 0b001, euclidean);
        assert_eq!(sign, 1);
        assert_eq!(result, 0b001);

        // 1 * e12 = e12
        let (sign, result) = basis_product(0, 0b011, euclidean);
        assert_eq!(sign, 1);
        assert_eq!(result, 0b011);

        // 1 * e123 = e123
        let (sign, result) = basis_product(0, 0b111, euclidean);
        assert_eq!(sign, 1);
        assert_eq!(result, 0b111);
    }

    #[test]
    fn basis_vectors_square_to_metric() {
        // e1 * e1 = 1 (Euclidean)
        let (sign, result) = basis_product(0b001, 0b001, euclidean);
        assert_eq!(sign, 1);
        assert_eq!(result, 0);

        // e2 * e2 = 1 (Euclidean)
        let (sign, result) = basis_product(0b010, 0b010, euclidean);
        assert_eq!(sign, 1);
        assert_eq!(result, 0);

        // e0 * e0 = 0 (PGA - degenerate)
        let (sign, result) = basis_product(0b001, 0b001, pga);
        assert_eq!(sign, 0);
        assert_eq!(result, 0);
    }

    #[test]
    fn vectors_anticommute() {
        // e1 * e2 = e12
        let (sign1, result1) = basis_product(0b001, 0b010, euclidean);
        assert_eq!(sign1, 1);
        assert_eq!(result1, 0b011);

        // e2 * e1 = -e12
        let (sign2, result2) = basis_product(0b010, 0b001, euclidean);
        assert_eq!(sign2, -1);
        assert_eq!(result2, 0b011);
    }

    #[test]
    fn bivector_products() {
        // e12 * e12 = e1 e2 e1 e2 = -e1 e1 e2 e2 = -1 * 1 * 1 = -1
        let (sign, result) = basis_product(0b011, 0b011, euclidean);
        assert_eq!(sign, -1);
        assert_eq!(result, 0);

        // e12 * e23 = e1 e2 e2 e3 = e1 e3 = e13
        // (after e2*e2 = 1 cancels)
        let (sign, result) = basis_product(0b011, 0b110, euclidean);
        assert_eq!(sign, 1);
        assert_eq!(result, 0b101); // e13

        // e23 * e12 = e2 e3 e1 e2 = -e1 e2 e3 e2 = -e1 e3 = -e13
        let (sign, result) = basis_product(0b110, 0b011, euclidean);
        assert_eq!(sign, -1);
        assert_eq!(result, 0b101);
    }

    #[test]
    fn trivector_products() {
        // e123 * e1 = e2 e3 (since e1 e1 = 1)
        let (sign, result) = basis_product(0b111, 0b001, euclidean);
        assert_eq!(result, 0b110); // e23

        // Sign: e123 * e1 = e1 e2 e3 e1
        // Move e1 past e3: -1, past e2: -1, past e1: contracts
        // = e2 e3
        assert_eq!(sign, 1);

        // e123 * e123 = -1 (for 3D Euclidean)
        let (sign, result) = basis_product(0b111, 0b111, euclidean);
        assert_eq!(sign, -1);
        assert_eq!(result, 0);
    }

    #[test]
    fn pga_degenerate_products() {
        // e0 * e1 = e01 (no contraction, just ordering)
        let (sign, result) = basis_product(0b001, 0b010, pga);
        assert_eq!(sign, 1);
        assert_eq!(result, 0b011);

        // e01 * e01 = e0 e1 e0 e1 = -e0 e0 e1 e1 = 0 * 1 = 0
        let (sign, result) = basis_product(0b011, 0b011, pga);
        assert_eq!(sign, 0);
        assert_eq!(result, 0);
    }

    #[test]
    fn minkowski_metric() {
        // Minkowski: e0 squares to -1, others to +1
        let minkowski = |i: usize| if i == 0 { -1 } else { 1 };

        // e0 * e0 = -1
        let (sign, result) = basis_product(0b001, 0b001, minkowski);
        assert_eq!(sign, -1);
        assert_eq!(result, 0);

        // e1 * e1 = +1
        let (sign, result) = basis_product(0b010, 0b010, minkowski);
        assert_eq!(sign, 1);
        assert_eq!(result, 0);
    }
}
