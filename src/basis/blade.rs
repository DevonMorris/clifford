//! Basis blade multiplication.
//!
//! This module computes the product of two basis blades, including both
//! the resulting blade index and the sign from reordering and metric application.
//!
//! # The Geometric Product of Basis Blades
//!
//! When multiplying two basis blades, the result is another basis blade
//! (possibly with a sign change). For example, in Euclidean 3D:
//!
//! - `e₁ * e₂ = e₁₂` (wedge, no shared vectors)
//! - `e₁ * e₁ = 1` (contraction, `e₁² = +1`)
//! - `e₁₂ * e₂ = e₁` (partial contraction)
//! - `e₁₂ * e₂₁ = -e₁₂ * e₁₂ = -(-1) = 1`
//!
//! # Sign Computation
//!
//! The sign comes from two sources:
//! 1. **Reordering**: Moving basis vectors to canonical order requires swaps.
//!    Each swap of adjacent vectors introduces a factor of `-1`.
//! 2. **Metric**: When a basis vector appears in both blades, it contracts
//!    using the metric: `e_i * e_i = metric(i)` which is `+1`, `-1`, or `0`.
//!
//! # Algorithm
//!
//! The algorithm counts the number of swaps needed to bring all matching
//! vectors together (using bit manipulation), then applies the metric for
//! each pair that cancels.

/// Computes the product of two basis blades.
///
/// Given blade indices `a` and `b`, computes the resulting blade index
/// and sign when they are multiplied under the geometric product.
///
/// # Arguments
///
/// * `a` - Index of the first basis blade (bitmask)
/// * `b` - Index of the second basis blade (bitmask)
/// * `metric` - Function returning the metric coefficient for each basis vector
///
/// # Returns
///
/// A tuple `(sign, result_index)` where:
/// - `sign` is `-1`, `0`, or `+1`
/// - `result_index` is the blade index of the product
///
/// A sign of `0` means the product is zero (from a null basis vector).
///
/// # Mathematical Background
///
/// For basis blades represented as ordered products of basis vectors,
/// the geometric product requires:
/// 1. Concatenating the two sequences of basis vectors
/// 2. Sorting to canonical order (counting swaps for sign)
/// 3. Canceling adjacent identical vectors using the metric
///
/// The XOR of indices gives the result (vectors present in exactly one blade).
/// The sign combines swap parity with metric contributions.
///
/// # Examples
///
/// ```
/// use clifford::basis::basis_product;
///
/// // Euclidean metric: all basis vectors square to +1
/// let euclidean = |_: usize| 1i8;
///
/// // e₁ * e₂ = e₁₂ (no sign change)
/// let (sign, result) = basis_product(0b01, 0b10, euclidean);
/// assert_eq!(sign, 1);
/// assert_eq!(result, 0b11);
///
/// // e₂ * e₁ = -e₁₂ (one swap needed)
/// let (sign, result) = basis_product(0b10, 0b01, euclidean);
/// assert_eq!(sign, -1);
/// assert_eq!(result, 0b11);
///
/// // e₁ * e₁ = 1 (contraction)
/// let (sign, result) = basis_product(0b01, 0b01, euclidean);
/// assert_eq!(sign, 1);
/// assert_eq!(result, 0b00);
/// ```
#[inline]
pub fn basis_product<F>(a: usize, b: usize, metric: F) -> (i8, usize)
where
    F: Fn(usize) -> i8,
{
    // The result index is always XOR (symmetric difference of basis vectors)
    let result_index = a ^ b;

    // Compute the sign from swaps and metric
    let sign = compute_sign(a, b, metric);

    (sign, result_index)
}

/// Computes the sign when multiplying basis blades a and b.
///
/// The sign has two components:
/// 1. Swap sign: (-1)^(number of swaps to reorder)
/// 2. Metric sign: product of metric(i) for all i in a ∩ b
///
/// # Algorithm
///
/// For each bit position `i` in `b` (going from low to high):
/// - Count how many bits in `a` are at positions higher than `i`
///   and will need to swap past `b[i]`
/// - If bit `i` is also set in `a`, apply the metric
///
/// The swap count uses a running mask to track which bits from `a`
/// are "to the right" of the current position in `b`.
#[inline]
fn compute_sign<F>(a: usize, b: usize, metric: F) -> i8
where
    F: Fn(usize) -> i8,
{
    // Bits that appear in both (will contract)
    let common = a & b;

    // Check for null vectors first (zero metric means zero product)
    let mut temp_common = common;
    while temp_common != 0 {
        let i = temp_common.trailing_zeros() as usize;
        if metric(i) == 0 {
            return 0;
        }
        temp_common &= temp_common - 1; // Clear lowest set bit
    }

    // Count swaps needed to bring matching vectors together
    let swaps = count_swaps(a, b);

    // Compute metric sign from contractions
    let metric_sign = compute_metric_sign(common, &metric);

    // Combine swap parity with metric sign
    if swaps.is_multiple_of(2) {
        metric_sign
    } else {
        -metric_sign
    }
}

/// Counts the number of swaps needed to reorder the product of blades a and b.
///
/// When we multiply blade `a` by blade `b`, we conceptually concatenate their
/// basis vectors and then sort to canonical order. Each swap introduces a `-1`.
///
/// # Algorithm
///
/// For each basis vector in `b` at position `j`, count how many basis vectors
/// in `a` are at positions greater than `j`. These must swap past `b[j]`.
#[inline]
fn count_swaps(a: usize, b: usize) -> usize {
    let mut swaps = 0;
    let mut b_remaining = b;

    while b_remaining != 0 {
        // Get the lowest set bit position in b
        let j = b_remaining.trailing_zeros() as usize;

        // Count bits in `a` that are at positions > j
        // These are the vectors from `a` that must swap past e_j
        let mask_above_j = !((1usize << (j + 1)) - 1); // bits above position j
        swaps += (a & mask_above_j).count_ones() as usize;

        // Clear the lowest set bit
        b_remaining &= b_remaining - 1;
    }

    swaps
}

/// Computes the sign contribution from the metric.
///
/// For each basis vector that appears in both blades (the common bits),
/// we get a factor of `metric(i)` when they contract.
#[inline]
fn compute_metric_sign<F>(common: usize, metric: F) -> i8
where
    F: Fn(usize) -> i8,
{
    let mut sign: i8 = 1;
    let mut remaining = common;

    while remaining != 0 {
        let i = remaining.trailing_zeros() as usize;
        sign *= metric(i);
        remaining &= remaining - 1; // Clear lowest set bit
    }

    sign
}

/// Checks if two basis blades anticommute.
///
/// Two basis blades anticommute if `a * b = -b * a`, which happens
/// when their grades are both odd.
///
/// # Arguments
///
/// * `a` - Index of the first basis blade
/// * `b` - Index of the second basis blade
///
/// # Returns
///
/// `true` if the blades anticommute, `false` if they commute.
///
/// # Example
///
/// ```
/// use clifford::basis::anticommutes;
///
/// // Two vectors anticommute
/// assert!(anticommutes(0b001, 0b010));
///
/// // Vector and bivector commute (in 3D Euclidean)
/// assert!(!anticommutes(0b001, 0b110));
///
/// // Two bivectors may commute or anticommute depending on overlap
/// // e₁₂ and e₂₃ share one vector, so they anticommute
/// // (But this function checks grade parity, not specific overlap)
/// ```
#[inline]
pub const fn anticommutes(a: usize, b: usize) -> bool {
    // Blades anticommute if the product of their grades is odd
    // This is a simplification; the actual anticommutation depends on the
    // overlap structure, but for grade-based analysis this is useful
    let grade_a = a.count_ones();
    let grade_b = b.count_ones();
    (grade_a * grade_b) % 2 == 1
}

#[cfg(test)]
mod tests {
    use super::*;
    use proptest::prelude::*;

    fn euclidean_metric(_i: usize) -> i8 {
        1
    }

    proptest! {
        /// The result index is always XOR of the inputs.
        #[test]
        fn result_is_xor(a in 0usize..64, b in 0usize..64) {
            let (_, result) = basis_product(a, b, euclidean_metric);
            prop_assert_eq!(result, a ^ b);
        }

        /// Scalar (index 0) is the identity: 1 * a = a with sign +1.
        #[test]
        fn scalar_identity(a in 0usize..64) {
            let (sign, result) = basis_product(0, a, euclidean_metric);
            prop_assert_eq!(sign, 1);
            prop_assert_eq!(result, a);

            let (sign, result) = basis_product(a, 0, euclidean_metric);
            prop_assert_eq!(sign, 1);
            prop_assert_eq!(result, a);
        }

        /// Associativity of sign: sign(a, bc) * sign(b, c) = sign(ab, c) * sign(a, b)
        /// (for non-zero products in Euclidean signature)
        #[test]
        fn sign_associativity(a in 0usize..16, b in 0usize..16, c in 0usize..16) {
            let (s_bc, bc) = basis_product(b, c, euclidean_metric);
            let (s_a_bc, _abc1) = basis_product(a, bc, euclidean_metric);

            let (s_ab, ab) = basis_product(a, b, euclidean_metric);
            let (s_ab_c, _abc2) = basis_product(ab, c, euclidean_metric);

            // Both paths should give same total sign
            prop_assert_eq!(
                (s_a_bc as i32) * (s_bc as i32),
                (s_ab_c as i32) * (s_ab as i32)
            );
        }
    }

    #[test]
    fn vector_products() {
        // e₁ * e₁ = 1
        let (sign, result) = basis_product(0b001, 0b001, euclidean_metric);
        assert_eq!(sign, 1);
        assert_eq!(result, 0);

        // e₁ * e₂ = e₁₂
        let (sign, result) = basis_product(0b001, 0b010, euclidean_metric);
        assert_eq!(sign, 1);
        assert_eq!(result, 0b011);

        // e₂ * e₁ = -e₁₂
        let (sign, result) = basis_product(0b010, 0b001, euclidean_metric);
        assert_eq!(sign, -1);
        assert_eq!(result, 0b011);
    }

    #[test]
    fn bivector_products() {
        // e₁₂ * e₁₂ = e₁e₂e₁e₂ = -e₁e₁e₂e₂ = -1
        let (sign, result) = basis_product(0b011, 0b011, euclidean_metric);
        assert_eq!(sign, -1);
        assert_eq!(result, 0);

        // e₁₂ * e₂ = e₁e₂e₂ = e₁
        let (sign, result) = basis_product(0b011, 0b010, euclidean_metric);
        assert_eq!(sign, 1);
        assert_eq!(result, 0b001);

        // e₂ * e₁₂ = e₂e₁e₂ = -e₁e₂e₂ = -e₁
        let (sign, result) = basis_product(0b010, 0b011, euclidean_metric);
        assert_eq!(sign, -1);
        assert_eq!(result, 0b001);
    }

    #[test]
    fn minkowski_metric() {
        // In Minkowski (1,3): e₀² = +1, e₁² = e₂² = e₃² = -1
        let minkowski = |i: usize| if i == 0 { 1 } else { -1 };

        // e₁ * e₁ = -1 (space-like)
        let (sign, result) = basis_product(0b0010, 0b0010, minkowski);
        assert_eq!(sign, -1);
        assert_eq!(result, 0);

        // e₀ * e₀ = +1 (time-like)
        let (sign, result) = basis_product(0b0001, 0b0001, minkowski);
        assert_eq!(sign, 1);
        assert_eq!(result, 0);
    }

    #[test]
    fn null_metric() {
        // PGA-like: e₀² = 0
        let pga_metric = |i: usize| if i == 0 { 0 } else { 1 };

        // e₀ * e₀ = 0
        let (sign, result) = basis_product(0b001, 0b001, pga_metric);
        assert_eq!(sign, 0);
        assert_eq!(result, 0);

        // e₀₁ * e₀ = 0 (e₀ appears in both, so contracts to 0)
        let (sign, _result) = basis_product(0b011, 0b001, pga_metric);
        assert_eq!(sign, 0);
    }

    #[test]
    fn swap_count_examples() {
        // e₁ * e₂: no swaps needed (1 comes before 2)
        assert_eq!(count_swaps(0b001, 0b010), 0);

        // e₂ * e₁: one swap needed (2 must pass 1)
        assert_eq!(count_swaps(0b010, 0b001), 1);

        // e₃ * e₁: one swap (3 passes 1)
        assert_eq!(count_swaps(0b100, 0b001), 1);

        // e₃ * e₁₂: two swaps (3 passes both 1 and 2)
        assert_eq!(count_swaps(0b100, 0b011), 2);

        // e₂₃ * e₁: one swap (neither 2 nor 3 needs to pass 1)
        // Wait: 2 is at position 1, 3 is at position 2, 1 is at position 0
        // For each bit in b (just bit 0), count bits in a above position 0
        // a = 0b110, above position 0 is bits 1 and 2, both set = 2 swaps
        assert_eq!(count_swaps(0b110, 0b001), 2);
    }
}
