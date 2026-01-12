//! Precomputed product tables for efficient lookup.
//!
//! Product tables precompute all basis blade products for an algebra,
//! enabling O(1) lookup during code generation. This is essential for
//! generating optimized product implementations.

use super::signature::Algebra;

/// Precomputed product table for a geometric algebra.
///
/// Stores the sign and result blade for all pairs of basis blades.
/// This enables O(1) lookup of any product during code generation.
///
/// # Example
///
/// ```
/// use clifford_codegen::algebra::{Algebra, ProductTable};
///
/// let algebra = Algebra::euclidean(3);
/// let table = ProductTable::new(&algebra);
///
/// // e1 * e2 = e12 with sign +1
/// let (sign, result) = table.geometric(1, 2);
/// assert_eq!(sign, 1);
/// assert_eq!(result, 3);
///
/// // e2 * e1 = -e12 (anticommutative)
/// let (sign, result) = table.geometric(2, 1);
/// assert_eq!(sign, -1);
/// assert_eq!(result, 3);
/// ```
#[derive(Clone, Debug)]
pub struct ProductTable {
    /// The algebra dimension.
    dim: usize,
    /// Signature (p, q, r): p positive, q negative, r degenerate.
    signature: (usize, usize, usize),
    /// Signs for geometric product: signs[a * n + b] = sign of e_a * e_b.
    signs: Vec<i8>,
    /// Result indices: results[a * n + b] = index of e_a * e_b.
    results: Vec<usize>,
}

impl ProductTable {
    /// Builds a product table for the given algebra.
    ///
    /// Precomputes all `n²` products where `n = 2^dim`.
    ///
    /// # Example
    ///
    /// ```
    /// use clifford_codegen::algebra::{Algebra, ProductTable};
    ///
    /// let algebra = Algebra::euclidean(3);
    /// let table = ProductTable::new(&algebra);
    ///
    /// // Table has 8*8 = 64 entries for 3D
    /// assert_eq!(table.dim(), 3);
    /// ```
    pub fn new(algebra: &Algebra) -> Self {
        let n = algebra.num_blades();
        let mut signs = vec![0i8; n * n];
        let mut results = vec![0usize; n * n];

        for a in 0..n {
            for b in 0..n {
                let (sign, result) = algebra.basis_product(a, b);
                signs[a * n + b] = sign;
                results[a * n + b] = result;
            }
        }

        Self {
            dim: algebra.dim(),
            signature: algebra.signature(),
            signs,
            results,
        }
    }

    /// Returns the algebra dimension.
    #[inline]
    pub fn dim(&self) -> usize {
        self.dim
    }

    /// Returns the number of blades (2^dim).
    #[inline]
    pub fn num_blades(&self) -> usize {
        1 << self.dim
    }

    /// Looks up the geometric product of two basis blades.
    ///
    /// # Returns
    ///
    /// A tuple `(sign, result)` where:
    /// - `sign` is the sign factor (-1, 0, or +1)
    /// - `result` is the blade index of the product
    ///
    /// # Panics
    ///
    /// Panics if `a` or `b` is out of bounds.
    ///
    /// # Example
    ///
    /// ```
    /// use clifford_codegen::algebra::{Algebra, ProductTable};
    ///
    /// let algebra = Algebra::euclidean(3);
    /// let table = ProductTable::new(&algebra);
    ///
    /// // e12 * e12 = -1 (bivector squares to -1)
    /// let (sign, result) = table.geometric(3, 3);
    /// assert_eq!(sign, -1);
    /// assert_eq!(result, 0);
    /// ```
    #[inline]
    pub fn geometric(&self, a: usize, b: usize) -> (i8, usize) {
        let n = self.num_blades();
        let idx = a * n + b;
        (self.signs[idx], self.results[idx])
    }

    /// Returns the sign of the product of two basis blades.
    #[inline]
    pub fn sign(&self, a: usize, b: usize) -> i8 {
        self.signs[a * self.num_blades() + b]
    }

    /// Returns the result blade of the product of two basis blades.
    #[inline]
    pub fn result(&self, a: usize, b: usize) -> usize {
        self.results[a * self.num_blades() + b]
    }

    /// Returns the metric value for a basis vector.
    ///
    /// - Returns `+1` for positive basis vectors (indices `0..p`)
    /// - Returns `-1` for negative basis vectors (indices `p..p+q`)
    /// - Returns `0` for degenerate basis vectors (indices `p+q..dim`)
    #[inline]
    #[allow(dead_code)]
    fn metric(&self, basis: usize) -> i8 {
        let (p, q, _r) = self.signature;
        if basis < p {
            1
        } else if basis < p + q {
            -1
        } else {
            0
        }
    }

    /// Returns the anti-metric value for a basis vector.
    ///
    /// The anti-metric "swaps" the degenerate subspace. In PGA (3,0,1):
    /// - Original metric: e₄ squares to 0
    /// - Anti-metric: e₄ squares to 1 (so products involving e₄² don't vanish)
    ///
    /// This is essential for correct antiproduct computation in degenerate algebras.
    #[inline]
    fn anti_metric(&self, basis: usize) -> i8 {
        let (p, q, _r) = self.signature;
        if basis < p {
            // Positive in original metric → positive in anti-metric
            1
        } else if basis < p + q {
            // Negative in original metric → negative in anti-metric
            -1
        } else {
            // Degenerate in original metric → positive in anti-metric
            // This is the key difference: e₀² = 0 becomes ē₀² = 1
            1
        }
    }

    /// Computes the metric contribution for the anti-metric product.
    ///
    /// For an overlap (shared basis vectors), returns the product of
    /// anti-metric values for each basis in the overlap.
    ///
    /// # Returns
    ///
    /// - The product of anti-metric values, or
    /// - 0 if the overlap represents a degenerate combination in anti-space
    fn anti_metric_contribution(&self, overlap: usize) -> i8 {
        let (_p, _q, r) = self.signature;

        // In the anti-metric, the degenerate subspace is the COMPLEMENT of
        // the original degenerate subspace. For PGA where e₄ is degenerate:
        // - Original: e₄² = 0
        // - Anti: (e₁₂₃)² maps to degenerate behavior
        //
        // The check: if ALL non-degenerate basis vectors are in the overlap
        // AND no degenerate basis vectors are in the overlap, the contribution
        // is 0 (this represents squaring the "anti-degenerate" subspace).
        //
        // More specifically: overlap is degenerate in anti-space if it consists
        // only of the non-degenerate part of the algebra.

        // Mask of non-degenerate basis vectors
        let non_degenerate_mask = (1usize << (self.dim - r)) - 1;

        // Check if overlap contains ONLY non-degenerate bases AND all of them
        if r > 0 && overlap == non_degenerate_mask {
            // This is the anti-degenerate case: e₁₂₃...eₚ₊ᵧ squared
            return 0;
        }

        // Otherwise, compute product of anti-metric values
        let mut contribution: i8 = 1;
        for basis in 0..self.dim {
            if (overlap >> basis) & 1 == 1 {
                let m = self.anti_metric(basis);
                // In anti-metric, degenerate bases become +1, so we never get 0 here
                contribution *= m;
            }
        }

        contribution
    }

    /// Computes the geometric product using the anti-metric.
    ///
    /// This is the dual of the geometric product, used for computing
    /// the antiproduct in degenerate algebras like PGA.
    ///
    /// # Returns
    ///
    /// A tuple `(sign, result)` where:
    /// - `sign` is the sign factor (-1, 0, or +1)
    /// - `result` is the blade index of the product
    pub fn geometric_anti(&self, a: usize, b: usize) -> (i8, usize) {
        let result_blade = a ^ b;
        let overlap = a & b;

        if overlap == 0 {
            // No metric contribution, just reordering sign
            return (self.exterior_sign(a, b), result_blade);
        }

        // Use anti-metric for overlap contribution
        let metric_sign = self.anti_metric_contribution(overlap);
        if metric_sign == 0 {
            return (0, 0);
        }

        // Compute reordering sign (same as regular geometric product)
        let reorder_sign = self.compute_reorder_sign(a, b);

        (reorder_sign * metric_sign, result_blade)
    }

    /// Computes the sign from reordering basis vectors in a product.
    ///
    /// When computing e_A * e_B, we need to count how many transpositions
    /// are needed to bring the basis vectors into canonical order.
    fn compute_reorder_sign(&self, a: usize, b: usize) -> i8 {
        // For each bit position, count how many bits in 'a' are above it
        // (need to swap past them when merging)
        let mut transpositions = 0;
        let mut b_remaining = b;

        while b_remaining != 0 {
            let lowest_b = b_remaining & b_remaining.wrapping_neg();
            let b_pos = lowest_b.trailing_zeros();

            // Count bits in 'a' that are above this position
            let a_above = a >> (b_pos + 1);
            transpositions += a_above.count_ones();

            b_remaining &= !lowest_b;
        }

        if transpositions % 2 == 0 { 1 } else { -1 }
    }

    /// Finds all products that contribute to a given result blade.
    ///
    /// Given sets of input blades A and B, returns all (a, b) pairs
    /// such that `a * b = result` (with non-zero sign).
    ///
    /// # Arguments
    ///
    /// * `a_blades` - Set of blade indices from the first operand
    /// * `b_blades` - Set of blade indices from the second operand
    /// * `result_blade` - The target result blade index
    ///
    /// # Returns
    ///
    /// Vector of `(sign, a_blade, b_blade)` tuples contributing to the result.
    ///
    /// # Example
    ///
    /// ```
    /// use clifford_codegen::algebra::{Algebra, ProductTable};
    ///
    /// let algebra = Algebra::euclidean(3);
    /// let table = ProductTable::new(&algebra);
    ///
    /// // Which vector * vector products contribute to e12?
    /// let vectors = vec![1, 2, 4]; // e1, e2, e3
    /// let contributions = table.product_contributions(&vectors, &vectors, 3);
    ///
    /// // e1 * e2 = +e12, e2 * e1 = -e12
    /// assert_eq!(contributions.len(), 2);
    /// assert!(contributions.contains(&(1, 1, 2)));  // e1 * e2 = +e12
    /// assert!(contributions.contains(&(-1, 2, 1))); // e2 * e1 = -e12
    /// ```
    pub fn product_contributions(
        &self,
        a_blades: &[usize],
        b_blades: &[usize],
        result_blade: usize,
    ) -> Vec<(i8, usize, usize)> {
        let mut contributions = Vec::new();

        for &a in a_blades {
            for &b in b_blades {
                let (sign, result) = self.geometric(a, b);
                if result == result_blade && sign != 0 {
                    contributions.push((sign, a, b));
                }
            }
        }

        contributions
    }

    /// Checks if a product contributes to a given grade.
    ///
    /// Returns true if any `a * b` product (for a in a_blades, b in b_blades)
    /// produces a blade of the target grade.
    pub fn has_contributions_to_grade(
        &self,
        a_blades: &[usize],
        b_blades: &[usize],
        target_grade: usize,
    ) -> bool {
        for &a in a_blades {
            for &b in b_blades {
                let (sign, result) = self.geometric(a, b);
                if sign != 0 && result.count_ones() as usize == target_grade {
                    return true;
                }
            }
        }
        false
    }

    /// Returns all result blades from products of a_blades × b_blades.
    ///
    /// # Returns
    ///
    /// Vector of (blade_index, contributions) pairs, sorted by blade index.
    /// Each contribution is (sign, a_blade, b_blade).
    #[allow(clippy::type_complexity)]
    pub fn all_products(
        &self,
        a_blades: &[usize],
        b_blades: &[usize],
    ) -> Vec<(usize, Vec<(i8, usize, usize)>)> {
        use std::collections::BTreeMap;

        let mut result_map: BTreeMap<usize, Vec<(i8, usize, usize)>> = BTreeMap::new();

        for &a in a_blades {
            for &b in b_blades {
                let (sign, result) = self.geometric(a, b);
                if sign != 0 {
                    result_map.entry(result).or_default().push((sign, a, b));
                }
            }
        }

        result_map.into_iter().collect()
    }

    /// Computes the right complement of a blade.
    ///
    /// The right complement maps a grade-k blade to a grade-(n-k) blade where n = dim.
    /// For blade index `i`, the complement index is `(2^dim - 1) XOR i`.
    ///
    /// The right complement is defined by: `u ∧ ū = I` (pseudoscalar)
    /// where ∧ is the exterior (wedge) product.
    ///
    /// # Returns
    ///
    /// A tuple `(sign, complement_index)` where sign is ±1 based on the
    /// permutation required to bring the blade and its complement into canonical
    /// order to form the pseudoscalar via the exterior product.
    ///
    /// # Example
    ///
    /// ```
    /// use clifford_codegen::algebra::{Algebra, ProductTable};
    ///
    /// let algebra = Algebra::euclidean(3);
    /// let table = ProductTable::new(&algebra);
    ///
    /// // In 3D: complement(e1) = e23 (with some sign)
    /// let (sign, result) = table.complement(1);
    /// assert_eq!(result, 6); // e23 = 0b110
    /// ```
    pub fn complement(&self, blade: usize) -> (i8, usize) {
        let pseudoscalar = self.num_blades() - 1; // All bits set
        let complement_blade = pseudoscalar ^ blade;

        // Sign is determined by: blade ∧ complement_blade = ±pseudoscalar
        // Since blade and complement_blade don't share any basis vectors,
        // the exterior product equals the geometric product for this pair.
        // We compute the sign by counting transpositions needed to put
        // the combined basis vectors into canonical order.
        let sign = self.exterior_sign(blade, complement_blade);

        (sign, complement_blade)
    }

    /// Computes the sign of the exterior product of two non-overlapping blades.
    ///
    /// For blades that don't share basis vectors, this is the sign from
    /// reordering the basis vectors into canonical order.
    fn exterior_sign(&self, a: usize, b: usize) -> i8 {
        debug_assert_eq!(a & b, 0, "blades must not overlap for exterior sign");

        // Count transpositions: for each bit in b, count how many bits in a
        // are to the right of it (i.e., have higher index but appear earlier)
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

    /// Computes the geometric antiproduct of two blades.
    ///
    /// The antiproduct is defined as:
    /// `a ⊛ b = ∁(∁a *_anti ∁b)`
    ///
    /// where `*_anti` is the geometric product using the **anti-metric** and
    /// `∁` is the complement.
    ///
    /// In degenerate algebras like PGA, the standard formula `∁(∁a * ∁b)` using
    /// the original metric fails because products involving the degenerate basis
    /// vanish. The anti-metric fixes this by "swapping" which directions are
    /// degenerate:
    /// - Original PGA metric: e₀² = 0
    /// - Anti-metric: ē₀² = 1 (the complement of e₀ becomes degenerate instead)
    ///
    /// # Returns
    ///
    /// A tuple `(sign, result)` where:
    /// - `sign` is the sign factor (-1, 0, or +1)
    /// - `result` is the blade index of the antiproduct
    ///
    /// Note: The sign may be 0 if the product vanishes.
    pub fn antiproduct(&self, a: usize, b: usize) -> (i8, usize) {
        // Get complements
        let (sign_ca, comp_a) = self.complement(a);
        let (sign_cb, comp_b) = self.complement(b);

        // Geometric product of complements using ANTI-METRIC
        let (sign_prod, prod) = self.geometric_anti(comp_a, comp_b);
        if sign_prod == 0 {
            return (0, 0);
        }

        // Complement of the product
        let (sign_result, result) = self.complement(prod);

        // Total sign
        let total_sign = sign_ca * sign_cb * sign_prod * sign_result;

        (total_sign, result)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn euclidean_3d_basic() {
        let algebra = Algebra::euclidean(3);
        let table = ProductTable::new(&algebra);

        assert_eq!(table.dim(), 3);
        assert_eq!(table.num_blades(), 8);

        // Scalar * anything = anything
        for b in 0..8 {
            let (sign, result) = table.geometric(0, b);
            assert_eq!(sign, 1);
            assert_eq!(result, b);
        }

        // Vectors square to +1
        assert_eq!(table.geometric(1, 1), (1, 0));
        assert_eq!(table.geometric(2, 2), (1, 0));
        assert_eq!(table.geometric(4, 4), (1, 0));

        // Vectors anticommute
        assert_eq!(table.geometric(1, 2), (1, 3));
        assert_eq!(table.geometric(2, 1), (-1, 3));
    }

    #[test]
    fn euclidean_3d_bivectors() {
        let algebra = Algebra::euclidean(3);
        let table = ProductTable::new(&algebra);

        // Bivectors square to -1
        assert_eq!(table.geometric(3, 3), (-1, 0)); // e12 * e12 = -1
        assert_eq!(table.geometric(5, 5), (-1, 0)); // e13 * e13 = -1
        assert_eq!(table.geometric(6, 6), (-1, 0)); // e23 * e23 = -1

        // e12 * e23 = e1 * e2 * e2 * e3 = e1 * e3 = e13
        let (sign, result) = table.geometric(3, 6);
        assert_eq!(result, 5); // e13
        assert_eq!(sign, 1);
    }

    #[test]
    fn euclidean_3d_pseudoscalar() {
        let algebra = Algebra::euclidean(3);
        let table = ProductTable::new(&algebra);

        // e123 * e123 = -1
        assert_eq!(table.geometric(7, 7), (-1, 0));

        // e1 * e23 = e123
        assert_eq!(table.geometric(1, 6), (1, 7));

        // e23 * e1 = e2 e3 e1 = -e2 e1 e3 = e1 e2 e3 = e123
        assert_eq!(table.geometric(6, 1), (1, 7));
    }

    #[test]
    fn pga_degenerate() {
        let algebra = Algebra::pga(3);
        let table = ProductTable::new(&algebra);

        // e4 is degenerate (index 8)
        assert_eq!(table.geometric(8, 8), (0, 0));

        // Products involving e4 * e4 are zero
        // e14 * e14 = e1 e4 e1 e4 = -e1 e1 e4 e4 = -1 * 0 = 0
        assert_eq!(table.geometric(9, 9), (0, 0));
    }

    #[test]
    fn product_contributions_vectors() {
        let algebra = Algebra::euclidean(3);
        let table = ProductTable::new(&algebra);

        let vectors = vec![1, 2, 4];

        // Contributions to e12 from vectors
        let contrib = table.product_contributions(&vectors, &vectors, 3);
        assert_eq!(contrib.len(), 2);
        assert!(contrib.contains(&(1, 1, 2)));
        assert!(contrib.contains(&(-1, 2, 1)));

        // Contributions to scalar from vectors
        let contrib = table.product_contributions(&vectors, &vectors, 0);
        assert_eq!(contrib.len(), 3); // e1*e1, e2*e2, e3*e3
    }

    #[test]
    fn all_products_vectors() {
        let algebra = Algebra::euclidean(3);
        let table = ProductTable::new(&algebra);

        let vectors = vec![1, 2, 4];
        let products = table.all_products(&vectors, &vectors);

        // vector * vector produces scalar + bivectors
        // Scalars: e1*e1, e2*e2, e3*e3
        // Bivectors: e1*e2, e2*e1, e1*e3, e3*e1, e2*e3, e3*e2

        // Check we get scalar and 3 bivectors
        let blades: Vec<usize> = products.iter().map(|(b, _)| *b).collect();
        assert!(blades.contains(&0)); // scalar
        assert!(blades.contains(&3)); // e12
        assert!(blades.contains(&5)); // e13
        assert!(blades.contains(&6)); // e23
    }

    #[test]
    fn has_contributions_to_grade_check() {
        let algebra = Algebra::euclidean(3);
        let table = ProductTable::new(&algebra);

        let vectors = vec![1, 2, 4];
        let bivectors = vec![3, 5, 6];

        // vector * vector contributes to grades 0 and 2
        assert!(table.has_contributions_to_grade(&vectors, &vectors, 0));
        assert!(table.has_contributions_to_grade(&vectors, &vectors, 2));
        assert!(!table.has_contributions_to_grade(&vectors, &vectors, 1));
        assert!(!table.has_contributions_to_grade(&vectors, &vectors, 3));

        // bivector * vector contributes to grades 1 and 3
        assert!(table.has_contributions_to_grade(&bivectors, &vectors, 1));
        assert!(table.has_contributions_to_grade(&bivectors, &vectors, 3));
        assert!(!table.has_contributions_to_grade(&bivectors, &vectors, 0));
        assert!(!table.has_contributions_to_grade(&bivectors, &vectors, 2));
    }

    #[test]
    fn pga_antiproduct_pseudoscalar_is_identity() {
        // PGA (3,0,1): e1, e2, e3 square to +1, e4 squares to 0
        // The PSEUDOSCALAR (e1234) acts as identity for antiproduct
        // (dual to scalar being identity for geometric product)
        let algebra = Algebra::pga(3);
        let table = ProductTable::new(&algebra);

        let pseudoscalar = 15; // e1234
        let e1 = 1;
        let e2 = 2;
        let e3 = 4;
        let e4 = 8;

        // Pseudoscalar ⊛ anything = ±anything (identity up to sign)
        let (sign, result) = table.antiproduct(pseudoscalar, e1);
        assert_eq!(result, e1, "pseudoscalar ⊛ e1 should give e1");
        assert_ne!(sign, 0);

        let (sign, result) = table.antiproduct(pseudoscalar, e2);
        assert_eq!(result, e2, "pseudoscalar ⊛ e2 should give e2");
        assert_ne!(sign, 0);

        let (sign, result) = table.antiproduct(pseudoscalar, e3);
        assert_eq!(result, e3, "pseudoscalar ⊛ e3 should give e3");
        assert_ne!(sign, 0);

        let (sign, result) = table.antiproduct(pseudoscalar, e4);
        assert_eq!(result, e4, "pseudoscalar ⊛ e4 should give e4");
        assert_ne!(sign, 0);
    }

    #[test]
    fn pga_antiproduct_scalar_behavior() {
        // In anti-space, scalar maps to "dual" behavior
        // scalar ⊛ non-degenerate-basis should NOT give 0
        // scalar ⊛ degenerate-basis CAN give 0 (anti-degenerate)
        let algebra = Algebra::pga(3);
        let table = ProductTable::new(&algebra);

        let scalar = 0;
        let e1 = 1; // non-degenerate
        let e4 = 8; // degenerate

        // scalar ⊛ e1 should be non-zero (the original bug fix)
        let (sign, result) = table.antiproduct(scalar, e1);
        assert_ne!(sign, 0, "scalar ⊛ e1 should not vanish in PGA");
        // The result should be the complement of e1 (e234)
        assert_eq!(result, 14, "scalar ⊛ e1 should give e234");

        // scalar ⊛ e4 gives 0 because e123 (complement of e4) is anti-degenerate
        let (sign, _result) = table.antiproduct(scalar, e4);
        assert_eq!(sign, 0, "scalar ⊛ e4 should vanish (anti-degenerate)");
    }

    #[test]
    fn pga_antiproduct_motor_point_nonzero() {
        // Motor components and point components should produce non-zero antiproducts
        let algebra = Algebra::pga(3);
        let table = ProductTable::new(&algebra);

        // Motor has scalar and e0i components (bivectors containing e4)
        // In 4D PGA: e14 = 9, e24 = 10, e34 = 12
        let e14 = 9;
        let e1 = 1;

        // e14 ⊛ e1 should be non-zero (translation affects points)
        let (sign, _result) = table.antiproduct(e14, e1);
        assert_ne!(sign, 0, "e14 ⊛ e1 should not vanish");
    }

    #[test]
    fn pga_anti_metric_contribution() {
        // Test the anti_metric_contribution directly
        let algebra = Algebra::pga(3);
        let table = ProductTable::new(&algebra);

        // In PGA (3,0,1), the anti-degenerate overlap is e123 (index 7)
        // which is all non-degenerate bases: bits 0, 1, 2
        assert_eq!(
            table.anti_metric_contribution(7),
            0,
            "e123 overlap should be degenerate in anti-space"
        );

        // Overlaps that include the degenerate basis (e4, bit 3) should NOT vanish
        // e4 alone (index 8)
        assert_ne!(
            table.anti_metric_contribution(8),
            0,
            "e4 overlap should not be degenerate in anti-space"
        );

        // e14 overlap (bits 0 and 3, index 9)
        assert_ne!(
            table.anti_metric_contribution(9),
            0,
            "e14 overlap should not be degenerate in anti-space"
        );

        // e234 overlap (bits 1, 2, 3, index 14)
        assert_ne!(
            table.anti_metric_contribution(14),
            0,
            "e234 overlap should not be degenerate in anti-space"
        );
    }

    #[test]
    fn euclidean_antiproduct_same_as_geometric() {
        // In Euclidean algebras (no degenerate directions), the antiproduct
        // should behave similarly to the geometric product (up to sign/complement)
        let algebra = Algebra::euclidean(3);
        let table = ProductTable::new(&algebra);

        // Antiscalar in 3D Euclidean is e123 (index 7)
        let antiscalar = 7;
        let e1 = 1;

        // e123 ⊛ e1 = e1 (identity property)
        assert_eq!(table.antiproduct(antiscalar, e1), (1, e1));
    }
}
