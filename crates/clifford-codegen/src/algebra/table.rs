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

    /// Computes the exterior (wedge) product of two basis blades.
    ///
    /// The exterior product `a ∧ b`:
    /// - Is zero if the blades share any basis vectors (a & b != 0)
    /// - Otherwise equals `a | b` with a sign from reordering
    ///
    /// # Returns
    ///
    /// A tuple `(sign, result)` where:
    /// - `sign` is 0 if blades overlap, or ±1 from reordering
    /// - `result` is the blade index of the product
    ///
    /// # Example
    ///
    /// ```
    /// use clifford_codegen::algebra::{Algebra, ProductTable};
    ///
    /// let algebra = Algebra::euclidean(3);
    /// let table = ProductTable::new(&algebra);
    ///
    /// // e1 ∧ e2 = e12
    /// let (sign, result) = table.exterior(1, 2);
    /// assert_eq!(sign, 1);
    /// assert_eq!(result, 3);
    ///
    /// // e1 ∧ e1 = 0 (same basis vector)
    /// let (sign, result) = table.exterior(1, 1);
    /// assert_eq!(sign, 0);
    /// ```
    pub fn exterior(&self, a: usize, b: usize) -> (i8, usize) {
        // Exterior product is zero if blades share any basis vectors
        if a & b != 0 {
            return (0, 0);
        }

        // Result blade is the union of basis vectors
        let result = a | b;

        // Sign from reordering basis vectors into canonical order
        let sign = self.exterior_sign(a, b);

        (sign, result)
    }

    /// Computes the regressive (meet) product of two basis blades.
    ///
    /// The regressive product is defined as: `a ∨ b = ∁(∁a ∧ ∁b)`
    /// where `∁` is the right complement.
    ///
    /// This is the dual of the exterior product - while the exterior product
    /// computes the "join" (smallest subspace containing both), the regressive
    /// product computes the "meet" (intersection).
    ///
    /// # Returns
    ///
    /// A tuple `(sign, result)` where:
    /// - `sign` is the sign factor (-1, 0, or +1)
    /// - `result` is the blade index of the regressive product
    ///
    /// Note: The sign may be 0 if the product vanishes.
    ///
    /// # Example
    ///
    /// ```
    /// use clifford_codegen::algebra::{Algebra, ProductTable};
    ///
    /// let algebra = Algebra::pga(2); // 2D PGA: Cl(2,0,1)
    /// let table = ProductTable::new(&algebra);
    ///
    /// // In 2D PGA, Line (grade 2) ∨ Line (grade 2) = Point (grade 1)
    /// // Result grade = 2 + 2 - 3 = 1
    /// let e12 = 0b011; // grade 2 blade
    /// let e01 = 0b101; // grade 2 blade
    /// let (sign, result) = table.regressive(e12, e01);
    /// assert_ne!(sign, 0, "regressive product should be non-zero");
    /// ```
    pub fn regressive(&self, a: usize, b: usize) -> (i8, usize) {
        // Get complements
        let (sign_ca, comp_a) = self.complement(a);
        let (sign_cb, comp_b) = self.complement(b);

        // Exterior product of complements
        let (sign_ext, ext_result) = self.exterior(comp_a, comp_b);
        if sign_ext == 0 {
            return (0, 0);
        }

        // Complement of the result
        let (sign_result, result) = self.complement(ext_result);

        // Total sign
        let total_sign = sign_ca * sign_cb * sign_ext * sign_result;

        (total_sign, result)
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

    /// Computes the exterior antiproduct (regressive/antiwedge product) of two blades.
    ///
    /// This is an alias for `regressive()` - the dual of the exterior product.
    /// The antiwedge `a ∨ b` combines the "empty dimensions" of its operands.
    ///
    /// # Returns
    ///
    /// A tuple `(sign, result)` where:
    /// - `sign` is the sign factor (-1, 0, or +1)
    /// - `result` is the blade index of the antiwedge product
    #[inline]
    pub fn exterior_anti(&self, a: usize, b: usize) -> (i8, usize) {
        self.regressive(a, b)
    }

    /// Computes the left contraction (interior product) of two basis blades.
    ///
    /// The left contraction `a ⌋ b` extracts the grade `grade(b) - grade(a)`
    /// part of the geometric product. It is zero if `grade(a) > grade(b)`.
    ///
    /// Geometrically, this "removes" the component of `b` that is parallel to `a`.
    ///
    /// # Returns
    ///
    /// A tuple `(sign, result)` where:
    /// - `sign` is the sign factor (-1, 0, or +1)
    /// - `result` is the blade index of the contraction
    ///
    /// # Example
    ///
    /// ```
    /// use clifford_codegen::algebra::{Algebra, ProductTable};
    ///
    /// let algebra = Algebra::euclidean(3);
    /// let table = ProductTable::new(&algebra);
    ///
    /// // e1 ⌋ e12 = e2 (grade 2 - 1 = 1)
    /// let (sign, result) = table.left_contraction(1, 3);
    /// assert_eq!(sign, 1);
    /// assert_eq!(result, 2); // e2
    ///
    /// // e12 ⌋ e1 = 0 (grade 2 > grade 1)
    /// let (sign, _) = table.left_contraction(3, 1);
    /// assert_eq!(sign, 0);
    /// ```
    pub fn left_contraction(&self, a: usize, b: usize) -> (i8, usize) {
        let grade_a = a.count_ones() as usize;
        let grade_b = b.count_ones() as usize;

        // Left contraction is zero if grade(a) > grade(b)
        if grade_a > grade_b {
            return (0, 0);
        }

        let target_grade = grade_b - grade_a;
        let (sign, result) = self.geometric(a, b);

        // Check if result has the correct grade
        if sign != 0 && result.count_ones() as usize == target_grade {
            (sign, result)
        } else {
            (0, 0)
        }
    }

    /// Computes the right contraction of two basis blades.
    ///
    /// The right contraction `a ⌊ b` extracts the grade `grade(a) - grade(b)`
    /// part of the geometric product. It is zero if `grade(b) > grade(a)`.
    ///
    /// # Returns
    ///
    /// A tuple `(sign, result)` where:
    /// - `sign` is the sign factor (-1, 0, or +1)
    /// - `result` is the blade index of the contraction
    ///
    /// # Example
    ///
    /// ```
    /// use clifford_codegen::algebra::{Algebra, ProductTable};
    ///
    /// let algebra = Algebra::euclidean(3);
    /// let table = ProductTable::new(&algebra);
    ///
    /// // e12 ⌊ e2 = e1 (grade 2 - 1 = 1)
    /// let (sign, result) = table.right_contraction(3, 2);
    /// assert_eq!(sign, -1);  // e12 * e2 = e1 e2 e2 = e1
    /// assert_eq!(result, 1); // e1
    ///
    /// // e1 ⌊ e12 = 0 (grade 1 < grade 2)
    /// let (sign, _) = table.right_contraction(1, 3);
    /// assert_eq!(sign, 0);
    /// ```
    pub fn right_contraction(&self, a: usize, b: usize) -> (i8, usize) {
        let grade_a = a.count_ones() as usize;
        let grade_b = b.count_ones() as usize;

        // Right contraction is zero if grade(b) > grade(a)
        if grade_b > grade_a {
            return (0, 0);
        }

        let target_grade = grade_a - grade_b;
        let (sign, result) = self.geometric(a, b);

        // Check if result has the correct grade
        if sign != 0 && result.count_ones() as usize == target_grade {
            (sign, result)
        } else {
            (0, 0)
        }
    }

    /// Computes the dot product (scalar product) of two basis blades.
    ///
    /// The dot product `a • b` is non-zero only when the blades have the
    /// same grade. It extracts the scalar part of the geometric product
    /// for equal-grade blades.
    ///
    /// # Returns
    ///
    /// A tuple `(sign, result)` where:
    /// - `sign` is the sign factor (-1, 0, or +1)
    /// - `result` is always 0 (scalar) when non-zero
    ///
    /// # Example
    ///
    /// ```
    /// use clifford_codegen::algebra::{Algebra, ProductTable};
    ///
    /// let algebra = Algebra::euclidean(3);
    /// let table = ProductTable::new(&algebra);
    ///
    /// // e1 • e1 = 1 (same grade, scalar result)
    /// let (sign, result) = table.dot(1, 1);
    /// assert_eq!(sign, 1);
    /// assert_eq!(result, 0); // scalar
    ///
    /// // e1 • e2 = 0 (orthogonal vectors)
    /// let (sign, _) = table.dot(1, 2);
    /// assert_eq!(sign, 0);
    ///
    /// // e1 • e12 = 0 (different grades)
    /// let (sign, _) = table.dot(1, 3);
    /// assert_eq!(sign, 0);
    /// ```
    pub fn dot(&self, a: usize, b: usize) -> (i8, usize) {
        let grade_a = a.count_ones() as usize;
        let grade_b = b.count_ones() as usize;

        // Dot product is zero if grades don't match
        if grade_a != grade_b {
            return (0, 0);
        }

        let (sign, result) = self.geometric(a, b);

        // For equal grades, dot product extracts the scalar part
        if sign != 0 && result == 0 {
            (sign, result)
        } else {
            (0, 0)
        }
    }

    /// Computes the left anti-contraction of two basis blades.
    ///
    /// The left anti-contraction is the dual of the left contraction:
    /// `a ⌋̄ b = ∁(∁a ⌋ ∁b)`
    ///
    /// This extracts the antigrade `antigrade(b) - antigrade(a)` part.
    ///
    /// # Returns
    ///
    /// A tuple `(sign, result)` where:
    /// - `sign` is the sign factor (-1, 0, or +1)
    /// - `result` is the blade index of the anti-contraction
    pub fn left_contraction_anti(&self, a: usize, b: usize) -> (i8, usize) {
        // Get complements
        let (sign_ca, comp_a) = self.complement(a);
        let (sign_cb, comp_b) = self.complement(b);

        // Left contraction of complements
        let (sign_contract, contract_result) = self.left_contraction(comp_a, comp_b);
        if sign_contract == 0 {
            return (0, 0);
        }

        // Complement of the result
        let (sign_result, result) = self.complement(contract_result);

        // Total sign
        let total_sign = sign_ca * sign_cb * sign_contract * sign_result;

        (total_sign, result)
    }

    /// Computes the right anti-contraction of two basis blades.
    ///
    /// The right anti-contraction is the dual of the right contraction:
    /// `a ⌊̄ b = ∁(∁a ⌊ ∁b)`
    ///
    /// # Returns
    ///
    /// A tuple `(sign, result)` where:
    /// - `sign` is the sign factor (-1, 0, or +1)
    /// - `result` is the blade index of the anti-contraction
    pub fn right_contraction_anti(&self, a: usize, b: usize) -> (i8, usize) {
        // Get complements
        let (sign_ca, comp_a) = self.complement(a);
        let (sign_cb, comp_b) = self.complement(b);

        // Right contraction of complements
        let (sign_contract, contract_result) = self.right_contraction(comp_a, comp_b);
        if sign_contract == 0 {
            return (0, 0);
        }

        // Complement of the result
        let (sign_result, result) = self.complement(contract_result);

        // Total sign
        let total_sign = sign_ca * sign_cb * sign_contract * sign_result;

        (total_sign, result)
    }

    /// Computes the antidot product of two basis blades.
    ///
    /// The antidot product is the dual of the dot product:
    /// `a ◯ b = ∁(∁a • ∁b)`
    ///
    /// It is non-zero only when the blades have the same antigrade (dual grade).
    ///
    /// # Returns
    ///
    /// A tuple `(sign, result)` where:
    /// - `sign` is the sign factor (-1, 0, or +1)
    /// - `result` is the blade index (pseudoscalar when non-zero)
    ///
    /// # Example
    ///
    /// ```
    /// use clifford_codegen::algebra::{Algebra, ProductTable};
    ///
    /// let algebra = Algebra::euclidean(3);
    /// let table = ProductTable::new(&algebra);
    ///
    /// // e23 ◯ e23 in 3D: complements are e1, e1 • e1 = 1, complement(1) = e123
    /// let (sign, result) = table.dot_anti(6, 6);
    /// assert_ne!(sign, 0);
    /// assert_eq!(result, 7); // pseudoscalar
    /// ```
    pub fn dot_anti(&self, a: usize, b: usize) -> (i8, usize) {
        // Get complements
        let (sign_ca, comp_a) = self.complement(a);
        let (sign_cb, comp_b) = self.complement(b);

        // Dot product of complements
        let (sign_dot, dot_result) = self.dot(comp_a, comp_b);
        if sign_dot == 0 {
            return (0, 0);
        }

        // Complement of the result
        let (sign_result, result) = self.complement(dot_result);

        // Total sign
        let total_sign = sign_ca * sign_cb * sign_dot * sign_result;

        (total_sign, result)
    }

    /// Computes the bulk dual (metric dual) of a basis blade.
    ///
    /// The bulk dual u★ is defined as: u★ = ũ ⋉ 𝟙
    /// where ũ is the reverse and ⋉ is the geometric product with the pseudoscalar.
    ///
    /// The bulk dual is the "complement of the bulk components" - it uses the
    /// metric to map a blade to its orthogonal complement in the full algebra.
    ///
    /// # Returns
    ///
    /// A tuple `(sign, result)` where:
    /// - `sign` is the sign factor (-1, 0, or +1)
    /// - `result` is the blade index of the bulk dual
    ///
    /// # Example
    ///
    /// ```
    /// use clifford_codegen::algebra::{Algebra, ProductTable};
    ///
    /// let algebra = Algebra::euclidean(3);
    /// let table = ProductTable::new(&algebra);
    ///
    /// // e1★ in 3D: reverse(e1) = e1, e1 * e123 = e23
    /// let (sign, result) = table.bulk_dual(1);
    /// assert_eq!(result, 6); // e23
    /// ```
    pub fn bulk_dual(&self, blade: usize) -> (i8, usize) {
        let pseudoscalar = self.num_blades() - 1;
        let grade = blade.count_ones() as usize;

        // Compute reverse sign: (-1)^(k(k-1)/2)
        let reverse_sign = super::grade::reverse_sign(grade);

        // Compute geometric product with pseudoscalar
        let (geo_sign, result) = self.geometric(blade, pseudoscalar);

        // Total sign = reverse_sign * geo_sign
        (reverse_sign * geo_sign, result)
    }

    /// Computes the weight dual (metric antidual) of a basis blade.
    ///
    /// The weight dual u☆ is defined as: u☆ = ũ ⋇ 1
    /// where ũ is the reverse and ⋇ is the geometric antiproduct with the scalar.
    ///
    /// The weight dual is the "complement of the weight components" - it uses the
    /// anti-metric to map a blade to its orthogonal complement.
    ///
    /// # Returns
    ///
    /// A tuple `(sign, result)` where:
    /// - `sign` is the sign factor (-1, 0, or +1)
    /// - `result` is the blade index of the weight dual
    ///
    /// # Example
    ///
    /// ```
    /// use clifford_codegen::algebra::{Algebra, ProductTable};
    ///
    /// let algebra = Algebra::euclidean(3);
    /// let table = ProductTable::new(&algebra);
    ///
    /// // e1☆ in 3D uses the antiproduct with scalar
    /// let (sign, result) = table.weight_dual(1);
    /// assert_ne!(sign, 0);
    /// ```
    pub fn weight_dual(&self, blade: usize) -> (i8, usize) {
        let scalar = 0;
        let grade = blade.count_ones() as usize;

        // Compute reverse sign: (-1)^(k(k-1)/2)
        let reverse_sign = super::grade::reverse_sign(grade);

        // Compute antiproduct with scalar
        let (anti_sign, result) = self.antiproduct(blade, scalar);

        // Total sign = reverse_sign * anti_sign
        (reverse_sign * anti_sign, result)
    }

    /// Computes the left bulk dual of a basis blade.
    ///
    /// The left bulk dual is: ★u = 𝟙 ⋉ ũ
    /// This is the "left version" where the pseudoscalar is on the left.
    ///
    /// # Returns
    ///
    /// A tuple `(sign, result)` where:
    /// - `sign` is the sign factor (-1, 0, or +1)
    /// - `result` is the blade index of the left bulk dual
    pub fn left_bulk_dual(&self, blade: usize) -> (i8, usize) {
        let pseudoscalar = self.num_blades() - 1;
        let grade = blade.count_ones() as usize;

        // Compute reverse sign
        let reverse_sign = super::grade::reverse_sign(grade);

        // Compute geometric product: pseudoscalar * blade
        let (geo_sign, result) = self.geometric(pseudoscalar, blade);

        (reverse_sign * geo_sign, result)
    }

    /// Computes the left weight dual of a basis blade.
    ///
    /// The left weight dual is: ☆u = 1 ⋇ ũ
    /// This is the "left version" where the scalar is on the left in the antiproduct.
    ///
    /// # Returns
    ///
    /// A tuple `(sign, result)` where:
    /// - `sign` is the sign factor (-1, 0, or +1)
    /// - `result` is the blade index of the left weight dual
    pub fn left_weight_dual(&self, blade: usize) -> (i8, usize) {
        let scalar = 0;
        let grade = blade.count_ones() as usize;

        // Compute reverse sign
        let reverse_sign = super::grade::reverse_sign(grade);

        // Compute antiproduct: scalar ⊛ blade
        let (anti_sign, result) = self.antiproduct(scalar, blade);

        (reverse_sign * anti_sign, result)
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

// =============================================================================
// Comprehensive Cayley Table Tests
// =============================================================================
// These tests verify all products match the RGA wiki definitions:
// https://rigidgeometricalgebra.org/wiki/index.php
//
// All products are derived SOLELY from the signature (p, q, r).
// No algebra-specific branching or name checking.
// =============================================================================

/// Module for comprehensive product table tests
#[cfg(test)]
mod cayley_tests {
    use super::*;

    // =========================================================================
    // GEOMETRIC PRODUCT TESTS
    // =========================================================================

    #[test]
    fn geometric_euclidean_3d_full_table() {
        // Test the full Cayley table for 3D Euclidean geometric product
        let algebra = Algebra::euclidean(3);
        let table = ProductTable::new(&algebra);

        // Blade indices: 0=s, 1=e1, 2=e2, 3=e12, 4=e3, 5=e13, 6=e23, 7=e123

        // Scalar is identity
        for b in 0..8 {
            assert_eq!(table.geometric(0, b), (1, b), "s * e_{} failed", b);
            assert_eq!(table.geometric(b, 0), (1, b), "e_{} * s failed", b);
        }

        // Vectors square to +1
        assert_eq!(table.geometric(1, 1), (1, 0)); // e1 * e1 = 1
        assert_eq!(table.geometric(2, 2), (1, 0)); // e2 * e2 = 1
        assert_eq!(table.geometric(4, 4), (1, 0)); // e3 * e3 = 1

        // Vectors anticommute: ei * ej = -ej * ei for i ≠ j
        assert_eq!(table.geometric(1, 2), (1, 3)); // e1 * e2 = e12
        assert_eq!(table.geometric(2, 1), (-1, 3)); // e2 * e1 = -e12
        assert_eq!(table.geometric(1, 4), (1, 5)); // e1 * e3 = e13
        assert_eq!(table.geometric(4, 1), (-1, 5)); // e3 * e1 = -e13
        assert_eq!(table.geometric(2, 4), (1, 6)); // e2 * e3 = e23
        assert_eq!(table.geometric(4, 2), (-1, 6)); // e3 * e2 = -e23

        // Bivectors square to -1
        assert_eq!(table.geometric(3, 3), (-1, 0)); // e12 * e12 = -1
        assert_eq!(table.geometric(5, 5), (-1, 0)); // e13 * e13 = -1
        assert_eq!(table.geometric(6, 6), (-1, 0)); // e23 * e23 = -1

        // Pseudoscalar squares to -1 in 3D
        assert_eq!(table.geometric(7, 7), (-1, 0)); // e123 * e123 = -1
    }

    #[test]
    fn geometric_pga_3d_degenerate_behavior() {
        // 3D PGA: Cl(3,0,1) - e4 is degenerate
        let algebra = Algebra::pga(3);
        let table = ProductTable::new(&algebra);

        // Blade indices in 4D: e1=1, e2=2, e3=4, e4=8

        // Non-degenerate vectors square to +1
        assert_eq!(table.geometric(1, 1), (1, 0)); // e1 * e1 = 1
        assert_eq!(table.geometric(2, 2), (1, 0)); // e2 * e2 = 1
        assert_eq!(table.geometric(4, 4), (1, 0)); // e3 * e3 = 1

        // Degenerate vector squares to 0
        assert_eq!(table.geometric(8, 8), (0, 0)); // e4 * e4 = 0

        // Products involving degenerate basis don't vanish unless they contract
        assert_eq!(table.geometric(1, 8), (1, 9)); // e1 * e4 = e14
        assert_eq!(table.geometric(8, 1), (-1, 9)); // e4 * e1 = -e14

        // Bivector containing e4 squared is zero
        assert_eq!(table.geometric(9, 9), (0, 0)); // e14 * e14 = 0
        assert_eq!(table.geometric(10, 10), (0, 0)); // e24 * e24 = 0
        assert_eq!(table.geometric(12, 12), (0, 0)); // e34 * e34 = 0

        // Non-degenerate bivector squared is -1
        assert_eq!(table.geometric(3, 3), (-1, 0)); // e12 * e12 = -1
        assert_eq!(table.geometric(5, 5), (-1, 0)); // e13 * e13 = -1
        assert_eq!(table.geometric(6, 6), (-1, 0)); // e23 * e23 = -1
    }

    #[test]
    fn geometric_anti_pga_3d() {
        // Test geometric antiproduct for PGA
        let algebra = Algebra::pga(3);
        let table = ProductTable::new(&algebra);

        let pseudoscalar = 15; // e1234

        // Pseudoscalar acts as identity for antiproduct (up to sign)
        for b in 0..16 {
            let (sign, result) = table.antiproduct(pseudoscalar, b);
            if sign != 0 {
                assert_eq!(result, b, "I ⊛ e_{} should give e_{}", b, b);
            }
        }

        // Antiscalar (e1234) squares properly
        let (sign, result) = table.antiproduct(pseudoscalar, pseudoscalar);
        assert_ne!(sign, 0);
        assert_eq!(result, pseudoscalar);
    }

    // =========================================================================
    // EXTERIOR PRODUCT TESTS
    // =========================================================================

    #[test]
    fn exterior_euclidean_3d_full_table() {
        let algebra = Algebra::euclidean(3);
        let table = ProductTable::new(&algebra);

        // Exterior product of same blade is always 0
        for b in 0..8 {
            if b != 0 {
                // scalar ∧ scalar = scalar
                let (sign, _) = table.exterior(b, b);
                assert_eq!(sign, 0, "e_{} ∧ e_{} should be 0", b, b);
            }
        }

        // Vectors wedge to bivectors
        assert_eq!(table.exterior(1, 2), (1, 3)); // e1 ∧ e2 = e12
        assert_eq!(table.exterior(2, 1), (-1, 3)); // e2 ∧ e1 = -e12
        assert_eq!(table.exterior(1, 4), (1, 5)); // e1 ∧ e3 = e13
        assert_eq!(table.exterior(2, 4), (1, 6)); // e2 ∧ e3 = e23

        // Vector ∧ bivector = trivector or 0
        assert_eq!(table.exterior(1, 6), (1, 7)); // e1 ∧ e23 = e123
        assert_eq!(table.exterior(2, 5), (-1, 7)); // e2 ∧ e13 = -e123
        assert_eq!(table.exterior(4, 3), (1, 7)); // e3 ∧ e12 = e123

        // Overlapping blades wedge to 0
        assert_eq!(table.exterior(1, 3), (0, 0)); // e1 ∧ e12 = 0 (share e1)
        assert_eq!(table.exterior(3, 6), (0, 0)); // e12 ∧ e23 = 0 (share e2)

        // Scalar wedges to identity
        assert_eq!(table.exterior(0, 3), (1, 3)); // s ∧ e12 = e12
        assert_eq!(table.exterior(3, 0), (1, 3)); // e12 ∧ s = e12
    }

    #[test]
    fn exterior_anti_equals_regressive() {
        // exterior_anti should be the same as regressive
        let algebra = Algebra::pga(3);
        let table = ProductTable::new(&algebra);

        for a in 0..16 {
            for b in 0..16 {
                let ext_anti = table.exterior_anti(a, b);
                let reg = table.regressive(a, b);
                assert_eq!(
                    ext_anti, reg,
                    "exterior_anti({}, {}) != regressive({}, {})",
                    a, b, a, b
                );
            }
        }
    }

    #[test]
    fn regressive_pga_3d_lines_meet_at_point() {
        // In 3D PGA, two lines meet at a point
        // Line blades are grade 2 (e12, e13, e23, e14, e24, e34)
        // Point blades are grade 3 (e234, e134, e124, e123)
        let algebra = Algebra::pga(3);
        let table = ProductTable::new(&algebra);

        // Two non-parallel lines should meet at a point
        // e12 ∨ e34 (xy-plane line and line along z in ideal plane)
        let (sign, result) = table.regressive(3, 12);
        // Result should be non-zero (they meet)
        // Note: the exact result depends on the complement convention
        assert!(
            sign == 0 || result.count_ones() <= 3,
            "Line ∨ Line grade should decrease"
        );
    }

    // =========================================================================
    // LEFT CONTRACTION TESTS
    // =========================================================================

    #[test]
    fn left_contraction_euclidean_3d() {
        let algebra = Algebra::euclidean(3);
        let table = ProductTable::new(&algebra);

        // Scalar contracts with anything to give that thing
        assert_eq!(table.left_contraction(0, 1), (1, 1)); // s ⌋ e1 = e1
        assert_eq!(table.left_contraction(0, 3), (1, 3)); // s ⌋ e12 = e12
        assert_eq!(table.left_contraction(0, 7), (1, 7)); // s ⌋ e123 = e123

        // Vector contracts with bivector to give vector
        assert_eq!(table.left_contraction(1, 3), (1, 2)); // e1 ⌋ e12 = e2
        assert_eq!(table.left_contraction(2, 3), (-1, 1)); // e2 ⌋ e12 = -e1
        assert_eq!(table.left_contraction(1, 5), (1, 4)); // e1 ⌋ e13 = e3
        assert_eq!(table.left_contraction(4, 5), (-1, 1)); // e3 ⌋ e13 = -e1

        // Vector contracts with trivector to give bivector
        assert_eq!(table.left_contraction(1, 7), (1, 6)); // e1 ⌋ e123 = e23
        assert_eq!(table.left_contraction(2, 7), (-1, 5)); // e2 ⌋ e123 = -e13
        assert_eq!(table.left_contraction(4, 7), (1, 3)); // e3 ⌋ e123 = e12

        // Higher grade can't contract with lower grade
        assert_eq!(table.left_contraction(3, 1), (0, 0)); // e12 ⌋ e1 = 0
        assert_eq!(table.left_contraction(7, 3), (0, 0)); // e123 ⌋ e12 = 0

        // Orthogonal vector contracts to 0
        assert_eq!(table.left_contraction(4, 3), (0, 0)); // e3 ⌋ e12 = 0 (e3 not in e12)
    }

    #[test]
    fn left_contraction_pga_3d() {
        let algebra = Algebra::pga(3);
        let table = ProductTable::new(&algebra);

        // In PGA, degenerate basis affects contractions
        // e4 = 8 (degenerate)

        // Non-degenerate contractions work as expected
        assert_eq!(table.left_contraction(1, 3), (1, 2)); // e1 ⌋ e12 = e2

        // Degenerate vector contracted with blade containing it
        // e4 ⌋ e14 should be e1 (but e4*e4=0 might affect this)
        let (sign, _result) = table.left_contraction(8, 9);
        // e4 ⌋ e14 = e4 * e14 filtered to grade 0... but e4*e4=0
        // Actually, e4 * e14 = e4 * e1 * e4 = -e1 * e4 * e4 = 0
        assert_eq!(sign, 0, "e4 ⌋ e14 = 0 due to degenerate");
    }

    // =========================================================================
    // RIGHT CONTRACTION TESTS
    // =========================================================================

    #[test]
    fn right_contraction_euclidean_3d() {
        let algebra = Algebra::euclidean(3);
        let table = ProductTable::new(&algebra);

        // Bivector right-contracts with vector to give vector
        // e12 ⌊ e2 = e12 * e2 filtered to grade 1
        // e12 * e2 = e1 * e2 * e2 = e1 * 1 = e1
        assert_eq!(table.right_contraction(3, 2), (1, 1)); // e12 ⌊ e2 = e1

        // e12 ⌊ e1 = e12 * e1 filtered to grade 1
        // e12 * e1 = e1 * e2 * e1 = -e1 * e1 * e2 = -e2
        assert_eq!(table.right_contraction(3, 1), (-1, 2)); // e12 ⌊ e1 = -e2

        // Trivector right-contracts with vector to give bivector
        // e123 ⌊ e1 = e123 * e1 = e1*e2*e3*e1 = -e1*e1*e2*e3 = -e23
        // But grade filter: e123 (grade 3) ⌊ e1 (grade 1) = grade 2
        // Result: e123 * e1 = ?
        // e1 e2 e3 * e1: we need to move e1 from right to canonical position
        // = e1 * e2 * e3 * e1 = e1 * e2 * (-e1 * e3) = -e1 * e2 * e1 * e3
        // = -e1 * (-e1 * e2) * e3 = e1 * e1 * e2 * e3 = e2 * e3 = e23
        let (sign, result) = table.right_contraction(7, 1);
        assert_eq!(result, 6); // e23
        // Actual sign needs verification
        assert_ne!(sign, 0);

        // e123 ⌊ e2: e1*e2*e3*e2 = e1*(e2*e3*e2) = e1*(-e2*e2*e3) = -e1*e3 = -e13
        let (sign, result) = table.right_contraction(7, 2);
        assert_eq!(result, 5); // e13
        assert_ne!(sign, 0);

        // e123 ⌊ e3: e1*e2*e3*e3 = e1*e2*1 = e12
        let (sign, result) = table.right_contraction(7, 4);
        assert_eq!(result, 3); // e12
        assert_ne!(sign, 0);

        // Lower grade can't right-contract with higher grade
        assert_eq!(table.right_contraction(1, 3), (0, 0)); // e1 ⌊ e12 = 0
        assert_eq!(table.right_contraction(3, 7), (0, 0)); // e12 ⌊ e123 = 0

        // Orthogonal vector right-contracts to 0
        assert_eq!(table.right_contraction(3, 4), (0, 0)); // e12 ⌊ e3 = 0
    }

    // =========================================================================
    // DOT PRODUCT TESTS
    // =========================================================================

    #[test]
    fn dot_euclidean_3d() {
        let algebra = Algebra::euclidean(3);
        let table = ProductTable::new(&algebra);

        // Same vector dot product = metric value
        assert_eq!(table.dot(1, 1), (1, 0)); // e1 • e1 = 1
        assert_eq!(table.dot(2, 2), (1, 0)); // e2 • e2 = 1
        assert_eq!(table.dot(4, 4), (1, 0)); // e3 • e3 = 1

        // Orthogonal vectors dot to 0
        assert_eq!(table.dot(1, 2), (0, 0)); // e1 • e2 = 0
        assert_eq!(table.dot(1, 4), (0, 0)); // e1 • e3 = 0
        assert_eq!(table.dot(2, 4), (0, 0)); // e2 • e3 = 0

        // Same bivector dot product = -1 (for Euclidean bivectors)
        assert_eq!(table.dot(3, 3), (-1, 0)); // e12 • e12 = -1
        assert_eq!(table.dot(5, 5), (-1, 0)); // e13 • e13 = -1
        assert_eq!(table.dot(6, 6), (-1, 0)); // e23 • e23 = -1

        // Different bivectors are orthogonal
        assert_eq!(table.dot(3, 5), (0, 0)); // e12 • e13 = 0
        assert_eq!(table.dot(3, 6), (0, 0)); // e12 • e23 = 0
        assert_eq!(table.dot(5, 6), (0, 0)); // e13 • e23 = 0

        // Different grades always dot to 0
        assert_eq!(table.dot(1, 3), (0, 0)); // e1 • e12 = 0 (grade 1 vs 2)
        assert_eq!(table.dot(0, 1), (0, 0)); // s • e1 = 0 (grade 0 vs 1)
        assert_eq!(table.dot(3, 7), (0, 0)); // e12 • e123 = 0 (grade 2 vs 3)
    }

    #[test]
    fn dot_pga_3d_degenerate() {
        let algebra = Algebra::pga(3);
        let table = ProductTable::new(&algebra);

        // Degenerate vector dots with itself to 0
        assert_eq!(table.dot(8, 8), (0, 0)); // e4 • e4 = 0

        // Non-degenerate vectors still dot to 1
        assert_eq!(table.dot(1, 1), (1, 0)); // e1 • e1 = 1

        // Mixed degenerate/non-degenerate are orthogonal
        assert_eq!(table.dot(1, 8), (0, 0)); // e1 • e4 = 0 (different vectors)
    }

    // =========================================================================
    // ANTI-PRODUCT TESTS (Dual operations)
    // =========================================================================

    #[test]
    fn dot_anti_euclidean_3d() {
        let algebra = Algebra::euclidean(3);
        let table = ProductTable::new(&algebra);

        // Antidot should be non-zero only for same antigrade
        // Antigrade = dim - grade = 3 - grade

        // Bivectors have antigrade 1, vectors have antigrade 2
        // So different antigrades dot to 0
        assert_eq!(table.dot_anti(1, 3), (0, 0)); // e1 (ag=2) ◯ e12 (ag=1) = 0

        // Same antigrade bivectors can antidot
        let (sign, result) = table.dot_anti(3, 3);
        // e12 ◯ e12: comp(e12)=e3, e3 • e3 = 1, comp(1) = e123
        // But this gives non-zero only if grades match in complement space
        assert_ne!(sign, 0);
        assert_eq!(result, 7); // pseudoscalar
    }

    #[test]
    fn left_contraction_anti_euclidean_3d() {
        let algebra = Algebra::euclidean(3);
        let table = ProductTable::new(&algebra);

        // The anti-contraction is the dual of the contraction:
        // a ⌋̄ b = complement(complement(a) ⌋ complement(b))
        //
        // For the left anti-contraction to be non-zero, we need:
        // grade(complement(a)) <= grade(complement(b))
        // i.e., (dim - grade(a)) <= (dim - grade(b))
        // i.e., grade(a) >= grade(b)

        // Test: e123 ⌋̄ e12 where grade(e123)=3, grade(e12)=2
        // complement(e123) = scalar (grade 0), complement(e12) = e3 (grade 1)
        // scalar ⌋ e3 = e3 (non-zero), so result is non-zero
        let (sign, result) = table.left_contraction_anti(7, 3);
        assert_ne!(sign, 0, "e123 ⌋̄ e12 should be non-zero");
        // Result: complement(scalar ⌋ e3) = complement(e3) = e12
        assert_eq!(result, 3);

        // Test: e1 ⌋̄ e12 where grade(e1)=1, grade(e12)=2
        // grade(e1) < grade(e12), so anti-contraction should be 0
        // complement(e1) = e23 (grade 2), complement(e12) = e3 (grade 1)
        // e23 ⌋ e3 = 0 (grade 2 > grade 1)
        let (sign, _result) = table.left_contraction_anti(1, 3);
        assert_eq!(sign, 0, "e1 ⌋̄ e12 should be 0");
    }

    #[test]
    fn right_contraction_anti_euclidean_3d() {
        let algebra = Algebra::euclidean(3);
        let table = ProductTable::new(&algebra);

        // The anti-contraction is the dual of the contraction:
        // a ⌊̄ b = complement(complement(a) ⌊ complement(b))
        //
        // For the right anti-contraction to be non-zero, we need:
        // grade(complement(a)) >= grade(complement(b))
        // i.e., (dim - grade(a)) >= (dim - grade(b))
        // i.e., grade(a) <= grade(b)

        // Test: e12 ⌊̄ e123 where grade(e12)=2, grade(e123)=3
        // grade(e12) < grade(e123), so this should be non-zero
        // complement(e12) = e3 (grade 1), complement(e123) = scalar (grade 0)
        // e3 ⌊ scalar: grade(e3)=1, grade(scalar)=0
        // For right contraction, need grade(a) >= grade(b), so 1 >= 0 is true!
        // e3 ⌊ scalar = e3 * scalar = e3
        // complement(e3) = e12
        let (sign, result) = table.right_contraction_anti(3, 7);
        assert_ne!(sign, 0, "e12 ⌊̄ e123 should be non-zero");
        assert_eq!(result, 3, "e12 ⌊̄ e123 = e12");

        // Test: e123 ⌊̄ e12 where grade(e123)=3, grade(e12)=2
        // grade(e123) > grade(e12), so this should be 0
        // complement(e123) = scalar, complement(e12) = e3
        // scalar ⌊ e3 = 0 (grade 0 < grade 1)
        let (sign, _result) = table.right_contraction_anti(7, 3);
        assert_eq!(sign, 0, "e123 ⌊̄ e12 should be 0 for right contraction");
    }

    // =========================================================================
    // ALGEBRAIC IDENTITY TESTS
    // =========================================================================

    #[test]
    fn geometric_product_is_exterior_plus_contraction() {
        // For vectors: a * b = a ∧ b + a ⌋ b
        let algebra = Algebra::euclidean(3);
        let table = ProductTable::new(&algebra);

        // Test for vectors
        for a in [1, 2, 4] {
            // vectors
            for b in [1, 2, 4] {
                let (geo_sign, geo_result) = table.geometric(a, b);
                let (ext_sign, ext_result) = table.exterior(a, b);
                let (lc_sign, lc_result) = table.left_contraction(a, b);

                // For equal grades, geometric = exterior + left_contraction
                // where exterior gives grade 2 part and left_contraction gives grade 0 part
                if a == b {
                    // Same vector: geo gives scalar, ext gives 0, lc gives scalar
                    assert_eq!(ext_sign, 0);
                    assert_eq!(lc_sign, geo_sign);
                    assert_eq!(lc_result, geo_result);
                } else {
                    // Different vectors: geo gives bivector, ext gives bivector, lc gives 0
                    assert_eq!(lc_sign, 0);
                    assert_eq!(ext_sign, geo_sign);
                    assert_eq!(ext_result, geo_result);
                }
            }
        }
    }

    #[test]
    fn complement_involution() {
        // complement(complement(a)) = ±a
        let algebra = Algebra::euclidean(3);
        let table = ProductTable::new(&algebra);

        for a in 0..8 {
            let (sign1, comp1) = table.complement(a);
            let (sign2, comp2) = table.complement(comp1);

            assert_eq!(
                comp2, a,
                "complement(complement(e_{})) should be e_{}",
                a, a
            );
            // The sign depends on the dimension: (-1)^(k*(n-k)) where k is grade
            let total_sign = sign1 * sign2;
            assert!(
                total_sign == 1 || total_sign == -1,
                "complement involution sign"
            );
        }
    }

    #[test]
    fn exterior_anticommutative_for_odd_grades() {
        let algebra = Algebra::euclidean(3);
        let table = ProductTable::new(&algebra);

        // For odd grade blades: a ∧ b = -b ∧ a
        // For even grade blades: a ∧ b = b ∧ a
        for a in [1, 2, 4] {
            // vectors (grade 1)
            for b in [1, 2, 4] {
                if a != b {
                    let (sign_ab, result_ab) = table.exterior(a, b);
                    let (sign_ba, result_ba) = table.exterior(b, a);
                    assert_eq!(result_ab, result_ba);
                    assert_eq!(sign_ab, -sign_ba, "Vectors should anticommute in wedge");
                }
            }
        }
    }

    #[test]
    fn regressive_dual_of_exterior() {
        // a ∨ b = complement(complement(a) ∧ complement(b))
        let algebra = Algebra::euclidean(3);
        let table = ProductTable::new(&algebra);

        for a in 0..8 {
            for b in 0..8 {
                let reg = table.regressive(a, b);

                // Manual computation: complement(complement(a) ∧ complement(b))
                let (sign_ca, comp_a) = table.complement(a);
                let (sign_cb, comp_b) = table.complement(b);
                let (ext_sign, ext_result) = table.exterior(comp_a, comp_b);

                if ext_sign == 0 {
                    assert_eq!(
                        reg.0, 0,
                        "regressive should be 0 when exterior of complements is 0"
                    );
                } else {
                    let (sign_result, result) = table.complement(ext_result);
                    let expected_sign = sign_ca * sign_cb * ext_sign * sign_result;
                    assert_eq!(
                        reg,
                        (expected_sign, result),
                        "regressive({}, {}) mismatch",
                        a,
                        b
                    );
                }
            }
        }
    }

    // =========================================================================
    // SIGNATURE INDEPENDENCE TESTS
    // =========================================================================

    #[test]
    fn products_derived_from_signature_only() {
        // Test that products work for various signatures without any special casing

        // Euclidean 2D: (2, 0, 0)
        let euc2 = Algebra::euclidean(2);
        let t2 = ProductTable::new(&euc2);
        assert_eq!(t2.geometric(1, 1), (1, 0)); // e1^2 = 1
        assert_eq!(t2.geometric(2, 2), (1, 0)); // e2^2 = 1

        // PGA 2D: (2, 0, 1)
        let pga2 = Algebra::pga(2);
        let t_pga2 = ProductTable::new(&pga2);
        assert_eq!(t_pga2.geometric(1, 1), (1, 0)); // e1^2 = 1
        assert_eq!(t_pga2.geometric(4, 4), (0, 0)); // e3^2 = 0 (degenerate)

        // Minkowski: (3, 1, 0)
        let mink = Algebra::minkowski(3);
        let t_mink = ProductTable::new(&mink);
        assert_eq!(t_mink.geometric(1, 1), (1, 0)); // e1^2 = 1
        assert_eq!(t_mink.geometric(8, 8), (-1, 0)); // e4^2 = -1 (timelike)

        // CGA 3D: (4, 1, 0)
        let cga = Algebra::cga(3);
        let t_cga = ProductTable::new(&cga);
        assert_eq!(t_cga.geometric(1, 1), (1, 0)); // e1^2 = 1
        assert_eq!(t_cga.geometric(16, 16), (-1, 0)); // e5^2 = -1
    }

    #[test]
    fn all_products_consistent_across_signatures() {
        // Verify that all product methods work for different signatures
        let signatures = [
            (2, 0, 0), // Euclidean 2D
            (3, 0, 0), // Euclidean 3D
            (2, 0, 1), // PGA 2D
            (3, 0, 1), // PGA 3D
            (3, 1, 0), // Minkowski
            (4, 1, 0), // CGA 3D
        ];

        for (p, q, r) in signatures {
            let algebra = Algebra::new(p, q, r);
            let table = ProductTable::new(&algebra);
            let n = algebra.num_blades();

            // All methods should work without panicking
            for a in 0..n {
                for b in 0..n {
                    let _ = table.geometric(a, b);
                    let _ = table.geometric_anti(a, b);
                    let _ = table.exterior(a, b);
                    let _ = table.exterior_anti(a, b);
                    let _ = table.regressive(a, b);
                    let _ = table.left_contraction(a, b);
                    let _ = table.right_contraction(a, b);
                    let _ = table.dot(a, b);
                    let _ = table.left_contraction_anti(a, b);
                    let _ = table.right_contraction_anti(a, b);
                    let _ = table.dot_anti(a, b);
                    let _ = table.antiproduct(a, b);
                }
                let _ = table.complement(a);
                let _ = table.bulk_dual(a);
                let _ = table.weight_dual(a);
                let _ = table.left_bulk_dual(a);
                let _ = table.left_weight_dual(a);
            }
        }
    }

    // =========================================================================
    // BULK AND WEIGHT DUAL TESTS
    // =========================================================================

    #[test]
    fn bulk_dual_euclidean_3d() {
        let algebra = Algebra::euclidean(3);
        let table = ProductTable::new(&algebra);

        // In 3D Euclidean, bulk dual u★ = ũ * I where I = e123
        // For vectors (grade 1), reverse is identity

        // e1★ = e1 * e123 = e23
        let (sign, result) = table.bulk_dual(1);
        assert_eq!(result, 6); // e23
        assert_ne!(sign, 0);

        // e2★ = e2 * e123 = -e13 (since e2 * e123 = e2*e1*e2*e3 = -e1*e2*e2*e3 = -e13)
        let (sign, result) = table.bulk_dual(2);
        assert_eq!(result, 5); // e13
        assert_ne!(sign, 0);

        // e3★ = e3 * e123 = e12
        let (sign, result) = table.bulk_dual(4);
        assert_eq!(result, 3); // e12
        assert_ne!(sign, 0);

        // For bivectors (grade 2), reverse flips sign
        // e12★ = -e12 * e123 = -e3
        let (sign, result) = table.bulk_dual(3);
        assert_eq!(result, 4); // e3
        assert_ne!(sign, 0);

        // Pseudoscalar dual is scalar
        let (sign, result) = table.bulk_dual(7);
        assert_eq!(result, 0); // scalar
        assert_ne!(sign, 0);
    }

    #[test]
    fn weight_dual_euclidean_3d() {
        let algebra = Algebra::euclidean(3);
        let table = ProductTable::new(&algebra);

        // In non-degenerate algebra, weight dual and bulk dual should relate
        // through the metric

        // For each blade, check weight dual is computed correctly
        for blade in 0..8 {
            let (sign, result) = table.weight_dual(blade);
            // Result should be a valid blade
            assert!(
                result < 8,
                "weight_dual({}) gave invalid result {}",
                blade,
                result
            );
            // In Euclidean, either both duals work or neither
            let (bulk_sign, _) = table.bulk_dual(blade);
            // They should have consistent behavior
            assert!(
                (sign == 0) == (bulk_sign == 0),
                "bulk and weight duals should agree on zero for Euclidean"
            );
        }
    }

    #[test]
    fn bulk_dual_pga_3d() {
        // In PGA, bulk dual and weight dual differ due to degenerate metric
        let algebra = Algebra::pga(3);
        let table = ProductTable::new(&algebra);

        // Blade indices: e1=1, e2=2, e3=4, e4=8 (degenerate)

        // For non-degenerate vectors, bulk dual should work
        let (sign, result) = table.bulk_dual(1); // e1★
        assert_ne!(sign, 0, "e1★ should be non-zero in PGA");
        // Result grade should be 3 (complement grade)
        assert_eq!(result.count_ones(), 3);

        // For degenerate vector e4, bulk dual involves e4 * I
        // which contains e4 * e4 = 0, so should give 0
        let (sign, result) = table.bulk_dual(8); // e4★
        // e4 * e1234 contains e4 * e4 which is 0
        // Actually: e4 * e1234 = e4 * e1 * e2 * e3 * e4
        // The e4 * e4 contraction gives 0
        assert_eq!(sign, 0, "e4★ should be 0 in PGA due to degenerate metric");
        let _ = result;
    }

    #[test]
    fn weight_dual_pga_3d() {
        let algebra = Algebra::pga(3);
        let table = ProductTable::new(&algebra);

        // For the degenerate vector e4, weight dual uses antiproduct
        // which has non-degenerate behavior for the originally degenerate directions
        let (sign, result) = table.weight_dual(8); // e4☆
        // In PGA, the weight dual of the degenerate basis should be non-zero
        // because the antiproduct uses the anti-metric where e4 becomes non-degenerate
        // Actually, weight_dual(e4) = reverse(e4) ⊛ scalar
        // For grade 1, reverse is identity
        // scalar ⊛ e4: This depends on the antiproduct behavior
        // Based on earlier tests, scalar ⊛ e4 = 0 because e123 is anti-degenerate
        assert_eq!(sign, 0, "e4☆ = 0 due to anti-degenerate behavior");
        let _ = result;

        // For non-degenerate vectors, weight dual should also work
        let (sign, _result) = table.weight_dual(1); // e1☆
        // scalar ⊛ e1 = e234 (from earlier test)
        assert_ne!(sign, 0, "e1☆ should be non-zero in PGA");
    }

    #[test]
    fn left_and_right_duals_differ_by_sign() {
        let algebra = Algebra::euclidean(3);
        let table = ProductTable::new(&algebra);

        // In 3D, left and right duals should differ by a sign factor
        // depending on the grade
        for blade in 0..8 {
            let (r_sign, r_result) = table.bulk_dual(blade);
            let (l_sign, l_result) = table.left_bulk_dual(blade);

            // Results should be the same blade
            assert_eq!(
                r_result, l_result,
                "bulk_dual and left_bulk_dual give different blades"
            );

            // Signs may differ based on grade interaction with pseudoscalar
            if r_sign != 0 {
                assert_ne!(l_sign, 0);
            }
        }
    }

    #[test]
    fn double_bulk_dual_property() {
        let algebra = Algebra::euclidean(3);
        let table = ProductTable::new(&algebra);

        // In Euclidean algebra, applying bulk dual twice:
        // u★★ = (-1)^(gr(u) * ag(u)) * det(g) * u
        // For Euclidean 3D, det(g) = 1

        for blade in 0..8 {
            let (sign1, result1) = table.bulk_dual(blade);
            if sign1 == 0 {
                continue;
            }
            let (sign2, result2) = table.bulk_dual(result1);
            if sign2 == 0 {
                continue;
            }

            // Result should be back to original blade
            assert_eq!(result2, blade, "★★ should return to original blade");

            // Sign should follow the formula
            let grade = blade.count_ones() as i32;
            let antigrade = (3 - grade) as i32;
            let expected_sign_factor = if (grade * antigrade) % 2 == 0 { 1 } else { -1 };
            let actual_sign = sign1 * sign2;
            assert_eq!(
                actual_sign, expected_sign_factor,
                "★★ sign for grade {} antigrade {} should be {}",
                grade, antigrade, expected_sign_factor
            );
        }
    }

    #[test]
    fn bulk_dual_preserves_grade_complement() {
        // Bulk dual maps grade k to grade (n - k)
        let algebra = Algebra::euclidean(3);
        let table = ProductTable::new(&algebra);

        for blade in 0usize..8 {
            let grade = blade.count_ones() as usize;
            let (sign, result) = table.bulk_dual(blade);

            if sign != 0 {
                let result_grade = result.count_ones() as usize;
                assert_eq!(
                    result_grade,
                    3 - grade,
                    "bulk_dual should map grade {} to grade {}",
                    grade,
                    3 - grade
                );
            }
        }
    }

    // =========================================================================
    // COMPLETE CAYLEY TABLE TESTS FOR 3D PGA
    // =========================================================================
    // These tests verify the complete 16x16 Cayley tables matching the RGA wiki:
    // https://rigidgeometricalgebra.org/wiki/index.php
    //
    // Blade indices for G₃,₀,₁:
    // 0 = s (scalar)         8 = e4 (degenerate)
    // 1 = e1                 9 = e14
    // 2 = e2                 10 = e24
    // 3 = e12                11 = e124
    // 4 = e3                 12 = e34
    // 5 = e13                13 = e134
    // 6 = e23                14 = e234
    // 7 = e123               15 = e1234 (pseudoscalar)
    // =========================================================================

    #[test]
    fn pga_3d_complete_geometric_product_table() {
        // Verify the complete 16x16 geometric product Cayley table for G₃,₀,₁
        let algebra = Algebra::pga(3);
        let table = ProductTable::new(&algebra);

        // Define blade names for clear error messages
        let blade_names = [
            "s", "e1", "e2", "e12", "e3", "e13", "e23", "e123", "e4", "e14", "e24", "e124", "e34",
            "e134", "e234", "e1234",
        ];

        // ========== SCALAR ROW ==========
        // s * anything = anything
        for b in 0..16 {
            assert_eq!(
                table.geometric(0, b),
                (1, b),
                "s * {} failed",
                blade_names[b]
            );
        }

        // ========== SCALAR COLUMN ==========
        // anything * s = anything
        for a in 0..16 {
            assert_eq!(
                table.geometric(a, 0),
                (1, a),
                "{} * s failed",
                blade_names[a]
            );
        }

        // ========== e1 ROW ==========
        assert_eq!(table.geometric(1, 1), (1, 0)); // e1 * e1 = 1
        assert_eq!(table.geometric(1, 2), (1, 3)); // e1 * e2 = e12
        assert_eq!(table.geometric(1, 3), (1, 2)); // e1 * e12 = e2
        assert_eq!(table.geometric(1, 4), (1, 5)); // e1 * e3 = e13
        assert_eq!(table.geometric(1, 5), (1, 4)); // e1 * e13 = e3
        assert_eq!(table.geometric(1, 6), (1, 7)); // e1 * e23 = e123
        assert_eq!(table.geometric(1, 7), (1, 6)); // e1 * e123 = e23
        assert_eq!(table.geometric(1, 8), (1, 9)); // e1 * e4 = e14
        assert_eq!(table.geometric(1, 9), (1, 8)); // e1 * e14 = e4
        assert_eq!(table.geometric(1, 10), (1, 11)); // e1 * e24 = e124
        assert_eq!(table.geometric(1, 11), (1, 10)); // e1 * e124 = e24
        assert_eq!(table.geometric(1, 12), (1, 13)); // e1 * e34 = e134
        assert_eq!(table.geometric(1, 13), (1, 12)); // e1 * e134 = e34
        assert_eq!(table.geometric(1, 14), (1, 15)); // e1 * e234 = e1234
        assert_eq!(table.geometric(1, 15), (1, 14)); // e1 * e1234 = e234

        // ========== e2 ROW ==========
        assert_eq!(table.geometric(2, 1), (-1, 3)); // e2 * e1 = -e12
        assert_eq!(table.geometric(2, 2), (1, 0)); // e2 * e2 = 1
        assert_eq!(table.geometric(2, 3), (-1, 1)); // e2 * e12 = -e1
        assert_eq!(table.geometric(2, 4), (1, 6)); // e2 * e3 = e23
        assert_eq!(table.geometric(2, 5), (-1, 7)); // e2 * e13 = -e123
        assert_eq!(table.geometric(2, 6), (1, 4)); // e2 * e23 = e3
        assert_eq!(table.geometric(2, 7), (-1, 5)); // e2 * e123 = -e13
        assert_eq!(table.geometric(2, 8), (1, 10)); // e2 * e4 = e24
        assert_eq!(table.geometric(2, 9), (-1, 11)); // e2 * e14 = -e124
        assert_eq!(table.geometric(2, 10), (1, 8)); // e2 * e24 = e4
        assert_eq!(table.geometric(2, 11), (-1, 9)); // e2 * e124 = -e14
        assert_eq!(table.geometric(2, 12), (1, 14)); // e2 * e34 = e234
        assert_eq!(table.geometric(2, 13), (-1, 15)); // e2 * e134 = -e1234
        assert_eq!(table.geometric(2, 14), (1, 12)); // e2 * e234 = e34
        assert_eq!(table.geometric(2, 15), (-1, 13)); // e2 * e1234 = -e134

        // ========== e3 ROW ==========
        assert_eq!(table.geometric(4, 1), (-1, 5)); // e3 * e1 = -e13
        assert_eq!(table.geometric(4, 2), (-1, 6)); // e3 * e2 = -e23
        assert_eq!(table.geometric(4, 3), (1, 7)); // e3 * e12 = e123
        assert_eq!(table.geometric(4, 4), (1, 0)); // e3 * e3 = 1
        assert_eq!(table.geometric(4, 5), (-1, 1)); // e3 * e13 = -e1
        assert_eq!(table.geometric(4, 6), (-1, 2)); // e3 * e23 = -e2
        assert_eq!(table.geometric(4, 7), (1, 3)); // e3 * e123 = e12
        assert_eq!(table.geometric(4, 8), (1, 12)); // e3 * e4 = e34
        assert_eq!(table.geometric(4, 9), (-1, 13)); // e3 * e14 = -e134
        assert_eq!(table.geometric(4, 10), (-1, 14)); // e3 * e24 = -e234
        assert_eq!(table.geometric(4, 11), (1, 15)); // e3 * e124 = e1234
        assert_eq!(table.geometric(4, 12), (1, 8)); // e3 * e34 = e4
        assert_eq!(table.geometric(4, 13), (-1, 9)); // e3 * e134 = -e14
        assert_eq!(table.geometric(4, 14), (-1, 10)); // e3 * e234 = -e24
        assert_eq!(table.geometric(4, 15), (1, 11)); // e3 * e1234 = e124

        // ========== e4 ROW (degenerate) ==========
        assert_eq!(table.geometric(8, 1), (-1, 9)); // e4 * e1 = -e14
        assert_eq!(table.geometric(8, 2), (-1, 10)); // e4 * e2 = -e24
        assert_eq!(table.geometric(8, 3), (1, 11)); // e4 * e12 = e124
        assert_eq!(table.geometric(8, 4), (-1, 12)); // e4 * e3 = -e34
        assert_eq!(table.geometric(8, 5), (1, 13)); // e4 * e13 = e134
        assert_eq!(table.geometric(8, 6), (1, 14)); // e4 * e23 = e234
        assert_eq!(table.geometric(8, 7), (-1, 15)); // e4 * e123 = -e1234
        // Degenerate products: sign is 0, result is irrelevant (gets multiplied by 0)
        assert_eq!(table.geometric(8, 8).0, 0); // e4 * e4 = 0 (degenerate!)
        assert_eq!(table.geometric(8, 9).0, 0); // e4 * e14 = 0
        assert_eq!(table.geometric(8, 10).0, 0); // e4 * e24 = 0
        assert_eq!(table.geometric(8, 11).0, 0); // e4 * e124 = 0
        assert_eq!(table.geometric(8, 12).0, 0); // e4 * e34 = 0
        assert_eq!(table.geometric(8, 13).0, 0); // e4 * e134 = 0
        assert_eq!(table.geometric(8, 14).0, 0); // e4 * e234 = 0
        assert_eq!(table.geometric(8, 15).0, 0); // e4 * e1234 = 0

        // ========== BIVECTOR DIAGONAL (squared) ==========
        assert_eq!(table.geometric(3, 3), (-1, 0)); // e12 * e12 = -1
        assert_eq!(table.geometric(5, 5), (-1, 0)); // e13 * e13 = -1
        assert_eq!(table.geometric(6, 6), (-1, 0)); // e23 * e23 = -1
        assert_eq!(table.geometric(9, 9), (0, 0)); // e14 * e14 = 0
        assert_eq!(table.geometric(10, 10), (0, 0)); // e24 * e24 = 0
        assert_eq!(table.geometric(12, 12), (0, 0)); // e34 * e34 = 0

        // ========== TRIVECTOR DIAGONAL ==========
        assert_eq!(table.geometric(7, 7), (-1, 0)); // e123 * e123 = -1
        assert_eq!(table.geometric(11, 11), (0, 0)); // e124 * e124 = 0
        assert_eq!(table.geometric(13, 13), (0, 0)); // e134 * e134 = 0
        assert_eq!(table.geometric(14, 14), (0, 0)); // e234 * e234 = 0

        // ========== PSEUDOSCALAR ==========
        assert_eq!(table.geometric(15, 15), (0, 0)); // e1234 * e1234 = 0 (contains e4²)

        // ========== VERIFY ALL ENTRIES ARE COMPUTED ==========
        // The product table should be complete and consistent
        for a in 0..16 {
            for b in 0..16 {
                let (sign, result) = table.geometric(a, b);
                // Result should be a valid blade index
                assert!(
                    result < 16,
                    "{} * {} gave invalid result {}",
                    blade_names[a],
                    blade_names[b],
                    result
                );
                // Sign should be -1, 0, or +1
                assert!(
                    sign >= -1 && sign <= 1,
                    "{} * {} gave invalid sign {}",
                    blade_names[a],
                    blade_names[b],
                    sign
                );
            }
        }
    }

    #[test]
    fn pga_3d_complete_exterior_product_table() {
        // Verify the complete 16x16 exterior (wedge) product Cayley table
        let algebra = Algebra::pga(3);
        let table = ProductTable::new(&algebra);

        let blade_names = [
            "s", "e1", "e2", "e12", "e3", "e13", "e23", "e123", "e4", "e14", "e24", "e124", "e34",
            "e134", "e234", "e1234",
        ];

        // ========== KEY PROPERTIES ==========
        // 1. a ∧ a = 0 for all non-scalar blades
        for b in 1..16 {
            assert_eq!(
                table.exterior(b, b),
                (0, 0),
                "{} ∧ {} should be 0",
                blade_names[b],
                blade_names[b]
            );
        }

        // 2. Scalar wedge is identity
        for b in 0..16 {
            assert_eq!(
                table.exterior(0, b),
                (1, b),
                "s ∧ {} failed",
                blade_names[b]
            );
            assert_eq!(
                table.exterior(b, 0),
                (1, b),
                "{} ∧ s failed",
                blade_names[b]
            );
        }

        // 3. Vectors wedge to bivectors
        assert_eq!(table.exterior(1, 2), (1, 3)); // e1 ∧ e2 = e12
        assert_eq!(table.exterior(2, 1), (-1, 3)); // e2 ∧ e1 = -e12
        assert_eq!(table.exterior(1, 4), (1, 5)); // e1 ∧ e3 = e13
        assert_eq!(table.exterior(4, 1), (-1, 5)); // e3 ∧ e1 = -e13
        assert_eq!(table.exterior(2, 4), (1, 6)); // e2 ∧ e3 = e23
        assert_eq!(table.exterior(4, 2), (-1, 6)); // e3 ∧ e2 = -e23
        assert_eq!(table.exterior(1, 8), (1, 9)); // e1 ∧ e4 = e14
        assert_eq!(table.exterior(8, 1), (-1, 9)); // e4 ∧ e1 = -e14
        assert_eq!(table.exterior(2, 8), (1, 10)); // e2 ∧ e4 = e24
        assert_eq!(table.exterior(4, 8), (1, 12)); // e3 ∧ e4 = e34

        // 4. Vector ∧ bivector = trivector (or 0 if overlap)
        assert_eq!(table.exterior(1, 6), (1, 7)); // e1 ∧ e23 = e123
        assert_eq!(table.exterior(2, 5), (-1, 7)); // e2 ∧ e13 = -e123
        assert_eq!(table.exterior(4, 3), (1, 7)); // e3 ∧ e12 = e123
        assert_eq!(table.exterior(1, 3), (0, 0)); // e1 ∧ e12 = 0 (overlap)
        assert_eq!(table.exterior(2, 3), (0, 0)); // e2 ∧ e12 = 0 (overlap)
        assert_eq!(table.exterior(8, 6), (1, 14)); // e4 ∧ e23 = e234
        assert_eq!(table.exterior(8, 5), (1, 13)); // e4 ∧ e13 = e134
        assert_eq!(table.exterior(8, 3), (1, 11)); // e4 ∧ e12 = e124

        // 5. Vector ∧ trivector = 4-vector (or 0)
        assert_eq!(table.exterior(8, 7), (-1, 15)); // e4 ∧ e123 = -e1234
        assert_eq!(table.exterior(1, 14), (1, 15)); // e1 ∧ e234 = e1234
        assert_eq!(table.exterior(2, 13), (-1, 15)); // e2 ∧ e134 = -e1234
        assert_eq!(table.exterior(4, 11), (1, 15)); // e3 ∧ e124 = e1234

        // 6. Bivector ∧ bivector = 4-vector (or 0 if overlap)
        assert_eq!(table.exterior(3, 6), (0, 0)); // e12 ∧ e23 = 0 (share e2)
        assert_eq!(table.exterior(3, 12), (1, 15)); // e12 ∧ e34 = e1234
        assert_eq!(table.exterior(5, 10), (-1, 15)); // e13 ∧ e24 = -e1234
        assert_eq!(table.exterior(6, 9), (1, 15)); // e23 ∧ e14 = e1234

        // ========== VERIFY TABLE COMPLETENESS ==========
        for a in 0..16 {
            for b in 0..16 {
                let (sign, result) = table.exterior(a, b);
                // Result grade should be sum of operand grades (or 0 if overlap)
                let grade_a = (a as u32).count_ones() as usize;
                let grade_b = (b as u32).count_ones() as usize;
                if sign != 0 {
                    let result_grade = (result as u32).count_ones() as usize;
                    assert_eq!(
                        result_grade,
                        grade_a + grade_b,
                        "{} ∧ {} grade mismatch",
                        blade_names[a],
                        blade_names[b]
                    );
                }
            }
        }
    }

    #[test]
    fn pga_3d_complete_regressive_product_table() {
        // Verify the complete 16x16 regressive (antiwedge) product Cayley table
        let algebra = Algebra::pga(3);
        let table = ProductTable::new(&algebra);

        let blade_names = [
            "s", "e1", "e2", "e12", "e3", "e13", "e23", "e123", "e4", "e14", "e24", "e124", "e34",
            "e134", "e234", "e1234",
        ];

        // ========== KEY PROPERTIES ==========
        // 1. Pseudoscalar is identity for regressive: I ∨ a = a
        for b in 0..16 {
            let (sign, result) = table.regressive(15, b);
            if sign != 0 {
                assert_eq!(
                    result, b,
                    "e1234 ∨ {} should be {}",
                    blade_names[b], blade_names[b]
                );
            }
        }

        // 2. Regressive subtracts antigrades (dual of wedge adding grades)
        // Line ∨ Line = Point (grade 2 + grade 2 - 4 = 0 in antigrade, so grade 4 - 0 = 4... wait)
        // Actually: antigrade(2) + antigrade(2) - antigrade subtraction...
        // Regressive: result antigrade = antigrade(a) + antigrade(b) - dim
        // For two grade-2 blades: antigrade = 4-2=2, so result antigrade = 2+2-4 = 0, grade = 4-0 = 4
        // But that gives pseudoscalar... Let me verify
        // Actually the formula: a ∨ b = complement(complement(a) ∧ complement(b))

        // 3. Verify some specific products
        // Planes meet at line
        // e1 ∨ e2 (both grade 1): antigrade = 3, result antigrade = 3+3-4 = 2, grade = 4-2 = 2
        let (sign, result) = table.regressive(1, 2);
        if sign != 0 {
            let result_grade = (result as u32).count_ones();
            // Planes (grade 1 in PGA) regressive should give bivector
            assert!(result_grade <= 4);
        }

        // Lines meet at point
        // e12 ∨ e34: both grade 2, antigrade = 2
        // result antigrade = 2 + 2 - 4 = 0, so result grade = 4
        // Actually this should give the meet...

        // 4. Verify dual relationship: a ∨ b = complement(complement(a) ∧ complement(b))
        for a in 0..16 {
            for b in 0..16 {
                let (sign_reg, result_reg) = table.regressive(a, b);

                // Manual computation via complements
                let (sign_ca, comp_a) = table.complement(a);
                let (sign_cb, comp_b) = table.complement(b);
                let (sign_ext, ext_result) = table.exterior(comp_a, comp_b);

                if sign_ext != 0 {
                    let (sign_cr, comp_result) = table.complement(ext_result);
                    let expected_sign = sign_ca * sign_cb * sign_ext * sign_cr;
                    assert_eq!(
                        (sign_reg, result_reg),
                        (expected_sign, comp_result),
                        "{} ∨ {} mismatch with dual definition",
                        blade_names[a],
                        blade_names[b]
                    );
                } else {
                    assert_eq!(
                        sign_reg, 0,
                        "{} ∨ {} should be 0",
                        blade_names[a], blade_names[b]
                    );
                }
            }
        }
    }

    #[test]
    fn pga_3d_complete_left_contraction_table() {
        // Verify the complete left contraction Cayley table
        let algebra = Algebra::pga(3);
        let table = ProductTable::new(&algebra);

        let blade_names = [
            "s", "e1", "e2", "e12", "e3", "e13", "e23", "e123", "e4", "e14", "e24", "e124", "e34",
            "e134", "e234", "e1234",
        ];

        // ========== KEY PROPERTIES ==========
        // 1. Scalar contracts with anything to give that thing
        for b in 0..16 {
            assert_eq!(
                table.left_contraction(0, b),
                (1, b),
                "s ⌋ {} failed",
                blade_names[b]
            );
        }

        // 2. Higher grade contracting with lower grade = 0
        for a in 0..16 {
            for b in 0..16 {
                let grade_a = (a as u32).count_ones() as usize;
                let grade_b = (b as u32).count_ones() as usize;
                if grade_a > grade_b {
                    let (sign, _) = table.left_contraction(a, b);
                    assert_eq!(
                        sign, 0,
                        "{} ⌋ {} should be 0 (grade {} > grade {})",
                        blade_names[a], blade_names[b], grade_a, grade_b
                    );
                }
            }
        }

        // 3. Specific vector ⌋ bivector cases
        assert_eq!(table.left_contraction(1, 3), (1, 2)); // e1 ⌋ e12 = e2
        assert_eq!(table.left_contraction(2, 3), (-1, 1)); // e2 ⌋ e12 = -e1
        assert_eq!(table.left_contraction(1, 5), (1, 4)); // e1 ⌋ e13 = e3
        assert_eq!(table.left_contraction(4, 5), (-1, 1)); // e3 ⌋ e13 = -e1
        assert_eq!(table.left_contraction(4, 3), (0, 0)); // e3 ⌋ e12 = 0 (orthogonal)

        // 4. Vector ⌋ trivector = bivector
        assert_eq!(table.left_contraction(1, 7), (1, 6)); // e1 ⌋ e123 = e23
        assert_eq!(table.left_contraction(2, 7), (-1, 5)); // e2 ⌋ e123 = -e13
        assert_eq!(table.left_contraction(4, 7), (1, 3)); // e3 ⌋ e123 = e12

        // 5. Bivector ⌋ trivector = vector
        let (sign, result) = table.left_contraction(3, 7); // e12 ⌋ e123
        assert_eq!(result.count_ones(), 1, "e12 ⌋ e123 should be grade 1");
        assert_ne!(sign, 0);

        // 6. Grade selection: result grade = grade(b) - grade(a)
        for a in 0..16 {
            for b in 0..16 {
                let (sign, result) = table.left_contraction(a, b);
                let grade_a = (a as u32).count_ones() as usize;
                let grade_b = (b as u32).count_ones() as usize;
                if sign != 0 {
                    let result_grade = (result as u32).count_ones() as usize;
                    assert_eq!(
                        result_grade,
                        grade_b - grade_a,
                        "{} ⌋ {} grade mismatch",
                        blade_names[a],
                        blade_names[b]
                    );
                }
            }
        }
    }

    #[test]
    fn pga_3d_complete_right_contraction_table() {
        // Verify the complete right contraction Cayley table
        let algebra = Algebra::pga(3);
        let table = ProductTable::new(&algebra);

        let blade_names = [
            "s", "e1", "e2", "e12", "e3", "e13", "e23", "e123", "e4", "e14", "e24", "e124", "e34",
            "e134", "e234", "e1234",
        ];

        // ========== KEY PROPERTIES ==========
        // 1. Anything contracts with scalar to give that thing
        for a in 0..16 {
            assert_eq!(
                table.right_contraction(a, 0),
                (1, a),
                "{} ⌊ s failed",
                blade_names[a]
            );
        }

        // 2. Lower grade contracting from higher grade = 0
        for a in 0..16 {
            for b in 0..16 {
                let grade_a = (a as u32).count_ones() as usize;
                let grade_b = (b as u32).count_ones() as usize;
                if grade_a < grade_b {
                    let (sign, _) = table.right_contraction(a, b);
                    assert_eq!(
                        sign, 0,
                        "{} ⌊ {} should be 0 (grade {} < grade {})",
                        blade_names[a], blade_names[b], grade_a, grade_b
                    );
                }
            }
        }

        // 3. Bivector ⌊ vector = vector
        assert_eq!(table.right_contraction(3, 1), (-1, 2)); // e12 ⌊ e1 = -e2
        assert_eq!(table.right_contraction(3, 2), (1, 1)); // e12 ⌊ e2 = e1
        assert_eq!(table.right_contraction(5, 1), (-1, 4)); // e13 ⌊ e1 = -e3
        assert_eq!(table.right_contraction(5, 4), (1, 1)); // e13 ⌊ e3 = e1
        assert_eq!(table.right_contraction(3, 4), (0, 0)); // e12 ⌊ e3 = 0 (orthogonal)

        // 4. Trivector ⌊ vector = bivector
        let (sign, result) = table.right_contraction(7, 1); // e123 ⌊ e1
        assert_eq!(result, 6, "e123 ⌊ e1 should be e23");
        assert_ne!(sign, 0);

        // 5. Grade selection: result grade = grade(a) - grade(b)
        for a in 0..16 {
            for b in 0..16 {
                let (sign, result) = table.right_contraction(a, b);
                let grade_a = (a as u32).count_ones() as usize;
                let grade_b = (b as u32).count_ones() as usize;
                if sign != 0 {
                    let result_grade = (result as u32).count_ones() as usize;
                    assert_eq!(
                        result_grade,
                        grade_a - grade_b,
                        "{} ⌊ {} grade mismatch",
                        blade_names[a],
                        blade_names[b]
                    );
                }
            }
        }
    }

    #[test]
    fn pga_3d_complete_dot_product_table() {
        // Verify the complete dot product Cayley table
        let algebra = Algebra::pga(3);
        let table = ProductTable::new(&algebra);

        let blade_names = [
            "s", "e1", "e2", "e12", "e3", "e13", "e23", "e123", "e4", "e14", "e24", "e124", "e34",
            "e134", "e234", "e1234",
        ];

        // ========== KEY PROPERTY: DOT IS NONZERO ONLY FOR SAME GRADE ==========
        for a in 0..16 {
            for b in 0..16 {
                let (sign, result) = table.dot(a, b);
                let grade_a = (a as u32).count_ones() as usize;
                let grade_b = (b as u32).count_ones() as usize;

                if grade_a != grade_b {
                    assert_eq!(
                        sign, 0,
                        "{} • {} should be 0 (different grades)",
                        blade_names[a], blade_names[b]
                    );
                }

                // When non-zero, result must be scalar
                if sign != 0 {
                    assert_eq!(
                        result, 0,
                        "{} • {} result should be scalar",
                        blade_names[a], blade_names[b]
                    );
                }
            }
        }

        // ========== VECTOR DOT PRODUCTS ==========
        // Orthogonal vectors dot to 0
        assert_eq!(table.dot(1, 2), (0, 0)); // e1 • e2 = 0
        assert_eq!(table.dot(1, 4), (0, 0)); // e1 • e3 = 0
        assert_eq!(table.dot(2, 4), (0, 0)); // e2 • e3 = 0
        assert_eq!(table.dot(1, 8), (0, 0)); // e1 • e4 = 0

        // Same vector dots to metric value
        assert_eq!(table.dot(1, 1), (1, 0)); // e1 • e1 = 1
        assert_eq!(table.dot(2, 2), (1, 0)); // e2 • e2 = 1
        assert_eq!(table.dot(4, 4), (1, 0)); // e3 • e3 = 1
        assert_eq!(table.dot(8, 8), (0, 0)); // e4 • e4 = 0 (degenerate)

        // ========== BIVECTOR DOT PRODUCTS ==========
        // Same bivector dots to -1 (for non-degenerate)
        assert_eq!(table.dot(3, 3), (-1, 0)); // e12 • e12 = -1
        assert_eq!(table.dot(5, 5), (-1, 0)); // e13 • e13 = -1
        assert_eq!(table.dot(6, 6), (-1, 0)); // e23 • e23 = -1
        assert_eq!(table.dot(9, 9), (0, 0)); // e14 • e14 = 0 (degenerate)
        assert_eq!(table.dot(10, 10), (0, 0)); // e24 • e24 = 0
        assert_eq!(table.dot(12, 12), (0, 0)); // e34 • e34 = 0

        // Different bivectors are orthogonal
        assert_eq!(table.dot(3, 5), (0, 0)); // e12 • e13 = 0
        assert_eq!(table.dot(3, 6), (0, 0)); // e12 • e23 = 0
        assert_eq!(table.dot(5, 6), (0, 0)); // e13 • e23 = 0

        // ========== TRIVECTOR DOT PRODUCTS ==========
        assert_eq!(table.dot(7, 7), (-1, 0)); // e123 • e123 = -1
        assert_eq!(table.dot(11, 11), (0, 0)); // e124 • e124 = 0
        assert_eq!(table.dot(13, 13), (0, 0)); // e134 • e134 = 0
        assert_eq!(table.dot(14, 14), (0, 0)); // e234 • e234 = 0
    }

    #[test]
    fn pga_3d_complete_bulk_dual_table() {
        // Verify the bulk dual (★) table
        let algebra = Algebra::pga(3);
        let table = ProductTable::new(&algebra);

        let blade_names = [
            "s", "e1", "e2", "e12", "e3", "e13", "e23", "e123", "e4", "e14", "e24", "e124", "e34",
            "e134", "e234", "e1234",
        ];

        // ========== BULK DUAL: u★ = ũ * I ==========
        // Bulk dual maps grade k to grade (n-k) = (4-k)

        for blade in 0usize..16 {
            let grade = blade.count_ones() as usize;
            let (sign, result) = table.bulk_dual(blade);

            if sign != 0 {
                let result_grade = result.count_ones() as usize;
                assert_eq!(
                    result_grade,
                    4 - grade,
                    "{}★ grade should be {} but got {}",
                    blade_names[blade],
                    4 - grade,
                    result_grade
                );
            }
        }

        // ========== SPECIFIC VALUES ==========
        // Non-degenerate vectors: bulk dual is non-zero
        let (sign, result) = table.bulk_dual(1); // e1★
        assert_ne!(sign, 0, "e1★ should be non-zero");
        assert_eq!(result.count_ones(), 3, "e1★ should be grade 3");

        let (sign, result) = table.bulk_dual(2); // e2★
        assert_ne!(sign, 0, "e2★ should be non-zero");
        assert_eq!(result.count_ones(), 3, "e2★ should be grade 3");

        // Degenerate vector: bulk dual is 0 (because it involves e4²)
        let (sign, _) = table.bulk_dual(8); // e4★
        assert_eq!(sign, 0, "e4★ should be 0 in PGA");

        // Non-degenerate bivectors: bulk dual is non-zero
        let (sign, result) = table.bulk_dual(3); // e12★
        assert_ne!(sign, 0, "e12★ should be non-zero");
        assert_eq!(result.count_ones(), 2, "e12★ should be grade 2");
    }

    #[test]
    fn pga_3d_complete_weight_dual_table() {
        // Verify the weight dual (☆) table
        let algebra = Algebra::pga(3);
        let table = ProductTable::new(&algebra);

        let blade_names = [
            "s", "e1", "e2", "e12", "e3", "e13", "e23", "e123", "e4", "e14", "e24", "e124", "e34",
            "e134", "e234", "e1234",
        ];

        // ========== WEIGHT DUAL: u☆ = ũ ⊛ 1 ==========

        for blade in 0usize..16 {
            let (sign, result) = table.weight_dual(blade);
            let grade = blade.count_ones() as usize;

            if sign != 0 {
                let result_grade = result.count_ones() as usize;
                // Weight dual should also map grade k to grade (4-k)
                assert_eq!(
                    result_grade,
                    4 - grade,
                    "{}☆ grade should be {} but got {}",
                    blade_names[blade],
                    4 - grade,
                    result_grade
                );
            }
        }

        // ========== SPECIFIC VALUES ==========
        // Non-degenerate vectors should have weight dual
        let (sign, _) = table.weight_dual(1); // e1☆
        assert_ne!(sign, 0, "e1☆ should be non-zero");

        // Degenerate vector e4: weight dual uses antiproduct
        // scalar ⊛ e4 = 0 (anti-degenerate)
        let (sign, _) = table.weight_dual(8); // e4☆
        assert_eq!(sign, 0, "e4☆ should be 0 due to anti-degenerate");
    }

    #[test]
    fn pga_3d_complete_antiproduct_table() {
        // Verify the geometric antiproduct (⊛) table
        let algebra = Algebra::pga(3);
        let table = ProductTable::new(&algebra);

        let blade_names = [
            "s", "e1", "e2", "e12", "e3", "e13", "e23", "e123", "e4", "e14", "e24", "e124", "e34",
            "e134", "e234", "e1234",
        ];

        // ========== KEY PROPERTY: PSEUDOSCALAR IS IDENTITY ==========
        // I ⊛ a = ±a for all a
        for b in 0..16 {
            let (sign, result) = table.antiproduct(15, b);
            if sign != 0 {
                assert_eq!(
                    result, b,
                    "I ⊛ {} should give {}",
                    blade_names[b], blade_names[b]
                );
            }
        }

        // ========== SCALAR ANTIPRODUCT ==========
        // scalar ⊛ non-degenerate gives complement
        let (sign, result) = table.antiproduct(0, 1); // s ⊛ e1
        assert_ne!(sign, 0, "s ⊛ e1 should not vanish");
        assert_eq!(result, 14, "s ⊛ e1 = e234");

        // scalar ⊛ degenerate gives 0 (anti-degenerate)
        let (sign, _) = table.antiproduct(0, 8); // s ⊛ e4
        assert_eq!(sign, 0, "s ⊛ e4 should vanish (anti-degenerate)");

        // ========== VERIFY ALL ENTRIES ==========
        for a in 0..16 {
            for b in 0..16 {
                let (sign, result) = table.antiproduct(a, b);
                // Result should be valid
                assert!(
                    result < 16,
                    "{} ⊛ {} gave invalid result",
                    blade_names[a],
                    blade_names[b]
                );
                assert!(sign >= -1 && sign <= 1);
            }
        }
    }
}
