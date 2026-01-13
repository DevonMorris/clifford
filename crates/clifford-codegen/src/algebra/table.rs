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
    /// `a ⊛ b = ∁(∁a × ∁b)`
    ///
    /// where `×` is the **regular** geometric product and `∁` is the complement.
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

        // Geometric product of complements using REGULAR metric
        let (sign_prod, prod) = self.geometric(comp_a, comp_b);
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
    /// assert_eq!(sign, 1);   // e12 * e2 = e1 e2 e2 = e1 * (+1) = e1
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

    /// Computes the interior product (symmetric contraction) of two basis blades.
    ///
    /// The interior product extracts the grade `|grade(a) - grade(b)|` part
    /// of the geometric product. It's the symmetric version of the contractions.
    ///
    /// # Returns
    ///
    /// A tuple `(sign, result)` where:
    /// - `sign` is the sign factor (-1, 0, or +1)
    /// - `result` is the blade index of the interior product
    ///
    /// # Example
    ///
    /// ```
    /// use clifford_codegen::algebra::{Algebra, ProductTable};
    ///
    /// let algebra = Algebra::euclidean(3);
    /// let table = ProductTable::new(&algebra);
    ///
    /// // e12 · e2 = e1 (grade |2 - 1| = 1)
    /// let (sign, result) = table.interior(3, 2);
    /// assert_eq!(sign, 1);
    /// assert_eq!(result, 1); // e1
    ///
    /// // e1 · e2 = 0 (orthogonal, grade |1 - 1| = 0 but no scalar part)
    /// let (sign, _) = table.interior(1, 2);
    /// assert_eq!(sign, 0);
    /// ```
    pub fn interior(&self, a: usize, b: usize) -> (i8, usize) {
        let grade_a = a.count_ones() as usize;
        let grade_b = b.count_ones() as usize;
        let target_grade = grade_a.abs_diff(grade_b);

        let (sign, result) = self.geometric(a, b);

        if sign != 0 && result.count_ones() as usize == target_grade {
            (sign, result)
        } else {
            (0, 0)
        }
    }

    /// Computes the scalar product of two basis blades.
    ///
    /// The scalar product extracts only the grade-0 (scalar) part of the
    /// geometric product. This is equivalent to the dot product for same-grade
    /// blades, but returns zero for different grades.
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
    /// // e1 * e1 = 1 (scalar)
    /// let (sign, result) = table.scalar(1, 1);
    /// assert_eq!(sign, 1);
    /// assert_eq!(result, 0);
    ///
    /// // e1 * e2 has no scalar part
    /// let (sign, _) = table.scalar(1, 2);
    /// assert_eq!(sign, 0);
    /// ```
    pub fn scalar(&self, a: usize, b: usize) -> (i8, usize) {
        let (sign, result) = self.geometric(a, b);

        if sign != 0 && result == 0 {
            (sign, 0)
        } else {
            (0, 0)
        }
    }

    /// Computes the antiscalar product of two basis blades.
    ///
    /// The antiscalar product extracts only the grade-n (pseudoscalar) part
    /// of the antiproduct, where n is the dimension of the algebra.
    ///
    /// # Returns
    ///
    /// A tuple `(sign, result)` where:
    /// - `sign` is the sign factor (-1, 0, or +1)
    /// - `result` is the pseudoscalar blade index when non-zero
    ///
    /// # Example
    ///
    /// ```
    /// use clifford_codegen::algebra::{Algebra, ProductTable};
    ///
    /// let algebra = Algebra::pga(3);
    /// let table = ProductTable::new(&algebra);
    ///
    /// // Pseudoscalar index for 4D algebra (PGA 3D has 4 basis vectors)
    /// let pseudoscalar = (1 << 4) - 1; // 0b1111 = 15
    ///
    /// // e1234 ∨ e1234 = e1234 (antiscalar part)
    /// // In the antiproduct, pseudoscalar acts like scalar does in geometric product
    /// let (sign, result) = table.antiscalar(pseudoscalar, pseudoscalar);
    /// assert_eq!(sign, 1);
    /// assert_eq!(result, pseudoscalar);
    ///
    /// // 1 ∨ e1234 has no antiscalar part (result is scalar, not pseudoscalar)
    /// let (sign, _) = table.antiscalar(0, pseudoscalar);
    /// assert_eq!(sign, 0);
    /// ```
    pub fn antiscalar(&self, a: usize, b: usize) -> (i8, usize) {
        let dim = self.dim;
        let pseudoscalar = (1 << dim) - 1;

        let (sign, result) = self.antiproduct(a, b);

        if sign != 0 && result == pseudoscalar {
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
        // With the correct antiproduct formula (∁(∁a × ∁b) using regular geometric):
        // scalar ⊛ non-degenerate-basis: complement(e1234 × complement(ei)) vanishes
        //   because e1234 × e_ijk involves degenerate e0 in the metric
        // scalar ⊛ degenerate-basis (e4/e0): complement(e1234 × e123) is NON-ZERO
        //   because e1234 × e123 = ±e0 (only Euclidean bases in metric)
        let algebra = Algebra::pga(3);
        let table = ProductTable::new(&algebra);

        let scalar = 0;
        let e1 = 1; // non-degenerate
        let e4 = 8; // degenerate

        // scalar ⊛ e1 vanishes because the intermediate product involves e0
        let (sign, _result) = table.antiproduct(scalar, e1);
        assert_eq!(
            sign, 0,
            "scalar ⊛ e1 vanishes (intermediate uses degenerate e0)"
        );

        // scalar ⊛ e4 (e0) is NON-ZERO because e1234 × e123 uses only Euclidean bases
        let (sign, result) = table.antiproduct(scalar, e4);
        assert_ne!(sign, 0, "scalar ⊛ e4 should not vanish");
        assert_eq!(result, 7, "scalar ⊛ e4 should give e123");
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
// Complete 16×16 Cayley table tests for 3D PGA (G₃,₀,₁)
// Each test verifies all 256 entries of a product table.
// =============================================================================

#[cfg(test)]
mod full_cayley_tests {
    use crate::algebra::{Algebra, ProductTable};

    #[test]
    fn pga_3d_full_geometric_table() {
        let algebra = Algebra::pga(3);
        let table = ProductTable::new(&algebra);

        // Row 1
        assert_eq!(table.geometric(0, 0), (1, 0));
        assert_eq!(table.geometric(0, 1), (1, 1));
        assert_eq!(table.geometric(0, 2), (1, 2));
        assert_eq!(table.geometric(0, 3), (1, 3));
        assert_eq!(table.geometric(0, 4), (1, 4));
        assert_eq!(table.geometric(0, 5), (1, 5));
        assert_eq!(table.geometric(0, 6), (1, 6));
        assert_eq!(table.geometric(0, 7), (1, 7));
        assert_eq!(table.geometric(0, 8), (1, 8));
        assert_eq!(table.geometric(0, 9), (1, 9));
        assert_eq!(table.geometric(0, 10), (1, 10));
        assert_eq!(table.geometric(0, 11), (1, 11));
        assert_eq!(table.geometric(0, 12), (1, 12));
        assert_eq!(table.geometric(0, 13), (1, 13));
        assert_eq!(table.geometric(0, 14), (1, 14));
        assert_eq!(table.geometric(0, 15), (1, 15));
        // Row e1
        assert_eq!(table.geometric(1, 0), (1, 1));
        assert_eq!(table.geometric(1, 1), (1, 0));
        assert_eq!(table.geometric(1, 2), (1, 3));
        assert_eq!(table.geometric(1, 3), (1, 2));
        assert_eq!(table.geometric(1, 4), (1, 5));
        assert_eq!(table.geometric(1, 5), (1, 4));
        assert_eq!(table.geometric(1, 6), (1, 7));
        assert_eq!(table.geometric(1, 7), (1, 6));
        assert_eq!(table.geometric(1, 8), (1, 9));
        assert_eq!(table.geometric(1, 9), (1, 8));
        assert_eq!(table.geometric(1, 10), (1, 11));
        assert_eq!(table.geometric(1, 11), (1, 10));
        assert_eq!(table.geometric(1, 12), (1, 13));
        assert_eq!(table.geometric(1, 13), (1, 12));
        assert_eq!(table.geometric(1, 14), (1, 15));
        assert_eq!(table.geometric(1, 15), (1, 14));
        // Row e2
        assert_eq!(table.geometric(2, 0), (1, 2));
        assert_eq!(table.geometric(2, 1), (-1, 3));
        assert_eq!(table.geometric(2, 2), (1, 0));
        assert_eq!(table.geometric(2, 3), (-1, 1));
        assert_eq!(table.geometric(2, 4), (1, 6));
        assert_eq!(table.geometric(2, 5), (-1, 7));
        assert_eq!(table.geometric(2, 6), (1, 4));
        assert_eq!(table.geometric(2, 7), (-1, 5));
        assert_eq!(table.geometric(2, 8), (1, 10));
        assert_eq!(table.geometric(2, 9), (-1, 11));
        assert_eq!(table.geometric(2, 10), (1, 8));
        assert_eq!(table.geometric(2, 11), (-1, 9));
        assert_eq!(table.geometric(2, 12), (1, 14));
        assert_eq!(table.geometric(2, 13), (-1, 15));
        assert_eq!(table.geometric(2, 14), (1, 12));
        assert_eq!(table.geometric(2, 15), (-1, 13));
        // Row e12
        assert_eq!(table.geometric(3, 0), (1, 3));
        assert_eq!(table.geometric(3, 1), (-1, 2));
        assert_eq!(table.geometric(3, 2), (1, 1));
        assert_eq!(table.geometric(3, 3), (-1, 0));
        assert_eq!(table.geometric(3, 4), (1, 7));
        assert_eq!(table.geometric(3, 5), (-1, 6));
        assert_eq!(table.geometric(3, 6), (1, 5));
        assert_eq!(table.geometric(3, 7), (-1, 4));
        assert_eq!(table.geometric(3, 8), (1, 11));
        assert_eq!(table.geometric(3, 9), (-1, 10));
        assert_eq!(table.geometric(3, 10), (1, 9));
        assert_eq!(table.geometric(3, 11), (-1, 8));
        assert_eq!(table.geometric(3, 12), (1, 15));
        assert_eq!(table.geometric(3, 13), (-1, 14));
        assert_eq!(table.geometric(3, 14), (1, 13));
        assert_eq!(table.geometric(3, 15), (-1, 12));
        // Row e3
        assert_eq!(table.geometric(4, 0), (1, 4));
        assert_eq!(table.geometric(4, 1), (-1, 5));
        assert_eq!(table.geometric(4, 2), (-1, 6));
        assert_eq!(table.geometric(4, 3), (1, 7));
        assert_eq!(table.geometric(4, 4), (1, 0));
        assert_eq!(table.geometric(4, 5), (-1, 1));
        assert_eq!(table.geometric(4, 6), (-1, 2));
        assert_eq!(table.geometric(4, 7), (1, 3));
        assert_eq!(table.geometric(4, 8), (1, 12));
        assert_eq!(table.geometric(4, 9), (-1, 13));
        assert_eq!(table.geometric(4, 10), (-1, 14));
        assert_eq!(table.geometric(4, 11), (1, 15));
        assert_eq!(table.geometric(4, 12), (1, 8));
        assert_eq!(table.geometric(4, 13), (-1, 9));
        assert_eq!(table.geometric(4, 14), (-1, 10));
        assert_eq!(table.geometric(4, 15), (1, 11));
        // Row e13
        assert_eq!(table.geometric(5, 0), (1, 5));
        assert_eq!(table.geometric(5, 1), (-1, 4));
        assert_eq!(table.geometric(5, 2), (-1, 7));
        assert_eq!(table.geometric(5, 3), (1, 6));
        assert_eq!(table.geometric(5, 4), (1, 1));
        assert_eq!(table.geometric(5, 5), (-1, 0));
        assert_eq!(table.geometric(5, 6), (-1, 3));
        assert_eq!(table.geometric(5, 7), (1, 2));
        assert_eq!(table.geometric(5, 8), (1, 13));
        assert_eq!(table.geometric(5, 9), (-1, 12));
        assert_eq!(table.geometric(5, 10), (-1, 15));
        assert_eq!(table.geometric(5, 11), (1, 14));
        assert_eq!(table.geometric(5, 12), (1, 9));
        assert_eq!(table.geometric(5, 13), (-1, 8));
        assert_eq!(table.geometric(5, 14), (-1, 11));
        assert_eq!(table.geometric(5, 15), (1, 10));
        // Row e23
        assert_eq!(table.geometric(6, 0), (1, 6));
        assert_eq!(table.geometric(6, 1), (1, 7));
        assert_eq!(table.geometric(6, 2), (-1, 4));
        assert_eq!(table.geometric(6, 3), (-1, 5));
        assert_eq!(table.geometric(6, 4), (1, 2));
        assert_eq!(table.geometric(6, 5), (1, 3));
        assert_eq!(table.geometric(6, 6), (-1, 0));
        assert_eq!(table.geometric(6, 7), (-1, 1));
        assert_eq!(table.geometric(6, 8), (1, 14));
        assert_eq!(table.geometric(6, 9), (1, 15));
        assert_eq!(table.geometric(6, 10), (-1, 12));
        assert_eq!(table.geometric(6, 11), (-1, 13));
        assert_eq!(table.geometric(6, 12), (1, 10));
        assert_eq!(table.geometric(6, 13), (1, 11));
        assert_eq!(table.geometric(6, 14), (-1, 8));
        assert_eq!(table.geometric(6, 15), (-1, 9));
        // Row e123
        assert_eq!(table.geometric(7, 0), (1, 7));
        assert_eq!(table.geometric(7, 1), (1, 6));
        assert_eq!(table.geometric(7, 2), (-1, 5));
        assert_eq!(table.geometric(7, 3), (-1, 4));
        assert_eq!(table.geometric(7, 4), (1, 3));
        assert_eq!(table.geometric(7, 5), (1, 2));
        assert_eq!(table.geometric(7, 6), (-1, 1));
        assert_eq!(table.geometric(7, 7), (-1, 0));
        assert_eq!(table.geometric(7, 8), (1, 15));
        assert_eq!(table.geometric(7, 9), (1, 14));
        assert_eq!(table.geometric(7, 10), (-1, 13));
        assert_eq!(table.geometric(7, 11), (-1, 12));
        assert_eq!(table.geometric(7, 12), (1, 11));
        assert_eq!(table.geometric(7, 13), (1, 10));
        assert_eq!(table.geometric(7, 14), (-1, 9));
        assert_eq!(table.geometric(7, 15), (-1, 8));
        // Row e4
        assert_eq!(table.geometric(8, 0), (1, 8));
        assert_eq!(table.geometric(8, 1), (-1, 9));
        assert_eq!(table.geometric(8, 2), (-1, 10));
        assert_eq!(table.geometric(8, 3), (1, 11));
        assert_eq!(table.geometric(8, 4), (-1, 12));
        assert_eq!(table.geometric(8, 5), (1, 13));
        assert_eq!(table.geometric(8, 6), (1, 14));
        assert_eq!(table.geometric(8, 7), (-1, 15));
        assert_eq!(table.geometric(8, 8).0, 0);
        assert_eq!(table.geometric(8, 9).0, 0);
        assert_eq!(table.geometric(8, 10).0, 0);
        assert_eq!(table.geometric(8, 11).0, 0);
        assert_eq!(table.geometric(8, 12).0, 0);
        assert_eq!(table.geometric(8, 13).0, 0);
        assert_eq!(table.geometric(8, 14).0, 0);
        assert_eq!(table.geometric(8, 15).0, 0);
        // Row e14
        assert_eq!(table.geometric(9, 0), (1, 9));
        assert_eq!(table.geometric(9, 1), (-1, 8));
        assert_eq!(table.geometric(9, 2), (-1, 11));
        assert_eq!(table.geometric(9, 3), (1, 10));
        assert_eq!(table.geometric(9, 4), (-1, 13));
        assert_eq!(table.geometric(9, 5), (1, 12));
        assert_eq!(table.geometric(9, 6), (1, 15));
        assert_eq!(table.geometric(9, 7), (-1, 14));
        assert_eq!(table.geometric(9, 8).0, 0);
        assert_eq!(table.geometric(9, 9).0, 0);
        assert_eq!(table.geometric(9, 10).0, 0);
        assert_eq!(table.geometric(9, 11).0, 0);
        assert_eq!(table.geometric(9, 12).0, 0);
        assert_eq!(table.geometric(9, 13).0, 0);
        assert_eq!(table.geometric(9, 14).0, 0);
        assert_eq!(table.geometric(9, 15).0, 0);
        // Row e24
        assert_eq!(table.geometric(10, 0), (1, 10));
        assert_eq!(table.geometric(10, 1), (1, 11));
        assert_eq!(table.geometric(10, 2), (-1, 8));
        assert_eq!(table.geometric(10, 3), (-1, 9));
        assert_eq!(table.geometric(10, 4), (-1, 14));
        assert_eq!(table.geometric(10, 5), (-1, 15));
        assert_eq!(table.geometric(10, 6), (1, 12));
        assert_eq!(table.geometric(10, 7), (1, 13));
        assert_eq!(table.geometric(10, 8).0, 0);
        assert_eq!(table.geometric(10, 9).0, 0);
        assert_eq!(table.geometric(10, 10).0, 0);
        assert_eq!(table.geometric(10, 11).0, 0);
        assert_eq!(table.geometric(10, 12).0, 0);
        assert_eq!(table.geometric(10, 13).0, 0);
        assert_eq!(table.geometric(10, 14).0, 0);
        assert_eq!(table.geometric(10, 15).0, 0);
        // Row e124
        assert_eq!(table.geometric(11, 0), (1, 11));
        assert_eq!(table.geometric(11, 1), (1, 10));
        assert_eq!(table.geometric(11, 2), (-1, 9));
        assert_eq!(table.geometric(11, 3), (-1, 8));
        assert_eq!(table.geometric(11, 4), (-1, 15));
        assert_eq!(table.geometric(11, 5), (-1, 14));
        assert_eq!(table.geometric(11, 6), (1, 13));
        assert_eq!(table.geometric(11, 7), (1, 12));
        assert_eq!(table.geometric(11, 8).0, 0);
        assert_eq!(table.geometric(11, 9).0, 0);
        assert_eq!(table.geometric(11, 10).0, 0);
        assert_eq!(table.geometric(11, 11).0, 0);
        assert_eq!(table.geometric(11, 12).0, 0);
        assert_eq!(table.geometric(11, 13).0, 0);
        assert_eq!(table.geometric(11, 14).0, 0);
        assert_eq!(table.geometric(11, 15).0, 0);
        // Row e34
        assert_eq!(table.geometric(12, 0), (1, 12));
        assert_eq!(table.geometric(12, 1), (1, 13));
        assert_eq!(table.geometric(12, 2), (1, 14));
        assert_eq!(table.geometric(12, 3), (1, 15));
        assert_eq!(table.geometric(12, 4), (-1, 8));
        assert_eq!(table.geometric(12, 5), (-1, 9));
        assert_eq!(table.geometric(12, 6), (-1, 10));
        assert_eq!(table.geometric(12, 7), (-1, 11));
        assert_eq!(table.geometric(12, 8).0, 0);
        assert_eq!(table.geometric(12, 9).0, 0);
        assert_eq!(table.geometric(12, 10).0, 0);
        assert_eq!(table.geometric(12, 11).0, 0);
        assert_eq!(table.geometric(12, 12).0, 0);
        assert_eq!(table.geometric(12, 13).0, 0);
        assert_eq!(table.geometric(12, 14).0, 0);
        assert_eq!(table.geometric(12, 15).0, 0);
        // Row e134
        assert_eq!(table.geometric(13, 0), (1, 13));
        assert_eq!(table.geometric(13, 1), (1, 12));
        assert_eq!(table.geometric(13, 2), (1, 15));
        assert_eq!(table.geometric(13, 3), (1, 14));
        assert_eq!(table.geometric(13, 4), (-1, 9));
        assert_eq!(table.geometric(13, 5), (-1, 8));
        assert_eq!(table.geometric(13, 6), (-1, 11));
        assert_eq!(table.geometric(13, 7), (-1, 10));
        assert_eq!(table.geometric(13, 8).0, 0);
        assert_eq!(table.geometric(13, 9).0, 0);
        assert_eq!(table.geometric(13, 10).0, 0);
        assert_eq!(table.geometric(13, 11).0, 0);
        assert_eq!(table.geometric(13, 12).0, 0);
        assert_eq!(table.geometric(13, 13).0, 0);
        assert_eq!(table.geometric(13, 14).0, 0);
        assert_eq!(table.geometric(13, 15).0, 0);
        // Row e234
        assert_eq!(table.geometric(14, 0), (1, 14));
        assert_eq!(table.geometric(14, 1), (-1, 15));
        assert_eq!(table.geometric(14, 2), (1, 12));
        assert_eq!(table.geometric(14, 3), (-1, 13));
        assert_eq!(table.geometric(14, 4), (-1, 10));
        assert_eq!(table.geometric(14, 5), (1, 11));
        assert_eq!(table.geometric(14, 6), (-1, 8));
        assert_eq!(table.geometric(14, 7), (1, 9));
        assert_eq!(table.geometric(14, 8).0, 0);
        assert_eq!(table.geometric(14, 9).0, 0);
        assert_eq!(table.geometric(14, 10).0, 0);
        assert_eq!(table.geometric(14, 11).0, 0);
        assert_eq!(table.geometric(14, 12).0, 0);
        assert_eq!(table.geometric(14, 13).0, 0);
        assert_eq!(table.geometric(14, 14).0, 0);
        assert_eq!(table.geometric(14, 15).0, 0);
        // Row e1234
        assert_eq!(table.geometric(15, 0), (1, 15));
        assert_eq!(table.geometric(15, 1), (-1, 14));
        assert_eq!(table.geometric(15, 2), (1, 13));
        assert_eq!(table.geometric(15, 3), (-1, 12));
        assert_eq!(table.geometric(15, 4), (-1, 11));
        assert_eq!(table.geometric(15, 5), (1, 10));
        assert_eq!(table.geometric(15, 6), (-1, 9));
        assert_eq!(table.geometric(15, 7), (1, 8));
        assert_eq!(table.geometric(15, 8).0, 0);
        assert_eq!(table.geometric(15, 9).0, 0);
        assert_eq!(table.geometric(15, 10).0, 0);
        assert_eq!(table.geometric(15, 11).0, 0);
        assert_eq!(table.geometric(15, 12).0, 0);
        assert_eq!(table.geometric(15, 13).0, 0);
        assert_eq!(table.geometric(15, 14).0, 0);
        assert_eq!(table.geometric(15, 15).0, 0);
    }

    #[test]
    fn pga_3d_full_exterior_table() {
        let algebra = Algebra::pga(3);
        let table = ProductTable::new(&algebra);

        // Row 1
        assert_eq!(table.exterior(0, 0), (1, 0));
        assert_eq!(table.exterior(0, 1), (1, 1));
        assert_eq!(table.exterior(0, 2), (1, 2));
        assert_eq!(table.exterior(0, 3), (1, 3));
        assert_eq!(table.exterior(0, 4), (1, 4));
        assert_eq!(table.exterior(0, 5), (1, 5));
        assert_eq!(table.exterior(0, 6), (1, 6));
        assert_eq!(table.exterior(0, 7), (1, 7));
        assert_eq!(table.exterior(0, 8), (1, 8));
        assert_eq!(table.exterior(0, 9), (1, 9));
        assert_eq!(table.exterior(0, 10), (1, 10));
        assert_eq!(table.exterior(0, 11), (1, 11));
        assert_eq!(table.exterior(0, 12), (1, 12));
        assert_eq!(table.exterior(0, 13), (1, 13));
        assert_eq!(table.exterior(0, 14), (1, 14));
        assert_eq!(table.exterior(0, 15), (1, 15));
        // Row e1
        assert_eq!(table.exterior(1, 0), (1, 1));
        assert_eq!(table.exterior(1, 1), (0, 0));
        assert_eq!(table.exterior(1, 2), (1, 3));
        assert_eq!(table.exterior(1, 3), (0, 0));
        assert_eq!(table.exterior(1, 4), (1, 5));
        assert_eq!(table.exterior(1, 5), (0, 0));
        assert_eq!(table.exterior(1, 6), (1, 7));
        assert_eq!(table.exterior(1, 7), (0, 0));
        assert_eq!(table.exterior(1, 8), (1, 9));
        assert_eq!(table.exterior(1, 9), (0, 0));
        assert_eq!(table.exterior(1, 10), (1, 11));
        assert_eq!(table.exterior(1, 11), (0, 0));
        assert_eq!(table.exterior(1, 12), (1, 13));
        assert_eq!(table.exterior(1, 13), (0, 0));
        assert_eq!(table.exterior(1, 14), (1, 15));
        assert_eq!(table.exterior(1, 15), (0, 0));
        // Row e2
        assert_eq!(table.exterior(2, 0), (1, 2));
        assert_eq!(table.exterior(2, 1), (-1, 3));
        assert_eq!(table.exterior(2, 2), (0, 0));
        assert_eq!(table.exterior(2, 3), (0, 0));
        assert_eq!(table.exterior(2, 4), (1, 6));
        assert_eq!(table.exterior(2, 5), (-1, 7));
        assert_eq!(table.exterior(2, 6), (0, 0));
        assert_eq!(table.exterior(2, 7), (0, 0));
        assert_eq!(table.exterior(2, 8), (1, 10));
        assert_eq!(table.exterior(2, 9), (-1, 11));
        assert_eq!(table.exterior(2, 10), (0, 0));
        assert_eq!(table.exterior(2, 11), (0, 0));
        assert_eq!(table.exterior(2, 12), (1, 14));
        assert_eq!(table.exterior(2, 13), (-1, 15));
        assert_eq!(table.exterior(2, 14), (0, 0));
        assert_eq!(table.exterior(2, 15), (0, 0));
        // Row e12
        assert_eq!(table.exterior(3, 0), (1, 3));
        assert_eq!(table.exterior(3, 1), (0, 0));
        assert_eq!(table.exterior(3, 2), (0, 0));
        assert_eq!(table.exterior(3, 3), (0, 0));
        assert_eq!(table.exterior(3, 4), (1, 7));
        assert_eq!(table.exterior(3, 5), (0, 0));
        assert_eq!(table.exterior(3, 6), (0, 0));
        assert_eq!(table.exterior(3, 7), (0, 0));
        assert_eq!(table.exterior(3, 8), (1, 11));
        assert_eq!(table.exterior(3, 9), (0, 0));
        assert_eq!(table.exterior(3, 10), (0, 0));
        assert_eq!(table.exterior(3, 11), (0, 0));
        assert_eq!(table.exterior(3, 12), (1, 15));
        assert_eq!(table.exterior(3, 13), (0, 0));
        assert_eq!(table.exterior(3, 14), (0, 0));
        assert_eq!(table.exterior(3, 15), (0, 0));
        // Row e3
        assert_eq!(table.exterior(4, 0), (1, 4));
        assert_eq!(table.exterior(4, 1), (-1, 5));
        assert_eq!(table.exterior(4, 2), (-1, 6));
        assert_eq!(table.exterior(4, 3), (1, 7));
        assert_eq!(table.exterior(4, 4), (0, 0));
        assert_eq!(table.exterior(4, 5), (0, 0));
        assert_eq!(table.exterior(4, 6), (0, 0));
        assert_eq!(table.exterior(4, 7), (0, 0));
        assert_eq!(table.exterior(4, 8), (1, 12));
        assert_eq!(table.exterior(4, 9), (-1, 13));
        assert_eq!(table.exterior(4, 10), (-1, 14));
        assert_eq!(table.exterior(4, 11), (1, 15));
        assert_eq!(table.exterior(4, 12), (0, 0));
        assert_eq!(table.exterior(4, 13), (0, 0));
        assert_eq!(table.exterior(4, 14), (0, 0));
        assert_eq!(table.exterior(4, 15), (0, 0));
        // Row e13
        assert_eq!(table.exterior(5, 0), (1, 5));
        assert_eq!(table.exterior(5, 1), (0, 0));
        assert_eq!(table.exterior(5, 2), (-1, 7));
        assert_eq!(table.exterior(5, 3), (0, 0));
        assert_eq!(table.exterior(5, 4), (0, 0));
        assert_eq!(table.exterior(5, 5), (0, 0));
        assert_eq!(table.exterior(5, 6), (0, 0));
        assert_eq!(table.exterior(5, 7), (0, 0));
        assert_eq!(table.exterior(5, 8), (1, 13));
        assert_eq!(table.exterior(5, 9), (0, 0));
        assert_eq!(table.exterior(5, 10), (-1, 15));
        assert_eq!(table.exterior(5, 11), (0, 0));
        assert_eq!(table.exterior(5, 12), (0, 0));
        assert_eq!(table.exterior(5, 13), (0, 0));
        assert_eq!(table.exterior(5, 14), (0, 0));
        assert_eq!(table.exterior(5, 15), (0, 0));
        // Row e23
        assert_eq!(table.exterior(6, 0), (1, 6));
        assert_eq!(table.exterior(6, 1), (1, 7));
        assert_eq!(table.exterior(6, 2), (0, 0));
        assert_eq!(table.exterior(6, 3), (0, 0));
        assert_eq!(table.exterior(6, 4), (0, 0));
        assert_eq!(table.exterior(6, 5), (0, 0));
        assert_eq!(table.exterior(6, 6), (0, 0));
        assert_eq!(table.exterior(6, 7), (0, 0));
        assert_eq!(table.exterior(6, 8), (1, 14));
        assert_eq!(table.exterior(6, 9), (1, 15));
        assert_eq!(table.exterior(6, 10), (0, 0));
        assert_eq!(table.exterior(6, 11), (0, 0));
        assert_eq!(table.exterior(6, 12), (0, 0));
        assert_eq!(table.exterior(6, 13), (0, 0));
        assert_eq!(table.exterior(6, 14), (0, 0));
        assert_eq!(table.exterior(6, 15), (0, 0));
        // Row e123
        assert_eq!(table.exterior(7, 0), (1, 7));
        assert_eq!(table.exterior(7, 1), (0, 0));
        assert_eq!(table.exterior(7, 2), (0, 0));
        assert_eq!(table.exterior(7, 3), (0, 0));
        assert_eq!(table.exterior(7, 4), (0, 0));
        assert_eq!(table.exterior(7, 5), (0, 0));
        assert_eq!(table.exterior(7, 6), (0, 0));
        assert_eq!(table.exterior(7, 7), (0, 0));
        assert_eq!(table.exterior(7, 8), (1, 15));
        assert_eq!(table.exterior(7, 9), (0, 0));
        assert_eq!(table.exterior(7, 10), (0, 0));
        assert_eq!(table.exterior(7, 11), (0, 0));
        assert_eq!(table.exterior(7, 12), (0, 0));
        assert_eq!(table.exterior(7, 13), (0, 0));
        assert_eq!(table.exterior(7, 14), (0, 0));
        assert_eq!(table.exterior(7, 15), (0, 0));
        // Row e4
        assert_eq!(table.exterior(8, 0), (1, 8));
        assert_eq!(table.exterior(8, 1), (-1, 9));
        assert_eq!(table.exterior(8, 2), (-1, 10));
        assert_eq!(table.exterior(8, 3), (1, 11));
        assert_eq!(table.exterior(8, 4), (-1, 12));
        assert_eq!(table.exterior(8, 5), (1, 13));
        assert_eq!(table.exterior(8, 6), (1, 14));
        assert_eq!(table.exterior(8, 7), (-1, 15));
        assert_eq!(table.exterior(8, 8), (0, 0));
        assert_eq!(table.exterior(8, 9), (0, 0));
        assert_eq!(table.exterior(8, 10), (0, 0));
        assert_eq!(table.exterior(8, 11), (0, 0));
        assert_eq!(table.exterior(8, 12), (0, 0));
        assert_eq!(table.exterior(8, 13), (0, 0));
        assert_eq!(table.exterior(8, 14), (0, 0));
        assert_eq!(table.exterior(8, 15), (0, 0));
        // Row e14
        assert_eq!(table.exterior(9, 0), (1, 9));
        assert_eq!(table.exterior(9, 1), (0, 0));
        assert_eq!(table.exterior(9, 2), (-1, 11));
        assert_eq!(table.exterior(9, 3), (0, 0));
        assert_eq!(table.exterior(9, 4), (-1, 13));
        assert_eq!(table.exterior(9, 5), (0, 0));
        assert_eq!(table.exterior(9, 6), (1, 15));
        assert_eq!(table.exterior(9, 7), (0, 0));
        assert_eq!(table.exterior(9, 8), (0, 0));
        assert_eq!(table.exterior(9, 9), (0, 0));
        assert_eq!(table.exterior(9, 10), (0, 0));
        assert_eq!(table.exterior(9, 11), (0, 0));
        assert_eq!(table.exterior(9, 12), (0, 0));
        assert_eq!(table.exterior(9, 13), (0, 0));
        assert_eq!(table.exterior(9, 14), (0, 0));
        assert_eq!(table.exterior(9, 15), (0, 0));
        // Row e24
        assert_eq!(table.exterior(10, 0), (1, 10));
        assert_eq!(table.exterior(10, 1), (1, 11));
        assert_eq!(table.exterior(10, 2), (0, 0));
        assert_eq!(table.exterior(10, 3), (0, 0));
        assert_eq!(table.exterior(10, 4), (-1, 14));
        assert_eq!(table.exterior(10, 5), (-1, 15));
        assert_eq!(table.exterior(10, 6), (0, 0));
        assert_eq!(table.exterior(10, 7), (0, 0));
        assert_eq!(table.exterior(10, 8), (0, 0));
        assert_eq!(table.exterior(10, 9), (0, 0));
        assert_eq!(table.exterior(10, 10), (0, 0));
        assert_eq!(table.exterior(10, 11), (0, 0));
        assert_eq!(table.exterior(10, 12), (0, 0));
        assert_eq!(table.exterior(10, 13), (0, 0));
        assert_eq!(table.exterior(10, 14), (0, 0));
        assert_eq!(table.exterior(10, 15), (0, 0));
        // Row e124
        assert_eq!(table.exterior(11, 0), (1, 11));
        assert_eq!(table.exterior(11, 1), (0, 0));
        assert_eq!(table.exterior(11, 2), (0, 0));
        assert_eq!(table.exterior(11, 3), (0, 0));
        assert_eq!(table.exterior(11, 4), (-1, 15));
        assert_eq!(table.exterior(11, 5), (0, 0));
        assert_eq!(table.exterior(11, 6), (0, 0));
        assert_eq!(table.exterior(11, 7), (0, 0));
        assert_eq!(table.exterior(11, 8), (0, 0));
        assert_eq!(table.exterior(11, 9), (0, 0));
        assert_eq!(table.exterior(11, 10), (0, 0));
        assert_eq!(table.exterior(11, 11), (0, 0));
        assert_eq!(table.exterior(11, 12), (0, 0));
        assert_eq!(table.exterior(11, 13), (0, 0));
        assert_eq!(table.exterior(11, 14), (0, 0));
        assert_eq!(table.exterior(11, 15), (0, 0));
        // Row e34
        assert_eq!(table.exterior(12, 0), (1, 12));
        assert_eq!(table.exterior(12, 1), (1, 13));
        assert_eq!(table.exterior(12, 2), (1, 14));
        assert_eq!(table.exterior(12, 3), (1, 15));
        assert_eq!(table.exterior(12, 4), (0, 0));
        assert_eq!(table.exterior(12, 5), (0, 0));
        assert_eq!(table.exterior(12, 6), (0, 0));
        assert_eq!(table.exterior(12, 7), (0, 0));
        assert_eq!(table.exterior(12, 8), (0, 0));
        assert_eq!(table.exterior(12, 9), (0, 0));
        assert_eq!(table.exterior(12, 10), (0, 0));
        assert_eq!(table.exterior(12, 11), (0, 0));
        assert_eq!(table.exterior(12, 12), (0, 0));
        assert_eq!(table.exterior(12, 13), (0, 0));
        assert_eq!(table.exterior(12, 14), (0, 0));
        assert_eq!(table.exterior(12, 15), (0, 0));
        // Row e134
        assert_eq!(table.exterior(13, 0), (1, 13));
        assert_eq!(table.exterior(13, 1), (0, 0));
        assert_eq!(table.exterior(13, 2), (1, 15));
        assert_eq!(table.exterior(13, 3), (0, 0));
        assert_eq!(table.exterior(13, 4), (0, 0));
        assert_eq!(table.exterior(13, 5), (0, 0));
        assert_eq!(table.exterior(13, 6), (0, 0));
        assert_eq!(table.exterior(13, 7), (0, 0));
        assert_eq!(table.exterior(13, 8), (0, 0));
        assert_eq!(table.exterior(13, 9), (0, 0));
        assert_eq!(table.exterior(13, 10), (0, 0));
        assert_eq!(table.exterior(13, 11), (0, 0));
        assert_eq!(table.exterior(13, 12), (0, 0));
        assert_eq!(table.exterior(13, 13), (0, 0));
        assert_eq!(table.exterior(13, 14), (0, 0));
        assert_eq!(table.exterior(13, 15), (0, 0));
        // Row e234
        assert_eq!(table.exterior(14, 0), (1, 14));
        assert_eq!(table.exterior(14, 1), (-1, 15));
        assert_eq!(table.exterior(14, 2), (0, 0));
        assert_eq!(table.exterior(14, 3), (0, 0));
        assert_eq!(table.exterior(14, 4), (0, 0));
        assert_eq!(table.exterior(14, 5), (0, 0));
        assert_eq!(table.exterior(14, 6), (0, 0));
        assert_eq!(table.exterior(14, 7), (0, 0));
        assert_eq!(table.exterior(14, 8), (0, 0));
        assert_eq!(table.exterior(14, 9), (0, 0));
        assert_eq!(table.exterior(14, 10), (0, 0));
        assert_eq!(table.exterior(14, 11), (0, 0));
        assert_eq!(table.exterior(14, 12), (0, 0));
        assert_eq!(table.exterior(14, 13), (0, 0));
        assert_eq!(table.exterior(14, 14), (0, 0));
        assert_eq!(table.exterior(14, 15), (0, 0));
        // Row e1234
        assert_eq!(table.exterior(15, 0), (1, 15));
        assert_eq!(table.exterior(15, 1), (0, 0));
        assert_eq!(table.exterior(15, 2), (0, 0));
        assert_eq!(table.exterior(15, 3), (0, 0));
        assert_eq!(table.exterior(15, 4), (0, 0));
        assert_eq!(table.exterior(15, 5), (0, 0));
        assert_eq!(table.exterior(15, 6), (0, 0));
        assert_eq!(table.exterior(15, 7), (0, 0));
        assert_eq!(table.exterior(15, 8), (0, 0));
        assert_eq!(table.exterior(15, 9), (0, 0));
        assert_eq!(table.exterior(15, 10), (0, 0));
        assert_eq!(table.exterior(15, 11), (0, 0));
        assert_eq!(table.exterior(15, 12), (0, 0));
        assert_eq!(table.exterior(15, 13), (0, 0));
        assert_eq!(table.exterior(15, 14), (0, 0));
        assert_eq!(table.exterior(15, 15), (0, 0));
    }

    #[test]
    fn pga_3d_full_regressive_table() {
        let algebra = Algebra::pga(3);
        let table = ProductTable::new(&algebra);

        // Row 1
        assert_eq!(table.regressive(0, 0), (0, 0));
        assert_eq!(table.regressive(0, 1), (0, 0));
        assert_eq!(table.regressive(0, 2), (0, 0));
        assert_eq!(table.regressive(0, 3), (0, 0));
        assert_eq!(table.regressive(0, 4), (0, 0));
        assert_eq!(table.regressive(0, 5), (0, 0));
        assert_eq!(table.regressive(0, 6), (0, 0));
        assert_eq!(table.regressive(0, 7), (0, 0));
        assert_eq!(table.regressive(0, 8), (0, 0));
        assert_eq!(table.regressive(0, 9), (0, 0));
        assert_eq!(table.regressive(0, 10), (0, 0));
        assert_eq!(table.regressive(0, 11), (0, 0));
        assert_eq!(table.regressive(0, 12), (0, 0));
        assert_eq!(table.regressive(0, 13), (0, 0));
        assert_eq!(table.regressive(0, 14), (0, 0));
        assert_eq!(table.regressive(0, 15), (1, 0));
        // Row e1
        assert_eq!(table.regressive(1, 0), (0, 0));
        assert_eq!(table.regressive(1, 1), (0, 0));
        assert_eq!(table.regressive(1, 2), (0, 0));
        assert_eq!(table.regressive(1, 3), (0, 0));
        assert_eq!(table.regressive(1, 4), (0, 0));
        assert_eq!(table.regressive(1, 5), (0, 0));
        assert_eq!(table.regressive(1, 6), (0, 0));
        assert_eq!(table.regressive(1, 7), (0, 0));
        assert_eq!(table.regressive(1, 8), (0, 0));
        assert_eq!(table.regressive(1, 9), (0, 0));
        assert_eq!(table.regressive(1, 10), (0, 0));
        assert_eq!(table.regressive(1, 11), (0, 0));
        assert_eq!(table.regressive(1, 12), (0, 0));
        assert_eq!(table.regressive(1, 13), (0, 0));
        assert_eq!(table.regressive(1, 14), (1, 0));
        assert_eq!(table.regressive(1, 15), (-1, 1));
        // Row e2
        assert_eq!(table.regressive(2, 0), (0, 0));
        assert_eq!(table.regressive(2, 1), (0, 0));
        assert_eq!(table.regressive(2, 2), (0, 0));
        assert_eq!(table.regressive(2, 3), (0, 0));
        assert_eq!(table.regressive(2, 4), (0, 0));
        assert_eq!(table.regressive(2, 5), (0, 0));
        assert_eq!(table.regressive(2, 6), (0, 0));
        assert_eq!(table.regressive(2, 7), (0, 0));
        assert_eq!(table.regressive(2, 8), (0, 0));
        assert_eq!(table.regressive(2, 9), (0, 0));
        assert_eq!(table.regressive(2, 10), (0, 0));
        assert_eq!(table.regressive(2, 11), (0, 0));
        assert_eq!(table.regressive(2, 12), (0, 0));
        assert_eq!(table.regressive(2, 13), (-1, 0));
        assert_eq!(table.regressive(2, 14), (0, 0));
        assert_eq!(table.regressive(2, 15), (-1, 2));
        // Row e12
        assert_eq!(table.regressive(3, 0), (0, 0));
        assert_eq!(table.regressive(3, 1), (0, 0));
        assert_eq!(table.regressive(3, 2), (0, 0));
        assert_eq!(table.regressive(3, 3), (0, 0));
        assert_eq!(table.regressive(3, 4), (0, 0));
        assert_eq!(table.regressive(3, 5), (0, 0));
        assert_eq!(table.regressive(3, 6), (0, 0));
        assert_eq!(table.regressive(3, 7), (0, 0));
        assert_eq!(table.regressive(3, 8), (0, 0));
        assert_eq!(table.regressive(3, 9), (0, 0));
        assert_eq!(table.regressive(3, 10), (0, 0));
        assert_eq!(table.regressive(3, 11), (0, 0));
        assert_eq!(table.regressive(3, 12), (1, 0));
        assert_eq!(table.regressive(3, 13), (-1, 1));
        assert_eq!(table.regressive(3, 14), (-1, 2));
        assert_eq!(table.regressive(3, 15), (1, 3));
        // Row e3
        assert_eq!(table.regressive(4, 0), (0, 0));
        assert_eq!(table.regressive(4, 1), (0, 0));
        assert_eq!(table.regressive(4, 2), (0, 0));
        assert_eq!(table.regressive(4, 3), (0, 0));
        assert_eq!(table.regressive(4, 4), (0, 0));
        assert_eq!(table.regressive(4, 5), (0, 0));
        assert_eq!(table.regressive(4, 6), (0, 0));
        assert_eq!(table.regressive(4, 7), (0, 0));
        assert_eq!(table.regressive(4, 8), (0, 0));
        assert_eq!(table.regressive(4, 9), (0, 0));
        assert_eq!(table.regressive(4, 10), (0, 0));
        assert_eq!(table.regressive(4, 11), (1, 0));
        assert_eq!(table.regressive(4, 12), (0, 0));
        assert_eq!(table.regressive(4, 13), (0, 0));
        assert_eq!(table.regressive(4, 14), (0, 0));
        assert_eq!(table.regressive(4, 15), (-1, 4));
        // Row e13
        assert_eq!(table.regressive(5, 0), (0, 0));
        assert_eq!(table.regressive(5, 1), (0, 0));
        assert_eq!(table.regressive(5, 2), (0, 0));
        assert_eq!(table.regressive(5, 3), (0, 0));
        assert_eq!(table.regressive(5, 4), (0, 0));
        assert_eq!(table.regressive(5, 5), (0, 0));
        assert_eq!(table.regressive(5, 6), (0, 0));
        assert_eq!(table.regressive(5, 7), (0, 0));
        assert_eq!(table.regressive(5, 8), (0, 0));
        assert_eq!(table.regressive(5, 9), (0, 0));
        assert_eq!(table.regressive(5, 10), (-1, 0));
        assert_eq!(table.regressive(5, 11), (1, 1));
        assert_eq!(table.regressive(5, 12), (0, 0));
        assert_eq!(table.regressive(5, 13), (0, 0));
        assert_eq!(table.regressive(5, 14), (-1, 4));
        assert_eq!(table.regressive(5, 15), (1, 5));
        // Row e23
        assert_eq!(table.regressive(6, 0), (0, 0));
        assert_eq!(table.regressive(6, 1), (0, 0));
        assert_eq!(table.regressive(6, 2), (0, 0));
        assert_eq!(table.regressive(6, 3), (0, 0));
        assert_eq!(table.regressive(6, 4), (0, 0));
        assert_eq!(table.regressive(6, 5), (0, 0));
        assert_eq!(table.regressive(6, 6), (0, 0));
        assert_eq!(table.regressive(6, 7), (0, 0));
        assert_eq!(table.regressive(6, 8), (0, 0));
        assert_eq!(table.regressive(6, 9), (1, 0));
        assert_eq!(table.regressive(6, 10), (0, 0));
        assert_eq!(table.regressive(6, 11), (1, 2));
        assert_eq!(table.regressive(6, 12), (0, 0));
        assert_eq!(table.regressive(6, 13), (1, 4));
        assert_eq!(table.regressive(6, 14), (0, 0));
        assert_eq!(table.regressive(6, 15), (1, 6));
        // Row e123
        assert_eq!(table.regressive(7, 0), (0, 0));
        assert_eq!(table.regressive(7, 1), (0, 0));
        assert_eq!(table.regressive(7, 2), (0, 0));
        assert_eq!(table.regressive(7, 3), (0, 0));
        assert_eq!(table.regressive(7, 4), (0, 0));
        assert_eq!(table.regressive(7, 5), (0, 0));
        assert_eq!(table.regressive(7, 6), (0, 0));
        assert_eq!(table.regressive(7, 7), (0, 0));
        assert_eq!(table.regressive(7, 8), (1, 0));
        assert_eq!(table.regressive(7, 9), (-1, 1));
        assert_eq!(table.regressive(7, 10), (-1, 2));
        assert_eq!(table.regressive(7, 11), (1, 3));
        assert_eq!(table.regressive(7, 12), (-1, 4));
        assert_eq!(table.regressive(7, 13), (1, 5));
        assert_eq!(table.regressive(7, 14), (1, 6));
        assert_eq!(table.regressive(7, 15), (-1, 7));
        // Row e4
        assert_eq!(table.regressive(8, 0), (0, 0));
        assert_eq!(table.regressive(8, 1), (0, 0));
        assert_eq!(table.regressive(8, 2), (0, 0));
        assert_eq!(table.regressive(8, 3), (0, 0));
        assert_eq!(table.regressive(8, 4), (0, 0));
        assert_eq!(table.regressive(8, 5), (0, 0));
        assert_eq!(table.regressive(8, 6), (0, 0));
        assert_eq!(table.regressive(8, 7), (-1, 0));
        assert_eq!(table.regressive(8, 8), (0, 0));
        assert_eq!(table.regressive(8, 9), (0, 0));
        assert_eq!(table.regressive(8, 10), (0, 0));
        assert_eq!(table.regressive(8, 11), (0, 0));
        assert_eq!(table.regressive(8, 12), (0, 0));
        assert_eq!(table.regressive(8, 13), (0, 0));
        assert_eq!(table.regressive(8, 14), (0, 0));
        assert_eq!(table.regressive(8, 15), (-1, 8));
        // Row e14
        assert_eq!(table.regressive(9, 0), (0, 0));
        assert_eq!(table.regressive(9, 1), (0, 0));
        assert_eq!(table.regressive(9, 2), (0, 0));
        assert_eq!(table.regressive(9, 3), (0, 0));
        assert_eq!(table.regressive(9, 4), (0, 0));
        assert_eq!(table.regressive(9, 5), (0, 0));
        assert_eq!(table.regressive(9, 6), (1, 0));
        assert_eq!(table.regressive(9, 7), (-1, 1));
        assert_eq!(table.regressive(9, 8), (0, 0));
        assert_eq!(table.regressive(9, 9), (0, 0));
        assert_eq!(table.regressive(9, 10), (0, 0));
        assert_eq!(table.regressive(9, 11), (0, 0));
        assert_eq!(table.regressive(9, 12), (0, 0));
        assert_eq!(table.regressive(9, 13), (0, 0));
        assert_eq!(table.regressive(9, 14), (-1, 8));
        assert_eq!(table.regressive(9, 15), (1, 9));
        // Row e24
        assert_eq!(table.regressive(10, 0), (0, 0));
        assert_eq!(table.regressive(10, 1), (0, 0));
        assert_eq!(table.regressive(10, 2), (0, 0));
        assert_eq!(table.regressive(10, 3), (0, 0));
        assert_eq!(table.regressive(10, 4), (0, 0));
        assert_eq!(table.regressive(10, 5), (-1, 0));
        assert_eq!(table.regressive(10, 6), (0, 0));
        assert_eq!(table.regressive(10, 7), (-1, 2));
        assert_eq!(table.regressive(10, 8), (0, 0));
        assert_eq!(table.regressive(10, 9), (0, 0));
        assert_eq!(table.regressive(10, 10), (0, 0));
        assert_eq!(table.regressive(10, 11), (0, 0));
        assert_eq!(table.regressive(10, 12), (0, 0));
        assert_eq!(table.regressive(10, 13), (1, 8));
        assert_eq!(table.regressive(10, 14), (0, 0));
        assert_eq!(table.regressive(10, 15), (1, 10));
        // Row e124
        assert_eq!(table.regressive(11, 0), (0, 0));
        assert_eq!(table.regressive(11, 1), (0, 0));
        assert_eq!(table.regressive(11, 2), (0, 0));
        assert_eq!(table.regressive(11, 3), (0, 0));
        assert_eq!(table.regressive(11, 4), (-1, 0));
        assert_eq!(table.regressive(11, 5), (1, 1));
        assert_eq!(table.regressive(11, 6), (1, 2));
        assert_eq!(table.regressive(11, 7), (-1, 3));
        assert_eq!(table.regressive(11, 8), (0, 0));
        assert_eq!(table.regressive(11, 9), (0, 0));
        assert_eq!(table.regressive(11, 10), (0, 0));
        assert_eq!(table.regressive(11, 11), (0, 0));
        assert_eq!(table.regressive(11, 12), (-1, 8));
        assert_eq!(table.regressive(11, 13), (1, 9));
        assert_eq!(table.regressive(11, 14), (1, 10));
        assert_eq!(table.regressive(11, 15), (-1, 11));
        // Row e34
        assert_eq!(table.regressive(12, 0), (0, 0));
        assert_eq!(table.regressive(12, 1), (0, 0));
        assert_eq!(table.regressive(12, 2), (0, 0));
        assert_eq!(table.regressive(12, 3), (1, 0));
        assert_eq!(table.regressive(12, 4), (0, 0));
        assert_eq!(table.regressive(12, 5), (0, 0));
        assert_eq!(table.regressive(12, 6), (0, 0));
        assert_eq!(table.regressive(12, 7), (-1, 4));
        assert_eq!(table.regressive(12, 8), (0, 0));
        assert_eq!(table.regressive(12, 9), (0, 0));
        assert_eq!(table.regressive(12, 10), (0, 0));
        assert_eq!(table.regressive(12, 11), (-1, 8));
        assert_eq!(table.regressive(12, 12), (0, 0));
        assert_eq!(table.regressive(12, 13), (0, 0));
        assert_eq!(table.regressive(12, 14), (0, 0));
        assert_eq!(table.regressive(12, 15), (1, 12));
        // Row e134
        assert_eq!(table.regressive(13, 0), (0, 0));
        assert_eq!(table.regressive(13, 1), (0, 0));
        assert_eq!(table.regressive(13, 2), (1, 0));
        assert_eq!(table.regressive(13, 3), (-1, 1));
        assert_eq!(table.regressive(13, 4), (0, 0));
        assert_eq!(table.regressive(13, 5), (0, 0));
        assert_eq!(table.regressive(13, 6), (1, 4));
        assert_eq!(table.regressive(13, 7), (-1, 5));
        assert_eq!(table.regressive(13, 8), (0, 0));
        assert_eq!(table.regressive(13, 9), (0, 0));
        assert_eq!(table.regressive(13, 10), (1, 8));
        assert_eq!(table.regressive(13, 11), (-1, 9));
        assert_eq!(table.regressive(13, 12), (0, 0));
        assert_eq!(table.regressive(13, 13), (0, 0));
        assert_eq!(table.regressive(13, 14), (1, 12));
        assert_eq!(table.regressive(13, 15), (-1, 13));
        // Row e234
        assert_eq!(table.regressive(14, 0), (0, 0));
        assert_eq!(table.regressive(14, 1), (-1, 0));
        assert_eq!(table.regressive(14, 2), (0, 0));
        assert_eq!(table.regressive(14, 3), (-1, 2));
        assert_eq!(table.regressive(14, 4), (0, 0));
        assert_eq!(table.regressive(14, 5), (-1, 4));
        assert_eq!(table.regressive(14, 6), (0, 0));
        assert_eq!(table.regressive(14, 7), (-1, 6));
        assert_eq!(table.regressive(14, 8), (0, 0));
        assert_eq!(table.regressive(14, 9), (-1, 8));
        assert_eq!(table.regressive(14, 10), (0, 0));
        assert_eq!(table.regressive(14, 11), (-1, 10));
        assert_eq!(table.regressive(14, 12), (0, 0));
        assert_eq!(table.regressive(14, 13), (-1, 12));
        assert_eq!(table.regressive(14, 14), (0, 0));
        assert_eq!(table.regressive(14, 15), (-1, 14));
        // Row e1234
        assert_eq!(table.regressive(15, 0), (1, 0));
        assert_eq!(table.regressive(15, 1), (-1, 1));
        assert_eq!(table.regressive(15, 2), (-1, 2));
        assert_eq!(table.regressive(15, 3), (1, 3));
        assert_eq!(table.regressive(15, 4), (-1, 4));
        assert_eq!(table.regressive(15, 5), (1, 5));
        assert_eq!(table.regressive(15, 6), (1, 6));
        assert_eq!(table.regressive(15, 7), (-1, 7));
        assert_eq!(table.regressive(15, 8), (-1, 8));
        assert_eq!(table.regressive(15, 9), (1, 9));
        assert_eq!(table.regressive(15, 10), (1, 10));
        assert_eq!(table.regressive(15, 11), (-1, 11));
        assert_eq!(table.regressive(15, 12), (1, 12));
        assert_eq!(table.regressive(15, 13), (-1, 13));
        assert_eq!(table.regressive(15, 14), (-1, 14));
        assert_eq!(table.regressive(15, 15), (1, 15));
    }

    #[test]
    fn pga_3d_full_left_contraction_table() {
        let algebra = Algebra::pga(3);
        let table = ProductTable::new(&algebra);

        // Row 1
        assert_eq!(table.left_contraction(0, 0), (1, 0));
        assert_eq!(table.left_contraction(0, 1), (1, 1));
        assert_eq!(table.left_contraction(0, 2), (1, 2));
        assert_eq!(table.left_contraction(0, 3), (1, 3));
        assert_eq!(table.left_contraction(0, 4), (1, 4));
        assert_eq!(table.left_contraction(0, 5), (1, 5));
        assert_eq!(table.left_contraction(0, 6), (1, 6));
        assert_eq!(table.left_contraction(0, 7), (1, 7));
        assert_eq!(table.left_contraction(0, 8), (1, 8));
        assert_eq!(table.left_contraction(0, 9), (1, 9));
        assert_eq!(table.left_contraction(0, 10), (1, 10));
        assert_eq!(table.left_contraction(0, 11), (1, 11));
        assert_eq!(table.left_contraction(0, 12), (1, 12));
        assert_eq!(table.left_contraction(0, 13), (1, 13));
        assert_eq!(table.left_contraction(0, 14), (1, 14));
        assert_eq!(table.left_contraction(0, 15), (1, 15));
        // Row e1
        assert_eq!(table.left_contraction(1, 0), (0, 0));
        assert_eq!(table.left_contraction(1, 1), (1, 0));
        assert_eq!(table.left_contraction(1, 2), (0, 0));
        assert_eq!(table.left_contraction(1, 3), (1, 2));
        assert_eq!(table.left_contraction(1, 4), (0, 0));
        assert_eq!(table.left_contraction(1, 5), (1, 4));
        assert_eq!(table.left_contraction(1, 6), (0, 0));
        assert_eq!(table.left_contraction(1, 7), (1, 6));
        assert_eq!(table.left_contraction(1, 8), (0, 0));
        assert_eq!(table.left_contraction(1, 9), (1, 8));
        assert_eq!(table.left_contraction(1, 10), (0, 0));
        assert_eq!(table.left_contraction(1, 11), (1, 10));
        assert_eq!(table.left_contraction(1, 12), (0, 0));
        assert_eq!(table.left_contraction(1, 13), (1, 12));
        assert_eq!(table.left_contraction(1, 14), (0, 0));
        assert_eq!(table.left_contraction(1, 15), (1, 14));
        // Row e2
        assert_eq!(table.left_contraction(2, 0), (0, 0));
        assert_eq!(table.left_contraction(2, 1), (0, 0));
        assert_eq!(table.left_contraction(2, 2), (1, 0));
        assert_eq!(table.left_contraction(2, 3), (-1, 1));
        assert_eq!(table.left_contraction(2, 4), (0, 0));
        assert_eq!(table.left_contraction(2, 5), (0, 0));
        assert_eq!(table.left_contraction(2, 6), (1, 4));
        assert_eq!(table.left_contraction(2, 7), (-1, 5));
        assert_eq!(table.left_contraction(2, 8), (0, 0));
        assert_eq!(table.left_contraction(2, 9), (0, 0));
        assert_eq!(table.left_contraction(2, 10), (1, 8));
        assert_eq!(table.left_contraction(2, 11), (-1, 9));
        assert_eq!(table.left_contraction(2, 12), (0, 0));
        assert_eq!(table.left_contraction(2, 13), (0, 0));
        assert_eq!(table.left_contraction(2, 14), (1, 12));
        assert_eq!(table.left_contraction(2, 15), (-1, 13));
        // Row e12
        assert_eq!(table.left_contraction(3, 0), (0, 0));
        assert_eq!(table.left_contraction(3, 1), (0, 0));
        assert_eq!(table.left_contraction(3, 2), (0, 0));
        assert_eq!(table.left_contraction(3, 3), (-1, 0));
        assert_eq!(table.left_contraction(3, 4), (0, 0));
        assert_eq!(table.left_contraction(3, 5), (0, 0));
        assert_eq!(table.left_contraction(3, 6), (0, 0));
        assert_eq!(table.left_contraction(3, 7), (-1, 4));
        assert_eq!(table.left_contraction(3, 8), (0, 0));
        assert_eq!(table.left_contraction(3, 9), (0, 0));
        assert_eq!(table.left_contraction(3, 10), (0, 0));
        assert_eq!(table.left_contraction(3, 11), (-1, 8));
        assert_eq!(table.left_contraction(3, 12), (0, 0));
        assert_eq!(table.left_contraction(3, 13), (0, 0));
        assert_eq!(table.left_contraction(3, 14), (0, 0));
        assert_eq!(table.left_contraction(3, 15), (-1, 12));
        // Row e3
        assert_eq!(table.left_contraction(4, 0), (0, 0));
        assert_eq!(table.left_contraction(4, 1), (0, 0));
        assert_eq!(table.left_contraction(4, 2), (0, 0));
        assert_eq!(table.left_contraction(4, 3), (0, 0));
        assert_eq!(table.left_contraction(4, 4), (1, 0));
        assert_eq!(table.left_contraction(4, 5), (-1, 1));
        assert_eq!(table.left_contraction(4, 6), (-1, 2));
        assert_eq!(table.left_contraction(4, 7), (1, 3));
        assert_eq!(table.left_contraction(4, 8), (0, 0));
        assert_eq!(table.left_contraction(4, 9), (0, 0));
        assert_eq!(table.left_contraction(4, 10), (0, 0));
        assert_eq!(table.left_contraction(4, 11), (0, 0));
        assert_eq!(table.left_contraction(4, 12), (1, 8));
        assert_eq!(table.left_contraction(4, 13), (-1, 9));
        assert_eq!(table.left_contraction(4, 14), (-1, 10));
        assert_eq!(table.left_contraction(4, 15), (1, 11));
        // Row e13
        assert_eq!(table.left_contraction(5, 0), (0, 0));
        assert_eq!(table.left_contraction(5, 1), (0, 0));
        assert_eq!(table.left_contraction(5, 2), (0, 0));
        assert_eq!(table.left_contraction(5, 3), (0, 0));
        assert_eq!(table.left_contraction(5, 4), (0, 0));
        assert_eq!(table.left_contraction(5, 5), (-1, 0));
        assert_eq!(table.left_contraction(5, 6), (0, 0));
        assert_eq!(table.left_contraction(5, 7), (1, 2));
        assert_eq!(table.left_contraction(5, 8), (0, 0));
        assert_eq!(table.left_contraction(5, 9), (0, 0));
        assert_eq!(table.left_contraction(5, 10), (0, 0));
        assert_eq!(table.left_contraction(5, 11), (0, 0));
        assert_eq!(table.left_contraction(5, 12), (0, 0));
        assert_eq!(table.left_contraction(5, 13), (-1, 8));
        assert_eq!(table.left_contraction(5, 14), (0, 0));
        assert_eq!(table.left_contraction(5, 15), (1, 10));
        // Row e23
        assert_eq!(table.left_contraction(6, 0), (0, 0));
        assert_eq!(table.left_contraction(6, 1), (0, 0));
        assert_eq!(table.left_contraction(6, 2), (0, 0));
        assert_eq!(table.left_contraction(6, 3), (0, 0));
        assert_eq!(table.left_contraction(6, 4), (0, 0));
        assert_eq!(table.left_contraction(6, 5), (0, 0));
        assert_eq!(table.left_contraction(6, 6), (-1, 0));
        assert_eq!(table.left_contraction(6, 7), (-1, 1));
        assert_eq!(table.left_contraction(6, 8), (0, 0));
        assert_eq!(table.left_contraction(6, 9), (0, 0));
        assert_eq!(table.left_contraction(6, 10), (0, 0));
        assert_eq!(table.left_contraction(6, 11), (0, 0));
        assert_eq!(table.left_contraction(6, 12), (0, 0));
        assert_eq!(table.left_contraction(6, 13), (0, 0));
        assert_eq!(table.left_contraction(6, 14), (-1, 8));
        assert_eq!(table.left_contraction(6, 15), (-1, 9));
        // Row e123
        assert_eq!(table.left_contraction(7, 0), (0, 0));
        assert_eq!(table.left_contraction(7, 1), (0, 0));
        assert_eq!(table.left_contraction(7, 2), (0, 0));
        assert_eq!(table.left_contraction(7, 3), (0, 0));
        assert_eq!(table.left_contraction(7, 4), (0, 0));
        assert_eq!(table.left_contraction(7, 5), (0, 0));
        assert_eq!(table.left_contraction(7, 6), (0, 0));
        assert_eq!(table.left_contraction(7, 7), (-1, 0));
        assert_eq!(table.left_contraction(7, 8), (0, 0));
        assert_eq!(table.left_contraction(7, 9), (0, 0));
        assert_eq!(table.left_contraction(7, 10), (0, 0));
        assert_eq!(table.left_contraction(7, 11), (0, 0));
        assert_eq!(table.left_contraction(7, 12), (0, 0));
        assert_eq!(table.left_contraction(7, 13), (0, 0));
        assert_eq!(table.left_contraction(7, 14), (0, 0));
        assert_eq!(table.left_contraction(7, 15), (-1, 8));
        // Row e4
        assert_eq!(table.left_contraction(8, 0), (0, 0));
        assert_eq!(table.left_contraction(8, 1), (0, 0));
        assert_eq!(table.left_contraction(8, 2), (0, 0));
        assert_eq!(table.left_contraction(8, 3), (0, 0));
        assert_eq!(table.left_contraction(8, 4), (0, 0));
        assert_eq!(table.left_contraction(8, 5), (0, 0));
        assert_eq!(table.left_contraction(8, 6), (0, 0));
        assert_eq!(table.left_contraction(8, 7), (0, 0));
        assert_eq!(table.left_contraction(8, 8), (0, 0));
        assert_eq!(table.left_contraction(8, 9), (0, 0));
        assert_eq!(table.left_contraction(8, 10), (0, 0));
        assert_eq!(table.left_contraction(8, 11), (0, 0));
        assert_eq!(table.left_contraction(8, 12), (0, 0));
        assert_eq!(table.left_contraction(8, 13), (0, 0));
        assert_eq!(table.left_contraction(8, 14), (0, 0));
        assert_eq!(table.left_contraction(8, 15), (0, 0));
        // Row e14
        assert_eq!(table.left_contraction(9, 0), (0, 0));
        assert_eq!(table.left_contraction(9, 1), (0, 0));
        assert_eq!(table.left_contraction(9, 2), (0, 0));
        assert_eq!(table.left_contraction(9, 3), (0, 0));
        assert_eq!(table.left_contraction(9, 4), (0, 0));
        assert_eq!(table.left_contraction(9, 5), (0, 0));
        assert_eq!(table.left_contraction(9, 6), (0, 0));
        assert_eq!(table.left_contraction(9, 7), (0, 0));
        assert_eq!(table.left_contraction(9, 8), (0, 0));
        assert_eq!(table.left_contraction(9, 9), (0, 0));
        assert_eq!(table.left_contraction(9, 10), (0, 0));
        assert_eq!(table.left_contraction(9, 11), (0, 0));
        assert_eq!(table.left_contraction(9, 12), (0, 0));
        assert_eq!(table.left_contraction(9, 13), (0, 0));
        assert_eq!(table.left_contraction(9, 14), (0, 0));
        assert_eq!(table.left_contraction(9, 15), (0, 0));
        // Row e24
        assert_eq!(table.left_contraction(10, 0), (0, 0));
        assert_eq!(table.left_contraction(10, 1), (0, 0));
        assert_eq!(table.left_contraction(10, 2), (0, 0));
        assert_eq!(table.left_contraction(10, 3), (0, 0));
        assert_eq!(table.left_contraction(10, 4), (0, 0));
        assert_eq!(table.left_contraction(10, 5), (0, 0));
        assert_eq!(table.left_contraction(10, 6), (0, 0));
        assert_eq!(table.left_contraction(10, 7), (0, 0));
        assert_eq!(table.left_contraction(10, 8), (0, 0));
        assert_eq!(table.left_contraction(10, 9), (0, 0));
        assert_eq!(table.left_contraction(10, 10), (0, 0));
        assert_eq!(table.left_contraction(10, 11), (0, 0));
        assert_eq!(table.left_contraction(10, 12), (0, 0));
        assert_eq!(table.left_contraction(10, 13), (0, 0));
        assert_eq!(table.left_contraction(10, 14), (0, 0));
        assert_eq!(table.left_contraction(10, 15), (0, 0));
        // Row e124
        assert_eq!(table.left_contraction(11, 0), (0, 0));
        assert_eq!(table.left_contraction(11, 1), (0, 0));
        assert_eq!(table.left_contraction(11, 2), (0, 0));
        assert_eq!(table.left_contraction(11, 3), (0, 0));
        assert_eq!(table.left_contraction(11, 4), (0, 0));
        assert_eq!(table.left_contraction(11, 5), (0, 0));
        assert_eq!(table.left_contraction(11, 6), (0, 0));
        assert_eq!(table.left_contraction(11, 7), (0, 0));
        assert_eq!(table.left_contraction(11, 8), (0, 0));
        assert_eq!(table.left_contraction(11, 9), (0, 0));
        assert_eq!(table.left_contraction(11, 10), (0, 0));
        assert_eq!(table.left_contraction(11, 11), (0, 0));
        assert_eq!(table.left_contraction(11, 12), (0, 0));
        assert_eq!(table.left_contraction(11, 13), (0, 0));
        assert_eq!(table.left_contraction(11, 14), (0, 0));
        assert_eq!(table.left_contraction(11, 15), (0, 0));
        // Row e34
        assert_eq!(table.left_contraction(12, 0), (0, 0));
        assert_eq!(table.left_contraction(12, 1), (0, 0));
        assert_eq!(table.left_contraction(12, 2), (0, 0));
        assert_eq!(table.left_contraction(12, 3), (0, 0));
        assert_eq!(table.left_contraction(12, 4), (0, 0));
        assert_eq!(table.left_contraction(12, 5), (0, 0));
        assert_eq!(table.left_contraction(12, 6), (0, 0));
        assert_eq!(table.left_contraction(12, 7), (0, 0));
        assert_eq!(table.left_contraction(12, 8), (0, 0));
        assert_eq!(table.left_contraction(12, 9), (0, 0));
        assert_eq!(table.left_contraction(12, 10), (0, 0));
        assert_eq!(table.left_contraction(12, 11), (0, 0));
        assert_eq!(table.left_contraction(12, 12), (0, 0));
        assert_eq!(table.left_contraction(12, 13), (0, 0));
        assert_eq!(table.left_contraction(12, 14), (0, 0));
        assert_eq!(table.left_contraction(12, 15), (0, 0));
        // Row e134
        assert_eq!(table.left_contraction(13, 0), (0, 0));
        assert_eq!(table.left_contraction(13, 1), (0, 0));
        assert_eq!(table.left_contraction(13, 2), (0, 0));
        assert_eq!(table.left_contraction(13, 3), (0, 0));
        assert_eq!(table.left_contraction(13, 4), (0, 0));
        assert_eq!(table.left_contraction(13, 5), (0, 0));
        assert_eq!(table.left_contraction(13, 6), (0, 0));
        assert_eq!(table.left_contraction(13, 7), (0, 0));
        assert_eq!(table.left_contraction(13, 8), (0, 0));
        assert_eq!(table.left_contraction(13, 9), (0, 0));
        assert_eq!(table.left_contraction(13, 10), (0, 0));
        assert_eq!(table.left_contraction(13, 11), (0, 0));
        assert_eq!(table.left_contraction(13, 12), (0, 0));
        assert_eq!(table.left_contraction(13, 13), (0, 0));
        assert_eq!(table.left_contraction(13, 14), (0, 0));
        assert_eq!(table.left_contraction(13, 15), (0, 0));
        // Row e234
        assert_eq!(table.left_contraction(14, 0), (0, 0));
        assert_eq!(table.left_contraction(14, 1), (0, 0));
        assert_eq!(table.left_contraction(14, 2), (0, 0));
        assert_eq!(table.left_contraction(14, 3), (0, 0));
        assert_eq!(table.left_contraction(14, 4), (0, 0));
        assert_eq!(table.left_contraction(14, 5), (0, 0));
        assert_eq!(table.left_contraction(14, 6), (0, 0));
        assert_eq!(table.left_contraction(14, 7), (0, 0));
        assert_eq!(table.left_contraction(14, 8), (0, 0));
        assert_eq!(table.left_contraction(14, 9), (0, 0));
        assert_eq!(table.left_contraction(14, 10), (0, 0));
        assert_eq!(table.left_contraction(14, 11), (0, 0));
        assert_eq!(table.left_contraction(14, 12), (0, 0));
        assert_eq!(table.left_contraction(14, 13), (0, 0));
        assert_eq!(table.left_contraction(14, 14), (0, 0));
        assert_eq!(table.left_contraction(14, 15), (0, 0));
        // Row e1234
        assert_eq!(table.left_contraction(15, 0), (0, 0));
        assert_eq!(table.left_contraction(15, 1), (0, 0));
        assert_eq!(table.left_contraction(15, 2), (0, 0));
        assert_eq!(table.left_contraction(15, 3), (0, 0));
        assert_eq!(table.left_contraction(15, 4), (0, 0));
        assert_eq!(table.left_contraction(15, 5), (0, 0));
        assert_eq!(table.left_contraction(15, 6), (0, 0));
        assert_eq!(table.left_contraction(15, 7), (0, 0));
        assert_eq!(table.left_contraction(15, 8), (0, 0));
        assert_eq!(table.left_contraction(15, 9), (0, 0));
        assert_eq!(table.left_contraction(15, 10), (0, 0));
        assert_eq!(table.left_contraction(15, 11), (0, 0));
        assert_eq!(table.left_contraction(15, 12), (0, 0));
        assert_eq!(table.left_contraction(15, 13), (0, 0));
        assert_eq!(table.left_contraction(15, 14), (0, 0));
        assert_eq!(table.left_contraction(15, 15), (0, 0));
    }

    #[test]
    fn pga_3d_full_right_contraction_table() {
        let algebra = Algebra::pga(3);
        let table = ProductTable::new(&algebra);

        // Row 1
        assert_eq!(table.right_contraction(0, 0), (1, 0));
        assert_eq!(table.right_contraction(0, 1), (0, 0));
        assert_eq!(table.right_contraction(0, 2), (0, 0));
        assert_eq!(table.right_contraction(0, 3), (0, 0));
        assert_eq!(table.right_contraction(0, 4), (0, 0));
        assert_eq!(table.right_contraction(0, 5), (0, 0));
        assert_eq!(table.right_contraction(0, 6), (0, 0));
        assert_eq!(table.right_contraction(0, 7), (0, 0));
        assert_eq!(table.right_contraction(0, 8), (0, 0));
        assert_eq!(table.right_contraction(0, 9), (0, 0));
        assert_eq!(table.right_contraction(0, 10), (0, 0));
        assert_eq!(table.right_contraction(0, 11), (0, 0));
        assert_eq!(table.right_contraction(0, 12), (0, 0));
        assert_eq!(table.right_contraction(0, 13), (0, 0));
        assert_eq!(table.right_contraction(0, 14), (0, 0));
        assert_eq!(table.right_contraction(0, 15), (0, 0));
        // Row e1
        assert_eq!(table.right_contraction(1, 0), (1, 1));
        assert_eq!(table.right_contraction(1, 1), (1, 0));
        assert_eq!(table.right_contraction(1, 2), (0, 0));
        assert_eq!(table.right_contraction(1, 3), (0, 0));
        assert_eq!(table.right_contraction(1, 4), (0, 0));
        assert_eq!(table.right_contraction(1, 5), (0, 0));
        assert_eq!(table.right_contraction(1, 6), (0, 0));
        assert_eq!(table.right_contraction(1, 7), (0, 0));
        assert_eq!(table.right_contraction(1, 8), (0, 0));
        assert_eq!(table.right_contraction(1, 9), (0, 0));
        assert_eq!(table.right_contraction(1, 10), (0, 0));
        assert_eq!(table.right_contraction(1, 11), (0, 0));
        assert_eq!(table.right_contraction(1, 12), (0, 0));
        assert_eq!(table.right_contraction(1, 13), (0, 0));
        assert_eq!(table.right_contraction(1, 14), (0, 0));
        assert_eq!(table.right_contraction(1, 15), (0, 0));
        // Row e2
        assert_eq!(table.right_contraction(2, 0), (1, 2));
        assert_eq!(table.right_contraction(2, 1), (0, 0));
        assert_eq!(table.right_contraction(2, 2), (1, 0));
        assert_eq!(table.right_contraction(2, 3), (0, 0));
        assert_eq!(table.right_contraction(2, 4), (0, 0));
        assert_eq!(table.right_contraction(2, 5), (0, 0));
        assert_eq!(table.right_contraction(2, 6), (0, 0));
        assert_eq!(table.right_contraction(2, 7), (0, 0));
        assert_eq!(table.right_contraction(2, 8), (0, 0));
        assert_eq!(table.right_contraction(2, 9), (0, 0));
        assert_eq!(table.right_contraction(2, 10), (0, 0));
        assert_eq!(table.right_contraction(2, 11), (0, 0));
        assert_eq!(table.right_contraction(2, 12), (0, 0));
        assert_eq!(table.right_contraction(2, 13), (0, 0));
        assert_eq!(table.right_contraction(2, 14), (0, 0));
        assert_eq!(table.right_contraction(2, 15), (0, 0));
        // Row e12
        assert_eq!(table.right_contraction(3, 0), (1, 3));
        assert_eq!(table.right_contraction(3, 1), (-1, 2));
        assert_eq!(table.right_contraction(3, 2), (1, 1));
        assert_eq!(table.right_contraction(3, 3), (-1, 0));
        assert_eq!(table.right_contraction(3, 4), (0, 0));
        assert_eq!(table.right_contraction(3, 5), (0, 0));
        assert_eq!(table.right_contraction(3, 6), (0, 0));
        assert_eq!(table.right_contraction(3, 7), (0, 0));
        assert_eq!(table.right_contraction(3, 8), (0, 0));
        assert_eq!(table.right_contraction(3, 9), (0, 0));
        assert_eq!(table.right_contraction(3, 10), (0, 0));
        assert_eq!(table.right_contraction(3, 11), (0, 0));
        assert_eq!(table.right_contraction(3, 12), (0, 0));
        assert_eq!(table.right_contraction(3, 13), (0, 0));
        assert_eq!(table.right_contraction(3, 14), (0, 0));
        assert_eq!(table.right_contraction(3, 15), (0, 0));
        // Row e3
        assert_eq!(table.right_contraction(4, 0), (1, 4));
        assert_eq!(table.right_contraction(4, 1), (0, 0));
        assert_eq!(table.right_contraction(4, 2), (0, 0));
        assert_eq!(table.right_contraction(4, 3), (0, 0));
        assert_eq!(table.right_contraction(4, 4), (1, 0));
        assert_eq!(table.right_contraction(4, 5), (0, 0));
        assert_eq!(table.right_contraction(4, 6), (0, 0));
        assert_eq!(table.right_contraction(4, 7), (0, 0));
        assert_eq!(table.right_contraction(4, 8), (0, 0));
        assert_eq!(table.right_contraction(4, 9), (0, 0));
        assert_eq!(table.right_contraction(4, 10), (0, 0));
        assert_eq!(table.right_contraction(4, 11), (0, 0));
        assert_eq!(table.right_contraction(4, 12), (0, 0));
        assert_eq!(table.right_contraction(4, 13), (0, 0));
        assert_eq!(table.right_contraction(4, 14), (0, 0));
        assert_eq!(table.right_contraction(4, 15), (0, 0));
        // Row e13
        assert_eq!(table.right_contraction(5, 0), (1, 5));
        assert_eq!(table.right_contraction(5, 1), (-1, 4));
        assert_eq!(table.right_contraction(5, 2), (0, 0));
        assert_eq!(table.right_contraction(5, 3), (0, 0));
        assert_eq!(table.right_contraction(5, 4), (1, 1));
        assert_eq!(table.right_contraction(5, 5), (-1, 0));
        assert_eq!(table.right_contraction(5, 6), (0, 0));
        assert_eq!(table.right_contraction(5, 7), (0, 0));
        assert_eq!(table.right_contraction(5, 8), (0, 0));
        assert_eq!(table.right_contraction(5, 9), (0, 0));
        assert_eq!(table.right_contraction(5, 10), (0, 0));
        assert_eq!(table.right_contraction(5, 11), (0, 0));
        assert_eq!(table.right_contraction(5, 12), (0, 0));
        assert_eq!(table.right_contraction(5, 13), (0, 0));
        assert_eq!(table.right_contraction(5, 14), (0, 0));
        assert_eq!(table.right_contraction(5, 15), (0, 0));
        // Row e23
        assert_eq!(table.right_contraction(6, 0), (1, 6));
        assert_eq!(table.right_contraction(6, 1), (0, 0));
        assert_eq!(table.right_contraction(6, 2), (-1, 4));
        assert_eq!(table.right_contraction(6, 3), (0, 0));
        assert_eq!(table.right_contraction(6, 4), (1, 2));
        assert_eq!(table.right_contraction(6, 5), (0, 0));
        assert_eq!(table.right_contraction(6, 6), (-1, 0));
        assert_eq!(table.right_contraction(6, 7), (0, 0));
        assert_eq!(table.right_contraction(6, 8), (0, 0));
        assert_eq!(table.right_contraction(6, 9), (0, 0));
        assert_eq!(table.right_contraction(6, 10), (0, 0));
        assert_eq!(table.right_contraction(6, 11), (0, 0));
        assert_eq!(table.right_contraction(6, 12), (0, 0));
        assert_eq!(table.right_contraction(6, 13), (0, 0));
        assert_eq!(table.right_contraction(6, 14), (0, 0));
        assert_eq!(table.right_contraction(6, 15), (0, 0));
        // Row e123
        assert_eq!(table.right_contraction(7, 0), (1, 7));
        assert_eq!(table.right_contraction(7, 1), (1, 6));
        assert_eq!(table.right_contraction(7, 2), (-1, 5));
        assert_eq!(table.right_contraction(7, 3), (-1, 4));
        assert_eq!(table.right_contraction(7, 4), (1, 3));
        assert_eq!(table.right_contraction(7, 5), (1, 2));
        assert_eq!(table.right_contraction(7, 6), (-1, 1));
        assert_eq!(table.right_contraction(7, 7), (-1, 0));
        assert_eq!(table.right_contraction(7, 8), (0, 0));
        assert_eq!(table.right_contraction(7, 9), (0, 0));
        assert_eq!(table.right_contraction(7, 10), (0, 0));
        assert_eq!(table.right_contraction(7, 11), (0, 0));
        assert_eq!(table.right_contraction(7, 12), (0, 0));
        assert_eq!(table.right_contraction(7, 13), (0, 0));
        assert_eq!(table.right_contraction(7, 14), (0, 0));
        assert_eq!(table.right_contraction(7, 15), (0, 0));
        // Row e4
        assert_eq!(table.right_contraction(8, 0), (1, 8));
        assert_eq!(table.right_contraction(8, 1), (0, 0));
        assert_eq!(table.right_contraction(8, 2), (0, 0));
        assert_eq!(table.right_contraction(8, 3), (0, 0));
        assert_eq!(table.right_contraction(8, 4), (0, 0));
        assert_eq!(table.right_contraction(8, 5), (0, 0));
        assert_eq!(table.right_contraction(8, 6), (0, 0));
        assert_eq!(table.right_contraction(8, 7), (0, 0));
        assert_eq!(table.right_contraction(8, 8), (0, 0));
        assert_eq!(table.right_contraction(8, 9), (0, 0));
        assert_eq!(table.right_contraction(8, 10), (0, 0));
        assert_eq!(table.right_contraction(8, 11), (0, 0));
        assert_eq!(table.right_contraction(8, 12), (0, 0));
        assert_eq!(table.right_contraction(8, 13), (0, 0));
        assert_eq!(table.right_contraction(8, 14), (0, 0));
        assert_eq!(table.right_contraction(8, 15), (0, 0));
        // Row e14
        assert_eq!(table.right_contraction(9, 0), (1, 9));
        assert_eq!(table.right_contraction(9, 1), (-1, 8));
        assert_eq!(table.right_contraction(9, 2), (0, 0));
        assert_eq!(table.right_contraction(9, 3), (0, 0));
        assert_eq!(table.right_contraction(9, 4), (0, 0));
        assert_eq!(table.right_contraction(9, 5), (0, 0));
        assert_eq!(table.right_contraction(9, 6), (0, 0));
        assert_eq!(table.right_contraction(9, 7), (0, 0));
        assert_eq!(table.right_contraction(9, 8), (0, 0));
        assert_eq!(table.right_contraction(9, 9), (0, 0));
        assert_eq!(table.right_contraction(9, 10), (0, 0));
        assert_eq!(table.right_contraction(9, 11), (0, 0));
        assert_eq!(table.right_contraction(9, 12), (0, 0));
        assert_eq!(table.right_contraction(9, 13), (0, 0));
        assert_eq!(table.right_contraction(9, 14), (0, 0));
        assert_eq!(table.right_contraction(9, 15), (0, 0));
        // Row e24
        assert_eq!(table.right_contraction(10, 0), (1, 10));
        assert_eq!(table.right_contraction(10, 1), (0, 0));
        assert_eq!(table.right_contraction(10, 2), (-1, 8));
        assert_eq!(table.right_contraction(10, 3), (0, 0));
        assert_eq!(table.right_contraction(10, 4), (0, 0));
        assert_eq!(table.right_contraction(10, 5), (0, 0));
        assert_eq!(table.right_contraction(10, 6), (0, 0));
        assert_eq!(table.right_contraction(10, 7), (0, 0));
        assert_eq!(table.right_contraction(10, 8), (0, 0));
        assert_eq!(table.right_contraction(10, 9), (0, 0));
        assert_eq!(table.right_contraction(10, 10), (0, 0));
        assert_eq!(table.right_contraction(10, 11), (0, 0));
        assert_eq!(table.right_contraction(10, 12), (0, 0));
        assert_eq!(table.right_contraction(10, 13), (0, 0));
        assert_eq!(table.right_contraction(10, 14), (0, 0));
        assert_eq!(table.right_contraction(10, 15), (0, 0));
        // Row e124
        assert_eq!(table.right_contraction(11, 0), (1, 11));
        assert_eq!(table.right_contraction(11, 1), (1, 10));
        assert_eq!(table.right_contraction(11, 2), (-1, 9));
        assert_eq!(table.right_contraction(11, 3), (-1, 8));
        assert_eq!(table.right_contraction(11, 4), (0, 0));
        assert_eq!(table.right_contraction(11, 5), (0, 0));
        assert_eq!(table.right_contraction(11, 6), (0, 0));
        assert_eq!(table.right_contraction(11, 7), (0, 0));
        assert_eq!(table.right_contraction(11, 8), (0, 0));
        assert_eq!(table.right_contraction(11, 9), (0, 0));
        assert_eq!(table.right_contraction(11, 10), (0, 0));
        assert_eq!(table.right_contraction(11, 11), (0, 0));
        assert_eq!(table.right_contraction(11, 12), (0, 0));
        assert_eq!(table.right_contraction(11, 13), (0, 0));
        assert_eq!(table.right_contraction(11, 14), (0, 0));
        assert_eq!(table.right_contraction(11, 15), (0, 0));
        // Row e34
        assert_eq!(table.right_contraction(12, 0), (1, 12));
        assert_eq!(table.right_contraction(12, 1), (0, 0));
        assert_eq!(table.right_contraction(12, 2), (0, 0));
        assert_eq!(table.right_contraction(12, 3), (0, 0));
        assert_eq!(table.right_contraction(12, 4), (-1, 8));
        assert_eq!(table.right_contraction(12, 5), (0, 0));
        assert_eq!(table.right_contraction(12, 6), (0, 0));
        assert_eq!(table.right_contraction(12, 7), (0, 0));
        assert_eq!(table.right_contraction(12, 8), (0, 0));
        assert_eq!(table.right_contraction(12, 9), (0, 0));
        assert_eq!(table.right_contraction(12, 10), (0, 0));
        assert_eq!(table.right_contraction(12, 11), (0, 0));
        assert_eq!(table.right_contraction(12, 12), (0, 0));
        assert_eq!(table.right_contraction(12, 13), (0, 0));
        assert_eq!(table.right_contraction(12, 14), (0, 0));
        assert_eq!(table.right_contraction(12, 15), (0, 0));
        // Row e134
        assert_eq!(table.right_contraction(13, 0), (1, 13));
        assert_eq!(table.right_contraction(13, 1), (1, 12));
        assert_eq!(table.right_contraction(13, 2), (0, 0));
        assert_eq!(table.right_contraction(13, 3), (0, 0));
        assert_eq!(table.right_contraction(13, 4), (-1, 9));
        assert_eq!(table.right_contraction(13, 5), (-1, 8));
        assert_eq!(table.right_contraction(13, 6), (0, 0));
        assert_eq!(table.right_contraction(13, 7), (0, 0));
        assert_eq!(table.right_contraction(13, 8), (0, 0));
        assert_eq!(table.right_contraction(13, 9), (0, 0));
        assert_eq!(table.right_contraction(13, 10), (0, 0));
        assert_eq!(table.right_contraction(13, 11), (0, 0));
        assert_eq!(table.right_contraction(13, 12), (0, 0));
        assert_eq!(table.right_contraction(13, 13), (0, 0));
        assert_eq!(table.right_contraction(13, 14), (0, 0));
        assert_eq!(table.right_contraction(13, 15), (0, 0));
        // Row e234
        assert_eq!(table.right_contraction(14, 0), (1, 14));
        assert_eq!(table.right_contraction(14, 1), (0, 0));
        assert_eq!(table.right_contraction(14, 2), (1, 12));
        assert_eq!(table.right_contraction(14, 3), (0, 0));
        assert_eq!(table.right_contraction(14, 4), (-1, 10));
        assert_eq!(table.right_contraction(14, 5), (0, 0));
        assert_eq!(table.right_contraction(14, 6), (-1, 8));
        assert_eq!(table.right_contraction(14, 7), (0, 0));
        assert_eq!(table.right_contraction(14, 8), (0, 0));
        assert_eq!(table.right_contraction(14, 9), (0, 0));
        assert_eq!(table.right_contraction(14, 10), (0, 0));
        assert_eq!(table.right_contraction(14, 11), (0, 0));
        assert_eq!(table.right_contraction(14, 12), (0, 0));
        assert_eq!(table.right_contraction(14, 13), (0, 0));
        assert_eq!(table.right_contraction(14, 14), (0, 0));
        assert_eq!(table.right_contraction(14, 15), (0, 0));
        // Row e1234
        assert_eq!(table.right_contraction(15, 0), (1, 15));
        assert_eq!(table.right_contraction(15, 1), (-1, 14));
        assert_eq!(table.right_contraction(15, 2), (1, 13));
        assert_eq!(table.right_contraction(15, 3), (-1, 12));
        assert_eq!(table.right_contraction(15, 4), (-1, 11));
        assert_eq!(table.right_contraction(15, 5), (1, 10));
        assert_eq!(table.right_contraction(15, 6), (-1, 9));
        assert_eq!(table.right_contraction(15, 7), (1, 8));
        assert_eq!(table.right_contraction(15, 8), (0, 0));
        assert_eq!(table.right_contraction(15, 9), (0, 0));
        assert_eq!(table.right_contraction(15, 10), (0, 0));
        assert_eq!(table.right_contraction(15, 11), (0, 0));
        assert_eq!(table.right_contraction(15, 12), (0, 0));
        assert_eq!(table.right_contraction(15, 13), (0, 0));
        assert_eq!(table.right_contraction(15, 14), (0, 0));
        assert_eq!(table.right_contraction(15, 15), (0, 0));
    }

    #[test]
    fn pga_3d_full_dot_table() {
        let algebra = Algebra::pga(3);
        let table = ProductTable::new(&algebra);

        // Row 1
        assert_eq!(table.dot(0, 0), (1, 0));
        assert_eq!(table.dot(0, 1), (0, 0));
        assert_eq!(table.dot(0, 2), (0, 0));
        assert_eq!(table.dot(0, 3), (0, 0));
        assert_eq!(table.dot(0, 4), (0, 0));
        assert_eq!(table.dot(0, 5), (0, 0));
        assert_eq!(table.dot(0, 6), (0, 0));
        assert_eq!(table.dot(0, 7), (0, 0));
        assert_eq!(table.dot(0, 8), (0, 0));
        assert_eq!(table.dot(0, 9), (0, 0));
        assert_eq!(table.dot(0, 10), (0, 0));
        assert_eq!(table.dot(0, 11), (0, 0));
        assert_eq!(table.dot(0, 12), (0, 0));
        assert_eq!(table.dot(0, 13), (0, 0));
        assert_eq!(table.dot(0, 14), (0, 0));
        assert_eq!(table.dot(0, 15), (0, 0));
        // Row e1
        assert_eq!(table.dot(1, 0), (0, 0));
        assert_eq!(table.dot(1, 1), (1, 0));
        assert_eq!(table.dot(1, 2), (0, 0));
        assert_eq!(table.dot(1, 3), (0, 0));
        assert_eq!(table.dot(1, 4), (0, 0));
        assert_eq!(table.dot(1, 5), (0, 0));
        assert_eq!(table.dot(1, 6), (0, 0));
        assert_eq!(table.dot(1, 7), (0, 0));
        assert_eq!(table.dot(1, 8), (0, 0));
        assert_eq!(table.dot(1, 9), (0, 0));
        assert_eq!(table.dot(1, 10), (0, 0));
        assert_eq!(table.dot(1, 11), (0, 0));
        assert_eq!(table.dot(1, 12), (0, 0));
        assert_eq!(table.dot(1, 13), (0, 0));
        assert_eq!(table.dot(1, 14), (0, 0));
        assert_eq!(table.dot(1, 15), (0, 0));
        // Row e2
        assert_eq!(table.dot(2, 0), (0, 0));
        assert_eq!(table.dot(2, 1), (0, 0));
        assert_eq!(table.dot(2, 2), (1, 0));
        assert_eq!(table.dot(2, 3), (0, 0));
        assert_eq!(table.dot(2, 4), (0, 0));
        assert_eq!(table.dot(2, 5), (0, 0));
        assert_eq!(table.dot(2, 6), (0, 0));
        assert_eq!(table.dot(2, 7), (0, 0));
        assert_eq!(table.dot(2, 8), (0, 0));
        assert_eq!(table.dot(2, 9), (0, 0));
        assert_eq!(table.dot(2, 10), (0, 0));
        assert_eq!(table.dot(2, 11), (0, 0));
        assert_eq!(table.dot(2, 12), (0, 0));
        assert_eq!(table.dot(2, 13), (0, 0));
        assert_eq!(table.dot(2, 14), (0, 0));
        assert_eq!(table.dot(2, 15), (0, 0));
        // Row e12
        assert_eq!(table.dot(3, 0), (0, 0));
        assert_eq!(table.dot(3, 1), (0, 0));
        assert_eq!(table.dot(3, 2), (0, 0));
        assert_eq!(table.dot(3, 3), (-1, 0));
        assert_eq!(table.dot(3, 4), (0, 0));
        assert_eq!(table.dot(3, 5), (0, 0));
        assert_eq!(table.dot(3, 6), (0, 0));
        assert_eq!(table.dot(3, 7), (0, 0));
        assert_eq!(table.dot(3, 8), (0, 0));
        assert_eq!(table.dot(3, 9), (0, 0));
        assert_eq!(table.dot(3, 10), (0, 0));
        assert_eq!(table.dot(3, 11), (0, 0));
        assert_eq!(table.dot(3, 12), (0, 0));
        assert_eq!(table.dot(3, 13), (0, 0));
        assert_eq!(table.dot(3, 14), (0, 0));
        assert_eq!(table.dot(3, 15), (0, 0));
        // Row e3
        assert_eq!(table.dot(4, 0), (0, 0));
        assert_eq!(table.dot(4, 1), (0, 0));
        assert_eq!(table.dot(4, 2), (0, 0));
        assert_eq!(table.dot(4, 3), (0, 0));
        assert_eq!(table.dot(4, 4), (1, 0));
        assert_eq!(table.dot(4, 5), (0, 0));
        assert_eq!(table.dot(4, 6), (0, 0));
        assert_eq!(table.dot(4, 7), (0, 0));
        assert_eq!(table.dot(4, 8), (0, 0));
        assert_eq!(table.dot(4, 9), (0, 0));
        assert_eq!(table.dot(4, 10), (0, 0));
        assert_eq!(table.dot(4, 11), (0, 0));
        assert_eq!(table.dot(4, 12), (0, 0));
        assert_eq!(table.dot(4, 13), (0, 0));
        assert_eq!(table.dot(4, 14), (0, 0));
        assert_eq!(table.dot(4, 15), (0, 0));
        // Row e13
        assert_eq!(table.dot(5, 0), (0, 0));
        assert_eq!(table.dot(5, 1), (0, 0));
        assert_eq!(table.dot(5, 2), (0, 0));
        assert_eq!(table.dot(5, 3), (0, 0));
        assert_eq!(table.dot(5, 4), (0, 0));
        assert_eq!(table.dot(5, 5), (-1, 0));
        assert_eq!(table.dot(5, 6), (0, 0));
        assert_eq!(table.dot(5, 7), (0, 0));
        assert_eq!(table.dot(5, 8), (0, 0));
        assert_eq!(table.dot(5, 9), (0, 0));
        assert_eq!(table.dot(5, 10), (0, 0));
        assert_eq!(table.dot(5, 11), (0, 0));
        assert_eq!(table.dot(5, 12), (0, 0));
        assert_eq!(table.dot(5, 13), (0, 0));
        assert_eq!(table.dot(5, 14), (0, 0));
        assert_eq!(table.dot(5, 15), (0, 0));
        // Row e23
        assert_eq!(table.dot(6, 0), (0, 0));
        assert_eq!(table.dot(6, 1), (0, 0));
        assert_eq!(table.dot(6, 2), (0, 0));
        assert_eq!(table.dot(6, 3), (0, 0));
        assert_eq!(table.dot(6, 4), (0, 0));
        assert_eq!(table.dot(6, 5), (0, 0));
        assert_eq!(table.dot(6, 6), (-1, 0));
        assert_eq!(table.dot(6, 7), (0, 0));
        assert_eq!(table.dot(6, 8), (0, 0));
        assert_eq!(table.dot(6, 9), (0, 0));
        assert_eq!(table.dot(6, 10), (0, 0));
        assert_eq!(table.dot(6, 11), (0, 0));
        assert_eq!(table.dot(6, 12), (0, 0));
        assert_eq!(table.dot(6, 13), (0, 0));
        assert_eq!(table.dot(6, 14), (0, 0));
        assert_eq!(table.dot(6, 15), (0, 0));
        // Row e123
        assert_eq!(table.dot(7, 0), (0, 0));
        assert_eq!(table.dot(7, 1), (0, 0));
        assert_eq!(table.dot(7, 2), (0, 0));
        assert_eq!(table.dot(7, 3), (0, 0));
        assert_eq!(table.dot(7, 4), (0, 0));
        assert_eq!(table.dot(7, 5), (0, 0));
        assert_eq!(table.dot(7, 6), (0, 0));
        assert_eq!(table.dot(7, 7), (-1, 0));
        assert_eq!(table.dot(7, 8), (0, 0));
        assert_eq!(table.dot(7, 9), (0, 0));
        assert_eq!(table.dot(7, 10), (0, 0));
        assert_eq!(table.dot(7, 11), (0, 0));
        assert_eq!(table.dot(7, 12), (0, 0));
        assert_eq!(table.dot(7, 13), (0, 0));
        assert_eq!(table.dot(7, 14), (0, 0));
        assert_eq!(table.dot(7, 15), (0, 0));
        // Row e4
        assert_eq!(table.dot(8, 0), (0, 0));
        assert_eq!(table.dot(8, 1), (0, 0));
        assert_eq!(table.dot(8, 2), (0, 0));
        assert_eq!(table.dot(8, 3), (0, 0));
        assert_eq!(table.dot(8, 4), (0, 0));
        assert_eq!(table.dot(8, 5), (0, 0));
        assert_eq!(table.dot(8, 6), (0, 0));
        assert_eq!(table.dot(8, 7), (0, 0));
        assert_eq!(table.dot(8, 8), (0, 0));
        assert_eq!(table.dot(8, 9), (0, 0));
        assert_eq!(table.dot(8, 10), (0, 0));
        assert_eq!(table.dot(8, 11), (0, 0));
        assert_eq!(table.dot(8, 12), (0, 0));
        assert_eq!(table.dot(8, 13), (0, 0));
        assert_eq!(table.dot(8, 14), (0, 0));
        assert_eq!(table.dot(8, 15), (0, 0));
        // Row e14
        assert_eq!(table.dot(9, 0), (0, 0));
        assert_eq!(table.dot(9, 1), (0, 0));
        assert_eq!(table.dot(9, 2), (0, 0));
        assert_eq!(table.dot(9, 3), (0, 0));
        assert_eq!(table.dot(9, 4), (0, 0));
        assert_eq!(table.dot(9, 5), (0, 0));
        assert_eq!(table.dot(9, 6), (0, 0));
        assert_eq!(table.dot(9, 7), (0, 0));
        assert_eq!(table.dot(9, 8), (0, 0));
        assert_eq!(table.dot(9, 9), (0, 0));
        assert_eq!(table.dot(9, 10), (0, 0));
        assert_eq!(table.dot(9, 11), (0, 0));
        assert_eq!(table.dot(9, 12), (0, 0));
        assert_eq!(table.dot(9, 13), (0, 0));
        assert_eq!(table.dot(9, 14), (0, 0));
        assert_eq!(table.dot(9, 15), (0, 0));
        // Row e24
        assert_eq!(table.dot(10, 0), (0, 0));
        assert_eq!(table.dot(10, 1), (0, 0));
        assert_eq!(table.dot(10, 2), (0, 0));
        assert_eq!(table.dot(10, 3), (0, 0));
        assert_eq!(table.dot(10, 4), (0, 0));
        assert_eq!(table.dot(10, 5), (0, 0));
        assert_eq!(table.dot(10, 6), (0, 0));
        assert_eq!(table.dot(10, 7), (0, 0));
        assert_eq!(table.dot(10, 8), (0, 0));
        assert_eq!(table.dot(10, 9), (0, 0));
        assert_eq!(table.dot(10, 10), (0, 0));
        assert_eq!(table.dot(10, 11), (0, 0));
        assert_eq!(table.dot(10, 12), (0, 0));
        assert_eq!(table.dot(10, 13), (0, 0));
        assert_eq!(table.dot(10, 14), (0, 0));
        assert_eq!(table.dot(10, 15), (0, 0));
        // Row e124
        assert_eq!(table.dot(11, 0), (0, 0));
        assert_eq!(table.dot(11, 1), (0, 0));
        assert_eq!(table.dot(11, 2), (0, 0));
        assert_eq!(table.dot(11, 3), (0, 0));
        assert_eq!(table.dot(11, 4), (0, 0));
        assert_eq!(table.dot(11, 5), (0, 0));
        assert_eq!(table.dot(11, 6), (0, 0));
        assert_eq!(table.dot(11, 7), (0, 0));
        assert_eq!(table.dot(11, 8), (0, 0));
        assert_eq!(table.dot(11, 9), (0, 0));
        assert_eq!(table.dot(11, 10), (0, 0));
        assert_eq!(table.dot(11, 11), (0, 0));
        assert_eq!(table.dot(11, 12), (0, 0));
        assert_eq!(table.dot(11, 13), (0, 0));
        assert_eq!(table.dot(11, 14), (0, 0));
        assert_eq!(table.dot(11, 15), (0, 0));
        // Row e34
        assert_eq!(table.dot(12, 0), (0, 0));
        assert_eq!(table.dot(12, 1), (0, 0));
        assert_eq!(table.dot(12, 2), (0, 0));
        assert_eq!(table.dot(12, 3), (0, 0));
        assert_eq!(table.dot(12, 4), (0, 0));
        assert_eq!(table.dot(12, 5), (0, 0));
        assert_eq!(table.dot(12, 6), (0, 0));
        assert_eq!(table.dot(12, 7), (0, 0));
        assert_eq!(table.dot(12, 8), (0, 0));
        assert_eq!(table.dot(12, 9), (0, 0));
        assert_eq!(table.dot(12, 10), (0, 0));
        assert_eq!(table.dot(12, 11), (0, 0));
        assert_eq!(table.dot(12, 12), (0, 0));
        assert_eq!(table.dot(12, 13), (0, 0));
        assert_eq!(table.dot(12, 14), (0, 0));
        assert_eq!(table.dot(12, 15), (0, 0));
        // Row e134
        assert_eq!(table.dot(13, 0), (0, 0));
        assert_eq!(table.dot(13, 1), (0, 0));
        assert_eq!(table.dot(13, 2), (0, 0));
        assert_eq!(table.dot(13, 3), (0, 0));
        assert_eq!(table.dot(13, 4), (0, 0));
        assert_eq!(table.dot(13, 5), (0, 0));
        assert_eq!(table.dot(13, 6), (0, 0));
        assert_eq!(table.dot(13, 7), (0, 0));
        assert_eq!(table.dot(13, 8), (0, 0));
        assert_eq!(table.dot(13, 9), (0, 0));
        assert_eq!(table.dot(13, 10), (0, 0));
        assert_eq!(table.dot(13, 11), (0, 0));
        assert_eq!(table.dot(13, 12), (0, 0));
        assert_eq!(table.dot(13, 13), (0, 0));
        assert_eq!(table.dot(13, 14), (0, 0));
        assert_eq!(table.dot(13, 15), (0, 0));
        // Row e234
        assert_eq!(table.dot(14, 0), (0, 0));
        assert_eq!(table.dot(14, 1), (0, 0));
        assert_eq!(table.dot(14, 2), (0, 0));
        assert_eq!(table.dot(14, 3), (0, 0));
        assert_eq!(table.dot(14, 4), (0, 0));
        assert_eq!(table.dot(14, 5), (0, 0));
        assert_eq!(table.dot(14, 6), (0, 0));
        assert_eq!(table.dot(14, 7), (0, 0));
        assert_eq!(table.dot(14, 8), (0, 0));
        assert_eq!(table.dot(14, 9), (0, 0));
        assert_eq!(table.dot(14, 10), (0, 0));
        assert_eq!(table.dot(14, 11), (0, 0));
        assert_eq!(table.dot(14, 12), (0, 0));
        assert_eq!(table.dot(14, 13), (0, 0));
        assert_eq!(table.dot(14, 14), (0, 0));
        assert_eq!(table.dot(14, 15), (0, 0));
        // Row e1234
        assert_eq!(table.dot(15, 0), (0, 0));
        assert_eq!(table.dot(15, 1), (0, 0));
        assert_eq!(table.dot(15, 2), (0, 0));
        assert_eq!(table.dot(15, 3), (0, 0));
        assert_eq!(table.dot(15, 4), (0, 0));
        assert_eq!(table.dot(15, 5), (0, 0));
        assert_eq!(table.dot(15, 6), (0, 0));
        assert_eq!(table.dot(15, 7), (0, 0));
        assert_eq!(table.dot(15, 8), (0, 0));
        assert_eq!(table.dot(15, 9), (0, 0));
        assert_eq!(table.dot(15, 10), (0, 0));
        assert_eq!(table.dot(15, 11), (0, 0));
        assert_eq!(table.dot(15, 12), (0, 0));
        assert_eq!(table.dot(15, 13), (0, 0));
        assert_eq!(table.dot(15, 14), (0, 0));
        assert_eq!(table.dot(15, 15), (0, 0));
    }

    #[test]
    fn pga_3d_full_antiproduct_table() {
        let algebra = Algebra::pga(3);
        let table = ProductTable::new(&algebra);

        // Row scalar
        assert_eq!(table.antiproduct(0, 0), (0, 0));
        assert_eq!(table.antiproduct(0, 1), (0, 0));
        assert_eq!(table.antiproduct(0, 2), (0, 0));
        assert_eq!(table.antiproduct(0, 3), (0, 0));
        assert_eq!(table.antiproduct(0, 4), (0, 0));
        assert_eq!(table.antiproduct(0, 5), (0, 0));
        assert_eq!(table.antiproduct(0, 6), (0, 0));
        assert_eq!(table.antiproduct(0, 7), (0, 0));
        assert_eq!(table.antiproduct(0, 8), (1, 7));
        assert_eq!(table.antiproduct(0, 9), (-1, 6));
        assert_eq!(table.antiproduct(0, 10), (1, 5));
        assert_eq!(table.antiproduct(0, 11), (-1, 4));
        assert_eq!(table.antiproduct(0, 12), (-1, 3));
        assert_eq!(table.antiproduct(0, 13), (1, 2));
        assert_eq!(table.antiproduct(0, 14), (-1, 1));
        assert_eq!(table.antiproduct(0, 15), (1, 0));
        // Row e1
        assert_eq!(table.antiproduct(1, 0), (0, 0));
        assert_eq!(table.antiproduct(1, 1), (0, 0));
        assert_eq!(table.antiproduct(1, 2), (0, 0));
        assert_eq!(table.antiproduct(1, 3), (0, 0));
        assert_eq!(table.antiproduct(1, 4), (0, 0));
        assert_eq!(table.antiproduct(1, 5), (0, 0));
        assert_eq!(table.antiproduct(1, 6), (0, 0));
        assert_eq!(table.antiproduct(1, 7), (0, 0));
        assert_eq!(table.antiproduct(1, 8), (-1, 6));
        assert_eq!(table.antiproduct(1, 9), (1, 7));
        assert_eq!(table.antiproduct(1, 10), (1, 4));
        assert_eq!(table.antiproduct(1, 11), (-1, 5));
        assert_eq!(table.antiproduct(1, 12), (-1, 2));
        assert_eq!(table.antiproduct(1, 13), (1, 3));
        assert_eq!(table.antiproduct(1, 14), (1, 0));
        assert_eq!(table.antiproduct(1, 15), (-1, 1));
        // Row e2
        assert_eq!(table.antiproduct(2, 0), (0, 0));
        assert_eq!(table.antiproduct(2, 1), (0, 0));
        assert_eq!(table.antiproduct(2, 2), (0, 0));
        assert_eq!(table.antiproduct(2, 3), (0, 0));
        assert_eq!(table.antiproduct(2, 4), (0, 0));
        assert_eq!(table.antiproduct(2, 5), (0, 0));
        assert_eq!(table.antiproduct(2, 6), (0, 0));
        assert_eq!(table.antiproduct(2, 7), (0, 0));
        assert_eq!(table.antiproduct(2, 8), (1, 5));
        assert_eq!(table.antiproduct(2, 9), (-1, 4));
        assert_eq!(table.antiproduct(2, 10), (1, 7));
        assert_eq!(table.antiproduct(2, 11), (-1, 6));
        assert_eq!(table.antiproduct(2, 12), (1, 1));
        assert_eq!(table.antiproduct(2, 13), (-1, 0));
        assert_eq!(table.antiproduct(2, 14), (1, 3));
        assert_eq!(table.antiproduct(2, 15), (-1, 2));
        // Row e12
        assert_eq!(table.antiproduct(3, 0), (0, 0));
        assert_eq!(table.antiproduct(3, 1), (0, 0));
        assert_eq!(table.antiproduct(3, 2), (0, 0));
        assert_eq!(table.antiproduct(3, 3), (0, 0));
        assert_eq!(table.antiproduct(3, 4), (0, 0));
        assert_eq!(table.antiproduct(3, 5), (0, 0));
        assert_eq!(table.antiproduct(3, 6), (0, 0));
        assert_eq!(table.antiproduct(3, 7), (0, 0));
        assert_eq!(table.antiproduct(3, 8), (-1, 4));
        assert_eq!(table.antiproduct(3, 9), (1, 5));
        assert_eq!(table.antiproduct(3, 10), (1, 6));
        assert_eq!(table.antiproduct(3, 11), (-1, 7));
        assert_eq!(table.antiproduct(3, 12), (1, 0));
        assert_eq!(table.antiproduct(3, 13), (-1, 1));
        assert_eq!(table.antiproduct(3, 14), (-1, 2));
        assert_eq!(table.antiproduct(3, 15), (1, 3));
        // Row e3
        assert_eq!(table.antiproduct(4, 0), (0, 0));
        assert_eq!(table.antiproduct(4, 1), (0, 0));
        assert_eq!(table.antiproduct(4, 2), (0, 0));
        assert_eq!(table.antiproduct(4, 3), (0, 0));
        assert_eq!(table.antiproduct(4, 4), (0, 0));
        assert_eq!(table.antiproduct(4, 5), (0, 0));
        assert_eq!(table.antiproduct(4, 6), (0, 0));
        assert_eq!(table.antiproduct(4, 7), (0, 0));
        assert_eq!(table.antiproduct(4, 8), (-1, 3));
        assert_eq!(table.antiproduct(4, 9), (1, 2));
        assert_eq!(table.antiproduct(4, 10), (-1, 1));
        assert_eq!(table.antiproduct(4, 11), (1, 0));
        assert_eq!(table.antiproduct(4, 12), (1, 7));
        assert_eq!(table.antiproduct(4, 13), (-1, 6));
        assert_eq!(table.antiproduct(4, 14), (1, 5));
        assert_eq!(table.antiproduct(4, 15), (-1, 4));
        // Row e13
        assert_eq!(table.antiproduct(5, 0), (0, 0));
        assert_eq!(table.antiproduct(5, 1), (0, 0));
        assert_eq!(table.antiproduct(5, 2), (0, 0));
        assert_eq!(table.antiproduct(5, 3), (0, 0));
        assert_eq!(table.antiproduct(5, 4), (0, 0));
        assert_eq!(table.antiproduct(5, 5), (0, 0));
        assert_eq!(table.antiproduct(5, 6), (0, 0));
        assert_eq!(table.antiproduct(5, 7), (0, 0));
        assert_eq!(table.antiproduct(5, 8), (1, 2));
        assert_eq!(table.antiproduct(5, 9), (-1, 3));
        assert_eq!(table.antiproduct(5, 10), (-1, 0));
        assert_eq!(table.antiproduct(5, 11), (1, 1));
        assert_eq!(table.antiproduct(5, 12), (1, 6));
        assert_eq!(table.antiproduct(5, 13), (-1, 7));
        assert_eq!(table.antiproduct(5, 14), (-1, 4));
        assert_eq!(table.antiproduct(5, 15), (1, 5));
        // Row e23
        assert_eq!(table.antiproduct(6, 0), (0, 0));
        assert_eq!(table.antiproduct(6, 1), (0, 0));
        assert_eq!(table.antiproduct(6, 2), (0, 0));
        assert_eq!(table.antiproduct(6, 3), (0, 0));
        assert_eq!(table.antiproduct(6, 4), (0, 0));
        assert_eq!(table.antiproduct(6, 5), (0, 0));
        assert_eq!(table.antiproduct(6, 6), (0, 0));
        assert_eq!(table.antiproduct(6, 7), (0, 0));
        assert_eq!(table.antiproduct(6, 8), (-1, 1));
        assert_eq!(table.antiproduct(6, 9), (1, 0));
        assert_eq!(table.antiproduct(6, 10), (-1, 3));
        assert_eq!(table.antiproduct(6, 11), (1, 2));
        assert_eq!(table.antiproduct(6, 12), (-1, 5));
        assert_eq!(table.antiproduct(6, 13), (1, 4));
        assert_eq!(table.antiproduct(6, 14), (-1, 7));
        assert_eq!(table.antiproduct(6, 15), (1, 6));
        // Row e123
        assert_eq!(table.antiproduct(7, 0), (0, 0));
        assert_eq!(table.antiproduct(7, 1), (0, 0));
        assert_eq!(table.antiproduct(7, 2), (0, 0));
        assert_eq!(table.antiproduct(7, 3), (0, 0));
        assert_eq!(table.antiproduct(7, 4), (0, 0));
        assert_eq!(table.antiproduct(7, 5), (0, 0));
        assert_eq!(table.antiproduct(7, 6), (0, 0));
        assert_eq!(table.antiproduct(7, 7), (0, 0));
        assert_eq!(table.antiproduct(7, 8), (1, 0));
        assert_eq!(table.antiproduct(7, 9), (-1, 1));
        assert_eq!(table.antiproduct(7, 10), (-1, 2));
        assert_eq!(table.antiproduct(7, 11), (1, 3));
        assert_eq!(table.antiproduct(7, 12), (-1, 4));
        assert_eq!(table.antiproduct(7, 13), (1, 5));
        assert_eq!(table.antiproduct(7, 14), (1, 6));
        assert_eq!(table.antiproduct(7, 15), (-1, 7));
        // Row e4
        assert_eq!(table.antiproduct(8, 0), (-1, 7));
        assert_eq!(table.antiproduct(8, 1), (1, 6));
        assert_eq!(table.antiproduct(8, 2), (-1, 5));
        assert_eq!(table.antiproduct(8, 3), (1, 4));
        assert_eq!(table.antiproduct(8, 4), (1, 3));
        assert_eq!(table.antiproduct(8, 5), (-1, 2));
        assert_eq!(table.antiproduct(8, 6), (1, 1));
        assert_eq!(table.antiproduct(8, 7), (-1, 0));
        assert_eq!(table.antiproduct(8, 8), (-1, 15));
        assert_eq!(table.antiproduct(8, 9), (1, 14));
        assert_eq!(table.antiproduct(8, 10), (-1, 13));
        assert_eq!(table.antiproduct(8, 11), (1, 12));
        assert_eq!(table.antiproduct(8, 12), (1, 11));
        assert_eq!(table.antiproduct(8, 13), (-1, 10));
        assert_eq!(table.antiproduct(8, 14), (1, 9));
        assert_eq!(table.antiproduct(8, 15), (-1, 8));
        // Row e14
        assert_eq!(table.antiproduct(9, 0), (-1, 6));
        assert_eq!(table.antiproduct(9, 1), (1, 7));
        assert_eq!(table.antiproduct(9, 2), (1, 4));
        assert_eq!(table.antiproduct(9, 3), (-1, 5));
        assert_eq!(table.antiproduct(9, 4), (-1, 2));
        assert_eq!(table.antiproduct(9, 5), (1, 3));
        assert_eq!(table.antiproduct(9, 6), (1, 0));
        assert_eq!(table.antiproduct(9, 7), (-1, 1));
        assert_eq!(table.antiproduct(9, 8), (1, 14));
        assert_eq!(table.antiproduct(9, 9), (-1, 15));
        assert_eq!(table.antiproduct(9, 10), (-1, 12));
        assert_eq!(table.antiproduct(9, 11), (1, 13));
        assert_eq!(table.antiproduct(9, 12), (1, 10));
        assert_eq!(table.antiproduct(9, 13), (-1, 11));
        assert_eq!(table.antiproduct(9, 14), (-1, 8));
        assert_eq!(table.antiproduct(9, 15), (1, 9));
        // Row e24
        assert_eq!(table.antiproduct(10, 0), (1, 5));
        assert_eq!(table.antiproduct(10, 1), (-1, 4));
        assert_eq!(table.antiproduct(10, 2), (1, 7));
        assert_eq!(table.antiproduct(10, 3), (-1, 6));
        assert_eq!(table.antiproduct(10, 4), (1, 1));
        assert_eq!(table.antiproduct(10, 5), (-1, 0));
        assert_eq!(table.antiproduct(10, 6), (1, 3));
        assert_eq!(table.antiproduct(10, 7), (-1, 2));
        assert_eq!(table.antiproduct(10, 8), (-1, 13));
        assert_eq!(table.antiproduct(10, 9), (1, 12));
        assert_eq!(table.antiproduct(10, 10), (-1, 15));
        assert_eq!(table.antiproduct(10, 11), (1, 14));
        assert_eq!(table.antiproduct(10, 12), (-1, 9));
        assert_eq!(table.antiproduct(10, 13), (1, 8));
        assert_eq!(table.antiproduct(10, 14), (-1, 11));
        assert_eq!(table.antiproduct(10, 15), (1, 10));
        // Row e124
        assert_eq!(table.antiproduct(11, 0), (1, 4));
        assert_eq!(table.antiproduct(11, 1), (-1, 5));
        assert_eq!(table.antiproduct(11, 2), (-1, 6));
        assert_eq!(table.antiproduct(11, 3), (1, 7));
        assert_eq!(table.antiproduct(11, 4), (-1, 0));
        assert_eq!(table.antiproduct(11, 5), (1, 1));
        assert_eq!(table.antiproduct(11, 6), (1, 2));
        assert_eq!(table.antiproduct(11, 7), (-1, 3));
        assert_eq!(table.antiproduct(11, 8), (1, 12));
        assert_eq!(table.antiproduct(11, 9), (-1, 13));
        assert_eq!(table.antiproduct(11, 10), (-1, 14));
        assert_eq!(table.antiproduct(11, 11), (1, 15));
        assert_eq!(table.antiproduct(11, 12), (-1, 8));
        assert_eq!(table.antiproduct(11, 13), (1, 9));
        assert_eq!(table.antiproduct(11, 14), (1, 10));
        assert_eq!(table.antiproduct(11, 15), (-1, 11));
        // Row e34
        assert_eq!(table.antiproduct(12, 0), (-1, 3));
        assert_eq!(table.antiproduct(12, 1), (1, 2));
        assert_eq!(table.antiproduct(12, 2), (-1, 1));
        assert_eq!(table.antiproduct(12, 3), (1, 0));
        assert_eq!(table.antiproduct(12, 4), (1, 7));
        assert_eq!(table.antiproduct(12, 5), (-1, 6));
        assert_eq!(table.antiproduct(12, 6), (1, 5));
        assert_eq!(table.antiproduct(12, 7), (-1, 4));
        assert_eq!(table.antiproduct(12, 8), (1, 11));
        assert_eq!(table.antiproduct(12, 9), (-1, 10));
        assert_eq!(table.antiproduct(12, 10), (1, 9));
        assert_eq!(table.antiproduct(12, 11), (-1, 8));
        assert_eq!(table.antiproduct(12, 12), (-1, 15));
        assert_eq!(table.antiproduct(12, 13), (1, 14));
        assert_eq!(table.antiproduct(12, 14), (-1, 13));
        assert_eq!(table.antiproduct(12, 15), (1, 12));
        // Row e134
        assert_eq!(table.antiproduct(13, 0), (-1, 2));
        assert_eq!(table.antiproduct(13, 1), (1, 3));
        assert_eq!(table.antiproduct(13, 2), (1, 0));
        assert_eq!(table.antiproduct(13, 3), (-1, 1));
        assert_eq!(table.antiproduct(13, 4), (-1, 6));
        assert_eq!(table.antiproduct(13, 5), (1, 7));
        assert_eq!(table.antiproduct(13, 6), (1, 4));
        assert_eq!(table.antiproduct(13, 7), (-1, 5));
        assert_eq!(table.antiproduct(13, 8), (-1, 10));
        assert_eq!(table.antiproduct(13, 9), (1, 11));
        assert_eq!(table.antiproduct(13, 10), (1, 8));
        assert_eq!(table.antiproduct(13, 11), (-1, 9));
        assert_eq!(table.antiproduct(13, 12), (-1, 14));
        assert_eq!(table.antiproduct(13, 13), (1, 15));
        assert_eq!(table.antiproduct(13, 14), (1, 12));
        assert_eq!(table.antiproduct(13, 15), (-1, 13));
        // Row e234
        assert_eq!(table.antiproduct(14, 0), (1, 1));
        assert_eq!(table.antiproduct(14, 1), (-1, 0));
        assert_eq!(table.antiproduct(14, 2), (1, 3));
        assert_eq!(table.antiproduct(14, 3), (-1, 2));
        assert_eq!(table.antiproduct(14, 4), (1, 5));
        assert_eq!(table.antiproduct(14, 5), (-1, 4));
        assert_eq!(table.antiproduct(14, 6), (1, 7));
        assert_eq!(table.antiproduct(14, 7), (-1, 6));
        assert_eq!(table.antiproduct(14, 8), (1, 9));
        assert_eq!(table.antiproduct(14, 9), (-1, 8));
        assert_eq!(table.antiproduct(14, 10), (1, 11));
        assert_eq!(table.antiproduct(14, 11), (-1, 10));
        assert_eq!(table.antiproduct(14, 12), (1, 13));
        assert_eq!(table.antiproduct(14, 13), (-1, 12));
        assert_eq!(table.antiproduct(14, 14), (1, 15));
        assert_eq!(table.antiproduct(14, 15), (-1, 14));
        // Row e1234
        assert_eq!(table.antiproduct(15, 0), (1, 0));
        assert_eq!(table.antiproduct(15, 1), (-1, 1));
        assert_eq!(table.antiproduct(15, 2), (-1, 2));
        assert_eq!(table.antiproduct(15, 3), (1, 3));
        assert_eq!(table.antiproduct(15, 4), (-1, 4));
        assert_eq!(table.antiproduct(15, 5), (1, 5));
        assert_eq!(table.antiproduct(15, 6), (1, 6));
        assert_eq!(table.antiproduct(15, 7), (-1, 7));
        assert_eq!(table.antiproduct(15, 8), (-1, 8));
        assert_eq!(table.antiproduct(15, 9), (1, 9));
        assert_eq!(table.antiproduct(15, 10), (1, 10));
        assert_eq!(table.antiproduct(15, 11), (-1, 11));
        assert_eq!(table.antiproduct(15, 12), (1, 12));
        assert_eq!(table.antiproduct(15, 13), (-1, 13));
        assert_eq!(table.antiproduct(15, 14), (-1, 14));
        assert_eq!(table.antiproduct(15, 15), (1, 15));
    }
}
