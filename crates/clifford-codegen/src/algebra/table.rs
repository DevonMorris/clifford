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
}
