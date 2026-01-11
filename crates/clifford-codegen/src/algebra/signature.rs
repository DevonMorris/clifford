//! Algebra type representing a geometric algebra with a specific metric signature.
//!
//! A geometric algebra is defined by its **metric signature** (p, q, r):
//! - `p` basis vectors that square to `+1` (Euclidean)
//! - `q` basis vectors that square to `-1` (anti-Euclidean/Minkowski)
//! - `r` basis vectors that square to `0` (degenerate/null)
//!
//! Common algebras:
//! - **Euclidean 3D**: (3, 0, 0) - standard 3D geometry
//! - **PGA 3D**: (3, 0, 1) - projective geometric algebra with e₀² = 0
//! - **CGA 3D**: (4, 1, 0) - conformal geometric algebra with e₊² = 1, e₋² = -1
//! - **Minkowski**: (3, 1, 0) - spacetime algebra

use std::collections::HashMap;

use super::blade::Blade;
use super::grade::blades_of_grade;
use super::sign::basis_product;

/// A geometric algebra defined by its metric signature.
///
/// The algebra encapsulates:
/// - The metric signature (p, q, r)
/// - Optional custom names for basis vectors and blades
///
/// # Example
///
/// ```
/// use clifford_codegen::algebra::Algebra;
///
/// // Create a 3D Euclidean algebra
/// let alg = Algebra::euclidean(3);
/// assert_eq!(alg.dim(), 3);
/// assert_eq!(alg.num_blades(), 8);
///
/// // Metric: all basis vectors square to +1
/// assert_eq!(alg.metric(0), 1);
/// assert_eq!(alg.metric(1), 1);
/// assert_eq!(alg.metric(2), 1);
/// ```
#[derive(Clone, Debug)]
pub struct Algebra {
    /// Number of basis vectors squaring to +1.
    p: usize,
    /// Number of basis vectors squaring to -1.
    q: usize,
    /// Number of basis vectors squaring to 0.
    r: usize,
    /// Custom names for basis vectors (1-indexed in display).
    basis_names: Vec<String>,
    /// Custom names for blades (blade index -> name).
    blade_names: HashMap<usize, String>,
}

impl Algebra {
    /// Creates a new algebra with signature (p, q, r).
    ///
    /// Basis vectors are ordered: first p positive, then q negative, then r null.
    ///
    /// # Example
    ///
    /// ```
    /// use clifford_codegen::algebra::Algebra;
    ///
    /// // PGA: 3 Euclidean + 1 degenerate
    /// let pga = Algebra::new(3, 0, 1);
    /// assert_eq!(pga.dim(), 4);
    /// assert_eq!(pga.metric(0), 1);  // e1 squares to +1
    /// assert_eq!(pga.metric(3), 0);  // e4 squares to 0
    /// ```
    pub fn new(p: usize, q: usize, r: usize) -> Self {
        let dim = p + q + r;
        let basis_names = (1..=dim).map(|i| format!("e{}", i)).collect();
        Self {
            p,
            q,
            r,
            basis_names,
            blade_names: HashMap::new(),
        }
    }

    /// Creates a Euclidean algebra of dimension n.
    ///
    /// All basis vectors square to +1.
    ///
    /// # Example
    ///
    /// ```
    /// use clifford_codegen::algebra::Algebra;
    ///
    /// let e3 = Algebra::euclidean(3);
    /// assert_eq!(e3.signature(), (3, 0, 0));
    /// ```
    pub fn euclidean(n: usize) -> Self {
        Self::new(n, 0, 0)
    }

    /// Creates a Projective Geometric Algebra for n-dimensional space.
    ///
    /// PGA(n) has signature (n, 0, 1) where the extra degenerate basis
    /// represents the point at infinity.
    ///
    /// # Example
    ///
    /// ```
    /// use clifford_codegen::algebra::Algebra;
    ///
    /// let pga3 = Algebra::pga(3);
    /// assert_eq!(pga3.dim(), 4);
    /// assert_eq!(pga3.signature(), (3, 0, 1));
    /// ```
    pub fn pga(n: usize) -> Self {
        Self::new(n, 0, 1)
    }

    /// Creates a Conformal Geometric Algebra for n-dimensional space.
    ///
    /// CGA(n) has signature (n+1, 1, 0), adding two extra basis vectors
    /// e₊ (positive) and e₋ (negative) for the conformal model.
    ///
    /// # Example
    ///
    /// ```
    /// use clifford_codegen::algebra::Algebra;
    ///
    /// let cga3 = Algebra::cga(3);
    /// assert_eq!(cga3.dim(), 5);
    /// assert_eq!(cga3.signature(), (4, 1, 0));
    /// ```
    pub fn cga(n: usize) -> Self {
        Self::new(n + 1, 1, 0)
    }

    /// Creates a Minkowski algebra with n spatial dimensions.
    ///
    /// Minkowski(n) has signature (n, 1, 0), with n spatial dimensions
    /// and 1 time-like dimension (squares to -1).
    ///
    /// # Example
    ///
    /// ```
    /// use clifford_codegen::algebra::Algebra;
    ///
    /// let minkowski = Algebra::minkowski(3);
    /// assert_eq!(minkowski.dim(), 4);
    /// assert_eq!(minkowski.signature(), (3, 1, 0));
    /// ```
    pub fn minkowski(n: usize) -> Self {
        Self::new(n, 1, 0)
    }

    /// Returns the total dimension (number of basis vectors).
    #[inline]
    pub fn dim(&self) -> usize {
        self.p + self.q + self.r
    }

    /// Returns the signature (p, q, r).
    #[inline]
    pub fn signature(&self) -> (usize, usize, usize) {
        (self.p, self.q, self.r)
    }

    /// Returns the total number of blades (2^dim).
    #[inline]
    pub fn num_blades(&self) -> usize {
        1 << self.dim()
    }

    /// Returns the metric value for basis vector i.
    ///
    /// - Returns `+1` for indices `0..p` (positive)
    /// - Returns `-1` for indices `p..p+q` (negative)
    /// - Returns `0` for indices `p+q..p+q+r` (null/degenerate)
    ///
    /// # Example
    ///
    /// ```
    /// use clifford_codegen::algebra::Algebra;
    ///
    /// let cga = Algebra::cga(3);  // signature (4, 1, 0)
    /// assert_eq!(cga.metric(0), 1);   // e1 (positive)
    /// assert_eq!(cga.metric(3), 1);   // e+ (positive)
    /// assert_eq!(cga.metric(4), -1);  // e- (negative)
    /// ```
    #[inline]
    pub fn metric(&self, i: usize) -> i8 {
        if i < self.p {
            1
        } else if i < self.p + self.q {
            -1
        } else {
            0
        }
    }

    /// Computes the product of two basis blades.
    ///
    /// # Returns
    ///
    /// A tuple `(sign, result)` where sign ∈ {-1, 0, +1} and result
    /// is the blade index of the product.
    ///
    /// # Example
    ///
    /// ```
    /// use clifford_codegen::algebra::Algebra;
    ///
    /// let alg = Algebra::euclidean(3);
    ///
    /// // e1 * e2 = e12
    /// let (sign, result) = alg.basis_product(0b001, 0b010);
    /// assert_eq!(sign, 1);
    /// assert_eq!(result, 0b011);
    ///
    /// // e2 * e1 = -e12
    /// let (sign, result) = alg.basis_product(0b010, 0b001);
    /// assert_eq!(sign, -1);
    /// assert_eq!(result, 0b011);
    /// ```
    pub fn basis_product(&self, a: usize, b: usize) -> (i8, usize) {
        basis_product(a, b, |i| self.metric(i))
    }

    /// Returns all blades of a given grade.
    ///
    /// # Example
    ///
    /// ```
    /// use clifford_codegen::algebra::Algebra;
    ///
    /// let alg = Algebra::euclidean(3);
    /// let bivectors = alg.blades_of_grade(2);
    ///
    /// assert_eq!(bivectors.len(), 3);
    /// assert_eq!(bivectors[0].index(), 0b011); // e12
    /// assert_eq!(bivectors[1].index(), 0b101); // e13
    /// assert_eq!(bivectors[2].index(), 0b110); // e23
    /// ```
    pub fn blades_of_grade(&self, grade: usize) -> Vec<Blade> {
        blades_of_grade(self.dim(), grade)
            .into_iter()
            .map(Blade::from_index)
            .collect()
    }

    /// Returns all blades in the algebra.
    ///
    /// Blades are returned in index order (canonical ordering).
    pub fn all_blades(&self) -> Vec<Blade> {
        (0..self.num_blades()).map(Blade::from_index).collect()
    }

    /// Sets a custom name for a basis vector.
    ///
    /// # Arguments
    ///
    /// * `i` - The basis vector index (0-based)
    /// * `name` - The custom name
    ///
    /// # Example
    ///
    /// ```
    /// use clifford_codegen::algebra::Algebra;
    ///
    /// let mut alg = Algebra::euclidean(3);
    /// alg.set_basis_name(0, "x".to_string());
    /// alg.set_basis_name(1, "y".to_string());
    /// alg.set_basis_name(2, "z".to_string());
    ///
    /// assert_eq!(alg.basis_name(0), "x");
    /// assert_eq!(alg.basis_name(1), "y");
    /// assert_eq!(alg.basis_name(2), "z");
    /// ```
    pub fn set_basis_name(&mut self, i: usize, name: String) {
        if i < self.basis_names.len() {
            self.basis_names[i] = name;
        }
    }

    /// Returns the name of a basis vector.
    pub fn basis_name(&self, i: usize) -> &str {
        &self.basis_names[i]
    }

    /// Sets a custom name for a blade.
    ///
    /// # Example
    ///
    /// ```
    /// use clifford_codegen::algebra::{Algebra, Blade};
    ///
    /// let mut alg = Algebra::euclidean(3);
    /// alg.set_blade_name(Blade::from_index(0b011), "xy".to_string());
    ///
    /// assert_eq!(alg.blade_name(Blade::from_index(0b011)), "xy");
    /// ```
    pub fn set_blade_name(&mut self, blade: Blade, name: String) {
        self.blade_names.insert(blade.index(), name);
    }

    /// Returns the name for a blade.
    ///
    /// If a custom name was set, returns that. Otherwise, builds
    /// the name from basis vector names.
    ///
    /// # Example
    ///
    /// ```
    /// use clifford_codegen::algebra::{Algebra, Blade};
    ///
    /// let alg = Algebra::euclidean(3);
    ///
    /// // Default names
    /// assert_eq!(alg.blade_name(Blade::scalar()), "s");
    /// assert_eq!(alg.blade_name(Blade::basis(0)), "e1");
    /// assert_eq!(alg.blade_name(Blade::from_index(0b011)), "e1e2");
    /// ```
    pub fn blade_name(&self, blade: Blade) -> String {
        if let Some(name) = self.blade_names.get(&blade.index()) {
            return name.clone();
        }

        if blade.index() == 0 {
            return "s".to_string();
        }

        // Build from basis names
        blade
            .basis_vectors()
            .map(|i| self.basis_names[i].as_str())
            .collect::<Vec<_>>()
            .join("")
    }

    /// Returns the canonical index-based name for a blade.
    ///
    /// This produces names compatible with the TOML parser format:
    /// - Scalar: "s" (but scalars shouldn't be in blades section)
    /// - Basis vectors: "e1", "e2", etc. (1-indexed)
    /// - Higher blades: "e12", "e123", etc.
    ///
    /// # Example
    ///
    /// ```
    /// use clifford_codegen::algebra::{Algebra, Blade};
    ///
    /// let alg = Algebra::euclidean(3);
    ///
    /// assert_eq!(alg.blade_index_name(Blade::basis(0)), "e1");
    /// assert_eq!(alg.blade_index_name(Blade::basis(1)), "e2");
    /// assert_eq!(alg.blade_index_name(Blade::from_index(0b011)), "e12");
    /// assert_eq!(alg.blade_index_name(Blade::from_index(0b111)), "e123");
    /// ```
    pub fn blade_index_name(&self, blade: Blade) -> String {
        if blade.index() == 0 {
            return "s".to_string();
        }

        // Build from 1-indexed basis indices
        let indices: String = blade
            .basis_vectors()
            .map(|i| char::from_digit((i + 1) as u32, 10).unwrap())
            .collect();

        format!("e{}", indices)
    }
}

impl Default for Algebra {
    fn default() -> Self {
        Self::euclidean(3)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn euclidean_signature() {
        let alg = Algebra::euclidean(3);
        assert_eq!(alg.signature(), (3, 0, 0));
        assert_eq!(alg.dim(), 3);
        assert_eq!(alg.num_blades(), 8);

        for i in 0..3 {
            assert_eq!(alg.metric(i), 1);
        }
    }

    #[test]
    fn pga_signature() {
        let alg = Algebra::pga(3);
        assert_eq!(alg.signature(), (3, 0, 1));
        assert_eq!(alg.dim(), 4);
        assert_eq!(alg.num_blades(), 16);

        // First 3 square to +1, last one to 0
        assert_eq!(alg.metric(0), 1);
        assert_eq!(alg.metric(1), 1);
        assert_eq!(alg.metric(2), 1);
        assert_eq!(alg.metric(3), 0);
    }

    #[test]
    fn cga_signature() {
        let alg = Algebra::cga(3);
        assert_eq!(alg.signature(), (4, 1, 0));
        assert_eq!(alg.dim(), 5);
        assert_eq!(alg.num_blades(), 32);

        // First 4 square to +1, last one to -1
        assert_eq!(alg.metric(0), 1);
        assert_eq!(alg.metric(1), 1);
        assert_eq!(alg.metric(2), 1);
        assert_eq!(alg.metric(3), 1);
        assert_eq!(alg.metric(4), -1);
    }

    #[test]
    fn minkowski_signature() {
        let alg = Algebra::minkowski(3);
        assert_eq!(alg.signature(), (3, 1, 0));
        assert_eq!(alg.dim(), 4);

        assert_eq!(alg.metric(0), 1);
        assert_eq!(alg.metric(1), 1);
        assert_eq!(alg.metric(2), 1);
        assert_eq!(alg.metric(3), -1);
    }

    #[test]
    fn blades_by_grade() {
        let alg = Algebra::euclidean(3);

        let scalars = alg.blades_of_grade(0);
        assert_eq!(scalars.len(), 1);
        assert_eq!(scalars[0].index(), 0);

        let vectors = alg.blades_of_grade(1);
        assert_eq!(vectors.len(), 3);
        assert_eq!(vectors[0].index(), 1);
        assert_eq!(vectors[1].index(), 2);
        assert_eq!(vectors[2].index(), 4);

        let bivectors = alg.blades_of_grade(2);
        assert_eq!(bivectors.len(), 3);
        assert_eq!(bivectors[0].index(), 3);
        assert_eq!(bivectors[1].index(), 5);
        assert_eq!(bivectors[2].index(), 6);
    }

    #[test]
    fn custom_names() {
        let mut alg = Algebra::euclidean(3);
        alg.set_basis_name(0, "x".to_string());
        alg.set_basis_name(1, "y".to_string());
        alg.set_basis_name(2, "z".to_string());

        assert_eq!(alg.blade_name(Blade::basis(0)), "x");
        assert_eq!(alg.blade_name(Blade::basis(1)), "y");
        assert_eq!(alg.blade_name(Blade::basis(2)), "z");
        assert_eq!(alg.blade_name(Blade::from_index(0b011)), "xy");
        assert_eq!(alg.blade_name(Blade::from_index(0b111)), "xyz");
    }

    #[test]
    fn product_euclidean() {
        let alg = Algebra::euclidean(3);

        // e1 * e2 = e12
        let (sign, result) = alg.basis_product(1, 2);
        assert_eq!(sign, 1);
        assert_eq!(result, 3);

        // e2 * e1 = -e12
        let (sign, result) = alg.basis_product(2, 1);
        assert_eq!(sign, -1);
        assert_eq!(result, 3);

        // e1 * e1 = 1
        let (sign, result) = alg.basis_product(1, 1);
        assert_eq!(sign, 1);
        assert_eq!(result, 0);
    }

    #[test]
    fn product_pga_degenerate() {
        let alg = Algebra::pga(3);

        // e4 * e4 = 0 (degenerate)
        let (sign, result) = alg.basis_product(8, 8);
        assert_eq!(sign, 0);
        assert_eq!(result, 0);

        // e1 * e1 = 1 (non-degenerate)
        let (sign, result) = alg.basis_product(1, 1);
        assert_eq!(sign, 1);
        assert_eq!(result, 0);
    }

    #[test]
    fn product_minkowski() {
        let alg = Algebra::minkowski(3);

        // e4 * e4 = -1 (time-like)
        let (sign, result) = alg.basis_product(8, 8);
        assert_eq!(sign, -1);
        assert_eq!(result, 0);

        // e1 * e1 = +1 (space-like)
        let (sign, result) = alg.basis_product(1, 1);
        assert_eq!(sign, 1);
        assert_eq!(result, 0);
    }
}
