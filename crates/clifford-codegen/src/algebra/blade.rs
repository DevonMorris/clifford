//! Blade representation and operations.
//!
//! A blade is a basis element in a geometric algebra, represented by its
//! bitmask index where bit `i` indicates the presence of basis vector `eᵢ`.

use std::fmt;

/// Maximum supported dimension.
///
/// Limits the algebra to at most 6 basis vectors (64 blades).
pub const MAX_DIM: usize = 6;

/// A basis blade in a geometric algebra.
///
/// Blades are represented as bitmasks where bit `i` indicates the presence
/// of basis vector `eᵢ`. This representation is always canonical since
/// bits are inherently ordered.
///
/// # Examples
///
/// ```
/// use clifford_codegen::algebra::Blade;
///
/// // Scalar (grade 0)
/// let scalar = Blade::scalar();
/// assert_eq!(scalar.grade(), 0);
/// assert_eq!(scalar.index(), 0);
///
/// // Basis vector e₁
/// let e1 = Blade::basis(0);
/// assert_eq!(e1.grade(), 1);
/// assert_eq!(e1.index(), 1);
///
/// // Bivector e₁₂
/// let e12 = Blade::from_index(0b11);
/// assert_eq!(e12.grade(), 2);
/// ```
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub struct Blade {
    /// Bitmask representation where bit `i` = 1 means `eᵢ` is present.
    index: usize,
}

impl Blade {
    /// Creates a blade from its bitmask index.
    ///
    /// # Example
    ///
    /// ```
    /// use clifford_codegen::algebra::Blade;
    ///
    /// let e12 = Blade::from_index(0b11);
    /// assert_eq!(e12.index(), 3);
    /// ```
    #[inline]
    pub const fn from_index(index: usize) -> Self {
        Self { index }
    }

    /// Creates the scalar blade (grade 0, index 0).
    ///
    /// # Example
    ///
    /// ```
    /// use clifford_codegen::algebra::Blade;
    ///
    /// let s = Blade::scalar();
    /// assert_eq!(s.index(), 0);
    /// assert_eq!(s.grade(), 0);
    /// ```
    #[inline]
    pub const fn scalar() -> Self {
        Self { index: 0 }
    }

    /// Creates a basis vector blade.
    ///
    /// `Blade::basis(i)` creates blade `eᵢ₊₁` (using 0-based indexing internally).
    ///
    /// # Example
    ///
    /// ```
    /// use clifford_codegen::algebra::Blade;
    ///
    /// let e1 = Blade::basis(0);
    /// assert_eq!(e1.index(), 1);
    ///
    /// let e2 = Blade::basis(1);
    /// assert_eq!(e2.index(), 2);
    ///
    /// let e3 = Blade::basis(2);
    /// assert_eq!(e3.index(), 4);
    /// ```
    #[inline]
    pub const fn basis(i: usize) -> Self {
        Self { index: 1 << i }
    }

    /// Returns the blade's bitmask index.
    ///
    /// # Example
    ///
    /// ```
    /// use clifford_codegen::algebra::Blade;
    ///
    /// let e12 = Blade::from_index(3);
    /// assert_eq!(e12.index(), 3);
    /// ```
    #[inline]
    pub const fn index(&self) -> usize {
        self.index
    }

    /// Returns the grade (number of basis vectors in this blade).
    ///
    /// The grade equals the number of 1-bits in the index.
    ///
    /// # Example
    ///
    /// ```
    /// use clifford_codegen::algebra::Blade;
    ///
    /// assert_eq!(Blade::scalar().grade(), 0);
    /// assert_eq!(Blade::basis(0).grade(), 1);
    /// assert_eq!(Blade::from_index(0b111).grade(), 3);
    /// ```
    #[inline]
    pub const fn grade(&self) -> usize {
        self.index.count_ones() as usize
    }

    /// Checks if this blade contains basis vector `i`.
    ///
    /// # Example
    ///
    /// ```
    /// use clifford_codegen::algebra::Blade;
    ///
    /// let e12 = Blade::from_index(0b11);
    /// assert!(e12.contains(0)); // contains e1
    /// assert!(e12.contains(1)); // contains e2
    /// assert!(!e12.contains(2)); // does not contain e3
    /// ```
    #[inline]
    pub const fn contains(&self, i: usize) -> bool {
        (self.index >> i) & 1 == 1
    }

    /// Returns an iterator over the basis vector indices in this blade.
    ///
    /// Indices are yielded in ascending order.
    ///
    /// # Example
    ///
    /// ```
    /// use clifford_codegen::algebra::Blade;
    ///
    /// let e135 = Blade::from_index(0b10101);
    /// let indices: Vec<_> = e135.basis_vectors().collect();
    /// assert_eq!(indices, vec![0, 2, 4]);
    /// ```
    pub fn basis_vectors(&self) -> impl Iterator<Item = usize> + '_ {
        (0..MAX_DIM).filter(move |&i| self.contains(i))
    }

    /// Returns the number of basis vectors in this blade.
    ///
    /// This is equivalent to `self.grade()`.
    #[inline]
    pub const fn len(&self) -> usize {
        self.grade()
    }

    /// Returns true if this is the scalar blade.
    #[inline]
    pub const fn is_empty(&self) -> bool {
        self.index == 0
    }

    /// Computes the XOR of two blade indices.
    ///
    /// This gives the result blade of the geometric product (before
    /// considering sign).
    ///
    /// # Example
    ///
    /// ```
    /// use clifford_codegen::algebra::Blade;
    ///
    /// let e1 = Blade::basis(0);
    /// let e2 = Blade::basis(1);
    /// let e12 = e1.xor(e2);
    /// assert_eq!(e12.index(), 3);
    /// ```
    #[inline]
    pub const fn xor(&self, other: Self) -> Self {
        Self {
            index: self.index ^ other.index,
        }
    }

    /// Returns true if this blade shares any basis vectors with another.
    ///
    /// # Example
    ///
    /// ```
    /// use clifford_codegen::algebra::Blade;
    ///
    /// let e12 = Blade::from_index(0b11);
    /// let e23 = Blade::from_index(0b110);
    /// let e45 = Blade::from_index(0b11000);
    ///
    /// assert!(e12.overlaps(e23)); // share e2
    /// assert!(!e12.overlaps(e45)); // disjoint
    /// ```
    #[inline]
    pub const fn overlaps(&self, other: Self) -> bool {
        (self.index & other.index) != 0
    }

    /// Returns the blade name using standard notation.
    ///
    /// Uses 1-based indexing for display: `e1`, `e12`, `e123`, etc.
    ///
    /// # Example
    ///
    /// ```
    /// use clifford_codegen::algebra::Blade;
    ///
    /// assert_eq!(Blade::scalar().name(), "1");
    /// assert_eq!(Blade::basis(0).name(), "e1");
    /// assert_eq!(Blade::from_index(0b11).name(), "e12");
    /// ```
    pub fn name(&self) -> String {
        if self.index == 0 {
            return "1".to_string();
        }

        let mut name = String::from("e");
        for i in self.basis_vectors() {
            name.push_str(&(i + 1).to_string());
        }
        name
    }
}

impl fmt::Display for Blade {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}", self.name())
    }
}

impl Default for Blade {
    fn default() -> Self {
        Self::scalar()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn scalar_properties() {
        let s = Blade::scalar();
        assert_eq!(s.index(), 0);
        assert_eq!(s.grade(), 0);
        assert!(s.is_empty());
        assert_eq!(s.name(), "1");
    }

    #[test]
    fn basis_vector_properties() {
        for i in 0..MAX_DIM {
            let e = Blade::basis(i);
            assert_eq!(e.index(), 1 << i);
            assert_eq!(e.grade(), 1);
            assert!(e.contains(i));
            assert!(!e.is_empty());

            for j in 0..MAX_DIM {
                if j != i {
                    assert!(!e.contains(j));
                }
            }
        }
    }

    #[test]
    fn bivector_properties() {
        let e12 = Blade::from_index(0b11);
        assert_eq!(e12.grade(), 2);
        assert!(e12.contains(0));
        assert!(e12.contains(1));
        assert!(!e12.contains(2));
        assert_eq!(e12.name(), "e12");
    }

    #[test]
    fn xor_computes_product_blade() {
        let e1 = Blade::basis(0);
        let e2 = Blade::basis(1);
        let e12 = e1.xor(e2);
        assert_eq!(e12.index(), 0b11);

        // e12 * e2 = e1
        let e1_back = e12.xor(e2);
        assert_eq!(e1_back.index(), e1.index());
    }

    #[test]
    fn basis_vectors_iterator() {
        let e135 = Blade::from_index(0b10101);
        let indices: Vec<_> = e135.basis_vectors().collect();
        assert_eq!(indices, vec![0, 2, 4]);
    }

    #[test]
    fn overlaps_detection() {
        let e12 = Blade::from_index(0b11);
        let e23 = Blade::from_index(0b110);
        let e34 = Blade::from_index(0b1100);

        assert!(e12.overlaps(e23)); // share e2
        assert!(e23.overlaps(e34)); // share e3
        assert!(!e12.overlaps(e34)); // disjoint
    }

    #[test]
    fn ordering() {
        let blades: Vec<Blade> = (0..8).map(Blade::from_index).collect();

        // Verify they sort by index
        let mut sorted = blades.clone();
        sorted.sort();
        assert_eq!(blades, sorted);
    }
}
