# PRD-14.8: Geometric Constraint Engine

**Status**: Draft
**Parent**: [PRD-14: Geometric Algebra Code Generator](prd-14-codegen.md)
**Goal**: Implement the mathematical core for checking geometric constraints that define valid geometric entities

## Overview

According to [Rigid Geometric Algebra](https://rigidgeometricalgebra.org/wiki/index.php?title=Geometric_constraint), geometric entities in GA satisfy two fundamental constraints:

1. **Geometric Product Constraint**: `u ⊡ ũ = u • u`
   - The geometric product of an element with its reverse yields a scalar

2. **Geometric Antiproduct Constraint**: `u ⊟ ũ̃ = u ⊙ u`
   - The antiproduct of an element with its antireverse yields an antiscalar

This PRD implements the mathematical machinery to check whether a given grade combination satisfies these constraints.

## Mathematical Background

### Reverse Operation

The **reverse** of a k-blade flips the order of its basis vectors:

```
ẽ₁₂₃ = e₃₂₁ = (-1)^(3*2/2) e₁₂₃ = -e₁₂₃
```

The sign factor for grade k is:

```
reverse_sign(k) = (-1)^(k(k-1)/2)
```

| Grade | k(k-1)/2 | Sign |
|-------|----------|------|
| 0     | 0        | +1   |
| 1     | 0        | +1   |
| 2     | 1        | -1   |
| 3     | 3        | -1   |
| 4     | 6        | +1   |
| 5     | 10       | +1   |
| 6     | 15       | -1   |
| 7     | 21       | -1   |

Pattern: `++--++--...`

### Antireverse Operation

The **antireverse** is the dual of the reverse (or equivalently, the reverse of the dual). For an n-dimensional algebra, the antireverse of a k-blade has sign:

```
antireverse_sign(k, n) = reverse_sign(n - k)
```

This follows from the duality relationship: the dual of a k-blade is an (n-k)-blade.

### Geometric Antiproduct

The **geometric antiproduct** `a ⊟ b` is the dual of the geometric product of duals:

```
a ⊟ b = dual(dual(a) ⊡ dual(b))
```

For checking constraints, we need to verify that certain grade combinations produce only specific output grades.

## Constraint Checking Algorithm

### Geometric Constraint

For a multivector `u` with grades G = {g₁, g₂, ..., gₘ}, the constraint `u ⊡ ũ = scalar` holds if and only if:

For all pairs (gᵢ, gⱼ) where gᵢ, gⱼ ∈ G:
- The geometric product of a gᵢ-blade with a gⱼ-blade (with reverse sign applied to gⱼ) produces only grade-0 terms

This happens when:
1. gᵢ = gⱼ (same grade produces scalar in dot product)
2. The cross-grade terms cancel out

### Antiproduct Constraint

For the constraint `u ⊟ ũ̃ = antiscalar` (grade n), we check that:

For all pairs (gᵢ, gⱼ) where gᵢ, gⱼ ∈ G:
- The antiproduct of a gᵢ-blade with a gⱼ-blade (with antireverse sign) produces only grade-n terms

## Implementation

### Core Functions

```rust
/// Sign factor for the reverse of a k-blade.
///
/// The reverse flips the order of basis vectors in a blade,
/// introducing a sign of (-1)^(k(k-1)/2).
///
/// # Example
///
/// ```
/// use clifford_codegen::algebra::reverse_sign;
///
/// assert_eq!(reverse_sign(0), 1);  // scalar: +
/// assert_eq!(reverse_sign(1), 1);  // vector: +
/// assert_eq!(reverse_sign(2), -1); // bivector: -
/// assert_eq!(reverse_sign(3), -1); // trivector: -
/// assert_eq!(reverse_sign(4), 1);  // 4-vector: +
/// ```
pub const fn reverse_sign(grade: usize) -> i8 {
    // (-1)^(k(k-1)/2)
    let exponent = (grade * grade.saturating_sub(1)) / 2;
    if exponent % 2 == 0 { 1 } else { -1 }
}

/// Sign factor for the antireverse of a k-blade in an n-dimensional algebra.
///
/// The antireverse is the dual of the reverse, so for a k-blade in n dimensions,
/// the sign is the reverse sign of the dual grade (n - k).
///
/// # Example
///
/// ```
/// use clifford_codegen::algebra::antireverse_sign;
///
/// // In 3D algebra:
/// assert_eq!(antireverse_sign(0, 3), -1); // dual is grade 3
/// assert_eq!(antireverse_sign(1, 3), -1); // dual is grade 2
/// assert_eq!(antireverse_sign(2, 3), 1);  // dual is grade 1
/// assert_eq!(antireverse_sign(3, 3), 1);  // dual is grade 0
/// ```
pub const fn antireverse_sign(grade: usize, dim: usize) -> i8 {
    reverse_sign(dim - grade)
}
```

### Constraint Checking

```rust
use crate::algebra::{Algebra, ProductTable, blades_of_grades, grade};

/// Checks if a grade combination satisfies the geometric constraint.
///
/// The geometric constraint requires that `u ⊡ ũ` produces only a scalar
/// (grade 0). This is checked by verifying that for all blade pairs in the
/// type, the product with reverse gives only grade-0 contributions.
///
/// # Arguments
///
/// * `grades` - The grades present in the type (e.g., `[0, 2]` for a rotor)
/// * `algebra` - The algebra definition
/// * `table` - Precomputed product table
///
/// # Returns
///
/// `true` if the grade combination satisfies the geometric constraint.
pub fn satisfies_geometric_constraint(
    grades: &[usize],
    algebra: &Algebra,
    table: &ProductTable,
) -> bool {
    let blades = blades_of_grades(algebra.dim(), grades);

    // Check that u * reverse(u) produces only grade 0
    for &a in &blades {
        for &b in &blades {
            let rev_sign = reverse_sign(grade(b));
            let (prod_sign, result) = table.geometric(a, b);

            // If this product contributes with non-zero sign...
            if prod_sign != 0 && rev_sign != 0 {
                let result_grade = grade(result);
                // ...it must be grade 0
                if result_grade != 0 {
                    return false;
                }
            }
        }
    }

    true
}

/// Checks if a grade combination satisfies the antiproduct constraint.
///
/// The antiproduct constraint requires that `u ⊟ ũ̃` produces only an
/// antiscalar (grade n). This is checked using the duality relationship.
///
/// # Arguments
///
/// * `grades` - The grades present in the type
/// * `algebra` - The algebra definition
/// * `table` - Precomputed product table
///
/// # Returns
///
/// `true` if the grade combination satisfies the antiproduct constraint.
pub fn satisfies_antiproduct_constraint(
    grades: &[usize],
    algebra: &Algebra,
    table: &ProductTable,
) -> bool {
    let dim = algebra.dim();
    let blades = blades_of_grades(dim, grades);

    // For antiproduct: a ⊟ b = dual(dual(a) * dual(b))
    // The dual of a k-blade is an (n-k)-blade
    // We need the result to be grade n (antiscalar)

    for &a in &blades {
        for &b in &blades {
            // Compute dual blade indices
            let dual_a = (1 << dim) - 1 - a; // XOR with pseudoscalar
            let dual_b = (1 << dim) - 1 - b;

            // Antireverse of b
            let antirev_sign = antireverse_sign(grade(b), dim);

            // Product of duals
            let (prod_sign, dual_result) = table.geometric(dual_a, dual_b);

            if prod_sign != 0 && antirev_sign != 0 {
                // Dual of the result (back to original grading)
                let result = (1 << dim) - 1 - dual_result;
                let result_grade = grade(result);

                // Must be grade n (antiscalar = pseudoscalar)
                if result_grade != dim {
                    return false;
                }
            }
        }
    }

    true
}

/// Checks if a grade combination satisfies both geometric constraints.
///
/// A valid geometric entity must satisfy:
/// 1. `u ⊡ ũ = u • u` (scalar)
/// 2. `u ⊟ ũ̃ = u ⊙ u` (antiscalar)
pub fn satisfies_all_constraints(
    grades: &[usize],
    algebra: &Algebra,
    table: &ProductTable,
) -> bool {
    satisfies_geometric_constraint(grades, algebra, table)
        && satisfies_antiproduct_constraint(grades, algebra, table)
}
```

## Module Structure

```
crates/clifford-codegen/src/algebra/
├── mod.rs           # Add exports for new functions
├── grade.rs         # Add reverse_sign, antireverse_sign
└── constraints.rs   # NEW: satisfies_geometric_constraint, etc.
```

## Testing

### Property-Based Tests

```rust
#[cfg(test)]
mod tests {
    use super::*;
    use proptest::prelude::*;

    #[test]
    fn reverse_sign_pattern() {
        // Pattern: ++--++--...
        assert_eq!(reverse_sign(0), 1);
        assert_eq!(reverse_sign(1), 1);
        assert_eq!(reverse_sign(2), -1);
        assert_eq!(reverse_sign(3), -1);
        assert_eq!(reverse_sign(4), 1);
        assert_eq!(reverse_sign(5), 1);
        assert_eq!(reverse_sign(6), -1);
        assert_eq!(reverse_sign(7), -1);
    }

    #[test]
    fn reverse_is_involution() {
        // reverse(reverse(x)) = x, so sign^2 = 1
        for k in 0..10 {
            let sign = reverse_sign(k);
            assert!(sign == 1 || sign == -1);
        }
    }

    proptest! {
        #[test]
        fn antireverse_is_dual_of_reverse(k in 0usize..8, n in 0usize..8) {
            prop_assume!(k <= n);
            let dual_grade = n - k;
            prop_assert_eq!(
                antireverse_sign(k, n),
                reverse_sign(dual_grade)
            );
        }
    }

    #[test]
    fn euclidean_3d_known_entities() {
        let algebra = Algebra::euclidean(3);
        let table = ProductTable::new(&algebra);

        // Scalars satisfy constraints
        assert!(satisfies_all_constraints(&[0], &algebra, &table));

        // Vectors satisfy constraints
        assert!(satisfies_all_constraints(&[1], &algebra, &table));

        // Bivectors satisfy constraints
        assert!(satisfies_all_constraints(&[2], &algebra, &table));

        // Trivector (pseudoscalar) satisfies constraints
        assert!(satisfies_all_constraints(&[3], &algebra, &table));

        // Even subalgebra (rotor) satisfies constraints
        assert!(satisfies_all_constraints(&[0, 2], &algebra, &table));

        // Odd subalgebra satisfies constraints
        assert!(satisfies_all_constraints(&[1, 3], &algebra, &table));

        // Full multivector satisfies constraints
        assert!(satisfies_all_constraints(&[0, 1, 2, 3], &algebra, &table));
    }

    #[test]
    fn pga_3d_known_entities() {
        let algebra = Algebra::pga(3); // 4D with one degenerate basis
        let table = ProductTable::new(&algebra);

        // Motors (grade 0 + grade 2) should satisfy constraints
        assert!(satisfies_all_constraints(&[0, 2], &algebra, &table));

        // Flectors (grade 1 + grade 3) should satisfy constraints
        assert!(satisfies_all_constraints(&[1, 3], &algebra, &table));
    }
}
```

### Verification Against Known Algebras

Test that discovered entities match expectations for:
- Euclidean 2D and 3D
- PGA 2D and 3D
- CGA 2D and 3D

## Deliverables

- [ ] `reverse_sign(grade: usize) -> i8` function in `grade.rs`
- [ ] `antireverse_sign(grade: usize, dim: usize) -> i8` function in `grade.rs`
- [ ] `constraints.rs` module with constraint checking functions
- [ ] `satisfies_geometric_constraint(grades, algebra, table) -> bool`
- [ ] `satisfies_antiproduct_constraint(grades, algebra, table) -> bool`
- [ ] `satisfies_all_constraints(grades, algebra, table) -> bool`
- [ ] Property-based tests for sign functions
- [ ] Verification tests against known algebras

## Success Criteria

1. **Correctness**: Sign functions match mathematical definitions
2. **Completeness**: Handles all algebra types (Euclidean, PGA, CGA)
3. **Performance**: Constraint checking is efficient (O(n²) in number of blades)
4. **Testing**: All known entity types pass constraint checks

## Dependencies

- PRD-14.1 (Blade Algebra Engine) - uses `Algebra`, `ProductTable`, grade functions

## Next Steps

- PRD-14.9 uses this engine to discover valid grade combinations
- PRD-14.10 uses discovered entities for product inference

## References

- [Geometric Constraint - RGA Wiki](https://rigidgeometricalgebra.org/wiki/index.php?title=Geometric_constraint)
- [Geometric Antiproduct - RGA Wiki](https://rigidgeometricalgebra.org/wiki/index.php?title=Geometric_antiproduct)
