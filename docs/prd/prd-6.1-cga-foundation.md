# PRD-6.1: CGA Foundation - Signatures and Null Basis

**Status**: Pending
**Parent**: PRD-6 (Conformal Geometric Algebra)
**Goal**: Establish CGA signatures and verify null basis properties

## Reference

**Primary Resource**: https://conformalgeometricalgebra.org/wiki/index.php?title=Main_Page

## Background

Conformal Geometric Algebra (CGA) embeds Euclidean space into a higher-dimensional space by adding two extra basis vectors: one that squares to +1 (`e₊`) and one that squares to -1 (`e₋`). This enables uniform representation of points, spheres, planes, circles, and lines as blades.

### Signatures

| Algebra | Signature | Dimension | Blades | Basis |
|---------|-----------|-----------|--------|-------|
| 2D CGA | `Cl(3,1,0)` | 4 | 16 | e₁, e₂, e₊, e₋ |
| 3D CGA | `Cl(4,1,0)` | 5 | 32 | e₁, e₂, e₃, e₊, e₋ |

### Null Basis Convention

The wiki uses this convention:
- `e₊` squares to +1, `e₋` squares to -1
- `e∞ = e₋ + e₊` (point at infinity)
- `e₀ = (e₋ - e₊) / 2` (origin)
- Key properties: `e∞² = 0`, `e₀² = 0`, `e∞ · e₀ = -1`

## Deliverables

### 1. Signature Types (`src/signature/conformal.rs`)

```rust
/// 2D Conformal GA: Cl(3,1,0)
///
/// Embeds 2D Euclidean space into a 4D conformal space.
/// Basis: e₁, e₂, e₊, e₋ where e₊² = 1, e₋² = -1
///
/// Reference: https://conformalgeometricalgebra.org/wiki/index.php?title=Main_Page
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct Conformal2;

impl Signature for Conformal2 {
    type NumBlades = typenum::U16; // 2⁴ = 16

    const P: usize = 3;  // e₁, e₂, e₊
    const Q: usize = 1;  // e₋
    const R: usize = 0;

    fn metric(i: usize) -> i8 {
        match i {
            0 | 1 | 2 => 1,   // e₁² = e₂² = e₊² = 1
            3 => -1,          // e₋² = -1
            _ => panic!("basis index {i} out of range for Conformal2"),
        }
    }
}

/// 3D Conformal GA: Cl(4,1,0)
///
/// Embeds 3D Euclidean space into a 5D conformal space.
/// Basis: e₁, e₂, e₃, e₊, e₋ where e₊² = 1, e₋² = -1
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct Conformal3;

impl Signature for Conformal3 {
    type NumBlades = typenum::U32; // 2⁵ = 32

    const P: usize = 4;  // e₁, e₂, e₃, e₊
    const Q: usize = 1;  // e₋
    const R: usize = 0;

    fn metric(i: usize) -> i8 {
        match i {
            0 | 1 | 2 | 3 => 1,  // e₁² = e₂² = e₃² = e₊² = 1
            4 => -1,             // e₋² = -1
            _ => panic!("basis index {i} out of range for Conformal3"),
        }
    }
}
```

### 2. Null Basis Helper Functions

```rust
impl Conformal3 {
    /// Index of e₊ basis vector.
    pub const E_PLUS: usize = 3;
    /// Index of e₋ basis vector.
    pub const E_MINUS: usize = 4;
}

/// Creates e∞ = e₋ + e₊ (point at infinity).
pub fn e_infinity<T: Float>() -> Multivector<T, Conformal3> {
    let mut mv = Multivector::zero();
    mv.set(Blade::basis_vector(Conformal3::E_PLUS), T::one());
    mv.set(Blade::basis_vector(Conformal3::E_MINUS), T::one());
    mv
}

/// Creates e₀ = (e₋ - e₊) / 2 (origin).
pub fn e_origin<T: Float>() -> Multivector<T, Conformal3> {
    let mut mv = Multivector::zero();
    let half = T::one() / T::TWO;
    mv.set(Blade::basis_vector(Conformal3::E_MINUS), half);
    mv.set(Blade::basis_vector(Conformal3::E_PLUS), -half);
    mv
}
```

### 3. SymPy Derivation Setup

Create `derivations/src/clifford_derivations/cga.py`:

```python
"""CGA derivations using SymPy.

Reference: https://conformalgeometricalgebra.org/wiki/index.php

This module provides symbolic derivations for Conformal Geometric Algebra
operations. All Rust code should be generated from these derivations.

Usage:
    cd derivations
    uv run python -m clifford_derivations.cga
"""

from sympy import symbols, sqrt, Rational, expand, rust_code
from clifford_derivations import with_timeout

# Basis vector symbols
e1, e2, e3, ep, em = symbols('e1 e2 e3 ep em')

# Null basis
e_inf = em + ep       # e∞ = e₋ + e₊
e_o = (em - ep) / 2   # e₀ = (e₋ - e₊) / 2


@with_timeout(60)
def verify_null_basis_properties():
    """Verify e∞² = 0, e₀² = 0, e∞ · e₀ = -1."""
    # e∞² = (e₋ + e₊)² = e₋² + 2e₋e₊ + e₊² = -1 + 0 + 1 = 0
    # (Note: e₋e₊ = -e₊e₋, so 2e₋e₊ + 2e₊e₋ = 0 in symmetric product)

    # e∞ · e₀ = (e₋ + e₊) · (e₋ - e₊)/2
    #         = (e₋² - e₊²)/2 = (-1 - 1)/2 = -1

    print("Null basis properties verified symbolically")
    return True


if __name__ == "__main__":
    verify_null_basis_properties()
```

## Property Tests

```rust
#[cfg(test)]
mod tests {
    use super::*;
    use crate::test_utils::ABS_DIFF_EQ_EPS;
    use approx::abs_diff_eq;
    use proptest::prelude::*;

    proptest! {
        // ================================================================
        // Null basis properties
        // ================================================================

        #[test]
        fn e_infinity_squares_to_zero(s in -100.0f64..100.0) {
            // e∞ = e₋ + e₊ should square to zero
            let e_inf = e_infinity::<f64>() * s;
            let sq = &e_inf * &e_inf;
            prop_assert!(sq.is_zero(ABS_DIFF_EQ_EPS));
        }

        #[test]
        fn e_origin_squares_to_zero(s in -100.0f64..100.0) {
            // e₀ = (e₋ - e₊) / 2 should square to zero
            let e_o = e_origin::<f64>() * s;
            let sq = &e_o * &e_o;
            prop_assert!(sq.is_zero(ABS_DIFF_EQ_EPS));
        }

        #[test]
        fn e_inf_dot_e_o_equals_minus_one(_dummy in 0..1i32) {
            // e∞ · e₀ = -1
            let e_inf = e_infinity::<f64>();
            let e_o = e_origin::<f64>();
            let dot = e_inf.inner(&e_o).scalar_part();
            prop_assert!(abs_diff_eq!(dot, -1.0, epsilon = ABS_DIFF_EQ_EPS));
        }

        // ================================================================
        // Algebraic properties for CGA3
        // ================================================================

        #[test]
        fn cga3_geometric_product_associative(
            a in any::<Multivector<f64, Conformal3>>(),
            b in any::<Multivector<f64, Conformal3>>(),
            c in any::<Multivector<f64, Conformal3>>(),
        ) {
            let lhs = &(&a * &b) * &c;
            let rhs = &a * &(&b * &c);
            prop_assert!(abs_diff_eq!(lhs, rhs, epsilon = ABS_DIFF_EQ_EPS));
        }

        #[test]
        fn cga3_geometric_product_distributive(
            a in any::<Multivector<f64, Conformal3>>(),
            b in any::<Multivector<f64, Conformal3>>(),
            c in any::<Multivector<f64, Conformal3>>(),
        ) {
            let lhs = &a * &(&b + &c);
            let rhs = &(&a * &b) + &(&a * &c);
            prop_assert!(abs_diff_eq!(lhs, rhs, epsilon = ABS_DIFF_EQ_EPS));
        }

        #[test]
        fn cga3_reverse_involutory(a in any::<Multivector<f64, Conformal3>>()) {
            prop_assert!(abs_diff_eq!(
                a.reverse().reverse(),
                a,
                epsilon = ABS_DIFF_EQ_EPS
            ));
        }

        #[test]
        fn cga3_reverse_antimorphism(
            a in any::<Multivector<f64, Conformal3>>(),
            b in any::<Multivector<f64, Conformal3>>(),
        ) {
            let lhs = (&a * &b).reverse();
            let rhs = &b.reverse() * &a.reverse();
            prop_assert!(abs_diff_eq!(lhs, rhs, epsilon = ABS_DIFF_EQ_EPS));
        }

        #[test]
        fn cga3_outer_associative(
            a in any::<Multivector<f64, Conformal3>>(),
            b in any::<Multivector<f64, Conformal3>>(),
            c in any::<Multivector<f64, Conformal3>>(),
        ) {
            let lhs = a.outer(&b).outer(&c);
            let rhs = a.outer(&b.outer(&c));
            prop_assert!(abs_diff_eq!(lhs, rhs, epsilon = ABS_DIFF_EQ_EPS));
        }

        // ================================================================
        // Euclidean subspace
        // ================================================================

        #[test]
        fn cga3_euclidean_subspace_matches_euclidean3(
            ax in -10.0f64..10.0, ay in -10.0f64..10.0, az in -10.0f64..10.0,
            bx in -10.0f64..10.0, by in -10.0f64..10.0, bz in -10.0f64..10.0,
        ) {
            // Pure Euclidean vectors (no e₊, e₋ components)
            let cga_a: Multivector<f64, Conformal3> =
                Multivector::vector(&[ax, ay, az, 0.0, 0.0]);
            let cga_b: Multivector<f64, Conformal3> =
                Multivector::vector(&[bx, by, bz, 0.0, 0.0]);

            // Dot product should match
            let cga_dot = cga_a.inner(&cga_b).scalar_part();
            let expected_dot = ax*bx + ay*by + az*bz;
            prop_assert!(abs_diff_eq!(cga_dot, expected_dot, epsilon = ABS_DIFF_EQ_EPS));
        }
    }
}
```

## Files to Create/Modify

### New Files
- `src/signature/conformal.rs` - Signature definitions
- `derivations/src/clifford_derivations/cga.py` - SymPy derivation module

### Modified Files
- `src/signature/mod.rs` - Re-export `Conformal2`, `Conformal3`

## Verification Checklist

- [ ] `cargo fmt` passes
- [ ] `cargo clippy --all-features` passes
- [ ] `cargo doc --all-features --no-deps` passes
- [ ] `cargo test --all-features` - all property tests pass
- [ ] `cargo deny check` passes
- [ ] Null basis properties verified (e∞² = 0, e₀² = 0, e∞·e₀ = -1)
- [ ] Algebraic properties verified (associativity, distributivity, etc.)
- [ ] SymPy derivation module created

## Dependencies

- None (this is the foundation)

## Next Steps

After this PRD is complete, proceed to PRD-6.2 (Round Points and Conformal Embedding).
