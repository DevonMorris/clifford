# PRD: CGA 2D Basis Convention Alignment

**Status:** Draft
**Author:** Claude
**Created:** 2026-01-18

## Problem Statement

The CGA 2D implementation has an inconsistency between the declared metric signature and the conceptual model used in constructors and formulas.

### Current State

1. **TOML Signature** declares orthonormal basis:
   ```toml
   [signature]
   positive = ["e1", "e2", "e3"]  # e3² = +1
   negative = ["e4"]              # e4² = -1
   ```

2. **Field Names** use null basis semantics:
   - `o` (origin weight) for e3 coefficient
   - `i` (infinity weight) for e4 coefficient

3. **Constructors** embed points assuming null basis:
   ```rust
   // Currently: P = x·e₁ + y·e₂ + 1·e₃ + ½(x²+y²)·e₄
   fn from_euclidean(x, y) -> Self {
       Self::new_unchecked(x, y, T::one(), (x*x + y*y) / T::TWO)
   }
   ```

4. **Distance formula** manually computes a bilinear form that matches null basis inner product, bypassing the generated metric:
   ```rust
   // Computes: x₁x₂ + y₁y₂ - o₁i₂ - i₁o₂ (null basis inner product)
   // NOT: x₁x₂ + y₁y₂ + o₁o₂ - i₁i₂ (Cl(3,1) orthonormal inner product)
   ```

### The Conflict

In standard CGA literature:
- **Orthonormal basis**: {e₁, e₂, e₊, e₋} where e₊² = +1, e₋² = -1
- **Null basis**: e₀ = (e₋ - e₊)/2, e∞ = e₋ + e₊

The null basis vectors satisfy:
- e₀² = 0 (null)
- e∞² = 0 (null)
- e₀ · e∞ = -1

The current code conflates these two representations:
- Signature says e3² = +1, e4² = -1 (orthonormal)
- Formulas treat e3, e4 as if they were e₀, e∞ (null)

This works only because extension methods bypass the generated algebra operations. But it's conceptually incorrect and will break if users try to:
1. Use generated inner products directly
2. Verify null vector constraints (P·P should equal 0 for points)
3. Apply standard CGA formulas from textbooks

## Requirements

### 1. Correct Orthonormal Embedding

Constructors must embed Euclidean points using the **orthonormal basis** {e₁, e₂, e₊, e₋}:

For a point P at Euclidean coordinates (x, y):
```
P = x·e₁ + y·e₂ + ep_coeff·e₊ + em_coeff·e₋
```

Where the e₊ and e₋ coefficients are computed from the null basis embedding:
- Null basis: P = x·e₁ + y·e₂ + o·e₀ + i·e∞ with o=1, i=½(x²+y²)
- Substituting e₀ = (e₋ - e₊)/2, e∞ = e₋ + e₊:
  - e₊ coefficient = -o/2 + i = -½ + ½(x²+y²) = ½(x²+y² - 1)
  - e₋ coefficient = o/2 + i = ½ + ½(x²+y²) = ½(x²+y² + 1)

### 2. Rename Basis Vectors in TOML

Change the signature to use descriptive names:
```toml
[signature]
positive = ["e1", "e2", "ep"]  # ep = e₊, squares to +1
negative = ["em"]              # em = e₋, squares to -1
```

### 3. Update Field Naming with Documentation

Field names should indicate their role while documentation explains the basis relationship:

For RoundPoint:
```toml
field_map = [
  { name = "x", blade = "e1" },    # Euclidean x-coordinate
  { name = "y", blade = "e2" },    # Euclidean y-coordinate
  { name = "ep", blade = "ep" },   # e₊ coefficient (positive conformal basis)
  { name = "em", blade = "em" }    # e₋ coefficient (negative conformal basis)
]
```

### 4. Document Non-Semantic Fields

Add a `semantic = false` or similar flag to indicate fields that are NOT semantically named:

```toml
field_map = [
  { name = "x", blade = "e1" },
  { name = "y", blade = "e2" },
  { name = "ep", blade = "ep", semantic = false, doc = "Conformal basis coefficient, not a coordinate" },
  { name = "em", blade = "em", semantic = false, doc = "Conformal basis coefficient, not a coordinate" }
]
```

### 5. Update Constructors

**RoundPoint::from_euclidean(x, y)**:
```rust
pub fn from_euclidean(x: T, y: T) -> Self {
    let sq = x * x + y * y;
    let half = T::one() / T::TWO;
    // Null basis: o=1, i=½(x²+y²)
    // e₊ = -o/2 + i = (sq - 1)/2
    // e₋ = o/2 + i = (sq + 1)/2
    let ep_coeff = (sq - T::one()) * half;
    let em_coeff = (sq + T::one()) * half;
    Self::new_unchecked(x, y, ep_coeff, em_coeff)
}
```

**RoundPoint::origin()**:
```rust
pub fn origin() -> Self {
    // (0, 0): sq=0, ep=-½, em=½
    let half = T::one() / T::TWO;
    Self::new_unchecked(T::zero(), T::zero(), -half, half)
}
```

**RoundPoint::infinity()**:
```rust
pub fn infinity() -> Self {
    // e∞ = e₋ + e₊, normalized to unit e∞
    Self::new_unchecked(T::zero(), T::zero(), T::one(), T::one())
}
```

### 6. Update Coordinate Extraction

**to_euclidean()**: Extract (x, y) from orthonormal representation:
```rust
pub fn to_euclidean(&self) -> Option<(T, T)> {
    // Recover o = em - ep, i = (em + ep)/2 from orthonormal basis
    // Then x_eucl = x/o, y_eucl = y/o
    let o = self.em() - self.ep();  // o = origin weight
    if o.abs() < T::epsilon() {
        None
    } else {
        Some((self.x() / o, self.y() / o))
    }
}
```

### 7. Null Vector Property Tests

Add tests verifying the null vector constraint P·P = 0:

```rust
#[test]
fn test_point_is_null_vector() {
    use crate::ops::Dot;  // or Inner

    let p = RoundPoint::from_euclidean(3.0_f64, 4.0);
    let p_dot_p = p.dot(&p);

    // For a properly embedded point, P·P should be 0
    assert!(
        relative_eq!(p_dot_p, 0.0, epsilon = 1e-10),
        "Point should be a null vector: P·P = {} ≠ 0", p_dot_p
    );
}
```

### 8. Circle Center/Radius Tests

Add comprehensive tests that circles from three points have expected properties:

```rust
#[test]
fn test_circle_from_three_points_center_radius() {
    // Test case: circle centered at (2, 3) with radius 5
    let cx = 2.0_f64;
    let cy = 3.0;
    let r = 5.0;

    // Three points on the circle (at 0°, 120°, 240°)
    let p1 = RoundPoint::from_euclidean(cx + r, cy);
    let p2 = RoundPoint::from_euclidean(
        cx + r * (2.0 * PI / 3.0).cos(),
        cy + r * (2.0 * PI / 3.0).sin()
    );
    let p3 = RoundPoint::from_euclidean(
        cx + r * (4.0 * PI / 3.0).cos(),
        cy + r * (4.0 * PI / 3.0).sin()
    );

    let circle = Circle::from_three_points(&p1, &p2, &p3);

    let (extracted_cx, extracted_cy) = circle.center().unwrap();
    let extracted_r = circle.radius().unwrap();

    assert!(relative_eq!(extracted_cx, cx, epsilon = 1e-10));
    assert!(relative_eq!(extracted_cy, cy, epsilon = 1e-10));
    assert!(relative_eq!(extracted_r, r, epsilon = 1e-10));

    // Verify the three points lie on the circle (distance from center = radius)
    for p in [&p1, &p2, &p3] {
        let (px, py) = p.to_euclidean().unwrap();
        let dist = ((px - cx).powi(2) + (py - cy).powi(2)).sqrt();
        assert!(relative_eq!(dist, r, epsilon = 1e-10));
    }
}
```

### 9. Update Circle Center/Radius Extraction

The formulas for center() and radius() may need to change based on the new basis. Derive the correct formulas from:

```
Circle = P₁ ∧ P₂ ∧ P₃
```

where P₁, P₂, P₃ are in orthonormal representation.

### 10. TOML Documentation Updates

Add a documentation section explaining the convention:

```toml
# ============================================================
# Basis Convention Documentation
# ============================================================
#
# This algebra uses the ORTHONORMAL conformal basis:
#   e1, e2: Euclidean basis vectors (e1² = e2² = +1)
#   ep:     Positive conformal basis (ep² = +1)
#   em:     Negative conformal basis (em² = -1)
#
# The commonly-used NULL basis (e₀, e∞) relates to orthonormal as:
#   e₀ (origin)   = (em - ep) / 2
#   e∞ (infinity) = em + ep
#
# Point embedding: For Euclidean point (x, y)
#   Null basis:       P = x·e₁ + y·e₂ + 1·e₀ + ½(x²+y²)·e∞
#   Orthonormal:      P = x·e₁ + y·e₂ + ½(x²+y²-1)·eₚ + ½(x²+y²+1)·eₘ
#
# Field naming convention:
#   Fields x, y are SEMANTIC (directly represent Euclidean coordinates)
#   Fields ep, em are NON-SEMANTIC (conformal embedding coefficients)
#
# ============================================================
```

## Acceptance Criteria

1. **Signature correctness**: `e1² = e2² = ep² = +1`, `em² = -1`
2. **Null vector property**: `RoundPoint::from_euclidean(x, y).inner(self) ≈ 0`
3. **Coordinate roundtrip**: `from_euclidean(x, y).to_euclidean() == Some((x, y))`
4. **Circle properties**: `Circle::from_three_points` produces correct center and radius
5. **Distance formula**: `p1.distance(&p2)` matches Euclidean distance
6. **Documentation**: TOML explains orthonormal vs null basis relationship
7. **Non-semantic fields**: TOML marks ep/em fields as non-semantic

## Implementation Steps

1. **Update TOML** (algebras/conformal2.toml)
   - Rename e3/e4 to ep/em in signature
   - Update all field_map entries
   - Add documentation section explaining conventions
   - Mark non-semantic fields

2. **Regenerate algebra**
   ```bash
   cargo run -p clifford-codegen -- generate algebras/conformal2.toml
   ```

3. **Update extensions.rs**
   - Fix `from_euclidean()` to use correct orthonormal embedding
   - Fix `to_euclidean()` to extract from orthonormal form
   - Fix `origin()` and `infinity()`
   - Update `distance_squared()` to use generated inner product
   - Update `Circle::center()` and `Circle::radius()` formulas

4. **Add tests**
   - Null vector property tests
   - Circle center/radius tests with various configurations
   - Motor transformation tests
   - Distance formula tests

5. **Update module documentation**
   - Explain the basis convention choice
   - Link to standard CGA references

## Out of Scope

- Changing CGA 3D (conformal3.toml) - separate PRD if needed
- Adding a `NullBasisPoint` type that stores (o, i) directly
- Dual number support for automatic differentiation

## Open Questions

1. Should we provide helper methods `null_origin_weight()` and `null_infinity_weight()` that compute o and i from the orthonormal coefficients?

2. Should the generated inner product trait be used directly in `distance_squared()`, or should we continue with the manual formula for clarity?

3. Do we need a `semantic` field in the TOML schema, or is documentation sufficient?
