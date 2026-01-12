# PRD-17.6: Geometric Antiproduct for PGA Transformations

**Status**: Draft
**Parent**: [PRD-17](prd-17-codegen-products.md)
**Goal**: Add geometric antiproduct support to fix PGA motor/flector transformations

## Problem Statement

The generated sandwich products for PGA (Projective Geometric Algebra) don't correctly transform points, lines, and planes. For example, a pure translation motor doesn't translate points:

```rust
let motor = Motor::from_translation(10.0, 20.0, 30.0);
let p = Point::from_cartesian(1.0, 2.0, 3.0);
let result = motor.transform_point(&p);
// Expected: (11.0, 22.0, 33.0)
// Actual: (1.0, 2.0, 3.0) - no translation!
```

### Root Cause

In PGA with degenerate metric (e0² = 0), the standard geometric product sandwich `M * P * rev(M)` doesn't correctly handle translations because:

1. Translation components (e01, e02, e03) multiply with e0 (point's homogeneous coordinate)
2. Since e0² = 0, terms like `e01 * e0 = e0*e1*e0 = -e0²*e1 = 0`
3. Translation contributions vanish from spatial coordinates

### Solution: Geometric Antiproduct

According to [Rigid Geometric Algebra](https://rigidgeometricalgebra.org/wiki/index.php?title=Motor), PGA transformations should use the **geometric antiproduct** sandwich:

```
x' = Q ⊛ x ⊛ rev(Q)   // NOT: Q * x * rev(Q)
```

Where `⊛` is the geometric antiproduct (also written as ⟑ or antidot).

## Background: Antiproduct Theory

### Duality Relationship

The geometric antiproduct is **dual** to the geometric product. For elements a, b in Cl(p,q,r):

```
complement(a ⊙ b) = complement(a) ⊛ complement(b)
complement(a ⊛ b) = complement(a) ⊙ complement(b)
```

Where:
- `⊙` = geometric product
- `⊛` = geometric antiproduct
- `complement(x)` = dual of x (multiply by pseudoscalar)

### Metric for Antiproduct

The antiproduct uses "anti-basis vectors" with inverted metric:
- Basis vectors that square to +1 in geometric product → anti-vectors square to +1
- Basis vectors that square to -1 in geometric product → anti-vectors square to -1
- **Degenerate basis (e0² = 0) → anti-vector squares to non-zero!**

This is why antiproduct works for PGA: the degenerate direction becomes non-degenerate.

### Antiwedge Product

The **antiwedge** (exterior antiproduct) combines "empty dimensions":

```
antigrade(a ∨ b) = antigrade(a) + antigrade(b)
```

Where antigrade = n - grade for n-dimensional algebra.

De Morgan laws:
```
complement(a) ∨ complement(b) = complement(a ∧ b)
complement(a) ∧ complement(b) = complement(a ∨ b)
```

## Implementation Plan

### Phase 1: Multivector Antiproduct

Add antiproduct method to generic Multivector for testing:

```rust
impl<T: Float, S: Signature> Multivector<T, S> {
    /// Geometric antiproduct: a ⊛ b = complement(complement(a) * complement(b))
    ///
    /// The antiproduct is dual to the geometric product. In PGA, it's used
    /// for motor transformations because it correctly handles the degenerate
    /// direction.
    pub fn antiproduct(&self, other: &Self) -> Self {
        // Method 1: Via complements (dual)
        // a ⊛ b = dual(dual(a) * dual(b))
        self.dual().geometric(&other.dual()).dual()
    }

    /// Antiwedge product: a ∨ b = complement(complement(a) ∧ complement(b))
    pub fn antiwedge(&self, other: &Self) -> Self {
        self.dual().exterior(&other.dual()).dual()
    }
}
```

### Phase 2: Codegen Antiproduct Products

Add antiproduct generation to clifford-codegen:

```rust
pub enum ProductKind {
    Geometric,
    Exterior,
    LeftContraction,
    Regressive,
    Scalar,
    // New
    Antiproduct,   // Geometric antiproduct (⊛)
    Antiwedge,     // Exterior antiproduct (∨)
}
```

#### Grade Rules for Antiproduct

The antiproduct output grades follow from the duality:

```rust
fn antiproduct_grades(grade_a: usize, grade_b: usize, dim: usize) -> Vec<usize> {
    // antigrade(a ⊛ b) computed via geometric product of complements
    let antigrade_a = dim - grade_a;
    let antigrade_b = dim - grade_b;

    // Geometric product of complements gives grades from |ag_a - ag_b| to ag_a + ag_b
    geometric_grades(antigrade_a, antigrade_b)
        .iter()
        .map(|&g| dim - g)  // Convert antigrade back to grade
        .collect()
}
```

### Phase 3: Antiproduct Sandwich

Add antiproduct sandwich generation for versors:

```rust
/// Generates antiproduct sandwich: v ⊛ x ⊛ rev(v)
fn generate_antiproduct_sandwich(
    &self,
    versor_type: &TypeSpec,
    operand_type: &TypeSpec,
) -> TokenStream {
    // Similar to geometric sandwich but using antiproduct
}
```

### Phase 4: TOML Configuration

Update algebra TOML to specify which sandwich type to use:

```toml
[types.Motor]
grades = [0, 2, 4]
versor = true
# Specify antiproduct sandwich for PGA
sandwich = { product = "antiproduct", targets = ["Point", "Line", "Plane", "Motor", "Flector"] }

[types.Flector]
grades = [1, 3]
versor = true
sandwich = { product = "antiproduct", targets = ["Point", "Line", "Plane", "Motor", "Flector"] }
```

### Phase 5: Verification Tests

Add tests comparing antiproduct sandwich against known-correct results:

```rust
proptest! {
    #[test]
    fn motor_translation_works(dx in -10.0..10.0, dy in -10.0..10.0, dz in -10.0..10.0,
                               px in -10.0..10.0, py in -10.0..10.0, pz in -10.0..10.0) {
        let motor = Motor::from_translation(dx, dy, dz);
        let p = Point::from_cartesian(px, py, pz);
        let result = motor.transform_point(&p);

        prop_assert!(abs_diff_eq!(result.x(), px + dx, epsilon = ABS_DIFF_EQ_EPS));
        prop_assert!(abs_diff_eq!(result.y(), py + dy, epsilon = ABS_DIFF_EQ_EPS));
        prop_assert!(abs_diff_eq!(result.z(), pz + dz, epsilon = ABS_DIFF_EQ_EPS));
    }

    #[test]
    fn antiproduct_sandwich_matches_multivector(
        m in any::<Motor<f64>>(),
        p in any::<Point<f64>>()
    ) {
        let specialized = m.transform_point(&p);

        let mv_m: Multivector<f64, Projective3> = m.into();
        let mv_p: Multivector<f64, Projective3> = p.into();
        let mv_result = mv_m.antiproduct(&mv_p).antiproduct(&mv_m.reverse());

        let expected: Point<f64> = mv_result.into();
        prop_assert!(abs_diff_eq!(specialized, expected, epsilon = ABS_DIFF_EQ_EPS));
    }
}
```

## Deliverables

### Multivector Methods
- [ ] Add `Multivector::antiproduct()` method
- [ ] Add `Multivector::antiwedge()` method
- [ ] Add `Multivector::dual()` method (if not present)
- [ ] Add tests for antiproduct algebraic properties

### Codegen Updates
- [ ] Add `ProductKind::Antiproduct` variant
- [ ] Add `ProductKind::Antiwedge` variant
- [ ] Implement antiproduct grade computation
- [ ] Implement antiproduct term generation
- [ ] Add antiproduct sandwich generation
- [ ] Update TOML parser to support `sandwich.product` field

### PGA Fixes
- [ ] Update projective3.toml to use antiproduct sandwich
- [ ] Regenerate projective3 algebra
- [ ] Verify all nalgebra comparison tests pass
- [ ] Verify motor/flector transformation tests pass

### Documentation
- [ ] Document antiproduct in codegen
- [ ] Add antiproduct examples to Multivector docs
- [ ] Update CLAUDE.md with PGA transformation guidance

## Files Changed

| File | Action |
|------|--------|
| `src/algebra/multivector.rs` | Add antiproduct, antiwedge methods |
| `crates/clifford-codegen/src/codegen/products.rs` | Add antiproduct generation |
| `crates/clifford-codegen/src/spec/ir.rs` | Add sandwich product type |
| `crates/clifford-codegen/src/spec/parser.rs` | Parse sandwich.product field |
| `algebras/projective3.toml` | Update to antiproduct sandwich |
| `src/specialized/projective/dim3/generated/*` | Regenerate with antiproduct |

## Testing Strategy

1. **Unit tests**: Verify antiproduct algebraic identities
2. **Property tests**: antiproduct sandwich matches Multivector computation
3. **Nalgebra comparison**: Motor transformations match nalgebra Isometry3
4. **Concrete cases**: Known transformations (90° rotations, axis translations)

## Success Criteria

1. `Motor::transform_point()` correctly translates points
2. `Motor::transform_point()` correctly rotates points
3. All nalgebra comparison tests pass
4. Antiproduct has verification tests against Multivector
5. Generated code uses antiproduct for PGA versors

## References

- [Rigid Geometric Algebra - Motor](https://rigidgeometricalgebra.org/wiki/index.php?title=Motor)
- [Rigid Geometric Algebra - Geometric Antiproduct](https://rigidgeometricalgebra.org/wiki/index.php?title=Geometric_antiproduct)
- [Rigid Geometric Algebra - Antiwedge Product](https://rigidgeometricalgebra.org/wiki/index.php?title=Antiwedge_product)
- [Look, Ma, No Matrices!](https://enkimute.github.io/LookMaNoMatrices/) - PGA tutorial
