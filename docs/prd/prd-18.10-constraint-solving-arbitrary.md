# PRD-18.10: Constraint-Solving Arbitrary Generation

**Status**: Draft
**Parent**: [PRD-18: Constraint System Redesign](prd-18-constraint-redesign.md)

## Goal

Generate `Arbitrary` implementations that solve geometric constraints symbolically, using a minimal set of free variables.

## Problem

Current Arbitrary implementations either:
- Use random coefficients that violate constraints (Study condition, Plücker condition)
- Use factory methods which don't explore the full valid parameter space

For example, the generated `Arbitrary for Motor<T>` creates 8 random coefficients, but valid motors must satisfy the Study condition:
```
s·e₀₁₂₃ + e₂₃·e₀₁ + e₃₁·e₀₂ + e₁₂·e₀₃ = 0
```

Random coefficients almost never satisfy this constraint.

## Solution

Use Symbolica to:
1. Identify the constraint equations for each constrained type
2. Solve for dependent variables in terms of independent variables
3. Generate Arbitrary that picks random values for independent variables, computes dependent variables

### Example: Motor with Study Condition

The Study condition is linear in `e0123`, so we can solve:
```
e0123 = -(e23*e01 + e31*e02 + e12*e03) / s
```

Generated Arbitrary:
```rust
impl<T: Float + Debug + 'static> Arbitrary for Motor<T> {
    fn arbitrary_with(_: Self::Parameters) -> Self::Strategy {
        (
            (-10.0f64..10.0).prop_filter("non-zero s", |s| s.abs() > 0.1),
            -10.0f64..10.0,  // e23
            -10.0f64..10.0,  // e31
            -10.0f64..10.0,  // e12
            -10.0f64..10.0,  // e01
            -10.0f64..10.0,  // e02
            -10.0f64..10.0,  // e03
        )
            .prop_map(|(s, e23, e31, e12, e01, e02, e03)| {
                // Symbolically derived solution
                let e0123 = -(e23 * e01 + e31 * e02 + e12 * e03) / s;
                Motor::new_unchecked(
                    T::from_f64(s), T::from_f64(e23), T::from_f64(e31), T::from_f64(e12),
                    T::from_f64(e01), T::from_f64(e02), T::from_f64(e03), T::from_f64(e0123),
                )
            })
            .boxed()
    }
}
```

### Constraint Types

| Constraint | Solvability | Approach |
|------------|-------------|----------|
| Linear in one variable | Solvable | Solve for that variable |
| Quadratic | Filter-based | Generate random, filter valid |
| Multiple constraints | Chain solutions | Solve sequentially |

## Implementation

### New Files

- `crates/clifford-codegen/src/symbolic/constraint_solve.rs` - Solve constraints for dependent variables

### Modified Files

- `crates/clifford-codegen/src/codegen/traits.rs` - Update Arbitrary generation

### Tasks

1. **Add `solve_constraint_for_variable()` function**
   - Input: Symbolica `Atom` expression that must equal zero, target variable
   - Output: Expression for target variable in terms of other variables
   - Use Symbolica's equation solving capabilities

2. **Identify which variable to solve for**
   - Prefer highest-grade variable (e.g., pseudoscalar `e0123`)
   - Ensure the constraint is linear in that variable
   - Fall back to quadratic filter if not linear

3. **Generate prop_filter for denominator conditions**
   - Extract denominator from solution
   - Add filter to ensure denominator is non-zero (with epsilon margin)

4. **Generate prop_map with constraint solution**
   - Generate code for N-1 random variables
   - Compute dependent variable using derived formula

5. **Handle quadratic constraints**
   - For Plücker condition (quadratic), use prop_filter instead of solving
   - Generate: `prop_filter("satisfies_plucker", |l| l.satisfies_plucker_condition(eps))`

## Verification

After implementation:
```rust
proptest! {
    #[test]
    fn arbitrary_motor_satisfies_study_condition(m in any::<Motor<f64>>()) {
        let residual = m.s() * m.e0123() + m.e23() * m.e01() + m.e31() * m.e02() + m.e12() * m.e03();
        prop_assert!(relative_eq!(residual, 0.0, epsilon = RELATIVE_EQ_EPS, max_relative = RELATIVE_EQ_EPS));
    }
}
```

## Success Criteria

1. All constrained types have Arbitrary impls that satisfy constraints
2. No factory-method workarounds needed
3. Full parameter space explored (not just rotation+translation compositions)
4. Proptest verifies constraint satisfaction
