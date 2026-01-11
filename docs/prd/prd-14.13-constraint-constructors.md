# PRD-14.13: Constraint-Enforcing Constructors

**Status**: Draft
**Parent**: [PRD-14: Geometric Algebra Code Generator](prd-14-codegen.md)
**Goal**: Generate constructors that enforce geometric constraints by solving for dependent coefficients

## Overview

Types with geometric constraints (e.g., Motor, Bivector, Flector in PGA) currently have all coefficients as independent parameters. However, these coefficients are related by constraint equations. This PRD implements:

1. **`new`** - Constructor over a minimal subset of coefficients; computes dependent coefficients from the constraint
2. **`new_checked`** - Constructor over all coefficients; validates constraints within tolerance; returns `Result`
3. **`new_unchecked`** - Constructor over all coefficients; no validation (for performance-critical code and AD)

## Problem Statement

Consider the PGA Motor with constraint:
```
2*s*e0123 - 2*e12*e03 + 2*e13*e02 - 2*e23*e01 = 0
```

Currently, users must manually ensure this constraint holds. This is error-prone and the constraint relationship isn't obvious from the API.

### Goals

1. Provide a `new()` constructor that takes 7 independent coefficients and computes the 8th
2. Allow users to choose which coefficient is computed (configurable in TOML)
3. Provide `new_checked()` for validation when all coefficients are known
4. Maintain `new_unchecked()` for raw construction (AD, deserialization, internal ops)

## Constraint Solving

### Linear Constraints

All current geometric constraints are linear in the coefficients. For example:

**Motor**: `2*s*e0123 - 2*e12*e03 + 2*e13*e02 - 2*e23*e01 = 0`

This can be solved for any variable:
- Solve for `e0123`: `e0123 = (e12*e03 - e13*e02 + e23*e01) / s`
- Solve for `e01`: `e01 = (s*e0123 - e12*e03 + e13*e02) / e23`
- etc.

### Solving Strategy

1. Parse the constraint expression using the symbolic parser
2. Identify the "solved-for" variable from the TOML config
3. Use Symbolica to solve the constraint algebraically
4. Generate the computation expression for the dependent variable

### Degenerate Cases (Division by Zero)

Solving for a variable may introduce division. When the divisor is zero, the constraint equation becomes **degenerate** - the solved-for variable drops out entirely, meaning any value satisfies the constraint.

For example, with Motor constraint `2*s*e0123 - 2*e12*e03 + 2*e13*e02 - 2*e23*e01 = 0`:
- When `s ≠ 0`: Solve for `e0123 = (e12*e03 - e13*e02 + e23*e01) / s`
- When `s = 0`: Constraint becomes `-2*e12*e03 + 2*e13*e02 - 2*e23*e01 = 0`, which doesn't involve `e0123` at all

In the degenerate case, `e0123` can be **any value** and still satisfy the constraint. We choose a **canonical value** (typically zero) for determinism:

```rust
impl<T: Float> Motor<T> {
    /// Creates a Motor from independent coefficients.
    ///
    /// The `e0123` coefficient is computed from the constraint:
    /// `2*s*e0123 - 2*e12*e03 + 2*e13*e02 - 2*e23*e01 = 0`
    ///
    /// When `s` is zero (degenerate case), `e0123` is set to zero.
    pub fn new(s: T, e12: T, e13: T, e23: T, e01: T, e02: T, e03: T) -> Self {
        let e0123 = if s.abs() > T::epsilon() {
            (e12 * e03 - e13 * e02 + e23 * e01) / s
        } else {
            T::zero()  // Canonical choice for degenerate case
        };
        Self { s, e12, e13, e23, e01, e02, e03, e0123 }
    }
}
```

This makes `new()` total (never panics), which is important for ergonomics and AD compatibility.

For constraints without division (e.g., `a + b + c = 0` solved for `c = -a - b`), this issue doesn't arise.

## TOML Configuration

### Syntax

Add `solve_for` field to types with constraints:

```toml
[types.Motor]
grades = [0, 2, 4]
fields = ["s", "e12", "e13", "e23", "e01", "e02", "e03", "e0123"]
geometric_constraint = "2*s*e0123 - 2*e12*e03 + 2*e13*e02 - 2*e23*e01 = 0"
antiproduct_constraint = "2*s*e0123 - 2*e12*e03 + 2*e13*e02 - 2*e23*e01 = 0"
solve_for = "e0123"  # <-- New field
```

### Constraint Selection

When both `geometric_constraint` and `antiproduct_constraint` are present and equal (as they often are), only one constraint needs solving. When they differ, both must be satisfied, which may require solving for multiple variables or imposing multiple constraints.

For the initial implementation, we assume a single constraint equation is sufficient (the two constraints are often equivalent expressions).

### Default Behavior

If `solve_for` is not specified:
- If the type has a constraint, emit a warning and skip `new()` generation
- Always generate `new_checked()` and `new_unchecked()`

## Generated Code

### Example: Motor

```rust
/// 3D PGA Motor (rigid body transformation).
///
/// Motors have 8 coefficients but only 7 degrees of freedom due to
/// the geometric constraint: `2*s*e0123 - 2*e12*e03 + 2*e13*e02 - 2*e23*e01 = 0`
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Motor<T: Float> {
    s: T,
    e12: T,
    e13: T,
    e23: T,
    e01: T,
    e02: T,
    e03: T,
    e0123: T,
}

impl<T: Float> Motor<T> {
    /// Creates a Motor from 7 independent coefficients.
    ///
    /// The `e0123` coefficient is computed from the geometric constraint.
    /// When `s` is zero (degenerate case), `e0123` defaults to zero.
    #[inline]
    pub fn new(s: T, e12: T, e13: T, e23: T, e01: T, e02: T, e03: T) -> Self {
        let e0123 = if s.abs() > T::epsilon() {
            (e12 * e03 - e13 * e02 + e23 * e01) / s
        } else {
            T::zero()
        };
        Self { s, e12, e13, e23, e01, e02, e03, e0123 }
    }

    /// Creates a Motor from all 8 coefficients with constraint validation.
    ///
    /// Returns an error if the geometric constraint is not satisfied within
    /// the given tolerance.
    ///
    /// # Errors
    ///
    /// Returns `ConstraintError` if `|2*s*e0123 - 2*e12*e03 + 2*e13*e02 - 2*e23*e01| > tolerance`.
    #[inline]
    pub fn new_checked(
        s: T,
        e12: T,
        e13: T,
        e23: T,
        e01: T,
        e02: T,
        e03: T,
        e0123: T,
        tolerance: T,
    ) -> Result<Self, ConstraintError> {
        let residual = T::TWO * s * e0123
            - T::TWO * e12 * e03
            + T::TWO * e13 * e02
            - T::TWO * e23 * e01;

        if residual.abs() > tolerance {
            return Err(ConstraintError::new(
                "Motor",
                "2*s*e0123 - 2*e12*e03 + 2*e13*e02 - 2*e23*e01 = 0",
                residual,
            ));
        }

        Ok(Self { s, e12, e13, e23, e01, e02, e03, e0123 })
    }

    /// Creates a Motor from all 8 coefficients without validation.
    ///
    /// # Safety (Logical)
    ///
    /// Caller must ensure the geometric constraint is satisfied.
    /// Use this for performance-critical code, automatic differentiation,
    /// or when coefficients come from trusted sources (e.g., product operations).
    #[inline]
    pub fn new_unchecked(
        s: T,
        e12: T,
        e13: T,
        e23: T,
        e01: T,
        e02: T,
        e03: T,
        e0123: T,
    ) -> Self {
        Self { s, e12, e13, e23, e01, e02, e03, e0123 }
    }
}
```

### Example: Bivector (simpler constraint)

```rust
/// PGA Bivector (line representation).
///
/// Constraint: `-2*e12*e03 + 2*e13*e02 - 2*e23*e01 = 0`
impl<T: Float> Bivector<T> {
    /// Creates a Bivector from 5 independent coefficients.
    ///
    /// The `e03` coefficient is computed from the constraint.
    /// When `e12` is zero (degenerate case), `e03` defaults to zero.
    pub fn new(e12: T, e13: T, e23: T, e01: T, e02: T) -> Self {
        let e03 = if e12.abs() > T::epsilon() {
            (e13 * e02 - e23 * e01) / e12
        } else {
            T::zero()
        };
        Self { e12, e13, e23, e01, e02, e03 }
    }

    /// Creates a Bivector from all 6 coefficients with validation.
    pub fn new_checked(
        e12: T,
        e13: T,
        e23: T,
        e01: T,
        e02: T,
        e03: T,
        tolerance: T,
    ) -> Result<Self, ConstraintError> {
        let residual = -T::TWO * e12 * e03
            + T::TWO * e13 * e02
            - T::TWO * e23 * e01;

        if residual.abs() > tolerance {
            return Err(ConstraintError::new("Bivector", "...", residual));
        }

        Ok(Self { e12, e13, e23, e01, e02, e03 })
    }

    /// Creates a Bivector without validation.
    pub fn new_unchecked(e12: T, e13: T, e23: T, e01: T, e02: T, e03: T) -> Self {
        Self { e12, e13, e23, e01, e02, e03 }
    }
}
```

## Error Type

Generate a constraint error type in the module:

```rust
/// Error returned when a geometric constraint is not satisfied.
#[derive(Debug, Clone)]
pub struct ConstraintError {
    /// The type that failed the constraint check.
    pub type_name: &'static str,
    /// The constraint expression that was violated.
    pub constraint: &'static str,
    /// The residual value (how far from zero).
    pub residual: f64,
}

impl ConstraintError {
    /// Creates a new constraint error.
    pub fn new<T: Into<f64>>(
        type_name: &'static str,
        constraint: &'static str,
        residual: T,
    ) -> Self {
        Self {
            type_name,
            constraint,
            residual: residual.into(),
        }
    }
}

impl std::fmt::Display for ConstraintError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "{} constraint violated: {} (residual: {:.2e})",
            self.type_name, self.constraint, self.residual
        )
    }
}

impl std::error::Error for ConstraintError {}
```

## Implementation Plan

### Phase 1: Symbolic Constraint Solving

1. Extend `ConstraintParser` to extract the solved-for variable
2. Use Symbolica's solver to isolate the variable
3. Generate the computation expression as Rust code

### Phase 2: Code Generation

1. Add `solve_for` field to `RawTypeSpec` and `TypeSpec`
2. Update `TypeGenerator` to:
   - Generate `new()` with computed coefficient
   - Generate `new_checked()` with validation
   - Rename existing `new()` to `new_unchecked()`
3. Generate `ConstraintError` type

### Phase 3: Integration

1. Update existing TOML specs with `solve_for` fields
2. Regenerate all algebras
3. Update tests to use new constructors

## TOML Specification Updates

### projective3.toml

```toml
[types.Bivector]
grades = [2]
fields = ["e12", "e13", "e23", "e01", "e02", "e03"]
geometric_constraint = "-2*e12*e03 + 2*e13*e02 - 2*e23*e01 = 0"
solve_for = "e03"  # Solve for moment component

[types.Motor]
grades = [0, 2, 4]
fields = ["s", "e12", "e13", "e23", "e01", "e02", "e03", "e0123"]
geometric_constraint = "2*s*e0123 - 2*e12*e03 + 2*e13*e02 - 2*e23*e01 = 0"
solve_for = "e0123"  # Solve for pseudoscalar component

[types.Flector]
grades = [1, 3]
fields = ["e1", "e2", "e3", "e0", "e123", "e012", "e013", "e023"]
geometric_constraint = "-2*e1*e023 + 2*e2*e013 - 2*e3*e012 + 2*e0*e123 = 0"
solve_for = "e023"  # Solve for bulk trivector component
```

## Testing

### Unit Tests

```rust
#[cfg(test)]
mod tests {
    use super::*;
    use approx::abs_diff_eq;

    #[test]
    fn motor_new_satisfies_constraint() {
        let m = Motor::new(1.0, 0.1, 0.2, 0.3, 0.01, 0.02, 0.03);

        // Verify constraint holds
        let residual = 2.0 * m.s() * m.e0123()
            - 2.0 * m.e12() * m.e03()
            + 2.0 * m.e13() * m.e02()
            - 2.0 * m.e23() * m.e01();

        assert!(abs_diff_eq!(residual, 0.0, epsilon = 1e-10));
    }

    #[test]
    fn motor_new_checked_accepts_valid() {
        // Construct valid motor using new(), then verify new_checked accepts it
        let m1 = Motor::new(1.0, 0.1, 0.2, 0.3, 0.01, 0.02, 0.03);

        let m2 = Motor::new_checked(
            m1.s(), m1.e12(), m1.e13(), m1.e23(),
            m1.e01(), m1.e02(), m1.e03(), m1.e0123(),
            1e-10,
        );

        assert!(m2.is_ok());
    }

    #[test]
    fn motor_new_checked_rejects_invalid() {
        // Construct motor with violated constraint
        let result = Motor::new_checked(
            1.0, 0.1, 0.2, 0.3, 0.01, 0.02, 0.03,
            999.0,  // Wrong value for e0123
            1e-10,
        );

        assert!(result.is_err());
    }

    #[test]
    fn motor_new_degenerate_case_uses_canonical_value() {
        // s = 0 is degenerate; e0123 should default to zero
        let m = Motor::new(0.0, 0.1, 0.2, 0.3, 0.01, 0.02, 0.03);

        assert!(abs_diff_eq!(m.e0123(), 0.0, epsilon = 1e-10));

        // Verify constraint is still satisfied (0 = 0 when s = 0)
        let residual = 2.0 * m.s() * m.e0123()
            - 2.0 * m.e12() * m.e03()
            + 2.0 * m.e13() * m.e02()
            - 2.0 * m.e23() * m.e01();

        assert!(abs_diff_eq!(residual, 0.0, epsilon = 1e-10));
    }
}
```

### Property-Based Tests

```rust
proptest! {
    #[test]
    fn motor_new_always_satisfies_constraint(
        s in -10.0f64..10.0,   // Full range including zero (degenerate case)
        e12 in -10.0f64..10.0,
        e13 in -10.0f64..10.0,
        e23 in -10.0f64..10.0,
        e01 in -10.0f64..10.0,
        e02 in -10.0f64..10.0,
        e03 in -10.0f64..10.0,
    ) {
        let m = Motor::new(s, e12, e13, e23, e01, e02, e03);

        let residual = 2.0 * m.s() * m.e0123()
            - 2.0 * m.e12() * m.e03()
            + 2.0 * m.e13() * m.e02()
            - 2.0 * m.e23() * m.e01();

        prop_assert!(abs_diff_eq!(residual, 0.0, epsilon = 1e-10));
    }
}
```

## Arbitrary Implementations

For types with constraints, the `Arbitrary` implementation **must** use `new()` (not `new_unchecked()`) to guarantee all generated test values satisfy the constraint.

### Example: Motor Arbitrary

```rust
impl<T: Float + Debug> Arbitrary for Motor<T> {
    type Parameters = ();
    type Strategy = BoxedStrategy<Self>;

    fn arbitrary_with(_: Self::Parameters) -> Self::Strategy {
        // Generate 7 independent coefficients, let new() compute e0123
        (
            -10.0f64..10.0,  // s (degenerate when zero, handled by new())
            -10.0f64..10.0,  // e12
            -10.0f64..10.0,  // e13
            -10.0f64..10.0,  // e23
            -10.0f64..10.0,  // e01
            -10.0f64..10.0,  // e02
            -10.0f64..10.0,  // e03
        )
            .prop_map(|(s, e12, e13, e23, e01, e02, e03)| {
                Motor::new(
                    T::from_f64(s),
                    T::from_f64(e12),
                    T::from_f64(e13),
                    T::from_f64(e23),
                    T::from_f64(e01),
                    T::from_f64(e02),
                    T::from_f64(e03),
                )
            })
            .boxed()
    }
}
```

### Key Points

1. **Generate independent coefficients only** - Don't generate the solved-for coefficient
2. **Use `new()` not `new_unchecked()`** - This guarantees all test values satisfy the constraint
3. **Full range is safe** - Since `new()` is total (handles degenerate cases), all coefficient ranges are valid

### Bivector Example

```rust
impl<T: Float + Debug> Arbitrary for Bivector<T> {
    type Parameters = ();
    type Strategy = BoxedStrategy<Self>;

    fn arbitrary_with(_: Self::Parameters) -> Self::Strategy {
        (
            -10.0f64..10.0,  // e12 (degenerate when zero, handled by new())
            -10.0f64..10.0,  // e13
            -10.0f64..10.0,  // e23
            -10.0f64..10.0,  // e01
            -10.0f64..10.0,  // e02
            // e03 computed from constraint
        )
            .prop_map(|(e12, e13, e23, e01, e02)| {
                Bivector::new(
                    T::from_f64(e12),
                    T::from_f64(e13),
                    T::from_f64(e23),
                    T::from_f64(e01),
                    T::from_f64(e02),
                )
            })
            .boxed()
    }
}
```

## Deliverables

- [ ] Add `solve_for` field to TOML spec (`RawTypeSpec`, `TypeSpec`)
- [ ] Implement constraint solving in symbolic module
- [ ] Generate `new()` constructor with computed coefficient
- [ ] Generate `new_checked()` constructor with validation
- [ ] Rename existing constructor to `new_unchecked()`
- [ ] Generate `ConstraintError` type
- [ ] Update `Arbitrary` implementations to use `new()` with independent coefficients
- [ ] Update projective3.toml with `solve_for` fields
- [ ] Regenerate projective3 algebra
- [ ] Unit tests for all three constructors
- [ ] Property-based tests for constraint satisfaction

## Success Criteria

1. `new()` always produces elements that satisfy constraints
2. `new()` is total (never panics) - degenerate cases use canonical values
3. `new_checked()` correctly accepts/rejects based on tolerance
4. `new_unchecked()` compiles and works for all coefficient combinations
5. Generated documentation explains the constraint and constructor behavior
6. `Arbitrary` implementations generate only valid (constraint-satisfying) values

## Dependencies

- PRD-14.11 (Constraint Simplification) - constraint expression format
- PRD-14.12 (Symbolica Verification) - symbolic solver infrastructure

## Future Extensions

1. **Multiple constraints**: Handle types with multiple independent constraints
2. **Non-linear constraints**: Solve quadratic or higher-order constraints
3. **Alternative solutions**: Allow choosing between multiple valid solve targets
4. **Configurable canonical values**: Allow TOML to specify what value to use in degenerate cases (default: zero)
