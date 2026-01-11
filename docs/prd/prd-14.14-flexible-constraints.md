# PRD-14.14: Flexible Constraint System

**Status**: Draft
**Parent**: [PRD-14: Geometric Algebra Code Generator](prd-14-codegen.md)
**Goal**: Support automatic constraint dependency detection and user-defined constraints for accurate degree-of-freedom handling

## Overview

The current constraint system (PRD-14.13) assumes geometric and antiproduct constraints are independent, requiring separate `solve_for` fields for each. However, in 4D PGA, these constraints are mathematically equivalent (dependent), producing identical expressions. Additionally, some types require constraints beyond the automatic geometric/antiproduct constraints, such as unit norm for proper rigid body transformations.

This PRD addresses:
1. **Automatic dependency detection**: Determine if constraints are independent or equivalent
2. **User-defined constraints**: Allow specifying additional constraints like unit norm
3. **Correct DOF calculation**: Generate constructors with the correct number of independent parameters

## Problem Statement

### Identical Constraints in PGA

In 3D PGA (Cl(3,0,1)), the Motor type has:
- **Geometric constraint**: `u * ũ = scalar` → `2*s*e0123 - 2*e12*e03 + 2*e13*e02 - 2*e23*e01 = 0`
- **Antiproduct constraint**: `u ⊟ ũ̃ = antiscalar` → `2*s*e0123 - 2*e12*e03 + 2*e13*e02 - 2*e23*e01 = 0`

These are **identical expressions**. This is not a bug—it's a consequence of the algebra's structure where grade-2 elements are self-dual in 4D.

### DOF Mismatch

| Type | Components | Algebraic Constraints | Current DOF | Expected DOF |
|------|------------|----------------------|-------------|--------------|
| Motor | 8 | 1 (Study condition) | 7 | 6 (rigid body) |
| Bivector (Line) | 6 | 1 (Plücker condition) | 5 | 4 (line in 3D) |

The missing constraint is **unit norm** (`‖M‖ = 1`), which reduces:
- Motor: 7 → 6 DOF (true rigid body transformation)
- Line: 5 → 4 DOF (unique line representation)

### Current Limitations

1. The generator doesn't detect that `geometric_constraint == antiproduct_constraint`
2. No way to specify additional constraints like unit norm
3. Tests expect 2 independent constraints but only 1 exists

## Proposed Solution

### 1. Automatic Constraint Dependency Detection

Compare the normalized constraint expressions:

```rust
impl TypeSpec {
    /// Returns true if geometric and antiproduct constraints are independent.
    /// Returns false if they are equivalent (same expression) or if only one exists.
    pub fn has_independent_constraints(&self) -> bool {
        match (&self.geometric_constraint, &self.antiproduct_constraint) {
            (Some(gc), Some(ac)) => {
                // Normalize and compare expressions
                normalize_constraint(gc) != normalize_constraint(ac)
            }
            _ => false,
        }
    }
}

fn normalize_constraint(expr: &str) -> String {
    // 1. Parse into symbolic form
    // 2. Sort terms canonically
    // 3. Normalize coefficients
    // 4. Return canonical string representation
}
```

### 2. User-Defined Constraints

Add a new `constraints` field to the TOML spec for additional constraints:

```toml
[types.Motor]
grades = [0, 2, 4]
fields = ["s", "e12", "e13", "e23", "e01", "e02", "e03", "e0123"]
geometric_constraint = "2*s*e0123 - 2*e12*e03 + 2*e13*e02 - 2*e23*e01 = 0"
antiproduct_constraint = "2*s*e0123 - 2*e12*e03 + 2*e13*e02 - 2*e23*e01 = 0"
geometric_solve_for = "e0123"

# Additional user-defined constraints
[[types.Motor.constraints]]
name = "unit"
expression = "s*s + e12*e12 + e13*e13 + e23*e23 = 1"
solve_for = "s"        # Solve for s, making new() return Option<Self>
sign = "positive"      # Use positive square root
```

### 3. Constraint Model

Every constraint is treated uniformly:

| Has `solve_for`? | Behavior |
|------------------|----------|
| Yes | Compute that field in `new()`, reduce parameter count by 1 |
| No | Validate in `new_checked()`, optionally generate enforcement method |

Any constraint can have `solve_for`, including quadratic ones:

```toml
[[types.Motor.constraints]]
name = "unit"
expression = "s*s + e12*e12 + e13*e13 + e23*e23 = 1"
solve_for = "s"
sign = "positive"  # Convention for ± ambiguity (default: positive)
```

This generates:
```rust
pub fn new(e12: T, e13: T, e23: T, e01: T, e02: T, e03: T) -> Option<Self> {
    let sum_sq = e12 * e12 + e13 * e13 + e23 * e23;
    if sum_sq > T::one() {
        return None;  // No real solution
    }
    let s = (T::one() - sum_sq).sqrt();  // positive root
    // ... compute e0123 from Study condition ...
}
```

**Domain handling**: When solving produces restricted domains (e.g., `√(1-x²)` requires `|x| ≤ 1`), `new()` returns `Option<Self>` instead of `Self`.

If the user prefers normalization over domain-restricted construction, they can omit `solve_for` and use `enforce`:

```toml
[[types.Motor.constraints]]
name = "unit"
expression = "s*s + e12*e12 + e13*e13 + e23*e23 = 1"
enforce = "normalize"  # Generate normalize() method instead
```

### 4. Constructor Generation Rules

| Solvable Constraints | `new()` Params | Return Type |
|---------------------|----------------|-------------|
| 0 | All N fields | `Self` |
| K (all linear) | N-K fields | `Self` |
| K (any with domain restriction) | N-K fields | `Option<Self>` |

A constraint is "solvable" if it has a `solve_for` field. The generator:
1. Collects all solvable constraints
2. Reduces parameters by K (one per solvable constraint)
3. Returns `Option<Self>` if any constraint has domain restrictions (e.g., square roots)

## TOML Specification

### Extended Syntax

```toml
[types.Motor]
grades = [0, 2, 4]
fields = ["s", "e12", "e13", "e23", "e01", "e02", "e03", "e0123"]
description = "Rigid body transformation (rotation + translation)"

# Automatically derived constraints (from discovery)
geometric_constraint = "2*s*e0123 - 2*e12*e03 + 2*e13*e02 - 2*e23*e01 = 0"
antiproduct_constraint = "2*s*e0123 - 2*e12*e03 + 2*e13*e02 - 2*e23*e01 = 0"
geometric_solve_for = "e0123"
# antiproduct_solve_for not needed since constraints are dependent

# User-defined additional constraints
[[types.Motor.constraints]]
name = "unit"
description = "Unit motor for proper rigid transformations"
expression = "s*s + e12*e12 + e13*e13 + e23*e23 = 1"
# No solve_for → validated only
enforce = "normalize"  # Generate normalize() method

[types.Line]
grades = [2]
fields = ["e12", "e13", "e23", "e01", "e02", "e03"]
description = "Line in 3D space (Plücker coordinates)"

geometric_constraint = "-2*e12*e03 + 2*e13*e02 - 2*e23*e01 = 0"
antiproduct_constraint = "-2*e12*e03 + 2*e13*e02 - 2*e23*e01 = 0"
geometric_solve_for = "e03"

[[types.Line.constraints]]
name = "unit"
description = "Normalized line (direction has unit length)"
expression = "e12*e12 + e13*e13 + e23*e23 = 1"
enforce = "normalize"
```

### Constraint Block Schema

```toml
[[types.TypeName.constraints]]
name = "constraint_name"           # Required: identifier for this constraint
description = "Human readable"     # Optional: for documentation
expression = "a + b + c = 0"       # Required: constraint equation
solve_for = "field_name"           # Optional: if present, solve for this field in new()
sign = "positive"                  # Optional: for ± ambiguity, "positive" (default) or "negative"
enforce = "normalize"              # Optional: generate enforcement method (e.g., normalize())
```

### Additional Solvable Constraint Example

A user might want to add a custom solvable constraint:

```toml
[types.Line]
grades = [2]
fields = ["e12", "e13", "e23", "e01", "e02", "e03"]

# Standard Plücker constraint (auto-derived)
geometric_constraint = "-2*e12*e03 + 2*e13*e02 - 2*e23*e01 = 0"
geometric_solve_for = "e03"

# User adds a custom solvable constraint
[[types.Line.constraints]]
name = "fixed_moment_ratio"
description = "Application-specific: moment x equals moment y"
expression = "e01 - e02 = 0"
solve_for = "e02"  # Compute e02 from e01
```

This generates `new()` with 4 parameters (e12, e13, e23, e01), computing:
- `e03` from Plücker constraint
- `e02` from custom constraint

## Generated Code

### Motor with Proper Constraints

```rust
/// 3D PGA Motor (rigid body transformation).
///
/// Motors have 8 coefficients but represent transformations with 6 degrees
/// of freedom. This is achieved through two constraints:
///
/// 1. **Study condition** (linear): `2s·e₀₁₂₃ - 2e₁₂·e₀₃ + 2e₁₃·e₀₂ - 2e₂₃·e₀₁ = 0`
/// 2. **Unit norm** (quadratic): `s² + e₁₂² + e₁₃² + e₂₃² = 1`
///
/// The `new()` constructor enforces the Study condition automatically.
/// For unit motors, use `UnitMotor` or normalize after construction.
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
    /// The `e0123` coefficient is computed from the Study condition.
    /// This guarantees the motor is a valid geometric entity, but
    /// NOT necessarily a unit motor.
    ///
    /// For unit motors (proper rigid transformations), either:
    /// - Use `UnitMotor::from_*` constructors
    /// - Call `.normalize()` after construction
    #[inline]
    pub fn new(s: T, e12: T, e13: T, e23: T, e01: T, e02: T, e03: T) -> Self {
        let e0123 = if s.abs() > T::epsilon() {
            (e12 * e03 - e13 * e02 + e23 * e01) / s
        } else {
            T::zero()
        };
        Self { s, e12, e13, e23, e01, e02, e03, e0123 }
    }

    /// Creates a Motor with full validation of all constraints.
    ///
    /// Validates:
    /// - Study condition (linear): residual < tolerance
    /// - Unit norm (quadratic): |norm - 1| < tolerance (if check_unit = true)
    pub fn new_checked(
        s: T, e12: T, e13: T, e23: T,
        e01: T, e02: T, e03: T, e0123: T,
        tolerance: T,
        check_unit: bool,
    ) -> Result<Self, ConstraintError> {
        // Check Study condition
        let study_residual = T::TWO * s * e0123
            - T::TWO * e12 * e03
            + T::TWO * e13 * e02
            - T::TWO * e23 * e01;

        if study_residual.abs() > tolerance {
            return Err(ConstraintError::Study { residual: study_residual.to_f64().unwrap_or(0.0) });
        }

        // Optionally check unit norm
        if check_unit {
            let norm_sq = s * s + e12 * e12 + e13 * e13 + e23 * e23;
            let unit_residual = norm_sq - T::one();
            if unit_residual.abs() > tolerance {
                return Err(ConstraintError::Unit { residual: unit_residual.to_f64().unwrap_or(0.0) });
            }
        }

        Ok(Self { s, e12, e13, e23, e01, e02, e03, e0123 })
    }

    /// Creates a Motor without any validation.
    #[inline]
    pub fn new_unchecked(
        s: T, e12: T, e13: T, e23: T,
        e01: T, e02: T, e03: T, e0123: T,
    ) -> Self {
        Self { s, e12, e13, e23, e01, e02, e03, e0123 }
    }

    /// Returns the rotor norm squared: s² + e₁₂² + e₁₃² + e₂₃².
    #[inline]
    pub fn rotor_norm_squared(&self) -> T {
        self.s * self.s + self.e12 * self.e12 + self.e13 * self.e13 + self.e23 * self.e23
    }

    /// Returns the rotor norm.
    #[inline]
    pub fn rotor_norm(&self) -> T {
        self.rotor_norm_squared().sqrt()
    }

    /// Returns true if this motor satisfies the unit constraint within tolerance.
    #[inline]
    pub fn is_unit(&self, tolerance: T) -> bool {
        (self.rotor_norm_squared() - T::one()).abs() < tolerance
    }

    /// Normalizes this motor to have unit rotor norm.
    ///
    /// Returns `None` if the rotor part is zero.
    pub fn normalize(&self) -> Option<Self> {
        let norm = self.rotor_norm();
        if norm < T::epsilon() {
            return None;
        }
        let inv_norm = T::one() / norm;
        Some(Self::new_unchecked(
            self.s * inv_norm,
            self.e12 * inv_norm,
            self.e13 * inv_norm,
            self.e23 * inv_norm,
            self.e01 * inv_norm,
            self.e02 * inv_norm,
            self.e03 * inv_norm,
            self.e0123 * inv_norm,
        ))
    }
}
```

### Extended ConstraintError

```rust
/// Error returned when a geometric constraint is not satisfied.
#[derive(Debug, Clone)]
pub enum ConstraintError {
    /// Study condition violated (linear constraint).
    Study { residual: f64 },
    /// Unit norm violated (quadratic constraint).
    Unit { residual: f64 },
    /// Plücker condition violated (for lines).
    Plucker { residual: f64 },
    /// Generic constraint violation.
    Custom { name: &'static str, residual: f64 },
}

impl std::fmt::Display for ConstraintError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::Study { residual } => {
                write!(f, "Study condition violated (residual: {:.2e})", residual)
            }
            Self::Unit { residual } => {
                write!(f, "Unit norm violated (residual: {:.2e})", residual)
            }
            Self::Plucker { residual } => {
                write!(f, "Plücker condition violated (residual: {:.2e})", residual)
            }
            Self::Custom { name, residual } => {
                write!(f, "{} constraint violated (residual: {:.2e})", name, residual)
            }
        }
    }
}

impl std::error::Error for ConstraintError {}
```

## Implementation Plan

### Phase 1: Constraint Dependency Detection

1. Add `normalize_constraint()` function to canonicalize constraint expressions
2. Update `TypeSpec::has_independent_constraints()` to compare normalized forms
3. Update validation to not require two `solve_for` fields when constraints are dependent
4. Fix tests that expect 2 constraints when only 1 exists

### Phase 2: User-Defined Constraints

1. Add `constraints` array to `RawTypeSpec` in TOML parsing
2. Create `UserConstraint` struct with name, expression, type, etc.
3. Update `TypeSpec` to carry user constraints
4. Update validation to check user constraint syntax

### Phase 3: Code Generation Updates

1. Generate `ConstraintError` enum with variants for each constraint type
2. Generate `new_checked()` with optional constraint validation
3. Generate constraint-specific methods (`is_unit()`, `satisfies_study()`, etc.)
4. Generate `normalize()` method for types with unit constraint
5. Update documentation to explain all constraints

### Phase 4: Discovery Tool Updates

1. Update discovery to report constraint dependency
2. Suggest user constraints based on algebra structure:
   - 4D algebras → likely unit norm constraint
   - Grade-2 types → likely Plücker-like constraint
3. Generate TOML template with constraint placeholders

## Testing

### Unit Tests

```rust
#[test]
fn constraint_dependency_detection() {
    let spec = parse_spec("algebras/projective3.toml").unwrap();
    let motor = spec.types.iter().find(|t| t.name == "Motor").unwrap();

    // Constraints should be detected as dependent (identical)
    assert!(!motor.has_independent_constraints());
}

#[test]
fn motor_new_satisfies_study_condition() {
    let m = Motor::new(1.0, 0.1, 0.2, 0.3, 0.01, 0.02, 0.03);

    let residual = 2.0 * m.s() * m.e0123()
        - 2.0 * m.e12() * m.e03()
        + 2.0 * m.e13() * m.e02()
        - 2.0 * m.e23() * m.e01();

    assert!(abs_diff_eq!(residual, 0.0, epsilon = 1e-10));
}

#[test]
fn motor_new_checked_validates_unit() {
    // Create unit motor
    let m = Motor::new(1.0, 0.0, 0.0, 0.0, 0.1, 0.2, 0.3);

    // Should pass unit check (s=1, bivector=0 → norm=1)
    let result = Motor::new_checked(
        m.s(), m.e12(), m.e13(), m.e23(),
        m.e01(), m.e02(), m.e03(), m.e0123(),
        1e-10,
        true,  // check_unit
    );
    assert!(result.is_ok());

    // Non-unit should fail
    let result = Motor::new_checked(
        2.0, 0.0, 0.0, 0.0,  // s=2 → norm=2 ≠ 1
        0.1, 0.2, 0.3, 0.0,
        1e-10,
        true,
    );
    assert!(matches!(result, Err(ConstraintError::Unit { .. })));
}
```

### Property-Based Tests

```rust
proptest! {
    #[test]
    fn normalized_motor_is_unit(m in any::<Motor<f64>>()) {
        if let Some(normalized) = m.normalize() {
            prop_assert!(normalized.is_unit(1e-10));
        }
    }

    #[test]
    fn motor_new_always_satisfies_study(
        s in -10.0f64..10.0,
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

## Deliverables

- [ ] Implement `normalize_constraint()` for canonical comparison
- [ ] Update `has_independent_constraints()` to use normalization
- [ ] Add `constraints` array to TOML spec (`RawTypeSpec`)
- [ ] Create `UserConstraint` struct
- [ ] Generate `ConstraintError` enum with named variants
- [ ] Generate `new_checked()` with configurable validation
- [ ] Generate constraint helper methods (`is_unit()`, etc.)
- [ ] Generate `normalize()` for types with unit constraint
- [ ] Update discovery tool to detect constraint dependency
- [ ] Update projective3.toml with user constraints
- [ ] Documentation explaining constraint system

## Success Criteria

1. Identical constraints correctly detected as dependent
2. User can specify additional constraints in TOML
3. `new()` takes correct number of independent parameters
4. `new_checked()` validates all specified constraints
5. Unit motors can be created via normalization
6. Generated documentation explains all constraints and DOF

## Dependencies

- PRD-14.13 (Constraint-Enforcing Constructors) - base constraint solving
- PRD-14.11 (Constraint Simplification) - expression normalization utilities

## Known Limitations

### Constraint Ordering

The order of constraints in the TOML determines the order they are solved. If constraint B references a field solved by constraint A, then A must appear before B in the TOML:

```toml
# CORRECT: unit solved first (computes s), then study (uses s)
[[types.Motor.constraints]]
name = "unit"
expression = "s*s + e12*e12 + e13*e13 + e23*e23 = 1"
solve_for = "s"

[[types.Motor.constraints]]
name = "study"
expression = "2*s*e0123 - 2*e12*e03 + 2*e13*e02 - 2*e23*e01 = 0"
solve_for = "e0123"
```

**Future improvement**: Automatic dependency detection and topological ordering.

## Future Extensions

1. **Automatic unit type generation**: Generate `UnitMotor`, `UnitLine` wrappers
2. **Constraint algebra**: Verify constraint compatibility (e.g., Study + unit = 6 DOF)
3. **Non-linear solving**: Use numerical methods for quadratic constraints
4. **Constraint preservation tracking**: Automatically determine which operations preserve which constraints
5. **Automatic constraint ordering**: Detect dependencies between constraints and order them correctly
