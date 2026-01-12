# PRD-16.2: Versor Identification and Sandwich Products

**Status**: Draft
**Parent**: PRD-16
**Goal**: Automatically identify versor types and generate efficient sandwich product functions

## Background: What is a Versor?

In Geometric Algebra, a **versor** is an element that can be expressed as a geometric product of invertible vectors. Versors are fundamental for representing transformations.

### Mathematical Definition

A versor `V` satisfies the **versor constraint**:
```
V * rev(V) = scalar  (or pseudoscalar for odd versors)
```

Where `rev(V)` is the reverse (reversion) of `V`.

### Types of Versors

| Versor Type | Grades | Vector Count | Transformation |
|-------------|--------|--------------|----------------|
| **Unit Vector** | [1] | 1 (odd) | Reflection |
| **Rotor** (2D) | [0, 2] | 2 (even) | Rotation in plane |
| **Rotor** (3D) | [0, 2] | 2 (even) | Rotation about axis |
| **Motor** (PGA) | [0, 2, 4] | 2 or 4 (even) | Rigid motion |
| **Flector** (PGA) | [1, 3] | 1 or 3 (odd) | Reflection + translation |

### Versor Properties

1. **Grade Parity**: All grades in a versor have the same parity (all even or all odd)
2. **Invertibility**: `V⁻¹ = rev(V) / (V * rev(V))` for non-null versors
3. **Closure**: Product of two even versors is even; product of even and odd is odd

### The Sandwich Product

Versors transform elements via the **sandwich product**:
```
X' = V * X * rev(V)     (for unit versors)
X' = V * X * V⁻¹        (general case)
```

For **unit versors** (where `V * rev(V) = 1`), the sandwich simplifies to:
```
X' = V * X * rev(V)
```

This is the primary mechanism for applying transformations in GA.

## Problem Statement

### Current Implementation

The current `is_versor_type()` function in `products.rs` is inadequate:

```rust
fn is_versor_type(&self, ty: &TypeSpec) -> bool {
    ty.grades.contains(&0) && ty.grades.iter().all(|g| g % 2 == 0) && ty.name.contains("Rotor")
}
```

**Issues**:
1. Only detects even versors (misses flectors)
2. Requires "Rotor" in the name (misses Motors, custom names)
3. Doesn't verify the versor constraint `V * rev(V) = scalar`
4. Doesn't distinguish unit vs non-unit versors

### Missing Functionality

1. **No automatic versor identification**: Types must be manually marked
2. **No sandwich product generation**: Transformations must be manually implemented
3. **No versor constraint inference**: The Study condition is manually specified

### Relationship to Constraints

The versor constraint relates to the existing constraint system:

| Algebra | Versor Constraint | Current TOML |
|---------|-------------------|--------------|
| Euclidean 3D | `s² + xy² + xz² + yz² = (scalar)` | Not captured |
| PGA 3D (Motor) | `s·e₀₁₂₃ - e₂₃·e₀₁ - e₃₁·e₀₂ - e₁₂·e₀₃ = 0` | `geometric_constraint` |
| PGA 3D (Motor) | `s² + e₁₂² + e₁₃² + e₂₃² = (scalar)` | Not captured |

The **Study condition** in PGA is exactly the versor constraint for motors - it ensures `M * rev(M)` is a scalar (not scalar + pseudoscalar).

## Solution

### Phase 1: Versor Classification

#### 1.1 Define Versor Parity

```rust
/// Versor classification by parity.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum VersorParity {
    /// Even versor: grades 0, 2, 4, ... (rotors, motors)
    Even,
    /// Odd versor: grades 1, 3, 5, ... (reflectors, flectors)
    Odd,
}
```

#### 1.2 Versor Identification Algorithm

A type is a **potential versor** if:
1. All its grades have the same parity
2. The grade set is closed under the sandwich product (grades that appear when applying the versor to itself)

```rust
/// Determines if a type has versor-compatible grades.
///
/// Returns the parity if all grades have the same parity, None otherwise.
pub fn versor_parity(grades: &[usize]) -> Option<VersorParity> {
    if grades.is_empty() {
        return None;
    }

    let first_parity = grades[0] % 2;
    if grades.iter().all(|&g| g % 2 == first_parity) {
        Some(if first_parity == 0 {
            VersorParity::Even
        } else {
            VersorParity::Odd
        })
    } else {
        None
    }
}
```

#### 1.3 Symbolic Verification

To confirm a type is truly a versor, verify symbolically that `V * rev(V)` produces only scalar (grade 0) output:

```rust
/// Verifies the versor constraint symbolically.
///
/// A type is a versor if V * rev(V) produces only grade-0 (or grade-n for odd versors).
pub fn verify_versor_constraint(
    ty: &TypeSpec,
    algebra: &Algebra,
) -> VersorVerification {
    let parity = match versor_parity(&ty.grades) {
        Some(p) => p,
        None => return VersorVerification::NotVersor("mixed grade parity"),
    };

    // Compute V * rev(V) symbolically
    let v_rev_v = compute_symbolic_product(ty, ty, ProductKind::Geometric, |field| {
        // Apply reverse sign: (-1)^(k(k-1)/2) for grade k
        let k = field.grade;
        let sign = if (k * (k - 1) / 2) % 2 == 0 { 1 } else { -1 };
        (field.name.clone(), sign)
    });

    // Check which grades have non-zero output
    let output_grades: HashSet<usize> = v_rev_v.non_zero_grades();

    match parity {
        VersorParity::Even => {
            // Even versors should produce only scalar (grade 0)
            if output_grades == hashset!{0} {
                VersorVerification::Versor(VersorInfo {
                    parity: VersorParity::Even,
                    is_unit: false, // Need additional constraint for unit
                    constraint_expr: build_versor_constraint_expr(ty, &v_rev_v),
                })
            } else if output_grades.iter().all(|&g| g == 0 || g == algebra.dim()) {
                // Produces scalar + pseudoscalar - needs Study constraint
                VersorVerification::ConstrainedVersor {
                    parity: VersorParity::Even,
                    required_constraint: build_versor_constraint_expr(ty, &v_rev_v),
                }
            } else {
                VersorVerification::NotVersor("V * rev(V) produces non-scalar grades")
            }
        }
        VersorParity::Odd => {
            // Odd versors produce pseudoscalar (grade n) or scalar
            let n = algebra.dim();
            if output_grades.iter().all(|&g| g == 0 || g == n) {
                VersorVerification::Versor(VersorInfo {
                    parity: VersorParity::Odd,
                    is_unit: false,
                    constraint_expr: None,
                })
            } else {
                VersorVerification::NotVersor("V * rev(V) produces unexpected grades")
            }
        }
    }
}
```

### Phase 2: TOML Specification Extension

#### 2.1 Explicit Versor Marking

Allow explicit versor annotation in TOML:

```toml
[types.Motor]
grades = [0, 2, 4]
fields = ["s", "e12", "e13", "e23", "e01", "e02", "e03", "e0123"]
versor = true  # Explicit versor marking

# The Study condition is automatically recognized as the versor constraint
geometric_constraint = "2*s*e0123 - 2*e12*e03 + 2*e13*e02 - 2*e23*e01 = 0"
geometric_solve_for = "e0123"

[types.Flector]
grades = [1, 3]
fields = ["e1", "e2", "e3", "e0", "e012", "e013", "e023", "e123"]
versor = true  # Odd versor
```

#### 2.2 Sandwich Product Targets

Specify which types a versor can transform:

```toml
[types.Motor.sandwich]
# Motor can transform these types via sandwich product
targets = ["Point", "Line", "Plane", "Motor"]

# Or auto-detect (default behavior)
# targets = "auto"
```

#### 2.3 Raw TOML Structure

Update `raw.rs`:

```rust
#[derive(Debug, Deserialize)]
pub struct RawTypeSpec {
    pub grades: Vec<usize>,
    pub fields: Vec<String>,
    // ... existing fields ...

    /// Whether this type is a versor (for sandwich products).
    #[serde(default)]
    pub versor: bool,

    /// Types this versor can transform via sandwich product.
    #[serde(default)]
    pub sandwich_targets: Option<SandwichTargets>,
}

#[derive(Debug, Deserialize)]
#[serde(untagged)]
pub enum SandwichTargets {
    /// Specific list of target type names.
    List(Vec<String>),
    /// Auto-detect targets based on grade compatibility.
    Auto(String),  // "auto"
}
```

### Phase 3: Versor Information in IR

Update the intermediate representation:

```rust
/// Information about a versor type.
#[derive(Debug, Clone)]
pub struct VersorInfo {
    /// Versor parity (even or odd).
    pub parity: VersorParity,

    /// Whether this is a unit versor (V * rev(V) = 1).
    pub is_unit: bool,

    /// The versor constraint expression (if one is required).
    /// For motors, this is the Study condition.
    pub constraint: Option<String>,

    /// Types this versor can transform.
    pub sandwich_targets: Vec<String>,
}

/// Extended type specification with versor info.
#[derive(Debug, Clone)]
pub struct TypeSpec {
    pub name: String,
    pub grades: Vec<usize>,
    pub fields: Vec<FieldSpec>,
    pub constraints: Vec<UserConstraint>,

    /// Versor information (if this type is a versor).
    pub versor: Option<VersorInfo>,
}
```

### Phase 4: Sandwich Product Generation

#### 4.1 Product Generator Extension

Add sandwich product generation to `ProductGenerator`:

```rust
impl<'a> ProductGenerator<'a> {
    /// Generates all sandwich product functions.
    pub fn generate_sandwich_products(&self) -> TokenStream {
        let mut products = Vec::new();

        for ty in &self.spec.types {
            if let Some(ref versor_info) = ty.versor {
                for target_name in &versor_info.sandwich_targets {
                    if let Some(target) = self.find_type(target_name) {
                        let product = self.generate_sandwich_product(ty, target, versor_info);
                        products.push(product);
                    }
                }
            }
        }

        quote! { #(#products)* }
    }

    /// Generates a single sandwich product function.
    fn generate_sandwich_product(
        &self,
        versor: &TypeSpec,
        operand: &TypeSpec,
        info: &VersorInfo,
    ) -> TokenStream {
        let versor_name = format_ident!("{}", versor.name);
        let operand_name = format_ident!("{}", operand.name);
        let fn_name = format_ident!(
            "sandwich_{}_{}",
            versor.name.to_lowercase(),
            operand.name.to_lowercase()
        );

        // Compute sandwich product terms: V * X * rev(V)
        // This is NOT computed as two separate products - we expand fully
        let terms = self.compute_sandwich_terms(versor, operand);

        let output_type = self.determine_sandwich_output(versor, operand);
        let output_name = format_ident!("{}", output_type.name);

        let field_exprs: Vec<TokenStream> = output_type
            .fields
            .iter()
            .map(|field| {
                let expr = self.terms_to_expression(&terms, field.blade_index);
                expr
            })
            .collect();

        let doc = format!(
            "Computes the sandwich product `{} * {} * rev({})`.\n\n\
             Transforms a {} by the {} versor.",
            versor.name, operand.name, versor.name,
            operand.name, versor.name
        );

        quote! {
            #[doc = #doc]
            #[inline]
            pub fn #fn_name(v: &#versor_name<T>, x: &#operand_name<T>) -> #output_name<T> {
                #output_name::new(#(#field_exprs),*)
            }
        }
    }

    /// Computes sandwich product terms: V * X * rev(V).
    ///
    /// Instead of computing (V * X) then (* rev(V)), we expand the full
    /// trilinear product directly. This allows for more simplification
    /// and avoids intermediate storage.
    fn compute_sandwich_terms(
        &self,
        versor: &TypeSpec,
        operand: &TypeSpec,
    ) -> Vec<SandwichTerm> {
        let mut terms = Vec::new();

        // For each output blade:
        // Sum over all (v_i, x_j, rev(v)_k) where v_i * x_j * rev(v)_k contributes

        for v_field in &versor.fields {
            for x_field in &operand.fields {
                for w_field in &versor.fields {
                    // w is rev(v), so it has a sign based on grade
                    let w_grade = w_field.grade;
                    let rev_sign = if (w_grade * (w_grade - 1) / 2) % 2 == 0 { 1 } else { -1 };

                    // Compute v * x * w
                    let vx_blade = v_field.blade_index ^ x_field.blade_index;
                    let vx_sign = self.compute_product_sign(v_field.blade_index, x_field.blade_index);

                    let result_blade = vx_blade ^ w_field.blade_index;
                    let result_sign = vx_sign * self.compute_product_sign(vx_blade, w_field.blade_index) * rev_sign;

                    if result_sign != 0 {
                        terms.push(SandwichTerm {
                            output_blade: result_blade,
                            v_field: v_field.name.clone(),
                            x_field: x_field.name.clone(),
                            w_field: w_field.name.clone(),
                            sign: result_sign,
                        });
                    }
                }
            }
        }

        // Collect and simplify terms for each output blade
        self.simplify_sandwich_terms(terms)
    }
}
```

#### 4.2 Method Generation on Types

Also generate methods on the versor type itself:

```rust
impl<T: Float> Motor<T> {
    /// Transforms a point by this motor.
    ///
    /// Computes `self * point * self.reverse()`.
    #[inline]
    pub fn transform_point(&self, point: &Point<T>) -> Point<T> {
        sandwich_motor_point(self, point)
    }

    /// Transforms a line by this motor.
    #[inline]
    pub fn transform_line(&self, line: &Line<T>) -> Line<T> {
        sandwich_motor_line(self, line)
    }

    /// Transforms a plane by this motor.
    #[inline]
    pub fn transform_plane(&self, plane: &Plane<T>) -> Plane<T> {
        sandwich_motor_plane(self, plane)
    }
}
```

### Phase 5: Auto-Detection of Sandwich Targets

When `sandwich_targets = "auto"`, determine valid targets:

```rust
/// Determines which types can be transformed by a versor.
///
/// A type can be sandwiched if V * X * rev(V) produces the same grades as X.
fn infer_sandwich_targets(
    versor: &TypeSpec,
    all_types: &[TypeSpec],
    algebra: &Algebra,
) -> Vec<String> {
    let mut targets = Vec::new();

    for candidate in all_types {
        // Skip the versor itself (self-sandwich is composition)
        if candidate.name == versor.name {
            continue;
        }

        // Compute output grades of V * X * rev(V)
        let output_grades = compute_sandwich_output_grades(
            &versor.grades,
            &candidate.grades,
            algebra.dim(),
        );

        // Valid target if output grades match input grades
        if output_grades == candidate.grades.iter().cloned().collect::<HashSet<_>>() {
            targets.push(candidate.name.clone());
        }
    }

    targets
}

/// Computes the output grades of a sandwich product.
fn compute_sandwich_output_grades(
    versor_grades: &[usize],
    operand_grades: &[usize],
    dim: usize,
) -> HashSet<usize> {
    let mut output = HashSet::new();

    for &vg in versor_grades {
        for &xg in operand_grades {
            for &wg in versor_grades {
                // Compute possible output grades from v * x * w
                // Grade rules: |g1 - g2| <= g_out <= min(g1 + g2, 2n - g1 - g2)
                let vx_grades = product_grade_range(vg, xg, dim);
                for vxg in vx_grades {
                    let vxw_grades = product_grade_range(vxg, wg, dim);
                    output.extend(vxw_grades);
                }
            }
        }
    }

    output
}
```

## TOML Examples

### 3D Euclidean Rotor

```toml
[types.Rotor]
grades = [0, 2]
fields = ["s", "xy", "xz", "yz"]
versor = true

[types.Rotor.sandwich]
targets = ["Vector", "Bivector", "Rotor"]

[[types.Rotor.constraints]]
name = "unit"
expression = "s*s + xy*xy + xz*xz + yz*yz = 1"
solve_for = "s"
```

### PGA Motor

```toml
[types.Motor]
grades = [0, 2, 4]
fields = ["s", "e12", "e13", "e23", "e01", "e02", "e03", "e0123"]
versor = true

# Study condition - ensures M * rev(M) is scalar only
geometric_constraint = "2*s*e0123 - 2*e12*e03 + 2*e13*e02 - 2*e23*e01 = 0"
geometric_solve_for = "e0123"

[types.Motor.sandwich]
targets = ["Point", "Line", "Plane", "Motor"]
```

### PGA Flector (Odd Versor)

```toml
[types.Flector]
grades = [1, 3]
fields = ["e1", "e2", "e3", "e0", "e012", "e013", "e023", "e123"]
versor = true

[types.Flector.sandwich]
targets = ["Point", "Line", "Plane"]
```

## Testing Strategy

### Versor Identification Tests

```rust
#[test]
fn identify_even_versor() {
    let rotor_grades = vec![0, 2];
    assert_eq!(versor_parity(&rotor_grades), Some(VersorParity::Even));

    let motor_grades = vec![0, 2, 4];
    assert_eq!(versor_parity(&motor_grades), Some(VersorParity::Even));
}

#[test]
fn identify_odd_versor() {
    let vector_grades = vec![1];
    assert_eq!(versor_parity(&vector_grades), Some(VersorParity::Odd));

    let flector_grades = vec![1, 3];
    assert_eq!(versor_parity(&flector_grades), Some(VersorParity::Odd));
}

#[test]
fn reject_mixed_parity() {
    let mixed_grades = vec![0, 1, 2];  // Mixed parity
    assert_eq!(versor_parity(&mixed_grades), None);
}
```

### Sandwich Product Tests

```rust
proptest! {
    #[test]
    fn sandwich_preserves_grade(
        r in any::<UnitRotor<f64>>(),
        v in any::<Vector<f64>>()
    ) {
        let transformed = r.transform_vector(&v);
        // Vector transforms to vector
        // (implicitly tested by return type)
    }

    #[test]
    fn sandwich_preserves_norm(
        r in any::<UnitRotor<f64>>(),
        v in any::<Vector<f64>>()
    ) {
        let transformed = r.transform_vector(&v);
        prop_assert!(abs_diff_eq!(v.norm(), transformed.norm(), epsilon = ABS_DIFF_EQ_EPS));
    }

    #[test]
    fn motor_sandwich_is_rigid(
        m in any::<UnitMotor<f64>>(),
        p1 in any::<Point<f64>>(),
        p2 in any::<Point<f64>>()
    ) {
        let p1_t = m.transform_point(&p1);
        let p2_t = m.transform_point(&p2);

        // Distance between points is preserved
        let dist_before = distance(&p1, &p2);
        let dist_after = distance(&p1_t, &p2_t);
        prop_assert!(abs_diff_eq!(dist_before, dist_after, epsilon = ABS_DIFF_EQ_EPS));
    }
}
```

### Versor Constraint Tests

```rust
#[test]
fn motor_satisfies_study_condition() {
    // Motor from axis-angle should satisfy Study condition
    let axis = Vector::new(0.0, 0.0, 1.0);
    let angle = std::f64::consts::PI / 4.0;
    let motor = Motor::from_axis_angle(&axis, angle);

    // Study condition: s*e0123 - e23*e01 - e31*e02 - e12*e03 = 0
    let residual = motor.s() * motor.e0123()
        - motor.e23() * motor.e01()
        - motor.e31() * motor.e02()
        - motor.e12() * motor.e03();

    assert!(abs_diff_eq!(residual, 0.0, epsilon = 1e-10));
}
```

## Deliverables

- [ ] Add `VersorParity` enum in `algebra/versor.rs`
- [ ] Add `versor_parity()` function
- [ ] Add symbolic versor constraint verification
- [ ] Update `RawTypeSpec` with `versor` and `sandwich_targets` fields
- [ ] Update `TypeSpec` IR with `VersorInfo`
- [ ] Update parser to process versor information
- [ ] Implement `generate_sandwich_products()` in ProductGenerator
- [ ] Implement `compute_sandwich_terms()` for trilinear expansion
- [ ] Add sandwich target auto-detection
- [ ] Generate `transform_*` methods on versor types
- [ ] Add comprehensive tests

## Success Criteria

1. Rotors, Motors, Flectors correctly identified as versors
2. Sandwich products generated for all versor-target pairs
3. Generated sandwich products pass algebraic property tests
4. Study condition recognized as versor constraint

## Files Changed

| File | Action | Description |
|------|--------|-------------|
| `crates/clifford-codegen/src/algebra/mod.rs` | Update | Re-export versor module |
| `crates/clifford-codegen/src/algebra/versor.rs` | Create | Versor identification |
| `crates/clifford-codegen/src/spec/raw.rs` | Update | Add versor fields |
| `crates/clifford-codegen/src/spec/ir.rs` | Update | Add VersorInfo |
| `crates/clifford-codegen/src/spec/parser.rs` | Update | Parse versor info |
| `crates/clifford-codegen/src/codegen/products.rs` | Update | Sandwich generation |
| `algebras/euclidean3.toml` | Update | Add versor annotations |
| `algebras/projective3.toml` | Update | Add versor annotations |

## References

- [Versors in Geometric Algebra (Wikipedia)](https://en.wikipedia.org/wiki/Versor#Versors_in_geometric_algebra)
- [Geometric Algebra for Computer Science](https://geometricalgebra.org/) - Chapter on versors
- [ganja.js versor implementation](https://github.com/enkimute/ganja.js)
