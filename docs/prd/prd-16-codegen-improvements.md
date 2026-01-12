# PRD-16: Codegen Improvements - Blade Ordering, Versors, and Expression Simplification

**Status**: Draft
**Goal**: Fix blade ordering inconsistencies, add versor identification with sandwich product generation, and use Symbolica for expression simplification

## Sub-PRDs

| PRD | Status | Description |
|-----|--------|-------------|
| [PRD-16.1](prd-16.1-blade-ordering.md) | Draft | Blade Ordering Audit and Validation |
| [PRD-16.2](prd-16.2-versors.md) | Draft | Versor Identification and Sandwich Products |
| [PRD-16.3](prd-16.3-expression-simplification.md) | Draft | Expression Simplification with Symbolica |

## Problem Statement

The current code generator (`clifford-codegen`) has three significant issues affecting correctness and code quality:

### Issue 1: Blade Ordering vs Canonical Ordering

The generator uses bitmask indices for blade representation (`Blade::from_index()`), which is inherently canonical. However, there's a potential mismatch between:

1. **Canonical blade ordering**: Determined by bitmask index (e.g., `e12 = 0b011 = 3`, `e13 = 0b101 = 5`, `e23 = 0b110 = 6`)
2. **Type field ordering**: Order of fields in the TOML specification
3. **Constructor parameter ordering**: Order of parameters in `new()` and `new_unchecked()`

For example, in PGA3D the Bivector type has fields:
```toml
[types.Bivector]
fields = ["e12", "e13", "e23", "e01", "e02", "e03"]
```

The blade indices are:
- `e12` = `0b0011` = 3 (grade 2)
- `e13` = `0b0101` = 5 (grade 2)
- `e23` = `0b0110` = 6 (grade 2)
- `e01` = `0b1001` = 9 (grade 2, degenerate)
- `e02` = `0b1010` = 10 (grade 2, degenerate)
- `e03` = `0b1100` = 12 (grade 2, degenerate)

The current implementation must ensure field ordering in the TOML matches the canonical blade index ordering when computing products. If the parser doesn't enforce this or the product generator doesn't account for it, sign errors can occur.

### Issue 2: Missing Versor Identification and Sandwich Products

**Versors** are elements that can be expressed as geometric products of vectors. They're fundamental for transformations in GA:

- **Grade-0,2 versors** (rotors): Product of 2 vectors, perform rotations
- **Grade-0,2,4 versors** (motors in PGA): Product of 2 or 4 vectors, perform rigid motions
- **Flectors** (odd versors): Product of odd number of vectors, perform reflections

The current `is_versor_type()` check is insufficient:
```rust
fn is_versor_type(&self, ty: &TypeSpec) -> bool {
    ty.grades.contains(&0) && ty.grades.iter().all(|g| g % 2 == 0) && ty.name.contains("Rotor")
}
```

This misses:
- Motors (grades 0,2,4) without "Rotor" in name
- Flectors (odd grade versors)
- Any versor type not named "Rotor"

Sandwich products (`v * x * rev(v)`) are the primary mechanism for applying transformations, and should be generated for all versor types.

### Issue 3: Expression Explosion - No Simplification

The generated code contains many redundant terms that could be simplified:

1. **Zero terms**: Many products produce zero terms due to metric (e.g., `e0 * e0 = 0` in PGA)
2. **Constraint-implied simplifications**: Type constraints like the Study condition eliminate degrees of freedom
3. **Symbolic cancellation**: Terms that algebraically cancel aren't detected

Example of bloated output:
```rust
// Generated Motor * Motor (8 fields = 64 initial terms per output field)
Motor::new_unchecked(
    a.s() * b.s() - a.e12() * b.e12() - a.e13() * b.e13() - a.e23() * b.e23(),
    // ... 7 more fields with 6-8 terms each
)
```

With Study condition (`s*e0123 - e23*e01 - e31*e02 - e12*e03 = 0`), some terms could be eliminated or factored.

## Solution

### Phase 1: Blade Ordering Audit and Fix

#### 1.1 Parser Validation

Add validation in `spec/parser.rs` to ensure field ordering matches canonical blade ordering:

```rust
fn validate_field_ordering(ty: &TypeSpec) -> Result<(), ParseError> {
    let mut prev_index = 0;
    for field in &ty.fields {
        if field.blade_index < prev_index {
            return Err(ParseError::NonCanonicalFieldOrdering {
                type_name: ty.name.clone(),
                field: field.name.clone(),
                expected_after: prev_index,
                actual: field.blade_index,
            });
        }
        prev_index = field.blade_index;
    }
    Ok(())
}
```

#### 1.2 Product Generator Fix

Ensure `compute_terms()` uses blade indices consistently:

```rust
fn compute_terms(
    &self,
    type_a: &TypeSpec,
    type_b: &TypeSpec,
    result_blade: usize,  // This is a blade INDEX, not field position
    kind: ProductKind,
) -> Vec<ProductTerm> {
    // Use field.blade_index (canonical) not field position
    for field_a in &type_a.fields {
        let a_blade = field_a.blade_index;  // Canonical index
        // ...
    }
}
```

#### 1.3 Diagnostic Tool

Add a `verify` subcommand to check existing generated code against fresh computation:

```bash
cargo run --package clifford-codegen -- verify src/generated/projective3/
```

### Phase 2: Versor Identification and Sandwich Products

#### 2.1 Versor Identification Algorithm

A type is a versor if it can be expressed as a geometric product of grade-1 elements (vectors):

```rust
/// Versor parity (even vs odd number of vector factors).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum VersorParity {
    /// Even versor: product of even number of vectors (rotor, motor).
    /// Contains only even grades.
    Even,
    /// Odd versor: product of odd number of vectors (flector, reflector).
    /// Contains only odd grades.
    Odd,
}

/// Determines if a type is a versor and its parity.
///
/// A type is a versor if:
/// 1. All its grades have the same parity (all even or all odd)
/// 2. For even versors: grades ⊆ {0, 2, 4, ...}
/// 3. For odd versors: grades ⊆ {1, 3, 5, ...}
///
/// Additionally, the type should satisfy the versor constraint:
/// `v * reverse(v)` is a scalar (for even versors) or pseudoscalar (for odd versors).
pub fn identify_versor(ty: &TypeSpec, algebra: &Algebra) -> Option<VersorParity> {
    if ty.grades.is_empty() {
        return None;
    }

    let first_parity = ty.grades[0] % 2;
    if !ty.grades.iter().all(|g| g % 2 == first_parity) {
        return None;  // Mixed parity, not a versor
    }

    // Verify the versor constraint: v * rev(v) produces scalar/pseudoscalar only
    // This requires symbolic computation or specific checks per algebra

    Some(if first_parity == 0 {
        VersorParity::Even
    } else {
        VersorParity::Odd
    })
}
```

#### 2.2 TOML Specification Extension

Add explicit versor annotation in the spec:

```toml
[types.Motor]
grades = [0, 2, 4]
fields = ["s", "e12", "e13", "e23", "e01", "e02", "e03", "e0123"]
versor = true  # Explicit versor marking

# Constraints that versors satisfy
constraints = [
    { name = "study", expression = "s*e0123 - e23*e01 - e31*e02 - e12*e03 = 0", solve_for = "e0123" }
]

# Sandwich product targets - what types can this versor transform
[types.Motor.transforms]
Point = "Point"       # Motor * Point * rev(Motor) -> Point
Line = "Line"         # Motor * Line * rev(Motor) -> Line
Plane = "Plane"       # Motor * Plane * rev(Motor) -> Plane
```

#### 2.3 Sandwich Product Generation

Generate sandwich products for all versor-operand pairs:

```rust
/// Generates sandwich product functions for all versors.
fn generate_all_sandwich(&self) -> TokenStream {
    let mut products = Vec::new();

    for versor_type in self.spec.types.iter().filter(|t| self.is_versor(t)) {
        // Get transformation targets from spec or infer from grades
        let targets = self.get_sandwich_targets(versor_type);

        for operand_type in targets {
            if let Some(product) = self.generate_sandwich_product(versor_type, &operand_type) {
                products.push(product);
            }
        }
    }

    quote! { #(#products)* }
}

/// Computes sandwich product: v * x * rev(v).
/// Uses expanded formula, not naive three-way multiplication.
fn compute_sandwich_terms(
    &self,
    versor_type: &TypeSpec,
    operand_type: &TypeSpec,
    result_blade: usize,
) -> Vec<SandwichTerm> {
    // ... existing implementation with proper versor handling ...
}
```

### Phase 3: Symbolica Integration for Expression Simplification

#### 3.1 Symbolic Expression Building

Extend `symbolic/product.rs` to build and simplify product expressions:

```rust
use symbolica::atom::Atom;
use symbolica::evaluate::EvaluateOptions;

/// Builds a symbolic expression for a product output field.
pub fn build_symbolic_expression(
    &self,
    type_a: &TypeSpec,
    type_b: &TypeSpec,
    output_field: &FieldSpec,
    kind: ProductKind,
    constraints: &[ConstraintExpr],
) -> Atom {
    let a_symbols = self.create_field_symbols(type_a, "a");
    let b_symbols = self.create_field_symbols(type_b, "b");

    // Build raw expression from product table
    let raw_expr = self.compute_field(type_a, type_b, output_field.blade_index, kind, &a_symbols, &b_symbols);

    // Apply constraint substitutions
    let constrained_expr = self.apply_constraints(raw_expr, constraints);

    // Simplify using Symbolica
    self.simplify(constrained_expr)
}

/// Applies constraint-based substitutions.
///
/// For example, with Study condition `s*e0123 = e23*e01 + e31*e02 + e12*e03`,
/// we can substitute `e0123 = (e23*e01 + e31*e02 + e12*e03) / s` when `s != 0`.
fn apply_constraints(&self, expr: Atom, constraints: &[ConstraintExpr]) -> Atom {
    let mut result = expr;

    for constraint in constraints {
        if let Some(solve_for) = &constraint.solve_for {
            // Get the solved expression for the constrained field
            let solution = self.solver.solve(constraint, solve_for)?;
            // Substitute in the expression
            result = result.substitute(solve_for, &solution);
        }
    }

    result
}

/// Simplifies an expression using Symbolica.
fn simplify(&self, expr: Atom) -> Atom {
    // Expand, collect, and simplify
    expr.expand()
        .collect_terms()
        .factor_common()
}
```

#### 3.2 Expression to Rust Code

Convert simplified Symbolica expressions to Rust TokenStream:

```rust
/// Converts a Symbolica Atom to a Rust TokenStream.
fn atom_to_rust(&self, expr: &Atom, field_prefixes: &HashMap<String, &str>) -> TokenStream {
    match expr {
        Atom::Num(n) => {
            if n.is_zero() {
                quote! { T::zero() }
            } else if n.is_one() {
                quote! { T::one() }
            } else {
                let val = n.to_f64();
                quote! { T::from_f64(#val) }
            }
        }
        Atom::Var(name) => {
            // Parse field reference: "a_x" -> a.x()
            let (prefix, field) = parse_field_name(name);
            let prefix_ident = format_ident!("{}", prefix);
            let field_ident = format_ident!("{}", field);
            quote! { #prefix_ident.#field_ident() }
        }
        Atom::Mul(factors) => {
            let factor_exprs: Vec<_> = factors.iter()
                .map(|f| self.atom_to_rust(f, field_prefixes))
                .collect();
            quote! { #(#factor_exprs) *& * }
        }
        Atom::Add(terms) => {
            // Handle leading minus signs properly
            let mut result = Vec::new();
            for (i, term) in terms.iter().enumerate() {
                let term_expr = self.atom_to_rust(term, field_prefixes);
                if i == 0 {
                    result.push(term_expr);
                } else if is_negative(term) {
                    let positive = negate(term);
                    let positive_expr = self.atom_to_rust(&positive, field_prefixes);
                    result.push(quote! { - #positive_expr });
                } else {
                    result.push(quote! { + #term_expr });
                }
            }
            quote! { #(#result)* }
        }
        Atom::Neg(inner) => {
            let inner_expr = self.atom_to_rust(inner, field_prefixes);
            quote! { -(#inner_expr) }
        }
        // ... other cases
    }
}
```

#### 3.3 Constraint-Aware Generation Pipeline

```rust
impl<'a> ProductGenerator<'a> {
    /// Generates a product function with symbolic simplification.
    fn generate_geometric_from_entry_simplified(&self, entry: &ProductEntry) -> Option<TokenStream> {
        let type_a = self.find_type(&entry.lhs)?;
        let type_b = self.find_type(&entry.rhs)?;
        let output_type = self.find_type(&entry.output)?;

        // Collect all relevant constraints
        let constraints = self.collect_constraints(type_a, type_b, output_type);

        // Build symbolic expressions for each output field
        let symbolic = SymbolicProduct::new(self.algebra);
        let field_exprs: Vec<TokenStream> = output_type.fields.iter().map(|field| {
            let expr = symbolic.build_symbolic_expression(
                type_a, type_b, field,
                ProductKind::Geometric,
                &constraints,
            );
            self.atom_to_rust(&expr, &[("a", type_a), ("b", type_b)])
        }).collect();

        // Generate function with simplified expressions
        // ...
    }
}
```

## Deliverables

### Phase 1: Blade Ordering
- [ ] Add `validate_field_ordering()` in parser
- [ ] Add warning/error for non-canonical field ordering
- [ ] Add `verify` CLI subcommand
- [ ] Audit existing algebra specs for ordering issues
- [ ] Fix any discovered ordering bugs

### Phase 2: Versors and Sandwich Products
- [ ] Add `VersorParity` enum
- [ ] Implement `identify_versor()` function
- [ ] Add `versor` field to TOML spec format
- [ ] Add `transforms` section to TOML spec
- [ ] Generate sandwich products for all versor-target pairs
- [ ] Add sandwich products to spec inference (`products.sandwich` section)
- [ ] Document versor identification criteria

### Phase 3: Expression Simplification
- [ ] Extend `SymbolicProduct` with constraint-aware building
- [ ] Implement `apply_constraints()` for substitution
- [ ] Implement `simplify()` using Symbolica
- [ ] Add `atom_to_rust()` converter
- [ ] Update `ProductGenerator` to use simplified pipeline
- [ ] Add tests comparing simplified vs. unsimplified output
- [ ] Benchmark codegen time with/without simplification

## TOML Specification Changes

### New Fields

```toml
# Type-level versor annotation
[types.Motor]
versor = true                    # Explicitly mark as versor
versor_parity = "even"           # Optional: "even" or "odd" (inferred if omitted)

# Sandwich product targets
[types.Motor.transforms]
Point = "Point"
Line = "Line"
Plane = "Plane"
Bivector = "Bivector"

# Product output section (extended)
[products.sandwich]
Motor_Point = "Point"
Motor_Line = "Line"
Flector_Point = "Point"
```

### Constraint Enhancement

```toml
[types.Motor.constraints.study]
expression = "s*e0123 - e23*e01 - e31*e02 - e12*e03 = 0"
solve_for = "e0123"
# New: enable substitution in product expressions
substitute_in_products = true
```

## Testing Strategy

### Blade Ordering Tests

```rust
#[test]
fn field_ordering_matches_blade_indices() {
    let spec = parse_spec(include_str!("../../algebras/projective3.toml")).unwrap();

    for ty in &spec.types {
        let mut prev_index = 0;
        for field in &ty.fields {
            assert!(field.blade_index >= prev_index,
                "Type {} field {} has non-canonical ordering: {} < {}",
                ty.name, field.name, field.blade_index, prev_index);
            prev_index = field.blade_index;
        }
    }
}
```

### Versor Identification Tests

```rust
proptest! {
    #[test]
    fn versor_reverse_product_is_scalar(v in any::<Motor<f64>>()) {
        // For a versor v, v * rev(v) should be a scalar (grade 0)
        let result = geometric_motor_motor(&v, &v.reverse());
        // All non-scalar components should be approximately zero
        prop_assert!(abs_diff_eq!(result.e12(), 0.0, epsilon = 1e-10));
        prop_assert!(abs_diff_eq!(result.e13(), 0.0, epsilon = 1e-10));
        // ...
    }
}
```

### Expression Simplification Tests

```rust
#[test]
fn simplified_product_equals_unsimplified() {
    let spec = parse_spec(include_str!("../../algebras/projective3.toml")).unwrap();
    let algebra = Algebra::projective(3);

    // Generate both ways
    let generator_simple = ProductGenerator::new(&spec, &algebra, true);
    let generator_full = ProductGenerator::new(&spec, &algebra, false);

    // For random inputs, results should match
    proptest!(|(a: Motor<f64>, b: Motor<f64>)| {
        let simple_result = eval_generated(generator_simple.code(), a, b);
        let full_result = eval_generated(generator_full.code(), a, b);
        prop_assert!(abs_diff_eq!(simple_result, full_result, epsilon = 1e-10));
    });
}
```

## Success Criteria

1. **Blade Ordering**: All algebra specs pass validation; verify command finds no discrepancies
2. **Versors**: Motor, Rotor, Flector correctly identified; sandwich products generated
3. **Simplification**:
   - Product expressions have fewer terms than naive expansion
   - Generated code compiles and passes all tests
   - No performance regression in generated code

## Dependencies

- Symbolica crate (already in dependencies)
- Existing `symbolic/` module infrastructure

## Files Changed

| File | Action |
|------|--------|
| `crates/clifford-codegen/src/spec/parser.rs` | Update - Add field ordering validation |
| `crates/clifford-codegen/src/spec/raw.rs` | Update - Add versor/transforms fields |
| `crates/clifford-codegen/src/spec/ir.rs` | Update - Add VersorInfo to TypeSpec |
| `crates/clifford-codegen/src/algebra/versor.rs` | Create - Versor identification |
| `crates/clifford-codegen/src/codegen/products.rs` | Update - Simplification pipeline |
| `crates/clifford-codegen/src/symbolic/simplify.rs` | Create - Expression simplification |
| `crates/clifford-codegen/src/symbolic/to_rust.rs` | Create - Atom to TokenStream |
| `crates/clifford-codegen/src/main.rs` | Update - Add verify subcommand |
| `algebras/projective3.toml` | Update - Add versor annotations |
| `algebras/euclidean3.toml` | Update - Add versor annotations |

## References

- [Versors in Geometric Algebra (Wikipedia)](https://en.wikipedia.org/wiki/Versor#Versors_in_geometric_algebra)
- [Symbolica documentation](https://symbolica.io/)
- [Rigid Body Dynamics with Clifford Algebra](https://arxiv.org/abs/2201.13282) - Motor composition formulas
