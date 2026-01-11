# PRD-14.12: Symbolic Constraint Verification with Symbolica

## Overview

This PRD describes the integration of [Symbolica](https://symbolica.io), a high-performance computer algebra system for Rust, to symbolically verify that product outputs satisfy algebraic constraints.

## Problem Statement

Types in geometric algebra often have constraints expressed over their fields. For example:
- `Rotor` has constraint `s * s + xy * xy + xz * xz + yz * yz = 1` (unit norm)
- `NullPoint` in CGA has constraint `x * x + y * y - ep * em = 0` (null vector)

When products are computed (e.g., `Rotor * Rotor`), we need to verify that the output still satisfies the constraint. Currently, this verification happens at runtime (or is simply assumed). Symbolic verification allows us to prove at compile time that constraints are preserved.

## Goals

1. **Symbolic constraint parsing**: Parse constraint expressions from the TOML spec into Symbolica expressions
2. **Product output symbolization**: Represent product outputs as symbolic expressions
3. **Constraint verification**: Prove that if inputs satisfy constraints, outputs also satisfy constraints
4. **Integration with codegen**: Use verification results to generate appropriate code (panic/assert vs unchecked)

## Non-Goals

- Full computer algebra system integration (we use Symbolica specifically for constraint checking)
- Arbitrary constraint solving (we only verify specific constraint forms)
- Runtime symbolic computation (all symbolic work happens at codegen time)

## Design

### 1. Dependency Addition

Add Symbolica to `Cargo.toml`:

```toml
[dependencies]
symbolica = "1.0"
```

### 2. Constraint Expression Parser

Parse constraint strings into Symbolica atoms:

```rust
use symbolica::atom::Atom;
use symbolica::state::State;

pub struct ConstraintParser {
    state: State,
}

impl ConstraintParser {
    /// Parses a constraint expression like "s * s + xy * xy = 1"
    pub fn parse(&self, expr: &str, field_names: &[&str]) -> Result<ConstraintExpr, Error> {
        // Split on "=" to get lhs and rhs
        let (lhs, rhs) = parse_equality(expr)?;

        // Parse each side as a Symbolica expression
        let lhs_atom = self.parse_expression(lhs, field_names)?;
        let rhs_atom = self.parse_expression(rhs, field_names)?;

        Ok(ConstraintExpr { lhs: lhs_atom, rhs: rhs_atom })
    }

    /// Parses a mathematical expression into a Symbolica Atom
    fn parse_expression(&self, expr: &str, field_names: &[&str]) -> Result<Atom, Error> {
        // Replace field names with Symbolica variables
        // Handle operators: +, -, *, /, ^
        // Handle numeric literals
        todo!()
    }
}
```

### 3. Product Symbolic Evaluation

Given the symbolic expressions for product inputs and the product formula, compute the symbolic output:

```rust
pub struct SymbolicProduct {
    /// The algebra for blade computations
    algebra: Algebra,
    /// Product table
    table: ProductTable,
}

impl SymbolicProduct {
    /// Computes the symbolic output of a product
    pub fn compute(&self, type_a: &TypeSpec, type_b: &TypeSpec, product_kind: ProductKind)
        -> Vec<SymbolicField>
    {
        // For each output field, compute the symbolic expression
        // by applying the product formula to symbolic input fields
        todo!()
    }
}
```

### 4. Constraint Verification

Given symbolic input constraints and symbolic product output, verify the output satisfies the constraint:

```rust
pub struct ConstraintVerifier;

impl ConstraintVerifier {
    /// Verifies that if inputs satisfy their constraints,
    /// the output satisfies its constraint.
    ///
    /// Returns:
    /// - `Verified::Always` - Constraint always holds
    /// - `Verified::Conditional(cond)` - Constraint holds if condition is met
    /// - `Verified::Never(counterexample)` - Constraint can be violated
    pub fn verify(
        &self,
        input_constraints: &[ConstraintExpr],
        output_constraint: &ConstraintExpr,
        output_fields: &[SymbolicField],
    ) -> VerificationResult {
        // Substitute output fields into constraint
        // Simplify under assumption that input constraints hold
        // Check if result simplifies to true
        todo!()
    }
}

pub enum VerificationResult {
    /// Constraint is always satisfied when inputs satisfy their constraints
    Always,
    /// Constraint is satisfied under additional condition
    Conditional(Atom),
    /// Constraint can be violated (with counterexample)
    Never(String),
}
```

### 5. Example Verification

For Rotor * Rotor in Euclidean 3D:

**Input constraint**: `a.s^2 + a.xy^2 + a.xz^2 + a.yz^2 = 1` (for both rotors)

**Product formula** (simplified):
```
c.s  = a.s*b.s - a.xy*b.xy - a.xz*b.xz - a.yz*b.yz
c.xy = a.s*b.xy + a.xy*b.s + a.xz*b.yz - a.yz*b.xz
c.xz = a.s*b.xz + a.xz*b.s - a.xy*b.yz + a.yz*b.xy
c.yz = a.s*b.yz + a.yz*b.s + a.xy*b.xz - a.xz*b.xy
```

**Verification**:
Substitute into output constraint and simplify:
```
c.s^2 + c.xy^2 + c.xz^2 + c.yz^2 = 1
```

This should simplify to `1 = 1` when input constraints are applied.

### 6. Integration with Codegen

Update the product generator to use verification results:

```rust
impl ProductGenerator {
    fn generate_product_body(&self, entry: &ProductEntry) -> TokenStream {
        // Get constraints from input and output types
        let input_constraints = self.get_constraints(&entry.lhs, &entry.rhs);
        let output_constraint = self.get_constraint(&entry.output);

        // Verify constraint preservation
        let result = self.verifier.verify(&input_constraints, &output_constraint, &output_fields);

        match result {
            VerificationResult::Always => {
                // No runtime check needed
                quote! { #output_type::new(#(#field_exprs),*) }
            }
            VerificationResult::Conditional(cond) => {
                // Add debug assertion
                quote! {
                    let result = #output_type::new(#(#field_exprs),*);
                    debug_assert!(result.satisfies_constraint(), "constraint violated: {}", #cond);
                    result
                }
            }
            VerificationResult::Never(reason) => {
                // Compile-time error or always check
                quote! {
                    let result = #output_type::new(#(#field_exprs),*);
                    assert!(result.satisfies_constraint(), "constraint not guaranteed: {}", #reason);
                    result
                }
            }
        }
    }
}
```

## Implementation Plan

### Phase 1: Core Integration
1. Add Symbolica dependency
2. Implement constraint expression parser
3. Add basic symbolic expression support

### Phase 2: Product Symbolization
1. Implement symbolic product computation
2. Generate symbolic expressions for all product formulas
3. Add tests comparing symbolic and numeric results

### Phase 3: Verification
1. Implement constraint verification logic
2. Handle special cases (zero constraints, inequality constraints)
3. Add verification results to ProductEntry

### Phase 4: Codegen Integration
1. Update ProductGenerator to use verification
2. Generate appropriate runtime checks
3. Add verification status to generated code comments

## Testing Strategy

1. **Unit tests**: Verify parser handles common constraint forms
2. **Integration tests**: Verify known products (Rotor*Rotor) pass verification
3. **Property tests**: Random symbolic expressions simplify correctly
4. **Regression tests**: Known constraint violations are detected

## Files to Create/Modify

### New Files
- `crates/clifford-codegen/src/symbolic/mod.rs` - Module root
- `crates/clifford-codegen/src/symbolic/parser.rs` - Constraint parser
- `crates/clifford-codegen/src/symbolic/product.rs` - Symbolic products
- `crates/clifford-codegen/src/symbolic/verify.rs` - Verification logic

### Modified Files
- `crates/clifford-codegen/Cargo.toml` - Add Symbolica dependency
- `crates/clifford-codegen/src/lib.rs` - Export symbolic module
- `crates/clifford-codegen/src/codegen/products.rs` - Integrate verification

## Success Criteria

1. Symbolica successfully parses all constraint expressions in test TOML files
2. Rotor*Rotor verification passes for Euclidean 2D/3D
3. Verification correctly identifies when constraints may be violated
4. Generated code includes appropriate runtime checks based on verification

## Dependencies

- Symbolica 1.0 (https://crates.io/crates/symbolica)

## References

- [Symbolica Documentation](https://docs.rs/symbolica)
- [Symbolica Blog Post](https://symbolica.io/posts/stable_release/)
- [PRD-14.11: Constraint Simplification](prd-14.11-constraint-simplification.md)
