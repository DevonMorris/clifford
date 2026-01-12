# PRD-16.3: Expression Simplification with Symbolica

**Status**: Draft
**Parent**: PRD-16
**Goal**: Use Symbolica to simplify generated product expressions before code emission

## Problem Statement

### Current State

The code generator produces product expressions by accumulating all non-zero terms from the product table. For complex products like Motor × Motor in PGA, this leads to expression explosion:

```rust
// Current generated code (simplified example)
Motor::new_unchecked(
    a.s() * b.s() - a.e12() * b.e12() - a.e13() * b.e13() - a.e23() * b.e23()
        + a.e01() * b.e01() * ZERO - a.e02() * b.e02() * ZERO - ...  // Many zero terms!
        + a.e0123() * b.e0123() * ZERO * ZERO,  // More redundant zeros
    // ... 7 more output fields with similar bloat
)
```

### Issues

1. **Zero terms**: Many terms produce zero due to degenerate metric (e₀² = 0 in PGA)
2. **Constraint redundancy**: Type constraints like Study condition could eliminate computed fields
3. **No algebraic simplification**: Terms that cancel aren't detected
4. **Duplicated computation**: Common subexpressions aren't factored out

### Example: Motor × Motor

For PGA Motor (8 fields) × Motor (8 fields):
- **Naive**: 8 × 8 = 64 input pairs per output field
- **After metric**: ~40 non-zero terms per field (degenerate e₀ eliminates some)
- **After simplification**: ~12-20 terms per field (algebraic cancellation)
- **With constraint substitution**: ~8-12 terms (Study condition eliminates e₀₁₂₃ computation)

### The Symbolic Module

The codebase already has a `symbolic/` module with:
- `parser.rs`: Parses constraint strings into Symbolica `Atom`s
- `product.rs`: Computes symbolic products (partially implemented)
- `solver.rs`: Solves constraints for specified variables
- `verify.rs`: Verifies constraints are satisfied

However, **simplification is not yet integrated** into the code generation pipeline.

## Solution

### Phase 1: Symbolic Expression Building

#### 1.1 Extend SymbolicProduct

Complete the `SymbolicProduct` implementation to build full expressions:

```rust
use symbolica::atom::{Atom, AtomCore};

/// Computes symbolic product expressions.
pub struct SymbolicProduct<'a> {
    algebra: &'a Algebra,
    /// Cache of field symbols for each type.
    type_symbols: HashMap<String, Vec<(String, Atom)>>,
}

impl<'a> SymbolicProduct<'a> {
    /// Creates symbols for a type's fields with a prefix.
    ///
    /// For Motor with prefix "a", creates: a_s, a_e12, a_e13, ...
    pub fn create_field_symbols(&mut self, ty: &TypeSpec, prefix: &str) -> Vec<(String, Atom)> {
        let key = format!("{}_{}", prefix, ty.name);

        if let Some(cached) = self.type_symbols.get(&key) {
            return cached.clone();
        }

        let symbols: Vec<(String, Atom)> = ty.fields.iter().map(|field| {
            let name = format!("{}_{}", prefix, field.name);
            let atom = Atom::parse(&name).unwrap();
            (field.name.clone(), atom)
        }).collect();

        self.type_symbols.insert(key, symbols.clone());
        symbols
    }

    /// Builds the symbolic expression for a product output field.
    pub fn build_field_expression(
        &mut self,
        type_a: &TypeSpec,
        type_b: &TypeSpec,
        output_blade: usize,
        kind: ProductKind,
    ) -> Atom {
        let a_symbols = self.create_field_symbols(type_a, "a");
        let b_symbols = self.create_field_symbols(type_b, "b");

        let mut result = Atom::num(0);

        // Accumulate all terms that contribute to output_blade
        for (a_field, a_sym) in &a_symbols {
            let a_blade = type_a.field_blade_index(a_field);

            for (b_field, b_sym) in &b_symbols {
                let b_blade = type_b.field_blade_index(b_field);

                // Check if this pair contributes to output_blade
                if let Some(sign) = self.compute_contribution(a_blade, b_blade, output_blade, kind) {
                    if sign != 0 {
                        let term = a_sym * b_sym;
                        let signed_term = if sign > 0 { term } else { -term };
                        result = result + signed_term;
                    }
                }
            }
        }

        result
    }

    /// Computes whether (a_blade, b_blade) contributes to output_blade.
    /// Returns the sign (+1, -1) or None if no contribution.
    fn compute_contribution(
        &self,
        a_blade: usize,
        b_blade: usize,
        output_blade: usize,
        kind: ProductKind,
    ) -> Option<i32> {
        match kind {
            ProductKind::Geometric => {
                let result_blade = a_blade ^ b_blade;
                if result_blade == output_blade {
                    Some(self.algebra.product_sign(a_blade, b_blade))
                } else {
                    None
                }
            }
            ProductKind::Outer => {
                // Outer: only non-zero if blades don't share basis vectors
                if a_blade & b_blade == 0 {
                    let result_blade = a_blade | b_blade;
                    if result_blade == output_blade {
                        Some(self.algebra.product_sign(a_blade, b_blade))
                    } else {
                        None
                    }
                } else {
                    None
                }
            }
            ProductKind::LeftContraction => {
                // Left contraction: grade(result) = grade(b) - grade(a)
                // ... similar logic ...
                todo!()
            }
            // ... other product kinds ...
        }
    }
}
```

### Phase 2: Expression Simplification

#### 2.1 Basic Simplification

Use Symbolica's built-in simplification:

```rust
impl<'a> SymbolicProduct<'a> {
    /// Simplifies an expression using Symbolica.
    pub fn simplify(&self, expr: Atom) -> Atom {
        // Expand all products and collect like terms
        expr.expand().collect()
    }
}
```

#### 2.2 Constraint-Based Simplification

Apply type constraints to eliminate terms:

```rust
/// Applies constraint substitutions to simplify expressions.
pub struct ConstraintSimplifier<'a> {
    /// The constraint solver for computing substitutions.
    solver: ConstraintSolver,
    /// Cached solutions for solve_for constraints.
    solutions: HashMap<String, Atom>,
    /// The algebra spec for constraint access.
    spec: &'a AlgebraSpec,
}

impl<'a> ConstraintSimplifier<'a> {
    /// Pre-computes all solve_for substitutions.
    pub fn new(spec: &'a AlgebraSpec) -> Self {
        let solver = ConstraintSolver::new();
        let mut solutions = HashMap::new();

        // For each type with constraints, solve for the solve_for field
        for ty in &spec.types {
            for constraint in &ty.constraints {
                if let Some(ref solve_for) = constraint.solve_for {
                    // Parse constraint and solve
                    if let Ok(result) = solver.solve(&constraint.expression, solve_for) {
                        // Convert solution to Symbolica Atom
                        let solution_atom = self.parse_solution(&result);
                        let var_name = format!("{}_{}", ty.name.to_lowercase(), solve_for);
                        solutions.insert(var_name, solution_atom);
                    }
                }
            }
        }

        Self { solver, solutions, spec }
    }

    /// Applies all known substitutions to an expression.
    pub fn apply(&self, expr: &Atom) -> Atom {
        let mut result = expr.clone();

        for (var_name, solution) in &self.solutions {
            // Create the variable to substitute
            let var = Atom::parse(var_name).unwrap();
            // Substitute var -> solution
            result = result.substitute(&var, solution);
        }

        // Simplify after substitution
        result.expand().collect()
    }

    /// Parses a SolveResult into a Symbolica Atom.
    fn parse_solution(&self, result: &SolveResult) -> Atom {
        let parser = ConstraintParser::new();

        match result.solution_type {
            SolutionType::Linear => {
                let numerator = Atom::parse(&result.numerator).unwrap_or(Atom::num(0));
                if let Some(ref divisor) = result.divisor {
                    let divisor_atom = Atom::parse(divisor).unwrap_or(Atom::num(1));
                    numerator / divisor_atom
                } else {
                    numerator
                }
            }
            SolutionType::Quadratic => {
                // For quadratic, we can't substitute directly (would need sqrt)
                // Instead, skip this substitution for code generation
                // The constraint is still used for validation
                Atom::parse(&format!("sqrt({})", result.numerator))
                    .unwrap_or(Atom::num(0))
            }
        }
    }
}
```

#### 2.3 Zero Elimination

Eliminate terms that are zero due to degenerate metric:

```rust
impl<'a> SymbolicProduct<'a> {
    /// Eliminates terms that are zero due to degenerate metric.
    ///
    /// In PGA, e₀² = 0, so any term with e₀*e₀ (or e₀₁*e₀₁, etc.) is zero.
    pub fn eliminate_zeros(&self, expr: Atom) -> Atom {
        // For each degenerate basis vector, e_i² = 0
        // This means any blade containing e_i twice squares to zero

        let mut result = expr;

        // Find all degenerate basis vectors
        for (i, metric) in self.algebra.basis_metrics().iter().enumerate() {
            if *metric == 0 {
                // e_i is degenerate
                // Create substitution rule: e_i * e_i -> 0
                let basis_name = self.algebra.basis_name(i);
                let var = Atom::parse(&format!("{}", basis_name)).unwrap();
                let squared = &var * &var;

                // This is a simplification hint - Symbolica should handle this
                // via the metric algebra rules
            }
        }

        result.expand()
    }
}
```

### Phase 3: Rust Code Generation

#### 3.1 Atom to TokenStream Conversion

Convert simplified Symbolica expressions to Rust code:

```rust
use proc_macro2::TokenStream;
use quote::quote;

/// Converts Symbolica atoms to Rust TokenStream.
pub struct AtomToRust<'a> {
    /// Map from symbol name to (prefix, field) for accessor generation.
    /// e.g., "a_s" -> ("a", "s") -> `a.s()`
    symbol_map: HashMap<String, (String, String)>,
    /// The types being operated on.
    types: &'a [&'a TypeSpec],
}

impl<'a> AtomToRust<'a> {
    /// Creates a new converter.
    pub fn new(types: &'a [&'a TypeSpec], prefixes: &[&str]) -> Self {
        let mut symbol_map = HashMap::new();

        for (ty, prefix) in types.iter().zip(prefixes.iter()) {
            for field in &ty.fields {
                let symbol_name = format!("{}_{}", prefix, field.name);
                symbol_map.insert(symbol_name, (prefix.to_string(), field.name.clone()));
            }
        }

        Self { symbol_map, types }
    }

    /// Converts a Symbolica Atom to Rust TokenStream.
    pub fn convert(&self, atom: &Atom) -> TokenStream {
        // Match on the atom type
        // Note: Symbolica's API may vary - this is conceptual

        match self.classify_atom(atom) {
            AtomKind::Number(n) => self.convert_number(n),
            AtomKind::Variable(name) => self.convert_variable(&name),
            AtomKind::Sum(terms) => self.convert_sum(&terms),
            AtomKind::Product(factors) => self.convert_product(&factors),
            AtomKind::Negation(inner) => self.convert_negation(&inner),
            AtomKind::Power(base, exp) => self.convert_power(&base, exp),
        }
    }

    /// Converts a number to Rust.
    fn convert_number(&self, n: i64) -> TokenStream {
        match n {
            0 => quote! { T::zero() },
            1 => quote! { T::one() },
            2 => quote! { T::TWO },
            -1 => quote! { -T::one() },
            -2 => quote! { -T::TWO },
            _ => {
                if n >= i8::MIN as i64 && n <= i8::MAX as i64 {
                    let n_i8 = n as i8;
                    quote! { T::from_i8(#n_i8) }
                } else {
                    let n_f64 = n as f64;
                    quote! { T::from_f64(#n_f64) }
                }
            }
        }
    }

    /// Converts a variable reference to Rust.
    fn convert_variable(&self, name: &str) -> TokenStream {
        if let Some((prefix, field)) = self.symbol_map.get(name) {
            let prefix_ident = format_ident!("{}", prefix);
            let field_ident = format_ident!("{}", field);
            quote! { #prefix_ident.#field_ident() }
        } else {
            // Unknown variable - emit as identifier (for constants)
            let ident = format_ident!("{}", name);
            quote! { #ident }
        }
    }

    /// Converts a sum to Rust.
    fn convert_sum(&self, terms: &[Atom]) -> TokenStream {
        if terms.is_empty() {
            return quote! { T::zero() };
        }

        // Convert first term
        let first = self.convert(&terms[0]);

        // Add remaining terms with proper signs
        let mut result = first;
        for term in &terms[1..] {
            let (is_negative, positive_term) = self.extract_sign(term);
            let term_expr = self.convert(&positive_term);

            if is_negative {
                result = quote! { #result - #term_expr };
            } else {
                result = quote! { #result + #term_expr };
            }
        }

        result
    }

    /// Converts a product to Rust.
    fn convert_product(&self, factors: &[Atom]) -> TokenStream {
        if factors.is_empty() {
            return quote! { T::one() };
        }

        let exprs: Vec<TokenStream> = factors.iter().map(|f| self.convert(f)).collect();

        // Chain multiplications
        let first = &exprs[0];
        let rest = &exprs[1..];

        if rest.is_empty() {
            first.clone()
        } else {
            quote! { #first #(* #rest)* }
        }
    }

    /// Extracts sign from a term.
    fn extract_sign(&self, atom: &Atom) -> (bool, Atom) {
        // Check if atom is a negation or has negative leading coefficient
        // Return (is_negative, positive_form)
        // Implementation depends on Symbolica's API
        (false, atom.clone())  // Placeholder
    }
}
```

### Phase 4: Integration with ProductGenerator

#### 4.1 Update ProductGenerator

Integrate simplification into the generation pipeline:

```rust
impl<'a> ProductGenerator<'a> {
    /// Generates a geometric product with simplification.
    pub fn generate_geometric_simplified(&self, entry: &ProductEntry) -> Option<TokenStream> {
        let type_a = self.find_type(&entry.lhs)?;
        let type_b = self.find_type(&entry.rhs)?;
        let output_type = self.find_type(&entry.output)?;

        // Build symbolic expressions for each output field
        let mut symbolic = SymbolicProduct::new(self.algebra);
        let simplifier = ConstraintSimplifier::new(self.spec);

        let field_exprs: Vec<TokenStream> = output_type.fields.iter().map(|field| {
            // Build raw symbolic expression
            let raw_expr = symbolic.build_field_expression(
                type_a,
                type_b,
                field.blade_index,
                ProductKind::Geometric,
            );

            // Apply constraint substitutions
            let constrained = simplifier.apply(&raw_expr);

            // Simplify
            let simplified = symbolic.simplify(constrained);

            // Convert to Rust
            let converter = AtomToRust::new(&[type_a, type_b], &["a", "b"]);
            converter.convert(&simplified)
        }).collect();

        // Generate the function
        let fn_name = format_ident!(
            "geometric_{}_{}",
            entry.lhs.to_lowercase(),
            entry.rhs.to_lowercase()
        );
        let a_type = format_ident!("{}", entry.lhs);
        let b_type = format_ident!("{}", entry.rhs);
        let out_type = format_ident!("{}", entry.output);

        let constructor = if output_type.has_constraints() {
            quote! { new_unchecked }
        } else {
            quote! { new }
        };

        Some(quote! {
            /// Computes the geometric product.
            #[inline]
            pub fn #fn_name(a: &#a_type<T>, b: &#b_type<T>) -> #out_type<T> {
                #out_type::#constructor(#(#field_exprs),*)
            }
        })
    }
}
```

### Phase 5: Configuration and Diagnostics

#### 5.1 Generation Options

Add configuration for simplification:

```toml
[options]
# Enable symbolic simplification (default: true)
simplify_expressions = true

# Apply constraint substitutions (default: true)
apply_constraints = true

# Debug: emit unsimplified expressions as comments
debug_expressions = false
```

#### 5.2 Diagnostic Output

Optionally emit comparison of simplified vs unsimplified:

```rust
// Generated with debug_expressions = true
/// Computes the geometric product.
///
/// # Simplification
///
/// Before: a_s * b_s - a_e12 * b_e12 - a_e13 * b_e13 - a_e23 * b_e23
///         + a_e01 * b_e01 * 0 - a_e02 * b_e02 * 0 + ...
/// After:  a_s * b_s - a_e12 * b_e12 - a_e13 * b_e13 - a_e23 * b_e23
/// Reduction: 64 terms -> 4 terms
#[inline]
pub fn geometric_motor_motor(a: &Motor<T>, b: &Motor<T>) -> Motor<T> {
    Motor::new_unchecked(
        a.s() * b.s() - a.e12() * b.e12() - a.e13() * b.e13() - a.e23() * b.e23(),
        // ...
    )
}
```

## Testing Strategy

### Correctness Tests

```rust
proptest! {
    #[test]
    fn simplified_equals_naive(
        a in any::<Motor<f64>>(),
        b in any::<Motor<f64>>()
    ) {
        let simplified = geometric_motor_motor(&a, &b);
        let naive = geometric_motor_motor_naive(&a, &b);

        prop_assert!(abs_diff_eq!(simplified.s(), naive.s(), epsilon = ABS_DIFF_EQ_EPS));
        prop_assert!(abs_diff_eq!(simplified.e12(), naive.e12(), epsilon = ABS_DIFF_EQ_EPS));
        // ... all fields
    }
}
```

### Simplification Tests

```rust
#[test]
fn zero_terms_eliminated() {
    let spec = parse_spec(include_str!("../../algebras/projective3.toml")).unwrap();
    let algebra = Algebra::projective(3, 0, 1);

    let mut symbolic = SymbolicProduct::new(&algebra);
    let motor = spec.types.iter().find(|t| t.name == "Motor").unwrap();

    // Build expression for s field of Motor * Motor
    let expr = symbolic.build_field_expression(motor, motor, 0, ProductKind::Geometric);
    let simplified = symbolic.simplify(expr);

    // Count terms - should be much less than 64
    let term_count = count_terms(&simplified);
    assert!(term_count <= 10, "Expected <= 10 terms, got {}", term_count);
}

#[test]
fn constraint_reduces_terms() {
    let spec = parse_spec(include_str!("../../algebras/projective3.toml")).unwrap();
    let simplifier = ConstraintSimplifier::new(&spec);

    // Expression that includes e0123
    let expr = Atom::parse("a_s * b_e0123 + a_e0123 * b_s").unwrap();

    // After applying Study condition, e0123 is replaced
    let simplified = simplifier.apply(&expr);

    // e0123 should not appear in the result
    let result_str = simplified.to_string();
    assert!(!result_str.contains("e0123"), "e0123 should be eliminated");
}
```

### Benchmarks

Compare codegen time and generated code quality:

```rust
#[bench]
fn bench_codegen_with_simplification(b: &mut Bencher) {
    let spec = parse_spec(include_str!("../../algebras/projective3.toml")).unwrap();

    b.iter(|| {
        let generator = ProductGenerator::new(&spec, true); // simplification on
        generator.generate_all()
    });
}

#[bench]
fn bench_codegen_without_simplification(b: &mut Bencher) {
    let spec = parse_spec(include_str!("../../algebras/projective3.toml")).unwrap();

    b.iter(|| {
        let generator = ProductGenerator::new(&spec, false); // simplification off
        generator.generate_all()
    });
}
```

## Expected Impact

### Term Reduction

| Product | Naive Terms | Simplified | Reduction |
|---------|-------------|------------|-----------|
| Motor × Motor (scalar field) | 64 | ~8 | 88% |
| Motor × Point | 32 | ~12 | 63% |
| Motor × Line | 48 | ~20 | 58% |
| Rotor × Vector | 8 | ~4 | 50% |

### Code Quality

Before:
```rust
a.s() * b.s() - a.e12() * b.e12() - a.e13() * b.e13() - a.e23() * b.e23()
    + a.e01() * b.e01() * T::zero() - a.e02() * b.e02() * T::zero()
    - a.e03() * b.e03() * T::zero() + a.e0123() * b.e0123() * T::zero()
```

After:
```rust
a.s() * b.s() - a.e12() * b.e12() - a.e13() * b.e13() - a.e23() * b.e23()
```

## Deliverables

- [ ] Complete `SymbolicProduct::build_field_expression()`
- [ ] Implement `SymbolicProduct::simplify()` using Symbolica
- [ ] Implement `ConstraintSimplifier` for constraint substitution
- [ ] Implement `AtomToRust` for code generation
- [ ] Integrate into `ProductGenerator`
- [ ] Add configuration options in TOML
- [ ] Add diagnostic/debug output
- [ ] Add correctness tests comparing simplified vs naive
- [ ] Add benchmarks

## Success Criteria

1. Generated expressions have fewer terms than naive expansion
2. All tests pass (simplified matches naive numerically)
3. Codegen time is acceptable (< 2x increase)
4. Generated code compiles without issues

## Files Changed

| File | Action | Description |
|------|--------|-------------|
| `crates/clifford-codegen/src/symbolic/product.rs` | Update | Complete symbolic product building |
| `crates/clifford-codegen/src/symbolic/simplify.rs` | Create | Expression simplification |
| `crates/clifford-codegen/src/symbolic/to_rust.rs` | Create | Atom to TokenStream conversion |
| `crates/clifford-codegen/src/symbolic/mod.rs` | Update | Re-export new modules |
| `crates/clifford-codegen/src/codegen/products.rs` | Update | Use simplified pipeline |
| `crates/clifford-codegen/src/spec/raw.rs` | Update | Add simplification options |

## Dependencies

- Symbolica crate (already in dependencies)
- proc-macro2 and quote (already in dependencies)

## References

- [Symbolica documentation](https://symbolica.io/docs/)
- [Computer Algebra in Clifford Algebras](https://arxiv.org/abs/math/0407364)
