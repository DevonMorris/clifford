# PRD-46: Groebner Basis Constraint Simplification

**Status**: Implemented
**Goal**: Use Symbolica's Groebner basis computation to simplify product expressions using input type constraints

## Problem Statement

### Current State

Generated product expressions contain many redundant terms because they don't account for the algebraic constraints of input types. For example:

**PGA Line × Line:**
- Lines (grade-2 bivectors) satisfy the **Plücker condition**: `e01*e23 + e02*e31 + e03*e12 = 0`
- The generated product naively computes all 36 blade pair combinations
- Many output terms could be eliminated by substituting the Plücker relation

**CGA Circle × Circle:**
- Circles (grade-3 trivectors) satisfy the **blade constraint**: `C ∧ C = 0`
- Additionally: `C * C̃ = scalar` (versor constraint)
- Generated products contain terms that vanish due to these constraints

### Example: Line Meet (Antiwedge)

For `line1.antiwedge(line2)` producing a point, the current generated code:

```rust
// Current: All 36 terms computed
Point::new(
    l1.e01() * l2.e23() + l1.e02() * l2.e31() + l1.e03() * l2.e12()  // direction1 · moment2
    + l1.e23() * l2.e01() + l1.e31() * l2.e02() + l1.e12() * l2.e03() // moment1 · direction1
    + ... // Many more terms
)
```

With Plücker constraint substitution:
```rust
// After: Redundant terms eliminated using constraint
Point::new(
    // Simplified expression with fewer terms
    ...
)
```

### The Constraint Gap

The constraint derivation system (PRD-14.8) already computes constraint expressions:
- `derive_geometric_constraint()` produces `u * ũ = scalar` constraint expressions
- `derive_antiproduct_constraint()` produces `u ⊟ ũ̃ = antiscalar` constraint expressions
- `derive_blade_constraint()` produces `B ∧ B = 0` for simple blades

**But these constraints are not applied during product simplification.** The `ConstraintSimplifier` currently returns empty substitutions.

## Solution: Groebner Basis Reduction

### Core Idea

1. **Collect constraints** as polynomial equations for input types
2. **Compute Groebner basis** of the constraint ideal
3. **Reduce product expressions** modulo the Groebner basis
4. **Generate code** from simplified expressions

### Why Groebner Bases?

Groebner bases provide the mathematically optimal way to simplify polynomials modulo an ideal:

1. **Complete reduction**: Any polynomial can be reduced to a unique canonical form
2. **Automatic**: No manual pattern matching or heuristics needed
3. **Sound**: Result is equivalent to original modulo constraints
4. **Available**: Symbolica provides `GroebnerBasis::new()` and reduction operations

### Mathematical Background

**Polynomial Ideal**: Given constraints `{f₁ = 0, f₂ = 0, ...}`, the ideal is `I = ⟨f₁, f₂, ...⟩`

**Groebner Basis**: A basis `G = {g₁, g₂, ...}` for `I` with special properties:
- Every polynomial in `I` reduces to 0 using `G`
- Reduction is confluent (unique normal form)

**Reduction**: Given expression `p` and basis `G`:
- `p mod G` = unique simplified form
- If `p ∈ I` (satisfies constraints), `p mod G = 0`
- Otherwise, `p mod G` has no further simplification possible

### Example: Plücker Constraint

For PGA lines with fields `{e01, e02, e03, e12, e31, e23}`:

**Constraint**: `e01*e23 + e02*e31 + e03*e12 = 0`

This is already a Groebner basis element (single polynomial).

**Reduction**: Any expression containing `e01*e23` can substitute:
```
e01*e23 → -(e02*e31 + e03*e12)
```

For two lines `a` and `b`:
```
a_e01*a_e23 → -(a_e02*a_e31 + a_e03*a_e12)
b_e01*b_e23 → -(b_e02*b_e31 + b_e03*b_e12)
```

## Implementation

### Phase 1: Constraint Collection

#### 1.1 Extend ConstraintDeriver for Product Inputs

```rust
/// Collects constraint polynomials for input types.
pub struct ProductConstraintCollector<'a> {
    deriver: ConstraintDeriver<'a>,
}

impl<'a> ProductConstraintCollector<'a> {
    /// Collects all constraint polynomials for a type.
    ///
    /// Returns polynomials that must equal zero for valid instances.
    pub fn collect_constraints(
        &self,
        ty: &TypeSpec,
        prefix: &str,  // e.g., "a" for left operand
    ) -> Vec<MultivariatePolynomial<RationalField, u16>> {
        let mut constraints = Vec::new();

        // 1. Geometric constraint: u * ũ = scalar (non-scalar terms = 0)
        if let Some(geo) = self.deriver.derive_geometric_constraint(ty, prefix) {
            for expr in geo.zero_expressions {
                constraints.push(self.atom_to_polynomial(&expr));
            }
        }

        // 2. Antiproduct constraint: u ⊟ ũ̃ = antiscalar
        if let Some(anti) = self.deriver.derive_antiproduct_constraint(ty, prefix) {
            for expr in anti.zero_expressions {
                constraints.push(self.atom_to_polynomial(&expr));
            }
        }

        // 3. Blade constraint: B ∧ B = 0 (for simple blades)
        if let Some(blade) = self.deriver.derive_blade_constraint(ty, prefix) {
            for expr in blade.zero_expressions {
                constraints.push(self.atom_to_polynomial(&expr));
            }
        }

        // 4. Null constraint: v * ṽ = 0 (for null vectors in CGA)
        if let Some(null) = self.deriver.derive_null_constraint(ty, prefix) {
            for expr in null.zero_expressions {
                constraints.push(self.atom_to_polynomial(&expr));
            }
        }

        constraints
    }

    /// Converts Symbolica Atom to MultivariatePolynomial.
    fn atom_to_polynomial(
        &self,
        atom: &Atom,
    ) -> MultivariatePolynomial<RationalField, u16> {
        // Convert symbolic expression to multivariate polynomial
        // over rationals for exact computation
        atom.to_polynomial(&RationalField::new(), None)
            .expect("constraint should be polynomial")
    }
}
```

### Phase 2: Groebner Basis Computation

#### 2.1 GroebnerSimplifier

```rust
use symbolica::poly::groebner::GroebnerBasis;
use symbolica::poly::polynomial::MultivariatePolynomial;
use symbolica::domains::rational::RationalField;

/// Simplifies expressions using Groebner basis reduction.
pub struct GroebnerSimplifier {
    /// The Groebner basis for the constraint ideal.
    basis: Option<GroebnerBasis<RationalField, u16>>,
    /// Variable ordering (affects basis and reduction).
    variables: Vec<String>,
}

impl GroebnerSimplifier {
    /// Creates a simplifier from constraint polynomials.
    ///
    /// # Arguments
    ///
    /// * `constraints` - Polynomials that equal zero (constraint ideal generators)
    /// * `use_grevlex` - Use grevlex ordering (faster) vs lex (better elimination)
    pub fn new(
        constraints: Vec<MultivariatePolynomial<RationalField, u16>>,
        use_grevlex: bool,
    ) -> Self {
        if constraints.is_empty() {
            return Self { basis: None, variables: Vec::new() };
        }

        // Collect all variables from constraints
        let variables = Self::collect_variables(&constraints);

        // Compute Groebner basis
        let basis = GroebnerBasis::new(&constraints, use_grevlex);

        Self {
            basis: Some(basis),
            variables,
        }
    }

    /// Reduces a polynomial modulo the Groebner basis.
    ///
    /// Returns the canonical form of the polynomial in the quotient ring.
    pub fn reduce(
        &self,
        poly: &MultivariatePolynomial<RationalField, u16>,
    ) -> MultivariatePolynomial<RationalField, u16> {
        match &self.basis {
            Some(basis) => {
                // Reduce polynomial by the Groebner basis
                basis.reduce(poly)
            }
            None => poly.clone(),
        }
    }

    /// Reduces a Symbolica Atom expression.
    pub fn reduce_atom(&self, atom: &Atom) -> Atom {
        // Convert to polynomial
        let poly = atom.to_polynomial(&RationalField::new(), None)
            .expect("expression should be polynomial");

        // Reduce
        let reduced = self.reduce(&poly);

        // Convert back to Atom
        reduced.to_atom()
    }

    /// Collects all variable names from polynomials.
    fn collect_variables(
        polys: &[MultivariatePolynomial<RationalField, u16>],
    ) -> Vec<String> {
        let mut vars = std::collections::HashSet::new();
        for poly in polys {
            for var in poly.get_variables() {
                vars.insert(var.to_string());
            }
        }
        vars.into_iter().collect()
    }
}
```

### Phase 3: Integration with Product Generation

#### 3.1 Update compute_product_expressions

```rust
impl<'a> TraitGenerator<'a> {
    /// Computes product expressions with Groebner simplification.
    fn compute_product_expressions(
        &self,
        type_a: &TypeSpec,
        type_b: &TypeSpec,
        output_type: &TypeSpec,
        kind: SymbolicProductKind,
    ) -> Vec<TokenStream> {
        let symbolic_product = SymbolicProduct::new(self.algebra);
        let expr_simplifier = ExpressionSimplifier::new();

        // Collect constraints for both input types
        let constraint_collector = ProductConstraintCollector::new(self.algebra);
        let mut constraints = Vec::new();
        constraints.extend(constraint_collector.collect_constraints(type_a, "self"));
        constraints.extend(constraint_collector.collect_constraints(type_b, "rhs"));

        // Compute Groebner basis for constraint ideal
        let groebner = GroebnerSimplifier::new(constraints, true); // grevlex for speed

        // Create symbolic field variables
        let a_symbols = symbolic_product.create_field_symbols(type_a, "self");
        let b_symbols = symbolic_product.create_field_symbols(type_b, "rhs");

        // Compute and simplify each output field
        let symbolic_fields = symbolic_product.compute(
            type_a, type_b, output_type, kind, &a_symbols, &b_symbols
        );

        let converter = AtomToRust::new(&[type_a, type_b], &["self", "rhs"]);

        symbolic_fields
            .iter()
            .map(|field| {
                // 1. Basic expansion
                let expanded = expr_simplifier.simplify(&field.expression);

                // 2. Groebner reduction (key step!)
                let reduced = groebner.reduce_atom(&expanded);

                // 3. Final simplification
                let simplified = expr_simplifier.simplify(&reduced);

                converter.convert(&simplified)
            })
            .collect()
    }
}
```

### Phase 4: Caching and Optimization

#### 4.1 Cache Groebner Bases

Computing Groebner bases can be expensive. Cache per type:

```rust
/// Caches Groebner bases for type pairs.
pub struct GroebnerCache {
    /// Cache key: (type_a, type_b) -> basis
    cache: HashMap<(String, String), Arc<GroebnerSimplifier>>,
}

impl GroebnerCache {
    pub fn get_or_compute(
        &mut self,
        type_a: &TypeSpec,
        type_b: &TypeSpec,
        collector: &ProductConstraintCollector,
    ) -> Arc<GroebnerSimplifier> {
        let key = (type_a.name.clone(), type_b.name.clone());

        if let Some(cached) = self.cache.get(&key) {
            return Arc::clone(cached);
        }

        // Collect and compute
        let mut constraints = Vec::new();
        constraints.extend(collector.collect_constraints(type_a, "a"));
        constraints.extend(collector.collect_constraints(type_b, "b"));

        let simplifier = Arc::new(GroebnerSimplifier::new(constraints, true));
        self.cache.insert(key, Arc::clone(&simplifier));
        simplifier
    }
}
```

#### 4.2 Variable Ordering Strategy

Variable ordering affects both basis computation time and result size:

```rust
impl GroebnerSimplifier {
    /// Creates simplifier with optimal variable ordering.
    ///
    /// Strategy: Order output fields last for better elimination.
    pub fn with_ordered_variables(
        constraints: Vec<MultivariatePolynomial<RationalField, u16>>,
        input_vars: &[String],   // Variables from input types
        output_vars: &[String],  // Variables in output expression
    ) -> Self {
        // For elimination ordering: input vars < output vars
        // This prioritizes eliminating output vars in favor of input
        let ordered_vars: Vec<_> = input_vars.iter()
            .chain(output_vars.iter())
            .cloned()
            .collect();

        // Create polynomials with specified variable order
        // (Implementation depends on Symbolica's API)

        Self::new(constraints, true)
    }
}
```

## Configuration

**No algebra TOML changes required.** Groebner simplification is an internal codegen optimization, not a user-facing configuration.

### CLI Flag (Optional)

```bash
# Disable Groebner simplification for debugging
clifford-codegen --no-groebner algebras/conformal3.toml
```

### Compile-Time Settings

All settings are constants in the codegen crate (see Phase 2 implementation). Users don't need to configure anything - the defaults are tuned for safety and performance.

## Testing Strategy

### Performance Benchmarks (Critical)

**Goal**: Prove that Groebner reduction improves runtime performance of generated code, not just reduces term count.

#### Benchmark Infrastructure

Generate two versions of each algebra for A/B comparison:

```rust
// In clifford-codegen, generate both versions for benchmarking
pub fn generate_for_benchmark(spec: &AlgebraSpec) -> BenchmarkOutput {
    BenchmarkOutput {
        // Standard generation with Groebner reduction
        optimized: generate_with_groebner(spec),
        // Naive generation without reduction
        naive: generate_without_groebner(spec),
    }
}
```

Create a benchmark harness that compares identical operations:

```rust
// benches/groebner_comparison.rs
use criterion::{black_box, criterion_group, criterion_main, Criterion, BenchmarkId};

// Import both versions
mod optimized {
    include!(concat!(env!("OUT_DIR"), "/optimized/projective3.rs"));
}
mod naive {
    include!(concat!(env!("OUT_DIR"), "/naive/projective3.rs"));
}

fn bench_line_antiwedge(c: &mut Criterion) {
    let mut group = c.benchmark_group("line_antiwedge");

    // Test data: random lines satisfying Plücker constraint
    let lines: Vec<(optimized::Line<f64>, optimized::Line<f64>)> =
        (0..1000).map(|_| (random_line(), random_line())).collect();

    group.bench_function("groebner_optimized", |b| {
        b.iter(|| {
            for (l1, l2) in &lines {
                black_box(l1.antiwedge(l2));
            }
        })
    });

    // Convert to naive types
    let naive_lines: Vec<(naive::Line<f64>, naive::Line<f64>)> =
        lines.iter().map(|(l1, l2)| (to_naive_line(l1), to_naive_line(l2))).collect();

    group.bench_function("naive", |b| {
        b.iter(|| {
            for (l1, l2) in &naive_lines {
                black_box(l1.antiwedge(l2));
            }
        })
    });

    group.finish();
}
```

#### Key Benchmarks to Run

| Benchmark | Operation | Constraint Used | Expected Improvement |
|-----------|-----------|-----------------|---------------------|
| `line_antiwedge` | Line ∨ Line → Point | Plücker | 20-40% faster |
| `line_geometric` | Line × Line → Motor | Plücker | 15-30% faster |
| `motor_compose` | Motor × Motor → Motor | Study condition | 10-25% faster |
| `circle_antiwedge` | Circle ∨ Circle (CGA) | Blade + Versor | 30-50% faster |
| `point_pair_product` | PointPair × PointPair | Null constraint | 20-35% faster |

#### Benchmark Metrics

```rust
/// Collected metrics for each benchmark.
#[derive(Debug)]
pub struct BenchmarkMetrics {
    /// Operation name
    pub name: String,
    /// Naive version: terms per output field
    pub naive_terms: Vec<usize>,
    /// Optimized version: terms per output field
    pub optimized_terms: Vec<usize>,
    /// Naive execution time (median, ns)
    pub naive_time_ns: u64,
    /// Optimized execution time (median, ns)
    pub optimized_time_ns: u64,
    /// Speedup factor (naive_time / optimized_time)
    pub speedup: f64,
    /// Term reduction percentage
    pub term_reduction_pct: f64,
}

impl BenchmarkMetrics {
    pub fn report(&self) {
        println!("{}:", self.name);
        println!("  Terms: {} → {} ({:.1}% reduction)",
            self.naive_terms.iter().sum::<usize>(),
            self.optimized_terms.iter().sum::<usize>(),
            self.term_reduction_pct);
        println!("  Time: {}ns → {}ns ({:.2}x speedup)",
            self.naive_time_ns,
            self.optimized_time_ns,
            self.speedup);
    }
}
```

#### Acceptance Criteria

The feature is only acceptable if benchmarks demonstrate:

1. **Positive speedup**: `optimized_time < naive_time` for all constrained type products
2. **Meaningful improvement**: At least 10% speedup for products with constraints
3. **No regression**: Unconstrained products (e.g., Vector × Vector) have < 5% overhead
4. **Correlation**: Term reduction correlates with runtime improvement

```rust
#[test]
fn groebner_provides_speedup() {
    let metrics = run_all_benchmarks();

    for m in &metrics {
        if m.term_reduction_pct > 10.0 {
            // If we reduced terms significantly, we should see speedup
            assert!(
                m.speedup > 1.05,
                "{}: {:.1}% term reduction but only {:.2}x speedup",
                m.name, m.term_reduction_pct, m.speedup
            );
        }
    }
}

#[test]
fn no_regression_for_unconstrained() {
    let metrics = run_all_benchmarks();

    for m in &metrics {
        if m.term_reduction_pct < 1.0 {
            // No reduction means no constraint applied - should not slow down
            assert!(
                m.speedup > 0.95,
                "{}: no term reduction and {:.2}x slowdown (regression)",
                m.name, m.speedup
            );
        }
    }
}
```

#### Static Analysis: Term Count Comparison

Before running runtime benchmarks, verify term reduction statically:

```rust
#[test]
fn symbolica_term_counts_reduced() {
    let spec = parse_spec(include_str!("algebras/projective3.toml")).unwrap();
    let algebra = Algebra::pga(3);

    let line = spec.types.iter().find(|t| t.name == "Line").unwrap();
    let point = spec.types.iter().find(|t| t.name == "Point").unwrap();

    // Compute Line ∨ Line → Point with and without Groebner
    let naive_terms = compute_product_terms_naive(&algebra, line, line, point, ProductKind::Antiwedge);
    let optimized_terms = compute_product_terms_groebner(&algebra, line, line, point, ProductKind::Antiwedge);

    let naive_total: usize = naive_terms.iter().sum();
    let optimized_total: usize = optimized_terms.iter().sum();

    println!("Line ∨ Line → Point:");
    println!("  Naive: {} terms", naive_total);
    println!("  Optimized: {} terms", optimized_total);
    println!("  Reduction: {:.1}%", 100.0 * (1.0 - optimized_total as f64 / naive_total as f64));

    assert!(
        optimized_total < naive_total,
        "Groebner should reduce term count: {} < {}",
        optimized_total, naive_total
    );
}
```

#### Codegen Time Overhead

Ensure Groebner computation doesn't excessively slow down code generation:

```rust
#[test]
fn codegen_time_acceptable() {
    let spec = parse_spec(include_str!("algebras/projective3.toml")).unwrap();

    let start_naive = Instant::now();
    let _naive = generate_without_groebner(&spec);
    let naive_duration = start_naive.elapsed();

    let start_opt = Instant::now();
    let _optimized = generate_with_groebner(&spec);
    let opt_duration = start_opt.elapsed();

    let overhead = opt_duration.as_secs_f64() / naive_duration.as_secs_f64();

    println!("Codegen time: {:.2}ms naive, {:.2}ms optimized ({:.1}x overhead)",
        naive_duration.as_secs_f64() * 1000.0,
        opt_duration.as_secs_f64() * 1000.0,
        overhead);

    // Codegen overhead should be < 3x (one-time cost, runtime benefit is ongoing)
    assert!(
        overhead < 3.0,
        "Codegen overhead too high: {:.1}x",
        overhead
    );
}
```

#### CI Integration

Add benchmark comparison to CI:

```yaml
# .github/workflows/benchmark.yml
name: Performance Benchmarks

on:
  pull_request:
    paths:
      - 'crates/clifford-codegen/src/symbolic/**'

jobs:
  benchmark:
    runs-on: ubuntu-latest
    steps:
      - uses: actions/checkout@v4

      - name: Run benchmarks
        run: cargo bench --bench groebner_comparison -- --save-baseline pr

      - name: Compare with main
        run: |
          git fetch origin main
          git checkout origin/main
          cargo bench --bench groebner_comparison -- --save-baseline main
          git checkout -
          cargo bench --bench groebner_comparison -- --baseline main --compare

      - name: Assert no regression
        run: cargo test --test benchmark_assertions
```

#### Example Expected Output

```
Running Groebner optimization benchmarks...

line_antiwedge:
  Terms: 36 → 24 (33.3% reduction)
  Time: 42ns → 31ns (1.35x speedup)  ✓

motor_compose:
  Terms: 64 → 48 (25.0% reduction)
  Time: 89ns → 71ns (1.25x speedup)  ✓

vector_geometric (no constraints):
  Terms: 9 → 9 (0.0% reduction)
  Time: 12ns → 12ns (1.00x speedup)  ✓ (no regression)

circle_antiwedge (CGA):
  Terms: 121 → 72 (40.5% reduction)
  Time: 156ns → 98ns (1.59x speedup)  ✓

Summary:
  4/4 benchmarks passed
  Average speedup for constrained types: 1.40x
  No regressions detected
```

### Correctness Tests

```rust
proptest! {
    /// Simplified product equals naive product numerically.
    #[test]
    fn groebner_simplified_matches_naive(
        l1 in arb_pga_line(),
        l2 in arb_pga_line(),
    ) {
        let simplified = l1.antiwedge(&l2);
        let naive = line_antiwedge_naive(&l1, &l2);

        prop_assert!(
            relative_eq!(simplified, naive, epsilon = RELATIVE_EQ_EPS, max_relative = RELATIVE_EQ_EPS),
            "Simplified: {:?}, Naive: {:?}", simplified, naive
        );
    }
}
```

### Reduction Verification

```rust
#[test]
fn symbolica_plucker_constraint_reduces_expression() {
    let algebra = Algebra::pga(3);
    let collector = ProductConstraintCollector::new(&algebra);

    // Get Line type spec
    let line_spec = /* ... */;

    let constraints = collector.collect_constraints(&line_spec, "l");
    let groebner = GroebnerSimplifier::new(constraints, true);

    // Expression containing Plücker terms
    let expr = parse_atom("l_e01 * l_e23 + l_e02 * l_e31 + l_e03 * l_e12");

    let reduced = groebner.reduce_atom(&expr);

    // Should reduce to 0 (Plücker constraint)
    assert_eq!(reduced, Atom::num(0), "Plücker constraint should reduce to zero");
}

#[test]
fn symbolica_groebner_reduces_product_terms() {
    let algebra = Algebra::pga(3);
    let spec = parse_spec(include_str!("algebras/projective3.toml")).unwrap();
    let line = spec.types.iter().find(|t| t.name == "Line").unwrap();

    // Count terms before and after Groebner reduction
    let symbolic = SymbolicProduct::new(&algebra);
    let a_symbols = symbolic.create_field_symbols(line, "a");
    let b_symbols = symbolic.create_field_symbols(line, "b");

    // Build antiwedge expression
    let raw_expr = /* compute line ∨ line */;
    let raw_terms = count_terms(&raw_expr);

    // Apply Groebner reduction
    let collector = ProductConstraintCollector::new(&algebra);
    let constraints = vec![
        collector.collect_constraints(line, "a"),
        collector.collect_constraints(line, "b"),
    ].concat();

    let groebner = GroebnerSimplifier::new(constraints, true);
    let reduced = groebner.reduce_atom(&raw_expr);
    let reduced_terms = count_terms(&reduced);

    assert!(
        reduced_terms < raw_terms,
        "Groebner should reduce term count: {} -> {}", raw_terms, reduced_terms
    );
}
```

### Benchmarks

```rust
#[bench]
fn bench_codegen_with_groebner(b: &mut Bencher) {
    let spec = parse_spec(include_str!("algebras/projective3.toml")).unwrap();

    b.iter(|| {
        generate_with_options(&spec, CodegenOptions {
            groebner_simplification: true,
            ..Default::default()
        })
    });
}

#[bench]
fn bench_codegen_without_groebner(b: &mut Bencher) {
    let spec = parse_spec(include_str!("algebras/projective3.toml")).unwrap();

    b.iter(|| {
        generate_with_options(&spec, CodegenOptions {
            groebner_simplification: false,
            ..Default::default()
        })
    });
}
```

## Expected Impact

### Term Reduction by Type

| Product | Type Constraints | Before | After | Reduction |
|---------|------------------|--------|-------|-----------|
| Line ∨ Line (PGA) | Plücker | ~36 | ~24 | 33% |
| Motor × Motor | Study condition | ~64 | ~48 | 25% |
| Circle ∨ Circle (CGA) | Blade + Versor | ~100+ | ~60 | 40% |
| PointPair × PointPair | Null constraint | ~81 | ~54 | 33% |

### Code Quality

**Before (naive):**
```rust
// Line meet - all 36 blade combinations
Point::new(
    l1.e01() * l2.e0123() + l1.e02() * l2.e0231() + ...  // 36 terms
)
```

**After (Groebner reduced):**
```rust
// Line meet - constraint-aware
Point::new(
    l1.e01() * l2.e23() - l1.e23() * l2.e01() + ...  // ~24 terms
)
```

## Implementation Phases

### Phase 1: Core Infrastructure ✅
- [x] `ProductConstraintCollector` to gather constraints per type
- [x] `GroebnerSimplifier` wrapper around Symbolica's `GroebnerBasis`
- [x] Atom ↔ MultivariatePolynomial conversion utilities

### Phase 2: Numerical Precision Safeguards ✅
- [x] `CoefficientValidator` for safe coefficient checking (integrated in `validate_coefficients`)
- [x] `ReductionResult` enum with `Reduced` and `Fallback` variants
- [x] Fallback mechanism when reduction fails validation
- [x] Term count validation (never increase terms)
- [ ] `StructureAnalyzer` for expression stability analysis (deferred - not needed for initial implementation)
- [ ] `NumericalVerifier` for high-precision equivalence testing (deferred - covered by property tests)

### Phase 3: Integration ✅
- [x] Update `compute_product_expressions()` to use safe Groebner reduction
- [x] Add `--no-groebner` CLI flag for disabling
- [x] `TraitsGenerator::with_options()` for explicit Groebner control
- [ ] Add debug assertion generation for compile-time verification (deferred)

### Phase 4: Optimization ✅
- [x] Implement `GroebnerCache` for type pair caching
- [x] GrevLex ordering for faster Groebner computation
- [x] Max constraints safeguard (`MAX_CONSTRAINT_POLYNOMIALS = 20`)

### Phase 5: Benchmarking ✅
- [x] Static term count reduction tests (verify optimization in codegen)
- [x] Runtime benchmarks via existing infrastructure (`benches/projective.rs`, etc.)
- [-] A/B dual codegen: Not needed - use git baseline comparison instead
- [-] CI benchmark workflow: Skipped - benchmarks too slow for CI

### Phase 6: Correctness Testing ✅
- [x] Property tests: reduction idempotence, semantic preservation
- [x] Reduction verification tests (constraints reduce to zero)
- [x] Coefficient validation (integrated in reduce_safe)
- [x] Edge case tests (zero, constants, single variables, many constraints, redundant constraints, nested products)

## Files Changed

| File | Action | Description |
|------|--------|-------------|
| `crates/clifford-codegen/src/symbolic/groebner.rs` | Created | `GroebnerSimplifier`, `ProductConstraintCollector`, `GroebnerCache`, `ReductionResult` (952 lines) |
| `crates/clifford-codegen/src/symbolic/mod.rs` | Updated | Export new Groebner types |
| `crates/clifford-codegen/src/codegen/traits.rs` | Updated | Integrate Groebner in `compute_product_expressions`, add `with_options` constructor |
| `crates/clifford-codegen/src/main.rs` | Updated | Add `--no-groebner` CLI flag |

### Deferred/Skipped Files
| File | Status | Reason |
|------|--------|--------|
| `src/symbolic/precision.rs` | Integrated | Coefficient validation in `GroebnerSimplifier::validate_coefficients` |
| `src/codegen/debug_verify.rs` | Skipped | Debug assertions not needed |
| `benches/groebner_comparison.rs` | Skipped | Use existing benchmarks with git baselines instead |
| `.github/workflows/benchmark.yml` | Skipped | Benchmarks too slow for CI |

## Dependencies

- **Symbolica** (existing): `symbolica::poly::groebner::GroebnerBasis`
- **Criterion** (existing, dev): For A/B performance benchmarks
- **rug** (optional, dev): GNU MPFR bindings for high-precision numerical verification
  - Only needed if `verify_numerically = true`
  - Can use pure-Rust alternative like `dashu` if MPFR unavailable
- **thiserror** (existing): For `CoefficientError` derive

## Risks and Mitigations

### Risk: Groebner Computation Time

Groebner basis computation can be expensive for large constraint systems.

**Mitigation:**
1. Cache bases per type pair
2. Use grevlex ordering (much faster than lex)
3. Add max constraint safeguard
4. Allow disabling via CLI for debugging

### Risk: Symbolica API Stability

Symbolica is actively developed; API may change.

**Mitigation:**
1. Wrap Symbolica calls in adapter layer
2. Pin Symbolica version in Cargo.toml
3. Add integration tests that fail on API breakage

### Risk: Numerical Precision

Groebner bases use exact rational arithmetic; conversion to float may introduce errors.

**This is the most critical risk and requires detailed analysis.**

#### Sources of Precision Loss

**1. Rational-to-Float Coefficient Conversion**

Groebner reduction operates over ℚ (rationals) but generated code uses `f32`/`f64`. Some rationals have no exact binary representation:

```
1/3 = 0.333... (infinite binary expansion)
1/5 = 0.2 (finite decimal, infinite binary)
1/7 = 0.142857... (repeating)
```

If reduction produces `(1/3) * a * b`, the generated code `a * b / 3.0` introduces representation error.

**2. Coefficient Explosion During Reduction**

Groebner basis computation can produce large intermediate coefficients:

```
Original: a*b + c*d
After reduction: (12345/6789) * a*b + (-9876/5432) * c*d
```

Large numerators/denominators may overflow or lose precision when converted to float.

**3. Expression Restructuring**

Algebraically equivalent expressions have different numerical properties:

| Original | Reduced | Numerical Issue |
|----------|---------|-----------------|
| `a² - b²` | `(a+b)(a-b)` | Catastrophic cancellation when a ≈ b |
| `a + b + c` | `a + (b + c)` | Different rounding accumulation |
| `2*a*b` | `a*b + a*b` | Double rounding vs single |

**4. Hidden Cancellation**

Reduction might eliminate terms that provided numerical stability:

```
Original: (a*b + ε) - (c*d + ε)  // ε terms stabilize near-cancellation
Reduced:  a*b - c*d              // lost stability term
```

**5. Division Introduction**

Groebner reduction might introduce divisions not present in the original:

```
Original: 2*a*b
Reduced:  a*b / (1/2)  →  a*b * 2.0  // OK
Reduced:  a*b / (3/7)  →  a*b * (7.0/3.0)  // precision loss
```

#### Mitigation Strategy: Coefficient Validation

**Principle: Only accept reductions with "safe" coefficients.**

```rust
/// Validates that a polynomial has numerically safe coefficients.
#[derive(Debug, Clone)]
pub struct CoefficientValidator {
    /// Maximum allowed numerator absolute value.
    max_numerator: i64,
    /// Maximum allowed denominator value.
    max_denominator: i64,
    /// Allowed "simple" denominators (exact in binary float).
    simple_denominators: HashSet<i64>,
}

impl Default for CoefficientValidator {
    fn default() -> Self {
        Self {
            max_numerator: 1_000_000,      // Fits in f32 mantissa
            max_denominator: 1_000_000,
            // Powers of 2 are exactly representable
            simple_denominators: [1, 2, 4, 8, 16, 32, 64, 128, 256, 512, 1024]
                .into_iter().collect(),
        }
    }
}

impl CoefficientValidator {
    /// Checks if a rational coefficient is safe for float conversion.
    pub fn is_safe(&self, numerator: i64, denominator: i64) -> bool {
        // Check bounds
        if numerator.abs() > self.max_numerator {
            return false;
        }
        if denominator > self.max_denominator {
            return false;
        }

        // Check if denominator is "simple" (power of 2)
        // Non-simple denominators introduce representation error
        if !self.simple_denominators.contains(&denominator) {
            // Allow if denominator divides numerator evenly
            if numerator % denominator != 0 {
                return false;
            }
        }

        true
    }

    /// Validates all coefficients in a polynomial.
    pub fn validate_polynomial(
        &self,
        poly: &MultivariatePolynomial<RationalField, u16>,
    ) -> Result<(), CoefficientError> {
        for (coef, _monomial) in poly.terms() {
            let (num, den) = coef.to_numerator_denominator();
            if !self.is_safe(num, den) {
                return Err(CoefficientError::UnsafeCoefficient {
                    numerator: num,
                    denominator: den,
                });
            }
        }
        Ok(())
    }
}

#[derive(Debug, thiserror::Error)]
pub enum CoefficientError {
    #[error("Unsafe coefficient {numerator}/{denominator}: may lose precision")]
    UnsafeCoefficient { numerator: i64, denominator: i64 },
}
```

#### Mitigation Strategy: Fallback on Validation Failure

If reduced expression has unsafe coefficients, fall back to original:

```rust
impl GroebnerSimplifier {
    /// Reduces with coefficient validation and fallback.
    pub fn reduce_safe(
        &self,
        expr: &Atom,
        validator: &CoefficientValidator,
    ) -> ReductionResult {
        let poly = expr.to_polynomial(&RationalField::new(), None)
            .expect("expression should be polynomial");

        let reduced = self.reduce(&poly);

        // Validate reduced coefficients
        match validator.validate_polynomial(&reduced) {
            Ok(()) => ReductionResult::Reduced(reduced.to_atom()),
            Err(e) => {
                // Log warning and return original
                log::warn!(
                    "Groebner reduction produced unsafe coefficients: {}. \
                     Falling back to original expression.",
                    e
                );
                ReductionResult::Fallback {
                    original: expr.clone(),
                    reason: e.to_string(),
                }
            }
        }
    }
}

pub enum ReductionResult {
    /// Successfully reduced with safe coefficients.
    Reduced(Atom),
    /// Fell back to original due to unsafe coefficients.
    Fallback { original: Atom, reason: String },
}
```

#### Mitigation Strategy: Expression Equivalence Verification

Verify numerical equivalence at codegen time using high-precision arithmetic:

```rust
use rug::{Float, Rational};  // GNU MPFR bindings for arbitrary precision

/// Verifies that two expressions are numerically equivalent.
pub struct NumericalVerifier {
    /// Precision in bits for verification (e.g., 256).
    precision: u32,
    /// Number of random test points.
    num_samples: usize,
    /// Relative tolerance for equivalence.
    tolerance: f64,
}

impl NumericalVerifier {
    /// Verifies original and reduced expressions are equivalent.
    pub fn verify_equivalence(
        &self,
        original: &Atom,
        reduced: &Atom,
        variables: &[String],
    ) -> VerificationResult {
        let mut rng = rand::thread_rng();

        for _ in 0..self.num_samples {
            // Generate random values for variables
            let values: HashMap<String, Float> = variables.iter()
                .map(|v| (v.clone(), self.random_value(&mut rng)))
                .collect();

            // Evaluate both expressions at high precision
            let original_value = self.evaluate(original, &values);
            let reduced_value = self.evaluate(reduced, &values);

            // Check relative error
            let error = self.relative_error(&original_value, &reduced_value);
            if error > self.tolerance {
                return VerificationResult::Failed {
                    sample: values.iter()
                        .map(|(k, v)| (k.clone(), v.to_f64()))
                        .collect(),
                    original_value: original_value.to_f64(),
                    reduced_value: reduced_value.to_f64(),
                    relative_error: error,
                };
            }
        }

        VerificationResult::Passed
    }

    fn random_value(&self, rng: &mut impl Rng) -> Float {
        // Use range that exercises both small and large values
        let exp: i32 = rng.gen_range(-10..10);
        let mantissa: f64 = rng.gen_range(-1.0..1.0);
        Float::with_val(self.precision, mantissa * 10f64.powi(exp))
    }

    fn evaluate(&self, expr: &Atom, values: &HashMap<String, Float>) -> Float {
        // Recursive evaluation at high precision
        // (Implementation uses Symbolica's evaluation or custom walker)
        todo!()
    }

    fn relative_error(&self, a: &Float, b: &Float) -> f64 {
        if a.is_zero() && b.is_zero() {
            return 0.0;
        }
        let diff = Float::with_val(self.precision, a - b).abs();
        let max_abs = a.clone().abs().max(&b.clone().abs());
        (diff / max_abs).to_f64()
    }
}

pub enum VerificationResult {
    Passed,
    Failed {
        sample: HashMap<String, f64>,
        original_value: f64,
        reduced_value: f64,
        relative_error: f64,
    },
}
```

#### Mitigation Strategy: Preserve Expression Structure

Reject reductions that introduce problematic patterns:

```rust
/// Analyzes expression structure for numerical stability.
pub struct StructureAnalyzer;

impl StructureAnalyzer {
    /// Checks if reduction introduced problematic patterns.
    pub fn check_stability(
        &self,
        original: &Atom,
        reduced: &Atom,
    ) -> StabilityCheck {
        let mut issues = Vec::new();

        // Check for introduced divisions
        if self.has_more_divisions(reduced, original) {
            issues.push(StabilityIssue::IntroducedDivision);
        }

        // Check for potential catastrophic cancellation
        // (subtraction of similar-magnitude terms not in original)
        if self.has_new_subtractions(reduced, original) {
            issues.push(StabilityIssue::PotentialCancellation);
        }

        // Check for deeper nesting (more rounding accumulation)
        let orig_depth = self.expression_depth(original);
        let red_depth = self.expression_depth(reduced);
        if red_depth > orig_depth + 2 {
            issues.push(StabilityIssue::IncreasedNesting {
                original: orig_depth,
                reduced: red_depth,
            });
        }

        if issues.is_empty() {
            StabilityCheck::Stable
        } else {
            StabilityCheck::Unstable(issues)
        }
    }

    fn has_more_divisions(&self, reduced: &Atom, original: &Atom) -> bool {
        self.count_divisions(reduced) > self.count_divisions(original)
    }

    fn count_divisions(&self, expr: &Atom) -> usize {
        // Walk expression tree counting Div nodes
        todo!()
    }

    fn has_new_subtractions(&self, reduced: &Atom, original: &Atom) -> bool {
        // Detect subtraction patterns not in original
        todo!()
    }

    fn expression_depth(&self, expr: &Atom) -> usize {
        // Compute AST depth
        todo!()
    }
}

pub enum StabilityCheck {
    Stable,
    Unstable(Vec<StabilityIssue>),
}

pub enum StabilityIssue {
    IntroducedDivision,
    PotentialCancellation,
    IncreasedNesting { original: usize, reduced: usize },
}
```

#### Integration: Safe Reduction Pipeline

Combine all mitigations into a safe reduction pipeline:

```rust
impl GroebnerSimplifier {
    /// Performs safe Groebner reduction with all validations.
    pub fn reduce_with_validation(
        &self,
        expr: &Atom,
        config: &SafeReductionConfig,
    ) -> SafeReductionResult {
        // Step 1: Perform reduction
        let reduced = self.reduce_atom(expr);

        // Step 2: Validate coefficients
        let coef_check = config.coefficient_validator.validate_atom(&reduced);
        if let Err(e) = coef_check {
            return SafeReductionResult::fallback(expr.clone(), format!("coefficient: {}", e));
        }

        // Step 3: Check structural stability
        let stability = config.structure_analyzer.check_stability(expr, &reduced);
        if let StabilityCheck::Unstable(issues) = stability {
            if config.reject_unstable {
                return SafeReductionResult::fallback(
                    expr.clone(),
                    format!("stability: {:?}", issues),
                );
            }
            // Log warning but continue
            log::warn!("Reduced expression has stability issues: {:?}", issues);
        }

        // Step 4: Numerical verification (optional, expensive)
        if config.verify_numerically {
            let variables = self.extract_variables(&reduced);
            match config.verifier.verify_equivalence(expr, &reduced, &variables) {
                VerificationResult::Passed => {}
                VerificationResult::Failed { relative_error, .. } => {
                    return SafeReductionResult::fallback(
                        expr.clone(),
                        format!("numerical verification failed: error={:.2e}", relative_error),
                    );
                }
            }
        }

        // All checks passed
        SafeReductionResult::Success {
            reduced,
            term_reduction: self.count_terms(expr) - self.count_terms(&reduced),
        }
    }
}

pub struct SafeReductionConfig {
    pub coefficient_validator: CoefficientValidator,
    pub structure_analyzer: StructureAnalyzer,
    pub verifier: NumericalVerifier,
    pub reject_unstable: bool,
    pub verify_numerically: bool,
}

pub enum SafeReductionResult {
    Success {
        reduced: Atom,
        term_reduction: i32,
    },
    Fallback {
        original: Atom,
        reason: String,
    },
}
```

#### Configuration

**No TOML spec changes required.** All precision settings are compile-time constants in the codegen crate:

```rust
// In clifford-codegen/src/symbolic/precision.rs

/// Default configuration - tuned for safety with good performance.
pub const DEFAULT_CONFIG: SafeReductionConfig = SafeReductionConfig {
    max_numerator: 1_000_000,
    max_denominator: 1_000_000,
    require_power_of_two_denominators: false,
    reject_unstable_patterns: true,
    // Numerical verification only in CI/release, not normal codegen
    verify_numerically: cfg!(feature = "verify-reductions"),
};
```

For users who need to adjust settings, they modify the codegen crate directly - this is an implementation detail, not a user-facing configuration.

#### Compile-Time Guarantees

For maximum safety, generate both versions and verify at compile time:

```rust
// Generated code with debug assertions
impl Line {
    #[inline]
    pub fn antiwedge(&self, other: &Line) -> Point {
        // Groebner-reduced expression
        let result = Point::new(
            self.e01() * other.e23() - self.e23() * other.e01() + ...
        );

        #[cfg(debug_assertions)]
        {
            // Naive expression for verification
            let naive = Point::new(
                self.e01() * other.e0123() + ...  // all 36 terms
            );
            debug_assert!(
                relative_eq!(result, naive, epsilon = 1e-10, max_relative = 1e-10),
                "Groebner reduction mismatch: reduced={:?}, naive={:?}",
                result, naive
            );
        }

        result
    }
}
```

## Success Criteria

**Performance (Required)**:
1. ≥10% runtime speedup for products involving constrained types (Line, Motor, Circle)
2. <5% regression for products without constraints (Vector, Scalar)
3. Codegen time overhead <3x

**Correctness (Required)**:
4. All property tests pass (simplified = naive numerically within tolerance)
5. No precision loss beyond 1e-10 relative error

**Quality (Required)**:
6. Term count reduction ≥20% for constrained type products
7. All coefficient validation checks pass (no unsafe rationals)

## References

- [Symbolica Groebner Documentation](https://symbolica.io/docs/polynomials.html)
- [Symbolica Rust API](https://docs.rs/symbolica/latest/symbolica/poly/groebner/index.html)
- [Groebner Bases - Wikipedia](https://en.wikipedia.org/wiki/Gr%C3%B6bner_basis)
- [PRD-14.8: Constraint Engine](prd-14.8-constraint-engine.md)
- [PRD-16.3: Expression Simplification](prd-16.3-expression-simplification.md)
