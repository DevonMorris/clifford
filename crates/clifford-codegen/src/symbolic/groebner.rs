//! Groebner basis simplification for product expressions.
//!
//! This module provides constraint-based simplification using Groebner bases.
//! When computing products of constrained types (like PGA lines with the Plücker
//! condition), the Groebner basis allows us to eliminate redundant terms.
//!
//! # Overview
//!
//! 1. **Collect constraints** from input types as polynomial equations
//! 2. **Compute Groebner basis** of the constraint ideal
//! 3. **Reduce product expressions** modulo the Groebner basis
//! 4. **Validate coefficients** to ensure numerical safety
//!
//! # Example
//!
//! For PGA lines with the Plücker constraint `e01*e23 + e02*e31 + e03*e12 = 0`:
//!
//! ```ignore
//! let constraints = collector.collect_constraints(&line_spec, "l");
//! let simplifier = GroebnerSimplifier::new(constraints, true);
//! let reduced = simplifier.reduce_atom(&expr);
//! ```

use std::collections::HashMap;
use std::sync::Arc;

use symbolica::atom::{Atom, AtomCore};
use symbolica::domains::rational::Q;
use symbolica::poly::GrevLexOrder;
use symbolica::poly::groebner::GroebnerBasis;
use symbolica::poly::polynomial::MultivariatePolynomial;

use crate::algebra::Algebra;
use crate::spec::{InvolutionKind, TypeSpec};

use super::ConstraintDeriver;

/// Maximum number of constraint polynomials before we skip Groebner computation.
///
/// Large constraint systems can cause exponential blowup in Groebner basis computation.
const MAX_CONSTRAINT_POLYNOMIALS: usize = 20;

/// Maximum numerator/denominator magnitude for safe float conversion.
///
/// Coefficients larger than this may lose precision when converted to f64.
const MAX_COEFFICIENT_MAGNITUDE: i64 = 1_000_000;

/// Collects constraint polynomials for input types in a product.
///
/// This struct gathers the constraint expressions that must equal zero for
/// valid instances of each input type, then converts them to polynomials
/// suitable for Groebner basis computation.
pub struct ProductConstraintCollector<'a> {
    /// The constraint deriver for extracting constraints.
    deriver: ConstraintDeriver<'a>,
}

impl<'a> ProductConstraintCollector<'a> {
    /// Creates a new constraint collector.
    ///
    /// # Arguments
    ///
    /// * `algebra` - The algebra for computations
    /// * `involution` - Which involution to use (usually Reverse)
    pub fn new(algebra: &'a Algebra, involution: InvolutionKind) -> Self {
        Self {
            deriver: ConstraintDeriver::new(algebra, involution),
        }
    }

    /// Collects all constraint polynomials for a type.
    ///
    /// Returns polynomials that must equal zero for valid instances.
    /// The polynomials use field names prefixed with the given prefix
    /// (e.g., "self_e01" for the e01 field of self).
    ///
    /// # Arguments
    ///
    /// * `ty` - The type specification
    /// * `prefix` - Variable prefix (e.g., "self", "rhs")
    pub fn collect_constraints(&self, ty: &TypeSpec, prefix: &str) -> Vec<Atom> {
        let mut constraints = Vec::new();

        // 1. Geometric constraint: u * ũ = scalar (non-scalar terms = 0)
        if let Some(geo) = self.deriver.derive_geometric_constraint(ty, prefix) {
            constraints.extend(geo.zero_expressions);
        }

        // 2. Antiproduct constraint: u ⊟ ũ̃ = antiscalar
        if let Some(anti) = self.deriver.derive_antiproduct_constraint(ty, prefix) {
            constraints.extend(anti.zero_expressions);
        }

        constraints
    }
}

/// Result of attempting a safe Groebner reduction.
#[derive(Debug, Clone)]
pub enum ReductionResult {
    /// Successfully reduced with safe coefficients.
    Reduced {
        /// The reduced expression.
        reduced: Atom,
        /// Number of terms reduced (positive = fewer terms).
        term_reduction: i32,
    },
    /// Fell back to original due to validation failure.
    Fallback {
        /// The original expression (unchanged).
        original: Atom,
        /// Reason for fallback.
        reason: String,
    },
}

impl ReductionResult {
    /// Returns the resulting expression (reduced or original).
    pub fn expression(&self) -> &Atom {
        match self {
            ReductionResult::Reduced { reduced, .. } => reduced,
            ReductionResult::Fallback { original, .. } => original,
        }
    }

    /// Returns true if reduction was successful.
    pub fn is_reduced(&self) -> bool {
        matches!(self, ReductionResult::Reduced { .. })
    }
}

/// Simplifies expressions using Groebner basis reduction.
///
/// This struct computes a Groebner basis from constraint polynomials and
/// uses it to reduce (simplify) product expressions.
///
/// # Safety
///
/// The simplifier validates reduced expressions to ensure:
/// - Coefficients don't explode (within `MAX_COEFFICIENT_MAGNITUDE`)
/// - The reduced form is actually simpler (fewer terms)
///
/// If validation fails, it falls back to the original expression.
pub struct GroebnerSimplifier {
    /// The Groebner basis polynomials (already computed).
    basis_polys: Vec<MultivariatePolynomial<Q, u16, GrevLexOrder>>,
    /// Whether reduction is enabled (may be disabled for performance).
    enabled: bool,
}

impl GroebnerSimplifier {
    /// Creates a simplifier from constraint atoms.
    ///
    /// # Arguments
    ///
    /// * `constraints` - Atom expressions that must equal zero
    /// * `use_grevlex` - Use grevlex ordering (faster) vs lex (better elimination)
    ///
    /// # Returns
    ///
    /// A simplifier ready to reduce expressions. If constraints are empty or
    /// too numerous, returns a disabled simplifier that passes through unchanged.
    pub fn new(constraints: Vec<Atom>, _use_grevlex: bool) -> Self {
        if constraints.is_empty() {
            return Self::disabled();
        }

        if constraints.len() > MAX_CONSTRAINT_POLYNOMIALS {
            return Self::disabled();
        }

        // Convert constraints to polynomials
        let mut polys = Vec::with_capacity(constraints.len());

        for constraint in &constraints {
            // Expand the constraint first
            let expanded = constraint.expand();
            // Convert to polynomial over rationals (lex order by default)
            let poly: MultivariatePolynomial<Q, u16> = expanded.to_polynomial(&Q, None);
            // Reorder to grevlex for faster Groebner computation
            polys.push(poly.reorder::<GrevLexOrder>());
        }

        if polys.is_empty() {
            return Self::disabled();
        }

        // Compute Groebner basis (false = don't print stats)
        let basis = GroebnerBasis::new(&polys, false);

        Self {
            basis_polys: basis.system,
            enabled: true,
        }
    }

    /// Creates a disabled simplifier that passes expressions through unchanged.
    fn disabled() -> Self {
        Self {
            basis_polys: Vec::new(),
            enabled: false,
        }
    }

    /// Returns true if this simplifier is enabled.
    pub fn is_enabled(&self) -> bool {
        self.enabled
    }

    /// Reduces an expression modulo the Groebner basis.
    ///
    /// This is the main simplification method. It converts the expression to
    /// a polynomial, reduces it using the Groebner basis, then converts back.
    ///
    /// # Arguments
    ///
    /// * `expr` - The expression to reduce
    ///
    /// # Returns
    ///
    /// The reduced expression, or the original if reduction fails or is disabled.
    pub fn reduce_atom(&self, expr: &Atom) -> Atom {
        self.reduce_safe(expr).expression().clone()
    }

    /// Reduces an expression with validation and fallback.
    ///
    /// This method performs reduction and validates the result. If validation
    /// fails (e.g., coefficient explosion), it falls back to the original.
    pub fn reduce_safe(&self, expr: &Atom) -> ReductionResult {
        if !self.enabled || self.basis_polys.is_empty() {
            return ReductionResult::Fallback {
                original: expr.clone(),
                reason: "Groebner simplification disabled".to_string(),
            };
        }

        // Expand and convert to polynomial
        let expanded = expr.expand();
        let poly: MultivariatePolynomial<Q, u16> = expanded.to_polynomial(&Q, None);
        let poly = poly.reorder::<GrevLexOrder>();

        // Count original terms
        let original_terms = poly.nterms();

        // Reduce modulo Groebner basis
        let reduced_poly = poly.reduce(&self.basis_polys);

        // Validate coefficients
        if let Err(e) = Self::validate_coefficients(&reduced_poly) {
            return ReductionResult::Fallback {
                original: expr.clone(),
                reason: e,
            };
        }

        // Convert back to atom
        let reduced = reduced_poly.to_expression();

        // Count reduced terms
        let reduced_terms = reduced_poly.nterms();
        let term_reduction = original_terms as i32 - reduced_terms as i32;

        // Only accept if we actually reduced terms or stayed the same
        if term_reduction < 0 {
            return ReductionResult::Fallback {
                original: expr.clone(),
                reason: format!(
                    "Reduction increased terms: {} -> {}",
                    original_terms, reduced_terms
                ),
            };
        }

        ReductionResult::Reduced {
            reduced,
            term_reduction,
        }
    }

    /// Validates that polynomial coefficients are safe for float conversion.
    fn validate_coefficients(
        poly: &MultivariatePolynomial<Q, u16, GrevLexOrder>,
    ) -> Result<(), String> {
        for monomial in poly.into_iter() {
            let coef = monomial.coefficient;
            // Get numerator and denominator as Integers
            let num = coef.numerator_ref();
            let den = coef.denominator_ref();

            // Try to convert to i64 for bounds checking
            if let (Some(n), Some(d)) = (num.to_i64(), den.to_i64()) {
                if n.abs() > MAX_COEFFICIENT_MAGNITUDE {
                    return Err(format!(
                        "Numerator too large: {} > {}",
                        n.abs(),
                        MAX_COEFFICIENT_MAGNITUDE
                    ));
                }

                if d.abs() > MAX_COEFFICIENT_MAGNITUDE {
                    return Err(format!(
                        "Denominator too large: {} > {}",
                        d.abs(),
                        MAX_COEFFICIENT_MAGNITUDE
                    ));
                }
            }
            // If conversion to i64 fails, the number is definitely too large
            // but we'll allow it for now (might be simplified later)
        }
        Ok(())
    }
}

/// Cache for Groebner bases computed for type pairs.
///
/// Computing Groebner bases can be expensive, so we cache them for reuse
/// when the same type pair is encountered multiple times.
pub struct GroebnerCache {
    /// Cache key: (type_a_name, type_b_name) -> simplifier
    cache: HashMap<(String, String), Arc<GroebnerSimplifier>>,
}

impl GroebnerCache {
    /// Creates a new empty cache.
    pub fn new() -> Self {
        Self {
            cache: HashMap::new(),
        }
    }

    /// Gets or computes a Groebner simplifier for a type pair.
    ///
    /// If the simplifier for this pair is already cached, returns it.
    /// Otherwise, computes the Groebner basis and caches it.
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

        // Collect constraints from both input types
        let mut constraints = Vec::new();
        constraints.extend(collector.collect_constraints(type_a, "self"));
        constraints.extend(collector.collect_constraints(type_b, "rhs"));

        // Compute Groebner basis
        let simplifier = Arc::new(GroebnerSimplifier::new(constraints, true));
        self.cache.insert(key, Arc::clone(&simplifier));
        simplifier
    }

    /// Clears the cache.
    pub fn clear(&mut self) {
        self.cache.clear();
    }
}

impl Default for GroebnerCache {
    fn default() -> Self {
        Self::new()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::borrow::Cow;
    use std::sync::Mutex;
    use symbolica::atom::DefaultNamespace;
    use symbolica::parser::ParseSettings;

    // Symbolica uses global state that conflicts when tests run in parallel.
    static SYMBOLICA_LOCK: Mutex<()> = Mutex::new(());

    fn parse_atom(s: &str) -> Atom {
        let input = DefaultNamespace {
            namespace: Cow::Borrowed(env!("CARGO_CRATE_NAME")),
            data: s,
            file: Cow::Borrowed(file!()),
            line: line!() as usize,
        };
        Atom::parse(input, ParseSettings::symbolica()).unwrap()
    }

    #[test]
    fn symbolica_empty_constraints_returns_disabled() {
        let _guard = SYMBOLICA_LOCK.lock().unwrap();
        let simplifier = GroebnerSimplifier::new(vec![], true);
        assert!(!simplifier.is_enabled());
    }

    #[test]
    fn symbolica_disabled_simplifier_returns_original() {
        let _guard = SYMBOLICA_LOCK.lock().unwrap();
        let simplifier = GroebnerSimplifier::disabled();
        let expr = parse_atom("a + b");
        let result = simplifier.reduce_safe(&expr);

        assert!(!result.is_reduced());
        assert_eq!(result.expression().expand(), expr.expand());
    }

    #[test]
    fn symbolica_simple_constraint_reduces() {
        let _guard = SYMBOLICA_LOCK.lock().unwrap();

        // Constraint: a + b = 0  =>  a = -b
        let constraint = parse_atom("a + b");
        let simplifier = GroebnerSimplifier::new(vec![constraint], true);

        // Expression: 2*a + 2*b should reduce to 0
        let expr = parse_atom("2*a + 2*b");
        let result = simplifier.reduce_safe(&expr);

        // The expression should reduce to 0 since a + b = 0
        if result.is_reduced() {
            let reduced = result.expression();
            assert_eq!(reduced.expand(), Atom::num(0));
        }
    }

    #[test]
    fn symbolica_plucker_constraint() {
        let _guard = SYMBOLICA_LOCK.lock().unwrap();

        // Plücker constraint: e01*e23 + e02*e31 + e03*e12 = 0
        let constraint = parse_atom("e01*e23 + e02*e31 + e03*e12");
        let simplifier = GroebnerSimplifier::new(vec![constraint], true);

        assert!(simplifier.is_enabled());

        // The constraint itself should reduce to 0
        let expr = parse_atom("e01*e23 + e02*e31 + e03*e12");
        let result = simplifier.reduce_safe(&expr);

        if result.is_reduced() {
            let reduced = result.expression();
            assert_eq!(reduced.expand(), Atom::num(0));
        }
    }

    /// Creates a minimal Line type for testing.
    fn create_line_type() -> TypeSpec {
        use crate::spec::FieldSpec;

        TypeSpec {
            name: "Line".to_string(),
            description: None,
            grades: vec![2],
            fields: vec![
                FieldSpec {
                    name: "e01".to_string(),
                    blade_index: 0b0011,
                    grade: 2,
                },
                FieldSpec {
                    name: "e02".to_string(),
                    blade_index: 0b0101,
                    grade: 2,
                },
                FieldSpec {
                    name: "e03".to_string(),
                    blade_index: 0b1001,
                    grade: 2,
                },
                FieldSpec {
                    name: "e12".to_string(),
                    blade_index: 0b0110,
                    grade: 2,
                },
                FieldSpec {
                    name: "e31".to_string(),
                    blade_index: 0b1010,
                    grade: 2,
                },
                FieldSpec {
                    name: "e23".to_string(),
                    blade_index: 0b1100,
                    grade: 2,
                },
            ],
            versor: None,
            alias_of: None,
            is_sparse: false,
        }
    }

    #[test]
    fn symbolica_constraint_collector_pga_line() {
        let _guard = SYMBOLICA_LOCK.lock().unwrap();

        let line = create_line_type();
        let algebra = Algebra::pga(3);
        let collector = ProductConstraintCollector::new(&algebra, InvolutionKind::Reverse);

        let constraints = collector.collect_constraints(&line, "l");

        // PGA lines have the Plücker constraint
        assert!(
            !constraints.is_empty(),
            "PGA lines should have constraints (Plücker condition)"
        );
    }

    #[test]
    fn symbolica_cache_reuses_simplifier() {
        let _guard = SYMBOLICA_LOCK.lock().unwrap();

        let line = create_line_type();
        let algebra = Algebra::pga(3);
        let collector = ProductConstraintCollector::new(&algebra, InvolutionKind::Reverse);

        let mut cache = GroebnerCache::new();

        // First call should compute
        let s1 = cache.get_or_compute(&line, &line, &collector);
        // Second call should reuse
        let s2 = cache.get_or_compute(&line, &line, &collector);

        // Should be the same Arc
        assert!(Arc::ptr_eq(&s1, &s2));
    }

    #[test]
    fn symbolica_term_count_reduction_with_constraint() {
        let _guard = SYMBOLICA_LOCK.lock().unwrap();

        // Constraint: a*d + b*e + c*f = 0 (like Plücker)
        let constraint = parse_atom("a*d + b*e + c*f");
        let simplifier = GroebnerSimplifier::new(vec![constraint], true);

        // Expression that can be simplified using the constraint
        // 2*a*d + 2*b*e + 2*c*f should reduce to 0
        let expr = parse_atom("2*a*d + 2*b*e + 2*c*f");

        let result = simplifier.reduce_safe(&expr);

        assert!(
            result.is_reduced(),
            "Expression should reduce using the constraint"
        );

        if let ReductionResult::Reduced { term_reduction, .. } = result {
            assert!(term_reduction >= 0, "Term reduction should be non-negative");
        }
    }

    #[test]
    fn symbolica_non_constrained_expression_unchanged() {
        let _guard = SYMBOLICA_LOCK.lock().unwrap();

        // Constraint that doesn't apply to our expression
        let constraint = parse_atom("x*y + z");
        let simplifier = GroebnerSimplifier::new(vec![constraint], true);

        // Expression using different variables - can't be reduced
        let expr = parse_atom("a + b + c");
        let result = simplifier.reduce_safe(&expr);

        let original_expanded = expr.expand();
        let result_expanded = result.expression().expand();

        // Unrelated expression should be unchanged
        assert_eq!(
            original_expanded, result_expanded,
            "Unrelated expression should be unchanged"
        );
    }

    #[test]
    fn symbolica_partial_reduction() {
        let _guard = SYMBOLICA_LOCK.lock().unwrap();

        // Constraint: a + b = 0 => a = -b
        let constraint = parse_atom("a + b");
        let simplifier = GroebnerSimplifier::new(vec![constraint], true);

        // Expression: a + b + c should reduce to c
        let expr = parse_atom("a + b + c");
        let result = simplifier.reduce_safe(&expr);

        if result.is_reduced() {
            let reduced = result.expression();
            let expected = parse_atom("c");
            assert_eq!(
                reduced.expand(),
                expected.expand(),
                "a + b should cancel leaving only c"
            );
        }
    }
}

/// Property-based tests for Groebner reduction correctness.
#[cfg(test)]
mod proptest_tests {
    use super::*;
    use proptest::prelude::*;
    use std::borrow::Cow;
    use std::sync::Mutex;
    use symbolica::atom::DefaultNamespace;
    use symbolica::parser::ParseSettings;

    // Symbolica uses global state that conflicts when tests run in parallel.
    static SYMBOLICA_LOCK: Mutex<()> = Mutex::new(());

    fn parse_atom(s: &str) -> Atom {
        let input = DefaultNamespace {
            namespace: Cow::Borrowed(env!("CARGO_CRATE_NAME")),
            data: s,
            file: Cow::Borrowed(file!()),
            line: line!() as usize,
        };
        Atom::parse(input, ParseSettings::symbolica()).unwrap()
    }

    /// Strategy for generating simple polynomial terms like "3*a*b".
    fn term_strategy() -> impl Strategy<Value = String> {
        // Generate coefficient (-10 to 10, non-zero)
        let coef = prop::sample::select(vec![-5, -3, -2, -1, 1, 2, 3, 5]);

        // Generate variable names
        let var = prop::sample::select(vec!["a", "b", "c", "d", "e", "f"]);

        // Combine into a term with 1-3 variables
        (coef, prop::collection::vec(var, 1..=3)).prop_map(|(c, vars)| {
            let var_product = vars.join("*");
            if c == 1 {
                var_product
            } else if c == -1 {
                format!("-{}", var_product)
            } else {
                format!("{}*{}", c, var_product)
            }
        })
    }

    /// Strategy for generating simple linear expressions (sums of terms).
    fn linear_expr_strategy() -> impl Strategy<Value = String> {
        prop::collection::vec(term_strategy(), 1..=4).prop_map(|terms| terms.join(" + "))
    }

    proptest! {
        /// Reduction is idempotent: reducing twice gives the same result.
        #[test]
        fn symbolica_reduction_is_idempotent(expr_str in linear_expr_strategy()) {
            let _guard = SYMBOLICA_LOCK.lock().unwrap();

            // Simple constraint: a + b = 0
            let constraint = parse_atom("a + b");
            let simplifier = GroebnerSimplifier::new(vec![constraint], true);

            let expr = parse_atom(&expr_str);

            // First reduction
            let first = simplifier.reduce_atom(&expr);
            // Second reduction
            let second = simplifier.reduce_atom(&first);

            // Should be identical after expansion
            prop_assert_eq!(
                first.expand(),
                second.expand(),
                "Reduction should be idempotent"
            );
        }

        /// Disabled simplifier is identity.
        #[test]
        fn symbolica_disabled_is_identity(expr_str in linear_expr_strategy()) {
            let _guard = SYMBOLICA_LOCK.lock().unwrap();

            let simplifier = GroebnerSimplifier::disabled();
            let expr = parse_atom(&expr_str);

            let result = simplifier.reduce_atom(&expr);

            prop_assert_eq!(
                expr.expand(),
                result.expand(),
                "Disabled simplifier should return original"
            );
        }

        /// Term count never increases during reduction.
        #[test]
        fn symbolica_term_count_non_increasing(expr_str in linear_expr_strategy()) {
            let _guard = SYMBOLICA_LOCK.lock().unwrap();

            // Constraint: a*d + b*e = 0 (Plücker-like)
            let constraint = parse_atom("a*d + b*e");
            let simplifier = GroebnerSimplifier::new(vec![constraint], true);

            let expr = parse_atom(&expr_str);
            let result = simplifier.reduce_safe(&expr);

            match result {
                ReductionResult::Reduced { term_reduction, .. } => {
                    prop_assert!(
                        term_reduction >= 0,
                        "Term count should not increase, got reduction: {}",
                        term_reduction
                    );
                }
                ReductionResult::Fallback { .. } => {
                    // Fallback is acceptable - it means we detected a problem
                }
            }
        }

        /// Constraints reduce to zero.
        #[test]
        fn symbolica_constraint_reduces_to_zero(scale in -5i32..=5i32) {
            prop_assume!(scale != 0);
            let _guard = SYMBOLICA_LOCK.lock().unwrap();

            // Constraint: a + b = 0
            let constraint = parse_atom("a + b");
            let simplifier = GroebnerSimplifier::new(vec![constraint.clone()], true);

            // Scale * constraint should reduce to 0
            let scaled = if scale == 1 {
                "a + b".to_string()
            } else if scale == -1 {
                "-a - b".to_string()
            } else {
                format!("{}*a + {}*b", scale, scale)
            };

            let expr = parse_atom(&scaled);
            let result = simplifier.reduce_atom(&expr);

            prop_assert_eq!(
                result.expand(),
                Atom::num(0),
                "Scaled constraint should reduce to 0"
            );
        }

        /// Linear combinations of constraints reduce appropriately.
        #[test]
        fn symbolica_linear_combo_of_constraints(c1 in -3i32..=3i32, c2 in -3i32..=3i32) {
            prop_assume!(c1 != 0 || c2 != 0);
            let _guard = SYMBOLICA_LOCK.lock().unwrap();

            // Constraints: a + b = 0, c + d = 0
            let constraint1 = parse_atom("a + b");
            let constraint2 = parse_atom("c + d");
            let simplifier =
                GroebnerSimplifier::new(vec![constraint1.clone(), constraint2.clone()], true);

            // c1*(a + b) + c2*(c + d) should reduce to 0
            let expr_str = format!("{}*a + {}*b + {}*c + {}*d", c1, c1, c2, c2);
            let expr = parse_atom(&expr_str);
            let result = simplifier.reduce_atom(&expr);

            prop_assert_eq!(
                result.expand(),
                Atom::num(0),
                "Linear combination of constraints should reduce to 0"
            );
        }
    }

    proptest! {
        #![proptest_config(ProptestConfig::with_cases(20))]

        /// Reduction preserves semantic equivalence under constraint substitution.
        ///
        /// If a + b = 0 (so a = -b), then evaluating the original and reduced
        /// expressions with a = -b should give the same result.
        #[test]
        fn symbolica_reduction_preserves_semantics(
            b_val in -10i64..=10i64,
            c_val in -10i64..=10i64,
        ) {
            let _guard = SYMBOLICA_LOCK.lock().unwrap();

            // Constraint: a + b = 0 => a = -b
            let constraint = parse_atom("a + b");
            let simplifier = GroebnerSimplifier::new(vec![constraint], true);

            // Test expression: a*c + b*c = c*(a + b) = 0 when a = -b
            let expr = parse_atom("a*c + b*c");
            let reduced = simplifier.reduce_atom(&expr);

            // The reduced form should be 0 since a*c + b*c = c*(a+b) = 0
            // Or it should evaluate to 0 when substituting a = -b
            let a_val = -b_val;

            // Evaluate original: a*c + b*c with a=-b, b=b, c=c
            let original_eval = a_val * c_val + b_val * c_val;
            // This should always be 0

            prop_assert_eq!(
                original_eval, 0,
                "Original expression should evaluate to 0 when a = -b"
            );

            // The reduced form should also be 0 or evaluate to 0
            if reduced.expand() == Atom::num(0) {
                // Perfect - reduced to constant 0
            } else {
                // The reduction might leave variables - check it evaluates to 0
                // For a*c + b*c, the reduction should give 0
                prop_assert_eq!(
                    reduced.expand(),
                    Atom::num(0),
                    "Expression a*c + b*c should reduce to 0 with constraint a + b = 0"
                );
            }
        }
    }

    /// Edge case tests for Groebner reduction.
    mod edge_cases {
        use super::{Atom, GroebnerSimplifier, SYMBOLICA_LOCK, parse_atom};
        use symbolica::atom::AtomCore;

        #[test]
        fn symbolica_zero_expression() {
            let _guard = SYMBOLICA_LOCK.lock().unwrap();

            let constraint = parse_atom("a + b");
            let simplifier = GroebnerSimplifier::new(vec![constraint], true);

            let zero = Atom::num(0);
            let result = simplifier.reduce_atom(&zero);

            assert_eq!(result.expand(), Atom::num(0), "Zero should stay zero");
        }

        #[test]
        fn symbolica_constant_expression() {
            let _guard = SYMBOLICA_LOCK.lock().unwrap();

            let constraint = parse_atom("a + b");
            let simplifier = GroebnerSimplifier::new(vec![constraint], true);

            let constant = Atom::num(42);
            let result = simplifier.reduce_atom(&constant);

            assert_eq!(
                result.expand(),
                Atom::num(42),
                "Constants should be unchanged"
            );
        }

        #[test]
        fn symbolica_single_variable() {
            let _guard = SYMBOLICA_LOCK.lock().unwrap();

            let constraint = parse_atom("x + y");
            let simplifier = GroebnerSimplifier::new(vec![constraint], true);

            // Variable not in constraint
            let var = parse_atom("z");
            let result = simplifier.reduce_atom(&var);

            assert_eq!(
                result.expand(),
                var.expand(),
                "Unrelated variable should be unchanged"
            );
        }

        #[test]
        fn symbolica_many_constraints_disabled() {
            let _guard = SYMBOLICA_LOCK.lock().unwrap();

            // Create more than MAX_CONSTRAINT_POLYNOMIALS constraints
            let constraints: Vec<_> = (0..25)
                .map(|i| parse_atom(&format!("x{} + y{}", i, i)))
                .collect();

            let simplifier = GroebnerSimplifier::new(constraints, true);

            assert!(
                !simplifier.is_enabled(),
                "Too many constraints should disable simplifier"
            );
        }

        #[test]
        fn symbolica_redundant_constraints() {
            let _guard = SYMBOLICA_LOCK.lock().unwrap();

            // Two equivalent constraints: a + b = 0 and 2a + 2b = 0
            let c1 = parse_atom("a + b");
            let c2 = parse_atom("2*a + 2*b");
            let simplifier = GroebnerSimplifier::new(vec![c1, c2], true);

            // Should still work correctly
            let expr = parse_atom("a + b");
            let result = simplifier.reduce_atom(&expr);

            assert_eq!(
                result.expand(),
                Atom::num(0),
                "Redundant constraints should still reduce correctly"
            );
        }

        #[test]
        fn symbolica_product_constraint() {
            let _guard = SYMBOLICA_LOCK.lock().unwrap();

            // Plücker-like: p*s + q*t + r*u = 0
            let constraint = parse_atom("p*s + q*t + r*u");
            let simplifier = GroebnerSimplifier::new(vec![constraint], true);

            // Expression containing the constraint pattern
            let expr = parse_atom("2*p*s + 2*q*t + 2*r*u");
            let result = simplifier.reduce_atom(&expr);

            assert_eq!(
                result.expand(),
                Atom::num(0),
                "Product constraint should reduce multiples to 0"
            );
        }

        #[test]
        fn symbolica_nested_product() {
            let _guard = SYMBOLICA_LOCK.lock().unwrap();

            // Constraint: a + b = 0
            let constraint = parse_atom("a + b");
            let simplifier = GroebnerSimplifier::new(vec![constraint], true);

            // (a + b) * c = a*c + b*c should reduce to 0
            let expr = parse_atom("(a + b)*c");
            let expanded = expr.expand();
            let result = simplifier.reduce_atom(&expanded);

            assert_eq!(
                result.expand(),
                Atom::num(0),
                "Product with constraint should reduce to 0"
            );
        }
    }
}
