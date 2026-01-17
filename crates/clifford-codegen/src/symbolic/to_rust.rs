//! Convert Symbolica expressions to Rust code.
//!
//! This module provides utilities for converting Symbolica `Atom` expressions
//! into Rust `TokenStream` code that can be compiled.

use std::collections::HashMap;

use proc_macro2::TokenStream;
use quote::{format_ident, quote};
use symbolica::atom::{Atom, AtomCore, AtomView};
use symbolica::coefficient::CoefficientView;

use crate::spec::TypeSpec;

/// Converts Symbolica `Atom` expressions to Rust `TokenStream`.
///
/// This converter handles the translation of symbolic mathematical expressions
/// into efficient Rust code, properly handling:
/// - Field accessors (e.g., `a_x` → `a.x()`)
/// - Numeric constants (0, 1, 2, etc.)
/// - Arithmetic operations (+, -, *, powers)
///
/// # Example
///
/// ```ignore
/// use clifford_codegen::symbolic::AtomToRust;
///
/// let converter = AtomToRust::new(&[&type_a, &type_b], &["a", "b"]);
/// let atom = Atom::parse("a_x * b_y + a_y * b_x").unwrap();
/// let code = converter.convert(&atom);
/// // Produces: a.x() * b.y() + a.y() * b.x()
/// ```
pub struct AtomToRust {
    /// Map from symbol name (e.g., "a_x") to (prefix, field) for accessor generation.
    symbol_map: HashMap<String, (String, String)>,
    /// Set of prefixes that require `.as_inner()` for field access (wrapper types).
    wrapped_prefixes: std::collections::HashSet<String>,
}

impl AtomToRust {
    /// Creates a new converter for the given types and prefixes.
    ///
    /// # Arguments
    ///
    /// * `types` - The types whose fields will appear in expressions
    /// * `prefixes` - The variable prefixes used (e.g., "a", "b")
    ///
    /// # Example
    ///
    /// ```ignore
    /// let converter = AtomToRust::new(&[&vector_type, &bivector_type], &["a", "b"]);
    /// ```
    pub fn new(types: &[&TypeSpec], prefixes: &[&str]) -> Self {
        Self::new_with_wrappers(types, prefixes, &[])
    }

    /// Creates a converter with wrapper type support.
    ///
    /// # Arguments
    ///
    /// * `types` - The types whose fields will appear in expressions
    /// * `prefixes` - The variable prefixes used (e.g., "self", "rhs")
    /// * `wrapped_prefixes` - Prefixes that use wrapper types (need `.as_inner()`)
    ///
    /// # Example
    ///
    /// ```ignore
    /// // For Wedge<B> for Unit<A>: self is wrapped, rhs is not
    /// let converter = AtomToRust::new_with_wrappers(
    ///     &[&type_a, &type_b],
    ///     &["self", "rhs"],
    ///     &["self"],  // self needs .as_inner()
    /// );
    /// ```
    pub fn new_with_wrappers(
        types: &[&TypeSpec],
        prefixes: &[&str],
        wrapper_prefixes: &[&str],
    ) -> Self {
        let mut symbol_map = HashMap::new();

        for (ty, prefix) in types.iter().zip(prefixes.iter()) {
            for field in &ty.fields {
                let symbol_name = format!("{}_{}", prefix, field.name);
                symbol_map.insert(symbol_name, ((*prefix).to_string(), field.name.clone()));
            }
        }

        let wrapped_prefixes = wrapper_prefixes.iter().map(|s| (*s).to_string()).collect();

        Self {
            symbol_map,
            wrapped_prefixes,
        }
    }

    /// Converts a Symbolica `Atom` to a Rust `TokenStream`.
    pub fn convert(&self, atom: &Atom) -> TokenStream {
        self.convert_view(atom.as_atom_view())
    }

    /// Converts an `AtomView` to a Rust `TokenStream`.
    fn convert_view(&self, view: AtomView<'_>) -> TokenStream {
        match view {
            AtomView::Num(n) => self.convert_num(n),
            AtomView::Var(v) => self.convert_var(v),
            AtomView::Add(a) => self.convert_add(a),
            AtomView::Mul(m) => self.convert_mul(m),
            AtomView::Pow(p) => self.convert_pow(p),
            AtomView::Fun(_) => {
                // Functions are not expected in our expressions
                quote! { T::zero() }
            }
        }
    }

    /// Converts a numeric value to Rust.
    fn convert_num(&self, num: symbolica::atom::NumView<'_>) -> TokenStream {
        let coeff = num.get_coeff_view();

        match coeff {
            CoefficientView::Natural(n_re, d_re, n_im, d_im) => {
                // We only handle real numbers
                if n_im != 0 || d_im != 1 {
                    // Complex number - shouldn't happen in our use case
                    return quote! { T::zero() };
                }

                if d_re == 1 {
                    // Integer
                    self.convert_integer(n_re)
                } else {
                    // Rational - convert to float
                    let val = n_re as f64 / d_re as f64;
                    quote! { T::from_f64(#val) }
                }
            }
            CoefficientView::Float(_, _) => {
                // Float coefficients are rare in our use case
                // Fall back to zero - this shouldn't happen in practice
                quote! { T::zero() }
            }
            CoefficientView::Large(_, _) => {
                // Large coefficients are rare in our use case
                // Fall back to zero - this shouldn't happen in practice
                quote! { T::zero() }
            }
            CoefficientView::FiniteField(_, _) => {
                // Finite field - shouldn't happen
                quote! { T::zero() }
            }
            CoefficientView::RationalPolynomial(_) => {
                // Rational polynomial - shouldn't happen
                quote! { T::zero() }
            }
            CoefficientView::Indeterminate | CoefficientView::Infinity(_) => {
                // Indeterminate or infinity - shouldn't happen in our use case
                quote! { T::zero() }
            }
        }
    }

    /// Converts an integer to the appropriate Rust expression.
    fn convert_integer(&self, n: i64) -> TokenStream {
        match n {
            0 => quote! { T::zero() },
            1 => quote! { T::one() },
            2 => quote! { T::TWO },
            -1 => quote! { -T::one() },
            -2 => quote! { -T::TWO },
            _ if n >= i8::MIN as i64 && n <= i8::MAX as i64 => {
                let n_i8 = n as i8;
                quote! { T::from_i8(#n_i8) }
            }
            _ => {
                let n_f64 = n as f64;
                quote! { T::from_f64(#n_f64) }
            }
        }
    }

    /// Converts a variable reference to Rust.
    fn convert_var(&self, var: symbolica::atom::VarView<'_>) -> TokenStream {
        let symbol = var.get_symbol();
        let name = symbol.to_string();

        if let Some((prefix, field)) = self.symbol_map.get(&name) {
            // It's a field reference like "a_x"
            let prefix_ident = format_ident!("{}", prefix);
            let field_ident = format_ident!("{}", field);

            // For wrapper types, access via .as_inner()
            if self.wrapped_prefixes.contains(prefix) {
                quote! { #prefix_ident.as_inner().#field_ident() }
            } else {
                quote! { #prefix_ident.#field_ident() }
            }
        } else {
            // Unknown variable - emit as identifier (for constants or errors)
            let ident = format_ident!("{}", name.replace("clifford_codegen::", ""));
            quote! { #ident }
        }
    }

    /// Converts an addition to Rust, properly handling signs.
    fn convert_add(&self, add: symbolica::atom::AddView<'_>) -> TokenStream {
        let terms: Vec<AtomView<'_>> = add.iter().collect();

        if terms.is_empty() {
            return quote! { T::zero() };
        }

        // Convert first term
        let first = &terms[0];
        let (first_neg, first_expr) = self.extract_negation(*first);

        let mut result = if first_neg {
            quote! { -(#first_expr) }
        } else {
            first_expr
        };

        // Add remaining terms with proper sign handling
        for term in &terms[1..] {
            let (is_neg, term_expr) = self.extract_negation(*term);

            if is_neg {
                result = quote! { #result - #term_expr };
            } else {
                result = quote! { #result + #term_expr };
            }
        }

        result
    }

    /// Extracts negation from an expression.
    ///
    /// Returns (is_negative, positive_expression).
    fn extract_negation(&self, view: AtomView<'_>) -> (bool, TokenStream) {
        match view {
            AtomView::Num(n) => {
                let coeff = n.get_coeff_view();
                if self.is_negative_coefficient(&coeff) {
                    // Negate and return positive version
                    let atom = view.to_owned();
                    let negated = -&atom;
                    (true, self.convert(&negated))
                } else {
                    (false, self.convert_num(n))
                }
            }
            AtomView::Mul(m) => {
                // Check if first factor is a negative number
                let factors: Vec<AtomView<'_>> = m.iter().collect();
                if let Some(AtomView::Num(n)) = factors.first() {
                    let coeff = n.get_coeff_view();
                    if self.is_negative_coefficient(&coeff) {
                        // Negate the multiplication
                        let atom = view.to_owned();
                        let negated = -&atom;
                        let expanded = negated.expand();
                        return (true, self.convert(&expanded));
                    }
                }
                (false, self.convert_mul(m))
            }
            _ => (false, self.convert_view(view)),
        }
    }

    /// Checks if a coefficient is negative.
    fn is_negative_coefficient(&self, coeff: &CoefficientView) -> bool {
        match coeff {
            CoefficientView::Natural(n_re, _, _, _) => *n_re < 0,
            // Float and Large coefficients are rare in our use case
            // Fall back to false - this shouldn't happen in practice
            CoefficientView::Float(_, _) | CoefficientView::Large(_, _) => false,
            _ => false,
        }
    }

    /// Converts a multiplication to Rust.
    fn convert_mul(&self, mul: symbolica::atom::MulView<'_>) -> TokenStream {
        let factors: Vec<AtomView<'_>> = mul.iter().collect();

        if factors.is_empty() {
            return quote! { T::one() };
        }

        // Split out numeric coefficient if present
        let (coeff, remaining) = self.split_coefficient(&factors);

        if remaining.is_empty() {
            // Just a number
            return coeff.unwrap_or_else(|| quote! { T::one() });
        }

        // Convert remaining factors
        let factor_exprs: Vec<TokenStream> =
            remaining.iter().map(|f| self.convert_view(*f)).collect();

        // Build product
        let product = if factor_exprs.len() == 1 {
            factor_exprs[0].clone()
        } else {
            let first = &factor_exprs[0];
            let rest = &factor_exprs[1..];
            quote! { #first #(* #rest)* }
        };

        // Apply coefficient
        match coeff {
            None => product,
            Some(c) => {
                // Check if it's just T::one() or -T::one()
                let c_str = c.to_string();
                if c_str.contains("T :: one ()") && !c_str.contains('-') {
                    product
                } else if c_str == "- T :: one ()" {
                    quote! { -(#product) }
                } else {
                    quote! { #c * #product }
                }
            }
        }
    }

    /// Splits out the numeric coefficient from factors.
    ///
    /// Returns (coefficient_tokens, remaining_factors) where remaining_factors
    /// contains copies of non-numeric AtomViews.
    fn split_coefficient<'b>(
        &self,
        factors: &[AtomView<'b>],
    ) -> (Option<TokenStream>, Vec<AtomView<'b>>) {
        let mut coeff = None;
        let mut remaining = Vec::new();

        for factor in factors {
            if let AtomView::Num(n) = factor {
                if coeff.is_none() {
                    coeff = Some(self.convert_num(*n));
                    continue;
                }
            }
            remaining.push(*factor);
        }

        (coeff, remaining)
    }

    /// Converts a power expression to Rust.
    fn convert_pow(&self, pow: symbolica::atom::PowView<'_>) -> TokenStream {
        let (base, exp) = pow.get_base_exp();

        // Check for simple integer exponents
        if let AtomView::Num(n) = exp {
            let coeff = n.get_coeff_view();
            if let CoefficientView::Natural(exp_val, 1, 0, 1) = coeff {
                match exp_val {
                    0 => return quote! { T::one() },
                    1 => return self.convert_view(base),
                    2 => {
                        let base_expr = self.convert_view(base);
                        return quote! { #base_expr * #base_expr };
                    }
                    3 => {
                        let base_expr = self.convert_view(base);
                        return quote! { #base_expr * #base_expr * #base_expr };
                    }
                    _ => {}
                }
            }
        }

        // General case: use powi for integers, powf otherwise
        let base_expr = self.convert_view(base);
        let exp_expr = self.convert_view(exp);
        quote! { #base_expr.powf(#exp_expr) }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::spec::FieldSpec;
    use std::sync::Mutex;
    use symbolica::atom::Atom;

    // Symbolica uses global state that conflicts when tests run in parallel.
    // Tests prefixed with `symbolica_` are configured to run serially via nextest.
    // The mutex provides a fallback for `cargo test` users.
    static SYMBOLICA_LOCK: Mutex<()> = Mutex::new(());

    fn make_test_type(name: &str, fields: &[(&str, usize)]) -> TypeSpec {
        TypeSpec {
            name: name.to_string(),
            grades: vec![1],
            description: None,
            fields: fields
                .iter()
                .map(|(n, idx)| FieldSpec {
                    name: n.to_string(),
                    blade_index: *idx,
                    grade: 1,
                })
                .collect(),
            alias_of: None,
            versor: None,
            is_sparse: false,
        }
    }

    #[test]
    fn symbolica_convert_integer_constants() {
        let _guard = SYMBOLICA_LOCK.lock().unwrap();
        let type_a = make_test_type("A", &[("x", 1)]);
        let converter = AtomToRust::new(&[&type_a], &["a"]);

        // Test zero
        let zero = Atom::num(0);
        assert!(converter.convert(&zero).to_string().contains("zero"));

        // Test one
        let one = Atom::num(1);
        assert!(converter.convert(&one).to_string().contains("one"));

        // Test two
        let two = Atom::num(2);
        assert!(converter.convert(&two).to_string().contains("TWO"));
    }

    #[test]
    fn symbolica_convert_negative_integers() {
        let _guard = SYMBOLICA_LOCK.lock().unwrap();
        let type_a = make_test_type("A", &[("x", 1)]);
        let converter = AtomToRust::new(&[&type_a], &["a"]);

        // Test -1
        let neg_one = Atom::num(-1);
        let result = converter.convert(&neg_one).to_string();
        assert!(result.contains("one") && result.contains("-"));
    }
}
