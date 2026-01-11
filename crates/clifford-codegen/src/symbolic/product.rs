//! Symbolic product computation.
//!
//! This module computes the symbolic output of product operations,
//! representing each output field as a symbolic expression of input fields.

use std::borrow::Cow;
use std::collections::HashMap;

use symbolica::atom::{Atom, DefaultNamespace};
use symbolica::parser::ParseSettings;

use crate::algebra::{Algebra, Blade, ProductTable, left_contraction_grade, outer_grade};
use crate::spec::TypeSpec;

/// The kind of product to compute symbolically.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ProductKind {
    /// Geometric product (full product).
    Geometric,
    /// Outer product (wedge, grade-raising).
    Outer,
    /// Left contraction (inner product).
    LeftContraction,
}

/// A symbolic field expression.
#[derive(Debug, Clone)]
pub struct SymbolicField {
    /// Field name in the output type.
    pub name: String,
    /// Blade index this field represents.
    pub blade_index: usize,
    /// Symbolic expression for this field's value.
    pub expression: Atom,
}

/// Computes symbolic product outputs.
///
/// Given input types with symbolic field values, computes the symbolic
/// expressions for each output field.
pub struct SymbolicProduct<'a> {
    /// The algebra for blade computations.
    #[allow(dead_code)]
    algebra: &'a Algebra,
    /// The product table.
    table: ProductTable,
}

impl<'a> SymbolicProduct<'a> {
    /// Creates a new symbolic product computer.
    pub fn new(algebra: &'a Algebra) -> Self {
        let table = ProductTable::new(algebra);
        Self { algebra, table }
    }

    /// Creates symbolic variables for a type's fields.
    ///
    /// Returns a map from field name to symbolic atom.
    pub fn create_field_symbols(&self, ty: &TypeSpec, prefix: &str) -> HashMap<String, Atom> {
        ty.fields
            .iter()
            .map(|f| {
                let symbol_name = format!("{}_{}", prefix, f.name);
                let input = DefaultNamespace {
                    namespace: Cow::Borrowed(env!("CARGO_CRATE_NAME")),
                    data: symbol_name.as_str(),
                    file: Cow::Borrowed(file!()),
                    line: line!() as usize,
                };
                let atom = Atom::parse(input, ParseSettings::symbolica()).unwrap();
                (f.name.clone(), atom)
            })
            .collect()
    }

    /// Computes the symbolic output of a product.
    ///
    /// # Arguments
    ///
    /// * `type_a` - Left operand type
    /// * `type_b` - Right operand type
    /// * `output_type` - Output type
    /// * `kind` - Product kind (geometric, outer, etc.)
    /// * `a_symbols` - Symbolic values for type_a fields
    /// * `b_symbols` - Symbolic values for type_b fields
    ///
    /// # Returns
    ///
    /// Symbolic expressions for each output field.
    pub fn compute(
        &self,
        type_a: &TypeSpec,
        type_b: &TypeSpec,
        output_type: &TypeSpec,
        kind: ProductKind,
        a_symbols: &HashMap<String, Atom>,
        b_symbols: &HashMap<String, Atom>,
    ) -> Vec<SymbolicField> {
        output_type
            .fields
            .iter()
            .map(|output_field| {
                let expr = self.compute_field(
                    type_a,
                    type_b,
                    output_field.blade_index,
                    kind,
                    a_symbols,
                    b_symbols,
                );
                SymbolicField {
                    name: output_field.name.clone(),
                    blade_index: output_field.blade_index,
                    expression: expr,
                }
            })
            .collect()
    }

    /// Computes the symbolic expression for a single output field.
    fn compute_field(
        &self,
        type_a: &TypeSpec,
        type_b: &TypeSpec,
        result_blade: usize,
        kind: ProductKind,
        a_symbols: &HashMap<String, Atom>,
        b_symbols: &HashMap<String, Atom>,
    ) -> Atom {
        let dim = self.table.dim();
        let mut terms: Vec<Atom> = Vec::new();

        for field_a in &type_a.fields {
            for field_b in &type_b.fields {
                let a_blade = field_a.blade_index;
                let b_blade = field_b.blade_index;

                let (sign, result) = self.table.geometric(a_blade, b_blade);

                if result != result_blade || sign == 0 {
                    continue;
                }

                // Filter based on product kind
                let include = match kind {
                    ProductKind::Geometric => true,
                    ProductKind::Outer => {
                        let a_grade = Blade::from_index(a_blade).grade();
                        let b_grade = Blade::from_index(b_blade).grade();
                        let result_grade = Blade::from_index(result_blade).grade();
                        outer_grade(a_grade, b_grade, dim)
                            .map(|g| g == result_grade)
                            .unwrap_or(false)
                    }
                    ProductKind::LeftContraction => {
                        let a_grade = Blade::from_index(a_blade).grade();
                        let b_grade = Blade::from_index(b_blade).grade();
                        let result_grade = Blade::from_index(result_blade).grade();
                        left_contraction_grade(a_grade, b_grade)
                            .map(|g| g == result_grade)
                            .unwrap_or(false)
                    }
                };

                if include {
                    let a_sym = a_symbols.get(&field_a.name).unwrap();
                    let b_sym = b_symbols.get(&field_b.name).unwrap();

                    // Create term: sign * a_field * b_field
                    let product = a_sym * b_sym;
                    let term = if sign > 0 { product } else { -product };
                    terms.push(term);
                }
            }
        }

        if terms.is_empty() {
            // Return symbolic zero
            Atom::num(0)
        } else {
            // Sum all terms
            terms.into_iter().reduce(|acc, t| acc + t).unwrap()
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::spec::parse_spec;

    // Note: This test creates symbols which can cause issues with Symbolica's
    // global state when running in parallel with other tests. Run with
    // --test-threads=1 or individually if it fails.
    #[test]
    #[ignore = "Symbolica global state conflicts when running in parallel"]
    fn create_symbols_for_vector() {
        let spec = parse_spec(
            r#"
            [algebra]
            name = "test"

            [signature]
            positive = ["e1", "e2", "e3"]

            [types.Vector]
            grades = [1]
            fields = ["x", "y", "z"]
            "#,
        )
        .unwrap();

        let algebra = Algebra::euclidean(3);
        let product = SymbolicProduct::new(&algebra);
        let vector = spec.types.iter().find(|t| t.name == "Vector").unwrap();

        let symbols = product.create_field_symbols(vector, "a");

        assert!(symbols.contains_key("x"));
        assert!(symbols.contains_key("y"));
        assert!(symbols.contains_key("z"));
    }
}
