//! Symbolic product computation.
//!
//! This module computes the symbolic output of product operations,
//! representing each output field as a symbolic expression of input fields.

use std::borrow::Cow;
use std::collections::HashMap;

use symbolica::atom::{Atom, DefaultNamespace};
use symbolica::parser::ParseSettings;

use crate::algebra::{Algebra, Blade, ProductTable};
use crate::spec::TypeSpec;

/// The kind of product to compute symbolically.
///
/// Product naming follows [Rigid Geometric Algebra](https://rigidgeometricalgebra.org/) conventions.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ProductKind {
    /// Geometric product (full product).
    Geometric,
    /// Geometric antiproduct (complement(complement(a) × complement(b))).
    /// Used for versor composition with antisandwich-based transformations.
    Antigeometric,
    /// Wedge product (∧, exterior, grade-raising).
    Wedge,
    /// Left contraction (a ⌋ b, grade gb - ga when ga <= gb).
    LeftContraction,
    /// Right contraction (a ⌊ b, grade ga - gb when gb <= ga).
    RightContraction,
    /// Antiwedge product (∨, regressive/meet).
    Antiwedge,
    /// Bulk contraction (a ∨ b★).
    BulkContraction,
    /// Weight contraction (a ∨ b☆).
    WeightContraction,
    /// Bulk expansion (a ∧ b★).
    BulkExpansion,
    /// Weight expansion (a ∧ b☆).
    WeightExpansion,
    /// Dot product (• metric inner, same-grade only, returns scalar).
    Dot,
    /// Antidot product (⊚ metric antiproduct inner, same-antigrade only, returns scalar).
    Antidot,
    /// Scalar product (grade-0 projection of geometric product).
    Scalar,
    /// Projection: target ∨ (self ∧ target☆).
    Project,
    /// Antiprojection: target ∧ (self ∨ target☆).
    Antiproject,
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
pub struct SymbolicProduct {
    /// The product table.
    table: ProductTable,
}

impl SymbolicProduct {
    /// Creates a new symbolic product computer.
    pub fn new(algebra: &Algebra) -> Self {
        let table = ProductTable::new(algebra);
        Self { table }
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
        let mut terms: Vec<Atom> = Vec::new();

        for field_a in &type_a.fields {
            for field_b in &type_b.fields {
                let a_blade = field_a.blade_index;
                let b_blade = field_b.blade_index;

                // Compute the product based on kind using ProductTable methods
                let (sign, result) = match kind {
                    ProductKind::Geometric => self.table.geometric(a_blade, b_blade),
                    ProductKind::Antigeometric => self.table.antiproduct(a_blade, b_blade),
                    ProductKind::Wedge => self.table.exterior(a_blade, b_blade),
                    ProductKind::LeftContraction => self.table.left_contraction(a_blade, b_blade),
                    ProductKind::RightContraction => self.table.right_contraction(a_blade, b_blade),
                    ProductKind::Antiwedge => self.table.regressive(a_blade, b_blade),
                    ProductKind::BulkContraction => self.table.bulk_contraction(a_blade, b_blade),
                    ProductKind::WeightContraction => {
                        self.table.weight_contraction(a_blade, b_blade)
                    }
                    ProductKind::BulkExpansion => self.table.bulk_expansion(a_blade, b_blade),
                    ProductKind::WeightExpansion => self.table.weight_expansion(a_blade, b_blade),
                    ProductKind::Dot => self.table.dot(a_blade, b_blade),
                    ProductKind::Antidot => self.table.antidot(a_blade, b_blade),
                    ProductKind::Scalar => {
                        // Scalar product: grade-0 projection of geometric product
                        let (s, r) = self.table.geometric(a_blade, b_blade);
                        // Only include if result is grade 0 (scalar blade index = 0)
                        let result_grade = Blade::from_index(r).grade();
                        if result_grade == 0 { (s, r) } else { (0, 0) }
                    }
                    ProductKind::Project => self.table.project(a_blade, b_blade),
                    ProductKind::Antiproject => self.table.antiproject(a_blade, b_blade),
                };

                if result != result_blade || sign == 0 {
                    continue;
                }

                let a_sym = a_symbols.get(&field_a.name).unwrap();
                let b_sym = b_symbols.get(&field_b.name).unwrap();

                // Create term: sign * a_field * b_field
                let product = a_sym * b_sym;
                let term = if sign > 0 { product } else { -product };
                terms.push(term);
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
    use std::sync::Mutex;

    // Symbolica uses global state that conflicts when tests run in parallel.
    // Tests prefixed with `symbolica_` are configured to run serially via nextest.
    // The mutex provides a fallback for `cargo test` users.
    static SYMBOLICA_LOCK: Mutex<()> = Mutex::new(());

    #[test]
    fn symbolica_create_symbols_for_vector() {
        let _guard = SYMBOLICA_LOCK.lock().unwrap();
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
