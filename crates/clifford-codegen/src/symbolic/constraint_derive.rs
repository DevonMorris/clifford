//! Constraint derivation from algebra structure.
//!
//! This module derives constraint expressions by computing `u * reverse(u)`
//! symbolically for a type and extracting the non-scalar terms. The geometric
//! constraint states that `u * reverse(u)` must produce only a scalar for
//! valid geometric entities.
//!
//! For example, for a rotor R = s + B (scalar + bivector), computing
//! `R * reverse(R)` produces scalar terms and potentially non-scalar terms.
//! The non-scalar terms must equal zero, giving us the constraint expression.

use std::collections::HashMap;

use symbolica::atom::Atom;

use crate::algebra::{Algebra, Blade, ProductTable, antireverse_sign, grade, reverse_sign};
use crate::spec::TypeSpec;

/// Result of deriving constraints for a type.
#[derive(Debug, Clone)]
pub struct DerivedConstraint {
    /// Name describing the constraint (e.g., "geometric", "study").
    pub name: String,

    /// Symbolic expressions that must equal zero.
    ///
    /// Each expression corresponds to a non-scalar grade in `u * reverse(u)`.
    /// All must be zero for the type to satisfy the geometric constraint.
    pub zero_expressions: Vec<Atom>,

    /// Grades of the non-scalar terms (for documentation).
    pub constraint_grades: Vec<usize>,
}

/// Derives constraint expressions from algebra structure.
pub struct ConstraintDeriver<'a> {
    /// The algebra for computations.
    algebra: &'a Algebra,
    /// Product table for efficient computation.
    table: ProductTable,
}

impl<'a> ConstraintDeriver<'a> {
    /// Creates a new constraint deriver.
    pub fn new(algebra: &'a Algebra) -> Self {
        let table = ProductTable::new(algebra);
        Self { algebra, table }
    }

    /// Derives the geometric constraint for a type.
    ///
    /// Computes `u * reverse(u)` symbolically and returns the non-scalar
    /// terms that must equal zero.
    ///
    /// # Returns
    ///
    /// - `Some(constraint)` if the type has non-trivial constraints
    /// - `None` if `u * reverse(u)` is purely scalar (no constraints needed)
    pub fn derive_geometric_constraint(
        &self,
        ty: &TypeSpec,
        field_prefix: &str,
    ) -> Option<DerivedConstraint> {
        let dim = self.algebra.dim();

        // Create symbolic variables for fields
        let symbols = self.create_symbols(ty, field_prefix);

        // Compute u * reverse(u) for each output grade
        let mut zero_expressions = Vec::new();
        let mut constraint_grades = Vec::new();

        // Check all grades from 1 to dim (skip scalar grade 0)
        for output_grade in 1..=dim {
            let expr = self.compute_product_reverse_at_grade(ty, output_grade, &symbols);

            // Check if the expression is non-trivial (not just zero)
            if !self.is_zero(&expr) {
                zero_expressions.push(expr);
                constraint_grades.push(output_grade);
            }
        }

        if zero_expressions.is_empty() {
            None
        } else {
            Some(DerivedConstraint {
                name: "geometric".to_string(),
                zero_expressions,
                constraint_grades,
            })
        }
    }

    /// Creates symbolic variables for a type's fields.
    fn create_symbols(&self, ty: &TypeSpec, prefix: &str) -> HashMap<usize, Atom> {
        use std::borrow::Cow;
        use symbolica::atom::DefaultNamespace;
        use symbolica::parser::ParseSettings;

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
                (f.blade_index, atom)
            })
            .collect()
    }

    /// Computes the expression for `u * reverse(u)` at a specific output grade.
    fn compute_product_reverse_at_grade(
        &self,
        ty: &TypeSpec,
        output_grade: usize,
        symbols: &HashMap<usize, Atom>,
    ) -> Atom {
        let mut terms: Vec<Atom> = Vec::new();

        // For each pair of fields (a, b), compute a * reverse(b)
        // The reverse of a blade B is: reverse_sign(grade(B)) * B
        for field_a in &ty.fields {
            for field_b in &ty.fields {
                let blade_a = field_a.blade_index;
                let blade_b = field_b.blade_index;
                let grade_b = Blade::from_index(blade_b).grade();

                // Get the product a * b
                let (product_sign, result_blade) = self.table.geometric(blade_a, blade_b);

                // Check if result is at the target grade
                let result_grade = grade(result_blade);
                if result_grade != output_grade || product_sign == 0 {
                    continue;
                }

                // Apply reverse sign to b
                let rev_sign = reverse_sign(grade_b);
                let total_sign = product_sign * rev_sign;

                // Get symbolic values
                let sym_a = symbols.get(&blade_a).unwrap();
                let sym_b = symbols.get(&blade_b).unwrap();

                // Create term: total_sign * a * b
                let product = sym_a * sym_b;
                let term = if total_sign > 0 { product } else { -product };
                terms.push(term);
            }
        }

        if terms.is_empty() {
            Atom::num(0)
        } else {
            terms.into_iter().reduce(|acc, t| acc + t).unwrap()
        }
    }

    /// Checks if an expression is zero.
    fn is_zero(&self, expr: &Atom) -> bool {
        // Check if the expression is the number 0
        *expr == Atom::num(0)
    }

    /// Derives the antiproduct constraint for a type.
    ///
    /// Computes `u ⊟ antireverse(u)` symbolically and returns the non-pseudoscalar
    /// terms that must equal zero.
    ///
    /// The antiproduct is defined as: `a ⊟ b = dual(dual(a) * dual(b))`
    /// where dual(x) = x * I⁻¹ (right complement with pseudoscalar).
    ///
    /// # Returns
    ///
    /// - `Some(constraint)` if the type has non-trivial antiproduct constraints
    /// - `None` if `u ⊟ antireverse(u)` is purely pseudoscalar (no constraints needed)
    pub fn derive_antiproduct_constraint(
        &self,
        ty: &TypeSpec,
        field_prefix: &str,
    ) -> Option<DerivedConstraint> {
        let dim = self.algebra.dim();

        // Create symbolic variables for fields
        let symbols = self.create_symbols(ty, field_prefix);

        // Compute u ⊟ antireverse(u) for each output grade
        let mut zero_expressions = Vec::new();
        let mut constraint_grades = Vec::new();

        // Check all grades from 0 to dim-1 (skip pseudoscalar grade dim)
        for output_grade in 0..dim {
            let expr = self.compute_antiproduct_antireverse_at_grade(ty, output_grade, &symbols);

            // Check if the expression is non-trivial (not just zero)
            if !self.is_zero(&expr) {
                zero_expressions.push(expr);
                constraint_grades.push(output_grade);
            }
        }

        if zero_expressions.is_empty() {
            None
        } else {
            Some(DerivedConstraint {
                name: "antiproduct".to_string(),
                zero_expressions,
                constraint_grades,
            })
        }
    }

    /// Derives all constraints for a type (both geometric and antiproduct).
    ///
    /// Returns the combined set of constraint expressions that must equal zero.
    /// Symbolically equivalent expressions are deduplicated using Symbolica.
    pub fn derive_all_constraints(
        &self,
        ty: &TypeSpec,
        field_prefix: &str,
    ) -> Option<DerivedConstraint> {
        let geometric = self.derive_geometric_constraint(ty, field_prefix);
        let antiproduct = self.derive_antiproduct_constraint(ty, field_prefix);

        match (geometric, antiproduct) {
            (None, None) => None,
            (Some(c), None) | (None, Some(c)) => Some(c),
            (Some(geo), Some(anti)) => {
                // Combine both constraints, deduplicating equivalent expressions
                let mut all_exprs = geo.zero_expressions;
                let mut all_grades = geo.constraint_grades;

                // Add antiproduct expressions if not symbolically equivalent to existing ones
                for (expr, grade) in anti
                    .zero_expressions
                    .into_iter()
                    .zip(anti.constraint_grades)
                {
                    if !self.is_equivalent_to_any(&expr, &all_exprs) {
                        all_exprs.push(expr);
                        all_grades.push(grade);
                    }
                }

                if all_exprs.is_empty() {
                    None
                } else {
                    Some(DerivedConstraint {
                        name: "combined".to_string(),
                        zero_expressions: all_exprs,
                        constraint_grades: all_grades,
                    })
                }
            }
        }
    }

    /// Checks if an expression is symbolically equivalent to any in a list.
    ///
    /// Two expressions are equivalent if their difference simplifies to zero.
    fn is_equivalent_to_any(&self, expr: &Atom, existing: &[Atom]) -> bool {
        use crate::symbolic::ExpressionSimplifier;
        let simplifier = ExpressionSimplifier::new();

        for other in existing {
            let diff = expr - other;
            let simplified = simplifier.simplify(&diff);
            if self.is_zero(&simplified) {
                return true;
            }
            // Also check if they're negatives of each other (same constraint)
            let neg_diff = expr + other;
            let neg_simplified = simplifier.simplify(&neg_diff);
            if self.is_zero(&neg_simplified) {
                return true;
            }
        }
        false
    }

    /// Computes the expression for `u ⊟ antireverse(u)` at a specific output grade.
    ///
    /// The antiproduct is: `a ⊟ b = dual(dual(a) * dual(b))`
    fn compute_antiproduct_antireverse_at_grade(
        &self,
        ty: &TypeSpec,
        output_grade: usize,
        symbols: &HashMap<usize, Atom>,
    ) -> Atom {
        let dim = self.algebra.dim();
        let pseudoscalar = (1 << dim) - 1; // All bits set
        let mut terms: Vec<Atom> = Vec::new();

        for field_a in &ty.fields {
            for field_b in &ty.fields {
                let blade_a = field_a.blade_index;
                let blade_b = field_b.blade_index;
                let grade_b = Blade::from_index(blade_b).grade();

                // Compute duals by XOR with pseudoscalar
                let dual_a = pseudoscalar ^ blade_a;
                let dual_b = pseudoscalar ^ blade_b;

                // Get the product of duals
                let (product_sign, dual_result) = self.table.geometric(dual_a, dual_b);

                // Convert back from dual (the antiproduct result)
                let result_blade = pseudoscalar ^ dual_result;

                // Check if result is at the target grade
                let result_grade = grade(result_blade);
                if result_grade != output_grade || product_sign == 0 {
                    continue;
                }

                // Apply antireverse sign to b (reverse sign of dual grade)
                let antirev_sign = antireverse_sign(grade_b, dim);
                let total_sign = product_sign * antirev_sign;

                // Get symbolic values
                let sym_a = symbols.get(&blade_a).unwrap();
                let sym_b = symbols.get(&blade_b).unwrap();

                // Create term: total_sign * a * b
                let product = sym_a * sym_b;
                let term = if total_sign > 0 { product } else { -product };
                terms.push(term);
            }
        }

        if terms.is_empty() {
            Atom::num(0)
        } else {
            terms.into_iter().reduce(|acc, t| acc + t).unwrap()
        }
    }

    /// Derives the scalar norm expression (scalar part of `u * reverse(u)`).
    ///
    /// This is useful for generating `norm_squared()` implementations.
    pub fn derive_norm_squared(&self, ty: &TypeSpec, field_prefix: &str) -> Atom {
        let symbols = self.create_symbols(ty, field_prefix);
        self.compute_product_reverse_at_grade(ty, 0, &symbols)
    }

    /// Derives the antiscalar norm expression (pseudoscalar part of `u ⊟ antireverse(u)`).
    ///
    /// This is useful for generating weight norm implementations in PGA.
    pub fn derive_antinorm_squared(&self, ty: &TypeSpec, field_prefix: &str) -> Atom {
        let dim = self.algebra.dim();
        let symbols = self.create_symbols(ty, field_prefix);
        self.compute_antiproduct_antireverse_at_grade(ty, dim, &symbols)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::spec::parse_spec;

    #[test]
    fn rotor_has_no_derived_constraints() {
        // Euclidean rotor: u * reverse(u) should be purely scalar
        // (no non-scalar constraints needed because cross-terms cancel)
        let spec = parse_spec(
            r#"
            [algebra]
            name = "test"

            [signature]
            positive = ["e1", "e2", "e3"]

            [types.Rotor]
            grades = [0, 2]
            fields = ["s", "xy", "xz", "yz"]
            "#,
        )
        .unwrap();

        let algebra = Algebra::euclidean(3);
        let deriver = ConstraintDeriver::new(&algebra);
        let rotor = spec.types.iter().find(|t| t.name == "Rotor").unwrap();

        let constraint = deriver.derive_geometric_constraint(rotor, "u");

        // For Euclidean rotors, cross-terms cancel due to anticommutativity
        // so there should be no derived constraint
        assert!(
            constraint.is_none(),
            "Euclidean rotor should have no derived geometric constraint"
        );
    }

    #[test]
    fn vector_norm_squared() {
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
        let deriver = ConstraintDeriver::new(&algebra);
        let vector = spec.types.iter().find(|t| t.name == "Vector").unwrap();

        let norm_sq = deriver.derive_norm_squared(vector, "v");

        // Should be x² + y² + z²
        let expr_str = format!("{}", norm_sq);
        assert!(expr_str.contains("v_x"), "norm² should contain v_x");
        assert!(expr_str.contains("v_y"), "norm² should contain v_y");
        assert!(expr_str.contains("v_z"), "norm² should contain v_z");
    }

    #[test]
    fn pga_motor_has_study_condition() {
        // PGA motor: grades [0, 2, 4] - should have Study condition constraint
        // The Study condition arises from the pseudoscalar term in motor * reverse(motor)
        let spec = parse_spec(
            r#"
            [algebra]
            name = "pga3"

            [signature]
            positive = ["e1", "e2", "e3"]
            zero = ["e0"]

            [types.Motor]
            grades = [0, 2, 4]
            fields = ["s", "e01", "e02", "e03", "e12", "e31", "e23", "e0123"]
            "#,
        )
        .unwrap();

        let algebra = Algebra::pga(3);
        let deriver = ConstraintDeriver::new(&algebra);
        let motor = spec.types.iter().find(|t| t.name == "Motor").unwrap();

        let constraint = deriver.derive_geometric_constraint(motor, "m");

        // PGA motors should have a derived constraint (the Study condition)
        assert!(
            constraint.is_some(),
            "PGA motor should have derived geometric constraint (Study condition)"
        );

        let constraint = constraint.unwrap();
        // The Study condition is on the pseudoscalar (grade 4) term
        assert!(
            constraint.constraint_grades.contains(&4),
            "Study condition should involve grade 4 (pseudoscalar)"
        );
    }

    #[test]
    fn pga_bivector_has_plucker_condition() {
        // PGA line (bivector): grade [2] - should have Plücker condition
        // In 4D PGA, bivectors don't automatically satisfy geometric constraint
        let spec = parse_spec(
            r#"
            [algebra]
            name = "pga3"

            [signature]
            positive = ["e1", "e2", "e3"]
            zero = ["e0"]

            [types.Line]
            grades = [2]
            fields = ["e01", "e02", "e03", "e12", "e31", "e23"]
            "#,
        )
        .unwrap();

        let algebra = Algebra::pga(3);
        let deriver = ConstraintDeriver::new(&algebra);
        let line = spec.types.iter().find(|t| t.name == "Line").unwrap();

        let constraint = deriver.derive_geometric_constraint(line, "l");

        // PGA bivectors should have a derived constraint (Plücker condition)
        // because disjoint bivectors commute in 4D
        assert!(
            constraint.is_some(),
            "PGA line should have derived geometric constraint (Plücker condition)"
        );
    }

    #[test]
    fn euclidean_rotor_all_constraints_deduplicated() {
        // For Euclidean rotors, both geometric and antiproduct constraints
        // should yield no constraints (cross-terms cancel in both cases)
        let spec = parse_spec(
            r#"
            [algebra]
            name = "test"

            [signature]
            positive = ["e1", "e2", "e3"]

            [types.Rotor]
            grades = [0, 2]
            fields = ["s", "xy", "xz", "yz"]
            "#,
        )
        .unwrap();

        let algebra = Algebra::euclidean(3);
        let deriver = ConstraintDeriver::new(&algebra);
        let rotor = spec.types.iter().find(|t| t.name == "Rotor").unwrap();

        let all_constraints = deriver.derive_all_constraints(rotor, "u");

        // Euclidean rotors satisfy both constraints trivially
        assert!(
            all_constraints.is_none(),
            "Euclidean rotor should have no combined constraints"
        );
    }

    #[test]
    fn pga_motor_all_constraints() {
        // For PGA motors, derive_all_constraints should combine geometric
        // and antiproduct constraints with deduplication
        let spec = parse_spec(
            r#"
            [algebra]
            name = "pga3"

            [signature]
            positive = ["e1", "e2", "e3"]
            zero = ["e0"]

            [types.Motor]
            grades = [0, 2, 4]
            fields = ["s", "e01", "e02", "e03", "e12", "e31", "e23", "e0123"]
            "#,
        )
        .unwrap();

        let algebra = Algebra::pga(3);
        let deriver = ConstraintDeriver::new(&algebra);
        let motor = spec.types.iter().find(|t| t.name == "Motor").unwrap();

        let all_constraints = deriver.derive_all_constraints(motor, "m");

        // PGA motors should have combined constraints
        assert!(
            all_constraints.is_some(),
            "PGA motor should have combined constraints"
        );

        let constraint = all_constraints.unwrap();
        // Print constraint info for debugging
        eprintln!(
            "PGA Motor combined constraint grades: {:?}",
            constraint.constraint_grades
        );
        eprintln!(
            "Number of constraint expressions: {}",
            constraint.zero_expressions.len()
        );
    }
}
