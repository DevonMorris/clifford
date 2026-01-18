//! Constraint derivation from algebra structure.
//!
//! This module derives constraint expressions by computing `u * involution(u)`
//! symbolically for a type and extracting the non-scalar terms. The geometric
//! constraint states that `u * involution(u)` must produce only a scalar for
//! valid geometric entities.
//!
//! The involution used depends on the algebra's configuration:
//! - `Reverse`: Most algebras (Euclidean, PGA) - versors satisfy `v * reverse(v) = scalar`
//! - `GradeInvolution`: Hyperbolic numbers - `z * involute(z) = a² - b²`
//! - `CliffordConjugate`: Complex numbers - `z * conjugate(z) = a² + b²`
//!
//! For example, for a rotor R = s + B (scalar + bivector), computing
//! `R * reverse(R)` produces scalar terms and potentially non-scalar terms.
//! The non-scalar terms must equal zero, giving us the constraint expression.

use std::collections::HashMap;

use symbolica::atom::Atom;

use crate::algebra::{
    Algebra, Blade, ProductTable, antireverse_sign, clifford_conjugate_sign, grade,
    grade_involution_sign, reverse_sign,
};
use crate::spec::{InvolutionKind, TypeSpec};

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
    /// Which involution to use for constraint derivation.
    involution: InvolutionKind,
}

impl<'a> ConstraintDeriver<'a> {
    /// Creates a new constraint deriver with the specified involution.
    pub fn new(algebra: &'a Algebra, involution: InvolutionKind) -> Self {
        let table = ProductTable::new(algebra);
        Self {
            algebra,
            table,
            involution,
        }
    }

    /// Derives the geometric constraint for a type.
    ///
    /// Computes `u * involution(u)` symbolically and returns the non-scalar
    /// terms that must equal zero. The involution used is determined by
    /// the algebra's configuration (reverse, grade involution, or Clifford conjugate).
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

        // Compute u * involution(u) for each output grade
        let mut zero_expressions = Vec::new();
        let mut constraint_grades = Vec::new();

        // Check all grades from 1 to dim (skip scalar grade 0)
        for output_grade in 1..=dim {
            let expr = self.compute_product_involution_at_grade(ty, output_grade, &symbols);

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

    /// Computes the expression for `u * involution(u)` at a specific output grade.
    fn compute_product_involution_at_grade(
        &self,
        ty: &TypeSpec,
        output_grade: usize,
        symbols: &HashMap<usize, Atom>,
    ) -> Atom {
        let mut terms: Vec<Atom> = Vec::new();

        // For each pair of fields (a, b), compute a * involution(b)
        // The involution of a blade B depends on the configured involution kind
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

                // Apply involution sign to b based on configured involution kind
                let inv_sign = match self.involution {
                    InvolutionKind::Reverse => reverse_sign(grade_b),
                    InvolutionKind::GradeInvolution => grade_involution_sign(grade_b),
                    InvolutionKind::CliffordConjugate => clifford_conjugate_sign(grade_b),
                };
                let total_sign = product_sign * inv_sign;

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
        self.compute_product_involution_at_grade(ty, 0, &symbols)
    }

    /// Derives the antiscalar norm expression (pseudoscalar part of `u ⊟ antireverse(u)`).
    ///
    /// This is useful for generating weight norm implementations in PGA.
    pub fn derive_antinorm_squared(&self, ty: &TypeSpec, field_prefix: &str) -> Atom {
        let dim = self.algebra.dim();
        let symbols = self.create_symbols(ty, field_prefix);
        self.compute_antiproduct_antireverse_at_grade(ty, dim, &symbols)
    }

    /// Derives the weight norm squared expression for PGA types.
    ///
    /// In PGA, the weight norm is the degenerate part of the element.
    /// For points, this is the w coordinate. For planes, it's the d coefficient.
    /// The weight norm squared is the sum of squares of all weight (degenerate) fields.
    ///
    /// Weight fields are those whose blade contains the degenerate basis vector (e0).
    pub fn derive_weight_norm_squared(&self, ty: &TypeSpec, field_prefix: &str) -> Atom {
        let symbols = self.create_symbols(ty, field_prefix);
        let (p, q, r) = self.algebra.signature();
        // Degenerate indices are p+q..p+q+r
        let degenerate_indices: Vec<_> = (p + q..p + q + r).collect();

        let mut terms: Vec<Atom> = Vec::new();

        for field in &ty.fields {
            // Check if this blade contains a degenerate basis vector
            let blade = field.blade_index;
            let is_weight = degenerate_indices
                .iter()
                .any(|&idx| (blade & (1 << idx)) != 0);

            if is_weight {
                if let Some(sym) = symbols.get(&blade) {
                    // Add sym²
                    terms.push(sym * sym);
                }
            }
        }

        if terms.is_empty() {
            Atom::num(0)
        } else {
            terms.into_iter().reduce(|acc, t| acc + t).unwrap()
        }
    }

    /// Derives the bulk norm squared expression for PGA types.
    ///
    /// In PGA, the bulk norm is the non-degenerate part of the element.
    /// The bulk norm squared is the sum of squares of all bulk (non-degenerate) fields.
    ///
    /// Bulk fields are those whose blade does NOT contain the degenerate basis vector (e0).
    pub fn derive_bulk_norm_squared(&self, ty: &TypeSpec, field_prefix: &str) -> Atom {
        let symbols = self.create_symbols(ty, field_prefix);
        let (p, q, r) = self.algebra.signature();
        // Degenerate indices are p+q..p+q+r
        let degenerate_indices: Vec<_> = (p + q..p + q + r).collect();

        let mut terms: Vec<Atom> = Vec::new();

        for field in &ty.fields {
            // Check if this blade does NOT contain any degenerate basis vectors
            let blade = field.blade_index;
            let is_bulk = degenerate_indices
                .iter()
                .all(|&idx| (blade & (1 << idx)) == 0);

            if is_bulk {
                if let Some(sym) = symbols.get(&blade) {
                    // Add sym²
                    terms.push(sym * sym);
                }
            }
        }

        if terms.is_empty() {
            Atom::num(0)
        } else {
            terms.into_iter().reduce(|acc, t| acc + t).unwrap()
        }
    }

    /// Derives individual weight component expressions for the Ideal wrapper.
    ///
    /// For Ideal elements (elements at infinity in PGA), all weight components
    /// must be zero. This returns a list of the symbolic expressions for each
    /// weight field, which can be used as constraints in Groebner simplification.
    pub fn derive_weight_components(&self, ty: &TypeSpec, field_prefix: &str) -> Vec<Atom> {
        let symbols = self.create_symbols(ty, field_prefix);
        let (p, q, r) = self.algebra.signature();
        // Degenerate indices are p+q..p+q+r
        let degenerate_indices: Vec<_> = (p + q..p + q + r).collect();

        let mut components = Vec::new();

        for field in &ty.fields {
            // Check if this blade contains a degenerate basis vector
            let blade = field.blade_index;
            let is_weight = degenerate_indices
                .iter()
                .any(|&idx| (blade & (1 << idx)) != 0);

            if is_weight {
                if let Some(sym) = symbols.get(&blade) {
                    components.push(sym.clone());
                }
            }
        }

        components
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::spec::parse_spec;

    #[test]
    fn symbolica_rotor_has_no_derived_constraints() {
        // Euclidean rotor: u * reverse(u) should be purely scalar
        // (no non-scalar constraints needed because cross-terms cancel)
        let spec = parse_spec(
            r#"
            [algebra]
            name = "test"
            complete = false

            [signature]
            positive = ["e1", "e2", "e3"]

            [types.Rotor]
            grades = [0, 2]
            field_map = [
                { name = "s", blade = "s" },
                { name = "xy", blade = "e12" },
                { name = "xz", blade = "e13" },
                { name = "yz", blade = "e23" }
            ]
            "#,
        )
        .unwrap();

        let algebra = Algebra::euclidean(3);
        let deriver = ConstraintDeriver::new(&algebra, InvolutionKind::Reverse);
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
    fn symbolica_vector_norm_squared() {
        let spec = parse_spec(
            r#"
            [algebra]
            name = "test"
            complete = false

            [signature]
            positive = ["e1", "e2", "e3"]

            [types.Vector]
            grades = [1]
            field_map = [
                { name = "x", blade = "e1" },
                { name = "y", blade = "e2" },
                { name = "z", blade = "e3" }
            ]
            "#,
        )
        .unwrap();

        let algebra = Algebra::euclidean(3);
        let deriver = ConstraintDeriver::new(&algebra, InvolutionKind::Reverse);
        let vector = spec.types.iter().find(|t| t.name == "Vector").unwrap();

        let norm_sq = deriver.derive_norm_squared(vector, "v");

        // Should be x² + y² + z²
        let expr_str = format!("{}", norm_sq);
        assert!(expr_str.contains("v_x"), "norm² should contain v_x");
        assert!(expr_str.contains("v_y"), "norm² should contain v_y");
        assert!(expr_str.contains("v_z"), "norm² should contain v_z");
    }

    #[test]
    fn symbolica_pga_motor_has_study_condition() {
        // PGA motor: grades [0, 2, 4] - should have Study condition constraint
        // The Study condition arises from the pseudoscalar term in motor * reverse(motor)
        let spec = parse_spec(
            r#"
            [algebra]
            name = "pga3"
            complete = false

            [signature]
            positive = ["e1", "e2", "e3"]
            zero = ["e0"]

            [types.Motor]
            grades = [0, 2, 4]
            field_map = [
                { name = "s", blade = "s" },
                { name = "e01", blade = "e14" },
                { name = "e02", blade = "e24" },
                { name = "e03", blade = "e34" },
                { name = "e12", blade = "e12" },
                { name = "e31", blade = "e13" },
                { name = "e23", blade = "e23" },
                { name = "e0123", blade = "e1234" }
            ]
            "#,
        )
        .unwrap();

        let algebra = Algebra::pga(3);
        let deriver = ConstraintDeriver::new(&algebra, InvolutionKind::Reverse);
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
    fn symbolica_pga_bivector_has_plucker_condition() {
        // PGA line (bivector): grade [2] - should have Plücker condition
        // In 4D PGA, bivectors don't automatically satisfy geometric constraint
        let spec = parse_spec(
            r#"
            [algebra]
            name = "pga3"
            complete = false

            [signature]
            positive = ["e1", "e2", "e3"]
            zero = ["e0"]

            [types.Line]
            grades = [2]
            field_map = [
                { name = "e01", blade = "e14" },
                { name = "e02", blade = "e24" },
                { name = "e03", blade = "e34" },
                { name = "e12", blade = "e12" },
                { name = "e31", blade = "e13" },
                { name = "e23", blade = "e23" }
            ]
            "#,
        )
        .unwrap();

        let algebra = Algebra::pga(3);
        let deriver = ConstraintDeriver::new(&algebra, InvolutionKind::Reverse);
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
    fn symbolica_euclidean_rotor_all_constraints_deduplicated() {
        // For Euclidean rotors, both geometric and antiproduct constraints
        // should yield no constraints (cross-terms cancel in both cases)
        let spec = parse_spec(
            r#"
            [algebra]
            name = "test"
            complete = false

            [signature]
            positive = ["e1", "e2", "e3"]

            [types.Rotor]
            grades = [0, 2]
            field_map = [
                { name = "s", blade = "s" },
                { name = "xy", blade = "e12" },
                { name = "xz", blade = "e13" },
                { name = "yz", blade = "e23" }
            ]
            "#,
        )
        .unwrap();

        let algebra = Algebra::euclidean(3);
        let deriver = ConstraintDeriver::new(&algebra, InvolutionKind::Reverse);
        let rotor = spec.types.iter().find(|t| t.name == "Rotor").unwrap();

        let all_constraints = deriver.derive_all_constraints(rotor, "u");

        // Euclidean rotors satisfy both constraints trivially
        assert!(
            all_constraints.is_none(),
            "Euclidean rotor should have no combined constraints"
        );
    }

    #[test]
    fn symbolica_pga_motor_all_constraints() {
        // For PGA motors, derive_all_constraints should combine geometric
        // and antiproduct constraints with deduplication
        let spec = parse_spec(
            r#"
            [algebra]
            name = "pga3"
            complete = false

            [signature]
            positive = ["e1", "e2", "e3"]
            zero = ["e0"]

            [types.Motor]
            grades = [0, 2, 4]
            field_map = [
                { name = "s", blade = "s" },
                { name = "e01", blade = "e14" },
                { name = "e02", blade = "e24" },
                { name = "e03", blade = "e34" },
                { name = "e12", blade = "e12" },
                { name = "e31", blade = "e13" },
                { name = "e23", blade = "e23" },
                { name = "e0123", blade = "e1234" }
            ]
            "#,
        )
        .unwrap();

        let algebra = Algebra::pga(3);
        let deriver = ConstraintDeriver::new(&algebra, InvolutionKind::Reverse);
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
