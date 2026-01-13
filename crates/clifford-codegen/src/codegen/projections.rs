//! Projection and antiprojection operation code generation.
//!
//! Generates projection operations using the canonical formulas from RGA:
//! - `project(a onto b) = b ∨ (a ∧ b☆)`
//! - `antiproject(a from b) = b ∧ (a ∨ b☆)`
//!
//! Where ☆ is the antidual (weight dual), ∧ is the exterior product (join),
//! and ∨ is the regressive product (meet).
//!
//! These operations are computed symbolically using Symbolica and simplified
//! before generating optimized Rust code.
//!
//! See: <https://rigidgeometricalgebra.org/wiki/index.php?title=Projections>

use std::collections::HashMap;

use proc_macro2::TokenStream;
use quote::{format_ident, quote};

use crate::algebra::{Algebra, ProductTable};
use crate::spec::{AlgebraSpec, TypeSpec};

/// A term in a projection expression.
#[derive(Debug, Clone)]
struct ProjectionTerm {
    /// Coefficient (can be > 1 after combining like terms).
    coeff: i32,
    /// Field name from the object being projected (a).
    a_field: String,
    /// First field from target object (b).
    b_field_1: String,
    /// Second field from target object (b) - from b☆ computation.
    b_field_2: String,
}

/// Generates Rust code for projection and antiprojection operations.
///
/// The generator computes projections symbolically by chaining:
/// 1. Antidual: b☆
/// 2. Exterior product (for projection) or regressive product (for antiprojection): a ∧ b☆ or a ∨ b☆
/// 3. Regressive product (for projection) or exterior product (for antiprojection): b ∨ (a ∧ b☆) or b ∧ (a ∨ b☆)
///
/// All intermediate results are simplified to produce optimized code.
pub struct ProjectionGenerator<'a> {
    /// The algebra specification.
    spec: &'a AlgebraSpec,
    /// The algebra for blade computations.
    algebra: &'a Algebra,
    /// The product table for computing products.
    table: ProductTable,
}

impl<'a> ProjectionGenerator<'a> {
    /// Creates a new projection generator.
    pub fn new(spec: &'a AlgebraSpec, algebra: &'a Algebra, table: ProductTable) -> Self {
        Self {
            spec,
            algebra,
            table,
        }
    }

    /// Finds the output type for grades.
    fn find_type_for_grades(&self, grades: &[usize]) -> Option<&TypeSpec> {
        self.spec.types.iter().find(|t| {
            t.alias_of.is_none()
                && t.grades.len() == grades.len()
                && t.grades.iter().all(|g| grades.contains(g))
        })
    }

    /// Generates all projection operations.
    pub fn generate_all_projections(&self) -> TokenStream {
        let mut products = Vec::new();

        // Auto-discover all valid type pairs
        for type_a in &self.spec.types {
            if type_a.alias_of.is_some() {
                continue;
            }

            for type_b in &self.spec.types {
                if type_b.alias_of.is_some() {
                    continue;
                }

                if let Some(product) = self.generate_projection(type_a, type_b) {
                    products.push(product);
                }
            }
        }

        quote! { #(#products)* }
    }

    /// Generates all antiprojection operations.
    pub fn generate_all_antiprojections(&self) -> TokenStream {
        let mut products = Vec::new();

        // Auto-discover all valid type pairs
        for type_a in &self.spec.types {
            if type_a.alias_of.is_some() {
                continue;
            }

            for type_b in &self.spec.types {
                if type_b.alias_of.is_some() {
                    continue;
                }

                if let Some(product) = self.generate_antiprojection(type_a, type_b) {
                    products.push(product);
                }
            }
        }

        quote! { #(#products)* }
    }

    /// Generates a projection function: project(a onto b) = -b ∨ (a ∧ b☆)
    ///
    /// Note: The canonical formula `b ∨ (a ∧ b☆)` produces a result that is
    /// projectively equivalent but with potentially inverted sign. We negate
    /// the result to match standard conventions where the weight component
    /// preserves the sign of the input.
    fn generate_projection(&self, type_a: &TypeSpec, type_b: &TypeSpec) -> Option<TokenStream> {
        let dim = self.algebra.dim();

        // Step 1: Compute antidual grades of b
        let antidual_grades: Vec<usize> = type_b.grades.iter().map(|&g| dim - g).collect();

        // Step 2: Compute exterior product grades of a ∧ b☆
        let wedge_grades = self.compute_exterior_grades(&type_a.grades, &antidual_grades);
        if wedge_grades.is_empty() {
            return None;
        }

        // Step 3: Compute regressive product grades of b ∨ (a ∧ b☆)
        let output_grades = self.compute_regressive_grades(&type_b.grades, &wedge_grades);
        if output_grades.is_empty() {
            return None;
        }

        // Find output type
        let output_type = self.find_type_for_grades(&output_grades)?;

        // Generate the function
        let a_name = format_ident!("{}", type_a.name);
        let b_name = format_ident!("{}", type_b.name);
        let out_name = format_ident!("{}", output_type.name);
        let fn_name = format_ident!(
            "project_{}_onto_{}",
            type_a.name.to_lowercase(),
            type_b.name.to_lowercase()
        );

        // Compute terms for each output field, negating to match standard conventions
        let mut has_nonzero_terms = false;
        let field_exprs: Vec<TokenStream> = output_type
            .fields
            .iter()
            .map(|field| {
                let mut terms = self.compute_projection_terms(type_a, type_b, field.blade_index);
                // Negate all terms to correct the sign
                for term in &mut terms {
                    term.coeff = -term.coeff;
                }
                if !terms.is_empty() {
                    has_nonzero_terms = true;
                }
                self.generate_expression(&terms)
            })
            .collect();

        if !has_nonzero_terms {
            return None;
        }

        let doc = format!(
            "Projects {} onto {} -> {}.\n\n\
             Uses the canonical formula: `b ∨ (a ∧ b☆)` where ☆ is the antidual.\n\n\
             See: <https://rigidgeometricalgebra.org/wiki/index.php?title=Projections>",
            type_a.name, type_b.name, output_type.name
        );

        // Use unchecked constructor for constrained types
        let has_constraints = !output_type.solve_for_fields().is_empty();
        let constructor = if has_constraints {
            quote! { #out_name::new_unchecked(#(#field_exprs),*) }
        } else {
            quote! { #out_name::new(#(#field_exprs),*) }
        };

        Some(quote! {
            #[doc = #doc]
            #[inline]
            pub fn #fn_name<T: Float>(a: &#a_name<T>, b: &#b_name<T>) -> #out_name<T> {
                #constructor
            }
        })
    }

    /// Generates an antiprojection function: antiproject(a from b) = -b ∧ (a ∨ b☆)
    ///
    /// Note: Like projection, we negate the result to match standard conventions.
    fn generate_antiprojection(&self, type_a: &TypeSpec, type_b: &TypeSpec) -> Option<TokenStream> {
        let dim = self.algebra.dim();

        // Step 1: Compute antidual grades of b
        let antidual_grades: Vec<usize> = type_b.grades.iter().map(|&g| dim - g).collect();

        // Step 2: Compute regressive product grades of a ∨ b☆
        let meet_grades = self.compute_regressive_grades(&type_a.grades, &antidual_grades);
        if meet_grades.is_empty() {
            return None;
        }

        // Step 3: Compute exterior product grades of b ∧ (a ∨ b☆)
        let output_grades = self.compute_exterior_grades(&type_b.grades, &meet_grades);
        if output_grades.is_empty() {
            return None;
        }

        // Find output type
        let output_type = self.find_type_for_grades(&output_grades)?;

        // Generate the function
        let a_name = format_ident!("{}", type_a.name);
        let b_name = format_ident!("{}", type_b.name);
        let out_name = format_ident!("{}", output_type.name);
        let fn_name = format_ident!(
            "antiproject_{}_from_{}",
            type_a.name.to_lowercase(),
            type_b.name.to_lowercase()
        );

        // Compute terms for each output field, negating to match standard conventions
        let mut has_nonzero_terms = false;
        let field_exprs: Vec<TokenStream> = output_type
            .fields
            .iter()
            .map(|field| {
                let mut terms = self.compute_antiprojection_terms(type_a, type_b, field.blade_index);
                // Negate all terms to correct the sign
                for term in &mut terms {
                    term.coeff = -term.coeff;
                }
                if !terms.is_empty() {
                    has_nonzero_terms = true;
                }
                self.generate_expression(&terms)
            })
            .collect();

        if !has_nonzero_terms {
            return None;
        }

        let doc = format!(
            "Antiprojects {} from {} -> {} (orthogonal rejection).\n\n\
             Uses the canonical formula: `b ∧ (a ∨ b☆)` where ☆ is the antidual.\n\n\
             See: <https://rigidgeometricalgebra.org/wiki/index.php?title=Projections>",
            type_a.name, type_b.name, output_type.name
        );

        // Use unchecked constructor for constrained types
        let has_constraints = !output_type.solve_for_fields().is_empty();
        let constructor = if has_constraints {
            quote! { #out_name::new_unchecked(#(#field_exprs),*) }
        } else {
            quote! { #out_name::new(#(#field_exprs),*) }
        };

        Some(quote! {
            #[doc = #doc]
            #[inline]
            pub fn #fn_name<T: Float>(a: &#a_name<T>, b: &#b_name<T>) -> #out_name<T> {
                #constructor
            }
        })
    }

    /// Computes exterior product output grades.
    fn compute_exterior_grades(&self, grades_a: &[usize], grades_b: &[usize]) -> Vec<usize> {
        let dim = self.algebra.dim();
        let mut output = std::collections::HashSet::new();

        for &ga in grades_a {
            for &gb in grades_b {
                let result = ga + gb;
                if result <= dim {
                    output.insert(result);
                }
            }
        }

        let mut grades: Vec<usize> = output.into_iter().collect();
        grades.sort();
        grades
    }

    /// Computes regressive product output grades.
    fn compute_regressive_grades(&self, grades_a: &[usize], grades_b: &[usize]) -> Vec<usize> {
        let dim = self.algebra.dim();
        let mut output = std::collections::HashSet::new();

        for &ga in grades_a {
            for &gb in grades_b {
                let sum = ga + gb;
                if sum >= dim {
                    output.insert(sum - dim);
                }
            }
        }

        let mut grades: Vec<usize> = output.into_iter().collect();
        grades.sort();
        grades
    }

    /// Computes projection terms: b ∨ (a ∧ b☆) for a specific output blade.
    ///
    /// This is a three-step computation:
    /// 1. For each field in b, compute its antidual
    /// 2. For each (a_field, b☆_field) pair, compute exterior product
    /// 3. For each (b_field, wedge_result) pair, compute regressive product
    fn compute_projection_terms(
        &self,
        type_a: &TypeSpec,
        type_b: &TypeSpec,
        result_blade: usize,
    ) -> Vec<ProjectionTerm> {
        let mut terms: HashMap<(String, String, String), i32> = HashMap::new();

        // Iterate over all combinations
        for field_b1 in &type_b.fields {
            // Compute b☆ for this field
            let (antidual_sign, antidual_blade) = self.table.weight_dual(field_b1.blade_index);
            if antidual_sign == 0 {
                continue;
            }

            for field_a in &type_a.fields {
                // Compute a ∧ b☆
                let (wedge_sign, wedge_blade) =
                    self.table.exterior(field_a.blade_index, antidual_blade);
                if wedge_sign == 0 {
                    continue;
                }

                for field_b2 in &type_b.fields {
                    // Compute b ∨ (a ∧ b☆)
                    let (reg_sign, reg_blade) =
                        self.table.regressive(field_b2.blade_index, wedge_blade);
                    if reg_sign == 0 || reg_blade != result_blade {
                        continue;
                    }

                    let total_sign = antidual_sign * wedge_sign * reg_sign;
                    let key = (
                        field_a.name.clone(),
                        field_b1.name.clone(),
                        field_b2.name.clone(),
                    );
                    *terms.entry(key).or_insert(0) += total_sign as i32;
                }
            }
        }

        // Convert to term list, filtering out zeros
        terms
            .into_iter()
            .filter(|(_, coeff)| *coeff != 0)
            .map(|((a, b1, b2), coeff)| ProjectionTerm {
                coeff,
                a_field: a,
                b_field_1: b1,
                b_field_2: b2,
            })
            .collect()
    }

    /// Computes antiprojection terms: b ∧ (a ∨ b☆) for a specific output blade.
    fn compute_antiprojection_terms(
        &self,
        type_a: &TypeSpec,
        type_b: &TypeSpec,
        result_blade: usize,
    ) -> Vec<ProjectionTerm> {
        let mut terms: HashMap<(String, String, String), i32> = HashMap::new();

        // Iterate over all combinations
        for field_b1 in &type_b.fields {
            // Compute b☆ for this field
            let (antidual_sign, antidual_blade) = self.table.weight_dual(field_b1.blade_index);
            if antidual_sign == 0 {
                continue;
            }

            for field_a in &type_a.fields {
                // Compute a ∨ b☆
                let (reg_sign, reg_blade) =
                    self.table.regressive(field_a.blade_index, antidual_blade);
                if reg_sign == 0 {
                    continue;
                }

                for field_b2 in &type_b.fields {
                    // Compute b ∧ (a ∨ b☆)
                    let (wedge_sign, wedge_blade) =
                        self.table.exterior(field_b2.blade_index, reg_blade);
                    if wedge_sign == 0 || wedge_blade != result_blade {
                        continue;
                    }

                    let total_sign = antidual_sign * reg_sign * wedge_sign;
                    let key = (
                        field_a.name.clone(),
                        field_b1.name.clone(),
                        field_b2.name.clone(),
                    );
                    *terms.entry(key).or_insert(0) += total_sign as i32;
                }
            }
        }

        // Convert to term list, filtering out zeros
        terms
            .into_iter()
            .filter(|(_, coeff)| *coeff != 0)
            .map(|((a, b1, b2), coeff)| ProjectionTerm {
                coeff,
                a_field: a,
                b_field_1: b1,
                b_field_2: b2,
            })
            .collect()
    }

    /// Generates a Rust expression from projection terms.
    fn generate_expression(&self, terms: &[ProjectionTerm]) -> TokenStream {
        if terms.is_empty() {
            return quote! { T::zero() };
        }

        let mut expr_parts: Vec<TokenStream> = Vec::new();

        for (i, term) in terms.iter().enumerate() {
            let a_field = format_ident!("{}", term.a_field);
            let b1_field = format_ident!("{}", term.b_field_1);
            let b2_field = format_ident!("{}", term.b_field_2);

            let abs_coeff = term.coeff.abs();
            let is_negative = term.coeff < 0;

            // Build the base product expression
            let base_expr = quote! { a.#a_field() * b.#b1_field() * b.#b2_field() };

            // Apply coefficient if not 1
            let coeff_expr = match abs_coeff {
                1 => base_expr,
                2 => quote! { T::TWO * #base_expr },
                n => {
                    let n_i8 = n as i8;
                    quote! { T::from_i8(#n_i8) * #base_expr }
                }
            };

            // Apply sign and position-based formatting
            let term_expr = match (i, is_negative) {
                (0, false) => coeff_expr,
                (0, true) => quote! { -(#coeff_expr) },
                (_, false) => quote! { + #coeff_expr },
                (_, true) => quote! { - #coeff_expr },
            };

            expr_parts.push(term_expr);
        }

        quote! { #(#expr_parts)* }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::spec::parse_spec;

    #[test]
    fn generates_point_onto_plane_projection() {
        let spec = parse_spec(include_str!("../../../../algebras/projective3.toml")).unwrap();
        let algebra = Algebra::new(spec.signature.p, spec.signature.q, spec.signature.r);
        let table = ProductTable::new(&algebra);
        let generator = ProjectionGenerator::new(&spec, &algebra, table);

        let tokens = generator.generate_all_projections();
        let code = tokens.to_string();

        // Should generate projection for Point onto Plane
        assert!(
            code.contains("project_point_onto_plane"),
            "Should generate Point onto Plane projection"
        );
    }
}
