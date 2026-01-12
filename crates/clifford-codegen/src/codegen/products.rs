//! Product code generation.
//!
//! This module provides the `ProductGenerator` for generating Rust code
//! for all binary products between algebra types.

use proc_macro2::TokenStream;
use quote::{format_ident, quote};

#[cfg(test)]
use crate::algebra::geometric_grades;
use crate::algebra::{Algebra, Blade, ProductTable, left_contraction_grade, outer_grade};
use crate::spec::{AlgebraSpec, ProductEntry, TypeSpec};
use crate::symbolic::{
    AtomToRust, ConstraintSimplifier, ExpressionSimplifier, ProductKind as SymbolicProductKind,
    SymbolicProduct,
};

/// The kind of product to generate.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[allow(dead_code)]
pub enum ProductKind {
    /// Geometric product (full product).
    Geometric,
    /// Outer product (wedge, grade-raising).
    Outer,
    /// Left contraction (inner product).
    LeftContraction,
    /// Regressive product (dual of outer).
    Regressive,
    /// Scalar product (grade-0 part of geometric).
    Scalar,
}

/// A term in a product expression.
#[derive(Debug, Clone)]
struct ProductTerm {
    /// Sign of this contribution (+1 or -1).
    sign: i8,
    /// Field name from type A.
    a_field: String,
    /// Field name from type B.
    b_field: String,
}

/// A term in a sandwich product expression (v * x * rev(v)).
#[derive(Debug, Clone)]
struct SandwichTerm {
    /// Sign of this contribution.
    sign: i8,
    /// First versor field (from v).
    v_field_1: String,
    /// Operand field (from x).
    x_field: String,
    /// Second versor field (from rev(v)).
    v_field_2: String,
}

/// Generates Rust code for binary products.
///
/// The generator produces functions for:
/// - Geometric products
/// - Outer products (wedge)
/// - Inner products (left contraction)
/// - Regressive products
/// - Sandwich products
/// - Scalar products
///
/// # Example
///
/// ```
/// use clifford_codegen::algebra::{Algebra, ProductTable};
/// use clifford_codegen::codegen::ProductGenerator;
/// use clifford_codegen::spec::parse_spec;
///
/// let spec = parse_spec(r#"
/// [algebra]
/// name = "euclidean2"
///
/// [signature]
/// positive = ["e1", "e2"]
///
/// [types.Scalar]
/// grades = [0]
/// fields = ["s"]
///
/// [types.Vector]
/// grades = [1]
/// fields = ["x", "y"]
///
/// [types.Bivector]
/// grades = [2]
/// fields = ["xy"]
///
/// [types.Rotor]
/// grades = [0, 2]
/// fields = ["s", "xy"]
///
/// [products.geometric]
/// Vector_Vector = "Rotor"
/// "#).unwrap();
///
/// let algebra = Algebra::euclidean(2);
/// let table = ProductTable::new(&algebra);
/// let generator = ProductGenerator::new(&spec, &algebra, table);
///
/// let tokens = generator.generate_products_file();
/// let code = tokens.to_string();
///
/// assert!(code.contains("geometric_vector_vector"));
/// ```
pub struct ProductGenerator<'a> {
    /// The algebra specification.
    spec: &'a AlgebraSpec,
    /// The algebra for blade computations.
    algebra: &'a Algebra,
    /// The product table for computing products.
    table: ProductTable,
}

impl<'a> ProductGenerator<'a> {
    /// Creates a new product generator.
    pub fn new(spec: &'a AlgebraSpec, algebra: &'a Algebra, table: ProductTable) -> Self {
        Self {
            spec,
            algebra,
            table,
        }
    }

    // ========================================================================
    // Type Resolution Helpers
    // ========================================================================

    /// Finds a TypeSpec by name.
    fn find_type(&self, name: &str) -> Option<&TypeSpec> {
        self.spec.types.iter().find(|t| t.name == name)
    }

    /// Generates the complete products.rs file.
    pub fn generate_products_file(&self) -> TokenStream {
        let header = self.generate_header();
        let imports = self.generate_imports();
        let geometric = self.generate_all_geometric();
        let outer = self.generate_all_outer();
        let inner = self.generate_all_inner();
        let scalar = self.generate_all_scalar();
        let sandwich = self.generate_all_sandwich();

        quote! {
            #header
            #imports

            // ============================================================
            // Geometric Products
            // ============================================================
            #geometric

            // ============================================================
            // Outer Products (Wedge)
            // ============================================================
            #outer

            // ============================================================
            // Inner Products (Left Contraction)
            // ============================================================
            #inner

            // ============================================================
            // Scalar Products
            // ============================================================
            #scalar

            // ============================================================
            // Sandwich Products
            // ============================================================
            #sandwich
        }
    }

    /// Generates the file header.
    fn generate_header(&self) -> TokenStream {
        let name = &self.spec.name;
        let header_doc = format!(
            r#"//! Product operations for {}.
//!
//! This file is auto-generated by clifford-codegen.
//! Do not edit manually."#,
            name
        );

        header_doc.parse().unwrap_or_else(|_| quote! {})
    }

    /// Generates import statements.
    fn generate_imports(&self) -> TokenStream {
        let type_names: Vec<_> = self
            .spec
            .types
            .iter()
            .filter(|t| t.alias_of.is_none())
            .map(|t| format_ident!("{}", t.name))
            .collect();

        quote! {
            use crate::scalar::Float;
            use super::types::{#(#type_names),*};
        }
    }

    /// Generates a constructor call for the given type.
    ///
    /// For types with solve_for fields (constrained types), uses `new_unchecked`
    /// to preserve the computed product values. For unconstrained types, uses `new`.
    fn generate_constructor_call(
        &self,
        ty: &TypeSpec,
        type_name: &proc_macro2::Ident,
        field_exprs: &[TokenStream],
    ) -> TokenStream {
        let has_constraints = !ty.solve_for_fields().is_empty();
        if has_constraints {
            quote! { #type_name::new_unchecked(#(#field_exprs),*) }
        } else {
            quote! { #type_name::new(#(#field_exprs),*) }
        }
    }

    // ========================================================================
    // Geometric Products
    // ========================================================================

    /// Generates all geometric product functions.
    fn generate_all_geometric(&self) -> TokenStream {
        // If no explicit products defined, generate nothing
        if self.spec.products.geometric.is_empty() {
            return quote! {};
        }

        let products: Vec<TokenStream> = self
            .spec
            .products
            .geometric
            .iter()
            .filter_map(|entry| self.generate_geometric_from_entry(entry))
            .collect();

        quote! { #(#products)* }
    }

    /// Generates a geometric product from a product entry.
    fn generate_geometric_from_entry(&self, entry: &ProductEntry) -> Option<TokenStream> {
        let type_a = self.find_type(&entry.lhs)?;
        let type_b = self.find_type(&entry.rhs)?;
        let output_type = self.find_type(&entry.output)?;

        let a_name = format_ident!("{}", entry.lhs);
        let b_name = format_ident!("{}", entry.rhs);
        let c_name = format_ident!("{}", entry.output);

        let fn_name = format_ident!(
            "geometric_{}_{}",
            entry.lhs.to_lowercase(),
            entry.rhs.to_lowercase()
        );

        // Use symbolic simplification for expression generation
        let field_exprs =
            self.generate_expression_symbolic(type_a, type_b, output_type, ProductKind::Geometric);

        let doc = format!(
            "Geometric product: {} * {} -> {}",
            entry.lhs, entry.rhs, entry.output
        );

        let constructor_call = self.generate_constructor_call(output_type, &c_name, &field_exprs);

        Some(quote! {
            #[doc = #doc]
            #[inline]
            pub fn #fn_name<T: Float>(a: &#a_name<T>, b: &#b_name<T>) -> #c_name<T> {
                #constructor_call
            }
        })
    }

    // ========================================================================
    // Outer Products
    // ========================================================================

    /// Generates all outer product functions.
    fn generate_all_outer(&self) -> TokenStream {
        // If no explicit products defined, generate nothing
        if self.spec.products.outer.is_empty() {
            return quote! {};
        }

        let products: Vec<TokenStream> = self
            .spec
            .products
            .outer
            .iter()
            .filter_map(|entry| self.generate_outer_from_entry(entry))
            .collect();

        quote! { #(#products)* }
    }

    /// Generates an outer product from a product entry.
    fn generate_outer_from_entry(&self, entry: &ProductEntry) -> Option<TokenStream> {
        let type_a = self.find_type(&entry.lhs)?;
        let type_b = self.find_type(&entry.rhs)?;
        let output_type = self.find_type(&entry.output)?;

        let a_name = format_ident!("{}", entry.lhs);
        let b_name = format_ident!("{}", entry.rhs);
        let c_name = format_ident!("{}", entry.output);

        let fn_name = format_ident!(
            "outer_{}_{}",
            entry.lhs.to_lowercase(),
            entry.rhs.to_lowercase()
        );

        // Use symbolic simplification for expression generation
        let field_exprs =
            self.generate_expression_symbolic(type_a, type_b, output_type, ProductKind::Outer);

        let doc = format!(
            "Outer product: {} ^ {} -> {}",
            entry.lhs, entry.rhs, entry.output
        );

        let constructor_call = self.generate_constructor_call(output_type, &c_name, &field_exprs);

        Some(quote! {
            #[doc = #doc]
            #[inline]
            pub fn #fn_name<T: Float>(a: &#a_name<T>, b: &#b_name<T>) -> #c_name<T> {
                #constructor_call
            }
        })
    }

    // ========================================================================
    // Inner Products (Left Contraction)
    // ========================================================================

    /// Generates all inner product (left contraction) functions.
    fn generate_all_inner(&self) -> TokenStream {
        // If no explicit products defined, generate nothing
        if self.spec.products.left_contraction.is_empty() {
            return quote! {};
        }

        let products: Vec<TokenStream> = self
            .spec
            .products
            .left_contraction
            .iter()
            .filter_map(|entry| self.generate_inner_from_entry(entry))
            .collect();

        quote! { #(#products)* }
    }

    /// Generates a left contraction from a product entry.
    fn generate_inner_from_entry(&self, entry: &ProductEntry) -> Option<TokenStream> {
        let type_a = self.find_type(&entry.lhs)?;
        let type_b = self.find_type(&entry.rhs)?;
        let output_type = self.find_type(&entry.output)?;

        let a_name = format_ident!("{}", entry.lhs);
        let b_name = format_ident!("{}", entry.rhs);
        let c_name = format_ident!("{}", entry.output);

        let fn_name = format_ident!(
            "left_contract_{}_{}",
            entry.lhs.to_lowercase(),
            entry.rhs.to_lowercase()
        );

        // Use symbolic simplification for expression generation
        let field_exprs = self.generate_expression_symbolic(
            type_a,
            type_b,
            output_type,
            ProductKind::LeftContraction,
        );

        let doc = format!(
            "Left contraction: {} | {} -> {}",
            entry.lhs, entry.rhs, entry.output
        );

        let constructor_call = self.generate_constructor_call(output_type, &c_name, &field_exprs);

        Some(quote! {
            #[doc = #doc]
            #[inline]
            pub fn #fn_name<T: Float>(a: &#a_name<T>, b: &#b_name<T>) -> #c_name<T> {
                #constructor_call
            }
        })
    }

    // ========================================================================
    // Scalar Products
    // ========================================================================

    /// Generates all scalar product functions.
    fn generate_all_scalar(&self) -> TokenStream {
        // If no explicit products defined, generate nothing
        if self.spec.products.scalar.is_empty() {
            return quote! {};
        }

        let products: Vec<TokenStream> = self
            .spec
            .products
            .scalar
            .iter()
            .filter_map(|entry| self.generate_scalar_from_entry(entry))
            .collect();

        quote! { #(#products)* }
    }

    /// Generates a scalar product from a product entry.
    fn generate_scalar_from_entry(&self, entry: &ProductEntry) -> Option<TokenStream> {
        let type_a = self.find_type(&entry.lhs)?;
        let type_b = self.find_type(&entry.rhs)?;

        let a_name = format_ident!("{}", entry.lhs);
        let b_name = format_ident!("{}", entry.rhs);

        let fn_name = format_ident!(
            "scalar_{}_{}",
            entry.lhs.to_lowercase(),
            entry.rhs.to_lowercase()
        );

        // Compute terms that produce grade 0 (scalar blade index = 0)
        let terms = self.compute_terms(type_a, type_b, 0, ProductKind::Geometric);
        if terms.is_empty() {
            return None;
        }
        let expr = self.generate_expression(&terms);

        let doc = format!(
            "Scalar product: {} * {} -> T (grade-0 part)",
            entry.lhs, entry.rhs
        );

        Some(quote! {
            #[doc = #doc]
            #[inline]
            pub fn #fn_name<T: Float>(a: &#a_name<T>, b: &#b_name<T>) -> T {
                #expr
            }
        })
    }

    // ========================================================================
    // Sandwich Products
    // ========================================================================

    /// Generates all sandwich product functions.
    ///
    /// Sandwich products are generated for types marked as versors in the TOML spec.
    /// For each versor type, we generate a sandwich function for each target type
    /// specified in its `sandwich.targets` list.
    fn generate_all_sandwich(&self) -> TokenStream {
        let mut products = Vec::new();

        for versor_type in &self.spec.types {
            // Skip non-versors and aliases
            if versor_type.alias_of.is_some() {
                continue;
            }

            if let Some(ref versor_spec) = versor_type.versor {
                // Get targets - either explicit or auto-detected
                let targets = if versor_spec.sandwich_targets.is_empty() {
                    // Auto-detect targets: all types that preserve grade under sandwich
                    self.infer_sandwich_targets(versor_type)
                } else {
                    versor_spec.sandwich_targets.clone()
                };

                for target_name in &targets {
                    if let Some(target_type) = self.find_type(target_name) {
                        if let Some(product) =
                            self.generate_sandwich_product_from_spec(versor_type, target_type)
                        {
                            products.push(product);
                        }
                    }
                }
            }
        }

        quote! { #(#products)* }
    }

    /// Infers valid sandwich targets for a versor type.
    ///
    /// A type is a valid target if the sandwich product V * X * rev(V) produces
    /// the same grades as X (grade-preserving transformation).
    fn infer_sandwich_targets(&self, versor_type: &TypeSpec) -> Vec<String> {
        let mut targets = Vec::new();
        let dim = self.algebra.dim();

        for candidate in &self.spec.types {
            // Skip aliases
            if candidate.alias_of.is_some() {
                continue;
            }

            // Check if sandwich preserves grades
            let output_grades =
                self.compute_sandwich_output_grades(&versor_type.grades, &candidate.grades, dim);

            // Valid target if output grades exactly match candidate grades
            let candidate_grades_set: std::collections::HashSet<usize> =
                candidate.grades.iter().copied().collect();
            if output_grades == candidate_grades_set {
                targets.push(candidate.name.clone());
            }
        }

        targets
    }

    /// Computes output grades of a sandwich product V * X * rev(V).
    fn compute_sandwich_output_grades(
        &self,
        versor_grades: &[usize],
        operand_grades: &[usize],
        dim: usize,
    ) -> std::collections::HashSet<usize> {
        use crate::algebra::geometric_grades;
        let mut output = std::collections::HashSet::new();

        for &vg in versor_grades {
            for &xg in operand_grades {
                for &wg in versor_grades {
                    // Compute possible output grades from v * x * w
                    let vx_grades = geometric_grades(vg, xg, dim);
                    for vxg in vx_grades {
                        let vxw_grades = geometric_grades(vxg, wg, dim);
                        output.extend(vxw_grades);
                    }
                }
            }
        }

        output
    }

    /// Generates a single sandwich product function from versor spec.
    fn generate_sandwich_product_from_spec(
        &self,
        versor_type: &TypeSpec,
        operand_type: &TypeSpec,
    ) -> Option<TokenStream> {
        self.generate_sandwich_product(versor_type, operand_type)
    }

    /// Checks if a type is a versor (e.g., Rotor).
    #[allow(dead_code)]
    fn is_versor_type(&self, ty: &TypeSpec) -> bool {
        ty.versor.is_some()
    }

    /// Generates a single sandwich product function.
    fn generate_sandwich_product(
        &self,
        versor_type: &TypeSpec,
        operand_type: &TypeSpec,
    ) -> Option<TokenStream> {
        let v_name = format_ident!("{}", versor_type.name);
        let x_name = format_ident!("{}", operand_type.name);

        let fn_name = format_ident!(
            "sandwich_{}_{}",
            versor_type.name.to_lowercase(),
            operand_type.name.to_lowercase()
        );

        let field_exprs: Vec<TokenStream> = operand_type
            .fields
            .iter()
            .map(|field| {
                let terms =
                    self.compute_sandwich_terms(versor_type, operand_type, field.blade_index);
                self.generate_sandwich_expression(&terms)
            })
            .collect();

        let doc = format!(
            "Sandwich product: {} * {} * rev({}) -> {}",
            versor_type.name, operand_type.name, versor_type.name, operand_type.name
        );

        let constructor_call = self.generate_constructor_call(operand_type, &x_name, &field_exprs);

        Some(quote! {
            #[doc = #doc]
            #[inline]
            pub fn #fn_name<T: Float>(v: &#v_name<T>, x: &#x_name<T>) -> #x_name<T> {
                #constructor_call
            }
        })
    }

    // ========================================================================
    // Term Computation
    // ========================================================================

    /// Computes output grades for a product.
    #[cfg(test)]
    fn compute_output_grades(
        &self,
        type_a: &TypeSpec,
        type_b: &TypeSpec,
        kind: ProductKind,
    ) -> Vec<usize> {
        let dim = self.algebra.dim();
        let mut grades = Vec::new();

        for &ga in &type_a.grades {
            for &gb in &type_b.grades {
                let result_grades = match kind {
                    ProductKind::Geometric => geometric_grades(ga, gb, dim),
                    ProductKind::Outer => {
                        if let Some(g) = outer_grade(ga, gb, dim) {
                            vec![g]
                        } else {
                            vec![]
                        }
                    }
                    ProductKind::LeftContraction => {
                        if let Some(g) = left_contraction_grade(ga, gb) {
                            vec![g]
                        } else {
                            vec![]
                        }
                    }
                    ProductKind::Regressive => {
                        // Regressive product: (a* ^ b*)* where * is dual
                        // Result grade = ga + gb - dim
                        let result = ga + gb;
                        if result >= dim {
                            vec![result - dim]
                        } else {
                            vec![]
                        }
                    }
                    ProductKind::Scalar => {
                        if ga == gb {
                            vec![0]
                        } else {
                            vec![]
                        }
                    }
                };

                for g in result_grades {
                    if !grades.contains(&g) {
                        grades.push(g);
                    }
                }
            }
        }

        grades.sort();
        grades
    }

    /// Computes product terms for a given output blade.
    fn compute_terms(
        &self,
        type_a: &TypeSpec,
        type_b: &TypeSpec,
        result_blade: usize,
        kind: ProductKind,
    ) -> Vec<ProductTerm> {
        let mut terms = Vec::new();

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
                        outer_grade(a_grade, b_grade, self.algebra.dim())
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
                    ProductKind::Regressive => {
                        // TODO: implement regressive product filtering
                        false
                    }
                    ProductKind::Scalar => result_blade == 0,
                };

                if include {
                    terms.push(ProductTerm {
                        sign,
                        a_field: field_a.name.clone(),
                        b_field: field_b.name.clone(),
                    });
                }
            }
        }

        terms
    }

    /// Computes sandwich product terms: v * x * rev(v) -> result_blade.
    fn compute_sandwich_terms(
        &self,
        versor_type: &TypeSpec,
        operand_type: &TypeSpec,
        result_blade: usize,
    ) -> Vec<SandwichTerm> {
        let mut terms = Vec::new();

        // For each combination: v_i * x_j * rev(v_k)
        for field_v1 in &versor_type.fields {
            for field_x in &operand_type.fields {
                for field_v2 in &versor_type.fields {
                    let v1_blade = field_v1.blade_index;
                    let x_blade = field_x.blade_index;
                    let v2_blade = field_v2.blade_index;

                    // Compute v_i * x_j
                    let (sign_vx, vx) = self.table.geometric(v1_blade, x_blade);
                    if sign_vx == 0 {
                        continue;
                    }

                    // Compute (v_i * x_j) * rev(v_k)
                    // rev(v_k) has sign (-1)^(k(k-1)/2) for grade k
                    let v2_grade = Blade::from_index(v2_blade).grade();
                    #[allow(clippy::manual_is_multiple_of)]
                    let rev_sign: i8 = if (v2_grade * v2_grade.saturating_sub(1) / 2) % 2 == 0 {
                        1
                    } else {
                        -1
                    };

                    let (sign_vxr, result) = self.table.geometric(vx, v2_blade);
                    if sign_vxr == 0 {
                        continue;
                    }

                    if result == result_blade {
                        let final_sign = sign_vx * sign_vxr * rev_sign;
                        terms.push(SandwichTerm {
                            sign: final_sign,
                            v_field_1: field_v1.name.clone(),
                            x_field: field_x.name.clone(),
                            v_field_2: field_v2.name.clone(),
                        });
                    }
                }
            }
        }

        self.simplify_sandwich_terms(terms)
    }

    /// Simplifies sandwich terms by combining like terms.
    fn simplify_sandwich_terms(&self, terms: Vec<SandwichTerm>) -> Vec<SandwichTerm> {
        use std::collections::HashMap;

        let mut combined: HashMap<(String, String, String), i8> = HashMap::new();

        for term in terms {
            let key = (
                term.v_field_1.clone(),
                term.x_field.clone(),
                term.v_field_2.clone(),
            );
            *combined.entry(key).or_insert(0) += term.sign;
        }

        combined
            .into_iter()
            .filter(|(_, sign)| *sign != 0)
            .map(|((v1, x, v2), sign)| SandwichTerm {
                sign, // Keep the accumulated coefficient (can be 2, 3, -2, etc.)
                v_field_1: v1,
                x_field: x,
                v_field_2: v2,
            })
            .collect()
    }

    // ========================================================================
    // Symbolic Expression Generation
    // ========================================================================

    /// Generates a simplified Rust expression using symbolic algebra.
    ///
    /// This method uses Symbolica to:
    /// 1. Build a symbolic expression for the product
    /// 2. Apply constraint substitutions (e.g., unit norm)
    /// 3. Simplify by expanding and collecting like terms
    /// 4. Convert the simplified expression to Rust code
    ///
    /// This produces more compact code than the term-based approach when
    /// there are opportunities for simplification (e.g., constraint substitution).
    fn generate_expression_symbolic(
        &self,
        type_a: &TypeSpec,
        type_b: &TypeSpec,
        output_type: &TypeSpec,
        kind: ProductKind,
    ) -> Vec<TokenStream> {
        let symbolic_product = SymbolicProduct::new(self.algebra);
        let expr_simplifier = ExpressionSimplifier::new();

        // Create constraint simplifier for input type constraints
        let constraint_simplifier = ConstraintSimplifier::new(&[type_a, type_b], &["a", "b"]);

        // Create symbolic field variables
        let a_symbols = symbolic_product.create_field_symbols(type_a, "a");
        let b_symbols = symbolic_product.create_field_symbols(type_b, "b");

        // Map ProductKind to SymbolicProductKind
        let symbolic_kind = match kind {
            ProductKind::Geometric => SymbolicProductKind::Geometric,
            ProductKind::Outer => SymbolicProductKind::Outer,
            ProductKind::LeftContraction => SymbolicProductKind::LeftContraction,
            _ => SymbolicProductKind::Geometric, // Fallback for other kinds
        };

        // Compute symbolic product
        let symbolic_fields = symbolic_product.compute(
            type_a,
            type_b,
            output_type,
            symbolic_kind,
            &a_symbols,
            &b_symbols,
        );

        // Create converter for Rust code generation
        let converter = AtomToRust::new(&[type_a, type_b], &["a", "b"]);

        // Apply constraint substitution, simplify, and convert each field expression
        symbolic_fields
            .iter()
            .map(|field| {
                // First apply constraint substitutions (e.g., s*s + xy*xy + ... = 1)
                let with_constraints = constraint_simplifier.apply(&field.expression);
                // Then simplify (expand and collect like terms)
                let simplified = expr_simplifier.simplify(&with_constraints);
                converter.convert(&simplified)
            })
            .collect()
    }

    // ========================================================================
    // Expression Generation
    // ========================================================================

    /// Generates a Rust expression from product terms.
    fn generate_expression(&self, terms: &[ProductTerm]) -> TokenStream {
        if terms.is_empty() {
            return quote! { T::zero() };
        }

        let mut expr_parts: Vec<TokenStream> = Vec::new();

        for (i, term) in terms.iter().enumerate() {
            let a_field = format_ident!("{}", term.a_field);
            let b_field = format_ident!("{}", term.b_field);

            let term_expr = match (i, term.sign) {
                (0, s) if s > 0 => quote! { a.#a_field() * b.#b_field() },
                (0, _) => quote! { -(a.#a_field() * b.#b_field()) },
                (_, s) if s > 0 => quote! { + a.#a_field() * b.#b_field() },
                (_, _) => quote! { - a.#a_field() * b.#b_field() },
            };

            expr_parts.push(term_expr);
        }

        quote! { #(#expr_parts)* }
    }

    /// Generates a Rust expression from sandwich terms.
    fn generate_sandwich_expression(&self, terms: &[SandwichTerm]) -> TokenStream {
        if terms.is_empty() {
            return quote! { T::zero() };
        }

        let mut expr_parts: Vec<TokenStream> = Vec::new();

        for (i, term) in terms.iter().enumerate() {
            let v1 = format_ident!("{}", term.v_field_1);
            let x_field = format_ident!("{}", term.x_field);
            let v2 = format_ident!("{}", term.v_field_2);

            let abs_coeff = term.sign.abs();
            let is_negative = term.sign < 0;

            // Build the base product expression
            let base_expr = quote! { v.#v1() * x.#x_field() * v.#v2() };

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
    fn generates_geometric_product() {
        let spec = parse_spec(include_str!("../../algebras/euclidean2.toml")).unwrap();
        let algebra = Algebra::euclidean(2);
        let table = ProductTable::new(&algebra);
        let generator = ProductGenerator::new(&spec, &algebra, table);

        let tokens = generator.generate_products_file();
        let code = tokens.to_string();

        assert!(code.contains("geometric_vector_vector"));
        assert!(code.contains("geometric_rotor_rotor"));
    }

    #[test]
    fn generates_outer_product() {
        let spec = parse_spec(include_str!("../../algebras/euclidean3.toml")).unwrap();
        let algebra = Algebra::euclidean(3);
        let table = ProductTable::new(&algebra);
        let generator = ProductGenerator::new(&spec, &algebra, table);

        let tokens = generator.generate_products_file();
        let code = tokens.to_string();

        assert!(code.contains("outer_vector_vector"));
    }

    #[test]
    fn skips_left_contraction_unless_specified() {
        // Left contraction is not auto-generated to avoid dead code
        // It will only be generated when explicitly specified in the spec
        let spec = parse_spec(include_str!("../../algebras/euclidean3.toml")).unwrap();
        let algebra = Algebra::euclidean(3);
        let table = ProductTable::new(&algebra);
        let generator = ProductGenerator::new(&spec, &algebra, table);

        let tokens = generator.generate_products_file();
        let code = tokens.to_string();

        // Left contraction should NOT be generated by default
        assert!(!code.contains("left_contract"));
    }

    #[test]
    fn generates_sandwich_for_marked_versors() {
        // Sandwich products are generated for types marked with versor = true
        // and explicitly listed targets in [types.TypeName.sandwich]
        let spec = parse_spec(include_str!("../../algebras/euclidean3.toml")).unwrap();
        let algebra = Algebra::euclidean(3);
        let table = ProductTable::new(&algebra);
        let generator = ProductGenerator::new(&spec, &algebra, table);

        let tokens = generator.generate_products_file();
        let code = tokens.to_string();

        // Rotor is marked as versor with targets Vector, Bivector, Rotor
        assert!(code.contains("sandwich_rotor_vector"));
        assert!(code.contains("sandwich_rotor_bivector"));
        assert!(code.contains("sandwich_rotor_rotor"));
    }

    #[test]
    fn skips_sandwich_without_versor_marking() {
        // Sandwich products are NOT generated for types without versor = true
        // Use euclidean2 which has no versor marking
        let spec = parse_spec(include_str!("../../algebras/euclidean2.toml")).unwrap();
        let algebra = Algebra::euclidean(2);
        let table = ProductTable::new(&algebra);
        let generator = ProductGenerator::new(&spec, &algebra, table);

        let tokens = generator.generate_products_file();
        let code = tokens.to_string();

        // Sandwich should NOT be generated (no versor marking)
        assert!(!code.contains("sandwich_"));
    }

    #[test]
    fn term_computation_vector_vector() {
        let spec = parse_spec(include_str!("../../algebras/euclidean3.toml")).unwrap();
        let algebra = Algebra::euclidean(3);
        let table = ProductTable::new(&algebra);
        let generator = ProductGenerator::new(&spec, &algebra, table);

        let vector = spec.types.iter().find(|t| t.name == "Vector").unwrap();

        // Scalar result from vector * vector
        let scalar_terms = generator.compute_terms(vector, vector, 0, ProductKind::Geometric);

        // Should have 3 terms: x*x + y*y + z*z
        assert_eq!(scalar_terms.len(), 3);
        assert!(scalar_terms.iter().all(|t| t.sign > 0));
    }

    #[test]
    fn output_grades_geometric() {
        let spec = parse_spec(include_str!("../../algebras/euclidean3.toml")).unwrap();
        let algebra = Algebra::euclidean(3);
        let table = ProductTable::new(&algebra);
        let generator = ProductGenerator::new(&spec, &algebra, table);

        let vector = spec.types.iter().find(|t| t.name == "Vector").unwrap();

        let grades = generator.compute_output_grades(vector, vector, ProductKind::Geometric);
        // Vector * Vector in 3D produces grades 0 and 2
        assert_eq!(grades, vec![0, 2]);
    }

    #[test]
    fn output_grades_outer() {
        let spec = parse_spec(include_str!("../../algebras/euclidean3.toml")).unwrap();
        let algebra = Algebra::euclidean(3);
        let table = ProductTable::new(&algebra);
        let generator = ProductGenerator::new(&spec, &algebra, table);

        let vector = spec.types.iter().find(|t| t.name == "Vector").unwrap();

        let grades = generator.compute_output_grades(vector, vector, ProductKind::Outer);
        // Vector ^ Vector produces grade 2
        assert_eq!(grades, vec![2]);
    }
}
