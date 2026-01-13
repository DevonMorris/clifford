//! Product code generation.
//!
//! This module provides the `ProductGenerator` for generating Rust code
//! for all binary products between algebra types.

use proc_macro2::TokenStream;
use quote::{format_ident, quote};

#[cfg(test)]
use crate::algebra::geometric_grades;
#[cfg(test)]
use crate::algebra::left_contraction_grade;
#[cfg(test)]
use crate::algebra::outer_grade;
use crate::algebra::{Algebra, Blade, ProductTable};
use crate::spec::{AlgebraSpec, ProductEntry, TypeSpec};
use crate::symbolic::{
    AtomToRust, ConstraintSimplifier, ExpressionSimplifier, ProductKind as SymbolicProductKind,
    SymbolicProduct,
};

use super::unary::UnaryGenerator;

/// The kind of product to generate.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[allow(dead_code)]
pub enum ProductKind {
    /// Geometric product (full product).
    Geometric,
    /// Exterior product (wedge, grade-raising).
    Exterior,
    /// Interior product (symmetric inner, grade |ga - gb|).
    Interior,
    /// Left contraction (A ⌋ B, grade gb - ga when ga <= gb).
    LeftContraction,
    /// Right contraction (A ⌊ B, grade ga - gb when gb <= ga).
    RightContraction,
    /// Regressive product (meet, dual of exterior, grade ga + gb - dim).
    Regressive,
    /// Scalar product (grade-0 part of geometric).
    Scalar,
    /// Geometric antiproduct (complement(complement(a) * complement(b))).
    Antigeometric,
    /// Antiscalar product (grade-n part of antigeometric).
    Antiscalar,
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
/// - Exterior products (wedge)
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
        let exterior = self.generate_all_exterior();
        let interior = self.generate_all_interior();
        let left_contraction = self.generate_all_inner(); // Left contraction (inner)
        let right_contraction = self.generate_all_right_contraction();
        let regressive = self.generate_all_regressive();
        let scalar = self.generate_all_scalar();
        let antigeometric = self.generate_all_antigeometric();
        let sandwich = self.generate_all_sandwich();
        let antisandwich = self.generate_all_antisandwich();

        // Generate unary operations
        let unary_gen = UnaryGenerator::new(self.spec);
        let unary = unary_gen.generate_all();

        quote! {
            #header
            #imports

            // ============================================================
            // Geometric Products
            // ============================================================
            #geometric

            // ============================================================
            // Exterior Products (Wedge)
            // ============================================================
            #exterior

            // ============================================================
            // Interior Products (Symmetric Inner)
            // ============================================================
            #interior

            // ============================================================
            // Left Contraction Products
            // ============================================================
            #left_contraction

            // ============================================================
            // Right Contraction Products
            // ============================================================
            #right_contraction

            // ============================================================
            // Regressive Products (Meet)
            // ============================================================
            #regressive

            // ============================================================
            // Scalar Products
            // ============================================================
            #scalar

            // ============================================================
            // Antigeometric Products
            // ============================================================
            #antigeometric

            // ============================================================
            // Sandwich Products (Geometric)
            // ============================================================
            #sandwich

            // ============================================================
            // Antisandwich Products (for PGA transformations)
            // ============================================================
            #antisandwich

            // ============================================================
            // Unary Operations (Reverse, Antireverse, Complement)
            // ============================================================
            #unary
        }
    }

    /// Generates the file header with comprehensive module documentation.
    fn generate_header(&self) -> TokenStream {
        let name = &self.spec.name;

        // Check if there are any versors to include sandwich documentation
        let has_versors = self.spec.types.iter().any(|t| t.versor.is_some());

        let sandwich_section = if has_versors {
            r#"
//!
//! # Sandwich Products
//!
//! For versor types (rotors, motors), sandwich products are provided:
//! - `sandwich_{versor}_{operand}(v, x)` computes `v × x × rev(v)`"#
        } else {
            ""
        };

        let header_doc = format!(
            r#"//! Product functions for the {} algebra.
//!
//! This module provides all algebraic products between types in the algebra.
//! Each function is named `{{product}}_{{lhs}}_{{rhs}}` where:
//! - `product` is one of: `geometric`, `exterior`, `left_contract`, `inner`, `scalar`
//! - `lhs` and `rhs` are the input type names in lowercase
//!
//! # Available Products
//!
//! | Product | Symbol | Description |
//! |---------|--------|-------------|
//! | `geometric_*` | `×` | Full geometric product |
//! | `exterior_*` | `∧` | Wedge/exterior product (grade sum) |
//! | `left_contract_*` | `⌋` | Left contraction |
//! | `inner_*` | `·` | Symmetric inner product (grade diff) |
//! | `scalar_*` | `⟨⟩₀` | Scalar (grade-0) product |{}
//!
//! This file is auto-generated by clifford-codegen. Do not edit manually."#,
            name, sandwich_section
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
    /// For types with constraints (solve_for fields), uses `new_unchecked` because:
    /// 1. Product outputs are mathematically correct as computed
    /// 2. Constraint solving would incorrectly modify the algebraic result
    /// 3. Constraints apply to *normalized* instances, not intermediate results
    ///
    /// For unconstrained types, uses the standard `new` constructor.
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
    // Exterior Products
    // ========================================================================

    /// Generates all exterior product functions.
    fn generate_all_exterior(&self) -> TokenStream {
        // If no explicit products defined, generate nothing
        if self.spec.products.exterior.is_empty() {
            return quote! {};
        }

        let products: Vec<TokenStream> = self
            .spec
            .products
            .exterior
            .iter()
            .filter_map(|entry| self.generate_exterior_from_entry(entry))
            .collect();

        quote! { #(#products)* }
    }

    /// Generates an exterior product from a product entry.
    fn generate_exterior_from_entry(&self, entry: &ProductEntry) -> Option<TokenStream> {
        let type_a = self.find_type(&entry.lhs)?;
        let type_b = self.find_type(&entry.rhs)?;
        let output_type = self.find_type(&entry.output)?;

        let a_name = format_ident!("{}", entry.lhs);
        let b_name = format_ident!("{}", entry.rhs);
        let c_name = format_ident!("{}", entry.output);

        let fn_name = format_ident!(
            "exterior_{}_{}",
            entry.lhs.to_lowercase(),
            entry.rhs.to_lowercase()
        );

        // Use symbolic simplification for expression generation
        let field_exprs =
            self.generate_expression_symbolic(type_a, type_b, output_type, ProductKind::Exterior);

        let doc = format!(
            "Exterior product: {} ^ {} -> {}",
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
    // Interior Products (Symmetric Inner)
    // ========================================================================

    /// Generates all interior (symmetric inner) product functions.
    fn generate_all_interior(&self) -> TokenStream {
        if self.spec.products.interior.is_empty() {
            return quote! {};
        }

        let products: Vec<TokenStream> = self
            .spec
            .products
            .interior
            .iter()
            .filter_map(|entry| self.generate_interior_from_entry(entry))
            .collect();

        quote! { #(#products)* }
    }

    /// Generates an interior product from a product entry.
    fn generate_interior_from_entry(&self, entry: &ProductEntry) -> Option<TokenStream> {
        let type_a = self.find_type(&entry.lhs)?;
        let type_b = self.find_type(&entry.rhs)?;
        let output_type = self.find_type(&entry.output)?;

        let a_name = format_ident!("{}", entry.lhs);
        let b_name = format_ident!("{}", entry.rhs);
        let c_name = format_ident!("{}", entry.output);

        let fn_name = format_ident!(
            "interior_{}_{}",
            entry.lhs.to_lowercase(),
            entry.rhs.to_lowercase()
        );

        let field_exprs =
            self.generate_expression_symbolic(type_a, type_b, output_type, ProductKind::Interior);

        let doc = format!(
            "Interior product (symmetric inner): {} · {} -> {}",
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
    // Right Contraction Products
    // ========================================================================

    /// Generates all right contraction product functions.
    fn generate_all_right_contraction(&self) -> TokenStream {
        if self.spec.products.right_contraction.is_empty() {
            return quote! {};
        }

        let products: Vec<TokenStream> = self
            .spec
            .products
            .right_contraction
            .iter()
            .filter_map(|entry| self.generate_right_contraction_from_entry(entry))
            .collect();

        quote! { #(#products)* }
    }

    /// Generates a right contraction from a product entry.
    fn generate_right_contraction_from_entry(&self, entry: &ProductEntry) -> Option<TokenStream> {
        let type_a = self.find_type(&entry.lhs)?;
        let type_b = self.find_type(&entry.rhs)?;
        let output_type = self.find_type(&entry.output)?;

        let a_name = format_ident!("{}", entry.lhs);
        let b_name = format_ident!("{}", entry.rhs);
        let c_name = format_ident!("{}", entry.output);

        let fn_name = format_ident!(
            "right_contract_{}_{}",
            entry.lhs.to_lowercase(),
            entry.rhs.to_lowercase()
        );

        let field_exprs = self.generate_expression_symbolic(
            type_a,
            type_b,
            output_type,
            ProductKind::RightContraction,
        );

        let doc = format!(
            "Right contraction: {} ⌊ {} -> {}",
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
    // Regressive Products (Meet)
    // ========================================================================

    /// Generates all regressive (meet) product functions.
    fn generate_all_regressive(&self) -> TokenStream {
        if self.spec.products.regressive.is_empty() {
            return quote! {};
        }

        let products: Vec<TokenStream> = self
            .spec
            .products
            .regressive
            .iter()
            .filter_map(|entry| self.generate_regressive_from_entry(entry))
            .collect();

        quote! { #(#products)* }
    }

    /// Generates a regressive product from a product entry.
    ///
    /// Uses term-based generation since regressive products use the complement-based
    /// formula `a ∨ b = ∁(∁a ∧ ∁b)` which requires the specialized table method.
    fn generate_regressive_from_entry(&self, entry: &ProductEntry) -> Option<TokenStream> {
        let type_a = self.find_type(&entry.lhs)?;
        let type_b = self.find_type(&entry.rhs)?;
        let output_type = self.find_type(&entry.output)?;

        let a_name = format_ident!("{}", entry.lhs);
        let b_name = format_ident!("{}", entry.rhs);
        let c_name = format_ident!("{}", entry.output);

        let fn_name = format_ident!(
            "regressive_{}_{}",
            entry.lhs.to_lowercase(),
            entry.rhs.to_lowercase()
        );

        // Generate field expressions using term-based computation
        // This uses table.regressive() which correctly computes ∁(∁a ∧ ∁b)
        let field_exprs: Vec<TokenStream> = output_type
            .fields
            .iter()
            .map(|field| {
                let terms =
                    self.compute_terms(type_a, type_b, field.blade_index, ProductKind::Regressive);
                self.generate_expression(&terms)
            })
            .collect();

        let doc = format!(
            "Regressive product (meet): {} ∨ {} -> {}",
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
    // Antigeometric Products
    // ========================================================================

    /// Generates all antigeometric product functions.
    ///
    /// The antigeometric product is defined as: a ⊛ b = complement(complement(a) * complement(b))
    fn generate_all_antigeometric(&self) -> TokenStream {
        // If no explicit products defined, generate nothing
        if self.spec.products.antigeometric.is_empty() {
            return quote! {};
        }

        let products: Vec<TokenStream> = self
            .spec
            .products
            .antigeometric
            .iter()
            .filter_map(|entry| self.generate_antigeometric_from_entry(entry))
            .collect();

        quote! { #(#products)* }
    }

    /// Generates an antigeometric product from a product entry.
    fn generate_antigeometric_from_entry(&self, entry: &ProductEntry) -> Option<TokenStream> {
        let type_a = self.find_type(&entry.lhs)?;
        let type_b = self.find_type(&entry.rhs)?;
        let output_type = self.find_type(&entry.output)?;

        let a_name = format_ident!("{}", entry.lhs);
        let b_name = format_ident!("{}", entry.rhs);
        let c_name = format_ident!("{}", entry.output);

        let fn_name = format_ident!(
            "antigeometric_{}_{}",
            entry.lhs.to_lowercase(),
            entry.rhs.to_lowercase()
        );

        // Use symbolic simplification for expression generation
        let field_exprs = self.generate_expression_symbolic(
            type_a,
            type_b,
            output_type,
            ProductKind::Antigeometric,
        );

        let doc = format!(
            "Antigeometric product: {} ⊛ {} -> {}\n\nDefined as complement(complement(a) * complement(b)).",
            entry.lhs, entry.rhs, entry.output
        );

        let constructor_call = self.generate_constructor_call(output_type, &c_name, &field_exprs);

        // Check if parameters are used in expressions
        let exprs_str = field_exprs
            .iter()
            .map(|e| e.to_string())
            .collect::<Vec<_>>()
            .join(" ");
        let a_used = exprs_str.contains("a .");
        let b_used = exprs_str.contains("b .");

        let param_a = if a_used {
            quote! { a: &#a_name<T> }
        } else {
            quote! { _a: &#a_name<T> }
        };
        let param_b = if b_used {
            quote! { b: &#b_name<T> }
        } else {
            quote! { _b: &#b_name<T> }
        };

        Some(quote! {
            #[doc = #doc]
            #[inline]
            pub fn #fn_name<T: Float>(#param_a, #param_b) -> #c_name<T> {
                #constructor_call
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

    /// Generates all antisandwich product functions.
    ///
    /// Antisandwich products use the geometric antiproduct instead of the
    /// geometric product. In PGA (Projective GA), antisandwich is required
    /// for correct motor transformations because it handles the degenerate
    /// metric (e0² = 0) properly.
    ///
    /// For each versor type, we generate an antisandwich function for each
    /// target type specified in its `sandwich.targets` list.
    fn generate_all_antisandwich(&self) -> TokenStream {
        let mut products = Vec::new();

        for versor_type in &self.spec.types {
            // Skip non-versors and aliases
            if versor_type.alias_of.is_some() {
                continue;
            }

            if let Some(ref versor_spec) = versor_type.versor {
                // Get targets - either explicit or auto-detected
                let targets = if versor_spec.sandwich_targets.is_empty() {
                    self.infer_sandwich_targets(versor_type)
                } else {
                    versor_spec.sandwich_targets.clone()
                };

                for target_name in &targets {
                    if let Some(target_type) = self.find_type(target_name) {
                        if let Some(product) =
                            self.generate_antisandwich_product(versor_type, target_type)
                        {
                            products.push(product);
                        }
                    }
                }
            }
        }

        quote! { #(#products)* }
    }

    /// Generates a single antisandwich product function.
    ///
    /// The antisandwich computes: v ⊛ x ⊛ antirev(v) where ⊛ is the geometric
    /// antiproduct.
    ///
    /// Returns None if the antiproduct produces no non-zero terms (which can
    /// happen in degenerate algebras like PGA where e0² = 0).
    fn generate_antisandwich_product(
        &self,
        versor_type: &TypeSpec,
        operand_type: &TypeSpec,
    ) -> Option<TokenStream> {
        let v_name = format_ident!("{}", versor_type.name);
        let x_name = format_ident!("{}", operand_type.name);

        let fn_name = format_ident!(
            "antisandwich_{}_{}",
            versor_type.name.to_lowercase(),
            operand_type.name.to_lowercase()
        );

        // Collect terms for each field, tracking if any field has non-empty terms
        let mut has_any_terms = false;
        let field_exprs: Vec<TokenStream> = operand_type
            .fields
            .iter()
            .map(|field| {
                let terms =
                    self.compute_antisandwich_terms(versor_type, operand_type, field.blade_index);
                if !terms.is_empty() {
                    has_any_terms = true;
                }
                self.generate_sandwich_expression(&terms)
            })
            .collect();

        // Skip generating function if no terms were computed (degenerate case)
        if !has_any_terms {
            return None;
        }

        let doc = format!(
            "Antisandwich product: {} ⊛ {} ⊛ antirev({}) -> {}\n\n\
             Uses the geometric antiproduct and antireverse. In PGA, use this\n\
             instead of the regular sandwich for correct motor transformations.",
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
                    ProductKind::Exterior => {
                        if let Some(g) = outer_grade(ga, gb, dim) {
                            vec![g]
                        } else {
                            vec![]
                        }
                    }
                    ProductKind::Interior => {
                        // Interior product: |ga - gb|
                        vec![ga.abs_diff(gb)]
                    }
                    ProductKind::LeftContraction => {
                        if let Some(g) = left_contraction_grade(ga, gb) {
                            vec![g]
                        } else {
                            vec![]
                        }
                    }
                    ProductKind::RightContraction => {
                        // Right contraction: ga - gb when gb <= ga
                        if gb <= ga { vec![ga - gb] } else { vec![] }
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
                    ProductKind::Antigeometric => {
                        // Antigeometric: same grade rules as geometric but with antigrades
                        // antigrade(a ⊛ b) = antigrade(a) + antigrade(b) ± 2k
                        // where antigrade = n - grade
                        // For now, use full geometric grades as the antiproduct produces similar structure
                        geometric_grades(ga, gb, dim)
                    }
                    ProductKind::Antiscalar => {
                        // Antiscalar: grade-n part of antigeometric
                        vec![dim]
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

                // Use specialized table methods for each product kind
                let (sign, result) = match kind {
                    ProductKind::Geometric => self.table.geometric(a_blade, b_blade),
                    ProductKind::Exterior => self.table.exterior(a_blade, b_blade),
                    ProductKind::Interior => self.table.interior(a_blade, b_blade),
                    ProductKind::LeftContraction => self.table.left_contraction(a_blade, b_blade),
                    ProductKind::RightContraction => self.table.right_contraction(a_blade, b_blade),
                    ProductKind::Regressive => self.table.regressive(a_blade, b_blade),
                    ProductKind::Scalar => self.table.scalar(a_blade, b_blade),
                    ProductKind::Antigeometric => self.table.antiproduct(a_blade, b_blade),
                    ProductKind::Antiscalar => self.table.antiscalar(a_blade, b_blade),
                };

                if result == result_blade && sign != 0 {
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

    /// Computes antisandwich product terms: v ⊛ x ⊛ antirev(v) -> result_blade.
    ///
    /// Uses the geometric antiproduct and antireverse instead of the
    /// geometric product and reverse. The antireverse is defined as:
    /// `antirev(x) = complement(reverse(complement(x)))`
    ///
    /// The antireverse sign for grade k in dimension n is: `(-1)^((n-k)(n-k-1)/2)`
    fn compute_antisandwich_terms(
        &self,
        versor_type: &TypeSpec,
        operand_type: &TypeSpec,
        result_blade: usize,
    ) -> Vec<SandwichTerm> {
        let mut terms = Vec::new();
        let dim = self.algebra.dim();

        // For each combination: v_i ⊛ x_j ⊛ antirev(v_k)
        for field_v1 in &versor_type.fields {
            for field_x in &operand_type.fields {
                for field_v2 in &versor_type.fields {
                    let v1_blade = field_v1.blade_index;
                    let x_blade = field_x.blade_index;
                    let v2_blade = field_v2.blade_index;

                    // Compute v_i ⊛ x_j using antiproduct
                    let (sign_vx, vx) = self.table.antiproduct(v1_blade, x_blade);
                    if sign_vx == 0 {
                        continue;
                    }

                    // Compute (v_i ⊛ x_j) ⊛ antirev(v_k)
                    // antirev(v_k) has sign (-1)^((n-k)(n-k-1)/2) for grade k in dimension n
                    let v2_grade = Blade::from_index(v2_blade).grade();
                    let antigrade = dim - v2_grade;
                    #[allow(clippy::manual_is_multiple_of)]
                    let antirev_sign: i8 = if (antigrade * antigrade.saturating_sub(1) / 2) % 2 == 0
                    {
                        1
                    } else {
                        -1
                    };

                    let (sign_vxr, result) = self.table.antiproduct(vx, v2_blade);
                    if sign_vxr == 0 {
                        continue;
                    }

                    if result == result_blade {
                        let final_sign = sign_vx * sign_vxr * antirev_sign;
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
            ProductKind::Exterior => SymbolicProductKind::Exterior,
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
                    quote! { T::from_i8(#n) * #base_expr }
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
        let spec = parse_spec(include_str!("../../../../algebras/euclidean2.toml")).unwrap();
        let algebra = Algebra::euclidean(2);
        let table = ProductTable::new(&algebra);
        let generator = ProductGenerator::new(&spec, &algebra, table);

        let tokens = generator.generate_products_file();
        let code = tokens.to_string();

        assert!(code.contains("geometric_vector_vector"));
        assert!(code.contains("geometric_rotor_rotor"));
    }

    #[test]
    fn generates_exterior_product() {
        let spec = parse_spec(include_str!("../../../../algebras/euclidean3.toml")).unwrap();
        let algebra = Algebra::euclidean(3);
        let table = ProductTable::new(&algebra);
        let generator = ProductGenerator::new(&spec, &algebra, table);

        let tokens = generator.generate_products_file();
        let code = tokens.to_string();

        assert!(code.contains("exterior_vector_vector"));
    }

    #[test]
    fn generates_left_contraction() {
        // Left contraction is auto-generated for all algebras
        // Left contraction a ⌋ b only gives non-zero when grade(a) <= grade(b)
        // Valid products: Vector ⌋ Bivector → Vector (grade 2-1=1)
        let spec = parse_spec(include_str!("../../../../algebras/euclidean3.toml")).unwrap();
        let algebra = Algebra::euclidean(3);
        let table = ProductTable::new(&algebra);
        let generator = ProductGenerator::new(&spec, &algebra, table);

        let tokens = generator.generate_products_file();
        let code = tokens.to_string();

        // Vector ⌋ Bivector → Vector is a valid left contraction
        assert!(code.contains("left_contract_vector_bivector"));
    }

    #[test]
    fn generates_sandwich_for_marked_versors() {
        // Sandwich products are generated for types marked with versor = true
        // and explicitly listed targets in [types.TypeName.sandwich]
        let spec = parse_spec(include_str!("../../../../algebras/euclidean3.toml")).unwrap();
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
        let spec = parse_spec(include_str!("../../../../algebras/euclidean2.toml")).unwrap();
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
        let spec = parse_spec(include_str!("../../../../algebras/euclidean3.toml")).unwrap();
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
        let spec = parse_spec(include_str!("../../../../algebras/euclidean3.toml")).unwrap();
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
        let spec = parse_spec(include_str!("../../../../algebras/euclidean3.toml")).unwrap();
        let algebra = Algebra::euclidean(3);
        let table = ProductTable::new(&algebra);
        let generator = ProductGenerator::new(&spec, &algebra, table);

        let vector = spec.types.iter().find(|t| t.name == "Vector").unwrap();

        let grades = generator.compute_output_grades(vector, vector, ProductKind::Exterior);
        // Vector ^ Vector produces grade 2
        assert_eq!(grades, vec![2]);
    }
}
