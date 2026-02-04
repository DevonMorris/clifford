//! Trait implementation generation.
//!
//! This module provides the `TraitsGenerator` for generating Rust trait
//! implementations including operators, approx traits, and Arbitrary.

use proc_macro2::TokenStream;
use quote::{format_ident, quote};

use crate::algebra::{Algebra, Blade, ProductTable};
use crate::spec::{AlgebraSpec, InvolutionKind, TypeSpec, WrapperKind};

/// Position of wrapper in a binary product.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum WrapperPosition {
    /// LHS is wrapped: `Trait<B> for Unit<A>`
    Lhs,
    /// RHS is wrapped: `Trait<Unit<B>> for A`
    Rhs,
    /// Both are wrapped: `Trait<Unit<B>> for Unit<A>`
    Both,
}

use crate::symbolic::{
    AtomToRust, ConstraintDeriver, ConstraintSolver, ExpressionSimplifier, GroebnerSimplifier,
    ProductConstraintCollector, ProductKind as SymbolicProductKind, SolutionType, SymbolicProduct,
};

/// Generates trait implementations for algebra types.
///
/// The generator produces:
/// - Arithmetic operators (Add, Sub, Neg, Mul)
/// - Scalar multiplication
/// - Geometric and outer product operators
/// - Approx traits (AbsDiffEq, RelativeEq, UlpsEq)
/// - Arbitrary implementations for proptest
///
/// # Example
///
/// ```
/// use clifford_codegen::algebra::{Algebra, ProductTable};
/// use clifford_codegen::codegen::TraitsGenerator;
/// use clifford_codegen::spec::parse_spec;
///
/// let spec = parse_spec(r#"
/// [algebra]
/// name = "euclidean2"
/// complete = false
///
/// [signature]
/// positive = ["e1", "e2"]
///
/// [types.Vector]
/// grades = [1]
/// field_map = [
///   { name = "x", blade = "e1" },
///   { name = "y", blade = "e2" }
/// ]
///
/// [types.Bivector]
/// grades = [2]
/// field_map = [
///   { name = "xy", blade = "e12" }
/// ]
/// "#).unwrap();
///
/// let algebra = Algebra::euclidean(2);
/// let table = ProductTable::new(&algebra);
/// let generator = TraitsGenerator::new(&spec, &algebra, table);
///
/// let (tokens, _tests) = generator.generate_traits_file();
/// let code = tokens.to_string();
///
/// assert!(code.contains("Add for Vector"));
/// ```
pub struct TraitsGenerator<'a> {
    /// The algebra specification.
    spec: &'a AlgebraSpec,
    /// The algebra for computations.
    algebra: &'a Algebra,
    /// The product table for term computation.
    table: ProductTable,
    /// Whether Groebner basis simplification is enabled.
    enable_groebner: bool,
}

impl<'a> TraitsGenerator<'a> {
    /// Creates a new traits generator with default options.
    ///
    /// Groebner simplification is enabled by default.
    pub fn new(spec: &'a AlgebraSpec, algebra: &'a Algebra, table: ProductTable) -> Self {
        Self::with_options(spec, algebra, table, true)
    }

    /// Creates a new traits generator with explicit Groebner option.
    ///
    /// # Arguments
    ///
    /// * `spec` - The algebra specification
    /// * `algebra` - The algebra for computations
    /// * `table` - The product table for term computation
    /// * `enable_groebner` - Whether to enable Groebner basis simplification
    pub fn with_options(
        spec: &'a AlgebraSpec,
        algebra: &'a Algebra,
        table: ProductTable,
        enable_groebner: bool,
    ) -> Self {
        Self {
            spec,
            algebra,
            table,
            enable_groebner,
        }
    }

    /// Generates the complete traits file.
    ///
    /// Returns a tuple of (TokenStream, String) where:
    /// - The TokenStream contains the main code (operators, traits, arbitrary)
    /// - The String contains pre-formatted verification tests that should be
    ///   appended after formatting the main code (rustfmt can't format inside
    ///   proptest! macros)
    pub fn generate_traits_file(&self) -> (TokenStream, String) {
        let header = self.generate_header();
        let imports = self.generate_imports();
        let ops = self.generate_all_ops();
        let product_traits = self.generate_all_product_traits();
        let normed = self.generate_all_normed();
        let approx = self.generate_all_approx();
        let arbitrary = self.generate_all_arbitrary();
        let verification_tests = self.generate_verification_tests_raw();

        let main_tokens = quote! {
            #header
            #imports

            // ============================================================
            // Operator Implementations
            // ============================================================
            #ops

            // ============================================================
            // Product Trait Implementations (clifford::ops)
            // ============================================================
            #product_traits

            // ============================================================
            // Normed Trait Implementations
            // ============================================================
            #normed

            // ============================================================
            // Approx Trait Implementations
            // ============================================================
            #approx

            // ============================================================
            // Arbitrary Implementations (for proptest)
            // ============================================================
            #arbitrary
        };

        (main_tokens, verification_tests)
    }

    /// Generates the file header.
    fn generate_header(&self) -> TokenStream {
        let name = &self.spec.name;
        let header_doc = format!(
            r#"//! Trait implementations for {}.
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

        // Check if any Versor trait impls will actually be generated
        // (versors need Mul output type to match a known type)
        let has_versor_impls = self.will_generate_versor_impls();

        let versor_import = if has_versor_impls {
            quote! { Versor, VersorInverse, InverseSandwich, InverseAntisandwich, }
        } else {
            quote! {}
        };

        // Project/Antiproject are only generated for degenerate algebras (PGA)
        let projection_import = if self.spec.signature.r > 0 {
            quote! { Project, Antiproject, }
        } else {
            quote! {}
        };

        // Import wrapper types based on algebra kind
        let is_degenerate = self.spec.signature.r > 0;
        let wrapper_import = if is_degenerate {
            quote! { use crate::wrappers::Unitized; }
        } else {
            quote! { use crate::wrappers::Unit; }
        };

        quote! {
            use crate::scalar::Float;
            #[allow(unused_imports)]
            use crate::ops::{
                Wedge, Antiwedge, LeftContract, RightContract,
                Sandwich, Antisandwich, Transform, ScalarProduct, BulkContract, WeightContract,
                BulkExpand, WeightExpand, Dot, Antidot, WeightDual,
                Reverse, Antireverse, Involute, RightComplement, #projection_import #versor_import
            };
            use super::types::{#(#type_names),*};
            #[allow(unused_imports)]
            #wrapper_import

            use std::ops::{Add, Sub, Neg, Mul};

            use approx::{AbsDiffEq, RelativeEq, UlpsEq};
        }
    }

    // ========================================================================
    // Operator Implementations
    // ========================================================================

    /// Generates all operator implementations.
    fn generate_all_ops(&self) -> TokenStream {
        let ops: Vec<TokenStream> = self
            .spec
            .types
            .iter()
            .filter(|t| t.alias_of.is_none())
            .flat_map(|ty| self.generate_ops_for_type(ty))
            .collect();

        quote! { #(#ops)* }
    }

    /// Generates operators for a single type.
    #[allow(clippy::vec_init_then_push)]
    fn generate_ops_for_type(&self, ty: &TypeSpec) -> Vec<TokenStream> {
        let mut impls = Vec::new();

        // Same-type operations
        impls.push(self.generate_add(ty));
        impls.push(self.generate_sub(ty));
        impls.push(self.generate_neg(ty));
        impls.push(self.generate_scalar_mul(ty));
        impls.push(self.generate_scalar_mul_reverse_f32(ty));
        impls.push(self.generate_scalar_mul_reverse_f64(ty));

        // Cross-type operations (geometric product) - only for explicit products
        for entry in &self.spec.products.geometric {
            // Only generate if lhs matches this type
            if entry.lhs == ty.name {
                if let Some(other) = self.find_type(&entry.rhs) {
                    if let Some(output_type) = self.find_type(&entry.output) {
                        impls.push(self.generate_geometric_mul_from_entry(
                            ty,
                            other,
                            output_type,
                            entry,
                        ));
                    }
                }
            }
        }

        impls
    }

    /// Finds a TypeSpec by name.
    fn find_type(&self, name: &str) -> Option<&TypeSpec> {
        self.spec.types.iter().find(|t| t.name == name)
    }

    /// Generates a constructor call.
    fn generate_constructor_call(_ty: &TypeSpec, field_exprs: &[TokenStream]) -> TokenStream {
        quote! { Self::new_unchecked(#(#field_exprs),*) }
    }

    /// Generates Add implementation.
    fn generate_add(&self, ty: &TypeSpec) -> TokenStream {
        let name = format_ident!("{}", ty.name);
        let field_adds: Vec<TokenStream> = ty
            .fields
            .iter()
            .map(|f| {
                let fname = format_ident!("{}", f.name);
                quote! { self.#fname() + rhs.#fname() }
            })
            .collect();

        let constructor_call = Self::generate_constructor_call(ty, &field_adds);

        quote! {
            impl<T: Float> Add for #name<T> {
                type Output = Self;

                #[inline]
                fn add(self, rhs: Self) -> Self {
                    #constructor_call
                }
            }
        }
    }

    /// Generates Sub implementation.
    fn generate_sub(&self, ty: &TypeSpec) -> TokenStream {
        let name = format_ident!("{}", ty.name);
        let field_subs: Vec<TokenStream> = ty
            .fields
            .iter()
            .map(|f| {
                let fname = format_ident!("{}", f.name);
                quote! { self.#fname() - rhs.#fname() }
            })
            .collect();

        let constructor_call = Self::generate_constructor_call(ty, &field_subs);

        quote! {
            impl<T: Float> Sub for #name<T> {
                type Output = Self;

                #[inline]
                fn sub(self, rhs: Self) -> Self {
                    #constructor_call
                }
            }
        }
    }

    /// Generates Neg implementation.
    fn generate_neg(&self, ty: &TypeSpec) -> TokenStream {
        let name = format_ident!("{}", ty.name);
        let field_negs: Vec<TokenStream> = ty
            .fields
            .iter()
            .map(|f| {
                let fname = format_ident!("{}", f.name);
                quote! { -self.#fname() }
            })
            .collect();

        let constructor_call = Self::generate_constructor_call(ty, &field_negs);

        quote! {
            impl<T: Float> Neg for #name<T> {
                type Output = Self;

                #[inline]
                fn neg(self) -> Self {
                    #constructor_call
                }
            }
        }
    }

    /// Generates scalar multiplication (Type * T).
    fn generate_scalar_mul(&self, ty: &TypeSpec) -> TokenStream {
        let name = format_ident!("{}", ty.name);

        quote! {
            impl<T: Float> Mul<T> for #name<T> {
                type Output = Self;

                #[inline]
                fn mul(self, scalar: T) -> Self {
                    self.scale(scalar)
                }
            }
        }
    }

    /// Generates reverse scalar multiplication for f32 (f32 * Type).
    fn generate_scalar_mul_reverse_f32(&self, ty: &TypeSpec) -> TokenStream {
        let name = format_ident!("{}", ty.name);

        quote! {
            impl Mul<#name<f32>> for f32 {
                type Output = #name<f32>;

                #[inline]
                fn mul(self, v: #name<f32>) -> #name<f32> {
                    v.scale(self)
                }
            }
        }
    }

    /// Generates reverse scalar multiplication for f64 (f64 * Type).
    fn generate_scalar_mul_reverse_f64(&self, ty: &TypeSpec) -> TokenStream {
        let name = format_ident!("{}", ty.name);

        quote! {
            impl Mul<#name<f64>> for f64 {
                type Output = #name<f64>;

                #[inline]
                fn mul(self, v: #name<f64>) -> #name<f64> {
                    v.scale(self)
                }
            }
        }
    }

    /// Generates geometric product (Type * Other -> Output) from a product entry.
    ///
    /// The formula is computed inline using symbolic simplification, rather than
    /// calling a separate function. This avoids generating geometric_* functions
    /// which are not type-safe in general.
    ///
    /// For self-complementary versors (like Motors in point-based PGA), this uses
    /// the antiproduct formula: `complement(complement(a) * complement(b))`.
    /// This ensures that `M1 * M2` composes transformations correctly.
    fn generate_geometric_mul_from_entry(
        &self,
        a: &TypeSpec,
        b: &TypeSpec,
        output: &TypeSpec,
        _entry: &crate::spec::ProductEntry,
    ) -> TokenStream {
        let a_name = format_ident!("{}", a.name);
        let b_name = format_ident!("{}", b.name);
        let out_name = format_ident!("{}", output.name);

        // Check if both types are self-complementary versors.
        // For such types (like Motor in 3D PGA), the standard geometric product
        // does NOT compose transformations correctly. We need to use the antiproduct:
        // complement(complement(a) * complement(b))
        let a_self_complement = a.versor.is_some()
            && self
                .find_complement_output_type(a)
                .map(|t| t == a.name)
                .unwrap_or(false);
        let b_self_complement = b.versor.is_some()
            && self
                .find_complement_output_type(b)
                .map(|t| t == b.name)
                .unwrap_or(false);

        // For self-complementary versors (like Motor in 3D PGA), use the antigeometric
        // product instead of the geometric product. This ensures M1 * M2 composes
        // transformations correctly in point-based PGA where antisandwich is used.
        let product_kind = if a_self_complement && b_self_complement {
            SymbolicProductKind::Antigeometric
        } else {
            SymbolicProductKind::Geometric
        };

        let field_exprs = self.compute_product_expressions(a, b, output, product_kind);
        let constructor_call = quote! { #out_name::new_unchecked(#(#field_exprs),*) };

        quote! {
            impl<T: Float> Mul<#b_name<T>> for #a_name<T> {
                type Output = #out_name<T>;

                #[inline]
                fn mul(self, rhs: #b_name<T>) -> #out_name<T> {
                    #constructor_call
                }
            }
        }
    }

    /// Computes product expressions using symbolic simplification.
    ///
    /// This is used by Mul operators to inline the formula instead of calling
    /// a separate function.
    ///
    /// The simplification pipeline is:
    /// 1. Apply constraint substitutions (e.g., s*s + xy*xy + ... = 1)
    /// 2. Expand and collect like terms
    /// 3. Apply Groebner basis reduction for constrained types (e.g., Lines with Plücker)
    fn compute_product_expressions(
        &self,
        type_a: &TypeSpec,
        type_b: &TypeSpec,
        output_type: &TypeSpec,
        kind: SymbolicProductKind,
    ) -> Vec<TokenStream> {
        self.compute_product_expressions_with_wrappers(
            type_a,
            None,
            type_b,
            None,
            output_type,
            kind,
        )
    }

    /// Computes product expressions with optional wrapper constraints.
    ///
    /// Like `compute_product_expressions`, but accepts optional wrapper kinds that
    /// add additional constraints to the Groebner basis (e.g., norm_squared = 1 for Unit).
    fn compute_product_expressions_with_wrappers(
        &self,
        type_a: &TypeSpec,
        wrapper_a: Option<WrapperKind>,
        type_b: &TypeSpec,
        wrapper_b: Option<WrapperKind>,
        output_type: &TypeSpec,
        kind: SymbolicProductKind,
    ) -> Vec<TokenStream> {
        let symbolic_product = SymbolicProduct::new(self.algebra);
        let expr_simplifier = ExpressionSimplifier::new();

        // Create Groebner simplifier for constraint-based reduction (with wrapper constraints)
        let groebner_simplifier =
            self.create_groebner_simplifier_with_wrappers(type_a, wrapper_a, type_b, wrapper_b);

        // Create symbolic field variables
        let a_symbols = symbolic_product.create_field_symbols(type_a, "self");
        let b_symbols = symbolic_product.create_field_symbols(type_b, "rhs");

        // Compute symbolic product
        let symbolic_fields =
            symbolic_product.compute(type_a, type_b, output_type, kind, &a_symbols, &b_symbols);

        // Determine which prefixes need .as_inner() access
        let mut wrapped_prefixes = Vec::new();
        if wrapper_a.is_some() {
            wrapped_prefixes.push("self");
        }
        if wrapper_b.is_some() {
            wrapped_prefixes.push("rhs");
        }

        // Create converter for Rust code generation (with wrapper access info)
        let converter =
            AtomToRust::new_with_wrappers(&[type_a, type_b], &["self", "rhs"], &wrapped_prefixes);

        // Apply simplification and convert each field expression
        symbolic_fields
            .iter()
            .map(|field| {
                // Simplify (expand and collect like terms)
                let simplified = expr_simplifier.simplify(&field.expression);
                // Apply Groebner reduction for constrained types
                let reduced = groebner_simplifier.reduce_atom(&simplified);
                converter.convert(&reduced)
            })
            .collect()
    }

    /// Creates a Groebner simplifier with optional wrapper constraints.
    ///
    /// When a wrapper kind is specified for an input type, the corresponding
    /// wrapper constraint is added to the Groebner basis (e.g., norm_squared = 1
    /// for Unit).
    fn create_groebner_simplifier_with_wrappers(
        &self,
        type_a: &TypeSpec,
        wrapper_a: Option<WrapperKind>,
        type_b: &TypeSpec,
        wrapper_b: Option<WrapperKind>,
    ) -> GroebnerSimplifier {
        // If Groebner is disabled, return a disabled simplifier
        if !self.enable_groebner {
            return GroebnerSimplifier::new(vec![], true);
        }

        let collector =
            ProductConstraintCollector::new(self.algebra, self.spec.norm.primary_involution);

        // Collect constraints from both input types, with optional wrapper constraints
        let mut constraints = Vec::new();
        if let Some(wrapper) = wrapper_a {
            constraints.extend(collector.collect_wrapper_constraints(type_a, wrapper, "self"));
        } else {
            constraints.extend(collector.collect_constraints(type_a, "self"));
        }
        if let Some(wrapper) = wrapper_b {
            constraints.extend(collector.collect_wrapper_constraints(type_b, wrapper, "rhs"));
        } else {
            constraints.extend(collector.collect_constraints(type_b, "rhs"));
        }

        // Create Groebner simplifier (uses grevlex ordering for faster computation)
        GroebnerSimplifier::new(constraints, true)
    }

    /// Creates a Groebner simplifier for sandwich products with wrapper constraints.
    ///
    /// For sandwich products, the versor appears twice (v × x × rev(v)), so its
    /// constraint is added once. The operand's constraint is also added if wrapped.
    fn create_groebner_simplifier_for_sandwich(
        &self,
        versor_type: &TypeSpec,
        wrapper_versor: Option<WrapperKind>,
        operand_type: &TypeSpec,
        wrapper_operand: Option<WrapperKind>,
    ) -> GroebnerSimplifier {
        // If Groebner is disabled, return a disabled simplifier
        if !self.enable_groebner {
            return GroebnerSimplifier::new(vec![], true);
        }

        let collector =
            ProductConstraintCollector::new(self.algebra, self.spec.norm.primary_involution);

        // Collect constraints
        let mut constraints = Vec::new();
        if let Some(wrapper) = wrapper_versor {
            constraints.extend(collector.collect_wrapper_constraints(versor_type, wrapper, "self"));
        } else {
            constraints.extend(collector.collect_constraints(versor_type, "self"));
        }
        if let Some(wrapper) = wrapper_operand {
            constraints.extend(collector.collect_wrapper_constraints(
                operand_type,
                wrapper,
                "operand",
            ));
        } else {
            constraints.extend(collector.collect_constraints(operand_type, "operand"));
        }

        GroebnerSimplifier::new(constraints, true)
    }

    /// Computes sandwich product expressions with wrapper constraint optimization.
    ///
    /// Uses symbolic computation and Groebner basis reduction to simplify
    /// the sandwich product when wrapper constraints are present.
    fn compute_sandwich_expressions_with_wrappers(
        &self,
        versor: &TypeSpec,
        wrapper_versor: Option<WrapperKind>,
        operand: &TypeSpec,
        wrapper_operand: Option<WrapperKind>,
        use_antiproduct: bool,
    ) -> Vec<TokenStream> {
        let symbolic_product = SymbolicProduct::new(self.algebra);
        let expr_simplifier = ExpressionSimplifier::new();

        // Create Groebner simplifier with wrapper constraints
        let groebner_simplifier = self.create_groebner_simplifier_for_sandwich(
            versor,
            wrapper_versor,
            operand,
            wrapper_operand,
        );

        // Create symbolic field variables
        let versor_symbols = symbolic_product.create_field_symbols(versor, "self");
        let operand_symbols = symbolic_product.create_field_symbols(operand, "operand");

        // Compute symbolic sandwich product
        let symbolic_fields = symbolic_product.compute_sandwich(
            versor,
            operand,
            &versor_symbols,
            &operand_symbols,
            use_antiproduct,
        );

        // Determine which prefixes need .as_inner() access
        let mut wrapped_prefixes = Vec::new();
        if wrapper_versor.is_some() {
            wrapped_prefixes.push("self");
        }
        if wrapper_operand.is_some() {
            wrapped_prefixes.push("operand");
        }

        // Create converter for Rust code generation
        let converter = AtomToRust::new_with_wrappers(
            &[versor, operand],
            &["self", "operand"],
            &wrapped_prefixes,
        );

        // Apply simplification and convert each field expression
        symbolic_fields
            .iter()
            .map(|field| {
                // Simplify (expand and collect)
                let simplified = expr_simplifier.simplify(&field.expression);
                // Apply Groebner reduction
                let reduced = groebner_simplifier.reduce_atom(&simplified);
                converter.convert(&reduced)
            })
            .collect()
    }

    /// Computes sandwich product expressions: v × x × rev(v).
    ///
    /// Returns TokenStream expressions for each field of the output type.
    fn compute_sandwich_expressions(
        &self,
        versor: &TypeSpec,
        operand: &TypeSpec,
    ) -> Vec<TokenStream> {
        operand
            .fields
            .iter()
            .map(|field| {
                let expr = self.compute_sandwich_field(versor, operand, field.blade_index, false);
                // Apply output field sign for non-canonical blade ordering
                if field.sign < 0 {
                    quote! { -(#expr) }
                } else {
                    expr
                }
            })
            .collect()
    }

    /// Computes antisandwich product expressions: v ⊛ x ⊛ antirev(v).
    ///
    /// Returns TokenStream expressions for each field of the output type.
    fn compute_antisandwich_expressions(
        &self,
        versor: &TypeSpec,
        operand: &TypeSpec,
    ) -> Vec<TokenStream> {
        operand
            .fields
            .iter()
            .map(|field| {
                let expr = self.compute_sandwich_field(versor, operand, field.blade_index, true);
                // Apply output field sign for non-canonical blade ordering
                if field.sign < 0 {
                    quote! { -(#expr) }
                } else {
                    expr
                }
            })
            .collect()
    }

    /// Computes a single sandwich field expression.
    ///
    /// If `use_antiproduct` is true, uses the antiproduct and antireverse (for antisandwich).
    fn compute_sandwich_field(
        &self,
        versor: &TypeSpec,
        operand: &TypeSpec,
        result_blade: usize,
        use_antiproduct: bool,
    ) -> TokenStream {
        let dim = self.algebra.dim();

        // Collect terms: for each combination v_i × x_j × rev(v_k) or v_i ⊛ x_j ⊛ antirev(v_k)
        let mut term_map: std::collections::HashMap<(String, String, String), i8> =
            std::collections::HashMap::new();

        for field_v1 in &versor.fields {
            for field_x in &operand.fields {
                for field_v2 in &versor.fields {
                    let v1_blade = field_v1.blade_index;
                    let x_blade = field_x.blade_index;
                    let v2_blade = field_v2.blade_index;

                    // Compute v_i × x_j (or v_i ⊛ x_j)
                    let (sign_vx, vx) = if use_antiproduct {
                        self.table.antiproduct(v1_blade, x_blade)
                    } else {
                        self.table.geometric(v1_blade, x_blade)
                    };
                    if sign_vx == 0 {
                        continue;
                    }

                    // Compute the reverse/antireverse sign
                    let v2_grade = Blade::from_index(v2_blade).grade();
                    let rev_sign: i8 = if use_antiproduct {
                        // Antireverse sign: (-1)^((n-k)(n-k-1)/2)
                        let antigrade = dim - v2_grade;
                        if (antigrade * antigrade.saturating_sub(1) / 2).is_multiple_of(2) {
                            1
                        } else {
                            -1
                        }
                    } else {
                        // Reverse sign: (-1)^(k(k-1)/2)
                        if (v2_grade * v2_grade.saturating_sub(1) / 2).is_multiple_of(2) {
                            1
                        } else {
                            -1
                        }
                    };

                    // Compute (v_i × x_j) × rev(v_k) (or (v_i ⊛ x_j) ⊛ antirev(v_k))
                    let (sign_vxr, result) = if use_antiproduct {
                        self.table.antiproduct(vx, v2_blade)
                    } else {
                        self.table.geometric(vx, v2_blade)
                    };
                    if sign_vxr == 0 {
                        continue;
                    }

                    if result == result_blade {
                        // Apply input field signs for non-canonical blade orderings.
                        // v1 and v2 are both from versor_type, x is from operand_type.
                        let input_sign = field_v1.sign * field_x.sign * field_v2.sign;
                        let final_sign = sign_vx * sign_vxr * rev_sign * input_sign;
                        let key = (
                            field_v1.name.clone(),
                            field_x.name.clone(),
                            field_v2.name.clone(),
                        );
                        *term_map.entry(key).or_insert(0) += final_sign;
                    }
                }
            }
        }

        // Convert to TokenStream
        if term_map.is_empty() {
            return quote! { T::zero() };
        }

        // Filter out zero coefficients and collect non-zero terms
        let mut terms: Vec<_> = term_map
            .into_iter()
            .filter(|(_, coeff)| *coeff != 0)
            .collect();

        if terms.is_empty() {
            return quote! { T::zero() };
        }

        // Sort for deterministic output
        terms.sort_by(|a, b| a.0.cmp(&b.0));

        let mut expr_parts: Vec<TokenStream> = Vec::new();
        for (i, ((v1, x, v2), coeff)) in terms.iter().enumerate() {
            let v1_ident = format_ident!("{}", v1);
            let x_ident = format_ident!("{}", x);
            let v2_ident = format_ident!("{}", v2);

            let abs_coeff = coeff.abs();
            let is_negative = *coeff < 0;

            // Build the base product expression
            let base_expr = quote! { self.#v1_ident() * operand.#x_ident() * self.#v2_ident() };

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

    /// Computes projection expressions: target ∨ (self ∧ target☆).
    ///
    /// Returns TokenStream expressions for each field of the output type.
    /// Only meaningful for degenerate algebras (PGA) where weight_dual is non-trivial.
    fn compute_project_expressions(
        &self,
        source: &TypeSpec,
        target: &TypeSpec,
        output: &TypeSpec,
    ) -> Vec<TokenStream> {
        output
            .fields
            .iter()
            .map(|field| self.compute_project_field(source, target, field.blade_index))
            .collect()
    }

    /// Computes antiprojection expressions: target ∧ (self ∨ target☆).
    ///
    /// Returns TokenStream expressions for each field of the output type.
    /// Only meaningful for degenerate algebras (PGA) where weight_dual is non-trivial.
    fn compute_antiproject_expressions(
        &self,
        source: &TypeSpec,
        target: &TypeSpec,
        output: &TypeSpec,
    ) -> Vec<TokenStream> {
        output
            .fields
            .iter()
            .map(|field| self.compute_antiproject_field(source, target, field.blade_index))
            .collect()
    }

    /// Computes a single projection field expression.
    ///
    /// For multi-blade types, the projection `B ∨ (A ∧ B☆)` expands to:
    /// `Σ_{i,j,k} coeff_a_i * coeff_b_j * coeff_b_k * (b_k ∨ (a_i ∧ b_j☆))`
    fn compute_project_field(
        &self,
        source: &TypeSpec,
        target: &TypeSpec,
        result_blade: usize,
    ) -> TokenStream {
        // Collect terms: for each combination a_i, b_j (dual), b_k (antiwedge)
        let mut term_map: std::collections::HashMap<(String, String, String), i8> =
            std::collections::HashMap::new();

        for field_a in &source.fields {
            for field_b_dual in &target.fields {
                for field_b_anti in &target.fields {
                    let a_blade = field_a.blade_index;
                    let b_dual_blade = field_b_dual.blade_index;
                    let b_anti_blade = field_b_anti.blade_index;

                    // Compute b_k ∨ (a_i ∧ b_j☆)
                    let (sign, result) =
                        self.table
                            .project_triple(a_blade, b_dual_blade, b_anti_blade);
                    if sign == 0 {
                        continue;
                    }

                    if result == result_blade {
                        let key = (
                            field_a.name.clone(),
                            field_b_dual.name.clone(),
                            field_b_anti.name.clone(),
                        );
                        *term_map.entry(key).or_insert(0) += sign;
                    }
                }
            }
        }

        self.triple_term_map_to_tokens(term_map, "target")
    }

    /// Computes a single antiprojection field expression.
    ///
    /// For multi-blade types, the antiprojection `B ∧ (A ∨ B☆)` expands to:
    /// `Σ_{i,j,k} coeff_a_i * coeff_b_j * coeff_b_k * (b_k ∧ (a_i ∨ b_j☆))`
    fn compute_antiproject_field(
        &self,
        source: &TypeSpec,
        target: &TypeSpec,
        result_blade: usize,
    ) -> TokenStream {
        // Collect terms: for each combination a_i, b_j (dual), b_k (wedge)
        let mut term_map: std::collections::HashMap<(String, String, String), i8> =
            std::collections::HashMap::new();

        for field_a in &source.fields {
            for field_b_dual in &target.fields {
                for field_b_wedge in &target.fields {
                    let a_blade = field_a.blade_index;
                    let b_dual_blade = field_b_dual.blade_index;
                    let b_wedge_blade = field_b_wedge.blade_index;

                    // Compute b_k ∧ (a_i ∨ b_j☆)
                    let (sign, result) =
                        self.table
                            .antiproject_triple(a_blade, b_dual_blade, b_wedge_blade);
                    if sign == 0 {
                        continue;
                    }

                    if result == result_blade {
                        let key = (
                            field_a.name.clone(),
                            field_b_dual.name.clone(),
                            field_b_wedge.name.clone(),
                        );
                        *term_map.entry(key).or_insert(0) += sign;
                    }
                }
            }
        }

        self.triple_term_map_to_tokens(term_map, "target")
    }

    /// Converts a triple term map to TokenStream.
    ///
    /// Used by project/antiproject field computations.
    fn triple_term_map_to_tokens(
        &self,
        term_map: std::collections::HashMap<(String, String, String), i8>,
        rhs_name: &str,
    ) -> TokenStream {
        if term_map.is_empty() {
            return quote! { T::zero() };
        }

        // Filter out zero coefficients and collect non-zero terms
        let mut terms: Vec<_> = term_map
            .into_iter()
            .filter(|(_, coeff)| *coeff != 0)
            .collect();

        if terms.is_empty() {
            return quote! { T::zero() };
        }

        // Sort for deterministic output
        terms.sort_by(|a, b| a.0.cmp(&b.0));

        let rhs_ident = format_ident!("{}", rhs_name);
        let mut expr_parts: Vec<TokenStream> = Vec::new();

        for (i, ((a, b1, b2), coeff)) in terms.iter().enumerate() {
            let a_ident = format_ident!("{}", a);
            let b1_ident = format_ident!("{}", b1);
            let b2_ident = format_ident!("{}", b2);

            let abs_coeff = coeff.abs();
            let is_negative = *coeff < 0;

            // Build the base product expression: self.a * target.b1 * target.b2
            let base_expr =
                quote! { self.#a_ident() * #rhs_ident.#b1_ident() * #rhs_ident.#b2_ident() };

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

    /// Checks if all expressions in a list are `T::zero()`.
    fn all_expressions_are_zero(exprs: &[TokenStream]) -> bool {
        exprs.iter().all(|e| {
            let s = e.to_string();
            // TokenStream string may have various spacing patterns
            s.contains("T :: zero ()") || s.contains("T::zero()") || s == "T :: zero ()"
        })
    }

    /// Computes scalar product expression (grade-0 projection of geometric product).
    ///
    /// Returns TokenStream expression for the scalar result.
    fn compute_scalar_product_expression(&self, a: &TypeSpec, b: &TypeSpec) -> TokenStream {
        let mut terms: Vec<TokenStream> = Vec::new();

        for field_a in &a.fields {
            for field_b in &b.fields {
                let (sign, result) = self
                    .table
                    .geometric(field_a.blade_index, field_b.blade_index);

                // Only include if result is grade 0 (scalar blade index = 0)
                let result_grade = Blade::from_index(result).grade();
                if result_grade != 0 || sign == 0 {
                    continue;
                }

                let a_ident = format_ident!("{}", field_a.name);
                let b_ident = format_ident!("{}", field_b.name);

                let term_expr = if terms.is_empty() {
                    if sign > 0 {
                        quote! { self.#a_ident() * rhs.#b_ident() }
                    } else {
                        quote! { -(self.#a_ident() * rhs.#b_ident()) }
                    }
                } else if sign > 0 {
                    quote! { + self.#a_ident() * rhs.#b_ident() }
                } else {
                    quote! { - self.#a_ident() * rhs.#b_ident() }
                };
                terms.push(term_expr);
            }
        }

        if terms.is_empty() {
            quote! { T::zero() }
        } else {
            quote! { #(#terms)* }
        }
    }

    /// Computes dot product expression (same-grade metric inner product).
    ///
    /// Returns TokenStream expression for the scalar result.
    fn compute_dot_expression(&self, a: &TypeSpec, b: &TypeSpec) -> TokenStream {
        let mut terms: Vec<TokenStream> = Vec::new();

        for field_a in &a.fields {
            for field_b in &b.fields {
                let (sign, _result) = self.table.dot(field_a.blade_index, field_b.blade_index);

                // Only include non-zero results
                if sign == 0 {
                    continue;
                }

                let a_ident = format_ident!("{}", field_a.name);
                let b_ident = format_ident!("{}", field_b.name);

                let term_expr = if terms.is_empty() {
                    if sign > 0 {
                        quote! { self.#a_ident() * rhs.#b_ident() }
                    } else {
                        quote! { -(self.#a_ident() * rhs.#b_ident()) }
                    }
                } else if sign > 0 {
                    quote! { + self.#a_ident() * rhs.#b_ident() }
                } else {
                    quote! { - self.#a_ident() * rhs.#b_ident() }
                };
                terms.push(term_expr);
            }
        }

        if terms.is_empty() {
            quote! { T::zero() }
        } else {
            quote! { #(#terms)* }
        }
    }

    /// Computes antidot product expression (same-antigrade metric anti-inner product).
    ///
    /// Returns TokenStream expression for the scalar result.
    fn compute_antidot_expression(&self, a: &TypeSpec, b: &TypeSpec) -> TokenStream {
        let mut terms: Vec<TokenStream> = Vec::new();

        for field_a in &a.fields {
            for field_b in &b.fields {
                let (sign, _result) = self.table.antidot(field_a.blade_index, field_b.blade_index);

                // Only include non-zero results
                if sign == 0 {
                    continue;
                }

                let a_ident = format_ident!("{}", field_a.name);
                let b_ident = format_ident!("{}", field_b.name);

                let term_expr = if terms.is_empty() {
                    if sign > 0 {
                        quote! { self.#a_ident() * rhs.#b_ident() }
                    } else {
                        quote! { -(self.#a_ident() * rhs.#b_ident()) }
                    }
                } else if sign > 0 {
                    quote! { + self.#a_ident() * rhs.#b_ident() }
                } else {
                    quote! { - self.#a_ident() * rhs.#b_ident() }
                };
                terms.push(term_expr);
            }
        }

        if terms.is_empty() {
            quote! { T::zero() }
        } else {
            quote! { #(#terms)* }
        }
    }

    // ========================================================================
    // Product Trait Implementations (clifford::ops)
    // ========================================================================

    /// Checks if a type represents a single-grade blade.
    ///
    /// Single-grade types have exactly one grade.
    /// Examples: Scalar (grade 0), Vector (grade 1), Bivector (grade 2), Circle (grade 3), etc.
    /// Counter-examples: Rotor (grades 0,2), Motor (grades 0,2,4), Flector (grades 1,3)
    ///
    /// Note: A single-grade blade can also be a versor (e.g., Vector is a reflector,
    /// Circle in CGA is an inversor). Being a blade and being a versor are orthogonal properties.
    fn is_single_grade_blade(&self, ty: &TypeSpec) -> bool {
        ty.grades.len() == 1
    }

    /// Checks if a type has non-zero bulk norm capability.
    ///
    /// A type has non-zero bulk norm if it has at least one field whose blade
    /// does NOT involve any degenerate basis vectors (those with metric == 0).
    ///
    /// This is used to filter out types like `DualUnit` in Cl(0,0,1) where all
    /// blades are purely degenerate, making bulk_norm structurally zero.
    fn has_nonzero_bulk_norm(&self, ty: &TypeSpec) -> bool {
        // Find indices of degenerate basis vectors (those with metric == 0)
        let degenerate_indices: Vec<usize> = self
            .spec
            .signature
            .basis
            .iter()
            .filter(|b| b.metric == 0)
            .map(|b| b.index)
            .collect();

        // If there are no degenerate basis vectors, all types have bulk norm
        if degenerate_indices.is_empty() {
            return true;
        }

        // Check if any field is a "bulk" field (no degenerate basis in its blade)
        ty.fields.iter().any(|f| {
            // Check if this field's blade involves any degenerate basis vector
            !degenerate_indices.iter().any(|&deg_idx| {
                // Check if bit `deg_idx` is set in the blade_index
                (f.blade_index >> deg_idx) & 1 == 1
            })
        })
    }

    /// Generates all product trait implementations.
    fn generate_all_product_traits(&self) -> TokenStream {
        let mut impls = Vec::new();

        // Note: GeometricProduct trait has been removed (PRD-24)
        // The geometric product cannot be type-safe for single-grade elements.
        // Use Dot for same-grade scalar products, or Mul operator for versors.

        // Product traits are only implemented for single-grade blade types.
        // Versors (Motor, Flector, Rotor) use Sandwich/Antisandwich instead.

        // Determine the appropriate wrapper kind for this algebra
        let is_degenerate = self.spec.signature.r > 0;
        let wrapper_kind = if is_degenerate {
            WrapperKind::Unitized // PGA uses unitized blades
        } else {
            WrapperKind::Unit // Euclidean uses unit normalization
        };

        // Wedge trait - single-grade types only
        for entry in &self.spec.products.wedge {
            if let (Some(a), Some(b), Some(out)) = (
                self.find_type(&entry.lhs),
                self.find_type(&entry.rhs),
                self.find_type(&entry.output),
            ) {
                if self.is_single_grade_blade(a) && self.is_single_grade_blade(b) {
                    impls.push(self.generate_wedge_trait(a, b, out, entry));

                    // Generate wrapper variants: Unit<A> ∧ B, A ∧ Unit<B>, Unit<A> ∧ Unit<B>
                    impls.push(self.generate_wrapper_product_trait(
                        a,
                        b,
                        out,
                        SymbolicProductKind::Wedge,
                        wrapper_kind,
                        WrapperPosition::Lhs,
                    ));
                    impls.push(self.generate_wrapper_product_trait(
                        a,
                        b,
                        out,
                        SymbolicProductKind::Wedge,
                        wrapper_kind,
                        WrapperPosition::Rhs,
                    ));
                    impls.push(self.generate_wrapper_product_trait(
                        a,
                        b,
                        out,
                        SymbolicProductKind::Wedge,
                        wrapper_kind,
                        WrapperPosition::Both,
                    ));
                }
            }
        }

        // Antiwedge trait - single-grade types only
        for entry in &self.spec.products.antiwedge {
            if let (Some(a), Some(b), Some(out)) = (
                self.find_type(&entry.lhs),
                self.find_type(&entry.rhs),
                self.find_type(&entry.output),
            ) {
                if self.is_single_grade_blade(a) && self.is_single_grade_blade(b) {
                    impls.push(self.generate_antiwedge_trait(a, b, out, entry));

                    // Generate wrapper variants
                    impls.push(self.generate_wrapper_product_trait(
                        a,
                        b,
                        out,
                        SymbolicProductKind::Antiwedge,
                        wrapper_kind,
                        WrapperPosition::Lhs,
                    ));
                    impls.push(self.generate_wrapper_product_trait(
                        a,
                        b,
                        out,
                        SymbolicProductKind::Antiwedge,
                        wrapper_kind,
                        WrapperPosition::Rhs,
                    ));
                    impls.push(self.generate_wrapper_product_trait(
                        a,
                        b,
                        out,
                        SymbolicProductKind::Antiwedge,
                        wrapper_kind,
                        WrapperPosition::Both,
                    ));
                }
            }
        }

        // LeftContract trait - single-grade types only
        for entry in &self.spec.products.left_contraction {
            if let (Some(a), Some(b), Some(out)) = (
                self.find_type(&entry.lhs),
                self.find_type(&entry.rhs),
                self.find_type(&entry.output),
            ) {
                if self.is_single_grade_blade(a) && self.is_single_grade_blade(b) {
                    impls.push(self.generate_left_contract_trait(a, b, out, entry));

                    // Generate wrapper variants
                    impls.push(self.generate_wrapper_product_trait(
                        a,
                        b,
                        out,
                        SymbolicProductKind::LeftContraction,
                        wrapper_kind,
                        WrapperPosition::Lhs,
                    ));
                    impls.push(self.generate_wrapper_product_trait(
                        a,
                        b,
                        out,
                        SymbolicProductKind::LeftContraction,
                        wrapper_kind,
                        WrapperPosition::Rhs,
                    ));
                    impls.push(self.generate_wrapper_product_trait(
                        a,
                        b,
                        out,
                        SymbolicProductKind::LeftContraction,
                        wrapper_kind,
                        WrapperPosition::Both,
                    ));
                }
            }
        }

        // RightContract trait - single-grade types only
        for entry in &self.spec.products.right_contraction {
            if let (Some(a), Some(b), Some(out)) = (
                self.find_type(&entry.lhs),
                self.find_type(&entry.rhs),
                self.find_type(&entry.output),
            ) {
                if self.is_single_grade_blade(a) && self.is_single_grade_blade(b) {
                    impls.push(self.generate_right_contract_trait(a, b, out, entry));

                    // Generate wrapper variants
                    impls.push(self.generate_wrapper_product_trait(
                        a,
                        b,
                        out,
                        SymbolicProductKind::RightContraction,
                        wrapper_kind,
                        WrapperPosition::Lhs,
                    ));
                    impls.push(self.generate_wrapper_product_trait(
                        a,
                        b,
                        out,
                        SymbolicProductKind::RightContraction,
                        wrapper_kind,
                        WrapperPosition::Rhs,
                    ));
                    impls.push(self.generate_wrapper_product_trait(
                        a,
                        b,
                        out,
                        SymbolicProductKind::RightContraction,
                        wrapper_kind,
                        WrapperPosition::Both,
                    ));
                }
            }
        }

        // Sandwich trait - generated from versor types
        for versor_type in &self.spec.types {
            if versor_type.alias_of.is_some() {
                continue;
            }
            if let Some(ref versor_spec) = versor_type.versor {
                let targets = if versor_spec.sandwich_targets.is_empty() {
                    // Infer valid targets: types where sandwich preserves grades
                    self.infer_sandwich_targets(versor_type)
                } else {
                    versor_spec.sandwich_targets.clone()
                };

                for target_name in &targets {
                    if let Some(target_type) = self.find_type(target_name) {
                        impls.push(
                            self.generate_sandwich_trait_from_versor(versor_type, target_type),
                        );

                        // Generate wrapper variants with Groebner optimization
                        impls.push(self.generate_wrapper_sandwich_trait(
                            versor_type,
                            target_type,
                            wrapper_kind,
                            WrapperPosition::Lhs,
                        ));
                        impls.push(self.generate_wrapper_sandwich_trait(
                            versor_type,
                            target_type,
                            wrapper_kind,
                            WrapperPosition::Rhs,
                        ));
                        impls.push(self.generate_wrapper_sandwich_trait(
                            versor_type,
                            target_type,
                            wrapper_kind,
                            WrapperPosition::Both,
                        ));
                    }
                }
            }
        }

        // Antisandwich trait - generated from versor types (uses same targets as sandwich)
        for versor_type in &self.spec.types {
            if versor_type.alias_of.is_some() {
                continue;
            }
            if let Some(ref versor_spec) = versor_type.versor {
                let targets = if versor_spec.sandwich_targets.is_empty() {
                    // Infer valid targets: types where sandwich preserves grades
                    self.infer_sandwich_targets(versor_type)
                } else {
                    versor_spec.sandwich_targets.clone()
                };

                for target_name in &targets {
                    if let Some(target_type) = self.find_type(target_name) {
                        impls.push(
                            self.generate_antisandwich_trait_from_versor(versor_type, target_type),
                        );

                        // Generate wrapper variants with Groebner optimization
                        impls.push(self.generate_wrapper_antisandwich_trait(
                            versor_type,
                            target_type,
                            wrapper_kind,
                            WrapperPosition::Lhs,
                        ));
                        impls.push(self.generate_wrapper_antisandwich_trait(
                            versor_type,
                            target_type,
                            wrapper_kind,
                            WrapperPosition::Rhs,
                        ));
                        impls.push(self.generate_wrapper_antisandwich_trait(
                            versor_type,
                            target_type,
                            wrapper_kind,
                            WrapperPosition::Both,
                        ));
                    }
                }
            }
        }

        // Transform trait - delegates to Sandwich (non-degenerate) or Antisandwich (degenerate)
        // based on whether the algebra has zero elements in its signature
        for versor_type in &self.spec.types {
            if versor_type.alias_of.is_some() {
                continue;
            }
            if let Some(ref versor_spec) = versor_type.versor {
                let targets = if versor_spec.sandwich_targets.is_empty() {
                    self.infer_sandwich_targets(versor_type)
                } else {
                    versor_spec.sandwich_targets.clone()
                };

                for target_name in &targets {
                    if let Some(target_type) = self.find_type(target_name) {
                        impls.push(
                            self.generate_transform_trait_from_versor(versor_type, target_type),
                        );

                        // Generate wrapper variants
                        impls.push(self.generate_wrapper_transform_trait(
                            versor_type,
                            target_type,
                            wrapper_kind,
                            WrapperPosition::Lhs,
                        ));
                        impls.push(self.generate_wrapper_transform_trait(
                            versor_type,
                            target_type,
                            wrapper_kind,
                            WrapperPosition::Rhs,
                        ));
                        impls.push(self.generate_wrapper_transform_trait(
                            versor_type,
                            target_type,
                            wrapper_kind,
                            WrapperPosition::Both,
                        ));
                    }
                }
            }
        }

        // InverseSandwich trait - generated from versor types
        // Uses inverse instead of reverse: v × x × v⁻¹
        for versor_type in &self.spec.types {
            if versor_type.alias_of.is_some() {
                continue;
            }
            if let Some(ref versor_spec) = versor_type.versor {
                let targets = if versor_spec.sandwich_targets.is_empty() {
                    self.infer_sandwich_targets(versor_type)
                } else {
                    versor_spec.sandwich_targets.clone()
                };

                for target_name in &targets {
                    if let Some(target_type) = self.find_type(target_name) {
                        impls.push(self.generate_inverse_sandwich_trait(versor_type, target_type));
                    }
                }
            }
        }

        // InverseSandwich trait - generated from explicit inverse_sandwich_targets
        // This allows non-versor types (like Circle in CGA) to perform inverse sandwiches
        for source_type in &self.spec.types {
            if source_type.alias_of.is_some() {
                continue;
            }
            // Skip versors (already handled above)
            if source_type.versor.is_some() {
                continue;
            }
            // Only process types with explicit inverse_sandwich_targets
            if source_type.inverse_sandwich_targets.is_empty() {
                continue;
            }

            for target_name in &source_type.inverse_sandwich_targets {
                if let Some(target_type) = self.find_type(target_name) {
                    impls.push(self.generate_inverse_sandwich_trait(source_type, target_type));
                }
            }
        }

        // InverseAntisandwich trait - generated from versor types
        // Uses inverse instead of antireverse: v ⊛ x ⊛ v⁻¹
        for versor_type in &self.spec.types {
            if versor_type.alias_of.is_some() {
                continue;
            }
            if let Some(ref versor_spec) = versor_type.versor {
                let targets = if versor_spec.sandwich_targets.is_empty() {
                    self.infer_sandwich_targets(versor_type)
                } else {
                    versor_spec.sandwich_targets.clone()
                };

                for target_name in &targets {
                    if let Some(target_type) = self.find_type(target_name) {
                        impls.push(
                            self.generate_inverse_antisandwich_trait(versor_type, target_type),
                        );
                    }
                }
            }
        }

        // InverseAntisandwich trait - generated from explicit inverse_sandwich_targets
        // This allows non-versor types (like Circle in CGA) to perform inverse antisandwiches
        for source_type in &self.spec.types {
            if source_type.alias_of.is_some() {
                continue;
            }
            // Skip versors (already handled above)
            if source_type.versor.is_some() {
                continue;
            }
            // Only process types with explicit inverse_sandwich_targets
            if source_type.inverse_sandwich_targets.is_empty() {
                continue;
            }

            for target_name in &source_type.inverse_sandwich_targets {
                if let Some(target_type) = self.find_type(target_name) {
                    impls.push(self.generate_inverse_antisandwich_trait(source_type, target_type));
                }
            }
        }

        // Versor trait - generated for versor types (Rotor, Motor, Flector)
        // Provides compose() method for versor composition
        impls.extend(self.generate_versor_traits());

        // ScalarProduct trait - single-grade types only
        for entry in &self.spec.products.scalar {
            if let (Some(a), Some(b), Some(out)) = (
                self.find_type(&entry.lhs),
                self.find_type(&entry.rhs),
                self.find_type(&entry.output),
            ) {
                if self.is_single_grade_blade(a) && self.is_single_grade_blade(b) {
                    impls.push(self.generate_scalar_product_trait(a, b, out, entry));

                    // Generate wrapper variants (scalar-returning)
                    impls.push(self.generate_wrapper_scalar_returning_product_trait(
                        a,
                        b,
                        SymbolicProductKind::Scalar,
                        wrapper_kind,
                        WrapperPosition::Lhs,
                    ));
                    impls.push(self.generate_wrapper_scalar_returning_product_trait(
                        a,
                        b,
                        SymbolicProductKind::Scalar,
                        wrapper_kind,
                        WrapperPosition::Rhs,
                    ));
                    impls.push(self.generate_wrapper_scalar_returning_product_trait(
                        a,
                        b,
                        SymbolicProductKind::Scalar,
                        wrapper_kind,
                        WrapperPosition::Both,
                    ));
                }
            }
        }

        // BulkContract trait - single-grade types only
        for entry in &self.spec.products.bulk_contraction {
            if let (Some(a), Some(b), Some(out)) = (
                self.find_type(&entry.lhs),
                self.find_type(&entry.rhs),
                self.find_type(&entry.output),
            ) {
                if self.is_single_grade_blade(a) && self.is_single_grade_blade(b) {
                    impls.push(self.generate_bulk_contract_trait(a, b, out, entry));

                    // Generate wrapper variants
                    impls.push(self.generate_wrapper_product_trait(
                        a,
                        b,
                        out,
                        SymbolicProductKind::BulkContraction,
                        wrapper_kind,
                        WrapperPosition::Lhs,
                    ));
                    impls.push(self.generate_wrapper_product_trait(
                        a,
                        b,
                        out,
                        SymbolicProductKind::BulkContraction,
                        wrapper_kind,
                        WrapperPosition::Rhs,
                    ));
                    impls.push(self.generate_wrapper_product_trait(
                        a,
                        b,
                        out,
                        SymbolicProductKind::BulkContraction,
                        wrapper_kind,
                        WrapperPosition::Both,
                    ));
                }
            }
        }

        // WeightContract trait - single-grade types only
        for entry in &self.spec.products.weight_contraction {
            if let (Some(a), Some(b), Some(out)) = (
                self.find_type(&entry.lhs),
                self.find_type(&entry.rhs),
                self.find_type(&entry.output),
            ) {
                if self.is_single_grade_blade(a) && self.is_single_grade_blade(b) {
                    impls.push(self.generate_weight_contract_trait(a, b, out, entry));

                    // Generate wrapper variants
                    impls.push(self.generate_wrapper_product_trait(
                        a,
                        b,
                        out,
                        SymbolicProductKind::WeightContraction,
                        wrapper_kind,
                        WrapperPosition::Lhs,
                    ));
                    impls.push(self.generate_wrapper_product_trait(
                        a,
                        b,
                        out,
                        SymbolicProductKind::WeightContraction,
                        wrapper_kind,
                        WrapperPosition::Rhs,
                    ));
                    impls.push(self.generate_wrapper_product_trait(
                        a,
                        b,
                        out,
                        SymbolicProductKind::WeightContraction,
                        wrapper_kind,
                        WrapperPosition::Both,
                    ));
                }
            }
        }

        // BulkExpand trait - single-grade types only
        for entry in &self.spec.products.bulk_expansion {
            if let (Some(a), Some(b), Some(out)) = (
                self.find_type(&entry.lhs),
                self.find_type(&entry.rhs),
                self.find_type(&entry.output),
            ) {
                if self.is_single_grade_blade(a) && self.is_single_grade_blade(b) {
                    impls.push(self.generate_bulk_expand_trait(a, b, out, entry));

                    // Generate wrapper variants
                    impls.push(self.generate_wrapper_product_trait(
                        a,
                        b,
                        out,
                        SymbolicProductKind::BulkExpansion,
                        wrapper_kind,
                        WrapperPosition::Lhs,
                    ));
                    impls.push(self.generate_wrapper_product_trait(
                        a,
                        b,
                        out,
                        SymbolicProductKind::BulkExpansion,
                        wrapper_kind,
                        WrapperPosition::Rhs,
                    ));
                    impls.push(self.generate_wrapper_product_trait(
                        a,
                        b,
                        out,
                        SymbolicProductKind::BulkExpansion,
                        wrapper_kind,
                        WrapperPosition::Both,
                    ));
                }
            }
        }

        // WeightExpand trait - single-grade types only
        for entry in &self.spec.products.weight_expansion {
            if let (Some(a), Some(b), Some(out)) = (
                self.find_type(&entry.lhs),
                self.find_type(&entry.rhs),
                self.find_type(&entry.output),
            ) {
                if self.is_single_grade_blade(a) && self.is_single_grade_blade(b) {
                    impls.push(self.generate_weight_expand_trait(a, b, out, entry));

                    // Generate wrapper variants
                    impls.push(self.generate_wrapper_product_trait(
                        a,
                        b,
                        out,
                        SymbolicProductKind::WeightExpansion,
                        wrapper_kind,
                        WrapperPosition::Lhs,
                    ));
                    impls.push(self.generate_wrapper_product_trait(
                        a,
                        b,
                        out,
                        SymbolicProductKind::WeightExpansion,
                        wrapper_kind,
                        WrapperPosition::Rhs,
                    ));
                    impls.push(self.generate_wrapper_product_trait(
                        a,
                        b,
                        out,
                        SymbolicProductKind::WeightExpansion,
                        wrapper_kind,
                        WrapperPosition::Both,
                    ));
                }
            }
        }

        // Note: Antigeometric trait has been removed (PRD-24)
        // The antigeometric product cannot be type-safe for single-grade elements.
        // Use Antidot for same-antigrade scalar products.

        // Dot trait (same-grade metric inner product, returns scalar)
        for entry in &self.spec.products.dot {
            if let (Some(a), Some(b)) = (self.find_type(&entry.lhs), self.find_type(&entry.rhs)) {
                impls.push(self.generate_dot_trait(a, b, entry));

                // Generate wrapper variants (scalar-returning)
                impls.push(self.generate_wrapper_scalar_returning_product_trait(
                    a,
                    b,
                    SymbolicProductKind::Dot,
                    wrapper_kind,
                    WrapperPosition::Lhs,
                ));
                impls.push(self.generate_wrapper_scalar_returning_product_trait(
                    a,
                    b,
                    SymbolicProductKind::Dot,
                    wrapper_kind,
                    WrapperPosition::Rhs,
                ));
                impls.push(self.generate_wrapper_scalar_returning_product_trait(
                    a,
                    b,
                    SymbolicProductKind::Dot,
                    wrapper_kind,
                    WrapperPosition::Both,
                ));
            }
        }

        // Antidot trait (same-antigrade metric inner product, returns scalar)
        for entry in &self.spec.products.antidot {
            if let (Some(a), Some(b)) = (self.find_type(&entry.lhs), self.find_type(&entry.rhs)) {
                impls.push(self.generate_antidot_trait(a, b, entry));

                // Generate wrapper variants (scalar-returning)
                impls.push(self.generate_wrapper_scalar_returning_product_trait(
                    a,
                    b,
                    SymbolicProductKind::Antidot,
                    wrapper_kind,
                    WrapperPosition::Lhs,
                ));
                impls.push(self.generate_wrapper_scalar_returning_product_trait(
                    a,
                    b,
                    SymbolicProductKind::Antidot,
                    wrapper_kind,
                    WrapperPosition::Rhs,
                ));
                impls.push(self.generate_wrapper_scalar_returning_product_trait(
                    a,
                    b,
                    SymbolicProductKind::Antidot,
                    wrapper_kind,
                    WrapperPosition::Both,
                ));
            }
        }

        // Project and Antiproject traits - only for degenerate algebras (PGA)
        // The weight dual formula used in these projections requires a degenerate metric.
        // In Euclidean algebras without a degenerate direction, weight_dual is trivially zero.
        if self.spec.signature.r > 0 {
            // Collect single-grade types (excluding Scalar which has grade 0)
            let single_grade_types: Vec<_> = self
                .spec
                .types
                .iter()
                .filter(|t| {
                    t.alias_of.is_none() && self.is_single_grade_blade(t) && !t.grades.contains(&0) // Exclude Scalar
                })
                .collect();

            // Project trait: target ∨ (self ∧ target☆)
            // Output has same grade as source
            for source in &single_grade_types {
                for target in &single_grade_types {
                    let source_grade = source.grades.first().copied().unwrap_or(0);

                    // Find output type with same grade as source
                    if let Some(output) = single_grade_types
                        .iter()
                        .find(|t| t.grades.contains(&source_grade))
                    {
                        // Compute expressions and check if any are non-zero
                        let field_exprs = self.compute_project_expressions(source, target, output);
                        if !Self::all_expressions_are_zero(&field_exprs) {
                            impls.push(self.generate_project_trait(source, target, output));

                            // Generate wrapper variants
                            impls.push(self.generate_wrapper_project_trait(
                                source,
                                target,
                                output,
                                wrapper_kind,
                                WrapperPosition::Lhs,
                            ));
                            impls.push(self.generate_wrapper_project_trait(
                                source,
                                target,
                                output,
                                wrapper_kind,
                                WrapperPosition::Rhs,
                            ));
                            impls.push(self.generate_wrapper_project_trait(
                                source,
                                target,
                                output,
                                wrapper_kind,
                                WrapperPosition::Both,
                            ));
                        }
                    }
                }
            }

            // Antiproject trait: target ∧ (self ∨ target☆)
            // Output has same grade as source
            for source in &single_grade_types {
                for target in &single_grade_types {
                    let source_grade = source.grades.first().copied().unwrap_or(0);

                    // Find output type with same grade as source
                    if let Some(output) = single_grade_types
                        .iter()
                        .find(|t| t.grades.contains(&source_grade))
                    {
                        // Compute expressions and check if any are non-zero
                        let field_exprs =
                            self.compute_antiproject_expressions(source, target, output);
                        if !Self::all_expressions_are_zero(&field_exprs) {
                            impls.push(self.generate_antiproject_trait(source, target, output));

                            // Generate wrapper variants
                            impls.push(self.generate_wrapper_antiproject_trait(
                                source,
                                target,
                                output,
                                wrapper_kind,
                                WrapperPosition::Lhs,
                            ));
                            impls.push(self.generate_wrapper_antiproject_trait(
                                source,
                                target,
                                output,
                                wrapper_kind,
                                WrapperPosition::Rhs,
                            ));
                            impls.push(self.generate_wrapper_antiproject_trait(
                                source,
                                target,
                                output,
                                wrapper_kind,
                                WrapperPosition::Both,
                            ));
                        }
                    }
                }
            }
        }

        // ====== Unary Operation Traits ======

        // Reverse trait - for all types
        for ty in &self.spec.types {
            if ty.alias_of.is_none() {
                impls.push(self.generate_reverse_trait(ty));
            }
        }

        // Antireverse trait - for all types
        for ty in &self.spec.types {
            if ty.alias_of.is_none() {
                impls.push(self.generate_antireverse_trait(ty));
            }
        }

        // Involute trait - for all types (algebra-specific involution for norm)
        for ty in &self.spec.types {
            if ty.alias_of.is_none() {
                impls.push(self.generate_involute_trait(ty));
            }
        }

        // RightComplement trait - only for types that have a complement function
        // (i.e., where the complement grades map to an existing type)
        for ty in &self.spec.types {
            if ty.alias_of.is_none() {
                if let Some(impl_tokens) = self.generate_right_complement_trait(ty) {
                    impls.push(impl_tokens);
                }
            }
        }

        // WeightDual trait - only for types where the weight dual grades map to an existing type
        for ty in &self.spec.types {
            if ty.alias_of.is_none() {
                if let Some(impl_tokens) = self.generate_weight_dual_trait(ty) {
                    impls.push(impl_tokens);
                }
            }
        }

        // VersorInverse trait - for versor types (Rotor, Motor, Flector)
        impls.extend(self.generate_versor_inverse_traits());

        quote! { #(#impls)* }
    }

    // Note: generate_geometric_product_trait has been removed (PRD-24)
    // The GeometricProduct trait cannot be type-safe for single-grade elements.

    /// Generates Dot trait impl (same-grade metric inner product, returns scalar).
    fn generate_dot_trait(
        &self,
        a: &TypeSpec,
        b: &TypeSpec,
        _entry: &crate::spec::ProductEntry,
    ) -> TokenStream {
        let a_name = format_ident!("{}", a.name);
        let b_name = format_ident!("{}", b.name);

        // Compute the dot product expression
        let expr = self.compute_dot_expression(a, b);

        quote! {
            impl<T: Float> Dot<#b_name<T>> for #a_name<T> {
                type Scalar = T;

                #[inline]
                fn dot(&self, rhs: &#b_name<T>) -> T {
                    #expr
                }
            }
        }
    }

    /// Generates Antidot trait impl (same-antigrade metric inner product, returns scalar).
    fn generate_antidot_trait(
        &self,
        a: &TypeSpec,
        b: &TypeSpec,
        _entry: &crate::spec::ProductEntry,
    ) -> TokenStream {
        let a_name = format_ident!("{}", a.name);
        let b_name = format_ident!("{}", b.name);

        // Compute the antidot product expression
        let expr = self.compute_antidot_expression(a, b);

        quote! {
            impl<T: Float> Antidot<#b_name<T>> for #a_name<T> {
                type Scalar = T;

                #[inline]
                fn antidot(&self, rhs: &#b_name<T>) -> T {
                    #expr
                }
            }
        }
    }

    /// Generates Project trait impl: target ∨ (self ∧ target☆).
    ///
    /// Only generated for degenerate algebras (PGA) where weight_dual is meaningful.
    fn generate_project_trait(
        &self,
        source: &TypeSpec,
        target: &TypeSpec,
        output: &TypeSpec,
    ) -> TokenStream {
        let source_name = format_ident!("{}", source.name);
        let target_name = format_ident!("{}", target.name);
        let out_name = format_ident!("{}", output.name);

        // Compute the projection expressions using the triple approach
        let field_exprs = self.compute_project_expressions(source, target, output);

        let constructor_call = quote! { #out_name::new_unchecked(#(#field_exprs),*) };

        quote! {
            impl<T: Float> Project<#target_name<T>> for #source_name<T> {
                type Output = #out_name<T>;

                #[inline]
                fn project(&self, target: &#target_name<T>) -> #out_name<T> {
                    #constructor_call
                }
            }
        }
    }

    /// Generates Antiproject trait impl: target ∧ (self ∨ target☆).
    ///
    /// Only generated for degenerate algebras (PGA) where weight_dual is meaningful.
    fn generate_antiproject_trait(
        &self,
        source: &TypeSpec,
        target: &TypeSpec,
        output: &TypeSpec,
    ) -> TokenStream {
        let source_name = format_ident!("{}", source.name);
        let target_name = format_ident!("{}", target.name);
        let out_name = format_ident!("{}", output.name);

        // Compute the antiprojection expressions using the triple approach
        let field_exprs = self.compute_antiproject_expressions(source, target, output);

        let constructor_call = quote! { #out_name::new_unchecked(#(#field_exprs),*) };

        quote! {
            impl<T: Float> Antiproject<#target_name<T>> for #source_name<T> {
                type Output = #out_name<T>;

                #[inline]
                fn antiproject(&self, target: &#target_name<T>) -> #out_name<T> {
                    #constructor_call
                }
            }
        }
    }

    /// Generates wrapper Project trait impl (e.g., `Project<Target> for Unit<Source>`).
    ///
    /// Project delegates via `.as_inner()` since it's a composite operation.
    fn generate_wrapper_project_trait(
        &self,
        source: &TypeSpec,
        target: &TypeSpec,
        output: &TypeSpec,
        wrapper_kind: WrapperKind,
        wrapper_pos: WrapperPosition,
    ) -> TokenStream {
        let source_name = format_ident!("{}", source.name);
        let target_name = format_ident!("{}", target.name);
        let out_name = format_ident!("{}", output.name);
        let wrapper_name = Self::wrapper_type_name(wrapper_kind);

        match wrapper_pos {
            WrapperPosition::Lhs => {
                quote! {
                    impl<T: Float> Project<#target_name<T>> for #wrapper_name<#source_name<T>> {
                        type Output = #out_name<T>;

                        #[inline]
                        fn project(&self, target: &#target_name<T>) -> #out_name<T> {
                            self.as_inner().project(target)
                        }
                    }
                }
            }
            WrapperPosition::Rhs => {
                quote! {
                    impl<T: Float> Project<#wrapper_name<#target_name<T>>> for #source_name<T> {
                        type Output = #out_name<T>;

                        #[inline]
                        fn project(&self, target: &#wrapper_name<#target_name<T>>) -> #out_name<T> {
                            self.project(target.as_inner())
                        }
                    }
                }
            }
            WrapperPosition::Both => {
                quote! {
                    impl<T: Float> Project<#wrapper_name<#target_name<T>>> for #wrapper_name<#source_name<T>> {
                        type Output = #out_name<T>;

                        #[inline]
                        fn project(&self, target: &#wrapper_name<#target_name<T>>) -> #out_name<T> {
                            self.as_inner().project(target.as_inner())
                        }
                    }
                }
            }
        }
    }

    /// Generates wrapper Antiproject trait impl (e.g., `Antiproject<Target> for Unit<Source>`).
    ///
    /// Antiproject delegates via `.as_inner()` since it's a composite operation.
    fn generate_wrapper_antiproject_trait(
        &self,
        source: &TypeSpec,
        target: &TypeSpec,
        output: &TypeSpec,
        wrapper_kind: WrapperKind,
        wrapper_pos: WrapperPosition,
    ) -> TokenStream {
        let source_name = format_ident!("{}", source.name);
        let target_name = format_ident!("{}", target.name);
        let out_name = format_ident!("{}", output.name);
        let wrapper_name = Self::wrapper_type_name(wrapper_kind);

        match wrapper_pos {
            WrapperPosition::Lhs => {
                quote! {
                    impl<T: Float> Antiproject<#target_name<T>> for #wrapper_name<#source_name<T>> {
                        type Output = #out_name<T>;

                        #[inline]
                        fn antiproject(&self, target: &#target_name<T>) -> #out_name<T> {
                            self.as_inner().antiproject(target)
                        }
                    }
                }
            }
            WrapperPosition::Rhs => {
                quote! {
                    impl<T: Float> Antiproject<#wrapper_name<#target_name<T>>> for #source_name<T> {
                        type Output = #out_name<T>;

                        #[inline]
                        fn antiproject(&self, target: &#wrapper_name<#target_name<T>>) -> #out_name<T> {
                            self.antiproject(target.as_inner())
                        }
                    }
                }
            }
            WrapperPosition::Both => {
                quote! {
                    impl<T: Float> Antiproject<#wrapper_name<#target_name<T>>> for #wrapper_name<#source_name<T>> {
                        type Output = #out_name<T>;

                        #[inline]
                        fn antiproject(&self, target: &#wrapper_name<#target_name<T>>) -> #out_name<T> {
                            self.as_inner().antiproject(target.as_inner())
                        }
                    }
                }
            }
        }
    }

    /// Generates Wedge trait impl.
    fn generate_wedge_trait(
        &self,
        a: &TypeSpec,
        b: &TypeSpec,
        output: &TypeSpec,
        _entry: &crate::spec::ProductEntry,
    ) -> TokenStream {
        let a_name = format_ident!("{}", a.name);
        let b_name = format_ident!("{}", b.name);
        let out_name = format_ident!("{}", output.name);

        // Compute the formula using symbolic machinery
        let field_exprs =
            self.compute_product_expressions(a, b, output, SymbolicProductKind::Wedge);

        // Generate constructor call
        let constructor_call = quote! { #out_name::new_unchecked(#(#field_exprs),*) };

        quote! {
            impl<T: Float> Wedge<#b_name<T>> for #a_name<T> {
                type Output = #out_name<T>;

                #[inline]
                fn wedge(&self, rhs: &#b_name<T>) -> #out_name<T> {
                    #constructor_call
                }
            }
        }
    }

    /// Generates Antiwedge trait impl.
    fn generate_antiwedge_trait(
        &self,
        a: &TypeSpec,
        b: &TypeSpec,
        output: &TypeSpec,
        _entry: &crate::spec::ProductEntry,
    ) -> TokenStream {
        let a_name = format_ident!("{}", a.name);
        let b_name = format_ident!("{}", b.name);
        let out_name = format_ident!("{}", output.name);

        // Compute the formula using symbolic machinery
        let field_exprs =
            self.compute_product_expressions(a, b, output, SymbolicProductKind::Antiwedge);

        // Generate constructor call
        let constructor_call = quote! { #out_name::new_unchecked(#(#field_exprs),*) };

        quote! {
            impl<T: Float> Antiwedge<#b_name<T>> for #a_name<T> {
                type Output = #out_name<T>;

                #[inline]
                fn antiwedge(&self, rhs: &#b_name<T>) -> #out_name<T> {
                    #constructor_call
                }
            }
        }
    }

    /// Generates a wrapper product trait impl (e.g., `Wedge<B> for Unit<A>`).
    ///
    /// This generates product implementations for wrapper types with
    /// constraint-based simplification via Groebner basis.
    fn generate_wrapper_product_trait(
        &self,
        a: &TypeSpec,
        b: &TypeSpec,
        output: &TypeSpec,
        kind: SymbolicProductKind,
        wrapper_kind: WrapperKind,
        wrapper_pos: WrapperPosition,
    ) -> TokenStream {
        let a_name = format_ident!("{}", a.name);
        let b_name = format_ident!("{}", b.name);
        let out_name = format_ident!("{}", output.name);
        let wrapper_name = Self::wrapper_type_name(wrapper_kind);

        // Determine wrapper constraints for each operand
        let (wrapper_a, wrapper_b) = match wrapper_pos {
            WrapperPosition::Lhs => (Some(wrapper_kind), None),
            WrapperPosition::Rhs => (None, Some(wrapper_kind)),
            WrapperPosition::Both => (Some(wrapper_kind), Some(wrapper_kind)),
        };

        // Compute the formula with wrapper constraints
        let field_exprs = self
            .compute_product_expressions_with_wrappers(a, wrapper_a, b, wrapper_b, output, kind);

        // Generate constructor call
        let constructor_call = quote! { #out_name::new_unchecked(#(#field_exprs),*) };

        // Generate the appropriate trait impl based on wrapper position
        // Note: #[allow(unused_variables)] is needed because wrapper constraints can
        // simplify products to constants, making some parameters unused.
        let trait_type = Self::product_trait_type_name(kind);
        let trait_method = Self::product_trait_method_name(kind);
        match wrapper_pos {
            WrapperPosition::Lhs => {
                quote! {
                    #[allow(unused_variables)]
                    impl<T: Float> #trait_type<#b_name<T>> for #wrapper_name<#a_name<T>> {
                        type Output = #out_name<T>;

                        #[inline]
                        fn #trait_method(&self, rhs: &#b_name<T>) -> #out_name<T> {
                            #constructor_call
                        }
                    }
                }
            }
            WrapperPosition::Rhs => {
                quote! {
                    #[allow(unused_variables)]
                    impl<T: Float> #trait_type<#wrapper_name<#b_name<T>>> for #a_name<T> {
                        type Output = #out_name<T>;

                        #[inline]
                        fn #trait_method(&self, rhs: &#wrapper_name<#b_name<T>>) -> #out_name<T> {
                            #constructor_call
                        }
                    }
                }
            }
            WrapperPosition::Both => {
                quote! {
                    #[allow(unused_variables)]
                    impl<T: Float> #trait_type<#wrapper_name<#b_name<T>>> for #wrapper_name<#a_name<T>> {
                        type Output = #out_name<T>;

                        #[inline]
                        fn #trait_method(&self, rhs: &#wrapper_name<#b_name<T>>) -> #out_name<T> {
                            #constructor_call
                        }
                    }
                }
            }
        }
    }

    /// Generates a wrapper product trait impl for scalar-returning products.
    ///
    /// This is used for ScalarProduct, Dot, and Antidot which return `T` directly
    /// with `type Scalar = T` instead of `type Output`.
    fn generate_wrapper_scalar_returning_product_trait(
        &self,
        a: &TypeSpec,
        b: &TypeSpec,
        kind: SymbolicProductKind,
        wrapper_kind: WrapperKind,
        wrapper_pos: WrapperPosition,
    ) -> TokenStream {
        let a_name = format_ident!("{}", a.name);
        let b_name = format_ident!("{}", b.name);
        let wrapper_name = Self::wrapper_type_name(wrapper_kind);

        // Determine wrapper constraints for each operand
        let (wrapper_a, wrapper_b) = match wrapper_pos {
            WrapperPosition::Lhs => (Some(wrapper_kind), None),
            WrapperPosition::Rhs => (None, Some(wrapper_kind)),
            WrapperPosition::Both => (Some(wrapper_kind), Some(wrapper_kind)),
        };

        // For scalar-returning products, we need to get the Scalar type
        let scalar_type = self
            .spec
            .types
            .iter()
            .find(|t| t.grades == vec![0] && t.alias_of.is_none())
            .expect("Scalar type must exist");

        // Compute the formula with wrapper constraints
        let field_exprs = self.compute_product_expressions_with_wrappers(
            a,
            wrapper_a,
            b,
            wrapper_b,
            scalar_type,
            kind,
        );

        // Scalar product returns the scalar directly, not a Scalar type
        let expr = if field_exprs.is_empty() {
            quote! { T::zero() }
        } else {
            field_exprs[0].clone()
        };

        // Generate the appropriate trait impl based on wrapper position
        let trait_type = Self::product_trait_type_name(kind);
        let trait_method = Self::product_trait_method_name(kind);
        match wrapper_pos {
            WrapperPosition::Lhs => {
                quote! {
                    #[allow(unused_variables)]
                    impl<T: Float> #trait_type<#b_name<T>> for #wrapper_name<#a_name<T>> {
                        type Scalar = T;

                        #[inline]
                        fn #trait_method(&self, rhs: &#b_name<T>) -> T {
                            #expr
                        }
                    }
                }
            }
            WrapperPosition::Rhs => {
                quote! {
                    #[allow(unused_variables)]
                    impl<T: Float> #trait_type<#wrapper_name<#b_name<T>>> for #a_name<T> {
                        type Scalar = T;

                        #[inline]
                        fn #trait_method(&self, rhs: &#wrapper_name<#b_name<T>>) -> T {
                            #expr
                        }
                    }
                }
            }
            WrapperPosition::Both => {
                quote! {
                    #[allow(unused_variables)]
                    impl<T: Float> #trait_type<#wrapper_name<#b_name<T>>> for #wrapper_name<#a_name<T>> {
                        type Scalar = T;

                        #[inline]
                        fn #trait_method(&self, rhs: &#wrapper_name<#b_name<T>>) -> T {
                            #expr
                        }
                    }
                }
            }
        }
    }

    /// Returns the Rust identifier for a wrapper type.
    fn wrapper_type_name(wrapper: WrapperKind) -> proc_macro2::Ident {
        let name = match wrapper {
            WrapperKind::Unit => "Unit",
            WrapperKind::Bulk => "Bulk",
            WrapperKind::Unitized => "Unitized",
            WrapperKind::Ideal => "Ideal",
            WrapperKind::Proper => "Proper",
            WrapperKind::Spacelike => "Spacelike",
            WrapperKind::Null => "Null",
        };
        format_ident!("{}", name)
    }

    /// Returns the Rust identifier for a product trait type (capitalized).
    fn product_trait_type_name(kind: SymbolicProductKind) -> proc_macro2::Ident {
        let name = match kind {
            SymbolicProductKind::Wedge => "Wedge",
            SymbolicProductKind::Antiwedge => "Antiwedge",
            SymbolicProductKind::LeftContraction => "LeftContract",
            SymbolicProductKind::RightContraction => "RightContract",
            SymbolicProductKind::Dot => "Dot",
            SymbolicProductKind::Antidot => "Antidot",
            SymbolicProductKind::Scalar => "ScalarProduct",
            SymbolicProductKind::BulkContraction => "BulkContract",
            SymbolicProductKind::WeightContraction => "WeightContract",
            SymbolicProductKind::BulkExpansion => "BulkExpand",
            SymbolicProductKind::WeightExpansion => "WeightExpand",
            _ => "UnknownProduct",
        };
        format_ident!("{}", name)
    }

    /// Returns the Rust identifier for a product trait method (lowercase).
    fn product_trait_method_name(kind: SymbolicProductKind) -> proc_macro2::Ident {
        let name = match kind {
            SymbolicProductKind::Wedge => "wedge",
            SymbolicProductKind::Antiwedge => "antiwedge",
            SymbolicProductKind::LeftContraction => "left_contract",
            SymbolicProductKind::RightContraction => "right_contract",
            SymbolicProductKind::Dot => "dot",
            SymbolicProductKind::Antidot => "antidot",
            SymbolicProductKind::Scalar => "scalar_product",
            SymbolicProductKind::BulkContraction => "bulk_contract",
            SymbolicProductKind::WeightContraction => "weight_contract",
            SymbolicProductKind::BulkExpansion => "bulk_expand",
            SymbolicProductKind::WeightExpansion => "weight_expand",
            _ => "unknown_product",
        };
        format_ident!("{}", name)
    }

    /// Generates LeftContract trait impl.
    fn generate_left_contract_trait(
        &self,
        a: &TypeSpec,
        b: &TypeSpec,
        output: &TypeSpec,
        _entry: &crate::spec::ProductEntry,
    ) -> TokenStream {
        let a_name = format_ident!("{}", a.name);
        let b_name = format_ident!("{}", b.name);
        let out_name = format_ident!("{}", output.name);

        // Compute the formula using symbolic machinery
        let field_exprs =
            self.compute_product_expressions(a, b, output, SymbolicProductKind::LeftContraction);

        // Generate constructor call
        let constructor_call = quote! { #out_name::new_unchecked(#(#field_exprs),*) };

        quote! {
            impl<T: Float> LeftContract<#b_name<T>> for #a_name<T> {
                type Output = #out_name<T>;

                #[inline]
                fn left_contract(&self, rhs: &#b_name<T>) -> #out_name<T> {
                    #constructor_call
                }
            }
        }
    }

    /// Generates RightContract trait impl.
    fn generate_right_contract_trait(
        &self,
        a: &TypeSpec,
        b: &TypeSpec,
        output: &TypeSpec,
        _entry: &crate::spec::ProductEntry,
    ) -> TokenStream {
        let a_name = format_ident!("{}", a.name);
        let b_name = format_ident!("{}", b.name);
        let out_name = format_ident!("{}", output.name);

        // Compute the formula using symbolic machinery
        let field_exprs =
            self.compute_product_expressions(a, b, output, SymbolicProductKind::RightContraction);

        // Generate constructor call
        let constructor_call = quote! { #out_name::new_unchecked(#(#field_exprs),*) };

        quote! {
            impl<T: Float> RightContract<#b_name<T>> for #a_name<T> {
                type Output = #out_name<T>;

                #[inline]
                fn right_contract(&self, rhs: &#b_name<T>) -> #out_name<T> {
                    #constructor_call
                }
            }
        }
    }

    /// Generates Sandwich trait impl from versor type.
    fn generate_sandwich_trait_from_versor(
        &self,
        versor: &TypeSpec,
        operand: &TypeSpec,
    ) -> TokenStream {
        let versor_name = format_ident!("{}", versor.name);
        let operand_name = format_ident!("{}", operand.name);

        // Compute the sandwich expression for each output field
        let field_exprs = self.compute_sandwich_expressions(versor, operand);

        // Generate constructor call
        let constructor_call = quote! { #operand_name::new_unchecked(#(#field_exprs),*) };

        // For sandwich, output is typically same type as operand
        // Allow unused variables for trivial sandwich products (e.g., Scalar on anything)
        quote! {
            #[allow(unused_variables)]
            impl<T: Float> Sandwich<#operand_name<T>> for #versor_name<T> {
                type Output = #operand_name<T>;

                #[inline]
                fn sandwich(&self, operand: &#operand_name<T>) -> #operand_name<T> {
                    #constructor_call
                }
            }
        }
    }

    /// Generates Antisandwich trait impl from versor type.
    fn generate_antisandwich_trait_from_versor(
        &self,
        versor: &TypeSpec,
        operand: &TypeSpec,
    ) -> TokenStream {
        let versor_name = format_ident!("{}", versor.name);
        let operand_name = format_ident!("{}", operand.name);

        // Compute the antisandwich expression for each output field
        let field_exprs = self.compute_antisandwich_expressions(versor, operand);

        // Generate constructor call
        let constructor_call = quote! { #operand_name::new_unchecked(#(#field_exprs),*) };

        // For antisandwich, output is typically same type as operand
        // Allow unused variables for trivial sandwich products (e.g., Scalar on anything)
        quote! {
            #[allow(unused_variables)]
            impl<T: Float> Antisandwich<#operand_name<T>> for #versor_name<T> {
                type Output = #operand_name<T>;

                #[inline]
                fn antisandwich(&self, operand: &#operand_name<T>) -> #operand_name<T> {
                    #constructor_call
                }
            }
        }
    }

    /// Generates InverseSandwich trait impl from versor type.
    ///
    /// InverseSandwich computes `v × x × v⁻¹` where `v⁻¹ = rev(v) / |v|²`.
    /// This is equivalent to `sandwich(x) / |v|²` and correctly handles non-unit versors.
    fn generate_inverse_sandwich_trait(
        &self,
        versor: &TypeSpec,
        operand: &TypeSpec,
    ) -> TokenStream {
        let versor_name = format_ident!("{}", versor.name);
        let operand_name = format_ident!("{}", operand.name);
        let is_degenerate = self.spec.signature.r > 0;

        // Compute the sandwich expression for each output field (same as regular sandwich)
        let field_exprs = self.compute_sandwich_expressions(versor, operand);

        // Scale each field by inv_norm_sq
        let scaled_fields: Vec<TokenStream> = field_exprs
            .iter()
            .map(|expr| quote! { (#expr) * inv_norm_sq })
            .collect();

        // Generate constructor call
        let constructor_call = quote! { #operand_name::new_unchecked(#(#scaled_fields),*) };

        // For degenerate algebras (PGA), use bulk_norm_squared
        // For non-degenerate algebras, use norm_squared
        let norm_computation = if is_degenerate {
            quote! {
                let norm_sq = <Self as crate::norm::DegenerateNormed>::bulk_norm_squared(self);
            }
        } else {
            quote! {
                let norm_sq = <Self as crate::norm::Normed>::norm_squared(self);
            }
        };

        // Allow unused variables for trivial sandwich products (e.g., Scalar on anything)
        quote! {
            #[allow(unused_variables)]
            impl<T: Float> InverseSandwich<#operand_name<T>> for #versor_name<T> {
                type Output = #operand_name<T>;

                #[inline]
                fn try_inverse_sandwich(&self, operand: &#operand_name<T>) -> Option<#operand_name<T>> {
                    #norm_computation
                    if norm_sq.abs() < T::epsilon() {
                        return None;
                    }
                    let inv_norm_sq = T::one() / norm_sq;
                    Some(#constructor_call)
                }
            }
        }
    }

    /// Generates InverseAntisandwich trait impl from versor type.
    ///
    /// InverseAntisandwich computes `v ⊛ x ⊛ v⁻¹` where `v⁻¹ = antirev(v) / |v|²`.
    /// This is equivalent to `antisandwich(x) / |v|²` and correctly handles non-unit versors.
    fn generate_inverse_antisandwich_trait(
        &self,
        versor: &TypeSpec,
        operand: &TypeSpec,
    ) -> TokenStream {
        let versor_name = format_ident!("{}", versor.name);
        let operand_name = format_ident!("{}", operand.name);
        let is_degenerate = self.spec.signature.r > 0;

        // Compute the antisandwich expression for each output field (same as regular antisandwich)
        let field_exprs = self.compute_antisandwich_expressions(versor, operand);

        // Scale each field by inv_norm_sq
        let scaled_fields: Vec<TokenStream> = field_exprs
            .iter()
            .map(|expr| quote! { (#expr) * inv_norm_sq })
            .collect();

        // Generate constructor call
        let constructor_call = quote! { #operand_name::new_unchecked(#(#scaled_fields),*) };

        // For degenerate algebras (PGA), use bulk_norm_squared
        // For non-degenerate algebras, use norm_squared
        let norm_computation = if is_degenerate {
            quote! {
                let norm_sq = <Self as crate::norm::DegenerateNormed>::bulk_norm_squared(self);
            }
        } else {
            quote! {
                let norm_sq = <Self as crate::norm::Normed>::norm_squared(self);
            }
        };

        // Allow unused variables for trivial sandwich products (e.g., Scalar on anything)
        quote! {
            #[allow(unused_variables)]
            impl<T: Float> InverseAntisandwich<#operand_name<T>> for #versor_name<T> {
                type Output = #operand_name<T>;

                #[inline]
                fn try_inverse_antisandwich(&self, operand: &#operand_name<T>) -> Option<#operand_name<T>> {
                    #norm_computation
                    if norm_sq.abs() < T::epsilon() {
                        return None;
                    }
                    let inv_norm_sq = T::one() / norm_sq;
                    Some(#constructor_call)
                }
            }
        }
    }

    /// Generates wrapper Sandwich trait impl with optimized computation.
    ///
    /// Uses Groebner basis reduction to simplify sandwich products when
    /// wrapper constraints are present (e.g., Unit<Motor> with |motor|² = 1).
    fn generate_wrapper_sandwich_trait(
        &self,
        versor: &TypeSpec,
        operand: &TypeSpec,
        wrapper_kind: WrapperKind,
        wrapper_pos: WrapperPosition,
    ) -> TokenStream {
        let versor_name = format_ident!("{}", versor.name);
        let operand_name = format_ident!("{}", operand.name);
        let wrapper_name = Self::wrapper_type_name(wrapper_kind);

        // Determine wrapper constraints
        let (wrapper_versor, wrapper_operand) = match wrapper_pos {
            WrapperPosition::Lhs => (Some(wrapper_kind), None),
            WrapperPosition::Rhs => (None, Some(wrapper_kind)),
            WrapperPosition::Both => (Some(wrapper_kind), Some(wrapper_kind)),
        };

        // Compute optimized expressions
        let field_exprs = self.compute_sandwich_expressions_with_wrappers(
            versor,
            wrapper_versor,
            operand,
            wrapper_operand,
            false, // sandwich uses geometric product
        );

        let constructor_call = quote! { #operand_name::new_unchecked(#(#field_exprs),*) };

        match wrapper_pos {
            WrapperPosition::Lhs => {
                quote! {
                    #[allow(unused_variables)]
                    impl<T: Float> Sandwich<#operand_name<T>> for #wrapper_name<#versor_name<T>> {
                        type Output = #operand_name<T>;

                        #[inline]
                        fn sandwich(&self, operand: &#operand_name<T>) -> #operand_name<T> {
                            #constructor_call
                        }
                    }
                }
            }
            WrapperPosition::Rhs => {
                quote! {
                    #[allow(unused_variables)]
                    impl<T: Float> Sandwich<#wrapper_name<#operand_name<T>>> for #versor_name<T> {
                        type Output = #operand_name<T>;

                        #[inline]
                        fn sandwich(&self, operand: &#wrapper_name<#operand_name<T>>) -> #operand_name<T> {
                            #constructor_call
                        }
                    }
                }
            }
            WrapperPosition::Both => {
                quote! {
                    #[allow(unused_variables)]
                    impl<T: Float> Sandwich<#wrapper_name<#operand_name<T>>> for #wrapper_name<#versor_name<T>> {
                        type Output = #operand_name<T>;

                        #[inline]
                        fn sandwich(&self, operand: &#wrapper_name<#operand_name<T>>) -> #operand_name<T> {
                            #constructor_call
                        }
                    }
                }
            }
        }
    }

    /// Generates wrapper Antisandwich trait impl with optimized computation.
    ///
    /// Uses Groebner basis reduction to simplify antisandwich products when
    /// wrapper constraints are present.
    fn generate_wrapper_antisandwich_trait(
        &self,
        versor: &TypeSpec,
        operand: &TypeSpec,
        wrapper_kind: WrapperKind,
        wrapper_pos: WrapperPosition,
    ) -> TokenStream {
        let versor_name = format_ident!("{}", versor.name);
        let operand_name = format_ident!("{}", operand.name);
        let wrapper_name = Self::wrapper_type_name(wrapper_kind);

        // Determine wrapper constraints
        let (wrapper_versor, wrapper_operand) = match wrapper_pos {
            WrapperPosition::Lhs => (Some(wrapper_kind), None),
            WrapperPosition::Rhs => (None, Some(wrapper_kind)),
            WrapperPosition::Both => (Some(wrapper_kind), Some(wrapper_kind)),
        };

        // Compute optimized expressions
        let field_exprs = self.compute_sandwich_expressions_with_wrappers(
            versor,
            wrapper_versor,
            operand,
            wrapper_operand,
            true, // antisandwich uses antiproduct
        );

        let constructor_call = quote! { #operand_name::new_unchecked(#(#field_exprs),*) };

        match wrapper_pos {
            WrapperPosition::Lhs => {
                quote! {
                    #[allow(unused_variables)]
                    impl<T: Float> Antisandwich<#operand_name<T>> for #wrapper_name<#versor_name<T>> {
                        type Output = #operand_name<T>;

                        #[inline]
                        fn antisandwich(&self, operand: &#operand_name<T>) -> #operand_name<T> {
                            #constructor_call
                        }
                    }
                }
            }
            WrapperPosition::Rhs => {
                quote! {
                    #[allow(unused_variables)]
                    impl<T: Float> Antisandwich<#wrapper_name<#operand_name<T>>> for #versor_name<T> {
                        type Output = #operand_name<T>;

                        #[inline]
                        fn antisandwich(&self, operand: &#wrapper_name<#operand_name<T>>) -> #operand_name<T> {
                            #constructor_call
                        }
                    }
                }
            }
            WrapperPosition::Both => {
                quote! {
                    #[allow(unused_variables)]
                    impl<T: Float> Antisandwich<#wrapper_name<#operand_name<T>>> for #wrapper_name<#versor_name<T>> {
                        type Output = #operand_name<T>;

                        #[inline]
                        fn antisandwich(&self, operand: &#wrapper_name<#operand_name<T>>) -> #operand_name<T> {
                            #constructor_call
                        }
                    }
                }
            }
        }
    }

    /// Generates wrapper Transform trait impl.
    ///
    /// Transform delegates to Sandwich (non-degenerate) or Antisandwich (degenerate).
    fn generate_wrapper_transform_trait(
        &self,
        versor: &TypeSpec,
        operand: &TypeSpec,
        wrapper_kind: WrapperKind,
        wrapper_pos: WrapperPosition,
    ) -> TokenStream {
        let versor_name = format_ident!("{}", versor.name);
        let operand_name = format_ident!("{}", operand.name);
        let wrapper_name = Self::wrapper_type_name(wrapper_kind);

        let is_degenerate = self.spec.signature.r > 0;
        let method_name = if is_degenerate {
            quote! { antisandwich }
        } else {
            quote! { sandwich }
        };

        match wrapper_pos {
            WrapperPosition::Lhs => {
                quote! {
                    impl<T: Float> Transform<#operand_name<T>> for #wrapper_name<#versor_name<T>> {
                        type Output = #operand_name<T>;

                        #[inline]
                        fn transform(&self, operand: &#operand_name<T>) -> #operand_name<T> {
                            self.#method_name(operand)
                        }
                    }
                }
            }
            WrapperPosition::Rhs => {
                quote! {
                    impl<T: Float> Transform<#wrapper_name<#operand_name<T>>> for #versor_name<T> {
                        type Output = #operand_name<T>;

                        #[inline]
                        fn transform(&self, operand: &#wrapper_name<#operand_name<T>>) -> #operand_name<T> {
                            self.#method_name(operand)
                        }
                    }
                }
            }
            WrapperPosition::Both => {
                quote! {
                    impl<T: Float> Transform<#wrapper_name<#operand_name<T>>> for #wrapper_name<#versor_name<T>> {
                        type Output = #operand_name<T>;

                        #[inline]
                        fn transform(&self, operand: &#wrapper_name<#operand_name<T>>) -> #operand_name<T> {
                            self.#method_name(operand)
                        }
                    }
                }
            }
        }
    }

    /// Generates Transform trait impl from versor type.
    ///
    /// The Transform trait delegates to either Sandwich or Antisandwich based on
    /// whether the algebra has degenerate elements (zero basis vectors):
    /// - Non-degenerate (r = 0): uses Sandwich
    /// - Degenerate (r > 0): uses Antisandwich
    fn generate_transform_trait_from_versor(
        &self,
        versor: &TypeSpec,
        operand: &TypeSpec,
    ) -> TokenStream {
        let versor_name = format_ident!("{}", versor.name);
        let operand_name = format_ident!("{}", operand.name);

        // Check if algebra is degenerate (has zero elements in signature)
        let is_degenerate = self.spec.signature.r > 0;

        let method_call = if is_degenerate {
            quote! { self.antisandwich(operand) }
        } else {
            quote! { self.sandwich(operand) }
        };

        quote! {
            impl<T: Float> Transform<#operand_name<T>> for #versor_name<T> {
                type Output = #operand_name<T>;

                #[inline]
                fn transform(&self, operand: &#operand_name<T>) -> #operand_name<T> {
                    #method_call
                }
            }
        }
    }

    /// Generates Versor trait impls for all versor×versor combinations.
    ///
    /// The Versor trait provides `compose()` for versor composition.
    /// This delegates to the Mul operator which already implements the
    /// correct product formula (geometric for Euclidean, antigeometric for PGA).
    ///
    /// Output types follow the algebraic rules:
    /// - Even × Even → Even (Motor × Motor → Motor)
    /// - Odd × Odd → Even (Flector × Flector → Motor)
    /// - Even × Odd → Odd (Motor × Flector → Flector)
    /// - Odd × Even → Odd (Flector × Motor → Flector)
    fn generate_versor_traits(&self) -> Vec<TokenStream> {
        let mut impls = Vec::new();

        // Find all versor types
        let versor_types: Vec<_> = self
            .spec
            .types
            .iter()
            .filter(|t| t.alias_of.is_none() && t.versor.is_some())
            .collect();

        // Generate impl for each pair of versors
        for lhs in &versor_types {
            for rhs in &versor_types {
                // Look up the output type from the geometric products
                // (the Mul operator output type)
                if let Some(output_type) = self.find_mul_output_type(&lhs.name, &rhs.name) {
                    let lhs_name = format_ident!("{}", lhs.name);
                    let rhs_name = format_ident!("{}", rhs.name);
                    let out_name = format_ident!("{}", output_type);

                    // Versor composition uses the Mul operator, which is:
                    // - Antigeometric product for self-complementary versors (Motor × Motor)
                    // - Geometric product for other types
                    // Since the Mul operator is already configured correctly, compose just delegates.
                    let compose_body = quote! {
                        *self * *other
                    };

                    impls.push(quote! {
                        impl<T: Float> Versor<#rhs_name<T>> for #lhs_name<T> {
                            type Output = #out_name<T>;

                            #[inline]
                            fn compose(&self, other: &#rhs_name<T>) -> #out_name<T> {
                                #compose_body
                            }
                        }
                    });
                }
            }
        }

        impls
    }

    /// Generates VersorInverse trait impls for all versor types and types with inverse_sandwich_targets.
    ///
    /// The VersorInverse trait provides `try_inverse()` for computing the
    /// multiplicative inverse of a versor: `V⁻¹ = rev(V) / |V|²`.
    ///
    /// For non-degenerate algebras (Euclidean), this uses the standard norm.
    /// For degenerate algebras (PGA), this uses the bulk norm.
    fn generate_versor_inverse_traits(&self) -> Vec<TokenStream> {
        let mut impls = Vec::new();
        let is_degenerate = self.spec.signature.r > 0;

        // Find all types that need VersorInverse:
        // 1. Versor types (Rotor, Motor, Flector)
        // 2. Non-versor types with explicit inverse_sandwich_targets (like Circle in CGA)
        let types_needing_inverse: Vec<_> = self
            .spec
            .types
            .iter()
            .filter(|t| {
                t.alias_of.is_none()
                    && (t.versor.is_some() || !t.inverse_sandwich_targets.is_empty())
            })
            .collect();

        for ty in types_needing_inverse {
            let type_name = format_ident!("{}", ty.name);

            // Compute scaled fields: rev(V) / norm_squared
            let scaled_fields: Vec<TokenStream> = ty
                .fields
                .iter()
                .map(|field| {
                    let field_name = format_ident!("{}", field.name);
                    let grade = field.grade;
                    // Reverse sign: (-1)^(k(k-1)/2)
                    if (grade * grade.saturating_sub(1) / 2).is_multiple_of(2) {
                        quote! { self.#field_name() * inv_norm_sq }
                    } else {
                        quote! { -self.#field_name() * inv_norm_sq }
                    }
                })
                .collect();

            // For degenerate algebras (PGA), use bulk_norm_squared
            // For non-degenerate algebras, use norm_squared
            let norm_computation = if is_degenerate {
                quote! {
                    let norm_sq = <Self as crate::norm::DegenerateNormed>::bulk_norm_squared(self);
                }
            } else {
                quote! {
                    let norm_sq = <Self as crate::norm::Normed>::norm_squared(self);
                }
            };

            impls.push(quote! {
                impl<T: Float> VersorInverse for #type_name<T> {
                    fn try_inverse(&self) -> Option<Self> {
                        #norm_computation
                        if norm_sq.abs() < T::epsilon() {
                            return None;
                        }
                        let inv_norm_sq = T::one() / norm_sq;
                        Some(Self::new_unchecked(#(#scaled_fields),*))
                    }
                }
            });
        }

        impls
    }

    /// Finds the output type for a Mul<Rhs> for Lhs operation.
    fn find_mul_output_type(&self, lhs: &str, rhs: &str) -> Option<String> {
        // Check geometric products for the output type
        for entry in &self.spec.products.geometric {
            if entry.lhs == lhs && entry.rhs == rhs {
                return Some(entry.output.clone());
            }
        }
        None
    }

    /// Checks if any Versor trait impls will be generated.
    ///
    /// Returns true if there's at least one pair of versor types where
    /// their geometric product produces a known output type.
    fn will_generate_versor_impls(&self) -> bool {
        let versor_types: Vec<_> = self
            .spec
            .types
            .iter()
            .filter(|t| t.alias_of.is_none() && t.versor.is_some())
            .collect();

        for lhs in &versor_types {
            for rhs in &versor_types {
                if self.find_mul_output_type(&lhs.name, &rhs.name).is_some() {
                    return true;
                }
            }
        }
        false
    }

    /// Generates ScalarProduct trait impl.
    fn generate_scalar_product_trait(
        &self,
        a: &TypeSpec,
        b: &TypeSpec,
        _output: &TypeSpec,
        _entry: &crate::spec::ProductEntry,
    ) -> TokenStream {
        let a_name = format_ident!("{}", a.name);
        let b_name = format_ident!("{}", b.name);

        // Compute the scalar product expression (grade-0 projection)
        let expr = self.compute_scalar_product_expression(a, b);

        quote! {
            impl<T: Float> ScalarProduct<#b_name<T>> for #a_name<T> {
                type Scalar = T;

                #[inline]
                fn scalar_product(&self, rhs: &#b_name<T>) -> T {
                    #expr
                }
            }
        }
    }

    /// Generates BulkContract trait impl.
    fn generate_bulk_contract_trait(
        &self,
        a: &TypeSpec,
        b: &TypeSpec,
        output: &TypeSpec,
        _entry: &crate::spec::ProductEntry,
    ) -> TokenStream {
        let a_name = format_ident!("{}", a.name);
        let b_name = format_ident!("{}", b.name);
        let out_name = format_ident!("{}", output.name);

        // Compute the formula using symbolic machinery
        let field_exprs =
            self.compute_product_expressions(a, b, output, SymbolicProductKind::BulkContraction);

        // Generate constructor call
        let constructor_call = quote! { #out_name::new_unchecked(#(#field_exprs),*) };

        quote! {
            impl<T: Float> BulkContract<#b_name<T>> for #a_name<T> {
                type Output = #out_name<T>;

                #[inline]
                fn bulk_contract(&self, rhs: &#b_name<T>) -> #out_name<T> {
                    #constructor_call
                }
            }
        }
    }

    /// Generates WeightContract trait impl.
    fn generate_weight_contract_trait(
        &self,
        a: &TypeSpec,
        b: &TypeSpec,
        output: &TypeSpec,
        _entry: &crate::spec::ProductEntry,
    ) -> TokenStream {
        let a_name = format_ident!("{}", a.name);
        let b_name = format_ident!("{}", b.name);
        let out_name = format_ident!("{}", output.name);

        // Compute the formula using symbolic machinery
        let field_exprs =
            self.compute_product_expressions(a, b, output, SymbolicProductKind::WeightContraction);

        // Generate constructor call
        let constructor_call = quote! { #out_name::new_unchecked(#(#field_exprs),*) };

        quote! {
            impl<T: Float> WeightContract<#b_name<T>> for #a_name<T> {
                type Output = #out_name<T>;

                #[inline]
                fn weight_contract(&self, rhs: &#b_name<T>) -> #out_name<T> {
                    #constructor_call
                }
            }
        }
    }

    /// Generates BulkExpand trait impl.
    fn generate_bulk_expand_trait(
        &self,
        a: &TypeSpec,
        b: &TypeSpec,
        output: &TypeSpec,
        _entry: &crate::spec::ProductEntry,
    ) -> TokenStream {
        let a_name = format_ident!("{}", a.name);
        let b_name = format_ident!("{}", b.name);
        let out_name = format_ident!("{}", output.name);

        // Compute the formula using symbolic machinery
        let field_exprs =
            self.compute_product_expressions(a, b, output, SymbolicProductKind::BulkExpansion);

        // Generate constructor call
        let constructor_call = quote! { #out_name::new_unchecked(#(#field_exprs),*) };

        quote! {
            impl<T: Float> BulkExpand<#b_name<T>> for #a_name<T> {
                type Output = #out_name<T>;

                #[inline]
                fn bulk_expand(&self, rhs: &#b_name<T>) -> #out_name<T> {
                    #constructor_call
                }
            }
        }
    }

    /// Generates WeightExpand trait impl.
    fn generate_weight_expand_trait(
        &self,
        a: &TypeSpec,
        b: &TypeSpec,
        output: &TypeSpec,
        _entry: &crate::spec::ProductEntry,
    ) -> TokenStream {
        let a_name = format_ident!("{}", a.name);
        let b_name = format_ident!("{}", b.name);
        let out_name = format_ident!("{}", output.name);

        // Compute the formula using symbolic machinery
        let field_exprs =
            self.compute_product_expressions(a, b, output, SymbolicProductKind::WeightExpansion);

        // Generate constructor call
        let constructor_call = quote! { #out_name::new_unchecked(#(#field_exprs),*) };

        quote! {
            impl<T: Float> WeightExpand<#b_name<T>> for #a_name<T> {
                type Output = #out_name<T>;

                #[inline]
                fn weight_expand(&self, rhs: &#b_name<T>) -> #out_name<T> {
                    #constructor_call
                }
            }
        }
    }

    // Note: generate_antigeometric_trait has been removed (PRD-24)
    // The Antigeometric trait cannot be type-safe for single-grade elements.

    // Note: generate_project_trait and generate_antiproject_trait have been removed.
    // The formulas b ∨ (a ∧ b☆) and b ∧ (a ∨ b☆) are NOT bilinear because 'b' appears
    // twice. The blade-pair product approach doesn't work for these compound operations.
    // Use Multivector::project() and Multivector::antiproject() instead.

    // ========================================================================
    // Sandwich Target Inference
    // ========================================================================

    /// Infers valid sandwich targets for a versor type.
    ///
    /// A type is a valid target if the sandwich product V * X * rev(V) produces
    /// the same grades as X (grade-preserving transformation).
    fn infer_sandwich_targets(&self, _versor_type: &TypeSpec) -> Vec<String> {
        // Versors have the grade-preserving property: V * X * rev(V) preserves the grade of X.
        // This means ANY type can be a valid sandwich target for a versor.
        // We include all non-alias types as valid targets.
        //
        // Note: This is a fundamental property of versors in geometric algebra.
        // A versor is a product of unit vectors, and conjugation by a versor
        // preserves grades (it only rotates/reflects within each grade subspace).
        self.spec
            .types
            .iter()
            .filter(|t| t.alias_of.is_none())
            .map(|t| t.name.clone())
            .collect()
    }

    // ========================================================================
    // Unary Operation Trait Implementations
    // ========================================================================

    /// Generates Reverse trait impl.
    fn generate_reverse_trait(&self, ty: &TypeSpec) -> TokenStream {
        let type_name = format_ident!("{}", ty.name);

        // Compute field expressions with reverse signs
        let field_exprs: Vec<TokenStream> = ty
            .fields
            .iter()
            .map(|field| {
                let field_name = format_ident!("{}", field.name);
                let grade = field.grade;
                // Reverse sign: (-1)^(k(k-1)/2)
                if (grade * grade.saturating_sub(1) / 2).is_multiple_of(2) {
                    quote! { self.#field_name() }
                } else {
                    quote! { -self.#field_name() }
                }
            })
            .collect();

        let constructor = quote! { Self::new_unchecked(#(#field_exprs),*) };

        quote! {
            impl<T: Float> Reverse for #type_name<T> {
                #[inline]
                fn reverse(&self) -> Self {
                    #constructor
                }
            }
        }
    }

    /// Generates Antireverse trait impl.
    fn generate_antireverse_trait(&self, ty: &TypeSpec) -> TokenStream {
        let type_name = format_ident!("{}", ty.name);
        let dim = self.algebra.dim();

        // Compute field expressions with antireverse signs
        let field_exprs: Vec<TokenStream> = ty
            .fields
            .iter()
            .map(|field| {
                let field_name = format_ident!("{}", field.name);
                let grade = field.grade;
                let antigrade = dim - grade;
                // Antireverse sign: (-1)^((n-k)(n-k-1)/2)
                if (antigrade * antigrade.saturating_sub(1) / 2).is_multiple_of(2) {
                    quote! { self.#field_name() }
                } else {
                    quote! { -self.#field_name() }
                }
            })
            .collect();

        let constructor = quote! { Self::new_unchecked(#(#field_exprs),*) };

        quote! {
            impl<T: Float> Antireverse for #type_name<T> {
                #[inline]
                fn antireverse(&self) -> Self {
                    #constructor
                }
            }
        }
    }

    /// Generates Involute trait impl (algebra-specific norm involution).
    ///
    /// The involution used depends on the algebra's `norm.primary_involution` setting:
    /// - Reverse: `(-1)^(k(k-1)/2)` for grade k
    /// - GradeInvolution: `(-1)^k` for grade k
    /// - CliffordConjugate: `(-1)^(k(k+1)/2)` for grade k
    fn generate_involute_trait(&self, ty: &TypeSpec) -> TokenStream {
        let type_name = format_ident!("{}", ty.name);
        let involution_kind = self.spec.norm.primary_involution;

        // Compute field expressions with appropriate involution signs
        let field_exprs: Vec<TokenStream> = ty
            .fields
            .iter()
            .map(|field| {
                let field_name = format_ident!("{}", field.name);
                let grade = field.grade;

                // Compute sign based on involution kind
                let is_positive = match involution_kind {
                    InvolutionKind::Reverse => {
                        // Reverse sign: (-1)^(k(k-1)/2)
                        (grade * grade.saturating_sub(1) / 2).is_multiple_of(2)
                    }
                    InvolutionKind::GradeInvolution => {
                        // Grade involution sign: (-1)^k
                        grade.is_multiple_of(2)
                    }
                    InvolutionKind::CliffordConjugate => {
                        // Clifford conjugate sign: (-1)^(k(k+1)/2)
                        (grade * (grade + 1) / 2).is_multiple_of(2)
                    }
                };

                if is_positive {
                    quote! { self.#field_name() }
                } else {
                    quote! { -self.#field_name() }
                }
            })
            .collect();

        let constructor = quote! { Self::new_unchecked(#(#field_exprs),*) };

        quote! {
            impl<T: Float> Involute for #type_name<T> {
                #[inline]
                fn involute(&self) -> Self {
                    #constructor
                }
            }
        }
    }

    /// Generates RightComplement trait impl.
    ///
    /// Returns None if the complement grades don't map to an existing type.
    fn generate_right_complement_trait(&self, ty: &TypeSpec) -> Option<TokenStream> {
        // Check if there's a matching output type for the complement
        let output_type_name = self.find_complement_output_type(ty)?;
        let output_type = self
            .spec
            .types
            .iter()
            .find(|t| t.name == output_type_name)?;

        let type_name = format_ident!("{}", ty.name);
        let out_name = format_ident!("{}", output_type_name);

        // Compute field expressions for complement
        let field_exprs: Vec<TokenStream> = output_type
            .fields
            .iter()
            .map(|out_field| {
                // Find the input field that complements to this output blade
                let out_blade = out_field.blade_index;

                for in_field in &ty.fields {
                    let (sign, comp_blade) = self.table.complement(in_field.blade_index);
                    if comp_blade == out_blade && sign != 0 {
                        let in_name = format_ident!("{}", in_field.name);
                        return if sign > 0 {
                            quote! { self.#in_name() }
                        } else {
                            quote! { -self.#in_name() }
                        };
                    }
                }
                // No input blade maps to this output blade
                quote! { T::zero() }
            })
            .collect();

        let constructor = quote! { #out_name::new_unchecked(#(#field_exprs),*) };

        Some(quote! {
            impl<T: Float> RightComplement for #type_name<T> {
                type Output = #out_name<T>;

                #[inline]
                fn right_complement(&self) -> #out_name<T> {
                    #constructor
                }
            }
        })
    }

    // Note: LeftComplement, BulkDual, and WeightDual trait generation removed
    // because the corresponding free functions don't exist in unary.rs yet.
    // These can be added in a future PR by extending unary.rs.

    /// Finds the output type for complement operations.
    ///
    /// The complement of a type with grades [g1, g2, ...] has grades [n-g1, n-g2, ...].
    /// Returns None if no type with matching grades exists.
    fn find_complement_output_type(&self, ty: &TypeSpec) -> Option<String> {
        let dim = self.algebra.dim();
        let complement_grades: Vec<usize> = ty.grades.iter().map(|g| dim - g).collect();

        // Find a type with matching grades
        for candidate in &self.spec.types {
            if candidate.alias_of.is_some() {
                continue;
            }
            let mut candidate_grades = candidate.grades.clone();
            candidate_grades.sort();
            let mut sorted_complement = complement_grades.clone();
            sorted_complement.sort();
            if candidate_grades == sorted_complement {
                return Some(candidate.name.clone());
            }
        }

        // No matching type found
        None
    }

    /// Generates WeightDual trait impl.
    ///
    /// Returns None if the weight dual grades don't map to an existing type.
    fn generate_weight_dual_trait(&self, ty: &TypeSpec) -> Option<TokenStream> {
        // Check if there's a matching output type for the weight dual
        let output_type_name = self.find_weight_dual_output_type(ty)?;
        let output_type = self
            .spec
            .types
            .iter()
            .find(|t| t.name == output_type_name)?;

        let type_name = format_ident!("{}", ty.name);
        let out_name = format_ident!("{}", output_type_name);

        // Compute field expressions for weight dual
        let field_exprs: Vec<TokenStream> = output_type
            .fields
            .iter()
            .map(|out_field| {
                // Find the input field that weight-duals to this output blade
                let out_blade = out_field.blade_index;

                for in_field in &ty.fields {
                    let (sign, dual_blade) = self.table.weight_dual(in_field.blade_index);
                    if dual_blade == out_blade && sign != 0 {
                        let in_name = format_ident!("{}", in_field.name);
                        return if sign > 0 {
                            quote! { self.#in_name() }
                        } else {
                            quote! { -self.#in_name() }
                        };
                    }
                }
                // No input blade maps to this output blade
                quote! { T::zero() }
            })
            .collect();

        let constructor = quote! { #out_name::new_unchecked(#(#field_exprs),*) };

        Some(quote! {
            impl<T: Float> WeightDual for #type_name<T> {
                type Output = #out_name<T>;

                #[inline]
                fn weight_dual(&self) -> #out_name<T> {
                    #constructor
                }
            }
        })
    }

    /// Finds the output type for weight dual operations.
    ///
    /// The weight dual of a type with grades [g1, g2, ...] has grades [n-g1, n-g2, ...].
    /// Returns None if no type with matching grades exists.
    fn find_weight_dual_output_type(&self, ty: &TypeSpec) -> Option<String> {
        let dim = self.algebra.dim();
        let dual_grades: Vec<usize> = ty.grades.iter().map(|g| dim - g).collect();

        // Find a type with matching grades
        for candidate in &self.spec.types {
            if candidate.alias_of.is_some() {
                continue;
            }
            let mut candidate_grades = candidate.grades.clone();
            candidate_grades.sort();
            let mut sorted_dual = dual_grades.clone();
            sorted_dual.sort();
            if candidate_grades == sorted_dual {
                return Some(candidate.name.clone());
            }
        }

        // No matching type found
        None
    }

    // ========================================================================
    // Normed Trait Implementations
    // ========================================================================

    /// Computes the metric sign for a blade squared: `blade * blade`.
    ///
    /// For a blade `eᵢeⱼ...` with grade k, the square involves:
    /// 1. Reordering sign: `(-1)^(k(k-1)/2)` (to bring pairs together)
    /// 2. Product of individual basis vector metrics: `Π(eᵢ²)`
    ///
    /// Returns the combined sign (-1, 0, or +1).
    fn compute_blade_metric_sign(&self, blade_index: usize, grade: usize) -> i8 {
        // Reordering sign: (-1)^(k(k-1)/2)
        let reorder_sign: i8 = if (grade * grade.saturating_sub(1) / 2).is_multiple_of(2) {
            1
        } else {
            -1
        };

        // Product of individual basis vector metrics
        let mut metric_product: i8 = 1;
        for basis in &self.spec.signature.basis {
            if (blade_index >> basis.index) & 1 == 1 {
                // This basis vector is in the blade
                metric_product *= basis.metric;
            }
        }

        reorder_sign * metric_product
    }

    /// Generates all Normed trait implementations.
    fn generate_all_normed(&self) -> TokenStream {
        let impls: Vec<TokenStream> = self
            .spec
            .types
            .iter()
            .filter(|t| t.alias_of.is_none())
            .map(|ty| self.generate_normed_impl(ty))
            .collect();

        // For PGA algebras (with degenerate basis), also generate DegenerateNormed
        let degenerate_impls: Vec<TokenStream> = if self.spec.signature.r > 0 {
            self.spec
                .types
                .iter()
                .filter(|t| t.alias_of.is_none())
                .filter_map(|ty| self.generate_degenerate_normed_impl(ty))
                .collect()
        } else {
            Vec::new()
        };

        quote! {
            #(#impls)*
            #(#degenerate_impls)*
        }
    }

    /// Generates `impl Normed for Type<T>`.
    ///
    /// The norm is computed as `scalar_part(A * involute(A))` where
    /// `involute` uses the algebra's canonical involution.
    ///
    /// Each field contributes `sign * field²` where:
    /// - `sign = inv_sign * metric_sign`
    /// - `inv_sign` depends on the involution kind and grade
    /// - `metric_sign` is the blade's metric signature (product of basis metrics)
    fn generate_normed_impl(&self, ty: &TypeSpec) -> TokenStream {
        let name = format_ident!("{}", ty.name);
        let involution_kind = self.spec.norm.primary_involution;

        // Generate norm_squared: sum of (sign * field²) for all fields
        // Fields with zero metric (degenerate basis vectors) contribute nothing
        let norm_squared_terms: Vec<TokenStream> = ty
            .fields
            .iter()
            .filter_map(|f| {
                let fname = format_ident!("{}", f.name);

                // Compute involution sign for this field's grade
                let inv_sign: i8 = match involution_kind {
                    InvolutionKind::Reverse => {
                        // Reverse sign: (-1)^(k(k-1)/2)
                        let k = f.grade;
                        if (k * k.saturating_sub(1) / 2) % 2 == 0 {
                            1
                        } else {
                            -1
                        }
                    }
                    InvolutionKind::GradeInvolution => {
                        // Grade involution sign: (-1)^k
                        if f.grade % 2 == 0 { 1 } else { -1 }
                    }
                    InvolutionKind::CliffordConjugate => {
                        // Clifford conjugate sign: (-1)^(k(k+1)/2)
                        let k = f.grade;
                        if (k * (k + 1) / 2) % 2 == 0 { 1 } else { -1 }
                    }
                };

                // Compute metric sign: the sign of blade² under geometric product
                // For blade with index b, this is the product of (eᵢ)² for all basis vectors
                // times the sign from reordering (which depends on grade)
                let metric_sign = self.compute_blade_metric_sign(f.blade_index, f.grade);

                // Combined sign for this term
                let total_sign = inv_sign * metric_sign;

                // Skip terms with zero metric (degenerate basis vectors contribute nothing)
                if total_sign == 0 {
                    None
                } else if total_sign > 0 {
                    Some(quote! { self.#fname() * self.#fname() })
                } else {
                    Some(quote! { -self.#fname() * self.#fname() })
                }
            })
            .collect();

        // Generate scale: multiply each field by factor
        let scale_fields: Vec<TokenStream> = ty
            .fields
            .iter()
            .map(|f| {
                let fname = format_ident!("{}", f.name);
                quote! { self.#fname() * factor }
            })
            .collect();

        // Handle edge case where type has no fields
        let norm_squared_expr = if norm_squared_terms.is_empty() {
            quote! { T::zero() }
        } else {
            quote! { #(#norm_squared_terms)+* }
        };

        quote! {
            impl<T: Float> crate::norm::Normed for #name<T> {
                type Scalar = T;

                #[inline]
                fn norm_squared(&self) -> T {
                    #norm_squared_expr
                }

                fn try_normalize(&self) -> Option<Self> {
                    let n = self.norm();
                    if n < T::epsilon() {
                        None
                    } else {
                        Some(self.scale(T::one() / n))
                    }
                }

                #[inline]
                fn scale(&self, factor: T) -> Self {
                    Self::new_unchecked(#(#scale_fields),*)
                }
            }
        }
    }

    /// Generates `impl DegenerateNormed for Type<T>` for PGA types.
    ///
    /// Returns None if the type has no bulk or weight components.
    fn generate_degenerate_normed_impl(&self, ty: &TypeSpec) -> Option<TokenStream> {
        let name = format_ident!("{}", ty.name);

        // Find indices of degenerate basis vectors (those with metric == 0)
        let degenerate_indices: Vec<usize> = self
            .spec
            .signature
            .basis
            .iter()
            .filter(|b| b.metric == 0)
            .map(|b| b.index)
            .collect();

        // Partition fields into bulk (no degenerate basis) and weight (has degenerate basis)
        let (bulk_fields, weight_fields): (Vec<_>, Vec<_>) = ty.fields.iter().partition(|f| {
            // Check if this field's blade involves any degenerate basis vector
            !degenerate_indices.iter().any(|&deg_idx| {
                // Check if bit `deg_idx` is set in the blade_index
                (f.blade_index >> deg_idx) & 1 == 1
            })
        });

        // Don't generate if there are no fields in either category
        if bulk_fields.is_empty() && weight_fields.is_empty() {
            return None;
        }

        // Generate bulk_norm_squared: sum of squares of bulk fields
        let bulk_norm_terms: Vec<TokenStream> = bulk_fields
            .iter()
            .map(|f| {
                let fname = format_ident!("{}", f.name);
                quote! { self.#fname() * self.#fname() }
            })
            .collect();

        // Generate weight_norm_squared: sum of squares of weight fields
        let weight_norm_terms: Vec<TokenStream> = weight_fields
            .iter()
            .map(|f| {
                let fname = format_ident!("{}", f.name);
                quote! { self.#fname() * self.#fname() }
            })
            .collect();

        let bulk_norm_expr = if bulk_norm_terms.is_empty() {
            quote! { T::zero() }
        } else {
            quote! { #(#bulk_norm_terms)+* }
        };

        let weight_norm_expr = if weight_norm_terms.is_empty() {
            quote! { T::zero() }
        } else {
            quote! { #(#weight_norm_terms)+* }
        };

        // Generate scale fields for try_unitize
        let scale_fields: Vec<TokenStream> = ty
            .fields
            .iter()
            .map(|f| {
                let fname = format_ident!("{}", f.name);
                quote! { self.#fname() * inv_w }
            })
            .collect();

        Some(quote! {
            impl<T: Float> crate::norm::DegenerateNormed for #name<T> {
                #[inline]
                fn bulk_norm_squared(&self) -> T {
                    #bulk_norm_expr
                }

                #[inline]
                fn weight_norm_squared(&self) -> T {
                    #weight_norm_expr
                }

                fn try_unitize(&self) -> Option<Self> {
                    let w = self.weight_norm();
                    if w < T::epsilon() {
                        None
                    } else {
                        let inv_w = T::one() / w;
                        Some(Self::new_unchecked(#(#scale_fields),*))
                    }
                }
            }
        })
    }

    // ========================================================================
    // Approx Trait Implementations
    // ========================================================================

    /// Generates all approx trait implementations.
    fn generate_all_approx(&self) -> TokenStream {
        let impls: Vec<TokenStream> = self
            .spec
            .types
            .iter()
            .filter(|t| t.alias_of.is_none())
            .map(|ty| self.generate_approx_impls(ty))
            .collect();

        quote! { #(#impls)* }
    }

    /// Generates approx trait implementations for a type.
    fn generate_approx_impls(&self, ty: &TypeSpec) -> TokenStream {
        let name = format_ident!("{}", ty.name);

        let abs_diff_checks: Vec<TokenStream> = ty
            .fields
            .iter()
            .map(|f| {
                let fname = format_ident!("{}", f.name);
                quote! { self.#fname().abs_diff_eq(&other.#fname(), epsilon) }
            })
            .collect();

        let relative_checks: Vec<TokenStream> = ty
            .fields
            .iter()
            .map(|f| {
                let fname = format_ident!("{}", f.name);
                quote! { self.#fname().relative_eq(&other.#fname(), epsilon, max_relative) }
            })
            .collect();

        let ulps_checks: Vec<TokenStream> = ty
            .fields
            .iter()
            .map(|f| {
                let fname = format_ident!("{}", f.name);
                quote! { self.#fname().ulps_eq(&other.#fname(), epsilon, max_ulps) }
            })
            .collect();

        quote! {
            impl<T: Float + AbsDiffEq<Epsilon = T>> AbsDiffEq for #name<T> {
                type Epsilon = T;

                fn default_epsilon() -> Self::Epsilon {
                    T::default_epsilon()
                }

                fn abs_diff_eq(&self, other: &Self, epsilon: Self::Epsilon) -> bool {
                    #(#abs_diff_checks)&&*
                }
            }

            impl<T: Float + RelativeEq<Epsilon = T>> RelativeEq for #name<T> {
                fn default_max_relative() -> Self::Epsilon {
                    T::default_max_relative()
                }

                fn relative_eq(
                    &self,
                    other: &Self,
                    epsilon: Self::Epsilon,
                    max_relative: Self::Epsilon,
                ) -> bool {
                    #(#relative_checks)&&*
                }
            }

            impl<T: Float + UlpsEq<Epsilon = T>> UlpsEq for #name<T> {
                fn default_max_ulps() -> u32 {
                    T::default_max_ulps()
                }

                fn ulps_eq(&self, other: &Self, epsilon: Self::Epsilon, max_ulps: u32) -> bool {
                    #(#ulps_checks)&&*
                }
            }
        }
    }

    // ========================================================================
    // Arbitrary Implementations
    // ========================================================================

    /// Generates all Arbitrary implementations.
    fn generate_all_arbitrary(&self) -> TokenStream {
        let impls: Vec<TokenStream> = self
            .spec
            .types
            .iter()
            .filter(|t| t.alias_of.is_none())
            .map(|ty| self.generate_arbitrary_impl(ty))
            .collect();

        quote! {
            #[cfg(any(test, feature = "proptest-support"))]
            #[allow(clippy::missing_docs_in_private_items)]
            mod arbitrary_impls {
                use super::*;
                use proptest::prelude::*;
                use proptest::strategy::BoxedStrategy;
                use std::fmt::Debug;

                #(#impls)*
            }
        }
    }

    /// Generates Arbitrary implementation for a type.
    ///
    /// For unconstrained types, generates random values for all fields.
    /// For constrained types (with geometric constraints like Study condition),
    /// solves the constraint to compute dependent variables from independent ones.
    fn generate_arbitrary_impl(&self, ty: &TypeSpec) -> TokenStream {
        // Check if this type has a derived constraint
        if let Some(constraint_impl) = self.try_generate_constrained_arbitrary(ty) {
            return constraint_impl;
        }

        // No constraint - generate simple random values for all fields
        self.generate_unconstrained_arbitrary(ty)
    }

    /// Attempts to generate a constraint-solving Arbitrary implementation.
    ///
    /// Returns `Some(TokenStream)` if the type has a derived constraint that can be solved,
    /// `None` otherwise.
    fn try_generate_constrained_arbitrary(&self, ty: &TypeSpec) -> Option<TokenStream> {
        // Derive constraints from algebra structure
        // Uses the algebra's configured involution (reverse, grade involution, or Clifford conjugate)
        let deriver = ConstraintDeriver::new(self.algebra, self.spec.norm.primary_involution);
        let constraint = deriver.derive_geometric_constraint(ty, "x")?;

        // Only handle single-constraint cases for now
        if constraint.zero_expressions.len() != 1 {
            return None;
        }

        let expr = &constraint.zero_expressions[0];

        // Convert Symbolica expression to string for the solver
        let expr_str = format!("{} = 0", expr);

        // Find the highest-grade field to solve for (typically pseudoscalar like e0123)
        let solve_for_field = ty.fields.iter().max_by_key(|f| f.grade)?;

        // Try to solve the constraint for this variable
        let solver = ConstraintSolver::new();
        let symbol_name = format!("x_{}", solve_for_field.name);
        let solution = solver.solve(&expr_str, &symbol_name).ok()?;

        // For quadratic constraints (like Plücker), use filter instead of solving
        if solution.solution_type == SolutionType::Quadratic {
            return self.generate_filtered_arbitrary(ty);
        }

        // Generate constraint-solving Arbitrary
        Some(self.generate_solving_arbitrary(ty, &solve_for_field.name, &solution))
    }

    /// Generates Arbitrary that solves a linear constraint.
    fn generate_solving_arbitrary(
        &self,
        ty: &TypeSpec,
        solve_for: &str,
        solution: &crate::symbolic::SolveResult,
    ) -> TokenStream {
        let name = format_ident!("{}", ty.name);
        let num_fields = ty.fields.len();

        // Find indices of free and dependent fields
        let solve_for_idx = ty.fields.iter().position(|f| f.name == solve_for).unwrap();
        let free_indices: Vec<usize> = (0..num_fields).filter(|&i| i != solve_for_idx).collect();

        // Proptest only supports tuples up to 12 elements.
        // For types with more free variables, use Vec-based approach.
        if free_indices.len() > 12 {
            return self.generate_vec_based_solving_arbitrary(ty, solve_for, solution);
        }

        // Generate ranges for free variables
        let range_tuple: Vec<TokenStream> = free_indices
            .iter()
            .map(|_| quote! { -100.0f64..100.0 })
            .collect();

        // Generate variable names for prop_map
        let prop_map_args: Vec<TokenStream> = free_indices
            .iter()
            .enumerate()
            .map(|(i, _)| {
                let var = format_ident!("x{}", i);
                quote! { #var }
            })
            .collect();

        // Build the solution expression
        let numerator_expr = self.convert_solution_to_tokens(&solution.numerator, ty);
        let solution_expr = if let Some(ref divisor) = solution.divisor {
            let divisor_expr = self.convert_solution_to_tokens(divisor, ty);
            quote! { (#numerator_expr) / (#divisor_expr) }
        } else {
            numerator_expr
        };

        // Build field initialization expressions
        let mut field_var_map: Vec<Option<usize>> = vec![None; num_fields];
        for (var_idx, &field_idx) in free_indices.iter().enumerate() {
            field_var_map[field_idx] = Some(var_idx);
        }

        let field_inits: Vec<TokenStream> = ty
            .fields
            .iter()
            .enumerate()
            .map(|(i, _f)| {
                if i == solve_for_idx {
                    quote! { T::from_f64(#solution_expr) }
                } else {
                    let var_idx = field_var_map[i].unwrap();
                    let var = format_ident!("x{}", var_idx);
                    quote! { T::from_f64(#var) }
                }
            })
            .collect();

        // Generate filter for divisor non-zero condition
        let filter_expr = if let Some(ref divisor) = solution.divisor {
            let divisor_var = self.find_divisor_variable(divisor, ty);
            if let Some(var_idx) = divisor_var {
                let var = format_ident!("x{}", var_idx);
                // Generate filter args with underscores for unused variables
                let filter_args: Vec<TokenStream> = free_indices
                    .iter()
                    .enumerate()
                    .map(|(i, _)| {
                        if i == var_idx {
                            let v = format_ident!("x{}", i);
                            quote! { #v }
                        } else {
                            let v = format_ident!("_x{}", i);
                            quote! { #v }
                        }
                    })
                    .collect();
                Some(quote! {
                    .prop_filter("non-zero divisor", |(#(#filter_args),*)| (#var).abs() > 0.1)
                })
            } else {
                None
            }
        } else {
            None
        };

        let filter_chain = filter_expr.unwrap_or_else(|| quote! {});

        if free_indices.len() == 1 {
            let var = format_ident!("x0");
            quote! {
                impl<T: Float + Debug + 'static> Arbitrary for #name<T> {
                    type Parameters = ();
                    type Strategy = BoxedStrategy<Self>;

                    fn arbitrary_with(_: Self::Parameters) -> Self::Strategy {
                        (-100.0f64..100.0)
                            #filter_chain
                            .prop_map(|#var| {
                                #name::new_unchecked(#(#field_inits),*)
                            })
                            .boxed()
                    }
                }
            }
        } else {
            quote! {
                impl<T: Float + Debug + 'static> Arbitrary for #name<T> {
                    type Parameters = ();
                    type Strategy = BoxedStrategy<Self>;

                    fn arbitrary_with(_: Self::Parameters) -> Self::Strategy {
                        (#(#range_tuple),*)
                            #filter_chain
                            .prop_map(|(#(#prop_map_args),*)| {
                                #name::new_unchecked(#(#field_inits),*)
                            })
                            .boxed()
                    }
                }
            }
        }
    }

    /// Generates Vec-based Arbitrary for constrained types with more than 12 free variables.
    fn generate_vec_based_solving_arbitrary(
        &self,
        ty: &TypeSpec,
        solve_for: &str,
        solution: &crate::symbolic::SolveResult,
    ) -> TokenStream {
        let name = format_ident!("{}", ty.name);
        let num_fields = ty.fields.len();

        // Find indices of free and dependent fields
        let solve_for_idx = ty.fields.iter().position(|f| f.name == solve_for).unwrap();
        let free_indices: Vec<usize> = (0..num_fields).filter(|&i| i != solve_for_idx).collect();
        let num_free = free_indices.len();

        // Build the solution expression using v[i] instead of x{i}
        let numerator_expr = self.convert_solution_to_vec_tokens(&solution.numerator, ty);
        let solution_expr = if let Some(ref divisor) = solution.divisor {
            let divisor_expr = self.convert_solution_to_vec_tokens(divisor, ty);
            quote! { (#numerator_expr) / (#divisor_expr) }
        } else {
            numerator_expr
        };

        // Build field initialization expressions
        let mut field_var_map: Vec<Option<usize>> = vec![None; num_fields];
        for (var_idx, &field_idx) in free_indices.iter().enumerate() {
            field_var_map[field_idx] = Some(var_idx);
        }

        let field_inits: Vec<TokenStream> = ty
            .fields
            .iter()
            .enumerate()
            .map(|(i, _f)| {
                if i == solve_for_idx {
                    quote! { T::from_f64(#solution_expr) }
                } else {
                    let var_idx = field_var_map[i].unwrap();
                    quote! { T::from_f64(v[#var_idx]) }
                }
            })
            .collect();

        // Generate filter for divisor non-zero condition
        let filter_expr = solution.divisor.as_ref().and_then(|divisor| {
            self.find_divisor_variable(divisor, ty).map(|var_idx| {
                quote! {
                    .prop_filter("non-zero divisor", |v| v[#var_idx].abs() > 0.1)
                }
            })
        });

        let filter_chain = filter_expr.unwrap_or_else(|| quote! {});

        quote! {
            impl<T: Float + Debug + 'static> Arbitrary for #name<T> {
                type Parameters = ();
                type Strategy = BoxedStrategy<Self>;

                fn arbitrary_with(_: Self::Parameters) -> Self::Strategy {
                    proptest::collection::vec(-100.0f64..100.0, #num_free)
                        #filter_chain
                        .prop_map(|v| {
                            #name::new_unchecked(#(#field_inits),*)
                        })
                        .boxed()
                }
            }
        }
    }

    /// Converts a solution expression string to TokenStream.
    ///
    /// The solution from `ConstraintSolver` uses variable names like "x_s", "x_e23", etc.
    /// (field names with the `x_` prefix from ConstraintDeriver).
    /// We need to convert these to the corresponding `x{i}` variables.
    fn convert_solution_to_tokens(&self, expr: &str, ty: &TypeSpec) -> TokenStream {
        let mut result = expr.to_string();

        // Find field indices for free variables (all except the solved-for one)
        let solve_for_field = ty.fields.iter().max_by_key(|f| f.grade).unwrap();

        let free_fields: Vec<_> = ty
            .fields
            .iter()
            .filter(|f| f.name != solve_for_field.name)
            .collect();

        // Replace "x_fieldname" with "x{i}" variables
        // Do longer names first to avoid partial replacements
        let mut sorted_fields: Vec<_> = free_fields.iter().enumerate().collect();
        sorted_fields.sort_by(|a, b| b.1.name.len().cmp(&a.1.name.len()));

        for (i, field) in sorted_fields {
            let field_pattern = format!("x_{}", field.name);
            let var_name = format!("x{}", i);
            result = result.replace(&field_pattern, &var_name);
        }

        // Parse and convert to TokenStream
        result.parse().unwrap_or_else(|_| quote! { T::zero() })
    }

    /// Converts a solution expression string to TokenStream using Vec indexing.
    ///
    /// Similar to `convert_solution_to_tokens` but uses `v[i]` instead of `x{i}`.
    /// This is used for types with more than 12 fields where we can't use tuples.
    fn convert_solution_to_vec_tokens(&self, expr: &str, ty: &TypeSpec) -> TokenStream {
        let mut result = expr.to_string();

        // Find field indices for free variables (all except the solved-for one)
        let solve_for_field = ty.fields.iter().max_by_key(|f| f.grade).unwrap();

        let free_fields: Vec<_> = ty
            .fields
            .iter()
            .filter(|f| f.name != solve_for_field.name)
            .collect();

        // Replace "x_fieldname" with "v[i]" variables
        // Do longer names first to avoid partial replacements
        let mut sorted_fields: Vec<_> = free_fields.iter().enumerate().collect();
        sorted_fields.sort_by(|a, b| b.1.name.len().cmp(&a.1.name.len()));

        for (i, field) in sorted_fields {
            let field_pattern = format!("x_{}", field.name);
            let var_name = format!("v[{}]", i);
            result = result.replace(&field_pattern, &var_name);
        }

        // Parse and convert to TokenStream
        result.parse().unwrap_or_else(|_| quote! { T::zero() })
    }

    /// Finds which free variable corresponds to the divisor.
    ///
    /// The divisor uses the `x_fieldname` format from ConstraintDeriver.
    fn find_divisor_variable(&self, divisor: &str, ty: &TypeSpec) -> Option<usize> {
        let solve_for_field = ty.fields.iter().max_by_key(|f| f.grade)?;

        let free_fields: Vec<_> = ty
            .fields
            .iter()
            .filter(|f| f.name != solve_for_field.name)
            .collect();

        for (i, field) in free_fields.iter().enumerate() {
            // Check for "x_fieldname" pattern
            let pattern = format!("x_{}", field.name);
            if divisor.contains(&pattern) {
                return Some(i);
            }
        }
        None
    }

    /// Generates Arbitrary with a filter for quadratic constraints.
    ///
    /// For constraints like the Plücker condition that are quadratic and can't
    /// be easily solved, we generate random values and filter.
    fn generate_filtered_arbitrary(&self, ty: &TypeSpec) -> Option<TokenStream> {
        let name = format_ident!("{}", ty.name);
        let num_fields = ty.fields.len();

        // Generate tuple of ranges
        let range_tuple: Vec<TokenStream> =
            (0..num_fields).map(|_| quote! { -10.0f64..10.0 }).collect();

        let prop_map_args: Vec<TokenStream> = (0..num_fields)
            .map(|i| {
                let var = format_ident!("x{}", i);
                quote! { #var }
            })
            .collect();

        let field_inits: Vec<TokenStream> = (0..num_fields)
            .map(|i| {
                let var = format_ident!("x{}", i);
                quote! { T::from_f64(#var) }
            })
            .collect();

        Some(quote! {
            impl<T: Float + Debug + 'static> Arbitrary for #name<T> {
                type Parameters = ();
                type Strategy = BoxedStrategy<Self>;

                fn arbitrary_with(_: Self::Parameters) -> Self::Strategy {
                    (#(#range_tuple),*)
                        .prop_map(|(#(#prop_map_args),*)| {
                            #name::new_unchecked(#(#field_inits),*)
                        })
                        .boxed()
                }
            }
        })
    }

    /// Generates simple Arbitrary for unconstrained types.
    fn generate_unconstrained_arbitrary(&self, ty: &TypeSpec) -> TokenStream {
        let name = format_ident!("{}", ty.name);
        let num_fields = ty.fields.len();

        // Proptest only supports tuples up to 12 elements.
        // For types with more fields, use Vec-based approach.
        if num_fields > 12 {
            return self.generate_vec_based_arbitrary(ty);
        }

        // Generate tuple of ranges
        let range_tuple: Vec<TokenStream> = (0..num_fields)
            .map(|_| quote! { -100.0f64..100.0 })
            .collect();

        let field_inits: Vec<TokenStream> = (0..num_fields)
            .map(|i| {
                let var = format_ident!("x{}", i);
                quote! { T::from_f64(#var) }
            })
            .collect();

        if num_fields == 1 {
            quote! {
                impl<T: Float + Debug + 'static> Arbitrary for #name<T> {
                    type Parameters = ();
                    type Strategy = BoxedStrategy<Self>;

                    fn arbitrary_with(_: Self::Parameters) -> Self::Strategy {
                        (-100.0f64..100.0)
                            .prop_map(|x0| {
                                #name::new_unchecked(#(#field_inits),*)
                            })
                            .boxed()
                    }
                }
            }
        } else {
            let prop_map_args: Vec<TokenStream> = (0..num_fields)
                .map(|i| {
                    let var = format_ident!("x{}", i);
                    quote! { #var }
                })
                .collect();

            quote! {
                impl<T: Float + Debug + 'static> Arbitrary for #name<T> {
                    type Parameters = ();
                    type Strategy = BoxedStrategy<Self>;

                    fn arbitrary_with(_: Self::Parameters) -> Self::Strategy {
                        (#(#range_tuple),*)
                            .prop_map(|(#(#prop_map_args),*)| {
                                #name::new_unchecked(#(#field_inits),*)
                            })
                            .boxed()
                    }
                }
            }
        }
    }

    /// Generates Arbitrary using Vec for types with more than 12 fields.
    ///
    /// Proptest only implements Strategy for tuples up to 12 elements,
    /// so we use `prop::collection::vec` for larger types.
    fn generate_vec_based_arbitrary(&self, ty: &TypeSpec) -> TokenStream {
        let name = format_ident!("{}", ty.name);
        let num_fields = ty.fields.len();

        let field_inits: Vec<TokenStream> = (0..num_fields)
            .map(|i| {
                quote! { T::from_f64(v[#i]) }
            })
            .collect();

        quote! {
            impl<T: Float + Debug + 'static> Arbitrary for #name<T> {
                type Parameters = ();
                type Strategy = BoxedStrategy<Self>;

                fn arbitrary_with(_: Self::Parameters) -> Self::Strategy {
                    proptest::collection::vec(-100.0f64..100.0, #num_fields)
                        .prop_map(|v| {
                            #name::new_unchecked(#(#field_inits),*)
                        })
                        .boxed()
                }
            }
        }
    }

    // ========================================================================
    // Verification Tests
    // ========================================================================

    /// Generates verification tests as a pre-formatted string.
    ///
    /// This returns a raw string instead of TokenStream because rustfmt
    /// cannot format code inside proptest! macros. By generating the tests
    /// as pre-formatted strings, we preserve the formatting.
    fn generate_verification_tests_raw(&self) -> String {
        let signature_name = self.generate_signature_name();
        let add_sub_tests = self.generate_add_sub_verification_tests_raw();
        let exterior_tests = self.generate_exterior_verification_tests_raw();
        let bulk_contraction_tests = self.generate_bulk_contraction_verification_tests_raw();
        let weight_contraction_tests = self.generate_weight_contraction_verification_tests_raw();
        let bulk_expansion_tests = self.generate_bulk_expansion_verification_tests_raw();
        let weight_expansion_tests = self.generate_weight_expansion_verification_tests_raw();
        let de_morgan_tests = self.generate_de_morgan_verification_tests_raw();
        // Note: Project/Antiproject are only idempotent for unitized blades.
        // The formula `b ∨ (a ∧ b☆)` scales with |b|² each application.
        // Tests use Unit<T> for Euclidean or Unitized<T> for projective algebras.
        let project_idempotency_tests = self.generate_project_idempotency_tests_raw();
        let antiproject_idempotency_tests = self.generate_antiproject_idempotency_tests_raw();
        let wrapper_equivalence_tests = self.generate_wrapper_equivalence_tests_raw();
        let is_degenerate = self.spec.signature.r > 0;

        // For PGA, we need both Unitized (for blades), Bulk (for versors), and Unit (for tests)
        // All algebras need Unit for wrapper equivalence tests
        let wrapper_imports = if is_degenerate {
            "Unit, Unitized, Bulk"
        } else {
            "Unit"
        };

        format!(
            r#"
// ============================================================
// Verification Tests (compare against Multivector)
// ============================================================

#[cfg(test)]
#[allow(clippy::missing_docs_in_private_items)]
mod verification_tests {{
    use super::*;
    use crate::algebra::Multivector;
    use crate::signature::{sig};
    #[allow(unused_imports)]
    use crate::wrappers::{{{wrapper_imports}}};
    #[allow(unused_imports)]
    use crate::norm::{{Normed, DegenerateNormed}};
    use approx::relative_eq;
    use proptest::prelude::*;

    /// Relative epsilon for floating-point comparisons in verification tests.
    /// Using relative comparison handles varying magnitudes better than absolute.
    const REL_EPSILON: f64 = 1e-10;
{add_sub}{exterior}{bulk_contraction}{weight_contraction}{bulk_expansion}{weight_expansion}{de_morgan}{project_idempotency}{antiproject_idempotency}{wrapper_equivalence}}}
"#,
            sig = signature_name,
            wrapper_imports = wrapper_imports,
            add_sub = add_sub_tests,
            exterior = exterior_tests,
            bulk_contraction = bulk_contraction_tests,
            weight_contraction = weight_contraction_tests,
            bulk_expansion = bulk_expansion_tests,
            weight_expansion = weight_expansion_tests,
            de_morgan = de_morgan_tests,
            project_idempotency = project_idempotency_tests,
            antiproject_idempotency = antiproject_idempotency_tests,
            wrapper_equivalence = wrapper_equivalence_tests,
        )
    }

    /// Generates the signature type name for this algebra.
    ///
    /// Uses the signature tuple (p, q, r) to determine the signature type,
    /// not the algebra name. This ensures generic handling of all algebras.
    fn generate_signature_name(&self) -> proc_macro2::Ident {
        let sig = &self.spec.signature;
        let (p, q, r) = (sig.p, sig.q, sig.r);

        // Derive signature type from (p, q, r)
        let sig_name = match (p, q, r) {
            // Euclidean: Cl(n, 0, 0)
            (2, 0, 0) => "Euclidean2",
            (3, 0, 0) => "Euclidean3",

            // Projective (PGA): Cl(n, 0, 1)
            (2, 0, 1) => "Projective2",
            (3, 0, 1) => "Projective3",

            // Conformal (CGA): Cl(n+1, 1, 0) - note: uses p+q=4/5 convention
            // CGA 2D: Cl(3, 1, 0) - 2D + 2 extra dimensions
            // CGA 3D: Cl(4, 1, 0) - 3D + 2 extra dimensions
            (4, 1, 0) => "Conformal3",

            // Minkowski spacetime: Cl(1, 3, 0)
            (1, 3, 0) => "Minkowski4",

            // Generic: use Cl{p}_{q}_{r} format
            _ => {
                return format_ident!("Cl{}_{}_{}", p, q, r);
            }
        };
        format_ident!("{}", sig_name)
    }

    /// Generates add/sub verification tests for each type as a formatted string.
    fn generate_add_sub_verification_tests_raw(&self) -> String {
        let signature_name = self.generate_signature_name();
        self.spec
            .types
            .iter()
            .filter(|t| t.alias_of.is_none())
            .map(|ty| {
                let name = &ty.name;
                let name_lower = ty.name.to_lowercase();

                format!(
                    r#"
    proptest! {{
        #[test]
        fn {name_lower}_add_matches_multivector(a in any::<{name}<f64>>(), b in any::<{name}<f64>>()) {{
            let mv_a: Multivector<f64, {sig}> = a.into();
            let mv_b: Multivector<f64, {sig}> = b.into();

            let specialized_result = a + b;
            let generic_result = mv_a + mv_b;

            let specialized_mv: Multivector<f64, {sig}> = specialized_result.into();
            prop_assert!(
                relative_eq!(specialized_mv, generic_result, epsilon = REL_EPSILON, max_relative = REL_EPSILON),
                "Add mismatch: specialized={{:?}}, generic={{:?}}",
                specialized_mv, generic_result
            );
        }}

        #[test]
        fn {name_lower}_sub_matches_multivector(a in any::<{name}<f64>>(), b in any::<{name}<f64>>()) {{
            let mv_a: Multivector<f64, {sig}> = a.into();
            let mv_b: Multivector<f64, {sig}> = b.into();

            let specialized_result = a - b;
            let generic_result = mv_a - mv_b;

            let specialized_mv: Multivector<f64, {sig}> = specialized_result.into();
            prop_assert!(
                relative_eq!(specialized_mv, generic_result, epsilon = REL_EPSILON, max_relative = REL_EPSILON),
                "Sub mismatch: specialized={{:?}}, generic={{:?}}",
                specialized_mv, generic_result
            );
        }}

        #[test]
        fn {name_lower}_neg_matches_multivector(a in any::<{name}<f64>>()) {{
            let mv_a: Multivector<f64, {sig}> = a.into();

            let specialized_result = -a;
            let generic_result = -mv_a;

            let specialized_mv: Multivector<f64, {sig}> = specialized_result.into();
            prop_assert!(
                relative_eq!(specialized_mv, generic_result, epsilon = REL_EPSILON, max_relative = REL_EPSILON),
                "Neg mismatch: specialized={{:?}}, generic={{:?}}",
                specialized_mv, generic_result
            );
        }}
    }}
"#,
                    name_lower = name_lower,
                    name = name,
                    sig = signature_name,
                )
            })
            .collect()
    }

    /// Generates exterior product verification tests as a formatted string.
    fn generate_exterior_verification_tests_raw(&self) -> String {
        // Only generate tests for products explicitly listed in the TOML
        if self.spec.products.wedge.is_empty() {
            return String::new();
        }

        let signature_name = self.generate_signature_name();
        self.spec
            .products
            .wedge
            .iter()
            // Only generate tests for single-grade types (where Wedge is implemented)
            .filter(|entry| {
                self.find_type(&entry.lhs)
                    .is_some_and(|t| self.is_single_grade_blade(t))
                    && self
                        .find_type(&entry.rhs)
                        .is_some_and(|t| self.is_single_grade_blade(t))
            })
            .map(|entry| {
                let lhs_lower = entry.lhs.to_lowercase();
                let rhs_lower = entry.rhs.to_lowercase();
                let out_lower = entry.output.to_lowercase();

                format!(
                    r#"
    proptest! {{
        #[test]
        fn wedge_{lhs_lower}_{rhs_lower}_{out_lower}_matches_multivector(a in any::<{lhs}<f64>>(), b in any::<{rhs}<f64>>()) {{
            use crate::ops::Wedge;
            let mv_a: Multivector<f64, {sig}> = a.into();
            let mv_b: Multivector<f64, {sig}> = b.into();

            let specialized_result: {out}<f64> = a.wedge(&b);
            let generic_result = mv_a.exterior(&mv_b);

            let specialized_mv: Multivector<f64, {sig}> = specialized_result.into();
            prop_assert!(
                relative_eq!(specialized_mv, generic_result, epsilon = REL_EPSILON, max_relative = REL_EPSILON),
                "Wedge product mismatch: specialized={{:?}}, generic={{:?}}",
                specialized_mv, generic_result
            );
        }}
    }}
"#,
                    lhs_lower = lhs_lower,
                    rhs_lower = rhs_lower,
                    out_lower = out_lower,
                    lhs = entry.lhs,
                    rhs = entry.rhs,
                    out = entry.output,
                    sig = signature_name,
                )
            })
            .collect()
    }

    /// Generates bulk contraction verification tests as a formatted string.
    fn generate_bulk_contraction_verification_tests_raw(&self) -> String {
        if self.spec.products.bulk_contraction.is_empty() {
            return String::new();
        }

        let signature_name = self.generate_signature_name();
        self.spec
            .products
            .bulk_contraction
            .iter()
            // Only generate tests for single-grade types (where BulkContract is implemented)
            .filter(|entry| {
                self.find_type(&entry.lhs)
                    .is_some_and(|t| self.is_single_grade_blade(t))
                    && self
                        .find_type(&entry.rhs)
                        .is_some_and(|t| self.is_single_grade_blade(t))
            })
            .map(|entry| {
                let lhs_lower = entry.lhs.to_lowercase();
                let rhs_lower = entry.rhs.to_lowercase();
                let out_lower = entry.output.to_lowercase();

                format!(
                    r#"
    proptest! {{
        #[test]
        fn bulk_contraction_{lhs_lower}_{rhs_lower}_{out_lower}_matches_multivector(a in any::<{lhs}<f64>>(), b in any::<{rhs}<f64>>()) {{
            use crate::ops::BulkContract;
            let mv_a: Multivector<f64, {sig}> = a.into();
            let mv_b: Multivector<f64, {sig}> = b.into();

            let specialized_result: {out}<f64> = a.bulk_contract(&b);
            let generic_result = mv_a.bulk_contraction(&mv_b);

            let specialized_mv: Multivector<f64, {sig}> = specialized_result.into();
            prop_assert!(
                relative_eq!(specialized_mv, generic_result, epsilon = REL_EPSILON, max_relative = REL_EPSILON),
                "Bulk contraction mismatch: specialized={{:?}}, generic={{:?}}",
                specialized_mv, generic_result
            );
        }}
    }}
"#,
                    lhs_lower = lhs_lower,
                    rhs_lower = rhs_lower,
                    out_lower = out_lower,
                    lhs = entry.lhs,
                    rhs = entry.rhs,
                    out = entry.output,
                    sig = signature_name,
                )
            })
            .collect()
    }

    /// Generates weight contraction verification tests as a formatted string.
    fn generate_weight_contraction_verification_tests_raw(&self) -> String {
        if self.spec.products.weight_contraction.is_empty() {
            return String::new();
        }

        let signature_name = self.generate_signature_name();
        self.spec
            .products
            .weight_contraction
            .iter()
            // Only generate tests for single-grade types (where WeightContract is implemented)
            .filter(|entry| {
                self.find_type(&entry.lhs)
                    .is_some_and(|t| self.is_single_grade_blade(t))
                    && self
                        .find_type(&entry.rhs)
                        .is_some_and(|t| self.is_single_grade_blade(t))
            })
            .map(|entry| {
                let lhs_lower = entry.lhs.to_lowercase();
                let rhs_lower = entry.rhs.to_lowercase();
                let out_lower = entry.output.to_lowercase();

                format!(
                    r#"
    proptest! {{
        #[test]
        fn weight_contraction_{lhs_lower}_{rhs_lower}_{out_lower}_matches_multivector(a in any::<{lhs}<f64>>(), b in any::<{rhs}<f64>>()) {{
            use crate::ops::WeightContract;
            let mv_a: Multivector<f64, {sig}> = a.into();
            let mv_b: Multivector<f64, {sig}> = b.into();

            let specialized_result: {out}<f64> = a.weight_contract(&b);
            let generic_result = mv_a.weight_contraction(&mv_b);

            let specialized_mv: Multivector<f64, {sig}> = specialized_result.into();
            prop_assert!(
                relative_eq!(specialized_mv, generic_result, epsilon = REL_EPSILON, max_relative = REL_EPSILON),
                "Weight contraction mismatch: specialized={{:?}}, generic={{:?}}",
                specialized_mv, generic_result
            );
        }}
    }}
"#,
                    lhs_lower = lhs_lower,
                    rhs_lower = rhs_lower,
                    out_lower = out_lower,
                    lhs = entry.lhs,
                    rhs = entry.rhs,
                    out = entry.output,
                    sig = signature_name,
                )
            })
            .collect()
    }

    /// Generates bulk expansion verification tests as a formatted string.
    fn generate_bulk_expansion_verification_tests_raw(&self) -> String {
        if self.spec.products.bulk_expansion.is_empty() {
            return String::new();
        }

        let signature_name = self.generate_signature_name();
        self.spec
            .products
            .bulk_expansion
            .iter()
            // Only generate tests for single-grade types (where BulkExpand is implemented)
            .filter(|entry| {
                self.find_type(&entry.lhs)
                    .is_some_and(|t| self.is_single_grade_blade(t))
                    && self
                        .find_type(&entry.rhs)
                        .is_some_and(|t| self.is_single_grade_blade(t))
            })
            .map(|entry| {
                let lhs_lower = entry.lhs.to_lowercase();
                let rhs_lower = entry.rhs.to_lowercase();
                let out_lower = entry.output.to_lowercase();

                format!(
                    r#"
    proptest! {{
        #[test]
        fn bulk_expansion_{lhs_lower}_{rhs_lower}_{out_lower}_matches_multivector(a in any::<{lhs}<f64>>(), b in any::<{rhs}<f64>>()) {{
            use crate::ops::BulkExpand;
            let mv_a: Multivector<f64, {sig}> = a.into();
            let mv_b: Multivector<f64, {sig}> = b.into();

            let specialized_result: {out}<f64> = a.bulk_expand(&b);
            let generic_result = mv_a.bulk_expansion(&mv_b);

            let specialized_mv: Multivector<f64, {sig}> = specialized_result.into();
            prop_assert!(
                relative_eq!(specialized_mv, generic_result, epsilon = REL_EPSILON, max_relative = REL_EPSILON),
                "Bulk expansion mismatch: specialized={{:?}}, generic={{:?}}",
                specialized_mv, generic_result
            );
        }}
    }}
"#,
                    lhs_lower = lhs_lower,
                    rhs_lower = rhs_lower,
                    out_lower = out_lower,
                    lhs = entry.lhs,
                    rhs = entry.rhs,
                    out = entry.output,
                    sig = signature_name,
                )
            })
            .collect()
    }

    /// Generates weight expansion verification tests as a formatted string.
    fn generate_weight_expansion_verification_tests_raw(&self) -> String {
        if self.spec.products.weight_expansion.is_empty() {
            return String::new();
        }

        let signature_name = self.generate_signature_name();
        self.spec
            .products
            .weight_expansion
            .iter()
            // Only generate tests for single-grade types (where WeightExpand is implemented)
            .filter(|entry| {
                self.find_type(&entry.lhs)
                    .is_some_and(|t| self.is_single_grade_blade(t))
                    && self
                        .find_type(&entry.rhs)
                        .is_some_and(|t| self.is_single_grade_blade(t))
            })
            .map(|entry| {
                let lhs_lower = entry.lhs.to_lowercase();
                let rhs_lower = entry.rhs.to_lowercase();
                let out_lower = entry.output.to_lowercase();

                format!(
                    r#"
    proptest! {{
        #[test]
        fn weight_expansion_{lhs_lower}_{rhs_lower}_{out_lower}_matches_multivector(a in any::<{lhs}<f64>>(), b in any::<{rhs}<f64>>()) {{
            use crate::ops::WeightExpand;
            let mv_a: Multivector<f64, {sig}> = a.into();
            let mv_b: Multivector<f64, {sig}> = b.into();

            let specialized_result: {out}<f64> = a.weight_expand(&b);
            let generic_result = mv_a.weight_expansion(&mv_b);

            let specialized_mv: Multivector<f64, {sig}> = specialized_result.into();
            prop_assert!(
                relative_eq!(specialized_mv, generic_result, epsilon = REL_EPSILON, max_relative = REL_EPSILON),
                "Weight expansion mismatch: specialized={{:?}}, generic={{:?}}",
                specialized_mv, generic_result
            );
        }}
    }}
"#,
                    lhs_lower = lhs_lower,
                    rhs_lower = rhs_lower,
                    out_lower = out_lower,
                    lhs = entry.lhs,
                    rhs = entry.rhs,
                    out = entry.output,
                    sig = signature_name,
                )
            })
            .collect()
    }

    // Note: Project and Antiproject verification tests are not generated because
    // the specialized types don't implement these traits (see comment above).

    /// Generates de-Morgan's law verification tests as a formatted string.
    ///
    /// Tests the fundamental identities from RGA:
    /// - complement(a * b) = complement(a) ⋇ complement(b)
    /// - complement(a ⋇ b) = complement(a) * complement(b)
    ///
    /// These tests verify that the complement and antiproduct operations
    /// on Multivector satisfy the de-Morgan duality laws.
    ///
    /// Note: Only generates tests for single-grade types. Mixed-grade types
    /// can have algebra-dependent sign factors in the De Morgan laws that
    /// make the simple identity not hold exactly.
    fn generate_de_morgan_verification_tests_raw(&self) -> String {
        let signature_name = self.generate_signature_name();

        // Generate tests for single-grade types only (non-alias)
        // Mixed-grade types can have sign factors in De Morgan laws
        self.spec
            .types
            .iter()
            .filter(|t| t.alias_of.is_none() && t.grades.len() == 1)
            .map(|ty| {
                let name = &ty.name;
                let name_lower = ty.name.to_lowercase();

                format!(
                    r#"
    proptest! {{
        /// De Morgan: complement(a * b) = complement(a) ⋇ complement(b)
        #[test]
        fn de_morgan_geometric_{name_lower}(a in any::<{name}<f64>>(), b in any::<{name}<f64>>()) {{
            let mv_a: Multivector<f64, {sig}> = a.into();
            let mv_b: Multivector<f64, {sig}> = b.into();

            // LHS: complement(a * b)
            let lhs = (mv_a * mv_b).complement();

            // RHS: complement(a) ⋇ complement(b)
            let rhs = mv_a.complement().antiproduct(&mv_b.complement());

            prop_assert!(
                relative_eq!(lhs, rhs, epsilon = REL_EPSILON, max_relative = REL_EPSILON),
                "De Morgan (geometric) failed: complement(a*b)={{:?}}, complement(a)⋇complement(b)={{:?}}",
                lhs, rhs
            );
        }}

        /// De Morgan: complement(a ⋇ b) = complement(a) * complement(b)
        #[test]
        fn de_morgan_antiproduct_{name_lower}(a in any::<{name}<f64>>(), b in any::<{name}<f64>>()) {{
            let mv_a: Multivector<f64, {sig}> = a.into();
            let mv_b: Multivector<f64, {sig}> = b.into();

            // LHS: complement(a ⋇ b)
            let lhs = mv_a.antiproduct(&mv_b).complement();

            // RHS: complement(a) * complement(b)
            let rhs = mv_a.complement() * mv_b.complement();

            prop_assert!(
                relative_eq!(lhs, rhs, epsilon = REL_EPSILON, max_relative = REL_EPSILON),
                "De Morgan (antiproduct) failed: complement(a⋇b)={{:?}}, complement(a)*complement(b)={{:?}}",
                lhs, rhs
            );
        }}
    }}
"#,
                    name_lower = name_lower,
                    name = name,
                    sig = signature_name,
                )
            })
            .collect()
    }

    /// Generates project idempotency tests with normalized targets.
    ///
    /// Tests that projection is idempotent when the target is normalized:
    /// `project(project(a, unit_b), unit_b) == project(a, unit_b)`
    ///
    /// Uses `Unit<T>` for Euclidean algebras (Normed) and `Unitized<T>` for
    /// projective algebras (DegenerateNormed).
    ///
    /// Only generates tests where grade(a) > grade(b), as projection is
    /// geometrically meaningful for projecting higher-grade onto lower-grade.
    fn generate_project_idempotency_tests_raw(&self) -> String {
        let sig = self.generate_signature_name().to_string();
        let is_degenerate = self.spec.signature.r > 0;
        let wrapper = if is_degenerate { "Unitized" } else { "Unit" };

        // Get all single-grade types (non-alias) with their grades
        let single_grade_types: Vec<_> = self
            .spec
            .types
            .iter()
            .filter(|t| t.alias_of.is_none() && self.is_single_grade_blade(t))
            .collect();

        // Generate tests only for grade(a) > grade(b)
        let mut result = String::new();
        for ty_a in &single_grade_types {
            for ty_b in &single_grade_types {
                // Get the grade of each type from its first field
                let grade_a = ty_a
                    .fields
                    .first()
                    .map(|f| (f.blade_index as u32).count_ones() as usize)
                    .unwrap_or(0);
                let grade_b = ty_b
                    .fields
                    .first()
                    .map(|f| (f.blade_index as u32).count_ones() as usize)
                    .unwrap_or(0);

                // Only test when projecting higher grade onto lower grade
                // Skip grade-0 (projecting scalar or onto scalar is meaningless)
                if grade_a <= grade_b || grade_a == 0 || grade_b == 0 {
                    continue;
                }

                let a_name = &ty_a.name;
                let b_name = &ty_b.name;
                let a_lower = ty_a.name.to_lowercase();
                let b_lower = ty_b.name.to_lowercase();

                result.push_str(&format!(
                    r#"
    proptest! {{
        /// Project idempotency with normalized target: project(project(a, unit_b), unit_b) == project(a, unit_b)
        #[test]
        fn project_idempotent_{a_lower}_{b_lower}(a in any::<{a_name}<f64>>(), unit_b in any::<{wrapper}<{b_name}<f64>>>()) {{
            let mv_a: Multivector<f64, {sig}> = a.into();
            let mv_b: Multivector<f64, {sig}> = unit_b.into_inner().into();

            let first = mv_a.project(&mv_b);
            let second = first.project(&mv_b);

            prop_assert!(
                relative_eq!(first, second, epsilon = REL_EPSILON, max_relative = REL_EPSILON),
                "Project idempotency failed: first={{:?}}, second={{:?}}",
                first, second
            );
        }}
    }}
"#,
                    a_lower = a_lower,
                    b_lower = b_lower,
                    a_name = a_name,
                    b_name = b_name,
                    wrapper = wrapper,
                    sig = sig,
                ));
            }
        }
        result
    }

    /// Generates antiproject idempotency tests with normalized targets.
    ///
    /// Tests that antiprojection is idempotent when the target is normalized:
    /// `antiproject(antiproject(a, unit_b), unit_b) == antiproject(a, unit_b)`
    ///
    /// Uses `Unit<T>` for Euclidean algebras (Normed) and `Unitized<T>` for
    /// projective algebras (DegenerateNormed).
    ///
    /// Only generates tests where grade(a) < grade(b), as antiprojection is
    /// geometrically meaningful for projecting lower-grade onto higher-grade.
    fn generate_antiproject_idempotency_tests_raw(&self) -> String {
        let sig = self.generate_signature_name().to_string();
        let is_degenerate = self.spec.signature.r > 0;
        let wrapper = if is_degenerate { "Unitized" } else { "Unit" };

        // Get all single-grade types (non-alias) with their grades
        let single_grade_types: Vec<_> = self
            .spec
            .types
            .iter()
            .filter(|t| t.alias_of.is_none() && self.is_single_grade_blade(t))
            .collect();

        // Generate tests only for grade(a) < grade(b)
        let mut result = String::new();
        for ty_a in &single_grade_types {
            for ty_b in &single_grade_types {
                // Get the grade of each type from its first field
                let grade_a = ty_a
                    .fields
                    .first()
                    .map(|f| (f.blade_index as u32).count_ones() as usize)
                    .unwrap_or(0);
                let grade_b = ty_b
                    .fields
                    .first()
                    .map(|f| (f.blade_index as u32).count_ones() as usize)
                    .unwrap_or(0);

                // Only test when projecting lower grade onto higher grade
                // Skip grade-0 sources (projecting scalar is meaningless)
                if grade_a >= grade_b || grade_a == 0 {
                    continue;
                }

                let a_name = &ty_a.name;
                let b_name = &ty_b.name;
                let a_lower = ty_a.name.to_lowercase();
                let b_lower = ty_b.name.to_lowercase();

                result.push_str(&format!(
                    r#"
    proptest! {{
        /// Antiproject idempotency with normalized target: antiproject(antiproject(a, unit_b), unit_b) == antiproject(a, unit_b)
        #[test]
        fn antiproject_idempotent_{a_lower}_{b_lower}(a in any::<{a_name}<f64>>(), unit_b in any::<{wrapper}<{b_name}<f64>>>()) {{
            let mv_a: Multivector<f64, {sig}> = a.into();
            let mv_b: Multivector<f64, {sig}> = unit_b.into_inner().into();

            let first = mv_a.antiproject(&mv_b);
            let second = first.antiproject(&mv_b);

            prop_assert!(
                relative_eq!(first, second, epsilon = REL_EPSILON, max_relative = REL_EPSILON),
                "Antiproject idempotency failed: first={{:?}}, second={{:?}}",
                first, second
            );
        }}
    }}
"#,
                    a_lower = a_lower,
                    b_lower = b_lower,
                    a_name = a_name,
                    b_name = b_name,
                    wrapper = wrapper,
                    sig = sig,
                ));
            }
        }
        result
    }

    /// Generates wrapper equivalence tests for Normed trait.
    ///
    /// For types that have `Normed` implementations, this generates tests verifying
    /// that `Unit<T>.method()` equals `T.method()` for all Normed methods.
    /// The wrapper implementations are optimized (e.g., `Unit<T>.norm() == 1`),
    /// but must produce semantically equivalent results.
    fn generate_wrapper_equivalence_tests_raw(&self) -> String {
        use crate::spec::InvolutionKind;

        let is_degenerate = self.spec.signature.r > 0;
        // Unit<T> only works for algebras with positive-definite norm:
        // - q == 0 (no negative squares in metric)
        // - r == 0 (no degenerate dimensions)
        // - primary_involution == Reverse (not GradeInvolution which is indefinite)
        let has_positive_definite_norm = self.spec.signature.q == 0
            && self.spec.signature.r == 0
            && self.spec.norm.primary_involution == InvolutionKind::Reverse;

        // Only generate Unit<T> tests for algebras with positive-definite norm
        let unit_tests: String = if has_positive_definite_norm {
            self.spec
                .types
                .iter()
                .filter(|t| t.alias_of.is_none())
                // Filter to types that are normalizable (have norm methods)
                .filter(|t| {
                    // All types with grades implement Normed
                    !t.grades.is_empty()
                })
            .map(|ty| {
                let name = &ty.name;
                let name_lower = ty.name.to_lowercase();

                format!(
                    r#"
    proptest! {{
        /// Unit<{name}>.norm() should equal inner's norm (both are 1.0).
        #[test]
        fn unit_{name_lower}_norm_matches_inner(u in any::<Unit<{name}<f64>>>()) {{
            // Use explicit trait syntax to specify the type
            let inner_norm = <{name}<f64> as Normed>::norm(u.as_inner());
            let wrapper_norm = <Unit<{name}<f64>> as Normed>::norm(&u);

            prop_assert!(
                relative_eq!(inner_norm, 1.0, epsilon = REL_EPSILON, max_relative = REL_EPSILON),
                "Inner norm should be 1.0, got {{}}", inner_norm
            );
            prop_assert!(
                relative_eq!(wrapper_norm, 1.0, epsilon = REL_EPSILON, max_relative = REL_EPSILON),
                "Wrapper norm should be 1.0, got {{}}", wrapper_norm
            );
            prop_assert!(
                relative_eq!(inner_norm, wrapper_norm, epsilon = REL_EPSILON, max_relative = REL_EPSILON),
                "Norms should match: {{}} vs {{}}", inner_norm, wrapper_norm
            );
        }}

        /// Unit<{name}>.norm_squared() should equal inner's norm_squared (both are 1.0).
        #[test]
        fn unit_{name_lower}_norm_squared_matches_inner(u in any::<Unit<{name}<f64>>>()) {{
            // Use explicit trait syntax to specify the type
            let inner_ns = <{name}<f64> as Normed>::norm_squared(u.as_inner());
            let wrapper_ns = <Unit<{name}<f64>> as Normed>::norm_squared(&u);

            prop_assert!(
                relative_eq!(inner_ns, 1.0, epsilon = REL_EPSILON, max_relative = REL_EPSILON),
                "Inner norm_squared should be 1.0, got {{}}", inner_ns
            );
            prop_assert!(
                relative_eq!(wrapper_ns, 1.0, epsilon = REL_EPSILON, max_relative = REL_EPSILON),
                "Wrapper norm_squared should be 1.0, got {{}}", wrapper_ns
            );
        }}
    }}
"#,
                    name = name,
                    name_lower = name_lower,
                )
            })
            .collect()
        } else {
            String::new()
        };

        // For degenerate algebras, also generate Bulk<T> tests for DegenerateNormed
        let bulk_tests: String = if is_degenerate {
            self.spec
                .types
                .iter()
                .filter(|t| t.alias_of.is_none())
                // Filter to types that are versors (have bulk_norm)
                .filter(|t| t.versor.is_some())
                // Filter to types that have non-zero bulk norm (exclude purely degenerate types)
                .filter(|t| self.has_nonzero_bulk_norm(t))
                .map(|ty| {
                    let name = &ty.name;
                    let name_lower = ty.name.to_lowercase();

                    format!(
                        r#"
    proptest! {{
        /// Bulk<{name}>.bulk_norm() should equal 1.0 (by definition of Bulk wrapper).
        #[test]
        fn bulk_{name_lower}_bulk_norm_matches_inner(b in any::<Bulk<{name}<f64>>>()) {{
            // Use explicit trait syntax to specify the type
            let inner_bulk = <{name}<f64> as DegenerateNormed>::bulk_norm(b.as_inner());
            let wrapper_bulk = <Bulk<{name}<f64>> as DegenerateNormed>::bulk_norm(&b);

            prop_assert!(
                relative_eq!(inner_bulk, 1.0, epsilon = REL_EPSILON, max_relative = REL_EPSILON),
                "Inner bulk_norm should be 1.0, got {{}}", inner_bulk
            );
            prop_assert!(
                relative_eq!(wrapper_bulk, 1.0, epsilon = REL_EPSILON, max_relative = REL_EPSILON),
                "Wrapper bulk_norm should be 1.0, got {{}}", wrapper_bulk
            );
        }}

        /// Bulk<{name}>.weight_norm() should match inner's weight_norm (delegation).
        #[test]
        fn bulk_{name_lower}_weight_norm_delegates(b in any::<Bulk<{name}<f64>>>()) {{
            // Use explicit trait syntax to specify the type
            let inner_weight = <{name}<f64> as DegenerateNormed>::weight_norm(b.as_inner());
            let wrapper_weight = <Bulk<{name}<f64>> as DegenerateNormed>::weight_norm(&b);

            prop_assert!(
                relative_eq!(inner_weight, wrapper_weight, epsilon = REL_EPSILON, max_relative = REL_EPSILON),
                "Weight norms should match: {{}} vs {{}}", inner_weight, wrapper_weight
            );
        }}
    }}
"#,
                        name = name,
                        name_lower = name_lower,
                    )
                })
                .collect()
        } else {
            String::new()
        };

        format!("{}{}", unit_tests, bulk_tests)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::spec::parse_spec;

    #[test]
    fn symbolica_generates_add_impl() {
        let spec = parse_spec(include_str!("../../algebras/euclidean3.toml")).unwrap();
        let algebra = Algebra::euclidean(3);
        let table = ProductTable::new(&algebra);
        let generator = TraitsGenerator::new(&spec, &algebra, table);

        let (tokens, _tests) = generator.generate_traits_file();
        let code = tokens.to_string();

        assert!(code.contains("impl < T : Float > Add for Vector"));
    }

    #[test]
    fn symbolica_generates_sub_impl() {
        let spec = parse_spec(include_str!("../../algebras/euclidean3.toml")).unwrap();
        let algebra = Algebra::euclidean(3);
        let table = ProductTable::new(&algebra);
        let generator = TraitsGenerator::new(&spec, &algebra, table);

        let (tokens, _tests) = generator.generate_traits_file();
        let code = tokens.to_string();

        assert!(code.contains("impl < T : Float > Sub for Vector"));
    }

    #[test]
    fn symbolica_generates_neg_impl() {
        let spec = parse_spec(include_str!("../../algebras/euclidean3.toml")).unwrap();
        let algebra = Algebra::euclidean(3);
        let table = ProductTable::new(&algebra);
        let generator = TraitsGenerator::new(&spec, &algebra, table);

        let (tokens, _tests) = generator.generate_traits_file();
        let code = tokens.to_string();

        assert!(code.contains("impl < T : Float > Neg for Vector"));
    }

    #[test]
    fn symbolica_generates_scalar_mul() {
        let spec = parse_spec(include_str!("../../algebras/euclidean3.toml")).unwrap();
        let algebra = Algebra::euclidean(3);
        let table = ProductTable::new(&algebra);
        let generator = TraitsGenerator::new(&spec, &algebra, table);

        let (tokens, _tests) = generator.generate_traits_file();
        let code = tokens.to_string();

        // Check that scalar mul is generated for Vector
        assert!(code.contains("Mul"));
        assert!(code.contains("for Vector"));
        assert!(code.contains("scale"));
    }

    #[test]
    fn symbolica_generates_geometric_mul() {
        let spec = parse_spec(include_str!("../../algebras/euclidean3.toml")).unwrap();
        let algebra = Algebra::euclidean(3);
        let table = ProductTable::new(&algebra);
        let generator = TraitsGenerator::new(&spec, &algebra, table);

        let (tokens, _tests) = generator.generate_traits_file();
        let code = tokens.to_string();

        // Vector * Vector should produce a geometric product via Mul trait
        assert!(
            code.contains("Mul") && code.contains("for Vector"),
            "Expected Mul trait impl for Vector"
        );
    }

    #[test]
    fn symbolica_generates_wedge() {
        let spec = parse_spec(include_str!("../../algebras/euclidean3.toml")).unwrap();
        let algebra = Algebra::euclidean(3);
        let table = ProductTable::new(&algebra);
        let generator = TraitsGenerator::new(&spec, &algebra, table);

        let (tokens, _tests) = generator.generate_traits_file();
        let code = tokens.to_string();

        // Vector wedge Vector should produce Bivector via Wedge trait
        assert!(
            code.contains("Wedge") && code.contains("for Vector"),
            "Expected Wedge trait impl for Vector, got:\n{}",
            &code[..2000.min(code.len())]
        );
    }

    #[test]
    fn symbolica_generates_approx_impls() {
        let spec = parse_spec(include_str!("../../algebras/euclidean3.toml")).unwrap();
        let algebra = Algebra::euclidean(3);
        let table = ProductTable::new(&algebra);
        let generator = TraitsGenerator::new(&spec, &algebra, table);

        let (tokens, _tests) = generator.generate_traits_file();
        let code = tokens.to_string();

        assert!(code.contains("AbsDiffEq for Vector"));
        assert!(code.contains("RelativeEq for Vector"));
        assert!(code.contains("UlpsEq for Vector"));
    }

    #[test]
    fn symbolica_generates_arbitrary_impls() {
        let spec = parse_spec(include_str!("../../algebras/euclidean3.toml")).unwrap();
        let algebra = Algebra::euclidean(3);
        let table = ProductTable::new(&algebra);
        let generator = TraitsGenerator::new(&spec, &algebra, table);

        let (tokens, _tests) = generator.generate_traits_file();
        let code = tokens.to_string();

        assert!(code.contains("impl < T : Float + Debug + 'static > Arbitrary for Vector"));
    }
}
