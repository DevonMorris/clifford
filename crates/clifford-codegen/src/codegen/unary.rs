//! Unary operation code generation.
//!
//! This module generates unary operations like reverse, antireverse, and complement.

use proc_macro2::TokenStream;
use quote::{format_ident, quote};

use crate::algebra::{Algebra, ProductTable};
use crate::spec::{AlgebraSpec, TypeSpec};

/// Generates Rust code for unary operations.
pub struct UnaryGenerator<'a> {
    /// The algebra specification.
    spec: &'a AlgebraSpec,
    /// Product table for complement computation.
    table: ProductTable,
}

impl<'a> UnaryGenerator<'a> {
    /// Creates a new unary generator.
    pub fn new(spec: &'a AlgebraSpec) -> Self {
        let algebra = Algebra::new(spec.signature.p, spec.signature.q, spec.signature.r);
        let table = ProductTable::new(&algebra);
        Self { spec, table }
    }

    /// Generates all unary operations.
    pub fn generate_all(&self) -> TokenStream {
        let reverse = self.generate_all_reverse();
        let antireverse = self.generate_all_antireverse();
        let complement = self.generate_all_complement();

        quote! {
            // ============================================================
            // Reverse Operations
            // ============================================================
            #reverse

            // ============================================================
            // Antireverse Operations
            // ============================================================
            #antireverse

            // ============================================================
            // Complement Operations
            // ============================================================
            #complement
        }
    }

    /// Generates all reverse operations.
    fn generate_all_reverse(&self) -> TokenStream {
        let ops: Vec<TokenStream> = self
            .spec
            .types
            .iter()
            .filter(|t| t.alias_of.is_none())
            .filter_map(|ty| self.generate_reverse(ty))
            .collect();

        quote! { #(#ops)* }
    }

    /// Generates a reverse operation for a type.
    fn generate_reverse(&self, ty: &TypeSpec) -> Option<TokenStream> {
        let type_name = format_ident!("{}", ty.name);
        let fn_name = format_ident!("reverse_{}", ty.name.to_lowercase());

        // Build field expressions with reverse signs
        let field_exprs: Vec<TokenStream> = ty
            .fields
            .iter()
            .map(|field| {
                let field_name = format_ident!("{}", field.name);
                let grade = field.grade;
                // Reverse sign: (-1)^(k(k-1)/2)
                if (grade * grade.saturating_sub(1) / 2).is_multiple_of(2) {
                    quote! { a.#field_name() }
                } else {
                    quote! { -a.#field_name() }
                }
            })
            .collect();

        let constructor = quote! { #type_name::new(#(#field_exprs),*) };

        let doc = format!(
            "Reverses the {} (negates grades where k(k-1)/2 is odd).",
            ty.name
        );

        Some(quote! {
            #[doc = #doc]
            #[inline]
            pub fn #fn_name<T: Float>(a: &#type_name<T>) -> #type_name<T> {
                #constructor
            }
        })
    }

    /// Generates all antireverse operations.
    fn generate_all_antireverse(&self) -> TokenStream {
        let ops: Vec<TokenStream> = self
            .spec
            .types
            .iter()
            .filter(|t| t.alias_of.is_none())
            .filter_map(|ty| self.generate_antireverse(ty))
            .collect();

        quote! { #(#ops)* }
    }

    /// Generates an antireverse operation for a type.
    fn generate_antireverse(&self, ty: &TypeSpec) -> Option<TokenStream> {
        let type_name = format_ident!("{}", ty.name);
        let fn_name = format_ident!("antireverse_{}", ty.name.to_lowercase());
        let dim = self.spec.signature.dim();

        // Build field expressions with antireverse signs
        let field_exprs: Vec<TokenStream> = ty
            .fields
            .iter()
            .map(|field| {
                let field_name = format_ident!("{}", field.name);
                let grade = field.grade;
                let antigrade = dim - grade;
                // Antireverse sign: (-1)^((n-k)(n-k-1)/2)
                if (antigrade * antigrade.saturating_sub(1) / 2).is_multiple_of(2) {
                    quote! { a.#field_name() }
                } else {
                    quote! { -a.#field_name() }
                }
            })
            .collect();

        let constructor = quote! { #type_name::new(#(#field_exprs),*) };

        let doc = format!(
            "Antireverses the {} (negates grades where (n-k)(n-k-1)/2 is odd).",
            ty.name
        );

        Some(quote! {
            #[doc = #doc]
            #[inline]
            pub fn #fn_name<T: Float>(a: &#type_name<T>) -> #type_name<T> {
                #constructor
            }
        })
    }

    /// Generates all complement operations.
    fn generate_all_complement(&self) -> TokenStream {
        let ops: Vec<TokenStream> = self
            .spec
            .types
            .iter()
            .filter(|t| t.alias_of.is_none())
            .filter_map(|ty| self.generate_complement(ty))
            .collect();

        quote! { #(#ops)* }
    }

    /// Generates a complement operation for a type.
    ///
    /// The complement maps each blade to its dual blade such that blade ^ complement = pseudoscalar.
    fn generate_complement(&self, ty: &TypeSpec) -> Option<TokenStream> {
        let dim = self.spec.signature.dim();

        // Find the output type (grades are complemented: k -> n-k)
        let mut output_grades: Vec<usize> = ty.grades.iter().map(|&g| dim - g).collect();
        output_grades.sort();

        // Find a type with matching grades (compare sorted since grade order may differ)
        let output_type = self.spec.types.iter().find(|t| {
            if t.alias_of.is_some() {
                return false;
            }
            let mut t_grades = t.grades.clone();
            t_grades.sort();
            t_grades == output_grades
        })?;

        let type_name = format_ident!("{}", ty.name);
        let output_name = format_ident!("{}", output_type.name);
        let fn_name = format_ident!("complement_{}", ty.name.to_lowercase());

        // Build field expressions using complement table
        let field_exprs: Vec<TokenStream> = output_type
            .fields
            .iter()
            .map(|out_field| {
                // Find which input field maps to this output blade
                let out_blade = out_field.blade_index;

                // Find the input blade that complements to this output blade
                let mut expr = quote! { T::zero() };
                for in_field in &ty.fields {
                    let (sign, comp_blade) = self.table.complement(in_field.blade_index);
                    if comp_blade == out_blade && sign != 0 {
                        let in_name = format_ident!("{}", in_field.name);
                        if sign > 0 {
                            expr = quote! { a.#in_name() };
                        } else {
                            expr = quote! { -a.#in_name() };
                        }
                        break;
                    }
                }
                expr
            })
            .collect();

        let constructor = quote! { #output_name::new(#(#field_exprs),*) };

        let doc = format!(
            "Computes the right complement of {} -> {}.",
            ty.name, output_type.name
        );

        Some(quote! {
            #[doc = #doc]
            #[inline]
            pub fn #fn_name<T: Float>(a: &#type_name<T>) -> #output_name<T> {
                #constructor
            }
        })
    }
}
