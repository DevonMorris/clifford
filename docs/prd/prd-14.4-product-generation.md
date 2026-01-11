# PRD-14.4: Product Code Generation

**Status**: Draft
**Parent**: [PRD-14: Geometric Algebra Code Generator](prd-14-codegen.md)
**Goal**: Generate optimized code for all binary products between types

## Overview

This phase generates:
1. Geometric product functions for all type pairs
2. Inner, outer, and regressive product functions
3. Sandwich product functions (for transformations)
4. Scalar product (returns just the scalar part)

Each product is generated as an optimized, inlined function with only the non-zero terms.

## Product Types

### Geometric Product

The fundamental product. All other products are derived from it.

```rust
/// Geometric product: A × B → C
///
/// Computes the full geometric product, keeping all resulting grades.
pub fn geometric_a_b<T: Float>(a: A<T>, b: B<T>) -> C<T>
```

### Outer Product (Wedge)

Grade-raising product: `grade(a ∧ b) = grade(a) + grade(b)`.

```rust
/// Outer product: A ∧ B → C
///
/// Keeps only terms where grades add.
pub fn outer_a_b<T: Float>(a: A<T>, b: B<T>) -> C<T>
```

### Inner Product (Left Contraction)

Grade-lowering product: `grade(a ⌋ b) = grade(b) - grade(a)` (when `grade(a) ≤ grade(b)`).

```rust
/// Left contraction: A ⌋ B → C
pub fn left_contract_a_b<T: Float>(a: A<T>, b: B<T>) -> C<T>
```

### Regressive Product

Dual of outer: `a ∨ b = (a* ∧ b*)*` where `*` is the dual.

```rust
/// Regressive product: A ∨ B → C
pub fn regressive_a_b<T: Float>(a: A<T>, b: B<T>) -> C<T>
```

### Sandwich Product

`sandwich(v, x) = v * x * reverse(v)` or `reverse(v) * x * v` depending on convention.

```rust
/// Sandwich product: V × X × V† → X
///
/// Applies V as a transformation to X.
pub fn sandwich_v_x<T: Float>(v: V<T>, x: X<T>) -> X<T>
```

### Scalar Product

Returns only the grade-0 part of the geometric product.

```rust
/// Scalar product: A × B → T
///
/// Returns the scalar part of the geometric product.
pub fn scalar_a_b<T: Float>(a: A<T>, b: B<T>) -> T
```

## Generation Strategy

### Step 1: Compute Contributing Terms

For each (input_type_A, input_type_B, output_blade) triple, find all pairs of input blades that contribute:

```rust
struct ProductTerm {
    /// Sign of this contribution (+1 or -1)
    sign: i8,
    /// Field name from type A
    a_field: String,
    /// Field name from type B
    b_field: String,
}

fn compute_terms(
    a_blades: &[usize],
    b_blades: &[usize],
    result_blade: usize,
    table: &ProductTable,
) -> Vec<ProductTerm> {
    let mut terms = Vec::new();

    for &a in a_blades {
        for &b in b_blades {
            let (sign, result) = table.geometric(a, b);
            if result == result_blade && sign != 0 {
                terms.push(ProductTerm {
                    sign,
                    a_field: blade_to_field(a),
                    b_field: blade_to_field(b),
                });
            }
        }
    }

    terms
}
```

### Step 2: Generate Expression

Convert terms to a Rust expression:

```rust
fn generate_expression(terms: &[ProductTerm]) -> TokenStream {
    if terms.is_empty() {
        return quote! { T::zero() };
    }

    let mut expr_parts: Vec<TokenStream> = Vec::new();

    for (i, term) in terms.iter().enumerate() {
        let a = format_ident!("{}", term.a_field);
        let b = format_ident!("{}", term.b_field);

        let term_expr = match (i, term.sign) {
            (0, 1)  => quote! { a.#a() * b.#b() },
            (0, -1) => quote! { -(a.#a() * b.#b()) },
            (_, 1)  => quote! { + a.#a() * b.#b() },
            (_, -1) => quote! { - a.#a() * b.#b() },
            _ => unreachable!(),
        };

        expr_parts.push(term_expr);
    }

    quote! { #(#expr_parts)* }
}
```

### Step 3: Determine Output Type

The output type depends on which grades are produced:

```rust
fn determine_output_type(
    a_grades: &[usize],
    b_grades: &[usize],
    product_kind: ProductKind,
    dim: usize,
    type_map: &HashMap<Vec<usize>, String>,
) -> String {
    let output_grades = match product_kind {
        ProductKind::Geometric => {
            compute_geometric_grades(a_grades, b_grades, dim)
        }
        ProductKind::Outer => {
            compute_outer_grades(a_grades, b_grades, dim)
        }
        ProductKind::Inner => {
            compute_inner_grades(a_grades, b_grades)
        }
        ProductKind::Regressive => {
            compute_regressive_grades(a_grades, b_grades, dim)
        }
    };

    // Find matching type
    type_map.get(&output_grades)
        .cloned()
        .unwrap_or_else(|| panic!("No type for grades {:?}", output_grades))
}
```

## Product Generator

```rust
pub struct ProductGenerator<'a> {
    spec: &'a AlgebraSpec,
    algebra: &'a Algebra,
    table: ProductTable,
}

impl<'a> ProductGenerator<'a> {
    pub fn generate_products_file(&self) -> TokenStream {
        let imports = self.generate_imports();
        let geometric = self.generate_all_geometric();
        let outer = self.generate_all_outer();
        let inner = self.generate_all_inner();
        let regressive = self.generate_all_regressive();
        let sandwich = self.generate_all_sandwich();
        let scalar = self.generate_all_scalar();

        quote! {
            //! Product operations for {algebra_name}.
            //!
            //! Auto-generated by clifford-codegen.

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
            // Regressive Products
            // ============================================================
            #regressive

            // ============================================================
            // Sandwich Products
            // ============================================================
            #sandwich

            // ============================================================
            // Scalar Products
            // ============================================================
            #scalar
        }
    }

    fn generate_geometric_product(
        &self,
        type_a: &TypeSpec,
        type_b: &TypeSpec,
    ) -> TokenStream {
        let a_name = format_ident!("{}", type_a.name);
        let b_name = format_ident!("{}", type_b.name);

        // Determine output type and blades
        let output_grades = compute_geometric_grades(&type_a.grades, &type_b.grades, self.algebra.dim());
        let output_type = self.find_type_for_grades(&output_grades)?;
        let c_name = format_ident!("{}", output_type.name);

        let fn_name = format_ident!(
            "geometric_{}_{}",
            type_a.name.to_lowercase(),
            type_b.name.to_lowercase()
        );

        // Generate each output field
        let field_exprs: Vec<TokenStream> = output_type.fields.iter().map(|field| {
            let terms = compute_terms(
                &self.blades_for_type(type_a),
                &self.blades_for_type(type_b),
                field.blade_index,
                &self.table,
            );
            generate_expression(&terms)
        }).collect();

        let field_names: Vec<_> = output_type.fields.iter()
            .map(|f| format_ident!("{}", f.name))
            .collect();

        let doc = format!(
            "Geometric product: {} × {} → {}\n\n\
            Computes the full geometric product.",
            type_a.name, type_b.name, output_type.name
        );

        quote! {
            #[doc = #doc]
            #[inline]
            pub fn #fn_name<T: Float>(a: #a_name<T>, b: #b_name<T>) -> #c_name<T> {
                #c_name::new(#(#field_exprs),*)
            }
        }
    }

    fn generate_sandwich_product(
        &self,
        versor_type: &TypeSpec,
        operand_type: &TypeSpec,
    ) -> TokenStream {
        let v_name = format_ident!("{}", versor_type.name);
        let x_name = format_ident!("{}", operand_type.name);

        let fn_name = format_ident!(
            "sandwich_{}_{}",
            versor_type.name.to_lowercase(),
            operand_type.name.to_lowercase()
        );

        // Sandwich product: v * x * reverse(v)
        // But we want to expand this to a single optimized formula

        // First, compute v * x (intermediate)
        // Then, compute (v * x) * reverse(v) (final)
        // Collect only terms that survive in the final operand type

        let field_exprs: Vec<TokenStream> = operand_type.fields.iter().map(|field| {
            let terms = self.compute_sandwich_terms(
                versor_type,
                operand_type,
                field.blade_index,
            );
            self.generate_sandwich_expression(&terms, versor_type)
        }).collect();

        let field_names: Vec<_> = operand_type.fields.iter()
            .map(|f| format_ident!("{}", f.name))
            .collect();

        let doc = format!(
            "Sandwich product: {} × {} × {}† → {}\n\n\
            Applies the versor transformation.",
            versor_type.name, operand_type.name, versor_type.name, operand_type.name
        );

        quote! {
            #[doc = #doc]
            #[inline]
            pub fn #fn_name<T: Float>(v: #v_name<T>, x: #x_name<T>) -> #x_name<T> {
                #x_name::new(#(#field_exprs),*)
            }
        }
    }

    /// Computes terms for sandwich product v * x * rev(v) -> result_blade
    fn compute_sandwich_terms(
        &self,
        versor_type: &TypeSpec,
        operand_type: &TypeSpec,
        result_blade: usize,
    ) -> Vec<SandwichTerm> {
        let v_blades = self.blades_for_type(versor_type);
        let x_blades = self.blades_for_type(operand_type);

        let mut terms = Vec::new();

        // For each combination v_i * x_j * rev(v_k)
        for &v_i in &v_blades {
            for &x_j in &x_blades {
                for &v_k in &v_blades {
                    // Compute (v_i * x_j) first
                    let (sign_vx, vx) = self.table.geometric(v_i, x_j);
                    if sign_vx == 0 { continue; }

                    // Then (v_i * x_j) * rev(v_k)
                    // rev(v_k) has sign based on grade
                    let grade_k = Blade::from_index(v_k).grade();
                    let rev_sign = if (grade_k * (grade_k.saturating_sub(1)) / 2) % 2 == 0 { 1 } else { -1 };

                    let (sign_vxr, result) = self.table.geometric(vx, v_k);
                    if sign_vxr == 0 { continue; }

                    if result == result_blade {
                        let final_sign = sign_vx * sign_vxr * rev_sign;
                        terms.push(SandwichTerm {
                            sign: final_sign,
                            v_field_1: blade_to_field(v_i),
                            x_field: blade_to_field(x_j),
                            v_field_2: blade_to_field(v_k),
                        });
                    }
                }
            }
        }

        // Combine like terms
        self.simplify_sandwich_terms(terms)
    }
}

struct SandwichTerm {
    sign: i8,
    v_field_1: String,
    x_field: String,
    v_field_2: String,
}
```

## Optimization: Common Subexpression Elimination

For sandwich products, extract common subexpressions:

```rust
fn optimize_sandwich(terms: &[SandwichTerm]) -> TokenStream {
    // Find terms that share v_field_1 and x_field
    // These can use a common intermediate variable

    let mut groups: HashMap<(String, String), Vec<&SandwichTerm>> = HashMap::new();
    for term in terms {
        let key = (term.v_field_1.clone(), term.x_field.clone());
        groups.entry(key).or_default().push(term);
    }

    let mut intermediates = Vec::new();
    let mut final_terms = Vec::new();

    for ((v1, x), group) in groups {
        if group.len() > 1 {
            // Create intermediate variable
            let inter_name = format_ident!("q_{}_{}", v1, x);
            let v1_ident = format_ident!("{}", v1);
            let x_ident = format_ident!("{}", x);

            intermediates.push(quote! {
                let #inter_name = v.#v1_ident() * x.#x_ident();
            });

            // Use intermediate in final expression
            for term in group {
                let v2 = format_ident!("{}", term.v_field_2);
                let sign = if term.sign > 0 { quote! { + } } else { quote! { - } };
                final_terms.push(quote! { #sign #inter_name * v.#v2() });
            }
        } else {
            // Direct computation
            let term = group[0];
            let v1_ident = format_ident!("{}", term.v_field_1);
            let x_ident = format_ident!("{}", term.x_field);
            let v2_ident = format_ident!("{}", term.v_field_2);
            let sign = if term.sign > 0 { quote! { + } } else { quote! { - } };
            final_terms.push(quote! { #sign v.#v1_ident() * x.#x_ident() * v.#v2_ident() });
        }
    }

    quote! {
        {
            #(#intermediates)*
            #(#final_terms)*
        }
    }
}
```

## Handling Output Type Selection

When the output type is ambiguous (product produces grades matching multiple types):

```rust
enum OutputTypeResolution {
    /// Exact match found
    Exact(TypeSpec),
    /// Multiple matches - use most specific
    Ambiguous(Vec<TypeSpec>),
    /// No match - would need new type
    NoMatch(Vec<usize>),
}

fn resolve_output_type(
    output_grades: &[usize],
    spec: &AlgebraSpec,
) -> OutputTypeResolution {
    let matches: Vec<_> = spec.types.iter()
        .filter(|t| t.grades == output_grades)
        .cloned()
        .collect();

    match matches.len() {
        0 => OutputTypeResolution::NoMatch(output_grades.to_vec()),
        1 => OutputTypeResolution::Exact(matches[0].clone()),
        _ => OutputTypeResolution::Ambiguous(matches),
    }
}
```

## Generated Examples

### Geometric: Vector × Vector → Even

```rust
/// Geometric product: Vector × Vector → Even
#[inline]
pub fn geometric_vector_vector<T: Float>(a: Vector<T>, b: Vector<T>) -> Even<T> {
    Even::new(
        // Scalar: a·b
        a.x()*b.x() + a.y()*b.y() + a.z()*b.z(),
        // e12
        a.x()*b.y() - a.y()*b.x(),
        // e13
        a.x()*b.z() - a.z()*b.x(),
        // e23
        a.y()*b.z() - a.z()*b.y(),
    )
}
```

### Outer: Vector ∧ Vector → Bivector

```rust
/// Outer product: Vector ∧ Vector → Bivector
#[inline]
pub fn outer_vector_vector<T: Float>(a: Vector<T>, b: Vector<T>) -> Bivector<T> {
    Bivector::new(
        a.x()*b.y() - a.y()*b.x(),  // e12
        a.x()*b.z() - a.z()*b.x(),  // e13
        a.y()*b.z() - a.z()*b.y(),  // e23
    )
}
```

### Sandwich: Rotor × Vector → Vector

```rust
/// Sandwich product: Rotor × Vector × Rotor† → Vector
#[inline]
pub fn sandwich_rotor_vector<T: Float>(r: Rotor<T>, v: Vector<T>) -> Vector<T> {
    // Optimized formula (not naive r * v * rev(r))
    let s = r.s();
    let xy = r.xy();
    let xz = r.xz();
    let yz = r.yz();

    // Intermediate: q = rev(r) * v
    let qx = s * v.x() - xy * v.y() - xz * v.z();
    let qy = s * v.y() + xy * v.x() - yz * v.z();
    let qz = s * v.z() + xz * v.x() + yz * v.y();
    let qt = -xy * v.z() + xz * v.y() - yz * v.x();

    // Result: q * r
    Vector::new(
        s * qx - xy * qy - xz * qz - yz * qt,
        s * qy + xy * qx + xz * qt - yz * qz,
        s * qz - xy * qt + xz * qx + yz * qy,
    )
}
```

## Testing

### Against Multivector

```rust
#[cfg(test)]
mod tests {
    use super::*;
    use proptest::prelude::*;
    use clifford::{Multivector, signature::Euclidean3};

    proptest! {
        #[test]
        fn geometric_vector_vector_matches_mv(
            a in any::<Vector<f64>>(),
            b in any::<Vector<f64>>(),
        ) {
            let gen_result = geometric_vector_vector(a, b);

            let mv_a = Multivector::<f64, Euclidean3>::from(a);
            let mv_b = Multivector::<f64, Euclidean3>::from(b);
            let mv_result = mv_a * mv_b;

            let mv_as_even = Even::try_from(mv_result).unwrap();

            prop_assert!(abs_diff_eq!(gen_result, mv_as_even, epsilon = 1e-10));
        }

        #[test]
        fn sandwich_rotor_vector_matches_mv(
            r in any::<UnitRotor<f64>>(),
            v in any::<Vector<f64>>(),
        ) {
            let gen_result = sandwich_rotor_vector(*r, v);

            let mv_r = Multivector::<f64, Euclidean3>::from(*r);
            let mv_v = Multivector::<f64, Euclidean3>::from(v);
            let mv_rev_r = mv_r.reverse();
            let mv_result = &(&mv_rev_r * &mv_v) * &mv_r;

            let mv_as_vec = Vector::try_from(mv_result).unwrap();

            prop_assert!(abs_diff_eq!(gen_result, mv_as_vec, epsilon = 1e-10));
        }
    }
}
```

### Algebraic Properties

```rust
proptest! {
    #[test]
    fn geometric_product_associative(
        a in any::<Rotor<f64>>(),
        b in any::<Rotor<f64>>(),
        c in any::<Rotor<f64>>(),
    ) {
        let ab_c = geometric_rotor_rotor(geometric_rotor_rotor(a, b), c);
        let a_bc = geometric_rotor_rotor(a, geometric_rotor_rotor(b, c));
        prop_assert!(abs_diff_eq!(ab_c, a_bc, epsilon = 1e-10));
    }

    #[test]
    fn outer_product_anticommutative(
        a in any::<Vector<f64>>(),
        b in any::<Vector<f64>>(),
    ) {
        let ab = outer_vector_vector(a, b);
        let ba = outer_vector_vector(b, a);
        prop_assert!(abs_diff_eq!(ab, -ba, epsilon = 1e-10));
    }

    #[test]
    fn sandwich_preserves_norm(
        r in any::<UnitRotor<f64>>(),
        v in any::<Vector<f64>>(),
    ) {
        let result = sandwich_rotor_vector(*r, v);
        prop_assert!(abs_diff_eq!(v.norm(), result.norm(), epsilon = 1e-10));
    }
}
```

## Deliverables

- [ ] `ProductGenerator` struct
- [ ] Geometric product generation
- [ ] Outer product generation
- [ ] Inner product (left contraction) generation
- [ ] Regressive product generation
- [ ] Sandwich product generation
- [ ] Scalar product generation
- [ ] Common subexpression optimization for sandwich
- [ ] Output type resolution
- [ ] Property tests against Multivector
- [ ] Algebraic property tests

## Success Criteria

1. All products compile and produce correct results
2. Products match `Multivector` implementation
3. Sandwich products are optimized (not naive 3-way product)
4. Algebraic properties (associativity, anticommutativity) verified
