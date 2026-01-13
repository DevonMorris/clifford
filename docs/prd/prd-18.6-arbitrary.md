# PRD-18.6: Bespoke Arbitrary Implementations

**Status**: Draft
**Parent**: [PRD-18](prd-18-constraint-redesign.md)
**Goal**: Generate Arbitrary implementations that solve inferred constraints

## Overview

For property-based testing, `Arbitrary` implementations should generate valid instances. For constrained types, this means solving the constraint during generation rather than filtering invalid values.

**Current approach (inefficient):**
```rust
impl Arbitrary for Motor<f64> {
    fn arbitrary_with(_: ()) -> BoxedStrategy<Self> {
        // Generate random values and filter out constraint violations
        (/* 8 random fields */)
            .prop_filter_map("constraint", |(...)| Motor::new(...))
            .boxed()
    }
}
```

**New approach (efficient):**
```rust
impl Arbitrary for Motor<f64> {
    fn arbitrary_with(_: ()) -> BoxedStrategy<Self> {
        // Generate 7 independent params, solve for e0123
        (/* 7 random fields */)
            .prop_map(|(s, e23, e31, e12, e01, e02, e03)| {
                let e0123 = -(e23 * e01 + e31 * e02 + e12 * e03) / s;
                Motor::new_unchecked(s, e23, e31, e12, e01, e02, e03, e0123)
            })
            .boxed()
    }
}
```

## Deliverables

### Modified Files (clifford-codegen)

#### `src/codegen/traits.rs`

Generate bespoke Arbitrary implementations:

```rust
fn generate_arbitrary(type_spec: &TypeSpec) -> TokenStream {
    let type_name = &type_spec.name;
    let fields = &type_spec.fields;

    if let Some(constraint) = &type_spec.inferred_constraint {
        generate_constrained_arbitrary(type_spec, constraint)
    } else {
        generate_unconstrained_arbitrary(type_spec)
    }
}

fn generate_unconstrained_arbitrary(type_spec: &TypeSpec) -> TokenStream {
    let type_name = &type_spec.name;
    let n_fields = type_spec.fields.len();

    // Generate tuple of n random f64 values
    let range_tuple = (0..n_fields).map(|_| quote!(-100.0f64..100.0));
    let field_names = type_spec.fields.iter().map(|f| &f.name);

    quote! {
        impl<T: Float + Debug + 'static> Arbitrary for #type_name<T> {
            type Parameters = ();
            type Strategy = BoxedStrategy<Self>;

            fn arbitrary_with(_: Self::Parameters) -> Self::Strategy {
                (#(#range_tuple),*)
                    .prop_map(|(#(#field_names),*)| {
                        #type_name::new_unchecked(
                            #(T::from_f64(#field_names)),*
                        )
                    })
                    .boxed()
            }
        }
    }
}

fn generate_constrained_arbitrary(
    type_spec: &TypeSpec,
    constraint: &InferredConstraint,
) -> TokenStream {
    let type_name = &type_spec.name;
    let fields = &type_spec.fields;

    // Determine which field to solve for
    let solve_for = determine_solve_for_field(constraint, fields);
    let independent_fields: Vec<_> = fields.iter()
        .filter(|f| f.name != solve_for)
        .collect();

    // Generate solution expression
    let solution_expr = generate_constraint_solution(constraint, &solve_for);

    // For quadratic constraints, filter for valid domain
    let needs_filter = constraint.constraint_type == ConstraintType::Quadratic;

    let n_independent = independent_fields.len();
    let range_tuple = (0..n_independent).map(|_| {
        if needs_filter {
            quote!((-0.5f64..0.5))  // Smaller range for quadratic
        } else {
            quote!((-10.0f64..10.0))
        }
    });
    let ind_names: Vec<_> = independent_fields.iter().map(|f| &f.name).collect();

    quote! {
        impl<T: Float + Debug + 'static> Arbitrary for #type_name<T> {
            type Parameters = ();
            type Strategy = BoxedStrategy<Self>;

            fn arbitrary_with(_: Self::Parameters) -> Self::Strategy {
                (#(#range_tuple),*)
                    .prop_filter("non-zero divisor", |(#(#ind_names),*)| {
                        // Filter out values that would cause division by zero
                        #divisor_check
                    })
                    .prop_map(|(#(#ind_names),*)| {
                        let #solve_for = #solution_expr;
                        #type_name::new_unchecked(
                            #(T::from_f64(#all_field_values)),*
                        )
                    })
                    .boxed()
            }
        }
    }
}
```

#### Generate Arbitrary for Wrapper Types

```rust
fn generate_wrapper_arbitrary() -> TokenStream {
    quote! {
        impl<T> Arbitrary for Unit<T>
        where
            T: Normed + Arbitrary + Debug + 'static,
            T::Strategy: 'static,
        {
            type Parameters = ();
            type Strategy = BoxedStrategy<Self>;

            fn arbitrary_with(_: Self::Parameters) -> Self::Strategy {
                any::<T>()
                    .prop_filter_map("normalizable", |t| Unit::try_new(t))
                    .boxed()
            }
        }

        impl<T> Arbitrary for Bulk<T>
        where
            T: DegenerateNormed + Arbitrary + Debug + 'static,
            T::Strategy: 'static,
        {
            type Parameters = ();
            type Strategy = BoxedStrategy<Self>;

            fn arbitrary_with(_: Self::Parameters) -> Self::Strategy {
                any::<T>()
                    .prop_filter_map("bulk-normalizable", |t| Bulk::try_new(t))
                    .boxed()
            }
        }

        // Similar for Ideal<T> and Proper<T>
    }
}
```

### Generated Code Example

**Motor (Study constraint - linear):**
```rust
impl<T: Float + Debug + 'static> Arbitrary for Motor<T> {
    type Parameters = ();
    type Strategy = BoxedStrategy<Self>;

    fn arbitrary_with(_: Self::Parameters) -> Self::Strategy {
        (
            (-10.0f64..10.0).prop_filter("non-zero", |s| s.abs() > 0.1),
            -10.0f64..10.0,
            -10.0f64..10.0,
            -10.0f64..10.0,
            -10.0f64..10.0,
            -10.0f64..10.0,
            -10.0f64..10.0,
        )
            .prop_map(|(s, e23, e31, e12, e01, e02, e03)| {
                // Solve Study condition: e0123 = -(e23*e01 + e31*e02 + e12*e03) / s
                let e0123 = -(e23 * e01 + e31 * e02 + e12 * e03) / s;
                Motor::new_unchecked(
                    T::from_f64(s),
                    T::from_f64(e23),
                    T::from_f64(e31),
                    T::from_f64(e12),
                    T::from_f64(e01),
                    T::from_f64(e02),
                    T::from_f64(e03),
                    T::from_f64(e0123),
                )
            })
            .boxed()
    }
}
```

**Rotor (Unit norm constraint - quadratic):**
```rust
impl<T: Float + Debug + 'static> Arbitrary for Rotor<T> {
    type Parameters = ();
    type Strategy = BoxedStrategy<Self>;

    fn arbitrary_with(_: Self::Parameters) -> Self::Strategy {
        (-0.5f64..0.5, -0.5f64..0.5, -0.5f64..0.5)
            .prop_filter_map("valid domain", |(xy, xz, yz)| {
                let sum_sq = xy * xy + xz * xz + yz * yz;
                if sum_sq > 1.0 { return None; }
                let s = (1.0 - sum_sq).sqrt();
                Some(Rotor::new_unchecked(
                    T::from_f64(s),
                    T::from_f64(xy),
                    T::from_f64(xz),
                    T::from_f64(yz),
                ))
            })
            .boxed()
    }
}
```

## Testing

```rust
proptest! {
    #[test]
    fn arbitrary_motor_satisfies_study(m in any::<Motor<f64>>()) {
        let residual = m.s() * m.e0123()
            + m.e23() * m.e01()
            + m.e31() * m.e02()
            + m.e12() * m.e03();
        prop_assert!(relative_eq!(residual, 0.0, epsilon = 1e-10, max_relative = 1e-10));
    }

    #[test]
    fn arbitrary_rotor_has_unit_norm(r in any::<Rotor<f64>>()) {
        let norm_sq = r.s() * r.s() + r.xy() * r.xy() + r.xz() * r.xz() + r.yz() * r.yz();
        prop_assert!(relative_eq!(norm_sq, 1.0, epsilon = 1e-10, max_relative = 1e-10));
    }

    #[test]
    fn arbitrary_unit_motor_has_unit_norm(m in any::<Unit<Motor<f64>>>()) {
        prop_assert!(relative_eq!(m.norm(), 1.0, epsilon = 1e-10, max_relative = 1e-10));
    }
}
```

## Success Criteria

1. Constrained types generate values that satisfy constraints
2. Linear constraints solved exactly (no filtering needed)
3. Quadratic constraints use filtered smaller range
4. Wrapper types (`Unit`, `Bulk`, etc.) have Arbitrary
5. All property tests pass reliably
6. No significant slowdown from constraint solving
