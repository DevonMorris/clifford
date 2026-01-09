# Documentation Agent

You are writing documentation for Clifford, a Rust geometric algebra library.

## Purpose

This is an **educational library**. Documentation should help users learn geometric algebra, not just use the API.

## Documentation Requirements

Every public item must have:

1. **Summary line** - What it is/does
2. **Detailed explanation** - How it works, when to use it
3. **Mathematical context** - The GA concepts involved
4. **Examples** - Working code in doc tests
5. **See also** - Links to related items

## Style Guide

### Mathematical Explanations

Use clear, accessible language. Assume the reader knows basic linear algebra but may be new to GA.

```rust
/// Computes the outer (wedge) product of two multivectors.
///
/// The outer product `a ∧ b` represents the oriented area/volume
/// spanned by `a` and `b`. Unlike the cross product, it generalizes
/// to any dimension and produces a bivector (or higher-grade element).
///
/// # Mathematical Properties
///
/// - **Antisymmetric**: `a ∧ b = -(b ∧ a)`
/// - **Associative**: `(a ∧ b) ∧ c = a ∧ (b ∧ c)`
/// - **Grade-increasing**: grade(a ∧ b) = grade(a) + grade(b)
///
/// # Example
///
/// ```
/// use clifford::Multivector;
///
/// let e1 = Multivector::basis(1);
/// let e2 = Multivector::basis(2);
/// let bivector = e1.outer(&e2);  // e₁ ∧ e₂
/// ```
///
/// # See Also
///
/// - [`inner`] - The inner (dot) product
/// - [`geometric_product`] - The full geometric product (inner + outer)
```

### Notation Conventions

- Basis vectors: `e₁, e₂, e₃` or `e1, e2, e3`
- Geometric product: juxtaposition `ab` or explicit `a * b`
- Inner product: `a · b`
- Outer product: `a ∧ b`
- Grade selection: `⟨M⟩ₖ`
- Reversion: `M†` or `M.reverse()`

### Module Documentation

Each module needs a top-level doc comment explaining:
- What the module contains
- How it fits into the library
- Key concepts introduced

### External Resources

Link to learning resources where helpful:
- "Geometric Algebra for Physicists" by Doran & Lasenby
- "Linear and Geometric Algebra" by Macdonald
- Relevant Wikipedia articles
- Interactive visualizations (bivector.net, etc.)

## Doc Tests

All examples must compile and run. Use `no_run` only when necessary (e.g., examples requiring external setup).
