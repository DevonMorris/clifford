# PRD-0: Architecture Overview

## Design Principles

### 1. Layered Abstraction

- **Foundation Layer**: Generic algebraic primitives with full generics
- **Specialized Layer**: Optimized concrete types for common dimensions (2D, 3D, 4D)
- **Application Layer**: High-level types (Rotors, Motors) with geometric meaning

### 2. Type System Strategy

- `Signature` trait encodes metric as `(P, Q, R)`:
  - P = basis vectors squaring to +1
  - Q = basis vectors squaring to -1
  - R = basis vectors squaring to 0 (degenerate)
- `Multivector<T, S>` generic over scalar `T` and signature `S`
- Const generics for compile-time dimension

### 3. Storage Strategy

Dense array representation with 2^N coefficients for dimensions ≤ 6 (64 elements max).

## Module Structure

```
src/
├── lib.rs                    # Crate root, re-exports
│
├── scalar/                   # Scalar type abstractions
│   ├── mod.rs
│   └── float.rs              # Float trait definition
│
├── signature/                # Metric signature system
│   ├── mod.rs
│   ├── metric.rs             # Signature trait
│   ├── euclidean.rs          # Euclidean signatures (n,0,0)
│   ├── pga.rs                # Projective GA (n,0,1)
│   └── cga.rs                # Conformal GA (n+1,1,0)
│
├── basis/                    # Basis blade system
│   ├── mod.rs
│   ├── blade.rs              # Basis product computation
│   ├── grade.rs              # Grade utilities
│   └── index.rs              # Canonical ordering
│
├── algebra/                  # Core algebraic types
│   ├── mod.rs
│   ├── multivector.rs        # Multivector<T, S>
│   ├── ops/
│   │   ├── mod.rs
│   │   ├── geometric.rs      # Geometric product
│   │   ├── inner.rs          # Inner product
│   │   ├── outer.rs          # Outer product
│   │   ├── regressive.rs     # Regressive product
│   │   └── sandwich.rs       # Sandwich product
│   ├── unary.rs              # Reverse, conjugate, dual
│   └── norms.rs              # Norms and normalization
│
├── specialized/              # Dimension-specific optimizations
│   ├── mod.rs
│   ├── euclidean::dim2/                 # 2D Euclidean GA
│   ├── euclidean::dim3/                 # 3D Euclidean GA
│   ├── peuclidean::dim3/                # 3D Projective GA
│   └── ceuclidean::dim3/                # 3D Conformal GA
│
└── transforms/               # High-level geometric types
    ├── mod.rs
    ├── rotor.rs              # Rotation
    └── motor.rs              # Rigid transformation
```

## Key Design Decisions

### Const Generics Approach

Use bounded arrays with `MAX_DIM = 6` (64 basis blades) rather than `generic_const_exprs` which is unstable.

```rust
pub const MAX_DIM: usize = 6;
pub const MAX_BLADES: usize = 1 << MAX_DIM; // 64

pub struct Multivector<T: Float, S: Signature> {
    coeffs: [T; MAX_BLADES],
    _signature: PhantomData<S>,
}
```

### Operator Overloading

- `*` for geometric product (most fundamental)
- Methods for other products: `.outer()`, `.inner()`, `.regressive()`
- `+`, `-` for addition/subtraction
- `/` for scalar division

### Generic Floats

All types are generic over a `Float` trait that abstracts f32/f64.

```rust
pub trait Float: Copy + Clone + Debug + ... {
    const ZERO: Self;
    const ONE: Self;
    fn sqrt(self) -> Self;
    fn sin(self) -> Self;
    fn cos(self) -> Self;
}
```
