# Clifford PRDs (Product Requirements Documents)

Implementation plan for the Clifford Geometric Algebra library.

## Overview

- **Dimensions**: Optimized 2D, 3D, 4D implementations
- **Signatures**: Full support (Euclidean, Minkowski, PGA, CGA, arbitrary)
- **Type System**: Compile-time const generics
- **Scalars**: Generic over float types (f32, f64)

## PRDs

| PRD | Status | Description |
|-----|--------|-------------|
| [PRD-0](prd-0-architecture.md) | Reference | Architecture overview |
| [PRD-1](prd-1-foundation.md) | Pending | Foundation: Float, Signature, basis |
| [PRD-2](prd-2-multivector.md) | Pending | Core Multivector, geometric product |
| [PRD-3](prd-3-products.md) | Pending | Derived products, grade operations |
| [PRD-4](prd-4-specialized.md) | Pending | Optimized 2D/3D types |
| [PRD-5](prd-5-pga.md) | Pending | Projective GA |
| [PRD-6](prd-6-cga.md) | Pending | Conformal GA, polish |

## Implementation Order

1. PRD-1: Foundation Layer
2. PRD-2: Core Multivector
3. PRD-3: Derived Products
4. PRD-4: Specialized Types
5. PRD-5: PGA
6. PRD-6: CGA + Polish

## Verification

Each PRD should pass:
- `cargo check` - no warnings
- `cargo test` - all proptest properties pass
- `cargo clippy` - no lints
- `cargo doc` - documentation builds
