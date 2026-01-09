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
| [PRD-1](prd-1-foundation.md) | **Complete** | Foundation: Float, Signature, Blade |
| [PRD-2](prd-2-multivector.md) | **Complete** | Core Multivector, geometric product |
| [PRD-3](prd-3-products.md) | **Complete** | Derived products, grade operations |
| [PRD-4](prd-4-specialized.md) | **Complete** | Optimized 2D/3D types with conversions |
| [PRD-5](prd-5-pga.md) | Pending | Projective GA (points, lines, planes, motors) |
| [PRD-6](prd-6-cga.md) | Pending | Conformal GA, documentation polish |
| [PRD-7](prd-7-traits.md) | Pending | Trait abstractions for generic algorithms |
| [PRD-8](prd-8-nalgebra.md) | Draft | nalgebra interoperability (feature-gated) |
| [PRD-9](prd-9-simd.md) | Draft | SIMD and simba integration (feature-gated) |

## Implementation Order

1. ~~PRD-1: Foundation Layer~~ **Done**
2. ~~PRD-2: Core Multivector~~ **Done**
3. ~~PRD-3: Derived Products~~ **Done**
4. ~~PRD-4: Specialized Types~~ **Done**
5. PRD-5: PGA
6. PRD-6: CGA + Polish
7. PRD-7: Trait Abstractions (can be done in parallel with PRD-5/6)
8. PRD-8: nalgebra Interop (can be done anytime, independent of other PRDs)
9. PRD-9: SIMD/simba (can be done anytime, independent of other PRDs)

## Verification

Each PRD should pass:
- `cargo check` - no warnings
- `cargo test` - all proptest properties pass
- `cargo clippy` - no lints
- `cargo doc` - documentation builds
- `cargo bench` - benchmarks captured in `benches/reports/`
