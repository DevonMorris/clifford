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
| [PRD-8](prd-8-nalgebra.md) | **Complete** | nalgebra interoperability (feature-gated) |
| [PRD-9](prd-9-simd.md) | Draft | SIMD and simba integration (feature-gated) |
| [PRD-10](prd-10-nalgebra-benchmarks.md) | **Complete** | nalgebra vs clifford benchmarks |
| [PRD-11](prd-11-autodiff.md) | Pending | Automatic differentiation via dual numbers |
| [PRD-12](prd-12-rerun.md) | **Complete** | Rerun visualization integration |
| [PRD-13](prd-13-release.md) | Draft | Release and publishing process |
| [PRD-14](prd-14-codegen.md) | In Progress | Geometric Algebra Code Generator |
| [PRD-15](prd-15-codegen-migration.md) | In Progress | Codegen Migration (Euclidean, Projective) |
| [PRD-16](prd-16-codegen-improvements.md) | Draft | Codegen Improvements: Blade Ordering, Versors, Simplification |
| [PRD-17](prd-17-codegen-products.md) | Draft | Codegen Product Completeness and Documentation |
| [PRD-18](prd-18-constraint-redesign.md) | In Progress | Constraint Redesign and Wrapper Types |
| [PRD-19](prd-19-symbolica-testing-strategy.md) | **Complete** | Symbolica Testing Strategy |
| [PRD-20](prd-20-projective-nalgebra-interop.md) | Draft | Projective-nalgebra Interoperability |
| [PRD-21](prd-21-2d-pga-dual-motors.md) | Draft | 2D PGA and Dual Motors |
| [PRD-22](prd-22-interior-products.md) | Draft | Interior Products and Product Naming Corrections |

## Implementation Order

1. ~~PRD-1: Foundation Layer~~ **Done**
2. ~~PRD-2: Core Multivector~~ **Done**
3. ~~PRD-3: Derived Products~~ **Done**
4. ~~PRD-4: Specialized Types~~ **Done**
5. PRD-5: PGA
6. PRD-6: CGA + Polish
7. PRD-7: Trait Abstractions (can be done in parallel with PRD-5/6)
8. ~~PRD-8: nalgebra Interop~~ **Done**
9. PRD-9: SIMD/simba (can be done anytime, independent of other PRDs)
10. ~~PRD-10: nalgebra Benchmarks~~ **Done**
11. PRD-11: Autodiff via Dual Numbers (independent, can be done anytime)

## Verification

Each PRD should pass:
- `cargo check` - no warnings
- `cargo test` - all proptest properties pass
- `cargo clippy` - no lints
- `cargo doc` - documentation builds
- `cargo bench` - benchmarks captured in `benches/reports/`
