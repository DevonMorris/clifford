# Clifford PRDs (Product Requirements Documents)

Implementation plan for the Clifford Geometric Algebra library.

## Overview

- **Dimensions**: Optimized 2D, 3D, 4D implementations
- **Signatures**: Full support (Euclidean, Minkowski, PGA, CGA, arbitrary)
- **Type System**: Code generation from TOML specifications
- **Scalars**: Generic over float types (f32, f64)

## PRDs

### Foundation & Core

| PRD | Status | Description |
|-----|--------|-------------|
| [PRD-0](prd-0-architecture.md) | Reference | Architecture overview |
| [PRD-1](prd-1-foundation.md) | **Complete** | Foundation: Float, Signature, Blade |
| [PRD-2](prd-2-multivector.md) | **Complete** | Core Multivector, geometric product |
| [PRD-3](prd-3-products.md) | **Complete** | Derived products, grade operations |
| [PRD-4](prd-4-specialized.md) | **Complete** | Optimized 2D/3D types with conversions |

### Geometric Algebras

| PRD | Status | Description |
|-----|--------|-------------|
| [PRD-5](prd-5-pga.md) | **Complete** | Projective GA (points, lines, planes, motors) |
| [PRD-6](prd-6-cga.md) | **Complete** | 3D Conformal GA Cl(4,1,0) |
| [PRD-50](prd-50-conformal2.md) | **Complete** | 2D Conformal GA Cl(3,1,0) |
| [PRD-32](prd-32-complex-numbers.md) | **Complete** | Complex numbers Cl(0,1,0) |
| [PRD-33](prd-33-dual-numbers.md) | **Complete** | Dual numbers Cl(0,0,1) |
| [PRD-36](prd-36-quaternions.md) | **Complete** | Quaternions Cl(0,2,0) |
| [PRD-37](prd-37-minkowski-plane.md) | **Complete** | Minkowski plane Cl(1,1,0) |
| [PRD-38](prd-38-dual-quaternions.md) | **Complete** | Dual quaternions Cl(0,2,1) |
| [PRD-39](prd-39-elliptic2.md) | **Complete** | 2D Elliptic geometry Cl(3,0,0) |
| [PRD-40](prd-40-hyperbolic-plane.md) | **Complete** | 2D Hyperbolic plane Cl(2,1,0) |
| [PRD-41](prd-41-minkowski3.md) | **Complete** | 3D Minkowski spacetime Cl(3,1,0) |

### Code Generation

| PRD | Status | Description |
|-----|--------|-------------|
| [PRD-14](prd-14-codegen.md) | **Complete** | Geometric Algebra Code Generator |
| [PRD-15](prd-15-codegen-migration.md) | **Complete** | Codegen Migration (Euclidean, Projective) |
| [PRD-16](prd-16-codegen-improvements.md) | **Complete** | Codegen Improvements: Blade Ordering, Versors |
| [PRD-17](prd-17-codegen-products.md) | **Complete** | Codegen Product Completeness |
| [PRD-23](prd-23-method-based-products.md) | **Complete** | Method-based Products |
| [PRD-24](prd-24-type-safe-products.md) | **Complete** | Type-safe Products |
| [PRD-25](prd-25-toml-cleanup.md) | **Complete** | TOML Cleanup |
| [PRD-26](prd-26-auto-versors.md) | **Complete** | Auto-identify Versors |
| [PRD-28](prd-28-build-script-regeneration.md) | **Complete** | Build Script Regeneration |
| [PRD-29](prd-29-semantic-field-names.md) | **Complete** | Semantic Field Names |

### Constraints & Norms

| PRD | Status | Description |
|-----|--------|-------------|
| [PRD-18](prd-18-constraint-redesign.md) | **Complete** | Constraint Redesign and Wrapper Types |
| [PRD-27](prd-27-auto-constraints.md) | **Complete** | Auto-infer Constraints |
| [PRD-34](prd-34-constraint-involution-fix.md) | **Complete** | Constraint Involution Fix |
| [PRD-35](prd-35-wrapper-trait-bounds.md) | **Complete** | Wrapper Trait Bounds |

### Products & Operations

| PRD | Status | Description |
|-----|--------|-------------|
| [PRD-22](prd-22-interior-products.md) | **Complete** | Interior Products |
| [PRD-30](prd-30-projections.md) | Draft | Projections and Antiprojections |

### Interoperability & Integration

| PRD | Status | Description |
|-----|--------|-------------|
| [PRD-7](prd-7-traits.md) | Draft | Trait abstractions for generic algorithms |
| [PRD-8](prd-8-nalgebra.md) | **Complete** | nalgebra interoperability |
| [PRD-9](prd-9-simd.md) | Draft | SIMD and simba integration |
| [PRD-10](prd-10-nalgebra-benchmarks.md) | **Complete** | nalgebra vs clifford benchmarks |
| [PRD-11](prd-11-autodiff.md) | Draft | Automatic differentiation via dual numbers |
| [PRD-12](prd-12-rerun.md) | **Complete** | Rerun visualization integration |
| [PRD-19](prd-19-symbolica-testing-strategy.md) | **Complete** | Symbolica Testing Strategy |
| [PRD-20](prd-20-projective-nalgebra-interop.md) | Draft | Projective-nalgebra Interoperability |
| [PRD-21](prd-21-2d-pga-dual-motors.md) | Draft | 2D PGA and Dual Motors |

### Release & Documentation

| PRD | Status | Description |
|-----|--------|-------------|
| [PRD-13](prd-13-release.md) | Draft | Release and publishing process |

## Implemented Algebras

The following algebras are implemented via code generation:

| Algebra | Signature | Module | Description |
|---------|-----------|--------|-------------|
| Euclidean 2D | Cl(2,0,0) | `euclidean::dim2` | 2D vectors, bivectors, rotors |
| Euclidean 3D | Cl(3,0,0) | `euclidean::dim3` | 3D vectors, bivectors, rotors |
| Projective 2D | Cl(2,0,1) | `projective::dim2` | 2D PGA: points, lines, motors |
| Projective 3D | Cl(3,0,1) | `projective::dim3` | 3D PGA: points, lines, planes, motors |
| Conformal 2D | Cl(3,1,0) | `conformal::dim2` | 2D CGA: points, circles, conformal transforms |
| Conformal 3D | Cl(4,1,0) | `conformal::dim3` | 3D CGA: points, spheres, circles |
| Complex | Cl(0,1,0) | `complex` | Complex numbers |
| Dual | Cl(0,0,1) | `dual` | Dual numbers (autodiff) |
| Quaternion | Cl(0,2,0) | `quaternion` | Quaternions |
| Hyperbolic | Cl(1,0,0) | `hyperbolic` | Hyperbolic (split-complex) numbers |
| Hyperbolic 2D | Cl(2,1,0) | `hyperbolic::dim2` | 2D hyperbolic plane |
| Elliptic 2D | Cl(3,0,0) | `elliptic::dim2` | 2D elliptic geometry |
| Minkowski 2D | Cl(1,1,0) | `minkowski::dim2` | Minkowski plane |
| Minkowski 3D | Cl(1,2,0) | `minkowski::dim3` | 3D spacetime (1 time + 2 space) |
| Dual Quaternion | Cl(0,2,1) | `dualquat` | Dual quaternions |

## Verification

Each PRD should pass:
- `cargo check` - no warnings
- `cargo test` - all proptest properties pass
- `cargo clippy` - no lints
- `cargo doc` - documentation builds
