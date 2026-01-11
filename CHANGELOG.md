# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

## [0.1.0] - 2026-01-11

### Added

#### Core Algebra
- Generic `Multivector<T, S>` type supporting any metric signature (#13)
- Geometric product implementation with compile-time dimension checking (#13)
- Inner product (left contraction) and outer product (wedge) (#17)
- Grade selection, extraction, and involution operations (#17)
- Dual and reverse operations (#17)
- `Float` trait for generic scalar types (f32, f64) (#11)
- `Signature` trait for metric signatures with compile-time blade count (#11)
- Built-in signatures: `Euclidean2`, `Euclidean3`, `Euclidean4` (#11)
- `Blade` type for basis blade representation (#11)

#### Specialized 2D/3D Types
- `euclidean::dim2`: `Vector`, `Bivector`, `Rotor`, `Even`, `Multivector2` (#20)
- `euclidean::dim3`: `Vector`, `Bivector`, `Trivector`, `Rotor`, `Even`, `Multivector3` (#20)
- Rotor construction from angle (2D) and angle-plane (3D) (#20)
- Vector rotation via sandwich product (#20)
- Bidirectional conversions between specialized and generic types (#24, #26)
- Private fields with public accessors for all specialized types (#37)

#### Projective Geometric Algebra (PGA)
- `projective::dim2`: `Point`, `Line`, `Motor` for 2D rigid transforms (#36)
- `projective::dim3`: `Point`, `Line`, `Plane`, `Motor`, `Flector` for 3D rigid transforms (#36)
- Motor construction from rotations, translations, and axis-angle (#36)
- Motor composition for chaining transforms (#36)
- Point, line, and plane transformations via motors (#36)
- Meet and join operations for geometric intersections (#36)
- Distance and angle calculations (#36)

#### Integrations
- **nalgebra**: Bidirectional conversions for vectors, rotations, quaternions, isometries (#29)
  - Feature-gated: `nalgebra-0_32`, `nalgebra-0_33` (default), `nalgebra-0_34`
- **Rerun**: Visualization support for geometric types (#39)
  - `AsPosition` and `AsArrow` wrapper types for disambiguation
  - Conversions to Rerun's `Vec3D`, `Position3D`, `Transform3D`, etc.
  - Feature-gated: `rerun-0_28`
- **serde**: Serialization/deserialization support (enabled by default)
- **proptest**: Property-based testing strategies via `Arbitrary` trait (#23)
  - Wrapper types: `NonZeroVector`, `UnitVector`, `UnitRotor`, etc.
  - Feature-gated: `proptest-support` (enabled by default)

#### Examples
- `rerun_visualization`: Basic Rerun integration demo (#39)
- `rerun_animation`: Animated rotors and motors (#39)
- `rerun_bivector`: Bivector and Hodge dual visualization (#39)

#### Benchmarks
- Generic multivector benchmarks: geometric product, grade operations (#18)
- Specialized type benchmarks: rotor rotation, vector operations (#20)
- nalgebra comparison benchmarks: rotation, isometry transforms (#34)
- PGA benchmarks: motor composition, point/line transforms (#36)

#### Documentation
- Comprehensive rustdoc for all public APIs
- Mathematical explanations in doc comments
- PRD (Product Requirements Documents) for architecture planning
- Claude Code agents for development workflows

### Infrastructure
- GitHub Actions CI: check, test, clippy, fmt, doc, cargo-deny (#6)
- Greptile AI code review integration (#6)
- Branch protection requiring CI and review
- crates.io publish workflow triggered by version tags (#7)
- License compliance via cargo-deny (MIT, Apache-2.0, BSD, etc.) (#32)
