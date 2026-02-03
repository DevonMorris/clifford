# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

## [0.2.0] - 2026-02-02

### Added

#### Code Generation System
- **clifford-codegen**: Full-featured code generation tool for geometric algebras (#51)
  - Define algebras declaratively via TOML specification files
  - Automatic generation of types, traits, products, and conversions
  - Groebner basis constraint simplification for wrapper types (#87)
  - Semantic field names based on geometric meaning, not blade indices (#70)
  - Auto-regeneration on every build via cargo build script (#69)
  - Project and Antiproject traits for geometric projections (#73)
  - Interior products with RGA duality laws (#59)

#### New Algebras (15 total)
All algebras defined declaratively via codegen TOML specs:
- **Complex numbers** Cl(0,1,0): Real and imaginary components (#75)
- **Dual numbers** Cl(0,0,1): Real and infinitesimal components (#76)
- **Hyperbolic numbers** Cl(1,0,0): Split-complex numbers (#74)
- **Quaternions** Cl(0,2,0): 3D rotation algebra (#79)
- **Dual quaternions** Cl(0,2,1): Rigid body transformations (#81)
- **2D Minkowski** Cl(1,1,0): 1+1 spacetime (#80)
- **3D Minkowski** Cl(3,1,0): 3+1 spacetime (#84)
- **2D Hyperbolic** Cl(2,1,0): Hyperbolic plane geometry (#83)
- **2D Elliptic** Cl(3,0,0): Elliptic projective geometry (#82)
- **2D Conformal GA** Cl(3,1,0): Circles, inversions, Mobius transforms (#95)
- **3D Conformal GA** Cl(4,1,0): Spheres, circles, conformal transforms (#86)

#### Visualization Crate (clifford-viz)
- **New crate**: Interactive visualization demos at clifford-rs.dev (#93)
- **2D rendering**: egui/eframe integration with dark/light theme support
- **3D rendering**: three-d backend for native and WASM (#124, #127)
- **WASM deployment**: Live demos at https://clifford-rs.dev (#105)
- **Mobile support**: Responsive design for all screen sizes (#111, #112)

##### Visualization Demos
- **Euclidean 2D**: Rotor rotation with angle interpolation (#94)
- **Euclidean 3D**: Quaternion-based 3D rotations (#127)
- **Projective 2D**: Interactive point/line/motor manipulations (#99)
- **Projective 3D**: 3D rigid transformations and robot arm demo (#101, #128)
- **Conformal 2D**: Circle inversion, Mobius transforms, circle from 3 points (#104, #108, #119)
- **Conformal 3D**: 3D conformal transformations (#130)
- **Complex/Dual**: Number system visualizations (#116)

#### Infrastructure
- **Workspace restructure**: Main crate moved to crates/clifford (#91)
- **Build performance**: Symbolica compile time reduced from 313s to 43s (#92)
- **CI improvements**: Cache generated code to avoid rebuild conflicts (#97)
- **Visual testing**: Property-based visual invariant testing (#102)

### Changed
- Migrated euclidean and projective types to use codegen (#53, #54, #56)
- Extension methods now prefer generated traits over manual implementations (#71, #72)
- Product traits use method-based API for ergonomics (#62)
- Versor `inverse()` is now provided by the `VersorInverse` trait (alongside `try_inverse()`), removing bespoke implementations from extensions

### Fixed
- Point::distance_squared now correctly handles non-unitized homogeneous points
- Clippy compliance for all targets including tests
- Orthonormal basis conventions for 2D CGA (#115)
- Rotation direction in euclidean algebras (#68)
- Infinite line clipping in visualization (#103)

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
