<p align="center">
  <img src="assets/clifford.png" alt="Clifford" width="240">
</p>

<h1 align="center">Clifford</h1>

<p align="center">
  <strong>Geometric Algebra for Rust</strong>
</p>

<p align="center">
  <a href="https://crates.io/crates/clifford"><img src="https://img.shields.io/crates/v/clifford.svg" alt="crates.io"></a>
  <a href="https://docs.rs/clifford"><img src="https://docs.rs/clifford/badge.svg" alt="docs.rs"></a>
  <a href="https://github.com/DevonMorris/clifford/actions"><img src="https://github.com/DevonMorris/clifford/workflows/CI/badge.svg" alt="CI"></a>
</p>

---

A Rust library for Geometric Algebra (Clifford Algebra) with a focus on correctness, performance, and education.

Geometric Algebra (GA) provides a unified mathematical language for geometry. It extends traditional vector algebra with the geometric product, which combines the dot product and wedge product into a single operation.

## Features

- **Generic multivector types** over any metric signature
- **Geometric, inner, and outer products** with grade operations
- **Euclidean, Projective (PGA), and Conformal (CGA) algebras**
- **Optimized 2D/3D specialized types** for common use cases
- **Strong nalgebra interop** - seamless conversions to/from nalgebra vectors, rotations, isometries, and quaternions
- **Compile-time dimension checking**

## Installation

```toml
[dependencies]
clifford = "0.1"
```

## Quick Start

### Specialized 2D/3D Types

For common 2D and 3D Euclidean geometry, use the optimized specialized types:

```rust
use clifford::specialized::euclidean::{dim2, dim3};
use std::f64::consts::FRAC_PI_2;

// 2D: Rotate a vector by 90 degrees
let v = dim2::Vector::new(1.0, 0.0);
let rotor = dim2::Rotor::from_angle(FRAC_PI_2);
let rotated = rotor.rotate(v);
assert!((rotated.x).abs() < 1e-10);
assert!((rotated.y - 1.0).abs() < 1e-10);

// 3D: Rotate around an axis
let v = dim3::Vector::new(1.0, 0.0, 0.0);
let plane = dim3::Bivector::unit_xy(); // rotation in xy-plane (around z-axis)
let rotor = dim3::Rotor::from_angle_plane(FRAC_PI_2, plane);
let rotated = rotor.rotate(v);
assert!((rotated.x).abs() < 1e-10);
assert!((rotated.y - 1.0).abs() < 1e-10);
assert!((rotated.z).abs() < 1e-10);
```

### Generic Multivector API

The generic API works with any metric signature:

```rust
use clifford::algebra::Multivector;
use clifford::signature::Euclidean3;

// Create two vectors
let a: Multivector<f64, Euclidean3> = Multivector::vector(&[1.0, 2.0, 0.0]);
let b: Multivector<f64, Euclidean3> = Multivector::vector(&[0.0, 1.0, 0.0]);

// Geometric product combines dot and wedge products
let ab = &a * &b;

// Dot product is the scalar part: a·b = 1*0 + 2*1 + 0*0 = 2
assert!((ab.scalar_part() - 2.0).abs() < 1e-10);
```

## Module Structure

- [`algebra`](https://docs.rs/clifford/latest/clifford/algebra/) - Core algebraic types (`Multivector`)
- [`scalar`](https://docs.rs/clifford/latest/clifford/scalar/) - Floating-point scalar type abstraction
- [`signature`](https://docs.rs/clifford/latest/clifford/signature/) - Metric signatures defining the algebra
- [`basis`](https://docs.rs/clifford/latest/clifford/basis/) - Basis blade representation and utilities
- [`specialized`](https://docs.rs/clifford/latest/clifford/specialized/) - Optimized types for 2D/3D Euclidean and Projective geometry

## Cargo Features

- `serde` - Enable serialization/deserialization via serde
- `proptest-support` - Enable proptest strategies for property-based testing
- `nalgebra-0_32` - Enable conversions to/from nalgebra 0.32.x types
- `nalgebra-0_33` - Enable conversions to/from nalgebra 0.33.x types (default)
- `nalgebra-0_34` - Enable conversions to/from nalgebra 0.34.x types

Note: The nalgebra features are mutually exclusive. Enable only one.

## License

MIT
