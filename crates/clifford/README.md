<p align="center">
  <img src="https://raw.githubusercontent.com/DevonMorris/clifford/main/assets/clifford.png" alt="Clifford" width="240">
</p>

<h1 align="center">Clifford</h1>

<p align="center">
  <strong>Geometric Algebra for Rust</strong>
</p>

<p align="center">
  <a href="https://crates.io/crates/clifford"><img src="https://img.shields.io/crates/v/clifford.svg" alt="crates.io"></a>
  <a href="https://docs.rs/clifford"><img src="https://docs.rs/clifford/badge.svg" alt="docs.rs"></a>
  <a href="https://github.com/DevonMorris/clifford/actions"><img src="https://github.com/DevonMorris/clifford/workflows/CI/badge.svg" alt="CI"></a>
  <a href="https://github.com/DevonMorris/clifford/blob/main/LICENSE"><img src="https://img.shields.io/badge/license-MIT-blue.svg" alt="MIT License"></a>
</p>

<p align="center">
  <em>Rotations without gimbal lock. Rigid transforms in a single type. The geometry you wish you learned in school.</em>
</p>

---

**Clifford** is a Rust library for [Geometric Algebra](https://en.wikipedia.org/wiki/Geometric_algebra) (also known as Clifford Algebra), providing tools for 3D rotations, rigid body transformations, and computational geometry. If you're looking for an alternative to quaternions, rotation matrices, or homogeneous coordinates, Geometric Algebra offers a unified, intuitive approach.

## Why Geometric Algebra?

Geometric Algebra (GA) unifies vectors, complex numbers, quaternions, and more into a single elegant framework. Instead of juggling matrices, quaternions, and Euler angles, GA gives you:

- **Rotors**: Like quaternions but they generalize to any dimension
- **Motors**: Rotation + translation in one composable object (no more 4x4 matrices!)
- **Bivectors**: Oriented planes that represent rotations naturally
- **The geometric product**: One operation that subsumes dot and cross products

If you work with **robotics**, **computer graphics**, **physics simulations**, or **game development**, GA simplifies your math while eliminating edge cases like gimbal lock.

## Features

- **15 Specialized Algebras**: Complex, Dual, Quaternion, Dual Quaternion, 2D/3D Euclidean, 2D/3D Projective (PGA), 2D/3D Conformal (CGA), 2D/3D Minkowski, 2D Elliptic, Hyperbolic
- **Projective Geometric Algebra (PGA)** for rigid body transforms (points, lines, planes, motors)
- **Conformal Geometric Algebra (CGA)** for circles, spheres, and conformal transformations
- **Code Generation System** (`clifford-codegen`) for deriving optimized algebraic operations
- **Optimized 2D/3D types** with zero-cost abstractions
- **Interior Products and Projections** for computing orthogonal components
- **Seamless nalgebra integration** for existing codebases
- **Rerun visualization** (`clifford-viz`) for debugging and education
- **Generic over float types** (f32, f64)
- **Compile-time dimension checking**
- **Property-tested correctness** with proptest
- **WASM support** for browser-based applications

## Installation

```toml
[dependencies]
clifford = "0.3"
```

## Quick Start

### 3D Rotations with Rotors

Rotors are the GA equivalent of quaternions, but they arise naturally from the geometry:

```rust,ignore
use clifford::ops::Transform;
use clifford::specialized::euclidean::dim3::{Vector, Bivector, Rotor};
use std::f64::consts::FRAC_PI_2;

// Create a rotation of 90 degrees in the xy-plane (around the z-axis)
let plane = Bivector::unit_xy();
let rotor = Rotor::from_angle_plane(FRAC_PI_2, plane);

// Rotate a vector
let v = Vector::new(1.0, 0.0, 0.0);
let rotated = rotor.transform(&v);

// x-axis is now pointing along y-axis
assert!((rotated.x()).abs() < 1e-10);
assert!((rotated.y() - 1.0).abs() < 1e-10);
```

### Rigid Transforms with PGA Motors

Motors combine rotation and translation into a single object that composes beautifully:

```rust,ignore
use clifford::specialized::projective::dim3::{Motor, Point};

// Create a motor: translate by (1, 2, 3)
let motor = Motor::from_translation(1.0_f64, 2.0, 3.0);

// Transform a point at the origin
let p = Point::origin();
let transformed = motor.transform_point(&p);

// The origin moved to (1, 2, 3)
assert!((transformed.x() - 1.0).abs() < 1e-10);
assert!((transformed.y() - 2.0).abs() < 1e-10);
assert!((transformed.z() - 3.0).abs() < 1e-10);
```

Compose motors for complex transforms - rotation then translation, or vice versa:

Motors can also transform lines and planes, making them ideal for robotics and graphics.

### 2D Geometry

The same patterns work in 2D:

```rust,ignore
use clifford::ops::Transform;
use clifford::specialized::euclidean::dim2::{Vector, Rotor};
use std::f64::consts::FRAC_PI_2;

let v = Vector::new(1.0, 0.0);
let rotor = Rotor::from_angle(FRAC_PI_2);
let rotated = rotor.transform(&v);

assert!((rotated.x()).abs() < 1e-10);
assert!((rotated.y() - 1.0).abs() < 1e-10);
```

### Working with nalgebra

Already using nalgebra? Clifford integrates seamlessly:

```rust,ignore
use clifford::specialized::euclidean::dim3::{Vector, Rotor, Bivector};
use nalgebra::{Vector3, UnitQuaternion};

// From nalgebra
let na_vec = Vector3::new(1.0_f64, 2.0, 3.0);
let cliff_vec = Vector::from(na_vec);

// To nalgebra
let back_to_na: Vector3<f64> = cliff_vec.into();

// Quaternions too
let rotor = Rotor::from_angle_plane(0.5, Bivector::unit_xy());
let quat: UnitQuaternion<f64> = rotor.into();
```

## Visualization with Rerun

Debug your geometric algebra computations visually:

```bash
cargo run --example rerun_bivector --features rerun-0_28
```

This shows an animated bivector (parallelogram) with its Hodge dual vector, demonstrating how the wedge product works.

<p align="center">
  <em>See <code>examples/</code> for more visualization demos.</em>
</p>

## Module Structure

### Specialized Algebras

| Module | Signature | Description |
|--------|-----------|-------------|
| [`specialized::complex`](https://docs.rs/clifford/latest/clifford/specialized/complex/) | Cl(0,1) | Complex numbers (i^2 = -1) |
| [`specialized::dual`](https://docs.rs/clifford/latest/clifford/specialized/dual/) | Cl(0,0,1) | Dual numbers (epsilon^2 = 0) for automatic differentiation |
| [`specialized::quaternion`](https://docs.rs/clifford/latest/clifford/specialized/quaternion/) | Cl(0,2) | Quaternions for 3D rotations |
| [`specialized::dualquat`](https://docs.rs/clifford/latest/clifford/specialized/dualquat/) | - | Dual quaternions for rigid transforms |
| [`specialized::euclidean::dim2`](https://docs.rs/clifford/latest/clifford/specialized/euclidean/dim2/) | Cl(2,0) | 2D Euclidean: Vector, Bivector, Rotor |
| [`specialized::euclidean::dim3`](https://docs.rs/clifford/latest/clifford/specialized/euclidean/dim3/) | Cl(3,0) | 3D Euclidean: Vector, Bivector, Trivector, Rotor |
| [`specialized::projective::dim2`](https://docs.rs/clifford/latest/clifford/specialized/projective/dim2/) | Cl(2,0,1) | 2D PGA: Point, Line, Motor |
| [`specialized::projective::dim3`](https://docs.rs/clifford/latest/clifford/specialized/projective/dim3/) | Cl(3,0,1) | 3D PGA: Point, Line, Plane, Motor, Flector |
| [`specialized::conformal::dim2`](https://docs.rs/clifford/latest/clifford/specialized/conformal/dim2/) | Cl(3,1) | 2D CGA: Points, circles, conformal transforms |
| [`specialized::conformal::dim3`](https://docs.rs/clifford/latest/clifford/specialized/conformal/dim3/) | Cl(4,1) | 3D CGA: Points, spheres, circles, conformal transforms |
| [`specialized::minkowski::dim2`](https://docs.rs/clifford/latest/clifford/specialized/minkowski/dim2/) | Cl(1,1) | 2D Minkowski spacetime |
| [`specialized::minkowski::dim3`](https://docs.rs/clifford/latest/clifford/specialized/minkowski/dim3/) | Cl(1,2) | 3D Minkowski spacetime |
| [`specialized::elliptic::dim2`](https://docs.rs/clifford/latest/clifford/specialized/elliptic/dim2/) | Cl(2,0) + positive curvature | 2D elliptic/spherical geometry |
| [`specialized::hyperbolic`](https://docs.rs/clifford/latest/clifford/specialized/hyperbolic/) | - | Hyperbolic numbers (j^2 = +1) |

### Core Infrastructure

| Module | Description |
|--------|-------------|
| [`algebra`](https://docs.rs/clifford/latest/clifford/algebra/) | Generic multivector for any metric signature |
| [`signature`](https://docs.rs/clifford/latest/clifford/signature/) | Metric signatures (Euclidean, Minkowski, etc.) |
| [`ops`](https://docs.rs/clifford/latest/clifford/ops/) | Algebraic operations (Wedge, Antiwedge, Sandwich, Transform, etc.) |

## Cargo Features

| Feature | Description | Default |
|---------|-------------|---------|
| `serde` | Serialization/deserialization | Yes |
| `proptest-support` | Property-based testing strategies | Yes |
| `nalgebra-0_33` | nalgebra 0.33.x conversions | Yes |
| `nalgebra-0_32` | nalgebra 0.32.x conversions | No |
| `nalgebra-0_34` | nalgebra 0.34.x conversions | No |
| `rerun-0_28` | Rerun visualization integration | No |

**Note**: The nalgebra features are mutually exclusive. Enable only one.

## Performance

Clifford is designed for performance:

- Specialized types use fixed-size arrays (no heap allocation)
- Operations are inlined and optimized by LLVM
- Benchmarks show competitive performance with hand-written quaternion code

Run benchmarks yourself:

```bash
cargo bench
```

## Learning Resources

New to Geometric Algebra? These resources helped us:

- [Geometric Algebra Primer](http://www.jaapsuter.com/geometric-algebra.pdf) - Jaap Suter's excellent introduction
- [Let's Remove Quaternions from Every 3D Engine](https://marctenbosch.com/quaternions/) - Marc ten Bosch on why GA is better
- [Siggraph 2019 GA Course](https://www.youtube.com/watch?v=tX4H_ctggYo) - Video introduction
- [bivector.net](https://bivector.net/) - Community resources and tools

## Minimum Supported Rust Version

Clifford requires Rust 1.87.0 or later (edition 2024).

## License

MIT License

---

<p align="center">
  <strong>Ready to try a better way to do geometry?</strong><br>
  <code>cargo add clifford</code>
</p>
