# PRD-9: SIMD and simba Integration

**Status**: Draft
**Goal**: High-performance SIMD operations and compatibility with simba's scalar abstraction

## Motivation

Geometric algebra operations (geometric product, rotations, etc.) are inherently parallel and benefit significantly from SIMD vectorization. This PRD covers two related features:

1. **simba integration**: Compatibility with [simba](https://crates.io/crates/simba), nalgebra's scalar abstraction crate, enabling use of SIMD scalar types
2. **Direct SIMD optimization**: Vectorized implementations of core operations using portable SIMD

## Background

### simba Crate

[simba](https://lib.rs/crates/simba) provides SIMD algebra abstractions for Rust. Key features:
- Traits like `SimdValue`, `SimdRealField`, `SimdComplexField`
- SIMD type wrappers (`WideF32x4`, `WideF64x2`) via the `wide` feature
- Works on stable Rust (with `wide` feature) or nightly (with `portable_simd`)

### SIMD Options in Rust (2025)

| Approach | Stability | Platforms | Notes |
|----------|-----------|-----------|-------|
| `std::simd` | Nightly | All LLVM targets | Portable, pairs with `multiversion` |
| `wide` crate | Stable | x86, NEON, WASM | Mature, used by simba |
| Intrinsics | Stable (1.87+) | Platform-specific | Most intrinsics now safe |

## Deliverables

### 1. Feature Flags

```toml
[features]
default = []

# simba integration
simba = ["dep:simba"]

# Direct SIMD (pick backend)
simd-wide = ["dep:wide"]           # Stable, recommended
simd-portable = []                  # Nightly, uses std::simd

[dependencies]
simba = { version = "0.9", optional = true }
wide = { version = "0.7", optional = true }
```

### 2. simba Trait Compatibility

Make our `Float` trait interoperable with simba:

```rust
// src/scalar/simba.rs
#[cfg(feature = "simba")]
use simba::scalar::{RealField, SubsetOf};

/// When simba is enabled, Float is implemented for any type
/// that satisfies simba's RealField requirements.
#[cfg(feature = "simba")]
impl<T> Float for T
where
    T: RealField + Copy + Default + 'static,
{
    // Delegate to simba's methods
}
```

Alternatively, provide bridging traits:

```rust
/// Bridge trait for using simba types with clifford.
#[cfg(feature = "simba")]
pub trait SimbaFloat: simba::scalar::RealField + Float {}

#[cfg(feature = "simba")]
impl<T: simba::scalar::RealField + Float> SimbaFloat for T {}
```

### 3. SIMD Multivector Operations

Vectorize the geometric product and other core operations:

```rust
// src/algebra/simd.rs

/// SIMD-optimized geometric product for Euclidean3.
///
/// Processes 4 multivector pairs simultaneously using f32x4.
#[cfg(feature = "simd-wide")]
pub fn geometric_product_batch_e3(
    a: &[Multivector<f32, Euclidean3>; 4],
    b: &[Multivector<f32, Euclidean3>; 4],
) -> [Multivector<f32, Euclidean3>; 4] {
    // SIMD implementation
}

/// SIMD-optimized rotation for batches of vectors.
#[cfg(feature = "simd-wide")]
impl Rotor<f32> {
    pub fn rotate_batch(&self, vectors: &[Vector<f32>; 4]) -> [Vector<f32>; 4] {
        // SIMD sandwich product
    }
}
```

### 4. SIMD Specialized Types

SIMD-native specialized types for batch processing:

```rust
// src/specialized/euclidean/dim3/simd.rs

/// Four 3D vectors stored in SIMD layout (SoA - Structure of Arrays).
#[cfg(feature = "simd-wide")]
#[repr(C)]
pub struct Vector4 {
    pub x: wide::f32x4,
    pub y: wide::f32x4,
    pub z: wide::f32x4,
}

impl Vector4 {
    /// Dot product of 4 vector pairs simultaneously.
    pub fn dot(&self, other: &Self) -> wide::f32x4 {
        self.x * other.x + self.y * other.y + self.z * other.z
    }

    /// Cross product of 4 vector pairs simultaneously.
    pub fn cross(&self, other: &Self) -> Self {
        Self {
            x: self.y * other.z - self.z * other.y,
            y: self.z * other.x - self.x * other.z,
            z: self.x * other.y - self.y * other.x,
        }
    }

    /// Normalize all 4 vectors.
    pub fn normalized(&self) -> Self {
        let norm = (self.x * self.x + self.y * self.y + self.z * self.z).sqrt();
        Self {
            x: self.x / norm,
            y: self.y / norm,
            z: self.z / norm,
        }
    }
}

/// Four 3D rotors in SIMD layout.
#[cfg(feature = "simd-wide")]
#[repr(C)]
pub struct Rotor4 {
    pub s: wide::f32x4,
    pub xy: wide::f32x4,
    pub xz: wide::f32x4,
    pub yz: wide::f32x4,
}

impl Rotor4 {
    /// Apply 4 rotations to 4 vectors simultaneously.
    pub fn rotate(&self, v: &Vector4) -> Vector4 {
        // SIMD sandwich product: R v R̃
    }

    /// Compose 4 rotor pairs simultaneously.
    pub fn compose(&self, other: &Self) -> Self {
        // SIMD geometric product of even elements
    }

    /// Slerp 4 rotor pairs simultaneously.
    pub fn slerp(&self, other: &Self, t: wide::f32x4) -> Self {
        // SIMD slerp
    }
}
```

### 5. Portable SIMD (Nightly)

For users on nightly, provide `std::simd` implementations:

```rust
// src/specialized/euclidean/dim3/portable_simd.rs
#![cfg(feature = "simd-portable")]

use std::simd::{f32x4, f64x2, SimdFloat};

/// Four 3D vectors using portable SIMD.
#[repr(C)]
pub struct Vector4Portable {
    pub x: f32x4,
    pub y: f32x4,
    pub z: f32x4,
}

// Similar implementations as Vector4, but using std::simd
```

## Module Structure

```
src/
  scalar/
    mod.rs
    float.rs
    simba.rs            # simba trait bridges (feature-gated)
  algebra/
    simd.rs             # SIMD batch operations (feature-gated)
  specialized/
    euclidean/
      dim3/
        simd.rs         # SIMD Vector4, Rotor4, etc.
        portable_simd.rs # std::simd versions (nightly)
```

## Performance Targets

Based on typical SIMD speedups:

| Operation | Scalar | SIMD (4-wide) | Expected Speedup |
|-----------|--------|---------------|------------------|
| Vector dot | ~3 ops | ~3 ops (4 pairs) | ~4x |
| Vector cross | ~9 ops | ~9 ops (4 pairs) | ~4x |
| Rotor rotate | ~28 ops | ~28 ops (4 vectors) | ~3-4x |
| Geometric product (3D) | ~64 ops | ~64 ops (4 pairs) | ~3-4x |

Speedups may be less than 4x due to:
- Memory bandwidth limitations
- Horizontal operations (reductions)
- Data layout conversion overhead

## Testing

### Property-Based Tests

```rust
#[cfg(all(test, feature = "simd-wide"))]
mod simd_tests {
    use super::*;
    use crate::test_utils::ABS_DIFF_EQ_EPS_F32;
    use proptest::prelude::*;

    proptest! {
        #[test]
        fn simd_dot_matches_scalar(
            v1 in any::<[dim3::Vector<f32>; 4]>(),
            v2 in any::<[dim3::Vector<f32>; 4]>(),
        ) {
            let simd_v1 = Vector4::from(v1);
            let simd_v2 = Vector4::from(v2);
            let simd_result: [f32; 4] = simd_v1.dot(&simd_v2).into();

            for i in 0..4 {
                let scalar_result = v1[i].dot(v2[i]);
                prop_assert!(abs_diff_eq!(
                    simd_result[i], scalar_result,
                    epsilon = ABS_DIFF_EQ_EPS_F32
                ));
            }
        }

        #[test]
        fn simd_rotate_matches_scalar(
            r in any::<[dim3::arbitrary::UnitRotor<f32>; 4]>(),
            v in any::<[dim3::Vector<f32>; 4]>(),
        ) {
            let simd_r = Rotor4::from(r.map(|x| *x));
            let simd_v = Vector4::from(v);
            let simd_result = simd_r.rotate(&simd_v);
            let result_arr: [dim3::Vector<f32>; 4] = simd_result.into();

            for i in 0..4 {
                let scalar_result = r[i].rotate(v[i]);
                prop_assert!(abs_diff_eq!(
                    result_arr[i], scalar_result,
                    epsilon = ABS_DIFF_EQ_EPS_F32
                ));
            }
        }
    }
}
```

Note: `ABS_DIFF_EQ_EPS_F32` would be a new constant for f32 precision (larger than the f64 `ABS_DIFF_EQ_EPS` due to reduced precision).

### Benchmarks

Add SIMD benchmarks to `benches/specialized.rs`:

```rust
#[cfg(feature = "simd-wide")]
fn bench_simd_rotate(c: &mut Criterion) {
    let rotors: [Rotor<f32>; 4] = [...];
    let vectors: [Vector<f32>; 4] = [...];
    let simd_r = Rotor4::from(rotors);
    let simd_v = Vector4::from(vectors);

    c.bench_function("euclidean/dim3/simd/rotor_rotate_x4", |b| {
        b.iter(|| black_box(&simd_r).rotate(black_box(&simd_v)))
    });
}
```

## simba Integration Details

### Trait Relationships

```
simba::scalar::RealField
    ├── simba::scalar::ComplexField
    │   └── simba::scalar::Field
    │       └── simba::scalar::Ring
    │           └── simba::scalar::ClosedAdd + ClosedMul + ...
    └── (inherits from all the above)

clifford::scalar::Float
    ├── Copy, Clone, Default, PartialEq, PartialOrd
    ├── Add, Sub, Mul, Div, Neg (+ Assign variants)
    ├── Sum, Product
    ├── approx traits (AbsDiffEq, RelativeEq, UlpsEq)
    └── Custom methods (abs, sqrt, sin, cos, atan2, acos, from_*)
```

### Integration Strategy

Option A: **Blanket impl Float for RealField types**
```rust
#[cfg(feature = "simba")]
impl<T> Float for T where T: RealField + Copy + Default + ... { }
```
- Pro: Automatic support for all simba-compatible types
- Con: May conflict with existing f32/f64 impls, orphan rules

Option B: **Wrapper types**
```rust
#[cfg(feature = "simba")]
pub struct SimbaScalar<T>(pub T);

impl<T: RealField> Float for SimbaScalar<T> { }
```
- Pro: No conflicts, explicit opt-in
- Con: Extra wrapper, less ergonomic

Option C: **Separate trait hierarchy** (Recommended)
```rust
/// Marker trait for simba-compatible Float types.
#[cfg(feature = "simba")]
pub trait SimbaCompatible: Float + simba::scalar::RealField {}

#[cfg(feature = "simba")]
impl SimbaCompatible for f32 {}
#[cfg(feature = "simba")]
impl SimbaCompatible for f64 {}
```
- Pro: Clear separation, no orphan issues
- Con: Users must import additional trait

## Documentation

Each SIMD type should document:
1. Memory layout (SoA vs AoS)
2. Lane semantics (which scalar maps to which lane)
3. Performance characteristics
4. Conversion to/from scalar arrays

Example:
```rust
/// Four 3D vectors in Structure-of-Arrays (SoA) layout.
///
/// # Memory Layout
///
/// Rather than storing 4 vectors as `[Vector; 4]` (Array of Structures),
/// this type stores components contiguously: all x values together,
/// then all y values, then all z values. This enables efficient SIMD
/// operations across corresponding components.
///
/// ```text
/// AoS: [x0,y0,z0], [x1,y1,z1], [x2,y2,z2], [x3,y3,z3]
/// SoA: [x0,x1,x2,x3], [y0,y1,y2,y3], [z0,z1,z2,z3]
/// ```
///
/// # Example
///
/// ```
/// use clifford::specialized::euclidean::dim3::{Vector, simd::Vector4};
///
/// let vectors = [
///     Vector::new(1.0, 0.0, 0.0),
///     Vector::new(0.0, 1.0, 0.0),
///     Vector::new(0.0, 0.0, 1.0),
///     Vector::new(1.0, 1.0, 1.0),
/// ];
/// let simd_vectors = Vector4::from(vectors);
/// let dots = simd_vectors.dot(&simd_vectors); // [1.0, 1.0, 1.0, 3.0]
/// ```
pub struct Vector4 { ... }
```

## Verification

- [ ] `cargo check --features simba` passes
- [ ] `cargo check --features simd-wide` passes
- [ ] `cargo test --features simba` passes
- [ ] `cargo test --features simd-wide` passes
- [ ] `cargo bench --features simd-wide` shows expected speedups
- [ ] SIMD results match scalar results (property tests)
- [ ] All SIMD types have comprehensive rustdoc

## Future Considerations

- **f64x2 / f64x4**: Double-precision SIMD for higher accuracy
- **AVX-512**: 8-wide or 16-wide operations on supported hardware
- **GPU offload**: For very large batches, consider compute shaders
- **Auto-vectorization hints**: Use `#[repr(simd)]` and alignment attributes
- **Runtime dispatch**: Use `multiversion` crate for runtime CPU feature detection

## References

- [simba on crates.io](https://crates.io/crates/simba)
- [wide crate](https://crates.io/crates/wide)
- [The state of SIMD in Rust in 2025](https://shnatsel.medium.com/the-state-of-simd-in-rust-in-2025-32c263e5f53d)
- [std::simd documentation](https://doc.rust-lang.org/std/simd/index.html)
