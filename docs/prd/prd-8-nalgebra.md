# PRD-8: nalgebra Interoperability

**Status**: Draft
**Goal**: Seamless conversions between clifford types and nalgebra vectors/matrices

## Motivation

[nalgebra](https://crates.io/crates/nalgebra) is the most widely-used linear algebra library in the Rust ecosystem (~2.8M downloads/month). Many users will want to:

1. Use clifford for geometric algebra operations while keeping nalgebra for general linear algebra
2. Interface with existing codebases that use nalgebra
3. Leverage nalgebra's ecosystem (nalgebra-glm, nalgebra-sparse, etc.)

## Version Support Strategy

Support multiple nalgebra versions via feature flags to accommodate users who cannot upgrade:

| Feature Flag | nalgebra Version | Notes |
|--------------|------------------|-------|
| `nalgebra-0-33` | 0.33.x | LTS - many projects still use this |
| `nalgebra-0-34` | 0.34.x | Current stable (default if any nalgebra feature enabled) |

**Rationale**:
- Supporting 2 minor versions covers ~18 months of releases
- Users on older versions often have dependency constraints from other crates
- Each minor version may have breaking API changes

## Deliverables

### 1. Feature Flags in Cargo.toml

```toml
[features]
default = []

# nalgebra interop (pick one)
nalgebra-0-33 = ["dep:nalgebra-0-33"]
nalgebra-0-34 = ["dep:nalgebra-0-34"]

[dependencies]
nalgebra-0-33 = { package = "nalgebra", version = "0.33", optional = true }
nalgebra-0-34 = { package = "nalgebra", version = "0.34", optional = true }
```

### 2. Conversion Traits

Define internal traits to abstract over nalgebra versions:

```rust
// src/interop/mod.rs
#[cfg(any(feature = "nalgebra-0-33", feature = "nalgebra-0-34"))]
pub mod nalgebra;

// src/interop/nalgebra/mod.rs
//! Conversions between clifford and nalgebra types.
//!
//! Enable with feature `nalgebra-0-33` or `nalgebra-0-34`.
```

### 3. Vector Conversions

```rust
// 2D Vector conversions
impl<T: Float + na::Scalar> From<na::Vector2<T>> for dim2::Vector<T>;
impl<T: Float + na::Scalar> From<dim2::Vector<T>> for na::Vector2<T>;

// 3D Vector conversions
impl<T: Float + na::Scalar> From<na::Vector3<T>> for dim3::Vector<T>;
impl<T: Float + na::Scalar> From<dim3::Vector<T>> for na::Vector3<T>;

// Generic multivector vector-part extraction
impl<T: Float + na::Scalar, S: Signature> Multivector<T, S> {
    /// Extract vector components as nalgebra vector.
    /// Returns None if dimension doesn't match or non-vector grades present.
    pub fn try_to_na_vector<const D: usize>(&self) -> Option<na::SVector<T, D>>;
}
```

### 4. Matrix Conversions (for Rotors)

Rotors can be converted to/from rotation matrices:

```rust
// 2D: Rotor <-> Rotation2 / Matrix2
impl<T: Float + na::RealField> From<dim2::Rotor<T>> for na::Rotation2<T>;
impl<T: Float + na::RealField> From<na::Rotation2<T>> for dim2::Rotor<T>;
impl<T: Float + na::RealField> From<dim2::Rotor<T>> for na::Matrix2<T>;

// 3D: Rotor <-> Rotation3 / Matrix3 / UnitQuaternion
impl<T: Float + na::RealField> From<dim3::Rotor<T>> for na::Rotation3<T>;
impl<T: Float + na::RealField> From<na::Rotation3<T>> for dim3::Rotor<T>;
impl<T: Float + na::RealField> From<dim3::Rotor<T>> for na::Matrix3<T>;
impl<T: Float + na::RealField> From<dim3::Rotor<T>> for na::UnitQuaternion<T>;
impl<T: Float + na::RealField> From<na::UnitQuaternion<T>> for dim3::Rotor<T>;
```

### 5. Bivector ↔ Antisymmetric Matrix

3D bivectors correspond to antisymmetric matrices (useful for angular velocity):

```rust
impl<T: Float + na::Scalar> From<dim3::Bivector<T>> for na::Matrix3<T> {
    /// Convert bivector to antisymmetric matrix.
    ///
    /// The bivector (xy, xz, yz) maps to:
    /// ```text
    /// [  0  -xy -xz ]
    /// [ xy   0  -yz ]
    /// [ xz  yz   0  ]
    /// ```
    fn from(b: Bivector<T>) -> Self;
}

impl<T: Float + na::Scalar> TryFrom<na::Matrix3<T>> for dim3::Bivector<T> {
    type Error = ConversionError;

    /// Extract bivector from antisymmetric matrix.
    /// Fails if matrix is not antisymmetric within tolerance.
    fn try_from(m: Matrix3<T>) -> Result<Self, Self::Error>;
}
```

### 6. Point Conversions (for future PGA/CGA)

```rust
// When PGA is implemented:
impl<T: Float + na::Scalar> From<na::Point2<T>> for pga2d::Point<T>;
impl<T: Float + na::Scalar> From<pga2d::Point<T>> for na::Point2<T>;

impl<T: Float + na::Scalar> From<na::Point3<T>> for pga3d::Point<T>;
impl<T: Float + na::Scalar> From<pga3d::Point<T>> for na::Point3<T>;
```

## Module Structure

```
src/
  interop/
    mod.rs              # Feature-gated module exports
    nalgebra/
      mod.rs            # Common traits and re-exports
      v0_33.rs          # 0.33-specific implementations (if needed)
      v0_34.rs          # 0.34-specific implementations (if needed)
      vector.rs         # Vector conversions
      matrix.rs         # Matrix/rotation conversions
      specialized.rs    # dim2/dim3 type conversions
```

## Testing

### Property-Based Tests

```rust
proptest! {
    #[test]
    fn vec3_nalgebra_roundtrip(v in any::<dim3::Vector<f64>>()) {
        let na_v: na::Vector3<f64> = v.into();
        let back: dim3::Vector<f64> = na_v.into();
        prop_assert!(abs_diff_eq!(v, back, epsilon = ABS_DIFF_EQ_EPS));
    }

    #[test]
    fn rotor3_to_quaternion_preserves_rotation(
        r in any::<UnitRotor<f64>>(),
        v in any::<dim3::Vector<f64>>(),
    ) {
        // Rotate with rotor
        let rotated_ga = r.rotate(v);

        // Rotate with quaternion
        let q: na::UnitQuaternion<f64> = (*r).into();
        let na_v: na::Vector3<f64> = v.into();
        let rotated_na = q * na_v;

        let rotated_back: dim3::Vector<f64> = rotated_na.into();
        prop_assert!(abs_diff_eq!(rotated_ga, rotated_back, epsilon = ABS_DIFF_EQ_EPS));
    }

    #[test]
    fn rotor3_quaternion_roundtrip(r in any::<UnitRotor<f64>>()) {
        let q: na::UnitQuaternion<f64> = (*r).into();
        let back: dim3::Rotor<f64> = q.into();
        // Note: rotors and quaternions have double-cover, so we check rotation equivalence
        let v = dim3::Vector::new(1.0, 2.0, 3.0);
        prop_assert!(abs_diff_eq!(r.rotate(v), back.rotate(v), epsilon = ABS_DIFF_EQ_EPS));
    }

    #[test]
    fn bivector_antisymmetric_matrix_roundtrip(b in any::<dim3::Bivector<f64>>()) {
        let m: na::Matrix3<f64> = b.into();
        // Verify antisymmetric
        prop_assert!(abs_diff_eq!(m, -m.transpose(), epsilon = ABS_DIFF_EQ_EPS));
        let back = dim3::Bivector::try_from(m).unwrap();
        prop_assert!(abs_diff_eq!(b, back, epsilon = ABS_DIFF_EQ_EPS));
    }
}
```

### Version-Specific Tests

```rust
#[cfg(feature = "nalgebra-0-33")]
mod tests_v0_33 {
    // Any 0.33-specific behavior tests
}

#[cfg(feature = "nalgebra-0-34")]
mod tests_v0_34 {
    // Any 0.34-specific behavior tests
}
```

## Documentation

Each conversion should document:
1. The mathematical correspondence
2. Any sign conventions or ordering differences
3. Precision/normalization considerations

Example:
```rust
/// Convert a 3D rotor to a nalgebra unit quaternion.
///
/// # Mathematical Correspondence
///
/// A rotor `R = s + xy·e₁₂ + xz·e₁₃ + yz·e₂₃` maps to quaternion
/// `q = s + xz·i + yz·j + xy·k`.
///
/// Note the component reordering: GA bivector components use lexicographic
/// basis ordering (e₁₂, e₁₃, e₂₃), while quaternions use (i, j, k) which
/// corresponds to (e₂₃, e₃₁, e₁₂) in GA.
///
/// # Normalization
///
/// This conversion preserves the rotor's normalization. A unit rotor
/// produces a unit quaternion.
impl<T: Float + na::RealField> From<Rotor<T>> for na::UnitQuaternion<T> {
    fn from(r: Rotor<T>) -> Self { ... }
}
```

## Error Handling

```rust
/// Error when converting from nalgebra types to clifford types.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum NalgebraConversionError {
    /// Matrix is not antisymmetric (for bivector conversion).
    NotAntisymmetric,
    /// Matrix is not a valid rotation matrix.
    NotRotation,
    /// Quaternion is not normalized.
    NotUnit,
}
```

## CI Considerations

Test matrix should include both supported nalgebra versions:

```yaml
# In .github/workflows/ci.yml
jobs:
  test-nalgebra:
    strategy:
      matrix:
        nalgebra: [nalgebra-0-33, nalgebra-0-34]
    steps:
      - run: cargo test --features ${{ matrix.nalgebra }}
```

## Verification

- [ ] `cargo check --features nalgebra-0-33` passes
- [ ] `cargo check --features nalgebra-0-34` passes
- [ ] `cargo test --features nalgebra-0-33` passes
- [ ] `cargo test --features nalgebra-0-34` passes
- [ ] Enabling both features simultaneously produces compile error (mutually exclusive)
- [ ] All conversions have comprehensive rustdoc
- [ ] Property tests verify mathematical equivalence

## Future Considerations

- **nalgebra-glm**: Could add conversions for glm types (vec2, vec3, mat3, quat)
- **simba**: nalgebra's scalar abstraction - ensure our Float trait is compatible
- **Const generics**: As nalgebra moves more toward const generics, conversions may simplify

## Dependencies

After implementation, update CLAUDE.md to document the new optional dependencies.
