# PRD-4: Specialized 2D/3D Types

**Status**: Pending
**Goal**: Optimized implementations for common dimensions

## Deliverables

### 1. GA2D Module (`src/specialized/ga2d/`)

```rust
/// 2D scalar (grade 0)
#[derive(Clone, Copy, Debug, PartialEq)]
#[repr(transparent)]
pub struct Scalar<T: Float>(pub T);

/// 2D vector (grade 1): e1, e2
#[derive(Clone, Copy, Debug, PartialEq)]
#[repr(C)]
pub struct Vec2<T: Float> {
    pub x: T,  // e1
    pub y: T,  // e2
}

/// 2D bivector/pseudoscalar (grade 2): e12
#[derive(Clone, Copy, Debug, PartialEq)]
#[repr(transparent)]
pub struct Bivec2<T: Float>(pub T);

/// 2D rotor: scalar + bivector
#[derive(Clone, Copy, Debug, PartialEq)]
#[repr(C)]
pub struct Rotor2<T: Float> {
    pub s: T,  // scalar
    pub xy: T, // e12
}

impl<T: Float> Rotor2<T> {
    /// Create from angle (radians)
    pub fn from_angle(angle: T) -> Self;

    /// Apply rotation to vector
    pub fn rotate(&self, v: Vec2<T>) -> Vec2<T>;

    /// Compose rotations
    pub fn compose(&self, other: &Self) -> Self;
}
```

### 2. GA3D Module (`src/specialized/ga3d/`)

```rust
/// 3D scalar (grade 0)
#[derive(Clone, Copy, Debug, PartialEq)]
#[repr(transparent)]
pub struct Scalar<T: Float>(pub T);

/// 3D vector (grade 1): e1, e2, e3
#[derive(Clone, Copy, Debug, PartialEq)]
#[repr(C)]
pub struct Vec3<T: Float> {
    pub x: T,  // e1
    pub y: T,  // e2
    pub z: T,  // e3
}

/// 3D bivector (grade 2): e12, e13, e23
#[derive(Clone, Copy, Debug, PartialEq)]
#[repr(C)]
pub struct Bivec3<T: Float> {
    pub xy: T,  // e12
    pub xz: T,  // e13
    pub yz: T,  // e23
}

/// 3D trivector/pseudoscalar (grade 3): e123
#[derive(Clone, Copy, Debug, PartialEq)]
#[repr(transparent)]
pub struct Trivec3<T: Float>(pub T);

/// Full 3D multivector (8 components)
#[derive(Clone, Copy, Debug, PartialEq)]
#[repr(C)]
pub struct Multivector3<T: Float> {
    pub s: T,
    pub v: Vec3<T>,
    pub b: Bivec3<T>,
    pub t: T,
}

/// 3D rotor: scalar + bivector (even subalgebra)
#[derive(Clone, Copy, Debug, PartialEq)]
#[repr(C)]
pub struct Rotor3<T: Float> {
    pub s: T,
    pub b: Bivec3<T>,
}
```

### 3. Rotor Operations

```rust
impl<T: Float> Rotor3<T> {
    /// Create from angle and axis (unit vector)
    pub fn from_angle_axis(angle: T, axis: Vec3<T>) -> Self;

    /// Create from angle and plane (unit bivector)
    pub fn from_angle_plane(angle: T, plane: Bivec3<T>) -> Self;

    /// Create rotation from vector a to vector b
    pub fn from_vectors(a: Vec3<T>, b: Vec3<T>) -> Self;

    /// Apply rotation: R * v * R̃
    pub fn rotate(&self, v: Vec3<T>) -> Vec3<T>;

    /// Compose rotations: R2 * R1
    pub fn compose(&self, other: &Self) -> Self;

    /// Inverse rotation
    pub fn inverse(&self) -> Self;

    /// Spherical linear interpolation
    pub fn slerp(&self, other: &Self, t: T) -> Self;
}
```

### 4. Conversions

```rust
// Generic to specialized
impl<T: Float> From<Multivector<T, Euclidean3>> for Multivector3<T>;
impl<T: Float> From<Multivector3<T>> for Multivector<T, Euclidean3>;

// Component extraction
impl<T: Float> TryFrom<Multivector<T, Euclidean3>> for Vec3<T>;
impl<T: Float> TryFrom<Multivector<T, Euclidean3>> for Bivec3<T>;
```

## Files to Create

- `src/specialized/mod.rs`
- `src/specialized/ga2d/mod.rs`
- `src/specialized/ga2d/types.rs`
- `src/specialized/ga2d/ops.rs`
- `src/specialized/ga3d/mod.rs`
- `src/specialized/ga3d/types.rs`
- `src/specialized/ga3d/ops.rs`
- `src/transforms/mod.rs`
- `src/transforms/rotor.rs`

## Testing (proptest)

```rust
proptest! {
    #[test]
    fn specialized_matches_generic(
        a in arb_vec3::<f64>(),
        b in arb_vec3::<f64>(),
    ) {
        let spec = a.geometric(b);
        let gen_a = Multivector::<f64, Euclidean3>::from(a);
        let gen_b = Multivector::<f64, Euclidean3>::from(b);
        let gen = gen_a * gen_b;
        prop_assert!(Multivector3::from(gen).approx_eq(&spec, 1e-10));
    }

    #[test]
    fn rotor_preserves_norm(
        rotor in arb_unit_rotor3::<f64>(),
        v in arb_vec3::<f64>(),
    ) {
        let rotated = rotor.rotate(v);
        prop_assert!((v.norm() - rotated.norm()).abs() < 1e-10);
    }

    #[test]
    fn rotor_composition(
        r1 in arb_unit_rotor3::<f64>(),
        r2 in arb_unit_rotor3::<f64>(),
        v in arb_vec3::<f64>(),
    ) {
        let sequential = r2.rotate(r1.rotate(v));
        let composed = r2.compose(&r1).rotate(v);
        prop_assert!(sequential.approx_eq(&composed, 1e-10));
    }

    #[test]
    fn rotor_inverse(
        r in arb_unit_rotor3::<f64>(),
        v in arb_vec3::<f64>(),
    ) {
        let roundtrip = r.inverse().rotate(r.rotate(v));
        prop_assert!(roundtrip.approx_eq(&v, 1e-10));
    }
}
```

## Verification

- [ ] `cargo check` passes
- [ ] `cargo test` - specialized matches generic
- [ ] `cargo clippy` - no warnings
- [ ] All types have comprehensive rustdoc
