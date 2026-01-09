# PRD-4: Specialized 2D/3D Types

**Status**: Complete
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

Bidirectional conversions between specialized types and generic `Multivector`:

```rust
// GA2D conversions
impl<T: Float> From<Vec2<T>> for Multivector<T, Euclidean2>;
impl<T: Float> From<Bivec2<T>> for Multivector<T, Euclidean2>;
impl<T: Float> From<Rotor2<T>> for Multivector<T, Euclidean2>;

impl<T: Float> TryFrom<Multivector<T, Euclidean2>> for Vec2<T>;
impl<T: Float> TryFrom<Multivector<T, Euclidean2>> for Bivec2<T>;
impl<T: Float> TryFrom<Multivector<T, Euclidean2>> for Rotor2<T>;

// GA3D conversions
impl<T: Float> From<Vec3<T>> for Multivector<T, Euclidean3>;
impl<T: Float> From<Bivec3<T>> for Multivector<T, Euclidean3>;
impl<T: Float> From<Trivec3<T>> for Multivector<T, Euclidean3>;
impl<T: Float> From<Rotor3<T>> for Multivector<T, Euclidean3>;
impl<T: Float> From<Multivector3<T>> for Multivector<T, Euclidean3>;

impl<T: Float> TryFrom<Multivector<T, Euclidean3>> for Vec3<T>;
impl<T: Float> TryFrom<Multivector<T, Euclidean3>> for Bivec3<T>;
impl<T: Float> TryFrom<Multivector<T, Euclidean3>> for Trivec3<T>;
impl<T: Float> TryFrom<Multivector<T, Euclidean3>> for Rotor3<T>;
impl<T: Float> From<Multivector<T, Euclidean3>> for Multivector3<T>;
```

### 5. Conversion Consistency

**Critical requirement**: Operations must produce identical results regardless of whether performed on specialized or generic types. Round-trip conversions must preserve mathematical behavior.

```rust
// These must be equivalent:
let spec_result = vec3_a.wedge(&vec3_b);
let gen_result = Multivector::from(vec3_a).outer(&Multivector::from(vec3_b));
assert!(Bivec3::try_from(gen_result).unwrap().approx_eq(&spec_result));

// Round-trip must preserve operations:
let v: Vec3<f64> = /* ... */;
let roundtrip = Vec3::try_from(Multivector::from(v)).unwrap();
assert!(v.approx_eq(&roundtrip));
```

**Consistency properties to verify**:
- `specialized.op(other) ≈ Specialized::from(Generic::from(specialized).op(Generic::from(other)))`
- Round-trip conversions are identity (within floating-point tolerance)
- Grade extraction after conversion matches original specialized type
- Normalization is preserved through conversions

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

### Rotor Properties
```rust
proptest! {
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

### Conversion Consistency (Critical)
```rust
proptest! {
    // Round-trip conversions
    #[test]
    fn vec3_roundtrip(v in arb_vec3::<f64>()) {
        let generic = Multivector::<f64, Euclidean3>::from(v);
        let back = Vec3::try_from(generic).unwrap();
        prop_assert!(v.approx_eq(&back, 1e-10));
    }

    #[test]
    fn bivec3_roundtrip(b in arb_bivec3::<f64>()) {
        let generic = Multivector::<f64, Euclidean3>::from(b);
        let back = Bivec3::try_from(generic).unwrap();
        prop_assert!(b.approx_eq(&back, 1e-10));
    }

    #[test]
    fn rotor3_roundtrip(r in arb_unit_rotor3::<f64>()) {
        let generic = Multivector::<f64, Euclidean3>::from(r);
        let back = Rotor3::try_from(generic).unwrap();
        prop_assert!(r.approx_eq(&back, 1e-10));
    }

    // Operation consistency: specialized vs generic
    #[test]
    fn wedge_consistency(
        a in arb_vec3::<f64>(),
        b in arb_vec3::<f64>(),
    ) {
        let spec_result = a.wedge(&b);
        let gen_a = Multivector::<f64, Euclidean3>::from(a);
        let gen_b = Multivector::<f64, Euclidean3>::from(b);
        let gen_result = gen_a.outer(&gen_b);
        let gen_as_bivec = Bivec3::try_from(gen_result).unwrap();
        prop_assert!(spec_result.approx_eq(&gen_as_bivec, 1e-10));
    }

    #[test]
    fn dot_consistency(
        a in arb_vec3::<f64>(),
        b in arb_vec3::<f64>(),
    ) {
        let spec_result = a.dot(&b);
        let gen_a = Multivector::<f64, Euclidean3>::from(a);
        let gen_b = Multivector::<f64, Euclidean3>::from(b);
        let gen_result = gen_a.inner(&gen_b).scalar();
        prop_assert!((spec_result - gen_result).abs() < 1e-10);
    }

    #[test]
    fn geometric_consistency(
        a in arb_vec3::<f64>(),
        b in arb_vec3::<f64>(),
    ) {
        let spec_result = a.geometric(&b);
        let gen_a = Multivector::<f64, Euclidean3>::from(a);
        let gen_b = Multivector::<f64, Euclidean3>::from(b);
        let gen_result = gen_a * gen_b;
        let gen_as_mv3 = Multivector3::from(gen_result);
        prop_assert!(spec_result.approx_eq(&gen_as_mv3, 1e-10));
    }

    #[test]
    fn rotor_rotation_consistency(
        r in arb_unit_rotor3::<f64>(),
        v in arb_vec3::<f64>(),
    ) {
        // Specialized rotation
        let spec_result = r.rotate(v);

        // Generic sandwich product: R * v * R̃
        let gen_r = Multivector::<f64, Euclidean3>::from(r);
        let gen_v = Multivector::<f64, Euclidean3>::from(v);
        let gen_r_rev = gen_r.reverse();
        let gen_result = &(&gen_r * &gen_v) * &gen_r_rev;
        let gen_as_vec = Vec3::try_from(gen_result).unwrap();

        prop_assert!(spec_result.approx_eq(&gen_as_vec, 1e-10));
    }
}
```

## Verification

- [ ] `cargo check` passes
- [ ] `cargo test` - specialized matches generic
- [ ] `cargo clippy` - no warnings
- [ ] All types have comprehensive rustdoc
