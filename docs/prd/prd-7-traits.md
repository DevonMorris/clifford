# PRD-7: Geometric Algebra Trait Abstractions

**Status**: Pending
**Goal**: Trait-based abstractions for generic algorithms across different geometric algebras

## Background

Currently, specialized types like `Rotor2`, `Rotor3`, and (future) PGA motors share
common operations but have no formal relationship. This makes it impossible to write
generic algorithms that work across different geometric algebras.

By introducing traits that capture common GA patterns, we enable:
- Generic algorithms (animation, interpolation, transformation chains)
- Unified API across Euclidean, Projective, and Conformal GAs
- Type-safe algebra switching
- Educational clarity showing GA's unifying structure

## Deliverables

### 1. Core Blade Trait (`src/traits/blade.rs`)

```rust
/// A blade (simple k-vector) in a geometric algebra.
pub trait Blade: Sized + Copy + Clone {
    type Scalar: Float;

    /// The grade of this blade.
    fn grade(&self) -> usize;

    /// Squared norm (may be negative in non-Euclidean signatures).
    fn norm_squared(&self) -> Self::Scalar;

    /// Norm (magnitude).
    fn norm(&self) -> Self::Scalar {
        self.norm_squared().abs().sqrt()
    }

    /// Returns a normalized (unit) blade, if possible.
    fn normalize(&self) -> Option<Self>;

    /// Reverses the blade: negates based on grade.
    fn reverse(&self) -> Self;
}
```

### 2. Vector Trait (`src/traits/vector.rs`)

```rust
/// A vector (grade-1 element) in a geometric algebra.
pub trait Vector: Blade {
    /// Dot product (inner product of vectors).
    fn dot(&self, other: &Self) -> Self::Scalar;

    /// Squared magnitude (same as norm_squared for vectors).
    fn magnitude_squared(&self) -> Self::Scalar {
        self.dot(self)
    }
}
```

### 3. Rotor Trait (`src/traits/rotor.rs`)

```rust
/// A rotor (even-grade element representing rotation).
pub trait Rotor: Sized + Copy + Clone {
    type Scalar: Float;
    type Vector: Vector<Scalar = Self::Scalar>;

    /// Identity rotation (no-op).
    fn identity() -> Self;

    /// Applies this rotation to a vector.
    fn rotate(&self, v: Self::Vector) -> Self::Vector;

    /// Composes two rotations: self then other.
    fn compose(&self, other: Self) -> Self;

    /// Returns the inverse rotation.
    fn inverse(&self) -> Self;

    /// Spherical linear interpolation.
    fn slerp(&self, other: Self, t: Self::Scalar) -> Self;

    /// Squared norm.
    fn norm_squared(&self) -> Self::Scalar;

    /// Normalizes to unit rotor.
    fn normalize(&self) -> Self;
}
```

### 4. Motor Trait (`src/traits/motor.rs`)

```rust
/// A motor (rigid body transformation: rotation + translation).
///
/// Motors generalize rotors to include translations. In Euclidean GA,
/// a motor is just a rotor. In PGA, motors can represent any rigid motion.
pub trait Motor: Rotor {
    type Point;

    /// Creates a pure translation motor.
    fn from_translation(direction: Self::Vector, distance: Self::Scalar) -> Self;

    /// Applies this transformation to a point.
    fn transform_point(&self, p: Self::Point) -> Self::Point;

    /// Logarithm (for interpolation and analysis).
    fn log(&self) -> Self;

    /// Exponential (from bivector to motor).
    fn exp(bivector: Self) -> Self;
}
```

### 5. GeometricAlgebra Trait (`src/traits/algebra.rs`)

```rust
/// A geometric algebra with associated types for each grade.
pub trait GeometricAlgebra {
    type Scalar: Float;
    type Vector: Vector<Scalar = Self::Scalar>;
    type Bivector: Blade<Scalar = Self::Scalar>;
    type Rotor: Rotor<Scalar = Self::Scalar, Vector = Self::Vector>;

    /// Dimension of the algebra (number of basis vectors).
    const DIM: usize;

    /// Metric signature (p, q, r).
    fn signature() -> (usize, usize, usize);
}

/// Marker for algebras that support motors (PGA, etc.).
pub trait MotorAlgebra: GeometricAlgebra {
    type Point;
    type Motor: Motor<Scalar = Self::Scalar, Vector = Self::Vector, Point = Self::Point>;
}
```

### 6. Implementations

```rust
// 2D Euclidean
pub struct Ga2d;

impl GeometricAlgebra for Ga2d {
    type Scalar = f64; // or generic
    type Vector = Vec2<f64>;
    type Bivector = Bivec2<f64>;
    type Rotor = Rotor2<f64>;

    const DIM: usize = 2;

    fn signature() -> (usize, usize, usize) { (2, 0, 0) }
}

// 3D Euclidean
pub struct Ga3d;

impl GeometricAlgebra for Ga3d {
    type Scalar = f64;
    type Vector = Vec3<f64>;
    type Bivector = Bivec3<f64>;
    type Rotor = Rotor3<f64>;

    const DIM: usize = 3;

    fn signature() -> (usize, usize, usize) { (3, 0, 0) }
}

// 3D PGA (after PRD-5)
pub struct Pga3d;

impl GeometricAlgebra for Pga3d { ... }
impl MotorAlgebra for Pga3d {
    type Point = Point<f64>;
    type Motor = Motor<f64>;
}
```

### 7. Generic Algorithms (`src/algorithms/`)

```rust
/// Generic animation/interpolation
pub fn interpolate<R: Rotor>(start: R, end: R, t: R::Scalar) -> R {
    start.slerp(end, t)
}

/// Chain multiple transformations
pub fn chain<R: Rotor>(transforms: impl IntoIterator<Item = R>) -> R {
    transforms.into_iter().fold(R::identity(), |a, b| a.compose(b))
}

/// Animate along a path of keyframes
pub fn animate_path<R: Rotor>(
    keyframes: &[R],
    t: R::Scalar,
) -> R {
    // Find segment and interpolate
    ...
}

/// Smooth step between rotations
pub fn smooth_step<R: Rotor>(
    start: R,
    end: R,
    t: R::Scalar,
) -> R {
    let smooth_t = t * t * (R::Scalar::from_i8(3) - R::Scalar::TWO * t);
    start.slerp(end, smooth_t)
}
```

## Files to Create

- `src/traits/mod.rs`
- `src/traits/blade.rs`
- `src/traits/vector.rs`
- `src/traits/rotor.rs`
- `src/traits/motor.rs`
- `src/traits/algebra.rs`
- `src/algorithms/mod.rs`
- `src/algorithms/interpolation.rs`
- `src/algorithms/transforms.rs`

## Testing (proptest)

```rust
proptest! {
    #[test]
    fn generic_rotor_identity<R: Rotor>(v in arb_vector::<R::Vector>()) {
        let rotated = R::identity().rotate(v);
        prop_assert!(rotated.approx_eq(&v, 1e-10));
    }

    #[test]
    fn generic_rotor_inverse<R: Rotor>(
        r in arb_rotor::<R>(),
        v in arb_vector::<R::Vector>(),
    ) {
        let roundtrip = r.inverse().rotate(r.rotate(v));
        prop_assert!(roundtrip.approx_eq(&v, 1e-10));
    }

    #[test]
    fn generic_slerp_endpoints<R: Rotor>(
        r1 in arb_rotor::<R>(),
        r2 in arb_rotor::<R>(),
    ) {
        let at_0 = r1.slerp(r2, R::Scalar::ZERO);
        let at_1 = r1.slerp(r2, R::Scalar::ONE);
        prop_assert!(at_0.approx_eq(&r1, 1e-10));
        prop_assert!(at_1.approx_eq(&r2, 1e-10));
    }

    #[test]
    fn chain_is_compose<R: Rotor>(
        r1 in arb_rotor::<R>(),
        r2 in arb_rotor::<R>(),
        r3 in arb_rotor::<R>(),
        v in arb_vector::<R::Vector>(),
    ) {
        let chained = chain([r1, r2, r3]);
        let manual = r1.compose(r2).compose(r3);
        prop_assert!(chained.rotate(v).approx_eq(&manual.rotate(v), 1e-10));
    }
}
```

## Design Considerations

### Why Traits Over Naming Conventions?

| Aspect | Naming Convention | Trait-Based |
|--------|-------------------|-------------|
| Generic algorithms | Not possible | `fn animate<R: Rotor>(...)` |
| Compile-time safety | Runtime errors | Type errors |
| Discoverability | Read docs | IDE autocomplete |
| Future algebras | Copy patterns | Implement traits |

### Trait Complexity

Start minimal and expand:
1. `Rotor` trait first (most immediate value)
2. `Motor` trait when PGA lands
3. `GeometricAlgebra` for full algebra abstraction
4. `Blade` and `Vector` for completeness

### Performance

- All traits should be object-safe where possible
- Use `#[inline]` on trait methods
- Consider `const fn` where applicable
- No dynamic dispatch in hot paths

## Dependencies

- **Requires**: PRD-4 (specialized types to implement traits)
- **Enhances**: PRD-5 (PGA motors implement Motor trait)
- **Enhances**: PRD-6 (CGA versors implement traits)

## Verification

- [ ] `cargo check` passes
- [ ] `cargo test` - generic algorithms work with 2D and 3D
- [ ] `cargo clippy` - no warnings
- [ ] Trait documentation explains geometric meaning
- [ ] Examples show generic algorithm usage
