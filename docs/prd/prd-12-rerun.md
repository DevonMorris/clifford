# PRD-12: Rerun Visualization Integration

**Status**: Draft
**Goal**: Enable visualization of geometric algebra types via the Rerun SDK (logging only, no viewer bundled)

## Motivation

[Rerun](https://rerun.io/) is a powerful open-source logging and visualization tool for computer vision and robotics. Integrating clifford with Rerun allows users to:

1. Debug geometric algebra computations visually
2. Visualize transformations (rotors, motors) in action
3. Create educational materials showing GA operations
4. Build robotics/CV applications using GA primitives

## Design Principles

### SDK-Only Dependency

**Critical**: We only provide type conversions for logging data. The rerun viewer is NOT bundled.

```toml
[dependencies]
rerun-0-28 = {
    package = "rerun",
    version = "0.28",
    optional = true,
    default-features = false,
    features = ["sdk"]  # SDK only, no viewer
}
```

Users must install the rerun viewer separately (`pip install rerun-sdk` or `cargo install rerun-cli`).

### Feature-Gated Integration

Following the nalgebra pattern (PRD-8):

```toml
[features]
rerun-0_28 = ["dep:rerun-0-28"]

[dependencies]
rerun-0-28 = { package = "rerun", version = "0.28", optional = true, default-features = false, features = ["sdk"] }
```

**Version rationale**: Start with 0.28.x (current stable as of Jan 2026). Add newer versions as needed.

## Type Mappings

### Euclidean 3D

| Clifford Type | Rerun Type | Conversion |
|---------------|------------|------------|
| `dim3::Vector<f32>` | `rerun::Vec3D` | Direction/displacement |
| `dim3::Vector<f32>` | `rerun::Position3D` | Via `AsPosition` wrapper |
| `dim3::Bivector<f32>` | `rerun::Vec3D` | Via dual vector (rotation axis) |
| `dim3::Rotor<f32>` | `rerun::RotationQuat` | Quaternion extraction |
| `dim3::Rotor<f32>` | `rerun::Transform3D` | Pure rotation transform |

### Euclidean 2D

| Clifford Type | Rerun Type | Conversion |
|---------------|------------|------------|
| `dim2::Vector<f32>` | `rerun::Vec2D` | Direction |
| `dim2::Vector<f32>` | `rerun::Position2D` | Via `AsPosition` wrapper |
| `dim2::Rotor<f32>` | `rerun::Angle` | Rotation angle extraction |

### Projective 3D (PGA)

| Clifford Type | Rerun Type | Conversion |
|---------------|------------|------------|
| `pga3d::Point<f32>` | `rerun::Position3D` | Cartesian from homogeneous |
| `pga3d::Line<f32>` | `(Position3D, Vec3D)` | Point + direction |
| `pga3d::Plane<f32>` | `(Position3D, Vec3D)` | Point + normal |
| `pga3d::Motor<f32>` | `rerun::Transform3D` | Full rigid transform |

### Projective 2D (PGA)

| Clifford Type | Rerun Type | Conversion |
|---------------|------------|------------|
| `pga2d::Point<f32>` | `rerun::Position2D` | Cartesian from homogeneous |
| `pga2d::Line<f32>` | Line segment data | Direction + point |
| `pga2d::Motor<f32>` | `rerun::Transform3D` | 2D rigid transform (embedded in 3D) |

## Helper Types

### Position vs Direction Disambiguation

Vectors can represent either positions (points) or directions (arrows). Use wrapper types:

```rust
/// Wrapper to interpret a vector as a position/point.
pub struct AsPosition<T>(pub T);

/// Wrapper to interpret a vector as a direction/arrow.
pub struct AsArrow<T>(pub T);

// Default: vectors become Vec3D (directions)
impl From<euclidean::dim3::Vector<f32>> for rerun::Vec3D { ... }

// Explicit: use AsPosition for points
impl From<AsPosition<euclidean::dim3::Vector<f32>>> for rerun::Position3D { ... }
```

### Usage Example

```rust
use clifford::specialized::euclidean::dim3::{Vector, Rotor, Bivector};
use clifford::visualization::rerun::AsPosition;
use std::f32::consts::FRAC_PI_4;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    // User must have rerun viewer running: `rerun`
    let rec = rerun::RecordingStreamBuilder::new("clifford_demo")
        .connect_tcp()?;  // Connect to running viewer

    // Log a point (position)
    let point = Vector::new(1.0_f32, 2.0, 3.0);
    rec.log("point", &rerun::Points3D::new([AsPosition(point)]))?;

    // Log direction arrows
    let direction = Vector::new(1.0_f32, 0.0, 0.0);
    rec.log("arrow", &rerun::Arrows3D::from_vectors([direction]))?;

    // Log a rotation transform
    let rotor = Rotor::from_angle_plane(FRAC_PI_4, Bivector::unit_xy());
    rec.log("rotation", &rerun::Transform3D::from(rotor))?;

    // Log a PGA motor (rotation + translation)
    use clifford::specialized::projective::dim3::Motor;
    let motor = Motor::from_rotation_translation(
        rotor,
        Vector::new(1.0, 0.0, 0.0)
    );
    rec.log("motor_transform", &rerun::Transform3D::from(motor))?;

    Ok(())
}
```

## Module Structure

Following the nalgebra pattern, add `rerun.rs` files in each type module:

```
src/specialized/
  euclidean/
    dim2/
      rerun.rs          # Vector, Bivector, Rotor -> Rerun
    dim3/
      rerun.rs          # Vector, Bivector, Trivector, Rotor -> Rerun
  projective/
    dim2/
      rerun.rs          # Point, Line, Motor -> Rerun
    dim3/
      rerun.rs          # Point, Line, Plane, Motor -> Rerun
  visualization/
    mod.rs              # Re-exports and helper types (AsPosition, AsArrow)
```

## Implementation Phases

### Phase 1: Euclidean 3D Core

1. Add `rerun-0_28` feature to Cargo.toml
2. Create `src/specialized/euclidean/dim3/rerun.rs`:
   - `Vector<f32> -> Vec3D`
   - `Vector<f32> -> Position3D` (via `AsPosition`)
   - `Bivector<f32> -> Vec3D` (dual vector)
   - `Rotor<f32> -> RotationQuat`
   - `Rotor<f32> -> Transform3D`
3. Create `src/specialized/visualization/mod.rs` with `AsPosition`, `AsArrow`
4. Property tests verifying component extraction

### Phase 2: Euclidean 2D + Full 3D

1. Create `src/specialized/euclidean/dim2/rerun.rs`
2. Add `Trivector` visualization (as scalar indicator)
3. Add `Multivector3` component extraction

### Phase 3: Projective GA

1. Create `src/specialized/projective/dim3/rerun.rs`:
   - `Point -> Position3D`
   - `Line -> (Position3D, Vec3D)`
   - `Plane -> (Position3D, Vec3D)`
   - `Motor -> Transform3D`
2. Create `src/specialized/projective/dim2/rerun.rs`

## Scalar Type Handling

Rerun uses `f32` internally. Strategy:

1. **Rerun conversions only for f32**: `From<Vector<f32>> for rerun::Vec3D`
2. **Ergonomic f64 -> f32 on our types**: `Vector<f64>::to_f32() -> Vector<f32>`
3. **Users chain conversions**: `my_vector_f64.to_f32().into()`

```rust
// Direct From impl for f32 (Rerun's native type)
impl From<Vector<f32>> for rerun::Vec3D { ... }

// Ergonomic scalar type conversion on clifford types (not rerun-specific)
impl<T: Float> Vector<T> {
    /// Convert to f32 vector (useful for rerun, graphics APIs, etc.)
    pub fn to_f32(&self) -> Vector<f32> {
        Vector::new(
            self.x().to_f32().unwrap_or(f32::NAN),
            self.y().to_f32().unwrap_or(f32::NAN),
            self.z().to_f32().unwrap_or(f32::NAN),
        )
    }

    /// Convert to f64 vector.
    pub fn to_f64(&self) -> Vector<f64> {
        Vector::new(
            self.x().to_f64().unwrap_or(f64::NAN),
            self.y().to_f64().unwrap_or(f64::NAN),
            self.z().to_f64().unwrap_or(f64::NAN),
        )
    }
}

// Usage:
let v: Vector<f64> = Vector::new(1.0, 2.0, 3.0);
let rerun_v: rerun::Vec3D = v.to_f32().into();  // Ergonomic chaining
```

**Note**: The `to_f32()`/`to_f64()` methods are generally useful and should be added to all specialized types, not just for rerun support. They can be implemented in a separate PR or as part of Phase 1.

## Testing

### Property-Based Tests

```rust
#[cfg(all(test, feature = "rerun-0_28"))]
mod tests {
    use super::*;
    use proptest::prelude::*;
    use crate::test_utils::ABS_DIFF_EQ_EPS;

    proptest! {
        #[test]
        fn vector_to_vec3d_components(v in any::<Vector<f32>>()) {
            let rerun_v: rerun::Vec3D = v.into();
            prop_assert!((rerun_v.x() - v.x()).abs() < f32::EPSILON);
            prop_assert!((rerun_v.y() - v.y()).abs() < f32::EPSILON);
            prop_assert!((rerun_v.z() - v.z()).abs() < f32::EPSILON);
        }

        #[test]
        fn rotor_to_quaternion_valid(r in any::<UnitRotor<f32>>()) {
            let q: rerun::RotationQuat = (*r).into();
            // Verify unit quaternion
            let norm_sq = q.x()*q.x() + q.y()*q.y() + q.z()*q.z() + q.w()*q.w();
            prop_assert!((norm_sq - 1.0).abs() < 1e-5);
        }

        #[test]
        fn pga_point_to_position(p in any::<projective::dim3::Point<f32>>()) {
            let pos: rerun::Position3D = p.into();
            // Verify Cartesian extraction matches
            prop_assert!((pos.x() - p.x()).abs() < 1e-5);
            prop_assert!((pos.y() - p.y()).abs() < 1e-5);
            prop_assert!((pos.z() - p.z()).abs() < 1e-5);
        }
    }
}
```

## CI Configuration

```yaml
# .github/workflows/ci.yml
jobs:
  check-rerun:
    runs-on: ubuntu-latest
    steps:
      - uses: actions/checkout@v4
      - run: cargo check --features rerun-0_28
```

**Note**: Only compilation check, not full test suite. Tests run without the rerun feature; this just verifies the feature compiles.

## Documentation

Each conversion documents:
1. Geometric interpretation
2. Convention mapping (quaternion ordering, etc.)
3. Precision notes (f64 -> f32)

Example:
```rust
/// Converts a 3D rotor to a Rerun rotation quaternion.
///
/// # Convention
///
/// Clifford rotor `R = s + xy·e₁₂ + xz·e₁₃ + yz·e₂₃` maps to
/// quaternion `(x, y, z, w) = (yz, -xz, xy, s)`.
///
/// This follows right-hand rotation conventions.
///
/// # Example
///
/// ```rust
/// use clifford::specialized::euclidean::dim3::{Rotor, Bivector};
///
/// let rotor = Rotor::from_angle_plane(
///     std::f32::consts::FRAC_PI_2,
///     Bivector::unit_xy()
/// );
/// let quat: rerun::RotationQuat = rotor.into();
/// ```
impl From<Rotor<f32>> for rerun::RotationQuat { ... }
```

## Files to Create/Modify

### New Files
- `src/specialized/euclidean/dim3/rerun.rs`
- `src/specialized/euclidean/dim2/rerun.rs`
- `src/specialized/projective/dim3/rerun.rs`
- `src/specialized/projective/dim2/rerun.rs`
- `src/specialized/visualization/mod.rs`
- `docs/prd/prd-12-rerun.md`

### Modified Files
- `Cargo.toml` - add rerun-0_28 feature
- `src/specialized/euclidean/dim3/mod.rs` - include rerun module
- `src/specialized/euclidean/dim2/mod.rs` - include rerun module
- `src/specialized/projective/dim3/mod.rs` - include rerun module
- `src/specialized/projective/dim2/mod.rs` - include rerun module
- `src/specialized/mod.rs` - include visualization module
- `.github/workflows/ci.yml` - add rerun test job
- `CLAUDE.md` - document new feature

## Verification Checklist

- [ ] `cargo check --features rerun-0_28` passes
- [ ] `cargo doc --features rerun-0_28` builds without warnings
- [ ] `cargo clippy --features rerun-0_28` passes
- [ ] All conversions have comprehensive rustdoc
- [ ] Property tests verify component extraction
- [ ] Example code compiles

## Resolved Decisions

1. **Rerun version**: 0.28.x (current stable)
2. **f64 support**: Ergonomic `to_f32()` methods on clifford types, rerun `From` impls only for f32
3. **SDK-only**: Use `default-features = false, features = ["sdk"]` - no viewer bundled
4. **Bivector visualization**: Dual vector arrow (simple, intuitive) - mesh option deferred

## Future Considerations

- **rerun-0_29+**: Add version features as Rerun releases
- **Custom archetypes**: `GaBivector`, `GaMotor` for richer visualization
- **Animation helpers**: Rotor slerp, motor interpolation visualization
- **Batch logging**: Efficient conversion for large point clouds
