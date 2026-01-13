# PRD-21: 2D PGA Dual Motor Representation

## Overview

This PRD addresses a fundamental inconsistency between 2D and 3D PGA motor representations. In 3D PGA, motors are stored in "dual form" (including the pseudoscalar e0123), which allows `antisandwich_motor_point` to work directly. In 2D PGA, motors are stored in "standard form" (grades 0+2, without the pseudoscalar e012), requiring a manual dual-sandwich workaround for point transformations.

This PRD proposes aligning 2D PGA with 3D PGA by using the dual representation for motors.

## Problem Statement

### Current State

**3D PGA Motor** (grades 0+2+4):
- Identity: `e0123=1` (pseudoscalar)
- Translation (3,0,0): `e23=1.5, e0123=1`
- Rotation Z 90°: `e03=0.707, e0123=0.707`
- `antisandwich_motor_point` works correctly

**2D PGA Motor** (grades 0+2):
- Identity: `s=1`
- Translation (3,4): `s=1, e01=-1.5, e02=-2`
- Rotation 90°: `s=0.707, e12=-0.707`
- `antisandwich_motor_point` **does not work** (returns zeros or wrong values)

### Root Cause

The antisandwich formula requires the versor to include the pseudoscalar for the "anti-identity" to work:
- 3D: Pseudoscalar `e0123` is grade 4 (even) → included in Motor
- 2D: Pseudoscalar `e012` is grade 3 (odd) → NOT in Motor, IS in Flector

### Current Workaround

The 2D `Motor::transform_point` uses a manual "dual sandwich":
```rust
pub fn transform_point(&self, p: &Point<T>) -> Point<T> {
    // complement(P) to get a line
    let p_as_line = Line::new(p.e0(), p.e2(), -p.e1());
    // Apply motor sandwich to line
    let transformed_line = products::sandwich_motor_line(self, &p_as_line);
    // complement back to point
    Point::new(-transformed_line.e02(), transformed_line.e01(), transformed_line.e12())
}
```

This works but is:
- Inconsistent with 3D (which uses `antisandwich_motor_point`)
- Confusing (why the manual duality?)
- Not using generated products

## Investigation Results

Testing confirmed that if we convert a 2D Motor to its dual representation (stored as a Flector), then `antisandwich_flector_point` produces correct results:

| Motor (s, e12, e01, e02) | Dual Flector (e1, e2, e0, e012) | Mapping |
|--------------------------|--------------------------------|---------|
| s | e012 | Identity scalar ↔ pseudoscalar |
| e12 | e0 | Rotation bivector ↔ weight |
| e01 | e2 | Translation ↔ dual direction |
| e02 | -e1 | Translation ↔ dual direction (negated) |

With this mapping:
- Identity: `e012=1` → antisandwich preserves points ✓
- Rotation 90°: `e0=-0.707, e012=0.707` → (1,0) → (0,1) ✓
- Translation (3,4): `e1=2, e2=-1.5, e012=1` → (1,2) → (4,6) ✓

## Proposed Solution

**Swap the grade definitions of Motor and Flector in `algebras/projective2.toml`:**

### Before
```toml
[types.Motor]
grades = [0, 2]
fields = ["s", "e12", "e01", "e02"]

[types.Flector]
grades = [1, 3]
fields = ["e1", "e2", "e0", "e012"]
```

### After
```toml
[types.Motor]
grades = [1, 3]
description = "2D rigid transformation (rotation + translation) in dual form"
fields = ["e1", "e2", "e0", "e012"]
versor = true
sandwich = { targets = ["Point", "Line", "Motor", "Flector"] }

[types.Flector]
grades = [0, 2]
description = "2D reflection (improper isometry) in dual form"
fields = ["s", "e12", "e01", "e02"]
versor = true
sandwich = { targets = ["Point", "Line", "Motor", "Flector"] }
```

## Implementation Plan

### Phase 1: Update Algebra Definition
1. Modify `algebras/projective2.toml` to swap Motor/Flector grades
2. Regenerate code: `cargo run --package clifford-codegen -- generate algebras/projective2.toml -o src/specialized/projective/dim2/generated --force`

### Phase 2: Update Extension Methods
Update `src/specialized/projective/dim2/extensions.rs`:

#### Motor Factory Methods
```rust
impl<T: Float> Motor<T> {
    /// Identity motor in dual form: e012 = 1
    pub fn identity() -> Self {
        Self::new(T::zero(), T::zero(), T::zero(), T::one())
    }

    /// Rotation motor in dual form
    /// Standard: s=cos(θ/2), e12=-sin(θ/2)
    /// Dual: e0=-sin(θ/2), e012=cos(θ/2)
    pub fn from_rotation(angle: T) -> Self {
        let half = angle / T::TWO;
        Self::new(T::zero(), T::zero(), -half.sin(), half.cos())
    }

    /// Translation motor in dual form
    /// Standard: s=1, e01=-dx/2, e02=-dy/2
    /// Dual: e1=dy/2, e2=-dx/2, e012=1 (note sign mapping)
    pub fn from_translation(dx: T, dy: T) -> Self {
        Self::new(dy / T::TWO, -dx / T::TWO, T::zero(), T::one())
    }

    /// Transform point using antisandwich (now works directly!)
    pub fn transform_point(&self, p: &Point<T>) -> Point<T> {
        products::antisandwich_motor_point(self, p)
    }
}
```

#### Flector Factory Methods
```rust
impl<T: Float> Flector<T> {
    /// Reflection through a line (now uses standard form)
    pub fn from_line(line: &Line<T>) -> Self {
        let l = line.unitized();
        Self::new(T::zero(), l.e12(), l.e01(), l.e02())
    }

    /// Transform point - flectors now use standard form
    /// Need to determine correct product (sandwich or antisandwich)
    pub fn transform_point(&self, p: &Point<T>) -> Point<T> {
        // TBD: verify which product is correct for swapped representation
        products::sandwich_flector_point(self, p)
    }
}
```

### Phase 3: Update nalgebra Conversions
Update `src/specialized/projective/dim2/nalgebra.rs`:

```rust
impl<T: Float + na::RealField> From<na::Isometry2<T>> for Motor<T> {
    fn from(iso: na::Isometry2<T>) -> Self {
        // Convert to dual form
        let angle = iso.rotation.angle();
        let tx = iso.translation.x;
        let ty = iso.translation.y;

        let rotation = Motor::from_rotation(angle);
        let translation = Motor::from_translation(tx, ty);
        translation.compose(&rotation)
    }
}

impl<T: Float + na::RealField> From<Motor<T>> for na::Isometry2<T> {
    fn from(m: Motor<T>) -> Self {
        // Extract from dual form
        let angle = m.rotation_angle();
        let trans = m.translation();
        na::Isometry2::new(na::Vector2::new(trans.x(), trans.y()), angle)
    }
}
```

### Phase 4: Update Tests
1. Update all tests in `extensions.rs` to use new representation
2. Verify nalgebra interop tests still pass
3. Add property-based tests verifying antisandwich works correctly

### Phase 5: Update mod.rs Exports
Ensure wrapper type aliases are correctly exported:
- `BulkMotor<T>` - normalized motor (bulk_norm = 1)
- `UnitizedPoint<T>` - finite point (weight = 1)
- etc.

## Verification

### Must Pass
1. All existing tests (may need updates for new representation)
2. Motor identity preserves points
3. Motor rotation transforms correctly
4. Motor translation transforms correctly
5. Motor composition matches nalgebra isometry composition
6. Flector reflection works correctly
7. nalgebra round-trip conversions preserve transformations

### Property Tests
```rust
proptest! {
    #[test]
    fn motor_antisandwich_works(m in any::<BulkMotor<f64>>(), p in any::<UnitizedPoint<f64>>()) {
        // antisandwich should now work directly
        let result = products::antisandwich_motor_point(&m, &p);
        // Should match the old dual-sandwich approach
        // (verify equivalence during migration)
    }
}
```

## Risks and Mitigations

### Risk: Breaking Change for Users
**Mitigation**: This is an internal representation change. The public API (`Motor::from_rotation`, `Motor::transform_point`, etc.) remains the same. Users don't directly access the internal components.

### Risk: Flector Transformation Correctness
**Mitigation**: Carefully verify which product (sandwich vs antisandwich) is correct for the swapped Flector. May need investigation similar to what we did for Motor.

### Risk: nalgebra Conversion Correctness
**Mitigation**: Comprehensive property-based tests comparing clifford and nalgebra operations.

## Success Criteria

1. `Motor::transform_point` uses `antisandwich_motor_point` (no manual duality)
2. 2D and 3D PGA have consistent patterns for motor transformations
3. All existing tests pass (with updates for new representation)
4. nalgebra interop works correctly
5. Code is simpler and more maintainable

## References

- [Rigid Geometric Algebra Wiki](https://rigidgeometricalgebra.org/wiki/index.php?title=Main_Page)
- Investigation transcript showing dual representation discovery
- 3D PGA motor implementation (`src/specialized/projective/dim3/extensions.rs`)
