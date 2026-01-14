# PRD-29: Semantic Field Names for Algebra Types

**Status**: Draft
**Goal**: Replace technical basis-blade field names with semantically meaningful names that convey geometric intent

## Background

The algebra TOML files currently use raw basis blade names for type fields:

```toml
[types.Motor]
fields = ["s", "e23", "e31", "e12", "e01", "e02", "e03", "e0123"]

[types.Line]
fields = ["e01", "e02", "e03", "e23", "e31", "e12"]
```

These names are mathematically correct but:
- Don't convey geometric meaning to users
- Require understanding of blade indexing conventions
- Make code harder to read and maintain

## Problem Statement

### Issue 1: Cryptic Field Names
`motor.e23()` tells you nothing about what that component represents geometrically. Compare to `motor.bx()` (bivector x-component) or `motor.rx()` (rotation about x).

### Issue 2: Inconsistent Conventions
Some types use semantic names (`x`, `y`, `z` for vectors) while others use blade names (`e23`, `e31`, `e12` for bivectors).

### Issue 3: Learning Barrier
New users must learn blade indexing before they can use the API effectively. Semantic names allow geometric intuition to guide usage.

## Solution

Replace blade-based field names with semantic names that convey geometric meaning.

### Naming Conventions

| Prefix | Meaning | Examples |
|--------|---------|----------|
| `x`, `y`, `z` | Spatial coordinates | Point, Vector |
| `w` | Homogeneous weight | Point |
| `nx`, `ny`, `nz` | Normal components | Plane |
| `d` | Distance/offset | Plane |
| `dx`, `dy`, `dz` | Direction components | Line |
| `mx`, `my`, `mz` | Moment components | Line, Motor |
| `bx`, `by`, `bz` | Bivector components | Bivector, Rotor, Motor |
| `s` | Scalar component | Rotor, Motor |
| `ps` | Pseudoscalar component | Motor |

### Euclidean 2D (`euclidean2.toml`)

| Type | Current | Proposed |
|------|---------|----------|
| Vector | `x`, `y` | `x`, `y` (unchanged) |
| Bivector | `xy` | `b` (single component) |
| Rotor | `s`, `xy` | `s`, `b` |

### Euclidean 3D (`euclidean3.toml`)

| Type | Current | Proposed |
|------|---------|----------|
| Vector | `x`, `y`, `z` | `x`, `y`, `z` (unchanged) |
| Bivector | `xy`, `xz`, `yz` | `bz`, `by`, `bx` |
| Rotor | `s`, `xy`, `xz`, `yz` | `s`, `bz`, `by`, `bx` |
| Trivector | `xyz` | `ps` |

**Note**: Bivector naming follows the convention that `e12` (xy-plane) corresponds to rotation about z-axis, hence `bz`. This matches `bx = e23`, `by = e31`, `bz = e12`.

### Projective 2D (`projective2.toml`)

| Type | Current | Proposed |
|------|---------|----------|
| Scalar | `s` | `s` (unchanged) |
| Point | `e1`, `e2`, `e0` | `x`, `y`, `w` |
| Line | `e01`, `e02`, `e12` | `mx`, `my`, `d` |
| Trivector | `e012` | `ps` |
| Motor | `e1`, `e2`, `e0`, `e012` | `x`, `y`, `w`, `ps` |
| Flector | `s`, `e01`, `e02`, `e12` | `s`, `mx`, `my`, `d` |

**Note**: 2D lines are characterized by moment (`mx`, `my`) and perpendicular distance `d` from origin.

### Projective 3D (`projective3.toml`)

| Type | Current | Proposed |
|------|---------|----------|
| Scalar | `s` | `s` (unchanged) |
| Point | `e1`, `e2`, `e3`, `e0` | `x`, `y`, `z`, `w` |
| Line | `e01`, `e02`, `e03`, `e23`, `e31`, `e12` | `mx`, `my`, `mz`, `dx`, `dy`, `dz` |
| Plane | `e023`, `e031`, `e012`, `e123` | `nx`, `ny`, `nz`, `d` |
| Quadvector | `e0123` | `ps` |
| Motor | `s`, `e23`, `e31`, `e12`, `e01`, `e02`, `e03`, `e0123` | `s`, `bx`, `by`, `bz`, `tx`, `ty`, `tz`, `ps` |
| Flector | `e1`, `e2`, `e3`, `e0`, `e023`, `e031`, `e012`, `e123` | `px`, `py`, `pz`, `pw`, `nx`, `ny`, `nz`, `d` |

**Motor components**:
- `s`: Scalar (cos of half-angle for pure rotation)
- `bx`, `by`, `bz`: Bivector/rotation components (sin of half-angle times axis)
- `tx`, `ty`, `tz`: Translation bivector components
- `ps`: Pseudoscalar (Study condition: `s*ps + bx*tx + by*ty + bz*tz = 0`)

**Flector components**:
- `px`, `py`, `pz`, `pw`: Point part (reflection point)
- `nx`, `ny`, `nz`, `d`: Plane part (reflection plane)

**Line components** (Plücker coordinates):
- `dx`, `dy`, `dz`: Direction vector (weight)
- `mx`, `my`, `mz`: Moment vector (bulk)
- Constraint: `dx*mx + dy*my + dz*mz = 0`

## Implementation Plan

### Step 1: Update TOML Specifications
Modify the `fields` arrays in each algebra TOML file to use the new semantic names.

### Step 2: Regenerate All Algebras
Run a build to regenerate all algebra code with new field names.

### Step 3: Update Extensions
Update manual extension methods in `src/specialized/` that reference field names:
- `euclidean/dim2/extensions.rs`
- `euclidean/dim3/extensions.rs`
- `projective/dim2/extensions.rs`
- `projective/dim3/extensions.rs`

### Step 4: Update Nalgebra Conversions
Update nalgebra interop code that maps between clifford types and nalgebra types:
- `src/interop/nalgebra.rs` (or similar)
- Ensure `From`/`Into` implementations use new field names

### Step 5: Update Rerun Conversions
Update rerun visualization code that converts clifford types to rerun primitives:
- `src/interop/rerun.rs` (or similar)
- Ensure visualization helpers use new field names

### Step 6: Update Documentation
- Update any examples in documentation
- Add field name explanations to type-level docs

### Step 7: Update Tests
If any tests reference fields by name, update them.

## API Changes

This is a **breaking change** for users who access fields directly:

```rust
// Before
let rx = motor.e23();
let tx = motor.e01();

// After
let rx = motor.bx();
let tx = motor.tx();
```

## Alternatives Considered

### Alternative 1: RGA Wiki Notation
Use `vx`, `vy`, `vz` for direction and `mx`, `my`, `mz` for moment (as in RGA wiki).

**Rejected**: `v` prefix is ambiguous (vector? velocity?). `d` for direction is clearer.

### Alternative 2: Full Word Names
Use `rotation_x`, `translation_x`, `normal_x`, etc.

**Rejected**: Too verbose for frequently-accessed fields. Single letters with prefixes balance clarity and brevity.

### Alternative 3: Keep Blade Names
Keep `e23`, `e01`, etc. and rely on documentation.

**Rejected**: Defeats the goal of making the API self-documenting and accessible.

## Testing

1. `cargo build` - Verify regeneration succeeds
2. `cargo nextest run` - Verify all tests pass
3. `cargo doc --open` - Verify documentation renders correctly
4. Manual review of generated accessor methods

## Success Criteria

1. All TOML files use semantic field names
2. Generated code compiles and tests pass
3. Documentation clearly explains field naming convention
4. API is more intuitive for new users

## References

- [RGA Wiki - Motor](https://rigidgeometricalgebra.org/wiki/index.php?title=Motor)
- [RGA Wiki - Line](https://rigidgeometricalgebra.org/wiki/index.php?title=Line)
- [RGA Wiki - Plane](https://rigidgeometricalgebra.org/wiki/index.php?title=Plane)
- [RGA Wiki - Point](https://rigidgeometricalgebra.org/wiki/index.php?title=Point)
