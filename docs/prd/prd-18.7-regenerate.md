# PRD-18.7: Regenerate All Algebras

**Status**: Draft
**Parent**: [PRD-18](prd-18-constraint-redesign.md)
**Goal**: Regenerate all algebra implementations with the new constraint system

## Overview

After implementing PRD-18.1 through PRD-18.6, regenerate all algebras to:
1. Use the new `Normed` trait hierarchy
2. Use geometry-specific wrappers
3. Remove user-configured constraints
4. Use inferred constraints
5. Use simplified constructors
6. Use bespoke Arbitrary implementations

## Deliverables

### Updated TOML Files

Remove all constraint sections from:
- `algebras/euclidean2.toml`
- `algebras/euclidean3.toml`
- `algebras/projective3.toml`

### Regeneration Commands

```bash
# Regenerate all algebras
for toml in algebras/*.toml; do
    cargo run --package clifford-codegen -- generate "$toml" --force
done

# Format generated code
cargo fmt

# Verify no warnings
cargo clippy

# Run all tests
cargo test

# Build docs
cargo doc --no-deps
```

### Expected Generated Files

#### Euclidean 2D (`src/specialized/euclidean/dim2/generated/`)
- `types.rs` - Vector, Bivector, Rotor with `Normed` impl
- `traits.rs` - Arbitrary implementations
- `products.rs` - No changes
- `conversions.rs` - No changes

#### Euclidean 3D (`src/specialized/euclidean/dim3/generated/`)
- `types.rs` - Vector, Bivector, Trivector, Rotor, Even with `Normed` impl
- `traits.rs` - Arbitrary implementations
- `products.rs` - No changes
- `conversions.rs` - No changes

#### Projective 3D (`src/specialized/projective/dim3/generated/`)
- `types.rs` - Point, Line, Plane, Motor, Flector with `Normed` + `DegenerateNormed` impl
- `traits.rs` - Arbitrary implementations with constraint solving
- `products.rs` - No changes
- `conversions.rs` - No changes

### Type Alias Generation

Each algebra module should have type aliases:

```rust
// euclidean/dim3/generated/types.rs
pub type UnitVector<T> = Unit<Vector<T>>;
pub type UnitBivector<T> = Unit<Bivector<T>>;
pub type UnitRotor<T> = Unit<Rotor<T>>;

// projective/dim3/generated/types.rs
pub type BulkMotor<T> = Bulk<Motor<T>>;
pub type BulkFlector<T> = Bulk<Flector<T>>;
pub type BulkLine<T> = Bulk<Line<T>>;
pub type IdealPoint<T> = Ideal<Point<T>>;
pub type IdealPlane<T> = Ideal<Plane<T>>;
```

## Verification Checklist

### Build Verification
- [ ] `cargo build` succeeds
- [ ] `cargo clippy` has no warnings
- [ ] `cargo doc --no-deps` succeeds
- [ ] `cargo deny check` passes

### Test Verification
- [ ] `cargo test` passes all tests
- [ ] Property tests verify constraints
- [ ] Existing extension tests still work

### API Verification
- [ ] `new_unchecked()` works for all types
- [ ] `new_checked()` validates constraints
- [ ] `Normed` trait implemented for all types
- [ ] `DegenerateNormed` implemented for PGA types
- [ ] Type aliases exported correctly

### Documentation Verification
- [ ] Generated docs include constraint equations
- [ ] Constructor docs explain usage
- [ ] Type aliases documented

## Rollback Plan

If issues are found:
1. Git stash/revert generated changes
2. Fix codegen
3. Regenerate

## Success Criteria

1. All algebras regenerated successfully
2. No compilation errors or warnings
3. All tests pass
4. Documentation builds correctly
5. Type aliases available in public API
6. Migration guide examples work
