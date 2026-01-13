# PRD-17.5: Regenerate All Algebras

**Status**: Draft
**Parent**: [PRD-17: Codegen Product Completeness and Documentation](prd-17-codegen-products.md)
**Goal**: Regenerate all existing algebras with the improved codegen

## Problem Statement

After completing PRD-17.1 through PRD-17.4, all existing algebras need to be regenerated to:

1. Include new product types (interior, left/right contraction, scalar)
2. Use correct naming (exterior instead of outer)
3. Have correct field documentation
4. Update extension files to use generated products

This is the final integration step that applies all improvements to the codebase.

## Prerequisites

Before executing this PRD, the following must be complete:

- [x] PRD-17.4: Guidance Updates (already complete)
- [ ] PRD-17.1: Missing Product Generation
- [ ] PRD-17.2: Exterior Product Naming
- [ ] PRD-17.3: Field Documentation Fix

## Algebras to Regenerate

| Algebra | TOML File | Generated Output | Specialized Output |
|---------|-----------|------------------|-------------------|
| Euclidean 2D | `algebras/euclidean2.toml` | `src/generated/euclidean2/` | `src/specialized/euclidean/dim2/generated/` |
| Euclidean 3D | `algebras/euclidean3.toml` | `src/generated/euclidean3/` | `src/specialized/euclidean/dim3/generated/` |
| Projective 2D | `algebras/projective2.toml` | `src/generated/projective2/` | `src/specialized/projective/dim2/generated/` |
| Projective 3D | `algebras/projective3.toml` | `src/generated/projective3/` | `src/specialized/projective/dim3/generated/` |

## Execution Steps

### Step 1: Update TOML Specifications

Enable new products in each TOML file:

```toml
[products]
geometric = true
exterior = true       # renamed from outer
interior = true       # NEW
left_contraction = true  # NEW (was partial)
right_contraction = true # NEW
scalar = true         # NEW
```

### Step 2: Regenerate All Algebras

```bash
# Regenerate all algebras
for toml in algebras/*.toml; do
    echo "Regenerating $(basename $toml)..."
    cargo run --package clifford-codegen -- generate "$toml" --force
done
```

Or individually:

```bash
cargo run --package clifford-codegen -- generate algebras/euclidean2.toml --force
cargo run --package clifford-codegen -- generate algebras/euclidean3.toml --force
cargo run --package clifford-codegen -- generate algebras/projective2.toml --force
cargo run --package clifford-codegen -- generate algebras/projective3.toml --force
```

### Step 3: Update Extension Files

Replace manual formulas with generated products:

#### projective/dim3/extensions.rs

```rust
// Before
impl<T: Float> Point<T> {
    pub fn join(&self, other: &Point<T>) -> Line<T> {
        // Manual formula
        let e01 = self.e1() * other.e2() - self.e2() * other.e1();
        // ...
    }

    pub fn left_contract_plane(&self, plane: &Plane<T>) -> T {
        // Manual formula
        self.e1() * plane.e023() + self.e2() * plane.e031() + ...
    }
}

// After
use super::generated::products;

impl<T: Float> Point<T> {
    pub fn join(&self, other: &Point<T>) -> Line<T> {
        products::exterior_point_point(self, other)
    }

    pub fn left_contract_plane(&self, plane: &Plane<T>) -> T {
        products::left_contract_point_plane(self, plane)
    }
}
```

#### Similar updates for:
- `projective/dim2/extensions.rs`
- `euclidean/dim2/extensions.rs`
- `euclidean/dim3/extensions.rs`

### Step 4: Run Verification

```bash
# Format
cargo fmt

# Lint
cargo clippy

# Build
cargo build

# Test all
cargo test

# Test with each nalgebra version
cargo test --no-default-features --features "serde proptest-support nalgebra-0_32"
cargo test --no-default-features --features "serde proptest-support nalgebra-0_33"
cargo test --no-default-features --features "serde proptest-support nalgebra-0_34"

# Documentation
cargo doc --no-deps

# License check
cargo deny check
```

### Step 5: Run Benchmarks

Verify no performance regression:

```bash
# Run benchmarks
cargo bench

# Compare with previous results if available
```

### Step 6: Commit Changes

```bash
# Stage all changes
git add algebras/
git add src/generated/
git add src/specialized/*/generated/
git add src/specialized/*/extensions.rs

# Commit
git commit -m "feat(codegen): regenerate all algebras with PRD-17 improvements

- Add interior product generation
- Add left/right contraction generation
- Add scalar product generation
- Rename outer → exterior products
- Fix field documentation to use TOML names
- Update extensions to use generated products

PRD-17, PRD-17.1, PRD-17.2, PRD-17.3, PRD-17.5

Co-Authored-By: Claude Opus 4.5 <noreply@anthropic.com>"
```

## Extension Updates Checklist

### projective/dim3/extensions.rs

| Method | Current | After |
|--------|---------|-------|
| `Point::join` | Manual formula | `products::exterior_point_point` |
| `Point::left_contract_line` | Manual formula | `products::left_contract_point_line` |
| `Point::left_contract_plane` | Manual formula | `products::left_contract_point_plane` |
| `Line::meet_plane` | Manual formula | `products::regressive_line_plane` or custom |
| `Line::left_contract_plane` | Manual formula | `products::left_contract_line_plane` |
| `Plane::meet` | Manual formula | `products::regressive_plane_plane` or custom |

### projective/dim2/extensions.rs

| Method | Current | After |
|--------|---------|-------|
| `Point::join` | Manual formula | `products::exterior_point_point` |
| `Line::meet` | Manual formula | `products::regressive_line_line` |

### euclidean/dim3/extensions.rs

| Method | Current | After |
|--------|---------|-------|
| `Vector::cross` | Manual formula | May keep (cross is specific to 3D) |
| `Vector::outer` | If exists | `products::exterior_vector_vector` |

### euclidean/dim2/extensions.rs

| Method | Current | After |
|--------|---------|-------|
| `Vector::outer` | If exists | `products::exterior_vector_vector` |

## Backward Compatibility

### Deprecated Outer Functions

The regenerated code will include deprecated aliases:

```rust
#[deprecated(since = "0.2.0", note = "use exterior_point_point instead")]
pub fn outer_point_point<T: Float>(a: &Point<T>, b: &Point<T>) -> Line<T> {
    exterior_point_point(a, b)
}
```

### Public API

The public API in extension files (`join`, `meet`, etc.) remains unchanged. Only the internal implementation switches to generated products.

## Testing Strategy

### Verify Generated Products Match Manual

Before removing manual formulas, verify equivalence:

```rust
#[test]
fn generated_exterior_matches_manual_join() {
    let p1 = Point::from_cartesian(1.0, 2.0, 3.0);
    let p2 = Point::from_cartesian(4.0, 5.0, 6.0);

    let via_generated = products::exterior_point_point(&p1, &p2);
    let via_manual = p1.join(&p2);  // before update

    assert_eq!(via_generated, via_manual);
}
```

### Property Tests Still Pass

All existing property tests must continue to pass after regeneration.

### New Product Tests

New products should have verification tests added in PRD-17.1.

## Deliverables

- [ ] Update `euclidean2.toml` with new product configuration
- [ ] Update `euclidean3.toml` with new product configuration
- [ ] Update `projective2.toml` with new product configuration
- [ ] Update `projective3.toml` with new product configuration
- [ ] Regenerate `src/generated/euclidean2/`
- [ ] Regenerate `src/generated/euclidean3/`
- [ ] Regenerate `src/generated/projective2/`
- [ ] Regenerate `src/generated/projective3/`
- [ ] Regenerate `src/specialized/euclidean/dim2/generated/`
- [ ] Regenerate `src/specialized/euclidean/dim3/generated/`
- [ ] Regenerate `src/specialized/projective/dim2/generated/`
- [ ] Regenerate `src/specialized/projective/dim3/generated/`
- [ ] Update `projective/dim3/extensions.rs` to use generated products
- [ ] Update `projective/dim2/extensions.rs` to use generated products
- [ ] Update `euclidean/dim3/extensions.rs` to use generated products
- [ ] Update `euclidean/dim2/extensions.rs` to use generated products
- [ ] Verify all tests pass
- [ ] Verify benchmarks show no regression
- [ ] Verify `cargo doc` builds successfully
- [ ] Create PR and merge

## Success Criteria

1. **All algebras regenerated**: Fresh generation with all new features
2. **Extensions use generated products**: No manual formulas remain (except documented exceptions)
3. **All tests pass**: Including nalgebra feature variants
4. **Documentation correct**: Field docs match TOML names
5. **No performance regression**: Benchmarks comparable to before

## Files Changed

| File | Action |
|------|--------|
| `algebras/euclidean2.toml` | Update - Enable new products |
| `algebras/euclidean3.toml` | Update - Enable new products |
| `algebras/projective2.toml` | Update - Enable new products |
| `algebras/projective3.toml` | Update - Enable new products |
| `src/generated/*/` | Regenerate - All files |
| `src/specialized/*/generated/` | Regenerate - All files |
| `src/specialized/*/extensions.rs` | Update - Use generated products |

## Estimated Scope

- TOML updates: ~20 lines each (4 files)
- Extension updates: ~50-100 lines each (4 files)
- Regenerated code: Automatic
- Testing/verification: ~2-4 hours

Total: ~1-2 days of work after PRD-17.1-17.3 are complete.
