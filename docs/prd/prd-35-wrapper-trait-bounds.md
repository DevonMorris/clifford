# PRD-35: Add Type-Level Trait Bounds to Wrapper Types

**Status**: Complete
**Goal**: Add compile-time trait bounds to wrapper types and fix codegen to only generate valid type aliases

## Problem Statement

The wrapper types (`Unit`, `Bulk`, `Unitized`, `Ideal`, `Proper`, `Spacelike`, `Null`) currently have **no type-level trait bounds**. They only check trait requirements at the `impl` level. This allows the creation of type aliases that can never be used.

### Example: Invalid Aliases

Complex numbers (Cl(0,1,0)) only implement `Normed`, not `DegenerateNormed`. But codegen still generates:

```rust
// In src/specialized/complex/generated/types.rs
pub type BulkComplex<T> = crate::wrappers::Bulk<Complex<T>>;
pub type UnitizedComplex<T> = crate::wrappers::Unitized<Complex<T>>;
```

These are **dead code** because:
1. `Bulk<T>` requires `T: DegenerateNormed` in its methods
2. `Complex` only implements `Normed`, not `DegenerateNormed`
3. Any attempt to use `BulkComplex<f64>::try_new(...)` fails at compile time

The same issue affects:
- Hyperbolic numbers: `BulkHyperbolic`, `UnitizedHyperbolic`
- Dual numbers: `BulkDual`, `UnitizedDual`, `BulkDualUnit`, `UnitizedDualUnit`

### Root Cause

1. **No type-level bounds**: Wrapper structs are defined as:
   ```rust
   pub struct Bulk<T> { inner: T }  // No bound on T
   ```
   Instead of:
   ```rust
   pub struct Bulk<T: DegenerateNormed> { inner: T }
   ```

2. **Codegen blindly generates aliases**: The `generate_wrapper_aliases()` function generates `Bulk` and `Unitized` for ALL types, checking only if the algebra has degenerate basis vectors, not whether each type implements the required trait.

## Proposed Solution

### Part 1: Add Type-Level Trait Bounds to Wrappers

Add explicit trait bounds to each wrapper struct:

| Wrapper | Current | Proposed |
|---------|---------|----------|
| `Unit<T>` | `pub struct Unit<T>` | `pub struct Unit<T: Normed>` |
| `Bulk<T>` | `pub struct Bulk<T>` | `pub struct Bulk<T: DegenerateNormed>` |
| `Unitized<T>` | `pub struct Unitized<T>` | `pub struct Unitized<T: DegenerateNormed>` |
| `Ideal<T>` | `pub struct Ideal<T>` | `pub struct Ideal<T: DegenerateNormed>` |
| `Proper<T>` | `pub struct Proper<T>` | `pub struct Proper<T: IndefiniteNormed>` |
| `Spacelike<T>` | `pub struct Spacelike<T>` | `pub struct Spacelike<T: IndefiniteNormed>` |
| `Null<T>` | `pub struct Null<T>` | `pub struct Null<T: IndefiniteNormed>` |

Benefits:
- Invalid type aliases fail at compile time
- Better IDE support and error messages
- Self-documenting API requirements

### Part 2: Fix Codegen to Generate Only Valid Aliases

Update `generate_wrapper_aliases()` in `types.rs` to check whether each wrapper is valid for the algebra type:

| Wrapper | Required Trait | When to Generate |
|---------|----------------|------------------|
| `Unit<T>` | `Normed` | Non-degenerate algebras (Euclidean, hyperbolic, complex, etc.) |
| `Bulk<T>` | `DegenerateNormed` | Degenerate algebras only (PGA) |
| `Unitized<T>` | `DegenerateNormed` | Degenerate algebras only (PGA) |
| `Ideal<T>` | `DegenerateNormed` | Degenerate algebras only (PGA) |
| `Proper<T>` | `IndefiniteNormed` | Indefinite algebras only (Minkowski) |
| `Spacelike<T>` | `IndefiniteNormed` | Indefinite algebras only (Minkowski) |
| `Null<T>` | `IndefiniteNormed` | Indefinite algebras only (Minkowski) |

Decision logic based on algebra signature:
- **Non-degenerate, definite** (p > 0, q = 0, r = 0): `Unit` only (e.g., Euclidean)
- **Non-degenerate, indefinite** (p > 0, q > 0, r = 0): `Unit`, `Proper`, `Spacelike`, `Null` (e.g., Minkowski)
- **Degenerate** (r > 0): `Bulk`, `Unitized`, `Ideal` (e.g., PGA)
- **1D algebras** (hyperbolic, complex, dual): `Unit` only for single-grade types that are normalizable

### Part 3: Update Existing Impls to Match New Bounds

Ensure all `impl` blocks on wrapper types have matching trait bounds on `T`:

```rust
// Current (only on methods, not type)
impl<T> Unit<T> {
    pub fn try_new(inner: T) -> Option<Self>
    where T: Normed { ... }
}

// Proposed (bound on impl)
impl<T: Normed> Unit<T> {
    pub fn try_new(inner: T) -> Option<Self> { ... }
}
```

## Implementation Plan

### Step 1: Add Trait Bounds to Wrapper Structs

In `src/wrappers.rs`, update each wrapper struct definition:

```rust
// Before
#[derive(Clone, Copy)]
#[repr(transparent)]
pub struct Bulk<T> {
    inner: T,
}

// After
#[derive(Clone, Copy)]
#[repr(transparent)]
pub struct Bulk<T: DegenerateNormed> {
    inner: T,
}
```

### Step 2: Add Trait Bounds to All Impl Blocks

Update all `impl` blocks to have the trait bound on the impl rather than individual methods:

```rust
// Before
impl<T> Bulk<T> {
    pub fn try_new(inner: T) -> Option<Self>
    where
        T: DegenerateNormed,
        T::Scalar: Float,
    { ... }
}

// After
impl<T: DegenerateNormed> Bulk<T>
where
    T::Scalar: Float,
{
    pub fn try_new(inner: T) -> Option<Self> { ... }
}
```

### Step 3: Update Codegen Alias Generation

In `crates/clifford-codegen/src/codegen/types.rs`, update `generate_wrapper_aliases()`:

```rust
fn generate_wrapper_aliases(&self) -> TokenStream {
    let mut aliases = Vec::new();
    let sig = &self.spec.signature;

    // Determine which wrapper categories to generate
    let is_degenerate = sig.r > 0;  // Has zero-squaring basis vectors
    let is_indefinite = sig.p > 0 && sig.q > 0 && sig.r == 0;  // Non-degenerate but indefinite
    let is_definite = (sig.p > 0 || sig.q > 0) && sig.r == 0 && !(sig.p > 0 && sig.q > 0);

    for ty in &self.spec.types {
        if ty.alias_of.is_some() || ty.name == "Scalar" {
            continue;
        }

        let type_name = format_ident!("{}", ty.name);

        // Unit<T> for non-degenerate algebras (types that implement Normed only)
        if !is_degenerate {
            let should_generate_unit = ty.grades.len() == 1
                || (ty.grades.contains(&0) && ty.grades.contains(&2) && ty.grades.len() == 2);
            if should_generate_unit {
                // Generate Unit alias
            }
        }

        // Bulk<T>, Unitized<T>, Ideal<T> for degenerate algebras only
        if is_degenerate {
            // Generate Bulk, Unitized aliases
        }

        // Proper<T>, Spacelike<T>, Null<T> for indefinite algebras only
        if is_indefinite {
            // Generate Proper, Spacelike, Null aliases
        }
    }

    // ...
}
```

### Step 4: Regenerate All Algebras

Run `cargo build` to regenerate all algebras with corrected aliases.

### Step 5: Verify Invalid Aliases Are Gone

Check that:
- Complex: no `BulkComplex`, `UnitizedComplex` (only `UnitComplex` if appropriate)
- Hyperbolic: no `BulkHyperbolic`, `UnitizedHyperbolic`
- Dual: no `BulkDual`, `UnitizedDual`, etc.
- PGA: has `Bulk*`, `Unitized*` aliases
- Euclidean: has `Unit*` aliases

## Files to Modify

1. `src/wrappers.rs`
   - Add trait bounds to `Unit<T>`, `Bulk<T>`, `Unitized<T>`, `Ideal<T>`, `Proper<T>`, `Spacelike<T>`, `Null<T>` struct definitions
   - Update all `impl` blocks to have bounds on the impl

2. `crates/clifford-codegen/src/codegen/types.rs`
   - Update `generate_wrapper_aliases()` to conditionally generate aliases based on algebra signature

## Verification

1. `cargo build` - ensure generated code compiles
2. `cargo nextest run` - all tests pass
3. `cargo clippy` - no warnings
4. Verify:
   - `src/specialized/complex/generated/types.rs` has no `Bulk*` or `Unitized*` aliases
   - `src/specialized/hyperbolic/generated/types.rs` has no `Bulk*` or `Unitized*` aliases
   - `src/specialized/dual/generated/types.rs` has no `Bulk*` or `Unitized*` aliases
   - `src/specialized/projective/dim3/generated/types.rs` still has `Bulk*` and `Unitized*` aliases
   - `src/specialized/euclidean/dim3/generated/types.rs` still has `Unit*` aliases

## Success Criteria

1. All wrapper structs have type-level trait bounds
2. Invalid type aliases are no longer generated
3. All tests pass
4. No clippy warnings
5. Compile-time errors for invalid wrapper usage (e.g., `Bulk<Complex<f64>>`)

## Notes

- This is a breaking change if anyone was relying on the type aliases existing (even though they couldn't be used)
- The `PhantomData` pattern is not needed since the inner field provides the type
- Trait bounds on struct definitions are idiomatic Rust for types that fundamentally require traits
