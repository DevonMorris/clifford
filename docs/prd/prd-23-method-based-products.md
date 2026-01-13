# PRD-23: Generate Methods on Types Instead of Free Functions

**Status**: Draft
**Goal**: Improve API discoverability by generating products as methods on types rather than free functions

## Background

The current codegen generates products as free functions in a `products` module:

```rust
// Current: Free functions in products.rs
use crate::specialized::projective::dim3::products;

let line = products::exterior_point_point(&p1, &p2);
let motor = products::geometric_motor_motor(&m1, &m2);
let transformed = products::sandwich_motor_point(&motor, &point);
```

This approach has poor discoverability because:
1. Users must know the `products` module exists
2. Users must remember the exact function naming convention (`{product}_{lhs}_{rhs}`)
3. IDE autocomplete doesn't help - typing `point.` shows nothing useful
4. The API feels foreign compared to idiomatic Rust

## Problem Statement

### Issue 1: Poor Discoverability

When a user has a `Point` and wants to find what operations are available, typing `point.` in an IDE shows only the basic accessors and constructors. The rich algebra of products is hidden in a separate module.

Compare:
```rust
// Current: User must know about products module
use products::exterior_point_point;
let line = exterior_point_point(&p1, &p2);

// Desired: Discoverable via autocomplete
let line = p1.wedge(&p2);  // IDE shows this when typing "p1."
```

### Issue 2: Verbose Call Sites

Free functions require verbose, repetitive code:
```rust
// Current
let result = products::geometric_motor_motor(&m1, &m2);
let transformed = products::sandwich_motor_point(&motor, &point);

// Desired
let result = m1.geometric(&m2);
let transformed = motor.sandwich(&point);
```

### Issue 3: Inconsistent with Rust Idioms

Rust APIs typically use methods for operations. Free functions feel like C-style APIs rather than idiomatic Rust:
```rust
// Idiomatic Rust (what users expect)
let normalized = vector.normalize();
let length = vector.norm();
let result = a.dot(&b);

// Current clifford (feels foreign)
let line = products::exterior_point_point(&p1, &p2);
```

### Issue 4: Operator Overloading Disconnect

We implement `std::ops` traits for some operations (e.g., `Mul` for geometric product), but the explicit product functions are free functions. This creates an inconsistent experience:
```rust
// Works (via Mul trait)
let result = motor1 * motor2;

// But explicit call requires free function
let result = products::geometric_motor_motor(&m1, &m2);

// Would be more consistent as
let result = m1.geometric(&m2);
```

## Solution

### Phase 1: Generate Methods on Types

Generate `impl` blocks on each type with methods for all products where that type is the left-hand side.

#### 1.1 Method Naming Convention

| Product | Method Name | Example |
|---------|-------------|---------|
| `geometric` | `.geometric(&rhs)` | `motor.geometric(&other)` |
| `wedge` | `.wedge(&rhs)` | `point.wedge(&other)` |
| `antiwedge` | `.antiwedge(&rhs)` | `plane.antiwedge(&other)` |
| `inner` | `.inner(&rhs)` | `vector.inner(&other)` |
| `left_contract` | `.left_contract(&rhs)` | `point.left_contract(&plane)` |
| `right_contract` | `.right_contract(&rhs)` | `plane.right_contract(&point)` |
| `bulk_contraction` | `.bulk_contract(&rhs)` | `a.bulk_contract(&b)` |
| `weight_contraction` | `.weight_contract(&rhs)` | `a.weight_contract(&b)` |
| `scalar` | `.scalar_product(&rhs)` | `a.scalar_product(&b)` |
| `sandwich` | `.sandwich(&operand)` | `motor.sandwich(&point)` |
| `antisandwich` | `.antisandwich(&operand)` | `flector.antisandwich(&line)` |

#### 1.2 Generated Code Structure

```rust
// In generated/types.rs or generated/ops.rs

impl<T: Float> Point<T> {
    /// Wedge product: Point ∧ Point -> Line
    ///
    /// Computes the exterior (wedge) product of two points,
    /// yielding the line through both points.
    #[inline]
    pub fn wedge(&self, other: &Point<T>) -> Line<T> {
        Line::new_unchecked(
            // ... generated expression
        )
    }

    /// Wedge product: Point ∧ Line -> Plane
    #[inline]
    pub fn wedge(&self, other: &Line<T>) -> Plane<T> {
        // ... but wait, this conflicts with the above!
    }
}
```

#### 1.3 Handling Multiple RHS Types

The challenge is that `wedge` can take different types as the RHS:
- `Point::wedge(&Point)` -> `Line`
- `Point::wedge(&Line)` -> `Plane`

**Option A: Trait-based dispatch (recommended)**

```rust
/// Trait for types that can be wedged with Point
pub trait WedgeWith<Rhs> {
    type Output;
    fn wedge(&self, rhs: &Rhs) -> Self::Output;
}

impl<T: Float> WedgeWith<Point<T>> for Point<T> {
    type Output = Line<T>;
    fn wedge(&self, rhs: &Point<T>) -> Line<T> { ... }
}

impl<T: Float> WedgeWith<Line<T>> for Point<T> {
    type Output = Plane<T>;
    fn wedge(&self, rhs: &Line<T>) -> Plane<T> { ... }
}

// Usage: works with any RHS type
let line = p1.wedge(&p2);      // Point ∧ Point
let plane = p1.wedge(&line);   // Point ∧ Line
```

**Option B: Suffixed methods (simpler but verbose)**

```rust
impl<T: Float> Point<T> {
    pub fn wedge_point(&self, other: &Point<T>) -> Line<T> { ... }
    pub fn wedge_line(&self, other: &Line<T>) -> Plane<T> { ... }
}
```

**Option C: Operator overloading only**

Rely on `BitXor` (`^`) for wedge, `BitOr` (`|`) for antiwedge, `Mul` (`*`) for geometric. But this doesn't help with contractions or explicit naming.

**Recommendation**: Option A (trait-based) for primary products, keeping free functions available for cases where explicit type specification is needed.

### Phase 2: Keep Free Functions as Implementation

The free functions can remain as the implementation, with methods delegating to them:

```rust
// In generated/products.rs (unchanged)
pub fn wedge_point_point<T: Float>(a: &Point<T>, b: &Point<T>) -> Line<T> { ... }

// In generated/ops.rs (new)
impl<T: Float> WedgeWith<Point<T>> for Point<T> {
    type Output = Line<T>;
    #[inline]
    fn wedge(&self, rhs: &Point<T>) -> Line<T> {
        super::products::wedge_point_point(self, rhs)
    }
}
```

This allows:
- Backwards compatibility (free functions still work)
- Method-based API for discoverability
- Single source of truth for the actual computation

### Phase 3: Trait Hierarchy

Define a trait hierarchy for each product type:

```rust
// In clifford::ops (not generated, part of core library)

/// Wedge (exterior) product
pub trait Wedge<Rhs = Self> {
    type Output;
    fn wedge(&self, rhs: &Rhs) -> Self::Output;
}

/// Antiwedge (regressive) product
pub trait Antiwedge<Rhs = Self> {
    type Output;
    fn antiwedge(&self, rhs: &Rhs) -> Self::Output;
}

/// Geometric product
pub trait GeometricProduct<Rhs = Self> {
    type Output;
    fn geometric(&self, rhs: &Rhs) -> Self::Output;
}

/// Sandwich product: v × x × rev(v)
pub trait Sandwich<Operand> {
    type Output;
    fn sandwich(&self, operand: &Operand) -> Self::Output;
}

// etc.
```

Codegen then generates `impl` blocks for these traits.

### Phase 4: Update TOML Spec (Optional)

Allow TOML to specify method generation preferences:

```toml
[codegen]
generate_methods = true      # Generate impl blocks with methods
generate_free_functions = true  # Keep free functions (for backwards compat)
generate_traits = true       # Generate trait impls

[products.wedge]
Point_Point = "Line"
Point_Line = "Plane"
# ...
```

## Implementation Plan

### Step 1: Define Core Traits

Add to `src/ops.rs` (or `src/algebra/ops.rs`):
- `Wedge<Rhs>` trait
- `Antiwedge<Rhs>` trait
- `GeometricProduct<Rhs>` trait
- `Sandwich<Operand>` trait
- `LeftContract<Rhs>`, `RightContract<Rhs>` traits
- `BulkContract<Rhs>`, `WeightContract<Rhs>` traits

### Step 2: Update Codegen

Modify `crates/clifford-codegen/src/codegen/products.rs`:
1. Keep existing free function generation
2. Add trait impl generation for each product
3. Generate one `impl Trait<Rhs> for Lhs` block per product combination

### Step 3: Generate Trait Impls

For each product entry like `Point_Point = "Line"`, generate:

```rust
impl<T: Float> Wedge<Point<T>> for Point<T> {
    type Output = Line<T>;

    #[inline]
    fn wedge(&self, rhs: &Point<T>) -> Line<T> {
        wedge_point_point(self, rhs)
    }
}
```

### Step 4: Update Documentation

- Document the new method-based API
- Show examples using methods rather than free functions
- Keep free function documentation for advanced use cases

### Step 5: Regenerate All Algebras

```bash
for toml in algebras/*.toml; do
    cargo run --package clifford-codegen -- generate "$toml" --force
done
```

## Migration Path

### Backwards Compatibility

Free functions remain available, so existing code continues to work:
```rust
// Still works
let line = products::wedge_point_point(&p1, &p2);

// New preferred style
let line = p1.wedge(&p2);
```

### Deprecation (Future)

In a future version, consider:
1. Marking free functions as `#[doc(hidden)]` (still available but not in docs)
2. Eventually deprecating free functions if method API proves superior

## Example: Before and After

### Before (Current)

```rust
use clifford::specialized::projective::dim3::{Point, Line, Motor, products};

fn transform_line(motor: &Motor<f64>, p1: &Point<f64>, p2: &Point<f64>) -> Line<f64> {
    // Create line through two points
    let line = products::wedge_point_point(p1, p2);

    // Transform the line
    products::sandwich_motor_line(motor, &line)
}
```

### After (With Methods)

```rust
use clifford::specialized::projective::dim3::{Point, Line, Motor};
use clifford::ops::{Wedge, Sandwich};

fn transform_line(motor: &Motor<f64>, p1: &Point<f64>, p2: &Point<f64>) -> Line<f64> {
    // Create line through two points - discoverable via IDE!
    let line = p1.wedge(p2);

    // Transform the line - discoverable via IDE!
    motor.sandwich(&line)
}
```

## Testing

1. **Trait impl tests**: Verify trait methods produce same results as free functions
2. **Type inference tests**: Verify compiler correctly infers output types
3. **Documentation tests**: Update doc examples to use method syntax

## Success Criteria

1. All products available as methods on the LHS type
2. IDE autocomplete shows available products when typing `point.`
3. Free functions remain available for backwards compatibility
4. No performance regression (methods should inline to same code)
5. Documentation updated to prefer method syntax

## Open Questions

1. **Trait location**: Should traits live in `clifford::ops` or be generated per-algebra?
2. **Naming**: Should we use `wedge` or `exterior` for the method name?
3. **Self-products**: For `Point::wedge(&Point)`, should we also have `Point::wedge_self()` or rely on `p.wedge(&p)`?
4. **Operator symbols**: Should we implement `BitXor` (`^`) for wedge alongside the method?

## References

- [Rust API Guidelines - Naming](https://rust-lang.github.io/api-guidelines/naming.html)
- [nalgebra's approach](https://docs.rs/nalgebra) - uses methods like `.dot()`, `.cross()`
- Current codegen: `crates/clifford-codegen/src/codegen/products.rs`
