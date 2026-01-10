# PRD-6: Conformal Geometric Algebra (CGA)

**Status**: Pending
**Goal**: Conformal GA support for circles, spheres, and uniform geometric operations

## Reference

**Primary Resource**: https://conformalgeometricalgebra.org/wiki/index.php?title=Main_Page

This wiki provides comprehensive formulas for CGA operations, geometric constraints,
and transformation versors. All implementations should be derived from and verified
against this reference.

## Lessons Learned from PGA (PRD-5)

### 1. Symbolic Derivation is Mandatory

**Never hardcode algebraic formulas.** Use SymPy for all derivations:

```python
# derivations/src/clifford_derivations/cga.py
from sympy import symbols, expand, simplify, rust_code

@with_timeout(300)
def derive_translator_sandwich():
    """Derive T * P * T̃ for CGA translator."""
    # ... symbolic computation
    return {component: simplify(expr) for component, expr in result.items()}

def generate_rust_translator(result):
    """Generate Rust code using sympy's rust_code()."""
    for name, expr in result.items():
        print(f"let {name} = {rust_code(expr)};")
```

### 2. Property Tests Against nalgebra

Every CGA operation must have property tests comparing against nalgebra equivalents:

```rust
proptest! {
    #[test]
    fn translator_matches_nalgebra(
        dx in -10.0f64..10.0, dy in -10.0f64..10.0, dz in -10.0f64..10.0,
        px in -10.0f64..10.0, py in -10.0f64..10.0, pz in -10.0f64..10.0,
    ) {
        let cga_trans = Translator::new(dx, dy, dz);
        let cga_point = Point::new(px, py, pz);
        let cga_result = cga_trans.transform_point(&cga_point);

        let na_trans = na::Translation3::new(dx, dy, dz);
        let na_point = na::Point3::new(px, py, pz);
        let na_result = na_trans * na_point;

        prop_assert!(abs_diff_eq!(cga_result.x(), na_result.x, epsilon = ABS_DIFF_EQ_EPS));
        prop_assert!(abs_diff_eq!(cga_result.y(), na_result.y, epsilon = ABS_DIFF_EQ_EPS));
        prop_assert!(abs_diff_eq!(cga_result.z(), na_result.z, epsilon = ABS_DIFF_EQ_EPS));
    }
}
```

### 3. Geometric Constraints as Invariants

CGA elements must satisfy geometric constraints (Reference: wiki Geometric_constraint page):

| Element | Constraint | Meaning |
|---------|------------|---------|
| **Null Vector** | `x · x = 0` | Valid conformal point |
| **Round** | `X ∧ X = 0, X · X ≠ 0` | Sphere/circle with real radius |
| **Flat** | `X ∧ e∞ = 0` | Plane/line (contains point at infinity) |
| **Translator** | `T T̃ = 1` | Unit versor |
| **Dilator** | `D D̃ = 1` | Unit versor |

Add verification methods:
```rust
impl<T: Float> Point<T> {
    pub fn is_null(&self, epsilon: T) -> bool;
}

impl<T: Float> Sphere<T> {
    pub fn is_real(&self, epsilon: T) -> bool;  // radius² > 0
    pub fn is_imaginary(&self, epsilon: T) -> bool;  // radius² < 0
}
```

### 4. Document Transformation Conventions

Clearly document which transformation is applied first:

```rust
/// Composes two versors: `self` applied first, then `other`.
///
/// `a.compose(&b).transform(x) == b.transform(a.transform(x))`
///
/// In geometric product terms: returns `self * other`.
pub fn compose(&self, other: &Self) -> Self;
```

### 5. Test All Composition Methods Are Equivalent

```rust
#[test]
fn all_composition_methods_equivalent(/* params */) {
    // All four should give the same result:
    // 1. composed_versor.transform(p)
    // 2. v2.transform(v1.transform(p))
    // 3. (na_t1 * na_t2) * na_point
    // 4. na_t1 * (na_t2 * na_point)
}
```

### 6. Wrapper Types for Constrained Elements

```rust
/// A conformal point (null vector in CGA).
/// Invariant: self · self = 0
pub struct Point<T: Float> { /* ... */ }

/// A unit translator versor.
/// Invariant: T T̃ = 1
pub struct UnitTranslator<T: Float>(Translator<T>);

impl<T: Float + Debug> Arbitrary for UnitTranslator<T> { /* ... */ }
```

---

## Background

Conformal Geometric Algebra embeds Euclidean space into a higher-dimensional
space where circles, spheres, planes, and lines are represented uniformly as blades.
All transformations (rotations, translations, dilations, inversions) are versors.

### Signature

| Algebra | Signature | Dimension | Blades |
|---------|-----------|-----------|--------|
| 2D CGA | `Cl(3,1,0)` | 4 | 16 |
| 3D CGA | `Cl(4,1,0)` | 5 | 32 |

### Null Basis Convention

Following the wiki convention:
- `e₊` squares to +1, `e₋` squares to -1
- `e∞ = e₋ + e₊` (point at infinity)
- `e₀ = (e₋ - e₊) / 2` (origin)
- `e∞ · e₀ = -1`, `e∞² = e₀² = 0`

---

## Implementation Strategy

**Incremental approach with verification at each step:**

1. **Phase 1: Generic Multivector with CGA Signatures**
2. **Phase 2: Specialized 2D Types** (Point, Circle, Line, Translator, Rotor)
3. **Phase 3: Specialized 3D Types** (Point, Sphere, Plane, Circle, Line, Translator, Rotor, Dilator)
4. **Phase 4: Geometric Constraints as Invariants**

---

## Phase 1: CGA Signatures

### Deliverables

#### Signature Types (`src/signature/conformal.rs`)

```rust
/// 2D Conformal GA: Cl(3,1,0)
///
/// Embeds 2D Euclidean space into a 4D conformal space.
/// - Grade 1 (vectors): Conformal points (null vectors)
/// - Grade 2 (bivectors): Point pairs, circles, lines
/// - Grade 3 (trivectors): Circles, lines (dual representation)
///
/// Basis: e₁, e₂, e₊, e₋ where e₊² = 1, e₋² = -1
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct Conformal2;

impl Signature for Conformal2 {
    type NumBlades = typenum::U16; // 2⁴ = 16

    const P: usize = 3;  // e₁, e₂, e₊
    const Q: usize = 1;  // e₋
    const R: usize = 0;

    fn metric(i: usize) -> i8 {
        match i {
            0 | 1 | 2 => 1,   // e₁² = e₂² = e₊² = 1
            3 => -1,          // e₋² = -1
            _ => panic!("basis index {i} out of range for Conformal2"),
        }
    }
}

/// 3D Conformal GA: Cl(4,1,0)
///
/// Embeds 3D Euclidean space into a 5D conformal space.
/// - Grade 1: Conformal points
/// - Grade 2: Point pairs
/// - Grade 3: Circles, lines
/// - Grade 4: Spheres, planes
///
/// Basis: e₁, e₂, e₃, e₊, e₋ where e₊² = 1, e₋² = -1
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct Conformal3;

impl Signature for Conformal3 {
    type NumBlades = typenum::U32; // 2⁵ = 32

    const P: usize = 4;  // e₁, e₂, e₃, e₊
    const Q: usize = 1;  // e₋
    const R: usize = 0;

    fn metric(i: usize) -> i8 {
        match i {
            0 | 1 | 2 | 3 => 1,  // e₁² = e₂² = e₃² = e₊² = 1
            4 => -1,             // e₋² = -1
            _ => panic!("basis index {i} out of range for Conformal3"),
        }
    }
}
```

### Phase 1 Testing

```rust
proptest! {
    #[test]
    fn cga3_null_basis_properties(s in -100.0f64..100.0) {
        // e∞ = e₋ + e₊ should square to zero
        let e_inf = e_minus() + e_plus();
        let sq = &e_inf * &e_inf;
        prop_assert!(sq.is_zero(ABS_DIFF_EQ_EPS));

        // e₀ = (e₋ - e₊) / 2 should square to zero
        let e_o = (e_minus() - e_plus()) * 0.5;
        let sq = &e_o * &e_o;
        prop_assert!(sq.is_zero(ABS_DIFF_EQ_EPS));

        // e∞ · e₀ = -1
        let dot = e_inf.inner(&e_o).scalar_part();
        prop_assert!(abs_diff_eq!(dot, -1.0, epsilon = ABS_DIFF_EQ_EPS));
    }

    #[test]
    fn cga3_geometric_product_associative(
        a in any::<Multivector<f64, Conformal3>>(),
        b in any::<Multivector<f64, Conformal3>>(),
        c in any::<Multivector<f64, Conformal3>>(),
    ) {
        let lhs = &(&a * &b) * &c;
        let rhs = &a * &(&b * &c);
        prop_assert!(abs_diff_eq!(lhs, rhs, epsilon = ABS_DIFF_EQ_EPS));
    }
}
```

### Phase 1 Verification

- [ ] `cargo fmt` passes
- [ ] `cargo clippy --all-features` passes
- [ ] `cargo doc --all-features --no-deps` passes
- [ ] `cargo test --all-features` - all property tests pass
- [ ] `cargo deny check` passes

---

## Phase 2: Specialized 2D Types

### Module Structure

```
src/specialized/conformal/
├── mod.rs              # Module root
└── dim2/
    ├── mod.rs          # Public API
    ├── types.rs        # Point, Circle, Line, Translator, Rotor
    ├── ops.rs          # Operations (join, meet, transform)
    ├── conversions.rs  # To/from Multivector
    ├── arbitrary.rs    # Proptest support
    └── nalgebra.rs     # nalgebra conversions (feature-gated)
```

### Type Definitions

```rust
/// A point in 2D CGA (null vector).
///
/// Conformal embedding: P = x·e₁ + y·e₂ + ½(x² + y²)·e∞ + e₀
///
/// Reference: https://conformalgeometricalgebra.org/wiki/index.php?title=Point
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Point<T: Float> {
    pub e1: T,   // x-coordinate
    pub e2: T,   // y-coordinate
    pub ep: T,   // e₊ component
    pub em: T,   // e₋ component
}

/// A circle in 2D CGA (grade-3 trivector or dual sphere).
///
/// Reference: https://conformalgeometricalgebra.org/wiki/index.php?title=Circle
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Circle<T: Float> { /* ... */ }

/// A line in 2D CGA (flat circle through infinity).
///
/// Reference: https://conformalgeometricalgebra.org/wiki/index.php?title=Line
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Line<T: Float> { /* ... */ }

/// A translator in 2D CGA.
///
/// T = 1 - (d/2)·e∞ where d is the translation vector.
///
/// Reference: https://conformalgeometricalgebra.org/wiki/index.php?title=Translator
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Translator<T: Float> { /* ... */ }

/// A rotor in 2D CGA (rotation around a point).
///
/// Reference: https://conformalgeometricalgebra.org/wiki/index.php?title=Motor
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Rotor<T: Float> { /* ... */ }
```

### SymPy Derivations Required

Create `derivations/src/clifford_derivations/cga.py`:

```python
"""CGA derivations using SymPy.

Reference: https://conformalgeometricalgebra.org/wiki/index.php

Usage:
    cd derivations
    uv run python -m clifford_derivations.cga
"""

@with_timeout(300)
def derive_point_embedding():
    """Derive P = x·e₁ + y·e₂ + ½(x² + y²)·e∞ + e₀"""
    pass

@with_timeout(300)
def derive_translator_sandwich():
    """Derive T * P * T̃ for translator versor."""
    pass

@with_timeout(300)
def derive_circle_from_three_points():
    """Derive circle = P₁ ∧ P₂ ∧ P₃"""
    pass
```

### Phase 2 Verification

- [ ] `cargo fmt` passes
- [ ] `cargo clippy --all-features` passes
- [ ] `cargo doc --all-features --no-deps` passes
- [ ] `cargo test --all-features` - all property tests pass
- [ ] `cargo deny check` passes
- [ ] Operations verified against SymPy derivations
- [ ] Operations match nalgebra equivalents

---

## Phase 3: Specialized 3D Types

### Module Structure

```
src/specialized/conformal/dim3/
├── mod.rs          # Public API
├── types.rs        # Point, Sphere, Plane, Circle, Line, Translator, Rotor, Dilator
├── ops.rs          # Operations
├── conversions.rs  # To/from Multivector
├── arbitrary.rs    # Proptest support
└── nalgebra.rs     # nalgebra conversions
```

### Type Definitions

```rust
/// A point in 3D CGA (null vector).
///
/// Conformal embedding: P = x·e₁ + y·e₂ + z·e₃ + ½(x² + y² + z²)·e∞ + e₀
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Point<T: Float> {
    pub e1: T, pub e2: T, pub e3: T,
    pub ep: T, pub em: T,
}

/// A sphere in 3D CGA (grade-4 quadvector).
///
/// Dual sphere: S = c - ½r²·e∞ where c is center point, r is radius.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Sphere<T: Float> { /* ... */ }

/// A plane in 3D CGA (flat sphere through infinity).
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Plane<T: Float> { /* ... */ }

/// A circle in 3D CGA (intersection of two spheres).
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Circle<T: Float> { /* ... */ }

/// A line in 3D CGA (intersection of two planes).
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Line<T: Float> { /* ... */ }

/// A translator in 3D CGA.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Translator<T: Float> { /* ... */ }

/// A rotor in 3D CGA (rotation).
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Rotor<T: Float> { /* ... */ }

/// A dilator in 3D CGA (uniform scaling from a point).
///
/// D = cosh(λ/2) + sinh(λ/2)·e∞·e₀ for scale factor e^λ.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Dilator<T: Float> { /* ... */ }

/// An inversor in 3D CGA (inversion in a sphere).
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Inversor<T: Float> { /* ... */ }
```

### Key Operations

```rust
impl<T: Float> Point<T> {
    /// Distance to another point.
    pub fn distance(&self, other: &Point<T>) -> T {
        // d² = -2 · (P₁ · P₂)
    }
}

impl<T: Float> Sphere<T> {
    /// Meet of two spheres: their intersection circle.
    pub fn meet(&self, other: &Sphere<T>) -> Circle<T>;

    /// Check if a point lies on this sphere.
    pub fn contains(&self, p: &Point<T>, epsilon: T) -> bool {
        // P · S = 0 iff P on S
    }
}

impl<T: Float> Translator<T> {
    /// Compose two translators.
    pub fn compose(&self, other: &Self) -> Self;

    /// Transform a point.
    pub fn transform_point(&self, p: &Point<T>) -> Point<T>;

    /// Transform a sphere.
    pub fn transform_sphere(&self, s: &Sphere<T>) -> Sphere<T>;
}

impl<T: Float> Dilator<T> {
    /// Uniform scaling centered at origin.
    pub fn from_scale(factor: T) -> Self;

    /// Uniform scaling centered at a point.
    pub fn from_center_scale(center: &Point<T>, factor: T) -> Self;
}
```

### Phase 3 Verification

- [ ] `cargo fmt` passes
- [ ] `cargo clippy --all-features` passes
- [ ] `cargo doc --all-features --no-deps` passes
- [ ] `cargo test --all-features` - all property tests pass
- [ ] `cargo deny check` passes
- [ ] Operations verified against SymPy derivations
- [ ] Operations match nalgebra equivalents

---

## Phase 4: Geometric Constraints

### Constraint Verification Methods

```rust
impl<T: Float> Point<T> {
    /// Returns true if this is a valid conformal point (null vector).
    pub fn is_null(&self, epsilon: T) -> bool {
        self.norm_squared().abs() < epsilon
    }
}

impl<T: Float> Sphere<T> {
    /// Returns true if the sphere has real (positive) radius.
    pub fn is_real(&self, epsilon: T) -> bool;

    /// Returns true if the sphere has imaginary radius (point pair).
    pub fn is_imaginary(&self, epsilon: T) -> bool;
}

impl<T: Float> Circle<T> {
    /// Returns true if this is a real circle (not a point pair).
    pub fn is_real(&self, epsilon: T) -> bool;
}
```

### Property Tests for Invariants

```rust
proptest! {
    #[test]
    fn point_is_null(x in -100.0f64..100.0, y in -100.0f64..100.0, z in -100.0f64..100.0) {
        let p = Point::new(x, y, z);
        prop_assert!(p.is_null(ABS_DIFF_EQ_EPS));
    }

    #[test]
    fn translator_preserves_null(
        t in any::<Translator<f64>>(),
        p in any::<Point<f64>>(),
    ) {
        let transformed = t.transform_point(&p);
        prop_assert!(transformed.is_null(ABS_DIFF_EQ_EPS));
    }

    #[test]
    fn sphere_contains_its_points(
        cx in -10.0f64..10.0, cy in -10.0f64..10.0, cz in -10.0f64..10.0,
        r in 0.1f64..10.0,
        theta in 0.0f64..std::f64::consts::TAU,
        phi in 0.0f64..std::f64::consts::PI,
    ) {
        let sphere = Sphere::from_center_radius(cx, cy, cz, r);
        let point = Point::new(
            cx + r * phi.sin() * theta.cos(),
            cy + r * phi.sin() * theta.sin(),
            cz + r * phi.cos(),
        );
        prop_assert!(sphere.contains(&point, ABS_DIFF_EQ_EPS));
    }
}
```

### Phase 4 Verification

- [ ] All constraint methods implemented
- [ ] Property tests verify constraints preserved by operations
- [ ] Documentation references wiki geometric constraint page

---

## nalgebra Interoperability

### Conversions (feature-gated)

```rust
// Point conversions
impl<T: Float + na::Scalar> From<na::Point3<T>> for Point<T>;
impl<T: Float + na::Scalar> TryFrom<Point<T>> for na::Point3<T>;

// Sphere conversions (center + radius)
impl<T: Float + na::Scalar> From<Sphere<T>> for (na::Point3<T>, T);
impl<T: Float + na::Scalar> From<(na::Point3<T>, T)> for Sphere<T>;

// Translator conversions
impl<T: Float + na::RealField> From<na::Translation3<T>> for Translator<T>;
impl<T: Float + na::RealField> From<Translator<T>> for na::Translation3<T>;

// Rotor conversions
impl<T: Float + na::RealField> From<na::UnitQuaternion<T>> for Rotor<T>;
impl<T: Float + na::RealField> From<Rotor<T>> for na::UnitQuaternion<T>;
```

### Error Types

```rust
#[derive(Clone, Copy, Debug, PartialEq, Eq, Error)]
pub enum PointConversionError {
    #[error("conformal point is not a valid Euclidean point")]
    InvalidConformalPoint,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq, Error)]
pub enum SphereConversionError {
    #[error("sphere has imaginary radius")]
    ImaginaryRadius,
    #[error("sphere is degenerate (zero radius)")]
    DegenerateSphere,
}
```

---

## Benchmarks

### Key Benchmarks

```rust
// Point embedding/extraction
fn bench_point_embed(c: &mut Criterion);
fn bench_point_extract(c: &mut Criterion);

// Versor transformations
fn bench_translator_transform_point(c: &mut Criterion);
fn bench_dilator_transform_point(c: &mut Criterion);
fn bench_rotor_transform_point(c: &mut Criterion);

// Geometric operations
fn bench_sphere_meet_sphere(c: &mut Criterion);  // circle
fn bench_plane_meet_sphere(c: &mut Criterion);   // circle
fn bench_three_points_circle(c: &mut Criterion); // P₁ ∧ P₂ ∧ P₃

// Comparison with nalgebra
fn bench_translate_cga_vs_nalgebra(c: &mut Criterion);
fn bench_sphere_intersect_cga_vs_algebraic(c: &mut Criterion);
```

### Expected Performance Targets

| Operation | Target | Notes |
|-----------|--------|-------|
| Point embed | < 10 ns | x² + y² + z² computation |
| Point extract | < 5 ns | Division by weight |
| Translator transform | < 30 ns | Sandwich product |
| Sphere meet | < 20 ns | Outer product |
| CGA vs nalgebra | Comparable | CGA has embed/extract overhead |

---

## SymPy Derivation Package

Add to `derivations/src/clifford_derivations/`:

```
cga.py              # CGA-specific derivations
├── derive_point_embedding()
├── derive_euclidean_extraction()
├── derive_translator_sandwich()
├── derive_dilator_sandwich()
├── derive_rotor_sandwich()
├── derive_sphere_from_center_radius()
├── derive_plane_from_normal_distance()
├── derive_circle_from_three_points()
├── verify_null_constraint()
└── generate_rust_*()  # Rust code generation using rust_code()
```

---

## Summary

| Phase | Scope | Key Deliverables | Status |
|-------|-------|------------------|--------|
| 1 | Generic CGA | Signatures, null basis tests | Pending |
| 2 | 2D Specialized | Point, Circle, Line, Translator, Rotor | Pending |
| 3 | 3D Specialized | Point, Sphere, Plane, Circle, Line, Translator, Rotor, Dilator | Pending |
| 4 | Geometric Constraints | Null constraint, real/imaginary detection | Pending |

---

## Final Verification

- [ ] `cargo fmt` passes
- [ ] `cargo clippy --all-features` passes
- [ ] `cargo doc --all-features --no-deps` passes
- [ ] `cargo test --all-features` - all property tests pass
- [ ] `cargo deny check` passes
- [ ] All SymPy derivations verified
- [ ] nalgebra conversions tested with feature flags
- [ ] Benchmarks run and documented
- [ ] Wiki references in rustdoc

## Release Checklist

- [ ] Version bump in Cargo.toml
- [ ] CHANGELOG updated
- [ ] All PRDs marked complete
- [ ] `CARGO_REGISTRY_TOKEN` set
- [ ] Tag and push: `git tag v0.1.0 && git push --tags`
