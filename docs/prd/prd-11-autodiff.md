# PRD-11: Automatic Differentiation via Dual Numbers

**Status**: Pending
**Goal**: Enable automatic differentiation for all GA operations using dual numbers

## Motivation

Automatic differentiation (autodiff) is essential for:

1. **Optimization**: Finding optimal rotations, poses, or configurations
2. **Physics simulation**: Computing forces from potential energy functions
3. **Machine learning**: Backpropagation through GA operations
4. **Sensitivity analysis**: Understanding how outputs change with inputs
5. **Robotics**: Jacobians for inverse kinematics

Dual numbers provide forward-mode autodiff with minimal implementation complexity and excellent performance for low-dimensional derivatives.

## Background: Dual Numbers

A dual number has the form `a + bε` where `ε² = 0` (nilpotent). This algebraic property means:

```
f(a + bε) = f(a) + f'(a)·b·ε
```

For any analytic function, evaluating at a dual number automatically computes the derivative:

| Operation | Result |
|-----------|--------|
| `(a + bε) + (c + dε)` | `(a+c) + (b+d)ε` |
| `(a + bε) × (c + dε)` | `ac + (ad + bc)ε` |
| `sin(a + bε)` | `sin(a) + cos(a)·b·ε` |
| `exp(a + bε)` | `exp(a) + exp(a)·b·ε` |
| `sqrt(a + bε)` | `sqrt(a) + b/(2·sqrt(a))·ε` |

## Design

### 1. Dual Number Type

```rust
/// A dual number for forward-mode automatic differentiation.
///
/// Represents `real + dual·ε` where `ε² = 0`.
///
/// # Example
///
/// ```
/// use clifford::autodiff::Dual;
///
/// // Compute derivative of x² at x=3
/// let x = Dual::new(3.0, 1.0);  // x = 3 + 1ε (seed derivative = 1)
/// let y = x * x;                 // y = 9 + 6ε
///
/// assert_eq!(y.real(), 9.0);     // f(3) = 9
/// assert_eq!(y.dual(), 6.0);     // f'(3) = 6
/// ```
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Dual<T> {
    real: T,
    dual: T,
}

impl<T> Dual<T> {
    /// Create a new dual number.
    pub fn new(real: T, dual: T) -> Self;

    /// Create a constant (derivative = 0).
    pub fn constant(real: T) -> Self where T: num_traits::Zero;

    /// Create a variable (derivative = 1).
    pub fn variable(real: T) -> Self where T: num_traits::One;

    /// Get the real (primal) part.
    pub fn real(&self) -> T;

    /// Get the dual (derivative) part.
    pub fn dual(&self) -> T;
}
```

### 2. Implement Float for Dual\<T\>

The key insight: if `Dual<T>` implements our `Float` trait, then **all existing GA operations automatically support autodiff** with zero additional code.

```rust
impl<T: Float> num_traits::Float for Dual<T> {
    fn sin(self) -> Self {
        Dual::new(
            self.real.sin(),
            self.dual * self.real.cos(),
        )
    }

    fn cos(self) -> Self {
        Dual::new(
            self.real.cos(),
            -self.dual * self.real.sin(),
        )
    }

    fn sqrt(self) -> Self {
        let sqrt_real = self.real.sqrt();
        Dual::new(
            sqrt_real,
            self.dual / (T::from_i8(2) * sqrt_real),
        )
    }

    fn exp(self) -> Self {
        let exp_real = self.real.exp();
        Dual::new(exp_real, self.dual * exp_real)
    }

    fn ln(self) -> Self {
        Dual::new(
            self.real.ln(),
            self.dual / self.real,
        )
    }

    // ... all other Float methods
}

impl<T: Float> Float for Dual<T> {
    // Note: TWO and PI are associated constants on Float trait.
    // For Dual<T>, these must be defined carefully since T::zero() is a method.
    // The actual implementation will need to handle this appropriately,
    // potentially by making TWO/PI methods instead of constants for Dual.
    const TWO: Self = Dual { real: T::TWO, dual: T::ZERO };  // T::ZERO if available
    const PI: Self = Dual { real: T::PI, dual: T::ZERO };

    fn from_f64(value: f64) -> Self {
        Dual::constant(T::from_f64(value))
    }
    // ...
}

// Note: The above const definitions are illustrative. In practice, since T::zero()
// is a method (not a const), implementing Float for Dual<T> may require either:
// 1. Changing TWO/PI to methods in the Float trait, or
// 2. Using a different approach for the dual part initialization
```

### 3. Approx Traits for Dual Numbers

For testing and comparisons:

```rust
impl<T: Float + AbsDiffEq> AbsDiffEq for Dual<T> {
    type Epsilon = T::Epsilon;

    fn default_epsilon() -> Self::Epsilon {
        T::default_epsilon()
    }

    fn abs_diff_eq(&self, other: &Self, epsilon: Self::Epsilon) -> bool {
        self.real.abs_diff_eq(&other.real, epsilon) &&
        self.dual.abs_diff_eq(&other.dual, epsilon)
    }
}
```

### 4. Usage with GA Types

Once `Dual<f64>` implements `Float`, it works everywhere:

```rust
use clifford::autodiff::Dual;
use clifford::specialized::euclidean::dim3::{Vector, Rotor, Bivector};

// Derivative of rotation angle on output vector
fn rotate_x_component_derivative() {
    let angle = Dual::variable(0.5);  // θ = 0.5, dθ = 1
    let plane = Bivector::<Dual<f64>>::unit_xy();

    let rotor = Rotor::from_angle_plane(angle, plane);
    let v = Vector::new(
        Dual::constant(1.0),
        Dual::constant(0.0),
        Dual::constant(0.0),
    );

    let rotated = rotor.rotate(v);

    // rotated.x.real() = cos(0.5) ≈ 0.877
    // rotated.x.dual() = -sin(0.5) ≈ -0.479 (derivative w.r.t. angle)
}

// Gradient of distance function
fn distance_gradient() {
    // Compute ∂dist/∂x, ∂dist/∂y, ∂dist/∂z
    let px = Dual::new(1.0, 1.0);  // seed x
    let py = Dual::new(2.0, 0.0);
    let pz = Dual::new(3.0, 0.0);

    let p = Vector::new(px, py, pz);
    let dist = p.norm();

    // dist.dual() = ∂||p||/∂x = x/||p||
}
```

### 5. Jacobian Helper

For computing full Jacobians:

```rust
/// Compute the Jacobian of a function f: R^n -> R^m
pub fn jacobian<F, const N: usize, const M: usize>(
    f: F,
    x: [f64; N],
) -> [[f64; N]; M]
where
    F: Fn([Dual<f64>; N]) -> [Dual<f64>; M],
{
    let mut jac = [[0.0; N]; M];

    for i in 0..N {
        // Seed the i-th input
        let mut x_dual = x.map(Dual::constant);
        x_dual[i] = Dual::variable(x[i]);

        let y = f(x_dual);

        for j in 0..M {
            jac[j][i] = y[j].dual();
        }
    }

    jac
}
```

### 6. Multi-Dual Numbers (Gradient in One Pass)

For computing gradients ∇f = (∂f/∂x₁, ..., ∂f/∂xₙ) in a single forward pass:

```rust
/// Multi-dual number with N independent infinitesimals.
///
/// Represents `real + g₁ε₁ + g₂ε₂ + ... + gₙεₙ` where `εᵢεⱼ = 0` for all i,j.
///
/// This computes f(x) and all N partial derivatives in ONE forward pass,
/// vs N passes with standard dual numbers.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct MultiDual<T, const N: usize> {
    real: T,
    grad: [T; N],  // [∂f/∂x₁, ∂f/∂x₂, ..., ∂f/∂xₙ]
}

impl<T, const N: usize> MultiDual<T, N> {
    /// Create a constant (all derivatives = 0).
    pub fn constant(real: T) -> Self where T: num_traits::Zero + Copy;

    /// Create the i-th variable (∂xᵢ/∂xᵢ = 1, others = 0).
    pub fn variable(real: T, index: usize) -> Self where T: num_traits::Zero + num_traits::One + Copy;

    /// Get the real (primal) value.
    pub fn real(&self) -> T;

    /// Get the full gradient vector.
    pub fn grad(&self) -> [T; N];

    /// Get a specific partial derivative.
    pub fn partial(&self, index: usize) -> T;
}
```

**Arithmetic rules** (all εᵢεⱼ = 0):

```rust
impl<T: Float, const N: usize> Mul for MultiDual<T, N> {
    fn mul(self, rhs: Self) -> Self {
        // (a + Σaᵢεᵢ)(b + Σbᵢεᵢ) = ab + Σ(a·bᵢ + b·aᵢ)εᵢ
        let mut grad = [T::zero(); N];
        for i in 0..N {
            grad[i] = self.real * rhs.grad[i] + rhs.real * self.grad[i];
        }
        MultiDual { real: self.real * rhs.real, grad }
    }
}

impl<T: Float, const N: usize> num_traits::Float for MultiDual<T, N> {
    fn sin(self) -> Self {
        let cos_real = self.real.cos();
        MultiDual {
            real: self.real.sin(),
            grad: self.grad.map(|g| g * cos_real),
        }
    }
    // ... other functions apply chain rule to each component
}
```

**Usage:**

```rust
use std::f64::consts::FRAC_PI_2;
use clifford::autodiff::MultiDual;

// Compute gradient of f(x,y,z) = x²y + sin(z) at (1, 2, π/2)
fn f<T: Float>(x: T, y: T, z: T) -> T {
    x * x * y + z.sin()
}

let x = MultiDual::<f64, 3>::variable(1.0, 0);  // ∂/∂x
let y = MultiDual::<f64, 3>::variable(2.0, 1);  // ∂/∂y
let z = MultiDual::<f64, 3>::variable(FRAC_PI_2, 2);  // ∂/∂z

let result = f(x, y, z);

assert_eq!(result.real(), 3.0);           // f(1,2,π/2) = 1²·2 + 1 = 3
assert_eq!(result.grad(), [4.0, 1.0, 0.0]); // ∇f = [2xy, x², cos(z)] = [4, 1, 0]
```

### 7. Hyper-Dual Numbers (Efficient Second Derivatives)

Nested dual numbers `Dual<Dual<f64>>` work but waste computation (redundant cross-terms). Hyper-dual numbers are optimized for second derivatives:

```rust
/// Hyper-dual number for computing f, ∇f, and ∂²f/∂x∂y efficiently.
///
/// Represents `a + b₁ε₁ + b₂ε₂ + cε₁ε₂` where:
/// - `ε₁² = ε₂² = 0` (nilpotent)
/// - `ε₁ε₂ ≠ 0` (captures mixed partial)
///
/// In one evaluation, computes:
/// - f(x,y)
/// - ∂f/∂x
/// - ∂f/∂y
/// - ∂²f/∂x∂y (mixed partial / Hessian off-diagonal)
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct HyperDual<T> {
    real: T,    // f(x,y)
    eps1: T,    // ∂f/∂x
    eps2: T,    // ∂f/∂y
    eps12: T,   // ∂²f/∂x∂y
}

impl<T> HyperDual<T> {
    pub fn new(real: T, eps1: T, eps2: T, eps12: T) -> Self;

    /// Create a constant.
    pub fn constant(real: T) -> Self where T: num_traits::Zero;

    /// Create x variable: x + 1·ε₁ + 0·ε₂ + 0·ε₁ε₂
    pub fn var_x(real: T) -> Self where T: num_traits::Zero + num_traits::One;

    /// Create y variable: y + 0·ε₁ + 1·ε₂ + 0·ε₁ε₂
    pub fn var_y(real: T) -> Self where T: num_traits::Zero + num_traits::One;

    pub fn real(&self) -> T;
    pub fn df_dx(&self) -> T;      // ∂f/∂x
    pub fn df_dy(&self) -> T;      // ∂f/∂y
    pub fn d2f_dxdy(&self) -> T;   // ∂²f/∂x∂y
}
```

**Arithmetic rules:**

```rust
impl<T: Float> Mul for HyperDual<T> {
    fn mul(self, rhs: Self) -> Self {
        // (a + b₁ε₁ + b₂ε₂ + cε₁ε₂)(d + e₁ε₁ + e₂ε₂ + fε₁ε₂)
        // = ad + (ae₁ + b₁d)ε₁ + (ae₂ + b₂d)ε₂ + (af + cd + b₁e₂ + b₂e₁)ε₁ε₂
        HyperDual {
            real: self.real * rhs.real,
            eps1: self.real * rhs.eps1 + self.eps1 * rhs.real,
            eps2: self.real * rhs.eps2 + self.eps2 * rhs.real,
            eps12: self.real * rhs.eps12 + self.eps12 * rhs.real
                 + self.eps1 * rhs.eps2 + self.eps2 * rhs.eps1,
        }
    }
}

impl<T: Float> num_traits::Float for HyperDual<T> {
    fn sin(self) -> Self {
        let (s, c) = (self.real.sin(), self.real.cos());
        HyperDual {
            real: s,
            eps1: self.eps1 * c,
            eps2: self.eps2 * c,
            eps12: self.eps12 * c - self.eps1 * self.eps2 * s,
        }
    }
    // ...
}
```

**Usage:**

```rust
use clifford::autodiff::HyperDual;

// Compute f, ∂f/∂x, ∂f/∂y, ∂²f/∂x∂y for f(x,y) = x²y
let x = HyperDual::var_x(3.0);  // x = 3
let y = HyperDual::var_y(2.0);  // y = 2

let f = x * x * y;

assert_eq!(f.real(), 18.0);      // f(3,2) = 9·2 = 18
assert_eq!(f.df_dx(), 12.0);     // ∂f/∂x = 2xy = 12
assert_eq!(f.df_dy(), 9.0);      // ∂f/∂y = x² = 9
assert_eq!(f.d2f_dxdy(), 6.0);   // ∂²f/∂x∂y = 2x = 6
```

### 8. Generalized Hyper-Dual (Full Hessian)

For the complete Hessian matrix of f: Rⁿ → R:

```rust
/// Hyper-dual with N variables, computing full Hessian.
///
/// Stores: real + Σgᵢεᵢ + Σᵢ≤ⱼ hᵢⱼεᵢεⱼ
///
/// Components:
/// - real: f(x)
/// - grad[i]: ∂f/∂xᵢ
/// - hess[i][j]: ∂²f/∂xᵢ∂xⱼ (symmetric, stored as upper triangle)
#[derive(Debug, Clone, PartialEq)]
pub struct HyperDualN<T, const N: usize> {
    real: T,
    grad: [T; N],
    hess: [[T; N]; N],  // Symmetric, but store full for simplicity
}

impl<T, const N: usize> HyperDualN<T, N> {
    /// Create the i-th variable.
    pub fn variable(real: T, index: usize) -> Self;

    /// Get the gradient vector.
    pub fn gradient(&self) -> [T; N];

    /// Get the full Hessian matrix.
    pub fn hessian(&self) -> [[T; N]; N];
}
```

**Usage for Newton's method optimization:**

```rust
use clifford::autodiff::HyperDualN;

fn cost_function<T: Float>(params: [T; 3]) -> T {
    // Some cost to minimize
    let [x, y, z] = params;
    (x - T::one()).powi(2) + (y - T::from_i8(2)).powi(2) + (z + T::one()).powi(2)
}

// Compute gradient and Hessian at current point
let x = [
    HyperDualN::<f64, 3>::variable(0.5, 0),
    HyperDualN::<f64, 3>::variable(1.0, 1),
    HyperDualN::<f64, 3>::variable(0.0, 2),
];

let cost = cost_function(x);

let grad = cost.gradient();   // ∇f for gradient descent
let hess = cost.hessian();    // H for Newton's method: x_new = x - H⁻¹∇f
```

### 9. Nested Dual Numbers (Alternative)

For completeness, nested `Dual<Dual<T>>` also works but is less efficient:

```rust
// Second derivative via nesting (works but redundant)
type Dual2<T> = Dual<Dual<T>>;

let x: Dual2<f64> = Dual::new(
    Dual::new(2.0, 1.0),  // x = 2, dx = 1
    Dual::new(1.0, 0.0),  // for second derivative
);

let y = x.sin();

// y.real().real() = sin(2)
// y.real().dual() = cos(2)      (first derivative)
// y.dual().dual() = -sin(2)     (second derivative)
// y.dual().real() = cos(2)      (redundant, same as y.real().dual())
```

**Tradeoff**: Nested duals are simpler to implement but use 4 values instead of 3 for second derivatives, and 8 instead of 4 for third derivatives.

## Comparison of Approaches

| Type | Memory | One-Pass Computes | Best For |
|------|--------|-------------------|----------|
| `Dual<T>` | 2×T | f, ∂f/∂xᵢ (one variable) | Single derivatives, simple cases |
| `MultiDual<T, N>` | (N+1)×T | f, ∇f (all N partials) | Full gradient, small N |
| `HyperDual<T>` | 4×T | f, ∂f/∂x, ∂f/∂y, ∂²f/∂x∂y | Mixed partials, 2 variables |
| `HyperDualN<T, N>` | (1+N+N²)×T | f, ∇f, H (full Hessian) | Newton optimization |
| `Dual<Dual<T>>` | 4×T | f, f', f'' (single var) | Higher-order, single variable |

### When to Use What

1. **Single derivative** (∂f/∂x): Use `Dual<T>`
2. **Gradient** (∇f): Use `MultiDual<T, N>` for N ≤ ~8, else multiple `Dual` passes
3. **Single second derivative** (f''): Use `Dual<Dual<T>>` or `HyperDual` with x=y
4. **Mixed partial** (∂²f/∂x∂y): Use `HyperDual<T>`
5. **Full Hessian**: Use `HyperDualN<T, N>` for small N, else finite differences on gradient
6. **Many parameters** (N > 10-20): Consider reverse-mode AD (separate PRD)

## Module Structure

```
src/
  autodiff/
    mod.rs           # Module root, re-exports
    dual.rs          # Dual<T> - basic forward-mode AD
    multi_dual.rs    # MultiDual<T, N> - gradient in one pass
    hyper_dual.rs    # HyperDual<T> - second derivatives (2 vars)
    hyper_dual_n.rs  # HyperDualN<T, N> - full Hessian
    jacobian.rs      # Jacobian/Hessian helper functions
    traits.rs        # Common traits for AD types
```

## Feature Flag

```toml
[features]
autodiff = []  # No additional dependencies needed!
```

## Deliverables

### Phase 1: Basic Dual Numbers
1. **`Dual<T>` type** with all arithmetic operations
2. **`num_traits::Float` implementation** for `Dual<T>`
3. **`Float` trait implementation** for `Dual<T>`
4. **`approx` trait implementations** for testing
5. **Jacobian helper functions**
6. **Documentation with examples**

### Phase 2: Multi-Dual and Hyper-Dual
7. **`MultiDual<T, N>` type** for gradients in one pass
8. **`HyperDual<T>` type** for efficient mixed partials
9. **`HyperDualN<T, N>` type** for full Hessian computation
10. **Hessian helper functions**
11. **Property-based tests** verifying all derivative types

## Verification

### Phase 1
- [ ] `Dual<f64>` implements `Float`
- [ ] All trig functions have correct derivatives
- [ ] All GA operations work with `Dual<f64>` scalars
- [ ] Rotor derivatives are correct (verified against numerical diff)
- [ ] Jacobian helper produces correct results

### Phase 2
- [ ] `MultiDual<f64, N>` computes correct gradients
- [ ] `HyperDual<f64>` computes correct mixed partials
- [ ] `HyperDualN<f64, N>` computes correct full Hessian
- [ ] Results match nested dual numbers (for correctness check)
- [ ] Performance benchmarks show expected speedups

## Example Use Cases

### 1. Optimal Rotation Finding

```rust
// Find rotation that aligns v1 to v2
fn alignment_cost(angle: Dual<f64>) -> Dual<f64> {
    let rotor = Rotor::from_angle_plane(angle, plane);
    let rotated = rotor.rotate(v1);
    (rotated - v2).norm_squared()
}

// Gradient descent
let mut angle = 0.0;
for _ in 0..100 {
    let cost = alignment_cost(Dual::variable(angle));
    angle -= 0.1 * cost.dual();  // gradient step
}
```

### 2. Inverse Kinematics

```rust
// Jacobian of end-effector position w.r.t. joint angles
fn forward_kinematics(joints: [Dual<f64>; 3]) -> Vector<Dual<f64>> {
    let r1 = Rotor::from_angle_plane(joints[0], Bivector::unit_xy());
    let r2 = Rotor::from_angle_plane(joints[1], Bivector::unit_xz());
    let r3 = Rotor::from_angle_plane(joints[2], Bivector::unit_yz());

    let rotor = r1.compose(r2).compose(r3);
    // Use scale() for generic scalar types like Dual<f64>
    rotor.rotate(Vector::unit_x().scale(arm_length))
}

let jac = jacobian(forward_kinematics, current_joints);
```

### 3. Physics: Force from Potential

```rust
// Force = -∇U (negative gradient of potential)
fn potential(pos: Vector<Dual<f64>>) -> Dual<f64> {
    // Spring potential: U = k/2 * ||x - x0||²
    let displacement = pos - equilibrium;
    Dual::constant(0.5 * k) * displacement.norm_squared()
}

// Compute force components
let fx = -potential(Vector::new(Dual::variable(x), Dual::constant(y), Dual::constant(z))).dual();
let fy = -potential(Vector::new(Dual::constant(x), Dual::variable(y), Dual::constant(z))).dual();
let fz = -potential(Vector::new(Dual::constant(x), Dual::constant(y), Dual::variable(z))).dual();
```

### 4. Hessian of Rotor Cost Function (HyperDual)

```rust
use clifford::autodiff::HyperDual;
use clifford::specialized::euclidean::dim3::{Vector, Rotor, Bivector};

// Cost: squared distance after rotation
// Find ∂²cost/∂θ₁∂θ₂ for two rotation angles
fn rotation_cost(theta1: HyperDual<f64>, theta2: HyperDual<f64>) -> HyperDual<f64> {
    let r1 = Rotor::from_angle_plane(theta1, Bivector::<HyperDual<f64>>::unit_xy());
    let r2 = Rotor::from_angle_plane(theta2, Bivector::<HyperDual<f64>>::unit_xz());

    let combined = r1.compose(r2);
    let v = Vector::new(
        HyperDual::constant(1.0),
        HyperDual::constant(0.0),
        HyperDual::constant(0.0),
    );
    let target = Vector::new(
        HyperDual::constant(0.0),
        HyperDual::constant(1.0),
        HyperDual::constant(0.0),
    );

    let rotated = combined.rotate(v);
    (rotated - target).norm_squared()
}

let theta1 = HyperDual::var_x(0.5);
let theta2 = HyperDual::var_y(0.3);

let cost = rotation_cost(theta1, theta2);

// Get all second-order info for Newton optimization
let f = cost.real();           // cost value
let df_dtheta1 = cost.df_dx(); // ∂cost/∂θ₁
let df_dtheta2 = cost.df_dy(); // ∂cost/∂θ₂
let d2f = cost.d2f_dxdy();     // ∂²cost/∂θ₁∂θ₂ (Hessian off-diagonal)
```

### 5. Full Gradient with MultiDual

```rust
use clifford::autodiff::MultiDual;

// Compute full gradient of a 3-joint robot arm cost in one pass
fn arm_cost<T: Float>(joints: [T; 3], target: [T; 3]) -> T {
    let r1 = Rotor::from_angle_plane(joints[0], Bivector::<T>::unit_xy());
    let r2 = Rotor::from_angle_plane(joints[1], Bivector::<T>::unit_xz());
    let r3 = Rotor::from_angle_plane(joints[2], Bivector::<T>::unit_yz());

    let end_effector = r1.compose(r2).compose(r3)
        .rotate(Vector::new(T::one(), T::zero(), T::zero()));

    let target_vec = Vector::new(target[0], target[1], target[2]);
    (end_effector - target_vec).norm_squared()
}

// Create variables for all 3 joints
let joints = [
    MultiDual::<f64, 3>::variable(0.1, 0),
    MultiDual::<f64, 3>::variable(0.2, 1),
    MultiDual::<f64, 3>::variable(0.3, 2),
];
let target = [0.5, 0.5, 0.5].map(MultiDual::constant);

let cost = arm_cost(joints, target);

// Full gradient in ONE forward pass!
let gradient = cost.grad();  // [∂cost/∂θ₁, ∂cost/∂θ₂, ∂cost/∂θ₃]
```

## API Considerations for Autodiff Support

For autodiff types to work seamlessly with existing GA types, some API adjustments may be needed:

### 1. Scalar Multiplication Order

Currently, `Vector<T> * T` works for any `T: Float`, but `T * Vector<T>` only works for `f32` and `f64` due to Rust's orphan rules. For `Dual<f64> * Vector<Dual<f64>>`, users must write `vector * scalar` not `scalar * vector`.

**Solution:** The `scale(scalar: T) -> Self` method is now available on all Vector types:

```rust
// These are equivalent:
let v1 = vector * scalar;      // Works for any T: Float
let v2 = vector.scale(scalar); // Also works, more discoverable
```

### 2. Float Trait Constants

The `Float` trait has associated constants `TWO` and `PI`. For `Dual<T>`, these require the dual part to be zero, but `T::zero()` is a method, not a const.

**Options:**
- Change `TWO`/`PI` from `const` to methods: `fn two() -> Self`
- Add a `ZERO` associated constant to `Float` trait
- Use a macro or const fn approach (if Rust const generics allow)

### 3. Vector Field Access

Vectors have public fields (`v.x`, `v.y`, `v.z`), which works well with autodiff since you can write `v.x.real()` and `v.x.dual()` to extract value and derivative.

## Performance Considerations

| Type | Memory per Value | Mul Cost | Best When |
|------|------------------|----------|-----------|
| `Dual<T>` | 2×T | 4 muls, 2 adds | Few derivatives needed |
| `MultiDual<T, N>` | (N+1)×T | 2N muls, N adds | N ≤ ~8 |
| `HyperDual<T>` | 4×T | 10 muls, 6 adds | Mixed partials needed |
| `HyperDualN<T, N>` | O(N²)×T | O(N²) | Full Hessian, small N |

**Key points:**
- Forward-mode is efficient when **outputs >> inputs** (few parameters, many outputs)
- All types are `Copy` when `T: Copy` - no heap allocation
- For N > ~10-20 parameters, forward-mode becomes expensive; consider reverse-mode
- Hyper-dual types have higher per-operation cost but compute more derivatives per pass

## Future Considerations

- **Reverse-mode AD**: For ML/optimization with many parameters (requires tape-based approach, separate PRD)
- **Sparse Jacobians/Hessians**: Exploit sparsity structure in large systems
- **Integration with nalgebra**: Use `Dual<f64>` as nalgebra scalar type
- **GPU-friendly variants**: SIMD-optimized multi-dual for batch operations
- **Third-order derivatives**: `HyperHyperDual` or nested hyper-duals for f'''
- **Interval arithmetic combination**: Dual + interval for verified derivatives with bounds
