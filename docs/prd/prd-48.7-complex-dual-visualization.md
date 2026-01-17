# PRD-48.7: Complex and Dual Number Visualization

**Status**: Draft
**Parent**: PRD-48
**Depends on**: PRD-48.1
**Goal**: Interactive demos for complex numbers and dual numbers (automatic differentiation)

## Overview

Visualize complex number operations with domain coloring, and dual numbers with tangent line computation for automatic differentiation.

---

## Demo 1: Complex Domain Coloring

**File**: `examples/visualization/complex_domain.rs`

### Features

1. **Function Selection**
   - Dropdown for common functions: z², z³, 1/z, eᶻ, sin(z), etc.
   - Custom function input (optional)

2. **Domain Coloring Scheme**
   - Hue = argument (angle) of f(z)
   - Brightness = |f(z)| (with logarithmic scaling)
   - Contour lines for constant |f(z)|

3. **Interactive Exploration**
   - Pan and zoom
   - Hover to see z and f(z) values
   - Click to trace function evaluation

4. **Zero/Pole Visualization**
   - Zeros: all colors meet at a point
   - Poles: colors cycle infinitely

### Domain Coloring Explanation

```
For f(z) at each pixel:
┌────────────────────────────────────────┐
│                                        │
│  Color = HSL(arg(f(z)), 1.0, L)       │
│                                        │
│  where L = brightness from |f(z)|     │
│                                        │
│  Hue wheel:                           │
│       0° = red    (positive real)      │
│      90° = yellow (positive imag)      │
│     180° = cyan   (negative real)      │
│     270° = blue   (negative imag)      │
│                                        │
└────────────────────────────────────────┘

f(z) = z²:
┌────────────────────────────────────────┐
│ ██████████████████████████████████████ │
│ ██████████▓▓▓▓▓▓▓▓▓▓▓▓████████████████ │
│ ████████▓▓░░░░░░░░░░░░▓▓██████████████ │
│ ██████▓▓░░            ░░▓▓████████████ │
│ ████▓▓░░      ●        ░░▓▓██████████ │
│ ██████▓▓░░   zero     ░░▓▓████████████ │
│ ████████▓▓░░░░░░░░░░░░▓▓██████████████ │
│ ██████████▓▓▓▓▓▓▓▓▓▓▓▓████████████████ │
│ ██████████████████████████████████████ │
└────────────────────────────────────────┘
Colors cycle TWICE around zero (degree 2)
```

### Implementation

```rust
use clifford::specialized::complex::Complex;

pub struct ComplexDomainDemo {
    function: ComplexFunction,
    view_center: Complex<f64>,
    view_scale: f64,
    show_contours: bool,
    show_grid: bool,
}

#[derive(Clone, Copy)]
enum ComplexFunction {
    Square,      // z²
    Cube,        // z³
    Reciprocal,  // 1/z
    Exp,         // eᶻ
    Sin,         // sin(z)
    Cos,         // cos(z)
    Log,         // log(z)
    Sqrt,        // √z
}

impl ComplexFunction {
    fn evaluate(&self, z: Complex<f64>) -> Complex<f64> {
        match self {
            Self::Square => z * z,
            Self::Cube => z * z * z,
            Self::Reciprocal => Complex::one() / z,
            Self::Exp => z.exp(),
            Self::Sin => z.sin(),
            Self::Cos => z.cos(),
            Self::Log => z.ln(),
            Self::Sqrt => z.sqrt(),
        }
    }
}

fn domain_color(z: Complex<f64>) -> egui::Color32 {
    let hue = z.arg() / std::f64::consts::TAU; // 0 to 1
    let mag = z.norm();
    let lightness = 1.0 - 1.0 / (1.0 + mag.ln().abs() * 0.1);
    hsl_to_rgb(hue as f32, 1.0, lightness as f32)
}
```

### Visualization Layout

```
┌─────────────────────────────────────────────────────────────┐
│  Complex Domain Coloring                               [×]  │
├─────────────────────────────────────┬───────────────────────┤
│                                     │ Function              │
│  ████████████████████████████████   │ [z²        ▼]         │
│  ████████████████████████████████   │                       │
│  ████████████████████████████████   │ Available:            │
│  ████████████████████████████████   │  z², z³, 1/z, eᶻ     │
│  ████████████████████████████████   │  sin(z), cos(z)      │
│  ████████████████████████████████   │  log(z), √z          │
│  ████████████████████████████████   │                       │
│  ████████████████████████████████   │ View                  │
│                                     │ Center: 0.0 + 0.0i    │
│  Hover: z = 1.5 + 0.5i             │ Scale: [====●====]    │
│         f(z) = 2.0 + 1.5i          │                       │
│                                     │ [✓] Show contours     │
│                                     │ [✓] Show grid         │
└─────────────────────────────────────┴───────────────────────┘
```

---

## Demo 2: Mandelbrot/Julia Explorer

**File**: `examples/visualization/complex_fractal.rs`

### Features

1. **Mandelbrot Set**
   - Classic z² + c iteration
   - Color by escape time
   - Deep zoom capability

2. **Julia Sets**
   - Select c parameter
   - See corresponding Julia set
   - Animate c along path

3. **Connection Display**
   - Show c location in Mandelbrot
   - Show corresponding Julia set side-by-side

### Implementation

```rust
pub struct FractalDemo {
    mode: FractalMode,
    julia_c: Complex<f64>,
    max_iterations: usize,
    view: ViewState,
}

enum FractalMode {
    Mandelbrot,
    Julia,
    SideBySide,
}

fn mandelbrot_escape(c: Complex<f64>, max_iter: usize) -> Option<usize> {
    let mut z = Complex::zero();
    for i in 0..max_iter {
        z = z * z + c;
        if z.norm_squared() > 4.0 {
            return Some(i);
        }
    }
    None // In the set
}

fn julia_escape(z0: Complex<f64>, c: Complex<f64>, max_iter: usize) -> Option<usize> {
    let mut z = z0;
    for i in 0..max_iter {
        z = z * z + c;
        if z.norm_squared() > 4.0 {
            return Some(i);
        }
    }
    None
}
```

---

## Demo 3: Dual Number Autodiff

**File**: `examples/visualization/dual_autodiff.rs`

### Features

1. **Function Definition**
   - Select from common functions
   - See function formula

2. **Evaluation Point**
   - Slider to move x
   - Show f(x) value

3. **Tangent Line**
   - Computed via dual numbers
   - Tangent = f(x) + f'(x)(t - x)
   - Line updates in real-time

4. **Dual Number Display**
   - Show dual computation: f(x + ε)
   - Real part = f(x)
   - Dual part = f'(x)

### The Magic of Dual Numbers

```
Dual number: a + bε  where ε² = 0

For any smooth function f:
  f(x + ε) = f(x) + f'(x)·ε

Example: f(x) = x²
  f(x + ε) = (x + ε)²
           = x² + 2xε + ε²
           = x² + 2xε + 0
           = x² + 2x·ε
             ↑     ↑
           f(x)  f'(x)

Automatic differentiation without limits!
```

### Implementation

```rust
use clifford::specialized::dual::Dual;

pub struct DualAutodiffDemo {
    function: DiffFunction,
    x: f64,
    show_tangent: bool,
    show_secant: bool,
    secant_h: f64,
}

#[derive(Clone, Copy)]
enum DiffFunction {
    Square,      // x²
    Cube,        // x³
    Sin,         // sin(x)
    Exp,         // eˣ
    Log,         // ln(x)
    Sqrt,        // √x
    Custom,      // x³ - 2x + 1
}

impl DiffFunction {
    fn evaluate_dual(&self, x: Dual<f64>) -> Dual<f64> {
        match self {
            Self::Square => x * x,
            Self::Cube => x * x * x,
            Self::Sin => x.sin(),
            Self::Exp => x.exp(),
            Self::Log => x.ln(),
            Self::Sqrt => x.sqrt(),
            Self::Custom => x * x * x - Dual::from_real(2.0) * x + Dual::from_real(1.0),
        }
    }

    fn derivative_at(&self, x: f64) -> (f64, f64) {
        let dual_x = Dual::new(x, 1.0); // x + ε
        let result = self.evaluate_dual(dual_x);
        (result.real(), result.dual()) // (f(x), f'(x))
    }
}
```

### Visualization Layout

```
┌─────────────────────────────────────────────────────────────┐
│  Dual Numbers - Automatic Differentiation              [×]  │
├─────────────────────────────────────┬───────────────────────┤
│                                     │ Function              │
│       ╭                             │ f(x) = [sin(x)   ▼]   │
│      ╱ │                            │                       │
│     ╱  │                            │ Evaluation Point      │
│ ───●───┼────────────────────        │ x = [=====●=====]     │
│   ╱│   │                            │     1.57              │
│  ╱ │   │                            │                       │
│ ╱  │   │  f(x)                      │ ─────────────────     │
│    │   │                            │ Dual Computation:     │
│    │   │                            │                       │
│    │   │ tangent                    │ f(x + ε) = sin(x + ε) │
│    x                                │          = sin(x)·1   │
│                                     │          + cos(x)·ε   │
│                                     │                       │
│                                     │ Result:               │
│                                     │   f(x)  = 1.000       │
│                                     │   f'(x) = 0.000       │
│                                     │                       │
│                                     │ [✓] Show tangent      │
│                                     │ [✓] Show secant       │
│                                     │ h = [====●====]       │
└─────────────────────────────────────┴───────────────────────┘
```

### Comparison with Numerical Differentiation

```
┌────────────────────────────────────────┐
│ Method Comparison at x = 1.5           │
│                                        │
│ Dual number (exact):                   │
│   f'(x) = 0.070737...                 │
│                                        │
│ Numerical (h = 0.01):                  │
│   (f(x+h) - f(x))/h = 0.070712...     │
│   Error: 0.000025                      │
│                                        │
│ Numerical (h = 0.001):                 │
│   (f(x+h) - f(x))/h = 0.070735...     │
│   Error: 0.000002                      │
│                                        │
│ Dual numbers: exact, no truncation!    │
└────────────────────────────────────────┘
```

---

## Demo 4: Taylor Series Visualization

**File**: `examples/visualization/dual_taylor.rs`

### Features

1. **Taylor Expansion**
   - Use iterated dual numbers for higher derivatives
   - Show Taylor polynomial approximation

2. **Order Selection**
   - Slider for Taylor order (1-6)
   - See approximation improve

3. **Error Visualization**
   - Show difference between function and Taylor polynomial
   - Highlight convergence radius

---

## Implementation Tasks

### complex_domain.rs
1. [ ] Domain coloring shader/computation
2. [ ] Function evaluation for each pixel
3. [ ] Hue from argument, brightness from magnitude
4. [ ] Pan and zoom controls
5. [ ] Function selection dropdown
6. [ ] Contour line overlay (optional)
7. [ ] Hover info display

### complex_fractal.rs
1. [ ] Mandelbrot iteration and coloring
2. [ ] Julia set iteration and coloring
3. [ ] c parameter selection
4. [ ] Side-by-side view
5. [ ] Deep zoom support
6. [ ] Iteration count slider

### dual_autodiff.rs
1. [ ] Function curve plotting
2. [ ] Tangent line computation via dual numbers
3. [ ] x position slider
4. [ ] Dual computation display
5. [ ] Secant line comparison
6. [ ] Numerical derivative comparison

### dual_taylor.rs
1. [ ] Higher-order derivatives via nested duals
2. [ ] Taylor polynomial computation
3. [ ] Order selection slider
4. [ ] Approximation vs actual plot
5. [ ] Error visualization

## Verification

```bash
cargo run --example complex_domain --release
cargo run --example complex_fractal --release
cargo run --example dual_autodiff --release
cargo run --example dual_taylor --release
```

### Checklist
- [ ] Domain coloring shows correct color patterns
- [ ] Zeros have all colors meeting, poles have infinite cycling
- [ ] Mandelbrot set renders correctly
- [ ] Julia sets change with c parameter
- [ ] Dual number derivative matches analytical
- [ ] Tangent line is tangent to curve
- [ ] Taylor approximation improves with order
