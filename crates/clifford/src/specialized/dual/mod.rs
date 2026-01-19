//! Dual Numbers - Cl(0,0,1)
//!
//! The dual numbers are a 2-dimensional algebra over the reals
//! with basis `{1, ε}` where `ε² = 0` (nilpotent).
//!
//! # Structure
//!
//! A dual number `z = a + bε` consists of:
//! - `a`: real part (grade 0)
//! - `b`: dual part (grade 1)
//!
//! # Multiplication
//!
//! The multiplication rule mirrors the product rule for derivatives:
//! ```text
//! (a + bε)(c + dε) = ac + (ad + bc)ε
//! ```
//!
//! Note that `ε² = 0`, so the `bdε²` term vanishes.
//!
//! # Automatic Differentiation
//!
//! Dual numbers are the foundation of **forward-mode automatic differentiation**.
//! When you evaluate `f(a + ε)`, the result is `f(a) + f'(a)·ε`:
//!
//! ```text
//! f(x + ε) = f(x) + f'(x)·ε
//! ```
//!
//! For example, to compute the derivative of `f(x) = x²` at `x = 3`:
//! ```text
//! f(3 + ε) = (3 + ε)² = 9 + 6ε
//! ```
//! So `f(3) = 9` and `f'(3) = 6`.
//!
//! # Norm
//!
//! The dual numbers have a **degenerate norm**:
//! ```text
//! z * involute(z) = (a + bε)(a - bε) = a² - b²ε² = a²
//! ```
//!
//! The norm only depends on the real part! All pure dual numbers `bε`
//! have norm zero, making them zero divisors.
//!
//! # Zero Divisors
//!
//! Every pure dual number is a zero divisor:
//! ```text
//! (bε)(cε) = bc·ε² = 0
//! ```
//!
//! # References
//!
//! - [Dual number](https://en.wikipedia.org/wiki/Dual_number)
//! - [Automatic differentiation](https://en.wikipedia.org/wiki/Automatic_differentiation)

mod extensions;
mod generated;

pub use generated::types::*;
