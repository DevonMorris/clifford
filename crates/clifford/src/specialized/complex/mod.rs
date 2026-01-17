//! Complex Numbers - Cl(0,1,0)
//!
//! The complex numbers are a 2-dimensional algebra over the reals
//! with basis `{1, i}` where `i² = -1`.
//!
//! # Structure
//!
//! A complex number `z = a + bi` consists of:
//! - `a`: real part (grade 0)
//! - `b`: imaginary part (grade 1)
//!
//! # Norm
//!
//! The complex numbers use **Clifford conjugate** for their norm:
//! ```text
//! conjugate(a + bi) = a - bi    (Clifford conjugate negates odd grades)
//! z * conjugate(z) = (a + bi)(a - bi) = a² + b²
//! ```
//!
//! This gives a **positive-definite** norm (always non-negative):
//! - `|z|² = a² + b² >= 0`
//! - `|z|² = 0` if and only if `z = 0`
//!
//! # Division Algebra
//!
//! Unlike hyperbolic numbers, complex numbers have **no zero divisors**.
//! Every non-zero complex number has a multiplicative inverse:
//! ```text
//! z⁻¹ = conjugate(z) / |z|²
//! ```
//!
//! # Comparison with Hyperbolic Numbers
//!
//! | Property | Complex Cl(0,1,0) | Hyperbolic Cl(1,0,0) |
//! |----------|-------------------|----------------------|
//! | Unit squares to | `i² = -1` | `j² = +1` |
//! | Norm | `a² + b²` | `a² - b²` |
//! | Zero divisors | None | Yes |
//!
//! # References
//!
//! - [Complex number](https://en.wikipedia.org/wiki/Complex_number)
//! - [Clifford algebra](https://en.wikipedia.org/wiki/Clifford_algebra)

mod generated;

pub use generated::types::*;
