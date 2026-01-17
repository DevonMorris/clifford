//! Hyperbolic geometric algebras.
//!
//! This module provides types for algebras related to hyperbolic structures.
//!
//! # Available Algebras
//!
//! - **Root module**: Hyperbolic numbers (split-complex) Cl(1,0,0)
//! - **[`dim2`]**: 2D hyperbolic plane (Lobachevsky plane) Cl(2,1,0)
//!
//! # Hyperbolic Numbers - Cl(1,0,0)
//!
//! The hyperbolic numbers are a 2-dimensional algebra over the reals
//! with basis `{1, j}` where `j² = +1` (unlike complex numbers where `i² = -1`).
//!
//! A hyperbolic number `z = a + bj` consists of:
//! - `a`: real part (grade 0)
//! - `b`: hyperbolic part (grade 1)
//!
//! ## Norm and Causal Character
//!
//! The hyperbolic numbers use **grade involution** for their norm:
//! ```text
//! involute(a + bj) = a - bj    (grade involution negates odd grades)
//! z * involute(z) = (a + bj)(a - bj) = a² - b²
//! ```
//!
//! This gives an **indefinite** norm that can be positive, negative, or zero:
//! - **Timelike** (`a² > b²`): norm_squared > 0
//! - **Spacelike** (`a² < b²`): norm_squared < 0
//! - **Lightlike** (`a² = b²`): norm_squared = 0 (zero divisors!)
//!
//! ## Zero Divisors
//!
//! Unlike the complex numbers, hyperbolic numbers have **zero divisors**:
//! ```text
//! (a + aj)(a - aj) = a² - a² = 0
//! ```
//!
//! The elements `1 ± j` are idempotent (square to themselves) and
//! lightlike elements form the "light cone" in the hyperbolic plane.
//!
//! # References
//!
//! - [Split-complex number](https://en.wikipedia.org/wiki/Split-complex_number)
//! - [Hyperbolic number](https://en.wikipedia.org/wiki/Hyperbolic_number)
//! - [Hyperbolic geometry](https://en.wikipedia.org/wiki/Hyperbolic_geometry)

mod generated;

pub mod dim2;

pub use generated::types::*;
