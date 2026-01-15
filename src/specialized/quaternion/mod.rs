//! Quaternions - Cl(0,2,0)
//!
//! Quaternions are a 4-dimensional algebra over the reals commonly used
//! for representing 3D rotations without gimbal lock.
//!
//! # Structure
//!
//! A quaternion `q = w + xi + yj + zk` consists of:
//! - `w`: scalar/real part (grade 0)
//! - `x, y`: imaginary components (grade 1, maps to basis vectors e1, e2)
//! - `z`: bivector component (grade 2, maps to e12 = e1·e2)
//!
//! # Basis Mapping
//!
//! The quaternion basis maps to Cl(0,2,0) as:
//! - `1 → 1` (scalar)
//! - `i → e1` (basis vector)
//! - `j → e2` (basis vector)
//! - `k → e12 = e1·e2` (bivector)
//!
//! With `e1² = e2² = -1`, this gives:
//! - `i² = j² = k² = -1`
//! - `ij = k`, `jk = i`, `ki = j`
//! - `ji = -k`, `kj = -i`, `ik = -j`
//!
//! # Norm
//!
//! Quaternions use **Clifford conjugate** for their norm:
//! ```text
//! conjugate(w + xi + yj + zk) = w - xi - yj - zk
//! q * conjugate(q) = w² + x² + y² + z²
//! ```
//!
//! This gives a **positive-definite** norm.
//!
//! # Division Algebra
//!
//! Like complex numbers, quaternions form a **division algebra** with no
//! zero divisors. Every non-zero quaternion has a multiplicative inverse:
//! ```text
//! q⁻¹ = conjugate(q) / |q|²
//! ```
//!
//! # References
//!
//! - [Quaternion](https://en.wikipedia.org/wiki/Quaternion)
//! - [Clifford algebra](https://en.wikipedia.org/wiki/Clifford_algebra)

mod generated;

pub use generated::types::*;
