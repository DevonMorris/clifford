//! Conformal Geometric Algebras.
//!
//! This module provides specialized types for Conformal Geometric Algebras (CGA),
//! which embed Euclidean space into a higher-dimensional space with one additional
//! positive and one additional negative basis vector.
//!
//! # Available Algebras
//!
//! - **[`dim2`]**: 2D Conformal GA Cl(3,1,0) - embeds 2D Euclidean space
//! - **[`dim3`]**: 3D Conformal GA Cl(4,1,0) - embeds 3D Euclidean space
//!
//! # Conformal Model
//!
//! The conformal model represents geometric objects elegantly:
//!
//! | Grade | 2D CGA | 3D CGA |
//! |-------|--------|--------|
//! | 1 | Round points | Round points |
//! | 2 | Point pairs, flat points | Dipoles, flat points |
//! | 3 | Circles, lines | Circles, lines |
//! | 4 | Pseudoscalar | Spheres, planes |
//! | 5 | - | Pseudoscalar |
//!
//! # Key Properties
//!
//! - Points are represented as **null vectors** on the null cone
//! - The inner product of two points gives -(squared distance)/2
//! - Transformations (rotation, translation, dilation) are versors
//! - Intersections are computed with the wedge product
//!
//! # Reference
//!
//! See <https://conformalgeometricalgebra.org/wiki/index.php?title=Main_Page>
//! for detailed documentation on CGA operations and conventions.

pub mod dim2;
pub mod dim3;
