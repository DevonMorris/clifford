//! Conformal Geometric Algebras.
//!
//! This module provides specialized types for Conformal Geometric Algebras (CGA),
//! which embed Euclidean space into a higher-dimensional space with one additional
//! positive and one additional negative basis vector.
//!
//! # Available Algebras
//!
//! - **[`dim3`]**: 3D Conformal GA Cl(4,1,0) - embeds 3D Euclidean space
//!
//! # Conformal Model
//!
//! The conformal model represents geometric objects elegantly:
//!
//! | Grade | Geometric Object |
//! |-------|-----------------|
//! | 1 | Round points (null vectors are actual points) |
//! | 2 | Dipoles (point pairs), flat points |
//! | 3 | Circles, lines (circles through infinity) |
//! | 4 | Spheres, planes (spheres through infinity) |
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

pub mod dim3;
