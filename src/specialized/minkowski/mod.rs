//! Minkowski spacetime algebras.
//!
//! This module provides specialized types for Minkowski (spacetime) geometric algebras,
//! which have mixed signature (both positive and negative squaring basis vectors).
//!
//! # Available Algebras
//!
//! - **[`dim2`]**: Minkowski plane Cl(1,1,0) - 2D spacetime
//!
//! # Indefinite Metric
//!
//! Minkowski algebras have an **indefinite metric**, meaning:
//! - Vectors can be **spacelike** (v² > 0), **timelike** (v² < 0), or **null** (v² = 0)
//! - The norm can be positive, negative, or zero for non-zero elements
//! - Null vectors form the "light cone" separating spacelike and timelike regions

pub mod dim2;
