//! 2D Hyperbolic Projective Geometry - Cl(2,1,0)
//!
//! The hyperbolic plane (Lobachevsky plane) is a non-Euclidean geometry
//! with constant negative curvature. This algebra provides types for
//! working with hyperbolic geometry using the hyperboloid model.
//!
//! # Structure
//!
//! | Grade | Blades | Squares | Description |
//! |-------|--------|---------|-------------|
//! | 0 | 1 | +1 | Scalar |
//! | 1 | e₁, e₂, e₃ | +1, +1, -1 | Points (homogeneous coordinates) |
//! | 2 | e₁₂, e₁₃, e₂₃ | -1, +1, +1 | Lines (geodesics) |
//! | 3 | e₁₂₃ | +1 | Pseudoscalar |
//!
//! # Indefinite Metric
//!
//! Unlike Euclidean geometry, the hyperbolic plane has an **indefinite metric**.
//! The norm of elements can be:
//!
//! - **Timelike** (norm² < 0): Ordinary points in the hyperbolic plane
//! - **Null/Lightlike** (norm² = 0): Ideal points on the "absolute"
//! - **Spacelike** (norm² > 0): Ultra-ideal points beyond infinity
//!
//! # Point Classification
//!
//! Points in homogeneous coordinates (x, y, t) are classified by x² + y² - t²:
//!
//! | Condition | Type | Interpretation |
//! |-----------|------|----------------|
//! | x² + y² - t² < 0 | Timelike | Ordinary point in hyperbolic plane |
//! | x² + y² - t² = 0 | Null | Ideal point (on the absolute/boundary) |
//! | x² + y² - t² > 0 | Spacelike | Ultra-ideal point (outside the plane) |
//!
//! # Comparison with Other Geometries
//!
//! | Property | Euclidean | Elliptic | Hyperbolic |
//! |----------|-----------|----------|------------|
//! | Curvature | 0 | Positive | Negative |
//! | Parallel lines | 1 | 0 | Infinitely many |
//! | Triangle angles | = 180° | > 180° | < 180° |
//! | Signature | Cl(3,0,0) | Cl(3,0,0) | Cl(2,1,0) |
//!
//! # Applications
//!
//! - **Non-Euclidean geometry**: Study of negatively curved spaces
//! - **Special relativity**: Rapidity addition, Lorentz boosts
//! - **Hyperbolic tessellations**: Escher-like tilings, circle packings
//! - **Complex analysis**: Poincaré disk model, modular forms
//! - **Machine learning**: Hyperbolic embeddings for hierarchical data

mod generated;

pub use generated::types::*;
