//! Code generator for optimized geometric algebra types.
//!
//! This crate provides tools to generate Rust code for geometric algebra
//! operations. Given an algebra specification (signature + type definitions),
//! it generates optimized, strongly-typed code for all products and operations.
//!
//! # Architecture
//!
//! - [`algebra`]: Core blade algebra computations (signs, products, grades)
//! - [`discovery`]: Automatic entity discovery from geometric constraints
//! - [`spec`]: TOML specification parsing and validation
//! - [`codegen`]: Code generation from specifications
//! - [`symbolic`]: Symbolic constraint verification (requires `symbolic` feature)
//!
//! # Example
//!
//! ```
//! use clifford_codegen::algebra::{Algebra, Blade, ProductTable};
//!
//! // Create a 3D Euclidean algebra
//! let algebra = Algebra::euclidean(3);
//!
//! // Build product table
//! let table = ProductTable::new(&algebra);
//!
//! // Compute e1 * e2
//! let e1 = Blade::basis(0);
//! let e2 = Blade::basis(1);
//! let (sign, result) = table.geometric(e1.index(), e2.index());
//!
//! assert_eq!(sign, 1);
//! assert_eq!(result, 0b11); // e12
//! ```

pub mod algebra;
pub mod codegen;
pub mod discovery;
pub mod spec;
pub mod symbolic;
