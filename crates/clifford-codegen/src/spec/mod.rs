//! Algebra specification parsing and representation.
//!
//! This module provides types and functions for parsing TOML specifications
//! that describe geometric algebras and their types.
//!
//! # Specification Format
//!
//! Specifications are written in TOML with the following structure:
//!
//! ```toml
//! [algebra]
//! name = "euclidean3"                    # Required: identifier
//! module_path = "euclidean::dim3"        # Optional: Rust module path
//! description = "3D Euclidean GA"        # Optional: documentation
//!
//! [signature]
//! positive = ["e1", "e2", "e3"]          # Basis vectors squaring to +1
//! negative = []                           # Basis vectors squaring to -1
//! zero = []                               # Basis vectors squaring to 0
//!
//! [blades]
//! e12 = "xy"                              # Custom blade names
//! e13 = "xz"
//! e23 = "yz"
//!
//! [types.Vector]
//! grades = [1]                            # Which grades this type contains
//! fields = ["x", "y", "z"]               # Field names
//!
//! [types.Rotor]
//! grades = [0, 2]
//! fields = ["s", "xy", "xz", "yz"]
//! ```
//!
//! # Example
//!
//! ```
//! use clifford_codegen::spec::parse_spec;
//!
//! let spec = parse_spec(r#"
//! [algebra]
//! name = "euclidean2"
//!
//! [signature]
//! positive = ["e1", "e2"]
//!
//! [types.Vector]
//! grades = [1]
//! fields = ["x", "y"]
//! "#).unwrap();
//!
//! assert_eq!(spec.name, "euclidean2");
//! assert_eq!(spec.signature.dim(), 2);
//! ```

mod bundled;
mod error;
mod ir;
mod parser;
mod raw;

pub use bundled::{EUCLIDEAN2, EUCLIDEAN3};
pub use error::ParseError;
pub use ir::{
    AlgebraSpec, BasisVector, FieldSpec, InvolutionKind, NormSpec, ProductEntry, ProductsSpec,
    SignatureSpec, TypeSpec,
};
pub use parser::parse_spec;
