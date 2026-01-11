//! Symbolic constraint verification using Symbolica.
//!
//! This module provides tools for symbolically verifying that product operations
//! preserve algebraic constraints. It uses the Symbolica computer algebra system
//! to prove properties at codegen time.
//!
//! # Overview
//!
//! The module provides:
//!
//! - [`ConstraintParser`]: Parse constraint expressions from strings
//! - [`SymbolicProduct`]: Compute symbolic product outputs
//! - [`ConstraintVerifier`]: Verify that outputs satisfy constraints
//!
//! # Example
//!
//! ```ignore
//! use clifford_codegen::symbolic::{ConstraintParser, ConstraintVerifier};
//!
//! // Parse a unit norm constraint
//! let parser = ConstraintParser::new();
//! let constraint = parser.parse("s * s + xy * xy = 1", &["s", "xy"])?;
//!
//! // Verify that Rotor * Rotor satisfies the constraint
//! let verifier = ConstraintVerifier::new();
//! let result = verifier.verify(&input_constraints, &constraint, &output_fields);
//! ```

mod parser;
mod product;
mod verify;

pub use parser::{ConstraintExpr, ConstraintParser, ParseError};
pub use product::{ProductKind, SymbolicField, SymbolicProduct};
pub use verify::{ConstraintVerifier, VerificationResult};
