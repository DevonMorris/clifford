//! Code generation for geometric algebra types.
//!
//! This module provides code generators that produce Rust code from
//! algebra specifications. The generated code includes struct definitions,
//! accessors, constructors, and basic operations.
//!
//! # Example
//!
//! ```
//! use clifford_codegen::algebra::Algebra;
//! use clifford_codegen::codegen::TypeGenerator;
//! use clifford_codegen::spec::parse_spec;
//!
//! let spec = parse_spec(r#"
//! [algebra]
//! name = "test"
//!
//! [signature]
//! positive = ["e1", "e2"]
//!
//! [types.Vector]
//! grades = [1]
//! fields = ["x", "y"]
//! "#).unwrap();
//!
//! let algebra = Algebra::euclidean(2);
//! let generator = TypeGenerator::new(&spec, &algebra);
//!
//! let tokens = generator.generate_types_file();
//! let code = tokens.to_string();
//!
//! assert!(code.contains("pub struct Vector"));
//! ```

mod constraints;
mod format;
mod products;
mod traits;
mod types;

pub use constraints::ConstraintGenerator;
pub use format::format_tokens;
pub use products::ProductGenerator;
pub use traits::TraitsGenerator;
pub use types::TypeGenerator;
