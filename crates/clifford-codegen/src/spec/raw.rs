//! Raw TOML structures for deserialization.
//!
//! These types directly match the TOML file format and are used by serde
//! to deserialize specification files. After deserialization, they are
//! validated and converted to the IR types in [`super::ir`].

use std::collections::HashMap;

use serde::Deserialize;

/// Raw algebra specification (matches TOML structure).
#[derive(Debug, Deserialize)]
pub struct RawAlgebraSpec {
    /// Algebra metadata.
    pub algebra: RawAlgebraInfo,
    /// Metric signature.
    pub signature: RawSignature,
    /// Custom blade names.
    #[serde(default)]
    pub blades: HashMap<String, String>,
    /// Type definitions.
    #[serde(default)]
    pub types: HashMap<String, RawTypeSpec>,
    /// Product definitions.
    #[serde(default)]
    pub products: RawProducts,
    /// Generation options.
    #[serde(default)]
    pub options: RawOptions,
}

/// Raw algebra info section.
#[derive(Debug, Deserialize)]
pub struct RawAlgebraInfo {
    /// Algebra identifier.
    pub name: String,
    /// Rust module path.
    pub module_path: Option<String>,
    /// Documentation.
    pub description: Option<String>,
}

/// Raw signature section.
#[derive(Debug, Deserialize, Default)]
pub struct RawSignature {
    /// Positive-square basis vectors.
    #[serde(default)]
    pub positive: Vec<String>,
    /// Negative-square basis vectors.
    #[serde(default)]
    pub negative: Vec<String>,
    /// Zero-square (degenerate) basis vectors.
    #[serde(default)]
    pub zero: Vec<String>,
}

/// Raw type specification.
#[derive(Debug, Deserialize)]
pub struct RawTypeSpec {
    /// Grades contained.
    pub grades: Vec<usize>,
    /// Documentation.
    pub description: Option<String>,
    /// Custom field names.
    #[serde(default)]
    pub fields: Vec<String>,
    /// Type this aliases.
    pub alias_of: Option<String>,
    /// Constraints.
    #[serde(default)]
    pub constraints: HashMap<String, RawConstraint>,
}

/// Raw constraint specification.
#[derive(Debug, Deserialize, Default)]
pub struct RawConstraint {
    /// Condition expression.
    pub condition: Option<String>,
}

/// Raw products section.
///
/// Each product type maps "Lhs_Rhs" -> "Output".
/// For example: `Vector_Vector = "Rotor"`.
#[derive(Debug, Deserialize, Default)]
pub struct RawProducts {
    /// Geometric product entries.
    #[serde(default)]
    pub geometric: HashMap<String, String>,
    /// Outer (wedge) product entries.
    #[serde(default)]
    pub outer: HashMap<String, String>,
    /// Left contraction entries.
    #[serde(default)]
    pub left_contraction: HashMap<String, String>,
    /// Right contraction entries.
    #[serde(default)]
    pub right_contraction: HashMap<String, String>,
    /// Regressive product entries.
    #[serde(default)]
    pub regressive: HashMap<String, String>,
    /// Scalar product entries.
    #[serde(default)]
    pub scalar: HashMap<String, String>,
}

/// Raw generation options.
#[derive(Debug, Deserialize)]
pub struct RawOptions {
    /// Generate serde implementations.
    #[serde(default)]
    pub generate_serde: bool,
    /// Generate Arbitrary implementations.
    #[serde(default = "default_true")]
    pub generate_arbitrary: bool,
    /// Generate test modules.
    #[serde(default = "default_true")]
    pub generate_tests: bool,
}

impl Default for RawOptions {
    fn default() -> Self {
        Self {
            generate_serde: false,
            generate_arbitrary: true,
            generate_tests: true,
        }
    }
}

/// Default value helper for serde.
fn default_true() -> bool {
    true
}
