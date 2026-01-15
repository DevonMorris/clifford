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
    /// Norm configuration.
    #[serde(default)]
    pub norm: RawNormSpec,
    /// Custom blade names.
    #[serde(default)]
    pub blades: HashMap<String, String>,
    /// Type definitions.
    #[serde(default)]
    pub types: HashMap<String, RawTypeSpec>,
}

/// Raw norm configuration section.
#[derive(Debug, Deserialize, Default)]
pub struct RawNormSpec {
    /// Which involution produces the primary norm.
    ///
    /// Options: "reverse" (default), "grade_involution", "clifford_conjugate"
    pub primary_involution: Option<String>,
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
}
