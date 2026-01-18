//! Raw TOML structures for deserialization.
//!
//! These types directly match the TOML file format and are used by serde
//! to deserialize specification files. After deserialization, they are
//! validated and converted to the IR types in [`super::ir`].

use std::collections::HashMap;

use serde::Deserialize;

/// Default value for the `complete` flag (true).
///
/// Algebras are expected to be complete by default (all products have output types).
/// Set `complete = false` explicitly for intentionally incomplete algebras.
fn default_complete() -> bool {
    true
}

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
    /// Whether to enforce algebra completeness.
    ///
    /// When `true` (default), codegen will error if any product between defined types
    /// produces grades that don't match a defined output type.
    /// Set to `false` to allow partial algebras.
    #[serde(default = "default_complete")]
    pub complete: bool,
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
    /// Explicit field-to-blade mappings (required for all types).
    ///
    /// Each entry maps a field name to a specific blade. Blade names use the
    /// format "e123" where digits are 1-based basis vector indices.
    ///
    /// Non-canonical blade orderings (e.g., "e20" instead of "e02") are supported
    /// and will automatically apply the correct sign correction.
    ///
    /// Example:
    /// ```toml
    /// field_map = [
    ///   { name = "nx", blade = "e20" },  # e20 = -e02
    ///   { name = "ny", blade = "e01" },
    ///   { name = "d", blade = "e12" }
    /// ]
    /// ```
    #[serde(default)]
    pub field_map: Vec<RawFieldMapping>,
    /// Type this aliases.
    pub alias_of: Option<String>,
}

/// A single field-to-blade mapping in the TOML specification.
#[derive(Debug, Deserialize)]
pub struct RawFieldMapping {
    /// Field name (e.g., "x", "nx", "dist").
    pub name: String,
    /// Blade name (e.g., "e1", "e20", "e123").
    ///
    /// Non-canonical orderings like "e20" (instead of "e02") are supported
    /// and will compute the appropriate sign correction.
    pub blade: String,
}
