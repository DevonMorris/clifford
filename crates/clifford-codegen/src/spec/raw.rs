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
    /// Product definitions (ignored - products are auto-inferred from types).
    #[serde(default)]
    #[allow(dead_code)]
    products: RawProducts,
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
    /// Geometric product constraint (e.g., "2*s*e0123 - 2*e12*e03 + ... = 0").
    pub geometric_constraint: Option<String>,
    /// Antiproduct constraint (often same as geometric_constraint).
    pub antiproduct_constraint: Option<String>,
    /// Field to solve for when enforcing the geometric constraint.
    /// If constraints are dependent (same expression), only this field is used.
    pub geometric_solve_for: Option<String>,
    /// Field to solve for when enforcing the antiproduct constraint.
    /// Only required when antiproduct_constraint differs from geometric_constraint.
    pub antiproduct_solve_for: Option<String>,
    /// User-defined additional constraints.
    #[serde(default)]
    pub constraints: Vec<RawUserConstraint>,
    /// Whether this type is a versor (for sandwich products).
    ///
    /// Versors are elements that can transform other elements via the
    /// sandwich product: `X' = V * X * rev(V)`.
    #[serde(default)]
    pub versor: bool,
    /// Sandwich product configuration.
    ///
    /// Specifies which types this versor can transform via sandwich product.
    #[serde(default)]
    pub sandwich: Option<RawSandwichConfig>,
}

/// Raw sandwich product configuration.
#[derive(Debug, Deserialize, Clone, Default)]
pub struct RawSandwichConfig {
    /// Target types for sandwich product.
    ///
    /// Can be:
    /// - A list of type names: `["Vector", "Bivector"]`
    /// - Empty to auto-detect targets based on grade compatibility
    #[serde(default)]
    pub targets: Vec<String>,
}

/// Raw user-defined constraint.
#[derive(Debug, Deserialize, Clone)]
pub struct RawUserConstraint {
    /// Constraint name (e.g., "unit", "normalized").
    pub name: String,
    /// Documentation.
    pub description: Option<String>,
    /// Constraint expression (e.g., "s*s + e12*e12 + e13*e13 + e23*e23 = 1").
    pub expression: String,
    /// Field to solve for. If present, this constraint reduces the parameter count.
    pub solve_for: Option<String>,
    /// Sign convention for ± ambiguity: "positive" (default) or "negative".
    #[serde(default = "default_positive")]
    pub sign: String,
    /// Enforcement method to generate (e.g., "normalize").
    pub enforce: Option<String>,
}

/// Default sign convention.
fn default_positive() -> String {
    "positive".to_string()
}

/// Raw products section.
///
/// **DEPRECATED**: Products are now auto-inferred from types.
/// This struct exists only for backward compatibility with old TOML files.
/// Any products specified in TOML are ignored.
#[derive(Debug, Deserialize, Default)]
#[allow(dead_code)]
pub struct RawProducts {
    /// Geometric product entries (ignored).
    #[serde(default)]
    geometric: HashMap<String, String>,
    /// Outer (wedge) product entries (ignored).
    #[serde(default)]
    outer: HashMap<String, String>,
    /// Left contraction entries (ignored).
    #[serde(default)]
    left_contraction: HashMap<String, String>,
    /// Right contraction entries (ignored).
    #[serde(default)]
    right_contraction: HashMap<String, String>,
    /// Regressive product entries (ignored).
    #[serde(default)]
    regressive: HashMap<String, String>,
    /// Scalar product entries (ignored).
    #[serde(default)]
    scalar: HashMap<String, String>,
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
