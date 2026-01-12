//! Error types for specification parsing.
//!
//! These errors provide detailed information about what went wrong
//! during parsing, making it easy to fix specification files.

use thiserror::Error;

/// Errors that can occur when parsing an algebra specification.
#[derive(Debug, Error)]
pub enum ParseError {
    /// TOML syntax error.
    #[error("TOML parse error: {0}")]
    Toml(#[from] toml::de::Error),

    /// Algebra dimension exceeds the maximum supported (6).
    #[error("dimension {0} exceeds maximum of 6")]
    DimensionTooLarge(usize),

    /// Signature must have at least one basis vector.
    #[error("signature must have at least one basis vector")]
    EmptySignature,

    /// Duplicate basis vector name in signature.
    #[error("duplicate basis vector name: '{0}'")]
    DuplicateBasisName(String),

    /// Type references a grade that exceeds the algebra dimension.
    #[error("type '{type_name}' has invalid grade {grade} (max: {max})")]
    InvalidGrade {
        /// The type name.
        type_name: String,
        /// The invalid grade.
        grade: usize,
        /// Maximum valid grade.
        max: usize,
    },

    /// Type has wrong number of fields for its grades.
    #[error("type '{type_name}' has {got} fields but needs {expected}")]
    FieldCountMismatch {
        /// The type name.
        type_name: String,
        /// Expected number of fields.
        expected: usize,
        /// Actual number of fields.
        got: usize,
    },

    /// Duplicate field name in a type.
    #[error("type '{type_name}' has duplicate field name: '{field}'")]
    DuplicateFieldName {
        /// The type name.
        type_name: String,
        /// The duplicate field name.
        field: String,
    },

    /// Duplicate type name in specification.
    #[error("duplicate type name: '{0}'")]
    DuplicateTypeName(String),

    /// Invalid blade name in blades section.
    #[error("invalid blade name: '{name}' (expected format: e1, e12, e123, etc.)")]
    InvalidBladeName {
        /// The invalid blade name.
        name: String,
    },

    /// Blade index exceeds algebra dimension.
    #[error("blade '{name}' references basis index {index} but algebra only has {dim} dimensions")]
    BladeIndexOutOfBounds {
        /// The blade name.
        name: String,
        /// The invalid index.
        index: usize,
        /// Algebra dimension.
        dim: usize,
    },

    /// Unknown type reference (e.g., in alias_of).
    #[error("unknown type reference: '{0}'")]
    UnknownType(String),

    /// Type alias references itself.
    #[error("type '{type_name}' cannot alias itself")]
    SelfAlias {
        /// The type name.
        type_name: String,
    },

    /// Type alias forms a cycle.
    #[error("type alias cycle detected involving '{type_name}'")]
    AliasCycle {
        /// The type name.
        type_name: String,
    },
}
