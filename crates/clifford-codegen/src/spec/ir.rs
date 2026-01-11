//! Intermediate representation for parsed algebra specifications.
//!
//! These types represent the validated, processed form of a TOML specification.
//! They contain all information needed for code generation.

use std::collections::HashMap;

/// Parsed algebra specification.
///
/// This is the main output of the parser, containing all information
/// needed to generate code for an algebra.
#[derive(Debug, Clone)]
pub struct AlgebraSpec {
    /// Algebra identifier (e.g., "euclidean3").
    pub name: String,
    /// Rust module path (e.g., "euclidean::dim3").
    pub module_path: Option<String>,
    /// Documentation string.
    pub description: Option<String>,
    /// Metric signature.
    pub signature: SignatureSpec,
    /// Custom blade name mappings (blade index -> name).
    pub blade_names: HashMap<usize, String>,
    /// Type definitions.
    pub types: Vec<TypeSpec>,
    /// Code generation options.
    pub options: GenerationOptions,
}

/// Metric signature specification.
///
/// Defines the basis vectors and their metric (how they square).
#[derive(Debug, Clone)]
pub struct SignatureSpec {
    /// All basis vectors with their properties.
    pub basis: Vec<BasisVector>,
    /// Number of positive-square basis vectors.
    pub p: usize,
    /// Number of negative-square basis vectors.
    pub q: usize,
    /// Number of zero-square (degenerate) basis vectors.
    pub r: usize,
}

impl SignatureSpec {
    /// Total dimension (number of basis vectors).
    #[inline]
    pub fn dim(&self) -> usize {
        self.p + self.q + self.r
    }

    /// Total number of blades (2^dim).
    #[inline]
    pub fn num_blades(&self) -> usize {
        1 << self.dim()
    }
}

/// A single basis vector in the signature.
#[derive(Debug, Clone)]
pub struct BasisVector {
    /// Name of the basis vector (e.g., "e1", "x").
    pub name: String,
    /// Index in the algebra (0-based).
    pub index: usize,
    /// Metric value: +1, -1, or 0.
    pub metric: i8,
}

/// A type definition in the specification.
#[derive(Debug, Clone)]
pub struct TypeSpec {
    /// Type name (e.g., "Rotor", "Vector").
    pub name: String,
    /// Grades contained in this type.
    pub grades: Vec<usize>,
    /// Documentation string.
    pub description: Option<String>,
    /// Fields in order (matching blade layout).
    pub fields: Vec<FieldSpec>,
    /// If this type aliases another (same storage layout).
    pub alias_of: Option<String>,
    /// Constraints (unit, nonzero, etc.).
    pub constraints: Vec<ConstraintSpec>,
}

/// A field in a type.
#[derive(Debug, Clone)]
pub struct FieldSpec {
    /// Field name (e.g., "x", "xy").
    pub name: String,
    /// Blade index this field holds.
    pub blade_index: usize,
    /// Grade of the blade.
    pub grade: usize,
}

/// A constraint on a type (unit, nonzero, etc.).
#[derive(Debug, Clone)]
pub struct ConstraintSpec {
    /// Kind of constraint.
    pub kind: ConstraintKind,
    /// Generated wrapper type name (e.g., "UnitRotor").
    pub wrapper_name: String,
    /// Norm type for unit constraints.
    pub norm_type: Option<NormType>,
    /// Condition expression for custom constraints.
    pub condition: Option<String>,
    /// Safe constructors for the constrained type.
    pub constructors: Vec<ConstructorSpec>,
    /// Operations that preserve this constraint.
    pub preserving_ops: Vec<String>,
}

/// Kind of constraint.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum ConstraintKind {
    /// ‖x‖ = 1
    Unit,
    /// ‖x‖ ≠ 0
    NonZero,
    /// Context-dependent canonical form
    Normalized,
    /// x · x = 0 (null vector in CGA)
    Null,
    /// In ideal subspace
    Ideal,
}

impl ConstraintKind {
    /// Returns the prefix for the wrapper type name.
    pub fn wrapper_prefix(&self) -> &'static str {
        match self {
            ConstraintKind::Unit => "Unit",
            ConstraintKind::NonZero => "NonZero",
            ConstraintKind::Normalized => "Normalized",
            ConstraintKind::Null => "Null",
            ConstraintKind::Ideal => "Ideal",
        }
    }
}

/// Type of norm used for unit constraints.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum NormType {
    /// s² + x² + y² + z² + ...
    Euclidean,
    /// PGA motor norm (special handling for degenerate components).
    Motor,
    /// CGA-specific norm.
    Cga,
    /// Custom norm (specified by condition).
    Custom,
}

/// A constructor for a constrained type.
#[derive(Debug, Clone)]
pub struct ConstructorSpec {
    /// Constructor name (e.g., "identity", "from_angle").
    pub name: String,
    /// Parameters.
    pub params: Vec<ParamSpec>,
}

/// A parameter in a constructor.
#[derive(Debug, Clone)]
pub struct ParamSpec {
    /// Parameter name.
    pub name: String,
    /// Type annotation (e.g., "T", "Vector<T>").
    pub ty: String,
}

/// Code generation options.
#[derive(Debug, Clone, Default)]
pub struct GenerationOptions {
    /// Generate serde implementations.
    pub generate_serde: bool,
    /// Generate Arbitrary implementations for proptest.
    pub generate_arbitrary: bool,
    /// Generate test modules.
    pub generate_tests: bool,
}
