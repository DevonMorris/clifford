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
    /// Product specifications.
    pub products: ProductsSpec,
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
    /// Constraint expression (e.g., "s * s + xy * xy = 1" for unit rotors).
    pub constraint: Option<String>,
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

/// Product specifications for all product types.
#[derive(Debug, Clone, Default)]
pub struct ProductsSpec {
    /// Geometric product entries.
    pub geometric: Vec<ProductEntry>,
    /// Outer (wedge) product entries.
    pub outer: Vec<ProductEntry>,
    /// Left contraction entries.
    pub left_contraction: Vec<ProductEntry>,
    /// Right contraction entries.
    pub right_contraction: Vec<ProductEntry>,
    /// Regressive product entries.
    pub regressive: Vec<ProductEntry>,
    /// Scalar product entries.
    pub scalar: Vec<ProductEntry>,
}

/// A single product entry specifying lhs × rhs → output.
#[derive(Debug, Clone)]
pub struct ProductEntry {
    /// Left-hand side type name (e.g., "Vector", "UnitRotor").
    pub lhs: String,
    /// Right-hand side type name.
    pub rhs: String,
    /// Output type name.
    pub output: String,
    /// Whether the output is a constrained type (uses new_unchecked).
    pub output_constrained: bool,
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
