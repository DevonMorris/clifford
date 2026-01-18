//! Intermediate representation for parsed algebra specifications.
//!
//! These types represent the validated, processed form of a TOML specification.
//! They contain all information needed for code generation.

use std::collections::HashMap;

/// Involution kind for norm computation.
///
/// Specifies which involution the algebra uses for its canonical norm.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum InvolutionKind {
    /// Reverse: (-1)^(k(k-1)/2) for grade k.
    ///
    /// Default for most algebras (Euclidean, PGA, Minkowski).
    /// Used for versor operations regardless of norm involution.
    #[default]
    Reverse,
    /// Grade involution: (-1)^k for grade k.
    ///
    /// Used for split-complex/hyperbolic numbers.
    GradeInvolution,
    /// Clifford conjugate: composition of reverse and grade involution.
    ///
    /// Sign = (-1)^(k(k+1)/2) for grade k.
    CliffordConjugate,
}

/// Norm configuration for an algebra.
///
/// Specifies which involution produces the "primary" norm for the algebra.
/// The codegen generates `Involute` based on this setting.
#[derive(Debug, Clone, Default)]
pub struct NormSpec {
    /// Which involution to use for primary norm computation.
    ///
    /// This determines what `involute()` returns for types in this algebra.
    pub primary_involution: InvolutionKind,
}

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
    /// Norm configuration.
    pub norm: NormSpec,
    /// Custom blade name mappings (blade index -> name).
    pub blade_names: HashMap<usize, String>,
    /// Type definitions.
    pub types: Vec<TypeSpec>,
    /// Product specifications.
    pub products: ProductsSpec,
    /// Whether completeness checking is enabled.
    ///
    /// When `true`, the parser verified that all products between defined types
    /// have matching output types.
    pub complete: bool,
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

    /// Returns indices of basis vectors that square to +1 (positive metric).
    pub fn positive_indices(&self) -> impl Iterator<Item = usize> + '_ {
        self.basis.iter().filter(|b| b.metric == 1).map(|b| b.index)
    }

    /// Returns indices of basis vectors that square to -1 (negative metric).
    pub fn negative_indices(&self) -> impl Iterator<Item = usize> + '_ {
        self.basis
            .iter()
            .filter(|b| b.metric == -1)
            .map(|b| b.index)
    }

    /// Returns indices of basis vectors that square to 0 (degenerate/null).
    ///
    /// In PGA, this is typically the e0 basis (projective origin).
    pub fn degenerate_indices(&self) -> impl Iterator<Item = usize> + '_ {
        self.basis.iter().filter(|b| b.metric == 0).map(|b| b.index)
    }

    /// Returns true if this algebra has a degenerate metric (r > 0).
    ///
    /// Algebras with degenerate metrics (like PGA) have different normalization
    /// behavior and require bulk/weight decomposition.
    #[inline]
    pub fn is_degenerate(&self) -> bool {
        self.r > 0
    }

    /// Returns true if this algebra has an indefinite metric (q > 0).
    ///
    /// Algebras with indefinite metrics (like Minkowski) can have timelike,
    /// spacelike, and lightlike vectors.
    #[inline]
    pub fn is_indefinite(&self) -> bool {
        self.q > 0
    }

    /// Returns the signature type name derived from (p, q, r).
    ///
    /// This generates a generic name like `Cl3_0_1` instead of algebra-specific
    /// names like "Projective3".
    pub fn signature_type_name(&self) -> String {
        format!("Cl{}_{}{}", self.p, self.q, self.r)
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
    /// Versor information (if this type is a versor).
    ///
    /// If present, this type can transform other elements via the sandwich
    /// product: `X' = V * X * rev(V)`.
    pub versor: Option<VersorSpec>,
    /// Whether this type is sparse (uses only a subset of blades within its grades).
    ///
    /// Sparse types have explicit blade mappings and don't use all blades of their grades.
    /// For example, a Line in CGA uses only 6 of the 10 grade-3 blades.
    pub is_sparse: bool,
    /// Types that can be transformed via inverse sandwich product.
    ///
    /// This allows non-versor types (like Circle in CGA) to perform
    /// inverse sandwich transformations: `X' = T * X * T⁻¹`.
    ///
    /// For versors, this is typically empty (uses auto-inferred targets).
    /// For blades like Circle, this explicitly lists valid targets.
    pub inverse_sandwich_targets: Vec<String>,
}

/// Versor specification for a type.
///
/// Versors are elements that can transform other elements via the sandwich
/// product. This includes rotors (rotations), motors (rigid motions),
/// and flectors (reflections).
#[derive(Debug, Clone)]
pub struct VersorSpec {
    /// Whether this versor has unit norm (`V * rev(V) = 1`).
    ///
    /// Unit versors have simpler sandwich formulas since `V⁻¹ = rev(V)`.
    pub is_unit: bool,
    /// Types that this versor can transform via sandwich product.
    ///
    /// Empty means auto-detect based on grade compatibility.
    pub sandwich_targets: Vec<String>,
}

/// A field in a type.
#[derive(Debug, Clone)]
pub struct FieldSpec {
    /// Field name (e.g., "x", "xy").
    pub name: String,
    /// Blade index this field holds (canonical bitmask form).
    pub blade_index: usize,
    /// Grade of the blade.
    pub grade: usize,
    /// Sign relative to canonical blade ordering (+1 or -1).
    ///
    /// When a blade is specified in non-canonical order (e.g., "e20" instead of "e02"),
    /// this sign captures the parity of the permutation needed to reach canonical form.
    /// For example, e20 = -e02, so sign = -1.
    ///
    /// This sign is applied during product generation to ensure correct results.
    pub sign: i8,
}

/// Product specifications for all product types.
///
/// Product naming follows [Rigid Geometric Algebra](https://rigidgeometricalgebra.org/) conventions:
/// - `∧` = wedge (exterior product)
/// - `∨` = antiwedge (regressive product)
/// - `★` = dual (bulk dual)
/// - `☆` = antidual (weight dual)
/// - `•` = dot (metric inner product, same-grade only)
/// - `⊚` = antidot (metric antiproduct inner, same-antigrade only)
#[derive(Debug, Clone, Default)]
pub struct ProductsSpec {
    /// Geometric product entries (used for Mul operator on versors).
    pub geometric: Vec<ProductEntry>,
    /// Wedge product entries (∧, exterior, grade-raising).
    pub wedge: Vec<ProductEntry>,
    /// Left contraction entries.
    pub left_contraction: Vec<ProductEntry>,
    /// Right contraction entries.
    pub right_contraction: Vec<ProductEntry>,
    /// Antiwedge product entries (∨, regressive/meet).
    pub antiwedge: Vec<ProductEntry>,
    /// Scalar product entries.
    pub scalar: Vec<ProductEntry>,
    /// Antigeometric product entries (used for antisandwich computation).
    pub antigeometric: Vec<ProductEntry>,
    /// Antiscalar product entries.
    pub antiscalar: Vec<ProductEntry>,
    /// Bulk contraction entries (a ∨ b★).
    pub bulk_contraction: Vec<ProductEntry>,
    /// Weight contraction entries (a ∨ b☆).
    pub weight_contraction: Vec<ProductEntry>,
    /// Bulk expansion entries (a ∧ b★).
    pub bulk_expansion: Vec<ProductEntry>,
    /// Weight expansion entries (a ∧ b☆).
    pub weight_expansion: Vec<ProductEntry>,
    /// Dot product entries (•, metric inner, same-grade only, returns scalar).
    pub dot: Vec<ProductEntry>,
    /// Antidot product entries (⊚, metric antiproduct inner, same-antigrade only, returns scalar).
    pub antidot: Vec<ProductEntry>,
    /// Projection entries: b ∨ (a ∧ b☆).
    pub project: Vec<ProductEntry>,
    /// Antiprojection entries: b ∧ (a ∨ b☆).
    pub antiproject: Vec<ProductEntry>,
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

/// Wrapper constraint kinds for constraint simplification.
///
/// These represent the different normalization and constraint wrappers
/// that can be applied to geometric types. Each wrapper has specific
/// algebraic constraints that can be used during Groebner basis simplification.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum WrapperKind {
    /// Unit wrapper: `norm() == 1` (Euclidean norm).
    ///
    /// Constraint: `norm_squared - 1 = 0`
    Unit,
    /// Bulk wrapper: `bulk_norm() == 1` (PGA versors).
    ///
    /// Constraint: `bulk_norm_squared - 1 = 0`
    Bulk,
    /// Unitized wrapper: `weight_norm() == 1` (PGA standard form).
    ///
    /// Constraint: `weight_norm_squared - 1 = 0`
    Unitized,
    /// Ideal wrapper: `weight_norm() ≈ 0` (PGA elements at infinity).
    ///
    /// Constraint: each weight component = 0
    Ideal,
    /// Proper wrapper: timelike, `|norm²| == 1` (Minkowski 4-velocities).
    ///
    /// Constraint: `norm_squared - 1 = 0`
    Proper,
    /// Spacelike wrapper: spacelike, `|norm²| == 1` (Minkowski spatial).
    ///
    /// Constraint: `norm_squared + 1 = 0`
    Spacelike,
    /// Null wrapper: `norm_squared ≈ 0` (Minkowski lightlike).
    ///
    /// Constraint: `norm_squared = 0`
    Null,
}
