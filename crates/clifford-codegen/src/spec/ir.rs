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
    /// Constraints on this type (both algebraic and user-defined).
    ///
    /// Constraints with `solve_for` reduce the parameter count in constructors.
    /// Linear constraints are always solvable; quadratic constraints may require
    /// `Option<Self>` return type due to domain restrictions.
    pub constraints: Vec<UserConstraint>,
    /// Versor information (if this type is a versor).
    ///
    /// If present, this type can transform other elements via the sandwich
    /// product: `X' = V * X * rev(V)`.
    pub versor: Option<VersorSpec>,
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

/// Sign convention for constraints with ± ambiguity.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum SignConvention {
    /// Use positive root (default).
    #[default]
    Positive,
    /// Use negative root.
    Negative,
}

/// User-defined constraint on a type.
#[derive(Debug, Clone)]
pub struct UserConstraint {
    /// Constraint name (e.g., "unit", "normalized").
    pub name: String,
    /// Documentation.
    pub description: Option<String>,
    /// Constraint expression (e.g., "s*s + e12*e12 + e13*e13 + e23*e23 = 1").
    pub expression: String,
    /// Field to solve for. If present, this constraint reduces the parameter count.
    pub solve_for: Option<String>,
    /// Sign convention for ± ambiguity.
    pub sign: SignConvention,
    /// Enforcement method to generate (e.g., "normalize").
    pub enforce: Option<String>,
    /// Whether this constraint has domain restrictions (e.g., requires sqrt).
    /// If true, `new()` returns `Option<Self>`.
    pub has_domain_restriction: bool,
}

impl TypeSpec {
    /// Returns true if this type has any constraints.
    pub fn has_constraint(&self) -> bool {
        !self.constraints.is_empty()
    }

    /// Returns all solve_for fields that need to be computed.
    pub fn solve_for_fields(&self) -> Vec<&str> {
        self.constraints
            .iter()
            .filter_map(|c| c.solve_for.as_deref())
            .collect()
    }

    /// Returns true if any solvable constraint has domain restrictions.
    ///
    /// If true, `new()` should return `Option<Self>` instead of `Self`.
    pub fn has_domain_restrictions(&self) -> bool {
        self.constraints
            .iter()
            .any(|c| c.solve_for.is_some() && c.has_domain_restriction)
    }

    /// Returns all user-defined constraints that should be solved.
    pub fn solvable_user_constraints(&self) -> impl Iterator<Item = &UserConstraint> {
        self.constraints.iter().filter(|c| c.solve_for.is_some())
    }

    /// Returns all user-defined constraints with enforcement methods.
    pub fn enforceable_constraints(&self) -> impl Iterator<Item = &UserConstraint> {
        self.constraints.iter().filter(|c| c.enforce.is_some())
    }
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
