//! Discovered entity representation.

/// Represents a discovered geometric entity.
///
/// An entity is a valid grade combination with associated metadata.
/// Valid entities must satisfy two constraints:
///
/// 1. **Geometric Product Constraint**: `u * ũ = scalar`
/// 2. **Antiproduct Constraint**: `u ⊟ ũ̃ = antiscalar`
///
/// Each constraint may require a field constraint expression that must equal zero.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct DiscoveredEntity {
    /// Suggested name for this entity (e.g., "Entity_0_2").
    pub name: String,

    /// The grades present in this entity.
    pub grades: Vec<usize>,

    /// Field constraint for the geometric product constraint (`u * ũ = scalar`).
    ///
    /// If `None`, the entity automatically satisfies the geometric constraint.
    /// If `Some(expr)`, the expression must equal zero for the constraint to hold.
    /// Example: "e01*e23 + e02*e31 + e03*e12 = 0" for PGA bivectors.
    pub geometric_constraint: Option<String>,

    /// Field constraint for the antiproduct constraint (`u ⊟ ũ̃ = antiscalar`).
    ///
    /// If `None`, the entity automatically satisfies the antiproduct constraint.
    /// If `Some(expr)`, the expression must equal zero for the constraint to hold.
    pub antiproduct_constraint: Option<String>,
}

impl DiscoveredEntity {
    /// Creates a new discovered entity with minimal information.
    pub fn new(name: impl Into<String>, grades: Vec<usize>) -> Self {
        Self {
            name: name.into(),
            grades,
            geometric_constraint: None,
            antiproduct_constraint: None,
        }
    }

    /// Sets the geometric field constraint expression.
    pub fn with_geometric_constraint(mut self, constraint: impl Into<String>) -> Self {
        self.geometric_constraint = Some(constraint.into());
        self
    }

    /// Sets the antiproduct field constraint expression.
    pub fn with_antiproduct_constraint(mut self, constraint: impl Into<String>) -> Self {
        self.antiproduct_constraint = Some(constraint.into());
        self
    }

    /// Returns the total number of constraints this entity has.
    pub fn constraint_count(&self) -> usize {
        let mut count = 0;
        if self.geometric_constraint.is_some() {
            count += 1;
        }
        if self.antiproduct_constraint.is_some() {
            count += 1;
        }
        count
    }

    /// Returns true if this entity has any field constraints.
    pub fn has_constraints(&self) -> bool {
        self.geometric_constraint.is_some() || self.antiproduct_constraint.is_some()
    }
}
