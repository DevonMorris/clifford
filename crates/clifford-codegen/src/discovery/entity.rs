//! Discovered entity representation.

/// Represents a discovered geometric entity.
///
/// An entity is a valid grade combination with associated metadata.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct DiscoveredEntity {
    /// Suggested name for this entity (e.g., "Entity_0_2").
    pub name: String,

    /// The grades present in this entity.
    pub grades: Vec<usize>,

    /// Field constraint expression required for geometric constraint satisfaction.
    ///
    /// If `None`, the entity automatically satisfies geometric constraints.
    /// If `Some(expr)`, the expression must equal zero for the constraint to hold.
    /// Example: "e01*e23 + e02*e31 + e03*e12 = 0" for PGA bivectors.
    pub constraint: Option<String>,
}

impl DiscoveredEntity {
    /// Creates a new discovered entity with minimal information.
    pub fn new(name: impl Into<String>, grades: Vec<usize>) -> Self {
        Self {
            name: name.into(),
            grades,
            constraint: None,
        }
    }

    /// Sets the field constraint expression.
    pub fn with_constraint(mut self, constraint: impl Into<String>) -> Self {
        self.constraint = Some(constraint.into());
        self
    }
}
