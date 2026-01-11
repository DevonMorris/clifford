//! Discovered entity representation.

/// Represents a discovered geometric entity.
///
/// An entity is a valid grade combination with associated metadata
/// including a suggested name and information about applicable constraints.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct DiscoveredEntity {
    /// Suggested name for this entity (e.g., "Vector", "Rotor").
    pub name: String,

    /// The grades present in this entity.
    pub grades: Vec<usize>,

    /// Human-readable description.
    pub description: String,

    /// Whether this type can support a unit constraint.
    /// True if `u * reverse(u)` produces a scalar (norm exists).
    pub can_be_unit: bool,

    /// Whether this type can support a nonzero constraint.
    /// True if normalization is meaningful.
    pub can_be_nonzero: bool,
}

impl DiscoveredEntity {
    /// Creates a new discovered entity with minimal information.
    pub fn new(name: impl Into<String>, grades: Vec<usize>) -> Self {
        Self {
            name: name.into(),
            grades,
            description: String::new(),
            can_be_unit: false,
            can_be_nonzero: false,
        }
    }

    /// Sets the description.
    pub fn with_description(mut self, description: impl Into<String>) -> Self {
        self.description = description.into();
        self
    }

    /// Marks this entity as supporting unit constraint.
    pub fn with_unit_support(mut self) -> Self {
        self.can_be_unit = true;
        self
    }

    /// Marks this entity as supporting nonzero constraint.
    pub fn with_nonzero_support(mut self) -> Self {
        self.can_be_nonzero = true;
        self
    }
}
