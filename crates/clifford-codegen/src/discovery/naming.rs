//! Entity naming.
//!
//! This module provides simple naming for grade combinations.
//! All entities are named "Entity_" followed by their grades.

/// Generates a name for a grade combination.
///
/// Names are simply "Entity_" followed by the sorted grades joined by underscores.
/// No heuristics or special names - the user can rename as desired.
///
/// # Arguments
///
/// * `grades` - The grades present in the entity
/// * `_dim` - The dimension of the algebra (unused, kept for API compatibility)
///
/// # Example
///
/// ```
/// use clifford_codegen::discovery::suggest_name;
///
/// assert_eq!(suggest_name(&[0], 3), "Entity_0");
/// assert_eq!(suggest_name(&[1], 3), "Entity_1");
/// assert_eq!(suggest_name(&[0, 2], 3), "Entity_0_2");
/// assert_eq!(suggest_name(&[0, 1, 2, 3], 3), "Entity_0_1_2_3");
/// ```
pub fn suggest_name(grades: &[usize], _dim: usize) -> String {
    // Sort grades for consistent naming
    let mut sorted = grades.to_vec();
    sorted.sort();

    format!(
        "Entity_{}",
        sorted
            .iter()
            .map(|g| g.to_string())
            .collect::<Vec<_>>()
            .join("_")
    )
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn single_grades() {
        assert_eq!(suggest_name(&[0], 3), "Entity_0");
        assert_eq!(suggest_name(&[1], 3), "Entity_1");
        assert_eq!(suggest_name(&[2], 3), "Entity_2");
        assert_eq!(suggest_name(&[3], 3), "Entity_3");
    }

    #[test]
    fn multi_grades() {
        assert_eq!(suggest_name(&[0, 2], 3), "Entity_0_2");
        assert_eq!(suggest_name(&[1, 3], 3), "Entity_1_3");
        assert_eq!(suggest_name(&[0, 1, 2, 3], 3), "Entity_0_1_2_3");
    }

    #[test]
    fn handles_unsorted_input() {
        // Should work even if input isn't sorted
        assert_eq!(suggest_name(&[2, 0], 3), "Entity_0_2");
        assert_eq!(suggest_name(&[3, 1], 3), "Entity_1_3");
        assert_eq!(suggest_name(&[3, 1, 2, 0], 3), "Entity_0_1_2_3");
    }
}
