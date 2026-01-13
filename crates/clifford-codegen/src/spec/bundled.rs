//! Bundled algebra specifications.
//!
//! This module provides pre-defined specifications for common algebras.

/// 2D Euclidean algebra specification (TOML).
pub const EUCLIDEAN2: &str = include_str!("../../../../algebras/euclidean2.toml");

/// 3D Euclidean algebra specification (TOML).
pub const EUCLIDEAN3: &str = include_str!("../../../../algebras/euclidean3.toml");

#[cfg(test)]
mod tests {
    use super::*;
    use crate::spec::parse_spec;

    #[test]
    fn parse_euclidean2_spec() {
        let spec = parse_spec(EUCLIDEAN2).unwrap();

        assert_eq!(spec.name, "euclidean2");
        assert_eq!(spec.module_path, Some("euclidean::dim2".to_string()));
        assert_eq!(spec.signature.p, 2);
        assert_eq!(spec.signature.q, 0);
        assert_eq!(spec.signature.r, 0);
        assert_eq!(spec.signature.dim(), 2);
        assert_eq!(spec.signature.num_blades(), 4);

        // Check types
        assert!(spec.types.iter().any(|t| t.name == "Scalar"));
        assert!(spec.types.iter().any(|t| t.name == "Vector"));
        assert!(spec.types.iter().any(|t| t.name == "Bivector"));
        assert!(spec.types.iter().any(|t| t.name == "Rotor"));

        // Check Rotor
        let rotor = spec.types.iter().find(|t| t.name == "Rotor").unwrap();
        assert_eq!(rotor.grades, vec![0, 2]);
        assert_eq!(rotor.fields.len(), 2);

        // Check options
        assert!(spec.options.generate_serde);
        assert!(spec.options.generate_arbitrary);
    }

    #[test]
    fn parse_euclidean3_spec() {
        let spec = parse_spec(EUCLIDEAN3).unwrap();

        assert_eq!(spec.name, "euclidean3");
        assert_eq!(spec.module_path, Some("euclidean::dim3".to_string()));
        assert_eq!(spec.signature.p, 3);
        assert_eq!(spec.signature.q, 0);
        assert_eq!(spec.signature.r, 0);
        assert_eq!(spec.signature.dim(), 3);
        assert_eq!(spec.signature.num_blades(), 8);

        // Check types
        assert!(spec.types.iter().any(|t| t.name == "Scalar"));
        assert!(spec.types.iter().any(|t| t.name == "Vector"));
        assert!(spec.types.iter().any(|t| t.name == "Bivector"));
        assert!(spec.types.iter().any(|t| t.name == "Trivector"));
        assert!(spec.types.iter().any(|t| t.name == "Rotor"));

        // Check Vector
        let vector = spec.types.iter().find(|t| t.name == "Vector").unwrap();
        assert_eq!(vector.grades, vec![1]);
        assert_eq!(vector.fields.len(), 3);

        // Check Bivector
        let bivector = spec.types.iter().find(|t| t.name == "Bivector").unwrap();
        assert_eq!(bivector.grades, vec![2]);
        assert_eq!(bivector.fields.len(), 3);

        // Check Rotor
        let rotor = spec.types.iter().find(|t| t.name == "Rotor").unwrap();
        assert_eq!(rotor.grades, vec![0, 2]);
        assert_eq!(rotor.fields.len(), 4);
    }
}
