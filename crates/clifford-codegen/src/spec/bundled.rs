//! Bundled algebra specifications.
//!
//! This module provides pre-defined specifications for common algebras.

/// 2D Euclidean algebra specification (TOML).
pub const EUCLIDEAN2: &str = include_str!("../../algebras/euclidean2.toml");

/// 3D Euclidean algebra specification (TOML).
pub const EUCLIDEAN3: &str = include_str!("../../algebras/euclidean3.toml");

#[cfg(test)]
mod tests {
    use super::*;
    use crate::spec::{ConstraintKind, NormType, parse_spec};

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

        // Check blade names
        assert_eq!(spec.blade_names.get(&1), Some(&"x".to_string()));
        assert_eq!(spec.blade_names.get(&2), Some(&"y".to_string()));
        assert_eq!(spec.blade_names.get(&3), Some(&"xy".to_string()));

        // Check types
        assert!(spec.types.iter().any(|t| t.name == "Scalar"));
        assert!(spec.types.iter().any(|t| t.name == "Vector"));
        assert!(spec.types.iter().any(|t| t.name == "Bivector"));
        assert!(spec.types.iter().any(|t| t.name == "Rotor"));
        assert!(spec.types.iter().any(|t| t.name == "Even"));
        assert!(spec.types.iter().any(|t| t.name == "Full"));

        // Check Rotor
        let rotor = spec.types.iter().find(|t| t.name == "Rotor").unwrap();
        assert_eq!(rotor.grades, vec![0, 2]);
        assert_eq!(rotor.fields.len(), 2);

        // Check Rotor unit constraint
        let unit = rotor
            .constraints
            .iter()
            .find(|c| c.kind == ConstraintKind::Unit)
            .unwrap();
        assert_eq!(unit.norm_type, Some(NormType::Euclidean));
        assert_eq!(unit.wrapper_name, "UnitRotor");

        // Check Even aliases Rotor
        let even = spec.types.iter().find(|t| t.name == "Even").unwrap();
        assert_eq!(even.alias_of, Some("Rotor".to_string()));

        // Check options
        assert!(spec.options.generate_serde);
        assert!(spec.options.generate_arbitrary);
        assert!(spec.options.generate_tests);
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

        // Check blade names
        assert_eq!(spec.blade_names.get(&1), Some(&"x".to_string()));
        assert_eq!(spec.blade_names.get(&2), Some(&"y".to_string()));
        assert_eq!(spec.blade_names.get(&4), Some(&"z".to_string()));
        assert_eq!(spec.blade_names.get(&3), Some(&"xy".to_string()));
        assert_eq!(spec.blade_names.get(&5), Some(&"xz".to_string()));
        assert_eq!(spec.blade_names.get(&6), Some(&"yz".to_string()));
        assert_eq!(spec.blade_names.get(&7), Some(&"xyz".to_string()));

        // Check types
        assert!(spec.types.iter().any(|t| t.name == "Scalar"));
        assert!(spec.types.iter().any(|t| t.name == "Vector"));
        assert!(spec.types.iter().any(|t| t.name == "Bivector"));
        assert!(spec.types.iter().any(|t| t.name == "Trivector"));
        assert!(spec.types.iter().any(|t| t.name == "Rotor"));
        assert!(spec.types.iter().any(|t| t.name == "Even"));
        assert!(spec.types.iter().any(|t| t.name == "Odd"));
        assert!(spec.types.iter().any(|t| t.name == "Full"));

        // Check Vector
        let vector = spec.types.iter().find(|t| t.name == "Vector").unwrap();
        assert_eq!(vector.grades, vec![1]);
        assert_eq!(vector.fields.len(), 3);
        assert_eq!(vector.constraints.len(), 2); // unit + nonzero

        // Check Bivector
        let bivector = spec.types.iter().find(|t| t.name == "Bivector").unwrap();
        assert_eq!(bivector.grades, vec![2]);
        assert_eq!(bivector.fields.len(), 3);

        // Check Rotor
        let rotor = spec.types.iter().find(|t| t.name == "Rotor").unwrap();
        assert_eq!(rotor.grades, vec![0, 2]);
        assert_eq!(rotor.fields.len(), 4);

        // Check Rotor unit constraint constructors
        let unit = rotor
            .constraints
            .iter()
            .find(|c| c.kind == ConstraintKind::Unit)
            .unwrap();
        assert!(unit.constructors.iter().any(|c| c.name == "identity"));
        assert!(
            unit.constructors
                .iter()
                .any(|c| c.name == "from_angle_plane")
        );
        assert!(
            unit.constructors
                .iter()
                .any(|c| c.name == "from_angle_axis")
        );
        assert!(unit.constructors.iter().any(|c| c.name == "from_vectors"));

        // Check Odd type
        let odd = spec.types.iter().find(|t| t.name == "Odd").unwrap();
        assert_eq!(odd.grades, vec![1, 3]);
        assert_eq!(odd.fields.len(), 4);

        // Check Full type
        let full = spec.types.iter().find(|t| t.name == "Full").unwrap();
        assert_eq!(full.grades, vec![0, 1, 2, 3]);
        assert_eq!(full.fields.len(), 8);
    }
}
