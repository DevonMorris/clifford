# PRD-14.10: Product Output Inference

**Status**: Draft
**Parent**: [PRD-14: Geometric Algebra Code Generator](prd-14-codegen.md)
**Goal**: Automatically infer product output types and constraint preservation from discovered entities

## Overview

Given a set of discovered entities (PRD-14.9), this module:
1. Computes the output grades for each product type (geometric, outer, inner, etc.)
2. Matches output grades to known entity types
3. Identifies which operations preserve constraints (unit, nonzero)
4. Auto-populates the `[products.*]` sections of the TOML specification

This eliminates the need for manual product specification, reducing errors and effort.

## Product Grade Computation

### Grade Rules by Product Type

| Product | Input Grades | Output Grades | Notes |
|---------|--------------|---------------|-------|
| Geometric | (a, b) | \|a-b\| to a+b (step 2) | Full geometric product |
| Outer (∧) | (a, b) | a + b | Zero if > dim |
| Inner (·) | (a, b) | \|a - b\| | Symmetric inner |
| Left Contraction (⌋) | (a, b) | b - a | Zero if a > b |
| Right Contraction (⌊) | (a, b) | a - b | Zero if b > a |
| Scalar | (a, b) | 0 | Grade-0 extraction |
| Regressive (∨) | (a, b) | a + b - n | Dual of outer |

### Computing Output for Multi-Grade Types

For types with multiple grades, the output is the union of outputs for each grade pair:

```rust
use std::collections::HashSet;

/// Computes output grades for a product of two multi-grade types.
fn compute_product_output_grades(
    lhs_grades: &[usize],
    rhs_grades: &[usize],
    product_type: ProductType,
    dim: usize,
) -> Vec<usize> {
    let mut output_grades = HashSet::new();

    for &ga in lhs_grades {
        for &gb in rhs_grades {
            let grades = match product_type {
                ProductType::Geometric => {
                    geometric_grades(ga, gb, dim)
                }
                ProductType::Outer => {
                    outer_grade(ga, gb, dim).map(|g| vec![g]).unwrap_or_default()
                }
                ProductType::Inner => {
                    vec![inner_grade(ga, gb)]
                }
                ProductType::LeftContraction => {
                    left_contraction_grade(ga, gb).map(|g| vec![g]).unwrap_or_default()
                }
                ProductType::Scalar => {
                    if ga == gb { vec![0] } else { vec![] }
                }
                ProductType::Regressive => {
                    let sum = ga + gb;
                    if sum >= dim { vec![sum - dim] } else { vec![] }
                }
            };

            output_grades.extend(grades);
        }
    }

    let mut result: Vec<_> = output_grades.into_iter().collect();
    result.sort();
    result
}
```

## Entity Matching

Match computed output grades to known entities:

```rust
/// Finds the entity that matches the given output grades.
fn find_matching_entity<'a>(
    output_grades: &[usize],
    entities: &'a [DiscoveredEntity],
) -> Option<&'a DiscoveredEntity> {
    entities.iter().find(|e| e.grades == output_grades)
}

/// Result of product inference.
#[derive(Debug, Clone)]
pub enum ProductOutput {
    /// Output matches a known entity
    Entity(String),
    /// Output is a scalar (grade 0)
    Scalar,
    /// Output doesn't match any entity (grades not in type system)
    NoMatch(Vec<usize>),
    /// Product is always zero
    Zero,
}

/// Infers the output type for a product.
fn infer_product_output(
    lhs: &DiscoveredEntity,
    rhs: &DiscoveredEntity,
    product_type: ProductType,
    entities: &[DiscoveredEntity],
    dim: usize,
) -> ProductOutput {
    let output_grades = compute_product_output_grades(
        &lhs.grades,
        &rhs.grades,
        product_type,
        dim,
    );

    if output_grades.is_empty() {
        return ProductOutput::Zero;
    }

    if output_grades == vec![0] {
        // Check if we have a Scalar type
        if let Some(entity) = find_matching_entity(&[0], entities) {
            return ProductOutput::Entity(entity.name.clone());
        }
        return ProductOutput::Scalar;
    }

    match find_matching_entity(&output_grades, entities) {
        Some(entity) => ProductOutput::Entity(entity.name.clone()),
        None => ProductOutput::NoMatch(output_grades),
    }
}
```

## Complete Product Table Generation

Generate all products between all entity pairs:

```rust
/// A single product entry.
#[derive(Debug, Clone)]
pub struct ProductEntry {
    pub lhs: String,
    pub rhs: String,
    pub output: ProductOutput,
}

/// Generates all products for a given product type.
fn generate_product_table(
    entities: &[DiscoveredEntity],
    product_type: ProductType,
    dim: usize,
) -> Vec<ProductEntry> {
    let mut entries = Vec::new();

    for lhs in entities {
        for rhs in entities {
            let output = infer_product_output(lhs, rhs, product_type, entities, dim);

            // Skip zero products
            if matches!(output, ProductOutput::Zero) {
                continue;
            }

            entries.push(ProductEntry {
                lhs: lhs.name.clone(),
                rhs: rhs.name.clone(),
                output,
            });
        }
    }

    entries
}

/// Generates all product tables.
fn generate_all_product_tables(
    entities: &[DiscoveredEntity],
    dim: usize,
) -> HashMap<ProductType, Vec<ProductEntry>> {
    let mut tables = HashMap::new();

    for product_type in [
        ProductType::Geometric,
        ProductType::Outer,
        ProductType::Inner,
        ProductType::LeftContraction,
        ProductType::Scalar,
    ] {
        tables.insert(product_type, generate_product_table(entities, product_type, dim));
    }

    tables
}
```

## Constraint Preservation Analysis

Analyze which products preserve unit and nonzero constraints:

### Unit Constraint Preservation

```rust
/// Checks if a product preserves unit constraint.
///
/// A product `a * b` preserves unit constraint if:
/// 1. Both inputs are the same type
/// 2. Output is the same type
/// 3. The algebraic property ||a * b|| = ||a|| * ||b|| holds
fn product_preserves_unit(
    lhs: &DiscoveredEntity,
    rhs: &DiscoveredEntity,
    output: &ProductOutput,
) -> bool {
    // Same type composition (e.g., UnitRotor * UnitRotor = UnitRotor)
    if lhs.name == rhs.name {
        if let ProductOutput::Entity(ref out_name) = output {
            if out_name == &lhs.name {
                // Check if this is a known unit-preserving pattern
                // Rotors, motors, etc. preserve unit under composition
                return is_unit_preserving_type(&lhs.name);
            }
        }
    }

    false
}

fn is_unit_preserving_type(name: &str) -> bool {
    matches!(
        name.to_lowercase().as_str(),
        "rotor" | "motor" | "flector" | "even" | "unitrotor" | "unitmotor"
    )
}

/// Collects all unit-preserving operations for a type.
fn collect_unit_preserving_ops(
    entity: &DiscoveredEntity,
    product_tables: &HashMap<ProductType, Vec<ProductEntry>>,
) -> Vec<String> {
    let mut ops = Vec::new();

    // Check geometric product
    if let Some(entries) = product_tables.get(&ProductType::Geometric) {
        for entry in entries {
            if entry.lhs == entity.name && entry.rhs == entity.name {
                if product_preserves_unit(entity, entity, &entry.output) {
                    ops.push("compose".to_string());
                }
            }
        }
    }

    // Reverse always preserves unit
    if entity.can_be_unit {
        ops.push("reverse".to_string());
        ops.push("inverse".to_string());
    }

    ops
}
```

### NonZero Constraint Preservation

```rust
/// Checks if a product preserves nonzero constraint.
///
/// A product `a * b` preserves nonzero if neither input can make the result zero.
fn product_preserves_nonzero(
    lhs: &DiscoveredEntity,
    rhs: &DiscoveredEntity,
    output: &ProductOutput,
    algebra: &Algebra,
) -> bool {
    // In Euclidean algebras, products of nonzero elements are often nonzero
    // But in degenerate algebras, this isn't guaranteed

    // For now, only guarantee nonzero preservation for:
    // 1. Scalar multiplication by nonzero scalar
    // 2. Same-type products in non-degenerate algebras

    if algebra.is_degenerate() {
        return false;
    }

    matches!(output, ProductOutput::Entity(_))
}
```

## TOML Generation

Generate the `[products.*]` sections:

```rust
fn write_products_section(
    product_tables: &HashMap<ProductType, Vec<ProductEntry>>,
    output: &mut dyn Write,
) -> std::io::Result<()> {
    writeln!(output, "# ============================================================")?;
    writeln!(output, "# Products (auto-inferred)")?;
    writeln!(output, "# ============================================================")?;
    writeln!(output)?;

    // Geometric products
    writeln!(output, "[products.geometric]")?;
    if let Some(entries) = product_tables.get(&ProductType::Geometric) {
        for entry in entries {
            if let ProductOutput::Entity(ref out_name) = entry.output {
                writeln!(output, "{}_{} = \"{}\"", entry.lhs, entry.rhs, out_name)?;
            }
        }
    }
    writeln!(output)?;

    // Outer products
    writeln!(output, "[products.outer]")?;
    if let Some(entries) = product_tables.get(&ProductType::Outer) {
        for entry in entries {
            if let ProductOutput::Entity(ref out_name) = entry.output {
                writeln!(output, "{}_{} = \"{}\"", entry.lhs, entry.rhs, out_name)?;
            }
        }
    }
    writeln!(output)?;

    // Scalar products
    writeln!(output, "[products.scalar]")?;
    if let Some(entries) = product_tables.get(&ProductType::Scalar) {
        for entry in entries {
            match &entry.output {
                ProductOutput::Entity(name) => {
                    writeln!(output, "{}_{} = \"{}\"", entry.lhs, entry.rhs, name)?;
                }
                ProductOutput::Scalar => {
                    writeln!(output, "{}_{} = \"T\"", entry.lhs, entry.rhs)?;
                }
                _ => {}
            }
        }
    }
    writeln!(output)?;

    // Left contraction
    writeln!(output, "[products.left_contraction]")?;
    if let Some(entries) = product_tables.get(&ProductType::LeftContraction) {
        for entry in entries {
            if let ProductOutput::Entity(ref out_name) = entry.output {
                writeln!(output, "{}_{} = \"{}\"", entry.lhs, entry.rhs, out_name)?;
            }
        }
    }

    Ok(())
}
```

## Constraint Wrapper Suggestions

Add constraint suggestions based on analysis:

```rust
fn write_constraint_suggestions(
    entity: &DiscoveredEntity,
    unit_preserving_ops: &[String],
    output: &mut dyn Write,
) -> std::io::Result<()> {
    if entity.can_be_unit || entity.can_be_nonzero {
        writeln!(output)?;
        writeln!(output, "# Suggested constraints for {}:", entity.name)?;
    }

    if entity.can_be_unit {
        writeln!(output, "# [types.{}.constraints.unit]", entity.name)?;
        if !unit_preserving_ops.is_empty() {
            writeln!(
                output,
                "# preserving_ops = {:?}",
                unit_preserving_ops
            )?;
        }
    }

    if entity.can_be_nonzero {
        writeln!(output, "# [types.{}.constraints.nonzero]", entity.name)?;
    }

    Ok(())
}
```

## Complete Integration

Integrate with entity discovery for full template generation:

```rust
/// Complete template generation with products and constraints.
pub fn generate_complete_template(
    algebra: &Algebra,
    output: &mut dyn Write,
) -> std::io::Result<()> {
    let table = ProductTable::new(algebra);

    // Discover entities (PRD-14.9)
    let mut entities = discover_entities(algebra, &table);

    // Generate product tables
    let product_tables = generate_all_product_tables(&entities, algebra.dim());

    // Analyze constraint preservation
    for entity in &mut entities {
        entity.unit_preserving_ops = collect_unit_preserving_ops(entity, &product_tables);
    }

    // Write complete TOML
    write_header(algebra, output)?;
    write_signature(algebra, output)?;
    write_types_with_constraints(&entities, output)?;
    write_products_section(&product_tables, output)?;
    write_options(output)?;

    Ok(())
}
```

## Example Output

For 3D Euclidean algebra:

```toml
# Auto-generated by clifford-codegen discover
# Algebra: Cl(3,0,0) - 3D Euclidean

[algebra]
name = "euclidean3"
module_path = "euclidean::dim3"
description = "3D Euclidean Geometric Algebra"

[signature]
positive = ["e1", "e2", "e3"]
negative = []
zero = []

# ============================================================
# Types
# ============================================================

[types.Scalar]
grades = [0]
description = "Grade-0 scalar"

[types.Vector]
grades = [1]
description = "Grade-1 vector"

# Suggested constraints:
# [types.Vector.constraints.unit]
# [types.Vector.constraints.nonzero]

[types.Bivector]
grades = [2]
description = "Grade-2 bivector"

# Suggested constraints:
# [types.Bivector.constraints.unit]

[types.Trivector]
grades = [3]
description = "Grade-3 pseudoscalar"

[types.Rotor]
grades = [0, 2]
description = "Even subalgebra (rotation)"

# Suggested constraints:
# [types.Rotor.constraints.unit]
# preserving_ops = ["compose", "reverse", "inverse"]

[types.Odd]
grades = [1, 3]
description = "Odd subalgebra"

[types.Full]
grades = [0, 1, 2, 3]
description = "Full multivector"

# ============================================================
# Products (auto-inferred)
# ============================================================

[products.geometric]
Vector_Vector = "Rotor"
Bivector_Vector = "Odd"
Vector_Bivector = "Odd"
Rotor_Rotor = "Rotor"
Rotor_Vector = "Odd"
Vector_Rotor = "Odd"
# ... all valid products

[products.outer]
Vector_Vector = "Bivector"
Vector_Bivector = "Trivector"
Bivector_Vector = "Trivector"
Scalar_Vector = "Vector"
# ... all valid products

[products.scalar]
Vector_Vector = "T"
Bivector_Bivector = "T"
Rotor_Rotor = "T"
# ... all scalar products

[products.left_contraction]
Vector_Bivector = "Vector"
Vector_Trivector = "Bivector"
Bivector_Trivector = "Vector"
# ... all contractions

[options]
generate_serde = true
generate_arbitrary = true
generate_nalgebra = true
```

## Module Structure

```
crates/clifford-codegen/src/discovery/
├── mod.rs
├── entity.rs         # From PRD-14.9
├── naming.rs         # From PRD-14.9
├── products.rs       # NEW: Product inference
├── constraints.rs    # NEW: Constraint analysis
└── template.rs       # Full template generation
```

## Testing

### Unit Tests

```rust
#[test]
fn geometric_product_inference() {
    let algebra = Algebra::euclidean(3);
    let entities = vec![
        entity("Vector", &[1]),
        entity("Bivector", &[2]),
        entity("Rotor", &[0, 2]),
    ];

    let output = infer_product_output(
        &entities[0], // Vector
        &entities[0], // Vector
        ProductType::Geometric,
        &entities,
        3,
    );

    assert!(matches!(output, ProductOutput::Entity(ref n) if n == "Rotor"));
}

#[test]
fn outer_product_inference() {
    let algebra = Algebra::euclidean(3);
    let entities = vec![
        entity("Vector", &[1]),
        entity("Bivector", &[2]),
        entity("Trivector", &[3]),
    ];

    let output = infer_product_output(
        &entities[0], // Vector
        &entities[1], // Bivector
        ProductType::Outer,
        &entities,
        3,
    );

    assert!(matches!(output, ProductOutput::Entity(ref n) if n == "Trivector"));
}

#[test]
fn unit_preservation_detection() {
    let rotor = entity("Rotor", &[0, 2]);
    rotor.can_be_unit = true;

    let output = ProductOutput::Entity("Rotor".to_string());

    assert!(product_preserves_unit(&rotor, &rotor, &output));
}
```

### Integration Tests

```rust
#[test]
fn complete_euclidean3_template() {
    let algebra = Algebra::euclidean(3);
    let mut output = Vec::new();

    generate_complete_template(&algebra, &mut output).unwrap();

    let toml_str = String::from_utf8(output).unwrap();
    let parsed: toml::Value = toml::from_str(&toml_str).unwrap();

    // Verify products section exists
    assert!(parsed.get("products").is_some());
    let products = parsed["products"].as_table().unwrap();
    assert!(products.contains_key("geometric"));
    assert!(products.contains_key("outer"));

    // Verify specific products
    let geometric = products["geometric"].as_table().unwrap();
    assert_eq!(geometric["Vector_Vector"].as_str(), Some("Rotor"));
}
```

## Deliverables

- [ ] `ProductType` enum for all product types
- [ ] `ProductOutput` enum for inference results
- [ ] `compute_product_output_grades(lhs, rhs, product_type, dim)`
- [ ] `infer_product_output(lhs, rhs, product_type, entities, dim)`
- [ ] `generate_product_table(entities, product_type, dim)`
- [ ] `generate_all_product_tables(entities, dim)`
- [ ] `product_preserves_unit(lhs, rhs, output)`
- [ ] `collect_unit_preserving_ops(entity, product_tables)`
- [ ] `generate_complete_template(algebra, output)` integration
- [ ] TOML generation for all product sections
- [ ] Tests for all standard algebras

## Success Criteria

1. **Correctness**: Inferred products match manual specification
2. **Completeness**: All valid products are discovered
3. **Constraint Analysis**: Correctly identifies unit-preserving operations
4. **Usability**: Generated TOML is immediately usable

## Dependencies

- PRD-14.8 (Geometric Constraint Engine)
- PRD-14.9 (Entity Discovery)
- PRD-14.1 (Blade Algebra Engine)

## References

- [Geometric Product - RGA Wiki](https://rigidgeometricalgebra.org/wiki/index.php?title=Geometric_product)
- [Products Summary - RGA Wiki](https://rigidgeometricalgebra.org/wiki/index.php?title=Products)
