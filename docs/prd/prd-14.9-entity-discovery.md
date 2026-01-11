# PRD-14.9: Entity Discovery

**Status**: Draft
**Parent**: [PRD-14: Geometric Algebra Code Generator](prd-14-codegen.md)
**Goal**: Automatically discover valid geometric entities from an algebra signature and generate a TOML template

## Overview

Given an algebra signature (p, q, r), this module:
1. Enumerates all possible grade combinations (2^(n+1) subsets for n-dimensional algebra)
2. Filters to those satisfying both geometric constraints (PRD-14.8)
3. Names known entities using heuristics
4. Generates a TOML template for users to customize

This enables a workflow where users provide only a signature, and the generator discovers the valid types automatically.

## Entity Discovery Algorithm

### Grade Combination Enumeration

For an n-dimensional algebra with grades 0 through n, we enumerate all non-empty subsets of grades:

```rust
/// Enumerates all non-empty grade combinations for an n-dimensional algebra.
fn enumerate_grade_combinations(dim: usize) -> Vec<Vec<usize>> {
    let num_grades = dim + 1; // grades 0 through dim

    // 2^(num_grades) - 1 non-empty subsets
    (1..(1 << num_grades))
        .map(|mask| {
            (0..num_grades)
                .filter(|&g| (mask >> g) & 1 == 1)
                .collect()
        })
        .collect()
}
```

### Constraint Filtering

Filter to combinations satisfying both constraints:

```rust
fn discover_valid_combinations(
    algebra: &Algebra,
    table: &ProductTable,
) -> Vec<Vec<usize>> {
    enumerate_grade_combinations(algebra.dim())
        .into_iter()
        .filter(|grades| satisfies_all_constraints(grades, algebra, table))
        .collect()
}
```

### Entity Naming Heuristics

Well-known grade combinations get standard names:

| Grades | Name | Description |
|--------|------|-------------|
| `[0]` | Scalar | Grade-0 element |
| `[1]` | Vector | Grade-1 element |
| `[2]` | Bivector | Grade-2 element |
| `[3]` | Trivector | Grade-3 element |
| `[n]` | Pseudoscalar | Top-grade element |
| `[0, 2]` | Even / Rotor | Even subalgebra |
| `[1, 3]` | Odd | Odd subalgebra |
| `[0, 2, 4, ...]` | Even | All even grades |
| `[1, 3, 5, ...]` | Odd | All odd grades |
| `[0, 1, ..., n]` | Full | Full multivector |

Unknown combinations get placeholder names like `Entity_0_1_3`.

```rust
/// Suggests a name for a grade combination based on known patterns.
fn suggest_name(grades: &[usize], dim: usize) -> String {
    // Sort grades for consistent matching
    let mut sorted = grades.to_vec();
    sorted.sort();

    // Single grades
    if sorted.len() == 1 {
        return match sorted[0] {
            0 => "Scalar".to_string(),
            1 => "Vector".to_string(),
            2 => "Bivector".to_string(),
            3 => "Trivector".to_string(),
            4 => "Quadvector".to_string(),
            g if g == dim => format!("Pseudoscalar"),
            g => format!("Grade{}", g),
        };
    }

    // Full multivector
    if sorted.len() == dim + 1 && sorted == (0..=dim).collect::<Vec<_>>() {
        return "Full".to_string();
    }

    // Even subalgebra
    let even_grades: Vec<_> = (0..=dim).step_by(2).collect();
    if sorted == even_grades {
        return if dim <= 3 { "Rotor".to_string() } else { "Even".to_string() };
    }

    // Odd subalgebra
    let odd_grades: Vec<_> = (1..=dim).step_by(2).collect();
    if sorted == odd_grades {
        return "Odd".to_string();
    }

    // Specific patterns for PGA
    if sorted == vec![0, 2] && dim == 4 {
        return "Motor".to_string();
    }
    if sorted == vec![1, 3] && dim == 4 {
        return "Flector".to_string();
    }

    // Placeholder for unknown combinations
    format!("Entity_{}", sorted.iter().map(|g| g.to_string()).collect::<Vec<_>>().join("_"))
}
```

## Constraint Wrapper Suggestions

Analyze which operations preserve constraints to suggest wrappers:

### Unit Constraint

A type supports a `Unit` wrapper if:
- It has a well-defined norm (e.g., rotors, vectors)
- Some operations preserve unit norm (e.g., composition of unit rotors)

```rust
/// Checks if a type could have a meaningful unit constraint.
fn can_have_unit_constraint(grades: &[usize], algebra: &Algebra) -> bool {
    // Must be able to compute a norm (squared norm is scalar)
    // Check if self * reverse(self) gives a scalar
    satisfies_geometric_constraint(grades, algebra, &ProductTable::new(algebra))
}

/// Checks if binary operation preserves unit constraint.
fn preserves_unit(
    lhs_grades: &[usize],
    rhs_grades: &[usize],
    output_grades: &[usize],
    algebra: &Algebra,
) -> bool {
    // If both inputs are unit and output matches input type,
    // check algebraically if ||a * b|| = ||a|| * ||b|| = 1 * 1 = 1
    // This is true for rotor composition, etc.
    lhs_grades == output_grades && rhs_grades == output_grades
}
```

### NonZero Constraint

A type supports a `NonZero` wrapper if:
- It has a well-defined norm
- Normalization is a useful operation

```rust
/// All types with well-defined norms can have NonZero constraint.
fn can_have_nonzero_constraint(grades: &[usize], algebra: &Algebra) -> bool {
    can_have_unit_constraint(grades, algebra)
}
```

## TOML Template Generation

### Output Format

```toml
# Auto-discovered entities for Cl(3,0,0)
# Generated by clifford-codegen discover
#
# Review and customize:
# - Rename placeholder entities (Entity_*)
# - Add/remove types as needed
# - Customize descriptions
# - Configure constraints

[algebra]
name = "euclidean3"
module_path = "euclidean::dim3"
description = "3D Euclidean Geometric Algebra Cl(3,0,0)"

[signature]
positive = ["e1", "e2", "e3"]
negative = []
zero = []

# ============================================================
# Discovered Types (7 entities satisfy geometric constraints)
# ============================================================

[types.Scalar]
grades = [0]
description = "Grade-0 scalar"

[types.Vector]
grades = [1]
description = "Grade-1 vector"

# Suggested constraints based on norm preservation analysis:
# [types.Vector.constraints.unit]
# [types.Vector.constraints.nonzero]

[types.Bivector]
grades = [2]
description = "Grade-2 bivector (oriented plane)"

# Suggested constraints:
# [types.Bivector.constraints.unit]

[types.Trivector]
grades = [3]
description = "Grade-3 pseudoscalar"

[types.Rotor]
grades = [0, 2]
description = "Even subalgebra element (rotation)"

# Suggested constraints:
# [types.Rotor.constraints.unit]
# Unit rotor composition preserves unit constraint

[types.Odd]
grades = [1, 3]
description = "Odd subalgebra element"

[types.Full]
grades = [0, 1, 2, 3]
description = "Full multivector (all grades)"

# ============================================================
# Products (auto-populated by PRD-14.10)
# ============================================================

[products.geometric]
# Will be populated by product inference

[products.outer]
# Will be populated by product inference

[options]
generate_serde = true
generate_arbitrary = true
generate_nalgebra = true
```

### Generator Implementation

```rust
use std::io::Write;

/// Represents a discovered entity.
#[derive(Debug, Clone)]
pub struct DiscoveredEntity {
    pub name: String,
    pub grades: Vec<usize>,
    pub description: String,
    pub can_be_unit: bool,
    pub can_be_nonzero: bool,
    pub unit_preserving_ops: Vec<String>,
}

/// Discovers entities and generates a TOML template.
pub fn discover_and_generate(
    algebra: &Algebra,
    output: &mut dyn Write,
) -> std::io::Result<()> {
    let table = ProductTable::new(algebra);

    // Discover valid grade combinations
    let combinations = discover_valid_combinations(algebra, &table);

    // Create entities with names and metadata
    let entities: Vec<DiscoveredEntity> = combinations
        .into_iter()
        .map(|grades| {
            let name = suggest_name(&grades, algebra.dim());
            let can_be_unit = can_have_unit_constraint(&grades, algebra);
            let can_be_nonzero = can_have_nonzero_constraint(&grades, algebra);

            DiscoveredEntity {
                name,
                grades: grades.clone(),
                description: generate_description(&grades, algebra.dim()),
                can_be_unit,
                can_be_nonzero,
                unit_preserving_ops: vec![], // Populated by PRD-14.10
            }
        })
        .collect();

    // Generate TOML
    write_toml_template(algebra, &entities, output)
}

fn generate_description(grades: &[usize], dim: usize) -> String {
    if grades.len() == 1 {
        format!("Grade-{} element", grades[0])
    } else if grades == &(0..=dim).collect::<Vec<_>>()[..] {
        "Full multivector (all grades)".to_string()
    } else {
        format!("Grades {:?}", grades)
    }
}

fn write_toml_template(
    algebra: &Algebra,
    entities: &[DiscoveredEntity],
    output: &mut dyn Write,
) -> std::io::Result<()> {
    writeln!(output, "# Auto-discovered entities for {}", algebra.name())?;
    writeln!(output, "# Generated by clifford-codegen discover")?;
    writeln!(output)?;
    // ... rest of TOML generation
    Ok(())
}
```

## CLI Integration

Add a `discover` subcommand to the CLI:

```bash
# Discover entities for 3D Euclidean algebra
clifford-codegen discover euclidean 3

# Discover entities for 3D PGA
clifford-codegen discover pga 3

# Discover entities for 3D CGA
clifford-codegen discover cga 3

# Output to file
clifford-codegen discover euclidean 3 -o euclidean3.toml

# Verbose mode (show constraint checking)
clifford-codegen discover euclidean 3 -v
```

### CLI Implementation

```rust
#[derive(Parser)]
enum Commands {
    /// Discover valid entities from a signature
    Discover {
        /// Algebra type: euclidean, pga, cga, or custom
        algebra_type: String,

        /// Dimension (for standard algebras)
        dimension: Option<usize>,

        /// Output file (stdout if not specified)
        #[arg(short, long)]
        output: Option<PathBuf>,

        /// Custom signature: p,q,r (e.g., "3,0,1" for PGA)
        #[arg(long)]
        signature: Option<String>,

        /// Verbose output
        #[arg(short, long)]
        verbose: bool,
    },
}

fn handle_discover(cmd: DiscoverCommand) -> Result<()> {
    let algebra = match cmd.algebra_type.as_str() {
        "euclidean" => Algebra::euclidean(cmd.dimension.unwrap_or(3)),
        "pga" => Algebra::pga(cmd.dimension.unwrap_or(3)),
        "cga" => Algebra::cga(cmd.dimension.unwrap_or(3)),
        "custom" => {
            let sig = cmd.signature.ok_or("--signature required for custom")?;
            let parts: Vec<usize> = sig.split(',').map(|s| s.parse()).collect::<Result<_, _>>()?;
            Algebra::new(parts[0], parts[1], parts[2])
        }
        _ => return Err(anyhow!("Unknown algebra type")),
    };

    let mut output: Box<dyn Write> = match cmd.output {
        Some(path) => Box::new(File::create(path)?),
        None => Box::new(std::io::stdout()),
    };

    discover_and_generate(&algebra, &mut output)?;

    Ok(())
}
```

## Module Structure

```
crates/clifford-codegen/src/
├── algebra/
│   ├── constraints.rs   # From PRD-14.8
│   └── ...
├── discovery/           # NEW
│   ├── mod.rs
│   ├── entity.rs        # DiscoveredEntity type
│   ├── naming.rs        # Name suggestion heuristics
│   ├── constraints.rs   # Wrapper constraint analysis
│   └── template.rs      # TOML generation
└── cli/
    └── discover.rs      # CLI subcommand
```

## Testing

### Unit Tests

```rust
#[test]
fn discover_euclidean_3d() {
    let algebra = Algebra::euclidean(3);
    let table = ProductTable::new(&algebra);
    let entities = discover_valid_combinations(&algebra, &table);

    // Should find at least: scalar, vector, bivector, trivector, rotor, odd, full
    assert!(entities.len() >= 7);

    // Check specific combinations are found
    assert!(entities.contains(&vec![0]));       // Scalar
    assert!(entities.contains(&vec![1]));       // Vector
    assert!(entities.contains(&vec![2]));       // Bivector
    assert!(entities.contains(&vec![3]));       // Trivector
    assert!(entities.contains(&vec![0, 2]));    // Rotor
    assert!(entities.contains(&vec![1, 3]));    // Odd
    assert!(entities.contains(&vec![0, 1, 2, 3])); // Full
}

#[test]
fn naming_heuristics() {
    assert_eq!(suggest_name(&[0], 3), "Scalar");
    assert_eq!(suggest_name(&[1], 3), "Vector");
    assert_eq!(suggest_name(&[2], 3), "Bivector");
    assert_eq!(suggest_name(&[3], 3), "Trivector");
    assert_eq!(suggest_name(&[0, 2], 3), "Rotor");
    assert_eq!(suggest_name(&[1, 3], 3), "Odd");
    assert_eq!(suggest_name(&[0, 1, 2, 3], 3), "Full");

    // Unknown combination gets placeholder
    assert_eq!(suggest_name(&[0, 1], 3), "Entity_0_1");
}

#[test]
fn pga_3d_entities() {
    let algebra = Algebra::pga(3); // 4D
    let table = ProductTable::new(&algebra);
    let entities = discover_valid_combinations(&algebra, &table);

    // PGA should have motors [0, 2] and flectors [1, 3]
    assert!(entities.contains(&vec![0, 2])); // Motor
    assert!(entities.contains(&vec![1, 3])); // Flector
}
```

### Integration Tests

```rust
#[test]
fn generate_euclidean3_template() {
    let algebra = Algebra::euclidean(3);
    let mut output = Vec::new();

    discover_and_generate(&algebra, &mut output).unwrap();

    let toml_str = String::from_utf8(output).unwrap();

    // Verify TOML is valid
    let parsed: toml::Value = toml::from_str(&toml_str).unwrap();

    // Check structure
    assert!(parsed.get("algebra").is_some());
    assert!(parsed.get("signature").is_some());
    assert!(parsed.get("types").is_some());
}
```

## Deliverables

- [ ] `DiscoveredEntity` type with name, grades, constraints
- [ ] `enumerate_grade_combinations(dim) -> Vec<Vec<usize>>`
- [ ] `discover_valid_combinations(algebra, table) -> Vec<Vec<usize>>`
- [ ] `suggest_name(grades, dim) -> String` with heuristics
- [ ] `can_have_unit_constraint(grades, algebra) -> bool`
- [ ] `can_have_nonzero_constraint(grades, algebra) -> bool`
- [ ] `discover_and_generate(algebra, output)` for TOML generation
- [ ] CLI `discover` subcommand
- [ ] Tests for all standard algebras

## Success Criteria

1. **Correctness**: Discovers all valid entities for known algebras
2. **Naming**: Provides sensible default names for common entities
3. **Usability**: Generated TOML is valid and well-commented
4. **Extensibility**: Easy to add new naming patterns

## Dependencies

- PRD-14.8 (Geometric Constraint Engine)
- PRD-14.1 (Blade Algebra Engine)

## Next Steps

- PRD-14.10 uses discovered entities to infer product outputs
