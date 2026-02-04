# Clifford Codegen TOML Specification v1.0

This document provides a complete reference for the TOML specification format used by `clifford-codegen` to define geometric algebras.

## Table of Contents

1. [Overview](#overview)
2. [Required Sections](#required-sections)
   - [[algebra]](#algebra)
   - [[signature]](#signature)
   - [[types]](#types)
3. [Optional Sections](#optional-sections)
   - [[norm]](#norm)
   - [[blades]](#blades)
4. [Field Map Syntax](#field-map-syntax)
5. [Automatic Detection](#automatic-detection)
6. [Examples](#examples)
7. [Validation Rules](#validation-rules)

---

## Overview

Each algebra is defined in a single TOML file in the `algebras/` directory. The codegen tool parses these specifications and generates optimized Rust code for the algebra's types and operations.

### Minimal Example

```toml
[algebra]
name = "euclidean2"
module_path = "euclidean::dim2"

[signature]
positive = ["e1", "e2"]

[types.Vector]
grades = [1]
field_map = [
  { name = "x", blade = "e1" },
  { name = "y", blade = "e2" }
]
```

---

## Required Sections

### [algebra]

Core metadata about the algebra.

| Field | Type | Required | Default | Description |
|-------|------|----------|---------|-------------|
| `name` | string | Yes | - | Unique identifier (e.g., "euclidean3") |
| `module_path` | string | No | derived from name | Rust module path (e.g., "euclidean::dim3") |
| `description` | string | No | - | Documentation string for the algebra |
| `complete` | bool | No | true | If true, validates all products have matching output types |

**Example:**
```toml
[algebra]
name = "projective3"
module_path = "projective::dim3"
description = "3D Projective Geometric Algebra"
complete = false  # Some products don't have explicit output types
```

### [signature]

Defines the metric signature of the algebra. Determines how basis vectors square.

| Field | Type | Required | Default | Description |
|-------|------|----------|---------|-------------|
| `positive` | string[] | Conditional | [] | Basis vectors squaring to +1 |
| `negative` | string[] | No | [] | Basis vectors squaring to -1 |
| `zero` | string[] | No | [] | Degenerate basis vectors squaring to 0 |

**Constraint:** At least one of `positive`, `negative`, or `zero` must be non-empty.

**Indexing Rules:**
- Numeric names (e1, e2, e3): Index = number - 1 (e1 -> 0, e2 -> 1)
- Special case: e0 -> dim - 1 (PGA convention for degenerate basis)
- Non-numeric names (ep, em): Positional indexing (positive first, then negative, then zero)
- Total dimension (p + q + r) must be 1-6

**Examples:**
```toml
# Euclidean 3D: Cl(3,0,0)
[signature]
positive = ["e1", "e2", "e3"]
negative = []
zero = []

# Minkowski plane: Cl(1,1,0)
[signature]
positive = ["e2"]  # time
negative = ["e1"]  # space
zero = []

# 3D PGA: Cl(3,0,1)
[signature]
positive = ["e1", "e2", "e3"]
negative = []
zero = ["e0"]

# 2D CGA: Cl(3,1,0)
[signature]
positive = ["e1", "e2", "ep"]
negative = ["em"]
zero = []
```

### [types]

Defines algebraic types (elements of the algebra). Each type is defined as `[types.TypeName]`.

| Field | Type | Required | Default | Description |
|-------|------|----------|---------|-------------|
| `grades` | int[] | Yes | - | Array of grades this type contains |
| `description` | string | No | - | Documentation for this type |
| `field_map` | array | Conditional | - | Required unless `alias_of` is set |
| `alias_of` | string | No | - | Makes this type an alias of another |
| `inverse_sandwich_targets` | string[] | No | [] | For non-versor types supporting inverse sandwich |

**Example:**
```toml
[types.Scalar]
grades = [0]
description = "Scalar (grade-0 element)"
field_map = [
  { name = "s", blade = "s" }
]

[types.Vector]
grades = [1]
description = "3D vector"
field_map = [
  { name = "x", blade = "e1" },
  { name = "y", blade = "e2" },
  { name = "z", blade = "e3" }
]

[types.Rotor]
grades = [0, 2]
description = "3D rotor (unit versor for rotations)"
field_map = [
  { name = "s", blade = "s" },
  { name = "rz", blade = "e12" },
  { name = "ry", blade = "e13" },
  { name = "rx", blade = "e23" }
]
```

---

## Optional Sections

### [norm]

Controls which involution is used for norm computation.

| Field | Type | Required | Default | Description |
|-------|------|----------|---------|-------------|
| `primary_involution` | string | No | "reverse" | Involution for norm: "reverse", "grade_involution", or "clifford_conjugate" |

**When to use each involution:**
- `"reverse"`: Standard for Euclidean and projective algebras. `|x|^2 = x * reverse(x)`
- `"grade_involution"`: Rarely used directly for norms
- `"clifford_conjugate"`: Used for complex numbers and some other algebras

**Example:**
```toml
[norm]
# Complex numbers use Clifford conjugate: z * conj(z) = a^2 + b^2
primary_involution = "clifford_conjugate"
```

### [blades]

Custom display names for blades. Maps internal blade names to user-friendly names.

**Example:**
```toml
[blades]
# PGA blade naming
e4 = "e0"      # e4 is the degenerate basis
e14 = "e01"    # Bivector naming
e123 = "e123"  # Trivector
e1234 = "e0123"  # Quadvector (pseudoscalar)
```

This is primarily used for documentation and display purposes. The internal blade ordering uses numeric indices based on the signature.

---

## Field Map Syntax

Each field in a type's `field_map` has the following structure:

```toml
{ name = "field_name", blade = "blade_reference" }
```

### Blade References

| Reference | Meaning |
|-----------|---------|
| `"s"` | Scalar (grade 0) |
| `"e1"`, `"e2"`, etc. | Basis vectors |
| `"e12"`, `"e123"`, etc. | Higher-grade blades |
| `"e21"`, `"e31"`, etc. | Non-canonical ordering (auto-computes sign) |

**Non-canonical blade ordering:**
When you specify a blade like `"e21"` instead of `"e12"`, the codegen automatically computes the sign flip. This is useful when the natural field semantics require a negated blade.

**Example:**
```toml
# Line in 2D PGA
[types.Line]
grades = [2]
field_map = [
  { name = "d", blade = "e12" },
  { name = "nx", blade = "e23" },
  { name = "ny", blade = "e31" }  # e31 = -e13, auto-computed
]
```

### Semantic Field Naming

**IMPORTANT:** Field names should reflect what the field DOES, not which blade it corresponds to.

**Good naming:**
```toml
# Motor fields describe their geometric effect
field_map = [
  { name = "s", blade = "s" },
  { name = "tz", blade = "e12" },  # Controls z-translation
  { name = "ty", blade = "e13" },  # Controls y-translation
  { name = "tx", blade = "e14" },  # Controls x-translation
  { name = "rx", blade = "e23" },  # Controls x-rotation
  { name = "ry", blade = "e24" },  # Controls y-rotation
  { name = "rz", blade = "e34" },  # Controls z-rotation
  { name = "ps", blade = "e1234" } # Pseudoscalar
]
```

**Bad naming:**
```toml
# Don't name fields after blade indices
field_map = [
  { name = "e12", blade = "e12" },  # Unclear what it does
  { name = "e13", blade = "e13" },
]
```

---

## Automatic Detection

The codegen tool automatically detects several properties:

### Versor Detection

A type is automatically detected as a versor if all its grades have the same parity (all even OR all odd). Versors can transform other elements via sandwich products.

```toml
# Rotor has grades [0, 2] - all even -> versor
[types.Rotor]
grades = [0, 2]

# Motor has grades [0, 2, 4] - all even -> versor
[types.Motor]
grades = [0, 2, 4]

# Flector has grades [1, 3] - all odd -> versor
[types.Flector]
grades = [1, 3]

# Vector has grades [1] - odd -> versor (single-grade is trivially same parity)
[types.Vector]
grades = [1]
```

### Sparse Type Detection

A type is sparse if its field count is less than the expected blade count for its grades. Sparse types have constraints that reduce the degrees of freedom.

```toml
# CGA Line has grade [3] but only 6 of 10 grade-3 blades
[types.Line]
grades = [3]
field_map = [
  # Only 6 fields instead of 10
  { name = "vx", blade = "e124" },
  { name = "vy", blade = "e134" },
  { name = "vz", blade = "e234" },
  { name = "mx", blade = "e125" },
  { name = "my", blade = "e135" },
  { name = "mz", blade = "e235" }
]
```

---

## Examples

### Euclidean 3D (Complete)

```toml
[algebra]
name = "euclidean3"
module_path = "euclidean::dim3"
description = "3D Euclidean Geometric Algebra Cl(3,0,0)"
complete = false

[signature]
positive = ["e1", "e2", "e3"]
negative = []
zero = []

[norm]
primary_involution = "reverse"

[types.Scalar]
grades = [0]
description = "Scalar (grade-0 element)"
field_map = [
  { name = "s", blade = "s" }
]

[types.Vector]
grades = [1]
description = "3D vector"
field_map = [
  { name = "x", blade = "e1" },
  { name = "y", blade = "e2" },
  { name = "z", blade = "e3" }
]

[types.Bivector]
grades = [2]
description = "3D bivector"
field_map = [
  { name = "rz", blade = "e12" },
  { name = "ry", blade = "e13" },
  { name = "rx", blade = "e23" }
]

[types.Trivector]
grades = [3]
description = "3D trivector (pseudoscalar)"
field_map = [
  { name = "ps", blade = "e123" }
]

[types.Rotor]
grades = [0, 2]
description = "3D rotor (unit versor for rotations)"
field_map = [
  { name = "s", blade = "s" },
  { name = "rz", blade = "e12" },
  { name = "ry", blade = "e13" },
  { name = "rx", blade = "e23" }
]
```

### Complex Numbers

```toml
[algebra]
name = "complex"
module_path = "complex"
description = "Complex numbers Cl(0,1,0)"

[signature]
positive = []
negative = ["e1"]
zero = []

[norm]
primary_involution = "clifford_conjugate"

[types.Scalar]
grades = [0]
field_map = [
  { name = "s", blade = "s" }
]

[types.ImagUnit]
grades = [1]
description = "Pure imaginary unit (i where i^2 = -1)"
field_map = [
  { name = "i", blade = "e1" }
]

[types.Complex]
grades = [0, 1]
description = "Complex number: a + bi where i^2 = -1"
field_map = [
  { name = "real", blade = "s" },
  { name = "imag", blade = "e1" }
]
```

---

## Validation Rules

The codegen tool validates specifications against these rules:

1. **Algebra section required:** Must have `[algebra]` with at least `name`
2. **Signature required:** Must have `[signature]` with at least one non-empty basis array
3. **Dimension limit:** Total dimension (p + q + r) must be 1-6
4. **Type requirements:** Each type must have `grades` and either `field_map` or `alias_of`
5. **Field uniqueness:** Field names must be unique within a type
6. **Blade validity:** All blade references must be valid for the algebra's dimension
7. **Grade consistency:** Fields must match the declared grades
8. **Completeness check:** If `complete = true`, all products must have matching output types

Run validation with:
```bash
cargo run -p clifford-codegen -- validate path/to/algebra.toml
```
