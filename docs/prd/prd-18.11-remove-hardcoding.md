# PRD-18.11: Remove Algebra-Specific Hardcoding

**Status**: Draft
**Parent**: [PRD-18: Constraint System Redesign](prd-18-constraint-redesign.md)

## Goal

Eliminate all PGA/Euclidean/CGA-specific string processing and derive everything from signature analysis. The code generator should be completely generic - no shortcuts based on algebra names.

## Problem

The codebase has hardcoded assumptions about specific algebras:

### Critical Issues Found

#### 1. Duplicate `generate_signature_name()` Functions

**Files**: `traits.rs:815-841`, `conversions.rs:113-137`

```rust
// PROBLEM: Checks string prefixes to determine algebra type
fn generate_signature_name(&self) -> proc_macro2::Ident {
    let name = &self.spec.name;
    let sig_name = if name.starts_with("euclidean") {
        match self.algebra.dim() {
            2 => "Euclidean2",
            3 => "Euclidean3",
            _ => "Euclidean3",
        }
    } else if name.starts_with("pga") || name.starts_with("projective") {
        match self.algebra.dim() {
            3 => "Projective2",
            4 => "Projective3",
            _ => "Projective3",
        }
    } else if name.starts_with("cga") || name.starts_with("conformal") {
        // ...
    } else {
        "CustomSignature"
    };
    format_ident!("{}", sig_name)
}
```

#### 2. Hardcoded Constraint Formulas

**File**: `extensions.rs`

- Study condition (lines 1013-1026) - hardcoded formula
- Plücker condition (lines 365-375) - hardcoded formula

These should be generated from symbolic analysis.

#### 3. Hardcoded Bulk/Weight Norm Formulas

**File**: `generated/traits.rs:1110-1251`

Each type has manually-specified bulk/weight norm formulas. These should be derived from which basis elements contain degenerate indices.

#### 4. Manual Product Formulas

**File**: `extensions.rs`

- Motor composition (lines 903-958) - 64-term explicit formula
- Regressive products (lines 421-442, 694-718)
- Left contraction (lines 580-594)

These should use generated products from `products.rs`.

#### 5. PGA vs Euclidean Branching

**File**: `types.rs:625-633`

```rust
fn generate_wrapper_aliases(&self) -> TokenStream {
    let is_pga = self.spec.signature.r > 0;
    if is_pga {
        self.generate_pga_wrapper_aliases()
    } else {
        self.generate_euclidean_wrapper_aliases()
    }
}
```

While this uses signature, it still has separate code paths.

## Solution: Derive Everything from Signature

The signature `Cl(p, q, r)` provides:
- `p` - number of basis vectors that square to +1
- `q` - number of basis vectors that square to -1
- `r` - number of basis vectors that square to 0 (degenerate)

### Derivation Rules

| Property | Derivation |
|----------|------------|
| Norm trait | `r > 0` → DegenerateNormed, `q > 0` → IndefiniteNormed, else → Normed |
| Bulk components | Fields whose blade contains NO degenerate basis indices |
| Weight components | Fields whose blade contains ANY degenerate basis index |
| Constraint expressions | Compute `u * rev(u)` symbolically, extract non-scalar terms |
| Wrapper aliases | Based on which norm traits the type implements |
| Signature type name | From `(p, q, r)` tuple: `Cl{p}_{q}_{r}` or similar |

### Bulk/Weight Decomposition Algorithm

```rust
fn derive_bulk_weight_split(ty: &TypeSpec, signature: &Signature) -> (Vec<Field>, Vec<Field>) {
    let degenerate_indices: HashSet<usize> = signature.degenerate_indices().collect();

    let (bulk, weight): (Vec<_>, Vec<_>) = ty.fields.iter().partition(|field| {
        let blade = Blade::from_index(field.blade_index);
        // Bulk: blade contains NO degenerate basis vectors
        !blade.indices().any(|i| degenerate_indices.contains(&i))
    });

    (bulk, weight)
}
```

## Implementation

### New Files

- `crates/clifford-codegen/src/codegen/norms.rs` - Generate bulk/weight formulas from degenerate indices
- `crates/clifford-codegen/src/symbolic/norm_derive.rs` - Derive norm formulas from metric

### Modified Files

```
crates/clifford-codegen/src/
├── codegen/
│   ├── traits.rs      - Remove generate_signature_name(), derive from signature
│   ├── conversions.rs - Remove duplicate generate_signature_name()
│   └── types.rs       - Unify wrapper alias generation
├── algebra/
│   └── signature.rs   - Add degenerate_indices(), positive_indices(), negative_indices()

src/specialized/
├── projective/dim3/
│   └── extensions.rs  - Replace manual formulas with generated products
└── euclidean/dim3/
    └── generated/products.rs - Remove "In PGA" comments
```

### Tasks

1. **Add Signature query methods**
   ```rust
   impl Signature {
       fn degenerate_indices(&self) -> impl Iterator<Item = usize>;
       fn positive_indices(&self) -> impl Iterator<Item = usize>;
       fn negative_indices(&self) -> impl Iterator<Item = usize>;
   }
   ```

2. **Create `derive_bulk_weight_split()`**
   - Partition type fields by degenerate basis involvement
   - Generate `bulk_norm_squared()` as sum of squares of bulk components
   - Generate `weight_norm_squared()` as sum of squares of weight components

3. **Remove all string-based algebra detection**
   - Delete `name.starts_with("euclidean")` checks
   - Delete `name.starts_with("pga")` checks
   - Delete `name.starts_with("cga")` checks

4. **Generate signature type names from tuple**
   - `Cl(3,0,0)` → `Signature3_0_0` or `Euclidean3`
   - `Cl(3,0,1)` → `Signature3_0_1` or `Degenerate3_0_1`
   - Decision: Use signature tuple directly, no algebra names

5. **Unify wrapper alias generation**
   - If type implements `DegenerateNormed` → `Bulk<T>`, `Ideal<T>` aliases
   - If type implements `IndefiniteNormed` → `Proper<T>` aliases
   - If type implements only `Normed` → `Unit<T>` aliases

6. **Replace manual formulas in extensions.rs**
   - Use `products::geometric_motor_motor()` for composition
   - Use `products::regressive_*()` for meet operations
   - Use `products::left_contract_*()` for contractions

7. **Remove PGA-specific comments from generated code**
   - Delete `// In PGA, use this instead` comments
   - Documentation should be signature-based, not algebra-based

## Verification

After implementation:
```bash
# Verify no algebra-specific strings remain
grep -r "starts_with.*euclidean" crates/clifford-codegen/src/
grep -r "starts_with.*pga" crates/clifford-codegen/src/
grep -r "starts_with.*projective" crates/clifford-codegen/src/

# All should return empty
```

```bash
# Regenerate and verify
for toml in algebras/*.toml; do
    cargo run --package clifford-codegen -- generate "$toml" --force
done
cargo fmt && cargo clippy && cargo test
```

## Success Criteria

1. **No string-based algebra detection** - Zero instances of `name.starts_with("euclidean")` or similar
2. **Bulk/weight derived from signature** - Formulas generated from degenerate index analysis
3. **Signature names from tuple** - No hardcoded "Euclidean3", "Projective3" strings
4. **No manual product formulas** - extensions.rs uses generated products
5. **Generalizes to any algebra** - Same code handles `Cl(p,q,r)` for any valid signature
