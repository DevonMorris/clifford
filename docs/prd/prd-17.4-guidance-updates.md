# PRD-17.4: Guidance Updates for Generated Products

**Status**: Complete
**Parent**: [PRD-17: Codegen Product Completeness and Documentation](prd-17-codegen-products.md)
**Goal**: Update CLAUDE.md and agent files to direct users to use generated products

## Problem Statement

Before this PRD, the documentation and agent guidance did not explicitly instruct implementers to use generated products from `generated/products.rs` in extension files. This led to:

1. Manual algebraic formulas in extension files
2. Error-prone implementations (e.g., the `Plane::meet` parameter ordering bug)
3. Inconsistent code quality between generated and hand-written code

## Solution

Update all guidance documentation to explicitly direct users to:

1. Check for existing generated products before implementing
2. Use generated products in extension files
3. File issues for missing products rather than hand-rolling formulas

## Changes Made

### CLAUDE.md

Added "Using Generated Products in Extensions" section after "When to Use Codegen":

```markdown
#### Using Generated Products in Extensions

When implementing domain-specific methods in extension files (`extensions.rs`),
**always use generated products** instead of manual formulas:

```rust
// CORRECT: Use generated product from generated/products.rs
use super::generated::products;

impl<T: Float> Point<T> {
    pub fn join(&self, other: &Point<T>) -> Line<T> {
        // Use generated exterior product
        products::exterior_point_point(self, other)
    }
}

// WRONG: Manual formula (error-prone, not verified)
impl<T: Float> Point<T> {
    pub fn join(&self, other: &Point<T>) -> Line<T> {
        // Don't do this - prone to sign/ordering errors
        Line::new_unchecked(
            self.e1() * other.e2() - self.e2() * other.e1(),  // might be wrong!
            // ...
        )
    }
}
```

**Acceptable exceptions** (when manual formulas are okay):
1. The product combination isn't generated yet (file an issue to add it!)
2. Performance-critical code where the generated form is suboptimal
3. Geometric shortcuts (e.g., `distance()` using a specialized formula)

Even for exceptions, add a comment citing the mathematical source of the formula.

**Red flags in code review**:
- Manual product formulas in extension files
- Sign corrections or coefficient adjustments
- Comments like "// manual implementation because codegen doesn't..."

These indicate missing codegen features that should be added via PRD-17.
```

### .claude/agents/implement.md

Added "Using Generated Products in Extensions" section:

```markdown
### Using Generated Products in Extensions

When implementing domain-specific methods in extension files:

1. **Check if a generated product exists** in `generated/products.rs`
2. **Use the generated product** rather than deriving formulas manually
3. **If no suitable product exists**, file an issue for codegen enhancement (PRD-17)

```rust
// CORRECT: Use generated products
use super::generated::products;

impl<T: Float> Point<T> {
    pub fn join(&self, other: &Point<T>) -> Line<T> {
        products::exterior_point_point(self, other)
    }

    pub fn left_contract_plane(&self, plane: &Plane<T>) -> T {
        products::left_contract_point_plane(self, plane)
    }
}

// WRONG: Manual formulas
impl<T: Float> Point<T> {
    pub fn join(&self, other: &Point<T>) -> Line<T> {
        // Don't hand-roll - use products::exterior_point_point
        Line::new_unchecked(
            self.e1() * other.e2() - self.e2() * other.e1(),
            // ...
        )
    }
}
```

**Red flags during implementation**:
- Writing multi-term algebraic expressions manually
- Sign corrections or "magic" coefficients
- Copy-pasting formulas from papers without codegen

These indicate missing codegen features - file an issue rather than working around.
```

Also updated the "What Codegen Handles" section to list all product types:

```markdown
### What Codegen Handles

- **Products**: Geometric, exterior, interior, left/right contraction - all generated correctly
```

### .claude/agents/review.md

Added "Generated Products in Extensions" checklist:

```markdown
### Generated Products in Extensions

When reviewing extension files (`extensions.rs`):

- [ ] **Uses generated products** - Methods use `products::*` functions from `generated/products.rs`
- [ ] **No manual formulas** - Algebraic products aren't hand-rolled with multi-term expressions
- [ ] **Missing products documented** - If manual formula needed, comment explains why and cites source
- [ ] **Issue filed for gaps** - Missing codegen products have issues filed (PRD-17)

**Acceptable exceptions** for manual formulas:
1. Product combination not yet generated (must have issue filed)
2. Performance-critical path with documented benchmark justification
3. Geometric shortcut formula with cited mathematical source

**Examples of what to flag**:
```rust
// FLAG: Manual exterior product - should use products::exterior_point_point
let line = Line::new_unchecked(
    self.e1() * other.e2() - self.e2() * other.e1(),
    // ...
);

// FLAG: Manual left contraction - should use products::left_contract_*
let result = self.e1() * plane.e023() + self.e2() * plane.e031() + ...;

// OK: Uses generated product
let line = products::exterior_point_point(self, other);
```
```

## Verification

The guidance is now in place. To verify effectiveness:

1. **New implementations** should follow the pattern
2. **Code reviews** should flag manual formulas
3. **Missing products** should result in PRD-17 issues, not workarounds

## Deliverables

- [x] Update CLAUDE.md with "Using Generated Products in Extensions" section
- [x] Update CLAUDE.md with "Regenerating Algebras After Codegen Changes" section
- [x] Update implement.md agent with generated product guidance
- [x] Update implement.md agent with algebra regeneration guidance
- [x] Update implement.md "What Codegen Handles" section
- [x] Update review.md agent with generated products checklist
- [x] Update review.md agent with codegen change requirements checklist

## Success Criteria

1. **Guidance visible**: New implementers find the guidance
2. **Reviews catch issues**: Manual formulas are flagged in review
3. **Issues filed**: Missing products result in codegen enhancement requests

## Files Changed

| File | Status | Change |
|------|--------|--------|
| `CLAUDE.md` | Complete | Added "Using Generated Products in Extensions" |
| `.claude/agents/implement.md` | Complete | Added product usage guidance |
| `.claude/agents/review.md` | Complete | Added products checklist |

## Notes

This PRD was completed as part of the initial PRD-17 creation. The guidance is now in place and will be effective for future implementations.

When PRD-17.1 (missing products) and PRD-17.2 (exterior naming) are complete, the examples in the guidance will need to be updated to reflect the actual function names (e.g., `exterior_*` instead of `outer_*`).
