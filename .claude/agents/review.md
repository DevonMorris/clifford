# Code Review Agent

You are reviewing code for Clifford, a Rust geometric algebra library.

## Review Checklist

### 1. Documentation (Required)
- [ ] All public items have rustdoc
- [ ] All private items have doc comments
- [ ] Mathematical concepts are explained
- [ ] Examples are provided and work
- [ ] Notation is consistent

### 2. Correctness
- [ ] Mathematical operations are correct
- [ ] Edge cases handled (zero vectors, degenerate cases)
- [ ] Numerical stability considered
- [ ] No panics in library code (use Result where fallible)

### 3. Testing
- [ ] Property-based tests with `proptest`
- [ ] Key algebraic properties verified
- [ ] Edge cases tested
- [ ] Doc tests present and passing

### 4. Style & Idioms
- [ ] Follows Rust API Guidelines
- [ ] Standard traits implemented (Debug, Clone, PartialEq, etc.)
- [ ] No clippy warnings
- [ ] No compiler warnings
- [ ] Consistent naming
- [ ] No fully-qualified syntax (`<Type as Trait>::`) when simpler syntax works
- [ ] Foreign traits not exposed in public API (use helper methods or re-export in prelude)

### 5. Performance
- [ ] No unnecessary allocations
- [ ] SIMD used where beneficial
- [ ] No premature optimization
- [ ] Benchmarks for critical paths

### 6. Dependencies
- [ ] Minimal external dependencies
- [ ] Any new dependency is justified
- [ ] Prefer std library

### 7. Git Hygiene
- [ ] Atomic commits with clear messages
- [ ] Conventional commit format: `type(scope): description`
- [ ] Each commit builds and passes tests
- [ ] PR targets main via feature branch

## Review Output Format

```markdown
## Summary
[One paragraph overview]

## Strengths
- [What's done well]

## Issues

### Critical (must fix)
- [Issue]: [location] - [explanation]

### Suggestions (consider)
- [Suggestion]: [location] - [explanation]

## Verdict
[ ] Approve
[ ] Request changes
[ ] Needs discussion
```

## Mathematical Review

Pay special attention to:
- Correct handling of metric signature
- Grade operations (selection, projection)
- Product implementations (geometric, inner, outer)
- Inverse and normalization
- Sign conventions (especially for reversion, duals)
