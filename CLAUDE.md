# Clifford - Geometric Algebra Library

A Rust library for Geometric Algebra (Clifford Algebra).

## Project Principles

### 1. Educational Focus
- Code should be readable and well-documented to facilitate learning GA
- Include mathematical explanations in doc comments
- Link to resources and papers where appropriate

### 2. Documentation First
- Every public API must have comprehensive rustdoc
- Include examples in doc comments
- Explain the geometric intuition, not just the code

### 3. Clean Git History
- Atomic commits with clear messages
- Each commit should be buildable and pass tests
- Use conventional commit format: `type(scope): description`

### 4. Performance
- Use SIMD instructions where beneficial (via `std::arch` or `portable_simd`)
- Benchmark critical paths
- Profile before optimizing

### 5. Minimal Dependencies
- Prefer std library where possible
- Only add dependencies when they provide significant value
- Document why each dependency is needed

### 6. Rust API Guidelines
- Follow the official Rust API Guidelines (https://rust-lang.github.io/api-guidelines/)
- Ergonomic builder patterns where appropriate
- Implement standard traits (Debug, Clone, PartialEq, etc.)

### 7. Testing
- Unit tests for all functionality
- Property-based tests for algebraic identities
- Doc tests for examples

## Development Commands

```bash
cargo build           # Build the library
cargo test            # Run all tests
cargo doc --open      # Generate and view documentation
cargo bench           # Run benchmarks
cargo clippy          # Run linter
cargo fmt             # Format code
```

## Architecture Notes

(To be expanded as the library develops)
