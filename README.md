<p align="center">
  <img src="assets/clifford.png" alt="Clifford" width="240">
</p>

<h1 align="center">Clifford</h1>

<p align="center">
  <strong>Geometric Algebra for Rust</strong>
</p>

<p align="center">
  <a href="https://crates.io/crates/clifford"><img src="https://img.shields.io/crates/v/clifford.svg" alt="crates.io"></a>
  <a href="https://docs.rs/clifford"><img src="https://docs.rs/clifford/badge.svg" alt="docs.rs"></a>
  <a href="https://github.com/DevonMorris/clifford/actions"><img src="https://github.com/DevonMorris/clifford/workflows/CI/badge.svg" alt="CI"></a>
</p>

---

A Rust library for Geometric Algebra (Clifford Algebra) with a focus on correctness, performance, and education.

## Installation

```toml
[dependencies]
clifford = "0.1"
```

## Overview

Geometric Algebra unifies and extends linear algebra, complex numbers, quaternions, and more into a single coherent framework. This library provides:

- Generic multivector types over any metric signature
- Geometric, inner, and outer products
- Support for Euclidean, projective (PGA), and conformal (CGA) algebras
- Compile-time dimension checking

## License

MIT
