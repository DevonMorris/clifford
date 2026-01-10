"""Symbolic derivations for Clifford geometric algebra library.

This package provides SymPy-based tools for deriving and verifying
algebraic formulas used in the Clifford Rust library, particularly
for Projective Geometric Algebra (PGA).

Reference: https://rigidgeometricalgebra.org/wiki/index.php?title=Motor

Usage:
    cd derivations
    uv run python -m clifford_derivations.motor

Convention (matching clifford::specialized::projective::dim3):
    - Motor field order: s, e23, e31, e12, e01, e02, e03, e0123
    - compose(a, b) applies a first, then b (returns b * a in GP terms)
"""

from .motor import (
    derive_composition,
    derive_inverse,
    generate_rust_compose,
    generate_rust_inverse,
    to_rust,
    with_timeout,
    TimeoutError,
    gp_bivectors,
)

from .transform import (
    derive_transform_point,
    verify_special_cases,
    generate_rust_transform_point,
)

__all__ = [
    # Motor composition and inverse
    "derive_composition",
    "derive_inverse",
    "generate_rust_compose",
    "generate_rust_inverse",
    "to_rust",
    "with_timeout",
    "TimeoutError",
    "gp_bivectors",
    # Transform point
    "derive_transform_point",
    "verify_special_cases",
    "generate_rust_transform_point",
]
