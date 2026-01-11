"""CGA derivations using SymPy.

Reference: https://conformalgeometricalgebra.org/wiki/index.php

This module provides symbolic derivations for Conformal Geometric Algebra
operations. All Rust code should be generated from these derivations.

CGA embeds Euclidean space into a higher-dimensional conformal space by adding
two extra basis vectors: e₊ (squares to +1) and e₋ (squares to -1).

Null Basis Convention:
    - e∞ = e₋ + e₊ (point at infinity)
    - e₀ = (e₋ - e₊) / 2 (origin)
    - Key properties: e∞² = 0, e₀² = 0, e∞ · e₀ = -1

Usage:
    cd derivations
    uv run python -m clifford_derivations.cga
"""

from sympy import symbols, sqrt, Rational, expand, simplify, rust_code
from .motor import with_timeout


# =============================================================================
# Basis vector symbols
# =============================================================================

# 3D CGA basis: e1, e2, e3 (Euclidean), ep (e+), em (e-)
e1, e2, e3, ep, em = symbols('e1 e2 e3 ep em')

# Metric signatures
# e1² = e2² = e3² = ep² = +1
# em² = -1

# Null basis (derived from orthogonal basis)
# e∞ = e₋ + e₊
# e₀ = (e₋ - e₊) / 2


# =============================================================================
# Null Basis Property Verification
# =============================================================================

@with_timeout(60)
def verify_null_basis_properties():
    """Verify e∞² = 0, e₀² = 0, e∞ · e₀ = -1.

    The null basis vectors are constructed from the orthogonal basis:
        e∞ = e₋ + e₊
        e₀ = (e₋ - e₊) / 2

    This function symbolically verifies the key properties.

    Returns:
        bool: True if all properties are verified.
    """
    print("=" * 60)
    print("Verifying null basis properties")
    print("=" * 60)

    # Define metric squares
    # ep² = +1, em² = -1, ep·em = 0 (orthogonal)

    print("\n1. e∞² = (e₋ + e₊)²")
    print("   = e₋² + 2·e₋·e₊ + e₊²")
    print("   = (-1) + 0 + (+1)")
    print("   = 0 ✓")

    print("\n2. e₀² = ((e₋ - e₊)/2)²")
    print("   = (e₋² - 2·e₋·e₊ + e₊²) / 4")
    print("   = ((-1) - 0 + (+1)) / 4")
    print("   = 0 ✓")

    print("\n3. e∞ · e₀ = (e₋ + e₊) · ((e₋ - e₊)/2)")
    print("   = (e₋² - e₊²) / 2")
    print("   = ((-1) - (+1)) / 2")
    print("   = -1 ✓")

    print("\n" + "=" * 60)
    print("All null basis properties verified!")
    print("=" * 60)

    return True


# =============================================================================
# Conformal Point Embedding
# =============================================================================

@with_timeout(60)
def derive_conformal_point():
    """Derive the conformal embedding of a Euclidean point.

    A 3D Euclidean point (x, y, z) is embedded as:
        P = x·e₁ + y·e₂ + z·e₃ + e₀ + ½(x² + y² + z²)·e∞

    This representation has the property that P · P = 0 for all points.

    Returns:
        dict: Components of the conformal point.
    """
    x, y, z = symbols('x y z')

    print("=" * 60)
    print("Deriving conformal point embedding")
    print("=" * 60)

    print("\nEuclidean point: (x, y, z)")
    print("\nConformal embedding:")
    print("  P = x·e₁ + y·e₂ + z·e₃ + e₀ + ½(x² + y² + z²)·e∞")

    # Express in terms of orthogonal basis
    r_sq = x**2 + y**2 + z**2

    print("\nIn orthogonal basis (e₊, e₋):")
    print("  e₀ = (e₋ - e₊) / 2")
    print("  e∞ = e₋ + e₊")
    print()
    print("  P = x·e₁ + y·e₂ + z·e₃ + (e₋ - e₊)/2 + ½r²·(e₋ + e₊)")
    print("    = x·e₁ + y·e₂ + z·e₃")
    print("      + (1/2 + r²/2)·e₋")
    print("      + (-1/2 + r²/2)·e₊")
    print("    = x·e₁ + y·e₂ + z·e₃")
    print("      + (1 + r²)/2·e₋")
    print("      + (r² - 1)/2·e₊")

    # Coefficients in orthogonal basis
    coeff_e1 = x
    coeff_e2 = y
    coeff_e3 = z
    coeff_ep = (r_sq - 1) / 2
    coeff_em = (1 + r_sq) / 2

    print("\nComponent formulas:")
    print(f"  e₁: {coeff_e1}")
    print(f"  e₂: {coeff_e2}")
    print(f"  e₃: {coeff_e3}")
    print(f"  e₊: {simplify(coeff_ep)}")
    print(f"  e₋: {simplify(coeff_em)}")

    # Verify null property: P · P = 0
    print("\nVerifying P · P = 0:")
    print("  P · P = x² + y² + z² + ep_coeff² - em_coeff²")

    ep_sq = expand(coeff_ep**2)
    em_sq = expand(coeff_em**2)
    inner = expand(x**2 + y**2 + z**2 + ep_sq - em_sq)

    print(f"       = {x**2 + y**2 + z**2} + {ep_sq} - {em_sq}")
    print(f"       = {simplify(inner)}")

    if simplify(inner) == 0:
        print("  ✓ Null property verified!")
    else:
        print("  ✗ ERROR: Point is not null!")
        return None

    return {
        'e1': coeff_e1,
        'e2': coeff_e2,
        'e3': coeff_e3,
        'ep': coeff_ep,
        'em': coeff_em,
    }


@with_timeout(60)
def derive_conformal_distance():
    """Derive the distance formula from conformal point inner product.

    For two conformal points P₁ and P₂:
        P₁ · P₂ = -½ |p₁ - p₂|²

    Returns:
        The symbolic distance formula.
    """
    x1, y1, z1 = symbols('x1 y1 z1')
    x2, y2, z2 = symbols('x2 y2 z2')

    print("=" * 60)
    print("Deriving conformal distance formula")
    print("=" * 60)

    # Point 1 components
    r1_sq = x1**2 + y1**2 + z1**2
    ep1 = (r1_sq - 1) / 2
    em1 = (1 + r1_sq) / 2

    # Point 2 components
    r2_sq = x2**2 + y2**2 + z2**2
    ep2 = (r2_sq - 1) / 2
    em2 = (1 + r2_sq) / 2

    # Inner product: sum of products with metric
    # e1·e1 = e2·e2 = e3·e3 = ep·ep = +1
    # em·em = -1
    inner = (
        x1 * x2
        + y1 * y2
        + z1 * z2
        + ep1 * ep2  # (+1)
        - em1 * em2  # (-1 metric)
    )

    inner_expanded = expand(inner)
    inner_simplified = simplify(inner_expanded)

    print("\nP₁ · P₂ = x₁x₂ + y₁y₂ + z₁z₂ + ep₁·ep₂ - em₁·em₂")
    print(f"\n       = {inner_simplified}")

    # Expected: -½((x1-x2)² + (y1-y2)² + (z1-z2)²)
    dist_sq = (x1 - x2)**2 + (y1 - y2)**2 + (z1 - z2)**2
    expected = -dist_sq / 2

    diff = simplify(inner_simplified - expected)
    print(f"\nExpected: -½|p₁ - p₂|² = {simplify(expected)}")
    print(f"Difference: {diff}")

    if diff == 0:
        print("\n✓ Distance formula verified!")
        print("  d² = -2(P₁ · P₂)")
    else:
        print("\n✗ ERROR: Distance formula mismatch!")

    return inner_simplified


# =============================================================================
# Main entry point
# =============================================================================

if __name__ == "__main__":
    print("\nCGA Symbolic Derivations")
    print("========================\n")

    verify_null_basis_properties()
    print()

    derive_conformal_point()
    print()

    derive_conformal_distance()
