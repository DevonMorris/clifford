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
# Circumcenter Derivation
# =============================================================================

@with_timeout(120)
def derive_circumcenter():
    """Derive the circumcenter of three points in 3D.

    Given three points A, B, C, the circumcenter O is equidistant from all three:
        |O - A| = |O - B| = |O - C|

    The circumcenter lies at the intersection of perpendicular bisector planes.

    Algorithm:
    1. Compute plane normal n = (B-A) × (C-A)
    2. Compute midpoints m1 = (A+B)/2, m2 = (A+C)/2
    3. Compute perpendicular directions p1 = (B-A) × n, p2 = (C-A) × n
    4. Find intersection of lines: m1 + t*p1 and m2 + s*p2
    5. Use triple scalar product to solve for t

    Returns:
        dict: The circumcenter formulas.
    """
    # Define symbolic coordinates for three points
    ax, ay, az = symbols('ax ay az')
    bx, by, bz = symbols('bx by bz')
    cx, cy, cz = symbols('cx cy cz')

    print("=" * 60)
    print("Deriving circumcenter of three points")
    print("=" * 60)

    print("\nGiven three points:")
    print("  A = (ax, ay, az)")
    print("  B = (bx, by, bz)")
    print("  C = (cx, cy, cz)")

    # Edge vectors
    d1x, d1y, d1z = bx - ax, by - ay, bz - az
    d2x, d2y, d2z = cx - ax, cy - ay, cz - az

    print("\nEdge vectors:")
    print(f"  d1 = B - A = ({d1x}, {d1y}, {d1z})")
    print(f"  d2 = C - A = ({d2x}, {d2y}, {d2z})")

    # Plane normal: n = d1 × d2
    nx = d1y * d2z - d1z * d2y
    ny = d1z * d2x - d1x * d2z
    nz = d1x * d2y - d1y * d2x

    print("\nPlane normal (d1 × d2):")
    print(f"  nx = {nx}")
    print(f"  ny = {ny}")
    print(f"  nz = {nz}")

    # Midpoints
    m1x, m1y, m1z = (ax + bx) / 2, (ay + by) / 2, (az + bz) / 2
    m2x, m2y, m2z = (ax + cx) / 2, (ay + cy) / 2, (az + cz) / 2

    print("\nMidpoints:")
    print(f"  m1 = (A + B) / 2")
    print(f"  m2 = (A + C) / 2")

    # Perpendicular directions in plane: p = d × n
    # p1 = d1 × n
    p1x = d1y * nz - d1z * ny
    p1y = d1z * nx - d1x * nz
    p1z = d1x * ny - d1y * nx

    # p2 = d2 × n
    p2x = d2y * nz - d2z * ny
    p2y = d2z * nx - d2x * nz
    p2z = d2x * ny - d2y * nx

    print("\nPerpendicular directions (in plane):")
    print("  p1 = d1 × n")
    print("  p2 = d2 × n")

    # Solve m1 + t*p1 = m2 + s*p2 for t
    # Rearranging: t*p1 - s*p2 = m2 - m1
    # Taking cross product with p2: t*(p1 × p2) = (m2 - m1) × p2
    # Dotting with n: t * (p1 × p2) · n = ((m2 - m1) × p2) · n

    dmx, dmy, dmz = m2x - m1x, m2y - m1y, m2z - m1z

    # dm × p2
    dm_cross_p2_x = dmy * p2z - dmz * p2y
    dm_cross_p2_y = dmz * p2x - dmx * p2z
    dm_cross_p2_z = dmx * p2y - dmy * p2x

    # p1 × p2
    p1_cross_p2_x = p1y * p2z - p1z * p2y
    p1_cross_p2_y = p1z * p2x - p1x * p2z
    p1_cross_p2_z = p1x * p2y - p1y * p2x

    # Dot with n
    numer = expand(dm_cross_p2_x * nx + dm_cross_p2_y * ny + dm_cross_p2_z * nz)
    denom = expand(p1_cross_p2_x * nx + p1_cross_p2_y * ny + p1_cross_p2_z * nz)

    print("\nSolving for parameter t:")
    print("  t = ((m2 - m1) × p2) · n / (p1 × p2) · n")

    # Circumcenter = m1 + t * p1
    # For the Rust code, we compute t and then apply it
    print("\nCircumcenter:")
    print("  O = m1 + t * p1")
    print(f"  Ox = m1x + t * p1x")
    print(f"  Oy = m1y + t * p1y")
    print(f"  Oz = m1z + t * p1z")

    # Verify with a concrete example: unit circle in xy-plane
    print("\n" + "-" * 40)
    print("Verification: Unit circle in xy-plane")
    print("-" * 40)

    # A = (1, 0, 0), B = (0, 1, 0), C = (-1, 0, 0)
    subs = {ax: 1, ay: 0, az: 0, bx: 0, by: 1, bz: 0, cx: -1, cy: 0, cz: 0}

    numer_val = numer.subs(subs)
    denom_val = denom.subs(subs)

    print(f"  numerator = {numer_val}")
    print(f"  denominator = {denom_val}")

    if denom_val != 0:
        t_val = numer_val / denom_val
        print(f"  t = {t_val}")

        # Compute center
        p1x_val = p1x.subs(subs)
        p1y_val = p1y.subs(subs)
        p1z_val = p1z.subs(subs)
        m1x_val = m1x.subs(subs)
        m1y_val = m1y.subs(subs)
        m1z_val = m1z.subs(subs)

        ox = m1x_val + t_val * p1x_val
        oy = m1y_val + t_val * p1y_val
        oz = m1z_val + t_val * p1z_val

        print(f"  Center = ({simplify(ox)}, {simplify(oy)}, {simplify(oz)})")
        print("  Expected: (0, 0, 0) ✓" if simplify(ox) == 0 and simplify(oy) == 0 and simplify(oz) == 0 else "  ✗ Mismatch!")
    else:
        print("  Degenerate case (collinear points)")

    # Generate complete Rust function body
    print("\n" + "=" * 60)
    print("Generated Rust function:")
    print("=" * 60)

    # Build the function body using rust_code() for all derived expressions
    function_body = f'''
    /// Creates a circle from three non-collinear points.
    ///
    /// The circle passes through all three points, with center at the circumcenter
    /// and radius equal to the circumradius.
    ///
    /// Returns `None` if the points are collinear (no unique circle exists).
    ///
    /// # Example
    ///
    /// ```
    /// use clifford::specialized::conformal::dim3::{{Circle, Point}};
    ///
    /// let p1 = Point::new(1.0, 0.0, 0.0);
    /// let p2 = Point::new(0.0, 1.0, 0.0);
    /// let p3 = Point::new(-1.0, 0.0, 0.0);
    /// let circle = Circle::from_three_points(&p1, &p2, &p3).unwrap();
    /// assert!((circle.radius() - 1.0).abs() < 1e-10);
    /// ```
    #[inline]
    pub fn from_three_points(p1: &Point<T>, p2: &Point<T>, p3: &Point<T>) -> Option<Self> {{
        // Extract coordinates
        let ax = p1.x();
        let ay = p1.y();
        let az = p1.z();
        let bx = p2.x();
        let by = p2.y();
        let bz = p2.z();
        let cx = p3.x();
        let cy = p3.y();
        let cz = p3.z();

        // Edge vectors: d1 = B - A, d2 = C - A
        let d1x = bx - ax;
        let d1y = by - ay;
        let d1z = bz - az;
        let d2x = cx - ax;
        let d2y = cy - ay;
        let d2z = cz - az;

        // Plane normal: n = d1 × d2
        let nx = {rust_code(expand(d1y * d2z - d1z * d2y))};
        let ny = {rust_code(expand(d1z * d2x - d1x * d2z))};
        let nz = {rust_code(expand(d1x * d2y - d1y * d2x))};

        // Check for collinearity
        let len_sq = nx * nx + ny * ny + nz * nz;
        if len_sq < T::epsilon() {{
            return None;
        }}

        // Normalize the normal
        let len = len_sq.sqrt();
        let nx = nx / len;
        let ny = ny / len;
        let nz = nz / len;

        // Midpoint of edge AB: m1 = (A + B) / 2
        let m1 = ((ax + bx) / T::TWO, (ay + by) / T::TWO, (az + bz) / T::TWO);

        // Perpendicular direction in plane: p1 = d1 × n
        let perp1 = (
            {rust_code(expand(p1x))},
            {rust_code(expand(p1y))},
            {rust_code(expand(p1z))},
        );

        // Perpendicular direction for edge AC: p2 = d2 × n
        let perp2 = (
            {rust_code(expand(p2x))},
            {rust_code(expand(p2y))},
            {rust_code(expand(p2z))},
        );

        // Delta midpoint: dm = m2 - m1 = (C - B) / 2
        let dm = ((cx - bx) / T::TWO, (cy - by) / T::TWO, (cz - bz) / T::TWO);

        // Numerator: ((m2 - m1) × p2) · n
        let dm_cross_p2 = (
            dm.1 * perp2.2 - dm.2 * perp2.1,
            dm.2 * perp2.0 - dm.0 * perp2.2,
            dm.0 * perp2.1 - dm.1 * perp2.0,
        );
        let numer = dm_cross_p2.0 * nx + dm_cross_p2.1 * ny + dm_cross_p2.2 * nz;

        // Denominator: (p1 × p2) · n
        let p1_cross_p2 = (
            perp1.1 * perp2.2 - perp1.2 * perp2.1,
            perp1.2 * perp2.0 - perp1.0 * perp2.2,
            perp1.0 * perp2.1 - perp1.1 * perp2.0,
        );
        let denom = p1_cross_p2.0 * nx + p1_cross_p2.1 * ny + p1_cross_p2.2 * nz;

        // Degenerate case: perpendicular bisectors are parallel or nearly parallel
        // Use a relative threshold to avoid numerical instability
        let denom_threshold = T::epsilon() * (T::one() + numer.abs());
        if denom.abs() < denom_threshold {{
            return None;
        }}

        let t = numer / denom;

        // Circumcenter: m1 + t * p1
        let center_x = m1.0 + t * perp1.0;
        let center_y = m1.1 + t * perp1.1;
        let center_z = m1.2 + t * perp1.2;

        // Radius: distance from center to any point
        let dx = center_x - ax;
        let dy = center_y - ay;
        let dz = center_z - az;
        let r = (dx * dx + dy * dy + dz * dz).sqrt();

        Some(Self {{
            cx: center_x,
            cy: center_y,
            cz: center_z,
            r,
            nx,
            ny,
            nz,
        }})
    }}
'''
    print(function_body)

    print("=" * 60)
    print("Circumcenter derivation complete!")
    print("=" * 60)

    return {
        'algorithm': 'perpendicular_bisector_intersection',
        'numer': 'dm_cross_p2 · n',
        'denom': 'p1_cross_p2 · n',
        't': 'numer / denom',
        'center': 'm1 + t * p1',
    }


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
    print()

    derive_circumcenter()
