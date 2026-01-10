"""Motor composition derivation using PGA geometric product.

Reference: https://rigidgeometricalgebra.org/wiki/index.php?title=Motor

PGA Signature: (3,0,1) meaning e1²=e2²=e3²=1, e0²=0

Approach:
    Use symbolic geometric product of two motors with proper PGA basis element
    multiplication rules.

Run:
    cd derivations
    uv run python -m clifford_derivations.motor
"""

from sympy import symbols, Matrix, sqrt, Rational, expand, simplify, collect, Symbol, rust_code
import signal


class TimeoutError(Exception):
    """Raised when a computation times out."""
    pass


def timeout_handler(signum, frame):
    raise TimeoutError("Computation timed out")


def with_timeout(seconds: int = 30):
    """Decorator to add timeout to a function."""
    def decorator(func):
        def wrapper(*args, **kwargs):
            old_handler = signal.signal(signal.SIGALRM, timeout_handler)
            signal.alarm(seconds)
            try:
                result = func(*args, **kwargs)
            finally:
                signal.alarm(0)
                signal.signal(signal.SIGALRM, old_handler)
            return result
        return wrapper
    return decorator


def to_rust(expr):
    """Convert a sympy expression to Rust code using sympy's built-in rust_code printer."""
    return rust_code(expand(expr))


# ============================================================
# PGA Geometric Product Tables
# ============================================================
#
# Basis elements in our motor:
#   1     (scalar)
#   e23   (bivector)
#   e31   (bivector)
#   e12   (bivector)
#   e01   (bivector with e0)
#   e02   (bivector with e0)
#   e03   (bivector with e0)
#   e0123 (pseudoscalar)
#
# Signature: e1²=e2²=e3²=1, e0²=0
# Anticommuting: ei*ej = -ej*ei for i≠j

# Compute geometric product of any two motor basis elements.
# Returns (coefficient, result_basis) where result is coeff * result_basis

def gp_bivectors():
    """Compute geometric product table for bivector pairs.

    Each bivector squares to a scalar:
    - e23² = e2e3e2e3 = -e2e2e3e3 = -1
    - e31² = e3e1e3e1 = -e3e3e1e1 = -1
    - e12² = e1e2e1e2 = -e1e1e2e2 = -1

    Bivector products (cyclic):
    - e23 * e31 = e2e3e3e1 = e2e1 = -e12
    - e31 * e23 = e3e1e2e3 = -e3e2e1e3 = e3e2e3e1 = -e2e1 = e12
    - e31 * e12 = e3e1e1e2 = e3e2 = -e23
    - e12 * e31 = e1e2e3e1 = -e1e3e2e1 = e3e2 = -e23... wait

    Let me be more careful:
    e12 * e31 = (e1 e2)(e3 e1)
              = e1 e2 e3 e1
              = -e1 e2 e1 e3  (swap e3 and e1, minus sign)
              = -e1 (-e1 e2) e3  (swap e2 and e1 inside, minus sign)
              = e1 e1 e2 e3
              = 1 * e2 e3
              = e23

    So e12 * e31 = +e23 (not -e23!)
    And e31 * e12 = -e23 (opposite sign)

    e23 * e12 = (e2 e3)(e1 e2)
              = e2 e3 e1 e2
              = -e2 e1 e3 e2  (swap e3 and e1)
              = e1 e2 e3 e2   (swap e2 and e1)
              = -e1 e2 e2 e3  (swap e3 and e2)
              = -e1 e3
              = e31

    So e23 * e12 = +e31
    And e12 * e23 = -e31

    e31 * e23 = (e3 e1)(e2 e3)
              = e3 e1 e2 e3
              = -e3 e1 e3 e2  (swap e2 and e3)
              = e3 e3 e1 e2   (swap e1 and e3)
              = 1 * e1 e2
              = e12

    So e31 * e23 = +e12
    And e23 * e31 = -e12
    """
    # Returns {(a, b): (coeff, result)} for a*b = coeff*result
    # where a, b, result are basis element names
    return {
        # Squares
        ('e23', 'e23'): (-1, '1'),
        ('e31', 'e31'): (-1, '1'),
        ('e12', 'e12'): (-1, '1'),
        # Cross products (note: OPPOSITE signs from quaternions!)
        ('e23', 'e31'): (-1, 'e12'),
        ('e31', 'e23'): (+1, 'e12'),
        ('e23', 'e12'): (+1, 'e31'),
        ('e12', 'e23'): (-1, 'e31'),
        ('e31', 'e12'): (-1, 'e23'),
        ('e12', 'e31'): (+1, 'e23'),
    }


@with_timeout(300)  # 5 minute timeout
def derive_composition():
    """Derive motor composition formula via symbolic geometric product.

    Convention: compose(m1, m2) applies m1 first, then m2.
    Geometric product: result = m2 * m1

    In PGA, the sandwich product (A * B) * p * (A * B)̃ applies A first, then B.
    So compose(self, other) = self * other to get "self first, then other".

    In the generated Rust code, we assign:
        M1 = other, M2 = self
    So that the result M2 * M1 = self * other.

    Uses proper PGA multiplication rules (NOT quaternion rules).

    Returns dict mapping component names to sympy expressions.
    """
    print("=== Motor Composition via PGA Geometric Product ===")
    print()
    print("Reference: https://rigidgeometricalgebra.org/wiki/index.php?title=Motor")
    print()
    print("PGA Signature: e1²=e2²=e3²=1, e0²=0")
    print()
    print("Convention: (A * B) * p * (A * B)̃ applies A first, then B")
    print("So compose(self, other) = self * other")
    print()

    # Motor 1 components - using our field names
    # In the Rust code, M1 = other (so it becomes the "second" operand in GP order)
    s1 = Symbol('s1', real=True)
    e23_1, e31_1, e12_1 = symbols('e23_1 e31_1 e12_1', real=True)
    e01_1, e02_1, e03_1 = symbols('e01_1 e02_1 e03_1', real=True)
    e0123_1 = Symbol('e0123_1', real=True)

    # Motor 2 components
    # In the Rust code, M2 = self (so it becomes the "first" operand in GP order)
    s2 = Symbol('s2', real=True)
    e23_2, e31_2, e12_2 = symbols('e23_2 e31_2 e12_2', real=True)
    e01_2, e02_2, e03_2 = symbols('e01_2 e02_2 e03_2', real=True)
    e0123_2 = Symbol('e0123_2', real=True)

    print("Motor 1 (maps to `other` in Rust):")
    print(f"  s1, e23_1, e31_1, e12_1, e01_1, e02_1, e03_1, e0123_1")
    print()
    print("Motor 2 (maps to `self` in Rust):")
    print(f"  s2, e23_2, e31_2, e12_2, e01_2, e02_2, e03_2, e0123_2")
    print()

    # ============================================================
    # Compute (m2 * m1) using PGA geometric product rules
    # ============================================================
    #
    # Motor M = s + e23*a + e31*b + e12*c + e01*d + e02*e + e03*f + e0123*g
    #
    # We need to compute m2 * m1 and collect terms by basis element.
    #
    # Key PGA products:
    # - Euclidean bivector squares: e23² = e31² = e12² = -1
    # - Euclidean bivector products (OPPOSITE of quaternion!):
    #     e23*e31 = -e12,  e31*e23 = +e12
    #     e31*e12 = -e23,  e12*e31 = +e23
    #     e12*e23 = -e31,  e23*e12 = +e31
    # - Degenerate products (involving e0): e01² = e02² = e03² = 0
    # - Mixed products: e23*e01 = e0123, e01*e23 = -e0123, etc.

    print("Computing geometric product m2 * m1...")
    print()

    # ---- Scalar component ----
    # From: s2*s1, e23_2*e23_1*(-1), e31_2*e31_1*(-1), e12_2*e12_1*(-1)
    # Degenerate terms: e01_2*e01_1=0, e02_2*e02_1=0, e03_2*e03_1=0
    s_new = s2*s1 - e23_2*e23_1 - e31_2*e31_1 - e12_2*e12_1

    # ---- e23 component ----
    # From: s2*e23_1, e23_2*s1
    # From: e31_2*e12_1 (e31*e12 = -e23, so coefficient is -1)
    # From: e12_2*e31_1 (e12*e31 = +e23, so coefficient is +1)
    e23_new = s2*e23_1 + e23_2*s1 - e31_2*e12_1 + e12_2*e31_1

    # ---- e31 component ----
    # From: s2*e31_1, e31_2*s1
    # From: e12_2*e23_1 (e12*e23 = -e31, so coefficient is -1)
    # From: e23_2*e12_1 (e23*e12 = +e31, so coefficient is +1)
    e31_new = s2*e31_1 + e31_2*s1 - e12_2*e23_1 + e23_2*e12_1

    # ---- e12 component ----
    # From: s2*e12_1, e12_2*s1
    # From: e23_2*e31_1 (e23*e31 = -e12, so coefficient is -1)
    # From: e31_2*e23_1 (e31*e23 = +e12, so coefficient is +1)
    e12_new = s2*e12_1 + e12_2*s1 - e23_2*e31_1 + e31_2*e23_1

    # ---- e01 component ----
    # From: s2*e01_1, e01_2*s1
    # From: e23_2*e0123_1 → e23*e0123 = e23*e0*e1*e2*e3 = e0*e1 = e01 (need to verify sign)
    # Actually: e0123 = e0*e1*e2*e3, so e23*e0123 = e2*e3*e0*e1*e2*e3
    #         = e0*e2*e3*e1*e2*e3 (moving e0 left, no sign change with e2,e3)
    #         = -e0*e2*e1*e3*e2*e3 = ... this is getting complex
    # Let me use: e23 = e2∧e3, e0123 = e0∧e1∧e2∧e3
    # e23 * e0123 = (e2e3)(e0e1e2e3) = e2e3e0e1e2e3
    #             = e0 e2e3e1e2e3 (commute e0 past e2e3)
    #             = e0 * (-e2e1e3e2e3) (swap e3,e1)
    #             = e0 * e1e2e3e2e3 (swap e2,e1)
    #             = e0 * e1 * (-e2e2e3e3) (swap inner e3,e2)
    #             = e0 * e1 * (-1*1) = -e0e1 = -e01
    # So e23 * e0123 gives -e01 contribution
    # e0123_2 * e23_1: e0123 * e23 = e0e1e2e3 * e2e3 = e0e1 * e2e3e2e3 = e0e1 * (-1) = -e01
    # So e0123 * e23 = -e01
    #
    # For e31: e31*e0123 = e3e1 * e0e1e2e3 = e0 * e3e1e1e2e3 = e0 * e3e2e3 = -e0 * e2e3e3 = -e0e2 = -e02
    # e0123 * e31 = e0e1e2e3 * e3e1 = e0e1e2 * e3e3e1 = e0e1e2e1 = -e0e1e1e2 = -e0e2 = -e02
    #
    # For e12: e12*e0123 = e1e2 * e0e1e2e3 = e0 * e1e2e1e2e3 = e0 * (-e1e1e2e2)e3 = -e0e3 = -e03
    # e0123 * e12 = e0e1e2e3 * e1e2 = e0e1e2e3e1e2 = e0(-e1e1)e2e3e2 = -e0e2e3e2 = e0e2e2e3 = e0e3 = e03
    # Wait that doesn't match. Let me be more careful.
    #
    # e0123 * e12 = (e0e1e2e3)(e1e2)
    #             = e0e1e2e3e1e2
    #             = e0e1e2 * e3e1 * e2 (swap nothing yet)
    #             = e0e1e2 * (-e1e3) * e2 (swap e3,e1)
    #             = -e0e1e2e1e3e2
    #             = -e0e1 * e2e1 * e3e2 (regrouping)
    #             = -e0e1 * (-e1e2) * e3e2
    #             = e0e1e1e2e3e2
    #             = e0 * 1 * e2e3e2
    #             = e0 * (-e2e2e3)
    #             = e0 * (-1) * e3
    #             = -e0e3 = -e03
    # So e0123 * e12 = -e03
    #
    # Let's also do:
    # e31_2 * e02_1: e31 * e02 = (e3e1)(e0e2) = e0 * e3e1e2 = e0 * (-e1e3e2) = e0 * e1e2e3 = e0123? No wait...
    # e31 * e02 = e3e1e0e2 = e0e3e1e2 = e0 * e3e1e2
    # e3e1e2: swap e1,e3 → -e1e3e2, swap e3,e2 → -e1(-e2e3) = e1e2e3
    # So e31 * e02 = e0 * e1e2e3 = e0e1e2e3 = e0123
    # Similarly for other cross-terms.
    #
    # For e01, the contributions are:
    # s2 * e01_1 (→ e01)
    # e01_2 * s1 (→ e01)
    # e12_2 * e02_1: e12 * e02 = e1e2e0e2 = e0e1e2e2 = e0e1 = e01
    # e02_2 * e12_1: e02 * e12 = e0e2e1e2 = e0 * e2e1e2 = e0 * (-e1e2e2) = -e0e1 = -e01
    # e31_2 * e03_1: e31 * e03 = e3e1e0e3 = e0e3e1e3 = e0 * e3e1e3 = e0 * (-e1e3e3) = -e0e1 = -e01
    # e03_2 * e31_1: e03 * e31 = e0e3e3e1 = e0e1 = e01
    # e23_2 * e0123_1: e23 * e0123 = -e01 (derived above)
    # e0123_2 * e23_1: e0123 * e23 = -e01 (derived above)
    e01_new = (s2*e01_1 + e01_2*s1
               + e12_2*e02_1 - e02_2*e12_1
               - e31_2*e03_1 + e03_2*e31_1
               - e23_2*e0123_1 - e0123_2*e23_1)

    # ---- e02 component ----
    # Similar analysis...
    # s2 * e02_1 (→ e02)
    # e02_2 * s1 (→ e02)
    # e23_2 * e03_1: e23 * e03 = e2e3e0e3 = e0e2e3e3 = e0e2 = e02
    # e03_2 * e23_1: e03 * e23 = e0e3e2e3 = e0 * e3e2e3 = e0 * (-e2e3e3) = -e0e2 = -e02
    # e12_2 * e01_1: e12 * e01 = e1e2e0e1 = e0e1e2e1 = e0 * e1e2e1 = e0 * (-e1e1e2) = -e0e2 = -e02
    # e01_2 * e12_1: e01 * e12 = e0e1e1e2 = e0e2 = e02
    # e31_2 * e0123_1: e31 * e0123 = -e02 (derived above)
    # e0123_2 * e31_1: e0123 * e31 = -e02 (derived above)
    e02_new = (s2*e02_1 + e02_2*s1
               + e23_2*e03_1 - e03_2*e23_1
               - e12_2*e01_1 + e01_2*e12_1
               - e31_2*e0123_1 - e0123_2*e31_1)

    # ---- e03 component ----
    # s2 * e03_1 (→ e03)
    # e03_2 * s1 (→ e03)
    # e31_2 * e01_1: e31 * e01 = e3e1e0e1 = e0e3e1e1 = e0e3 = e03
    # e01_2 * e31_1: e01 * e31 = e0e1e3e1 = e0 * e1e3e1 = e0 * (-e1e1e3) = -e0e3 = -e03
    # e23_2 * e02_1: e23 * e02 = e2e3e0e2 = e0e2e3e2 = e0 * e2e3e2 = e0 * (-e2e2e3) = -e0e3 = -e03
    # e02_2 * e23_1: e02 * e23 = e0e2e2e3 = e0e3 = e03
    # e12_2 * e0123_1: e12 * e0123 = -e03 (derived above)
    # e0123_2 * e12_1: e0123 * e12 = -e03 (derived above)
    e03_new = (s2*e03_1 + e03_2*s1
               + e31_2*e01_1 - e01_2*e31_1
               - e23_2*e02_1 + e02_2*e23_1
               - e12_2*e0123_1 - e0123_2*e12_1)

    # ---- e0123 component ----
    # s2 * e0123_1 (→ e0123)
    # e0123_2 * s1 (→ e0123)
    # e23_2 * e01_1: e23 * e01 = e2e3e0e1 = e0e2e3e1 = e0 * (-e1e2e3) = -e0e1e2e3 = -e0123? Wait:
    # e2e3e1 = -e2e1e3 = e1e2e3, so e23 * e01 = e0 * e1e2e3 = e0123
    # e01_2 * e23_1: e01 * e23 = e0e1e2e3 = e0123
    # e31_2 * e02_1: e31 * e02 = e0123 (derived above)
    # e02_2 * e31_1: e02 * e31 = e0e2e3e1 = e0 * e2e3e1 = e0 * (-e1e2e3) = -e0123? Wait:
    # e2e3e1: move e1 left: e2(-e1e3) = -e2e1e3 = e1e2e3
    # So e02 * e31 = e0 * e1e2e3 = e0123
    # e12_2 * e03_1: e12 * e03 = e1e2e0e3 = e0e1e2e3 = e0123
    # e03_2 * e12_1: e03 * e12 = e0e3e1e2 = e0 * e3e1e2 = e0 * (-e1e3e2) = e0 * e1e2e3 = e0123
    # Wait, that means all cross-terms give +e0123, but that can't be right for anticommutation.
    # Let me redo e01 * e23:
    # e01 * e23 = (e0e1)(e2e3) = e0e1e2e3 = e0123 ✓
    # e23 * e01 = (e2e3)(e0e1) = e2e3e0e1 = e0 * e2e3e1 = e0 * (-e2e1e3) = e0 * e1e2e3 = e0123
    # Hmm so they're both +e0123. That's because e0 anticommutes with e1,e2,e3 so moving it
    # to the front introduces a sign, but then the remaining indices are the same.
    #
    # Similarly:
    # e02 * e31 = e0e2 * e3e1 = e0e2e3e1 = e0 * e2e3e1 = e0 * (-e1e2e3 ... wait no:
    # e2e3e1: swap e3,e1 → e2(-e1e3) = -e2e1e3, swap e2,e1 → -(-e1e2)e3 = e1e2e3
    # So e02 * e31 = e0 * e1e2e3 = e0123 ✓
    # e31 * e02 = e3e1 * e0e2 = e0 * e3e1e2
    # e3e1e2: swap e1,e3 → -e1e3e2, swap e3,e2 → -e1(-e2e3) = e1e2e3
    # So e31 * e02 = e0 * e1e2e3 = e0123 ✓
    #
    # e03 * e12 = e0e3 * e1e2 = e0e3e1e2 = e0 * e3e1e2 = e0 * e1e2e3 = e0123 ✓
    # e12 * e03 = e1e2 * e0e3 = e0 * e1e2e3 = e0123 ✓
    #
    # So ALL cross-terms between translation bivectors and rotation bivectors give +e0123.
    # That's 6 terms: e23*e01, e01*e23, e31*e02, e02*e31, e12*e03, e03*e12
    e0123_new = (s2*e0123_1 + e0123_2*s1
                 + e23_2*e01_1 + e01_2*e23_1
                 + e31_2*e02_1 + e02_2*e31_1
                 + e12_2*e03_1 + e03_2*e12_1)

    print("Simplifying expressions (may take a few minutes)...")

    result = {
        's': simplify(s_new),
        'e23': simplify(e23_new),
        'e31': simplify(e31_new),
        'e12': simplify(e12_new),
        'e01': simplify(e01_new),
        'e02': simplify(e02_new),
        'e03': simplify(e03_new),
        'e0123': simplify(e0123_new),
    }

    print()
    print("=== Simplified Results ===")
    for name, expr in result.items():
        print(f"{name} = {expr}")

    return result


def generate_rust_compose(composition_result):
    """Generate Rust code for Motor::compose from derived sympy expressions."""
    print()
    print("=== Generated Rust Code for Motor::compose ===")
    print()

    # Variable substitution for Rust: e23_1 -> b23_1, etc.
    subs = {
        's1': 's1', 's2': 's2',
        'e23_1': 'b23_1', 'e31_1': 'b31_1', 'e12_1': 'b12_1',
        'e23_2': 'b23_2', 'e31_2': 'b31_2', 'e12_2': 'b12_2',
        'e01_1': 'd01_1', 'e02_1': 'd02_1', 'e03_1': 'd03_1',
        'e01_2': 'd01_2', 'e02_2': 'd02_2', 'e03_2': 'd03_2',
        'e0123_1': 'i1', 'e0123_2': 'i2',
    }

    def rust_expr(expr):
        """Convert sympy expression to Rust, substituting variable names."""
        s = to_rust(expr)
        for old, new in subs.items():
            s = s.replace(old, new)
        return s

    print('''    /// Composes two motors: apply `self` first, then `other`.
    ///
    /// Returns the geometric product `self * other`.
    ///
    /// In PGA, the sandwich product (A * B) * p * (A * B)̃ applies A first, then B.
    ///
    /// # Derivation
    ///
    /// Derived via PGA geometric product from:
    /// https://rigidgeometricalgebra.org/wiki/index.php?title=Motor
    ///
    /// See: `derivations/src/clifford_derivations/motor.py`
    #[inline]
    pub fn compose(&self, other: &Self) -> Self {
        // M1 = other, M2 = self, result = M2 * M1 = self * other
        let s1 = other.s;
        let b23_1 = other.e23;
        let b31_1 = other.e31;
        let b12_1 = other.e12;
        let d01_1 = other.e01;
        let d02_1 = other.e02;
        let d03_1 = other.e03;
        let i1 = other.e0123;

        let s2 = self.s;
        let b23_2 = self.e23;
        let b31_2 = self.e31;
        let b12_2 = self.e12;
        let d01_2 = self.e01;
        let d02_2 = self.e02;
        let d03_2 = self.e03;
        let i2 = self.e0123;

        // Geometric product: self * other = M2 * M1''')

    print(f"        let s = {rust_expr(composition_result['s'])};")
    print(f"        let e23 = {rust_expr(composition_result['e23'])};")
    print(f"        let e31 = {rust_expr(composition_result['e31'])};")
    print(f"        let e12 = {rust_expr(composition_result['e12'])};")
    print(f"        let e01 = {rust_expr(composition_result['e01'])};")
    print(f"        let e02 = {rust_expr(composition_result['e02'])};")
    print(f"        let e03 = {rust_expr(composition_result['e03'])};")
    print(f"        let e0123 = {rust_expr(composition_result['e0123'])};")

    print('''
        Self { s, e23, e31, e12, e01, e02, e03, e0123 }
    }''')


def derive_inverse():
    """Derive motor inverse formula.

    M⁻¹ = M̃ / ||M||²

    where M̃ is the reverse and ||M||² = s² + e23² + e31² + e12²
    """
    print()
    print("=== Deriving Motor::inverse ===")
    print()

    s = Symbol('s', real=True)
    e23, e31, e12 = symbols('e23 e31 e12', real=True)
    e01, e02, e03 = symbols('e01 e02 e03', real=True)
    e0123 = Symbol('e0123', real=True)

    # Norm squared (rotation part only contributes)
    norm_sq = s**2 + e23**2 + e31**2 + e12**2

    # Reverse: negate bivector and pseudoscalar parts
    # M̃ = s - e23 - e31 - e12 - e01 - e02 - e03 - e0123

    # Inverse = reverse / norm_sq
    result = {
        's': s / norm_sq,
        'e23': -e23 / norm_sq,
        'e31': -e31 / norm_sq,
        'e12': -e12 / norm_sq,
        'e01': -e01 / norm_sq,
        'e02': -e02 / norm_sq,
        'e03': -e03 / norm_sq,
        'e0123': -e0123 / norm_sq,
    }

    print(f"||M||² = {norm_sq}")
    print()
    print("M⁻¹ = M̃ / ||M||²:")
    for name, expr in result.items():
        print(f"  {name} = {expr}")

    return result, norm_sq


def generate_rust_inverse(inverse_result, norm_sq):
    """Generate Rust code for Motor::inverse from derived sympy expressions."""
    print()
    print("=== Generated Rust Code for Motor::inverse ===")
    print()

    # For inverse, we use self.X directly
    def rust_expr(expr):
        s = to_rust(expr)
        # Replace bare variable names with self.X
        for var in ['s', 'e23', 'e31', 'e12', 'e01', 'e02', 'e03', 'e0123']:
            # Only replace if it's a standalone variable (not part of another name)
            import re
            s = re.sub(r'\b' + var + r'\b', f'self.{var}', s)
        return s

    print('''    /// Returns the inverse of this motor: M⁻¹ = M̃ / ||M||².
    ///
    /// # Derivation
    ///
    /// See: `derivations/src/clifford_derivations/motor.py`
    #[inline]
    pub fn inverse(&self) -> Self {''')

    print(f"        let norm_sq = {rust_expr(norm_sq)};")
    print()
    print("        Self {")
    print(f"            s: {rust_expr(inverse_result['s'])},")
    print(f"            e23: {rust_expr(inverse_result['e23'])},")
    print(f"            e31: {rust_expr(inverse_result['e31'])},")
    print(f"            e12: {rust_expr(inverse_result['e12'])},")
    print(f"            e01: {rust_expr(inverse_result['e01'])},")
    print(f"            e02: {rust_expr(inverse_result['e02'])},")
    print(f"            e03: {rust_expr(inverse_result['e03'])},")
    print(f"            e0123: {rust_expr(inverse_result['e0123'])},")
    print("        }")
    print("    }")


if __name__ == "__main__":
    try:
        # Derive composition
        composition = derive_composition()
        generate_rust_compose(composition)

        # Derive inverse
        inverse, norm_sq = derive_inverse()
        generate_rust_inverse(inverse, norm_sq)

    except TimeoutError:
        print("Derivation timed out")
