# Clifford Derivations

Symbolic derivations for the Clifford geometric algebra library using SymPy.

## Reference

- [Rigid Geometric Algebra Wiki - Motor](https://rigidgeometricalgebra.org/wiki/index.php?title=Motor)
  - "Conversion from Motor to Matrix"
  - "Conversion from Matrix to Motor"

## Usage

```bash
# Run the derivations
cd derivations
uv run python -m clifford_derivations.motor
```

## Package Structure

```
src/clifford_derivations/
  __init__.py     # Package exports
  motor.py        # Motor composition derivation via matrix conversion
```

## Convention (Matching clifford::specialized::projective::dim3)

- **Motor field order**: s, e23, e31, e12, e01, e02, e03, e0123
- **Compose order**: `a.compose(&b)` applies `a` first, then `b` (returns `a * b` in GP terms)
- **Sandwich convention**: In PGA, `(A * B) * p * (A * B)̃` applies `A` first, then `B`

## Derivation Approach

Motor composition is derived by:

1. **Define PGA basis products**: Compute all basis element products using PGA signature (3,0,1)
2. **Symbolic GP expansion**: Compute the geometric product of two motors symbolically
3. **Simplify**: Use SymPy's `simplify()` to get the final expressions
4. **Generate Rust code**: Use `sympy.rust_code()` to generate the implementation

Key PGA bivector products (differ from quaternions!):
- `e23*e31 = -e12`, `e31*e23 = +e12`
- `e31*e12 = -e23`, `e12*e31 = +e23`
- `e12*e23 = -e31`, `e23*e12 = +e31`

## Critical Rule: Generate Rust from SymPy

**Never hardcode algebraic formulas.** We cannot trust manual algebra—we can trust SymPy.

All Rust code generation uses `sympy.rust_code()`:

```python
from sympy import symbols, expand, rust_code

x, y = symbols('x y')
expr = x**2 + 2*x*y + y**2
rust_output = rust_code(expand(expr))  # Uses sympy's rust_code printer
```

## Available Functions

### `derive_composition()`
Derive motor composition via PGA geometric product. Returns dict of sympy expressions.

### `derive_inverse()`
Derive motor inverse formula. Returns tuple of (result_dict, norm_sq_expr).

### `generate_rust_compose(composition_result)`
Generate Rust code for `Motor::compose` from sympy expressions.

### `generate_rust_inverse(inverse_result, norm_sq)`
Generate Rust code for `Motor::inverse` from sympy expressions.

### `to_rust(expr)`
Convert a sympy expression to Rust code using sympy's built-in `rust_code` printer.

### `gp_bivectors()`
Returns the geometric product table for PGA bivector pairs.

## Timeout Handling

Use `@with_timeout(seconds)` for potentially slow operations:

```python
from clifford_derivations import with_timeout

@with_timeout(60)
def my_derivation():
    ...
```
