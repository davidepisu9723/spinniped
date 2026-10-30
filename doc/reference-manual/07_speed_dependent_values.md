# Chapter 7 -- Speed-dependent values

## Table of contents

- [Record contract](#record-contract)
- [Allowed placements](#allowed-placements)
- [Resolution and validation](#resolution-and-validation)
- [Interpolation contract](#interpolation-contract)
- [Structural and load evaluation](#structural-and-load-evaluation)

## Record contract

`SpeedDependentValue` is a frozen, slotted dataclass in `spinniped.records`:

```python
SpeedDependentValue(
    speeds=(omega_0, omega_1, ...),
    values=(value_0, value_1, ...),
)
```

Both fields must be tuples. Mutable lists and other sequence types are
rejected at record construction rather than converted implicitly.

Despite the conventional lowercase variable names in the example, each
abscissa represents signed rotor angular speed $\Omega$ in rad/s, not forcing
angular frequency $\omega$. `speeds` contains numeric knots. `values` contains
numeric scalars or the same scalar and multivariate random references accepted
by other parameter fields.

The record has no interpolation or extrapolation selector. Its fixed contract
is piecewise-linear interpolation with no extrapolation.

## Allowed placements

The builder accepts the wrapper only where the tabulated scalar enters the
assembled quantity linearly:

| Record | Fields |
|---|---|
| `BearingProperty` | `kxx`, `kyy`, `kzz`, `kxy`, `kyx`, `cxx`, `cyy`, `czz`, `cxy`, `cyx` |
| `HarmonicLoad` | each component of `force` and `moment` |
| `RotatingLoad` | `force_magnitude` |

All other placements are rejected before sampling. Record constructors reject
tables in grids, coordinate systems, materials, shaft and disk properties,
load phases, and `Unbalance.mass_eccentricity`. Context-level validation covers
distribution parameters and `ModelDefinition.spin_axis`.

## Resolution and validation

Generic tuple traversal cannot safely distinguish a two-integer speed grid
from a random reference. The resolver therefore handles
`SpeedDependentValue` explicitly: it preserves `speeds` and recursively
resolves only the items in `values`. One distribution is sampled once per
realization through the normal sample cache, so repeated references retain
their established reuse semantics.

Record construction enforces all of the following local invariants:

- at least two knots;
- equal knot and ordinate lengths;
- finite real speeds and numeric ordinates, excluding Boolean values;
- strictly increasing, and therefore unique, speed knots.
- structurally valid random ordinate references.

When random ordinate references are resolved, constructing the resolved record
validates the resulting numeric values again. A `RotatingLoad` therefore
rejects known negative magnitude ordinates immediately and repeats that check
after random-reference resolution. The builder retains only checks that require
surrounding context, such as reference lookup and supported field placement.

## Interpolation contract

The dedicated interpolation helper accepts one resolved table, one nonempty
one-dimensional sequence of finite requested speeds, and a field path used in
diagnostics. It first validates the complete request against the closed knot
interval

$$
\Omega_0\leq\Omega_i\leq\Omega_{n-1},
$$

then applies one-dimensional linear interpolation. It raises if any requested
speed lies outside the interval; no partial interpolation, extrapolation, or
endpoint clamping occurs.

Signed speed is preserved. A table containing only nonnegative knots cannot
serve a reverse-speed operating point. Frequency-response forcing frequencies
are a separate independent variable and never select a table ordinate.

## Structural and load evaluation

For every tabulated bearing coefficient, the builder omits that coefficient
from the permanent base matrix and retains the resolved table, transformed
unit-coefficient matrix, affected global indices, and diagnostic path. If
$v_j(\Omega)$ is one interpolated coefficient and $\mathbf{B}_j$ its basis,
the solver constructs

$$
\Delta\mathbf{K}^{(s)}(\Omega_i)
=\sum_{j\in K}v_j^{(s)}(\Omega_i)\mathbf{B}_j^{(s)},
\qquad
\Delta\mathbf{C}^{(s)}(\Omega_i)
=\sum_{j\in C}v_j^{(s)}(\Omega_i)\mathbf{B}_j^{(s)}.
$$

Before the first numerical solve, all structural tables are interpolated over
the complete operating-speed request. The transient delta arrays have shape
`(samples, requested_speeds, ndof, ndof)` and are discarded after the
analysis; the permanent `BuiltModel` matrices retain shape
`(samples, ndof, ndof)`.

Declarative load components are evaluated from the resolved definition at the
frequency-response rotor `speed`, transformed into global coordinates, and
summed before the forcing-frequency loop. Loads do not modify the structural
base matrices or structural delta arrays. Constant and tabulated components
therefore implement

$$
\mathbf{F}^{(s)}(\Omega)=
\mathbf{F}^{(s)}_0+\Delta\mathbf{F}^{(s)}(\Omega).
$$

Frequency response exposes one scalar rotor `speed`, not a rotor-speed result
axis. A multi-speed CFD-load study invokes the analysis once per rotor speed;
within each call, all requested forcing frequencies reuse the same evaluated
load vector.
