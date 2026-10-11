# Chapter 3 -- Speed-dependent values

## Table of contents

- [Defining a table](#defining-a-table)
- [Supported fields](#supported-fields)
- [Interpolation and operating speed](#interpolation-and-operating-speed)
- [Combining speed dependence and randomness](#combining-speed-dependence-and-randomness)
- [Model and result dimensions](#model-and-result-dimensions)

`SpeedDependentValue` defines a scalar coefficient at selected rotor angular
speeds. It is useful for bearing coefficients or applied loads obtained from
tests, CFD, or another calculation at a finite set of operating points.

## Defining a table

Pass separate tuples of rotor speeds and corresponding values:

```python
from spinniped import BearingProperty, SpeedDependentValue

kxx = SpeedDependentValue(
    speeds=(0.0, 100.0, 200.0),       # rad/s
    values=(1.00e8, 0.92e8, 0.80e8),
)

bearing = BearingProperty(
    id=1,
    kxx=kxx,
    kyy=kxx,
    cxx=100.0,
    cyy=100.0,
)
```

The builder requires at least two points. `speeds` and `values` must be
one-dimensional, have equal lengths, contain finite real numbers after random
references are resolved, and list speeds in strictly increasing order.
Boolean values are rejected as mechanical numbers.

## Supported fields

Speed dependence is deliberately limited to quantities that enter a matrix or
load vector linearly:

- all stiffness coefficients `kxx`, `kyy`, `kzz`, `kxy`, and `kyx` of
  `BearingProperty`;
- all damping coefficients `cxx`, `cyy`, `czz`, `cxy`, and `cyx` of
  `BearingProperty`;
- each component of `HarmonicLoad.force` and `HarmonicLoad.moment`;
- `RotatingLoad.force_magnitude`.

For example, a tabulated force component and circular-force magnitude can be
defined as follows:

```python
from spinniped import HarmonicLoad, RotatingLoad

fluid_force = SpeedDependentValue(
    speeds=(50.0, 100.0, 150.0),
    values=(8.0, 19.0, 35.0),
)

fixed_direction = HarmonicLoad(
    id=10,
    grid=2,
    force=(fluid_force, 0.0, 0.0),
)

rotating = RotatingLoad(
    id=11,
    grid=2,
    force_magnitude=fluid_force,
)
```

Speed tables are rejected in geometry, material, shaft, disk, phase, spin-axis,
distribution-parameter, and unbalance fields. In particular,
`Unbalance.mass_eccentricity` remains constant; its force already follows the
explicit physical law $m e\Omega^2$.

## Interpolation and operating speed

Values are linearly interpolated at the signed rotor angular speed $\Omega$ in
radians per second. This is not the independent forcing frequency $\omega$.
Consequently, a frequency-response analysis evaluates each table at its one
scalar `speed`; it does not reevaluate the table while sweeping
`frequencies`.

Exact endpoints are accepted. Extrapolation and endpoint clamping are not:
every requested speed must lie in the closed interval from the first to the
last tabulated speed. Because speed is signed, an analysis containing negative
speeds requires a table whose range also contains those negative values.

The solver validates and interpolates every required structural table before
starting the first numerical solution. Campbell and unbalance-response sweeps
precompute values for their complete `speeds` arrays. Modal, frequency, and
time response evaluate one scalar operating `speed`. Selected declarative
frequency-response loads are also assembled at that operating speed before
the frequency sweep begins. An out-of-range value therefore fails before a
partially completed analysis can be returned.

## Combining speed dependence and randomness

Speed interpolation and random sampling are distinct operations, but the
ordinate values may contain ordinary random-distribution references:

```python
uncertain_kxx = SpeedDependentValue(
    speeds=(0.0, 100.0, 200.0),
    values=((21, 0), (21, 1), (21, 2)),
)
```

Here distribution `21` is a three-component multivariate distribution. For
each realization, the builder samples it once, resolves all three ordinates,
and retains one deterministic table for that realization. The solver then
interpolates that realized table at every requested speed. It never draws a
new random value during a speed loop.

Reusing one random reference has the same reuse semantics as elsewhere in the
API: every occurrence receives the same draw within a realization. The speed
knot locations themselves must be deterministic numeric values.

## Model and result dimensions

Adding speed dependence does not add a permanent speed axis to `BuiltModel`.
Its matrices keep shape `(samples, ndof, ndof)`. The builder stores constant
base matrices and compact linear bearing contributions; the solver evaluates
temporary deltas only at the operating speeds required by one analysis:

$$
\mathbf{K}^{(s)}(\Omega)=
\mathbf{K}^{(s)}_0+\Delta\mathbf{K}^{(s)}(\Omega),
\qquad
\mathbf{C}^{(s)}(\Omega)=
\mathbf{C}^{(s)}_0+\Delta\mathbf{C}^{(s)}(\Omega).
$$

Selected declarative loads follow the analogous linear superposition

$$
\mathbf{F}^{(s)}(\Omega)=
\mathbf{F}^{(s)}_0+\Delta\mathbf{F}^{(s)}(\Omega).
$$

`frequency_response` accepts one scalar rotor `speed`. To evaluate a CFD load
table across several rotor operating speeds, call that analysis once per rotor
speed; each call may still sweep any number of forcing `frequencies`.

The existing result layouts are unchanged. A Campbell frequency array, for
example, still has shape `(samples, speeds, modes)`, while frequency response
still has shape `(samples, frequencies, free_dofs)`.
