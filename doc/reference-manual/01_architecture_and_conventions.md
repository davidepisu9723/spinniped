# Chapter 1 -- Architecture and conventions

## Table of contents

- [Declarative and numerical layers](#declarative-and-numerical-layers)
- [Shared conventions](#shared-conventions)

## Declarative and numerical layers

The public model definition consists of immutable records:

```text
Grid, CoordinateSystem, Material
ShaftProperty, BearingProperty, DiskProperty
ShaftElement, BearingElement, DiskElement
HarmonicLoad, RotatingLoad, Unbalance
SpeedDependentValue, RandomDistribution
```

These records store identifiers, connectivity, and Python-valued parameters.
`SpeedDependentValue` is a nested value record rather than a member of a
`ModelDefinition` list. Definition records do not own or calculate element
matrices. `SpeedDependentValue`, `RandomDistribution`, `Grid`,
`CoordinateSystem`, `Material`, and `ShaftProperty` validate their local
contracts during record construction. `BearingProperty` does the same for its
constant, random, and tabulated coefficients, while `DiskProperty` validates
its nonnegative mass, inertia, and damping parameters. `ShaftElement` validates
its local identity and connectivity contract; `BearingElement` validates its
optional two-grid connectivity and frame reference as well.
`ModelBuilder` resolves
random parameters, validates table placement and field-specific constraints,
calculates constant local matrices and unit-coefficient bearing bases using the
functions in `spinniped.stiffness`, `spinniped.mass`, `spinniped.damping`, and
`spinniped.gyroscopic`, transforms them, and assembles the global matrices.

Load records are resolved and validated with the other definition records but
do not contribute to the structural matrices. They remain attached to each
resolved definition for analysis-specific force assembly by the steady
response solvers.

Declarative record collections use lists so a definition can be assembled
incrementally. Atomic value sequences use tuples and reject mutable lists.
This applies to speed-table sequences, coordinate-system vectors, harmonic
force and moment vectors, and the model spin axis. The builder converts the
completed definition into tuples and NumPy arrays in `BuiltModel`.

Speed-dependent bearing terms remain attached to the built model as compact
per-sample contributions. Solvers interpolate their coefficients over the
requested rotor-speed operating points before numerical solution. This keeps
the permanent global arrays sample-first and three-dimensional.

A deterministic build has one realization. A stochastic build has $n_s$
Monte Carlo realizations. Global numerical arrays therefore always retain a
leading realization index:

$$
\mathbf{K}^{(s)},\ \mathbf{M}^{(s)},\ \mathbf{C}^{(s)},\ \mathbf{G}^{(s)},
\qquad s=0,1,\ldots,n_s-1.
$$

In code their shape is `(samples, ndof, ndof)`, including when `samples == 1`.

## Shared conventions

Spinniped implementation chapters use six mechanical degrees of freedom per
grid, dense sample-first matrix arrays, coherent units, right-handed Cartesian
frames, and the global equation-of-motion convention documented throughout
this reference manual.

The public records are implemented in `spinniped.records`. Validation,
stochastic resolution, speed-table validation and interpolation, coordinate
conversion, and assembly are implemented in `spinniped.builder`. Numerical
kernels remain stateless in
`spinniped.stiffness`, `spinniped.mass`, `spinniped.damping`, and
`spinniped.gyroscopic`.
