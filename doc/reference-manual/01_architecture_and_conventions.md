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
```

These records store identifiers, connectivity, and Python-valued parameters.
They do not own or calculate element matrices. `ModelBuilder` resolves all
parameters, calculates local matrices using the functions in
`spinniped.stiffness`, `spinniped.mass`, `spinniped.damping`, and
`spinniped.gyroscopic`, transforms them, and assembles the global matrices.

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
stochastic resolution, coordinate conversion, and assembly are implemented in
`spinniped.builder`. Numerical kernels remain stateless in
`spinniped.stiffness`, `spinniped.mass`, `spinniped.damping`, and
`spinniped.gyroscopic`.
