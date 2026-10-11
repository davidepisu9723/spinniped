# Chapter 5 -- Built models

## Table of contents

- [Overview](#overview)
- [Speed-dependent contributions](#speed-dependent-contributions)

## Overview

`ModelBuilder.build` returns a `BuiltModel`. Deterministic and stochastic builds
use the same layout:

| Attribute | Shape or type | Meaning |
|---|---|---|
| `K`, `stiffness` | `(samples, ndof, ndof)` | Global base stiffness matrices |
| `M`, `mass` | `(samples, ndof, ndof)` | Global mass matrices |
| `C`, `damping` | `(samples, ndof, ndof)` | Global base viscous damping matrices |
| `G`, `gyroscopic` | `(samples, ndof, ndof)` | Global unit-speed gyroscopic matrices |
| `coordinates` | `(samples, grids, 3)` | Resolved global grid coordinates |
| `definitions` | tuple | Fully resolved deterministic definitions |
| `grid_ids` | tuple | User grid IDs in assembly order |
| `samples` | integer | Number of realizations |
| `ndof` | integer | Number of global degrees of freedom |
| `stochastic` | Boolean | Whether random sampling was used |
| `seed` | integer or `None` | Requested random seed |

Even a deterministic model has a leading sample dimension of length one. For
example, its base stiffness matrix is `model.K[0]`, not `model.K`.

Every grid has six degrees of freedom in this order:

```text
x, y, z, tx, ty, tz
```

If a grid has compact assembly index `p`, its global DOFs are `6*p` through
`6*p + 5`. User grid IDs may be non-consecutive; use `model.grid_ids` to see
their assembly order.

## Speed-dependent contributions

Without any `SpeedDependentValue`, `K` and `C` are the complete assembled
matrices. When bearing coefficients are tabulated, `K` and `C` instead contain
the constant base terms. The model also retains compact, implementation-owned
unit-coefficient contributions. The solver combines them at a requested rotor
speed:

$$
\mathbf{K}^{(s)}(\Omega)=
\mathbf{K}^{(s)}_0+\Delta\mathbf{K}^{(s)}(\Omega),
\qquad
\mathbf{C}^{(s)}(\Omega)=
\mathbf{C}^{(s)}_0+\Delta\mathbf{C}^{(s)}(\Omega).
$$

`M` remains the complete mass matrix, and `G` remains the complete unit-speed
gyroscopic matrix. The permanent matrix arrays do not acquire a speed
dimension. A solver call may temporarily evaluate deltas at all requested
speeds, but it neither mutates the base matrices nor accumulates a change from
the preceding speed.

The resulting equation at operating speed $\Omega$ is

$$
\mathbf{M}\ddot{\mathbf{q}}+
(\mathbf{C}(\Omega)+\Omega\mathbf{G})\dot{\mathbf{q}}+
\mathbf{K}(\Omega)\mathbf{q}=\mathbf{f}(\Omega).
$$

`G` is a unit-speed matrix and is multiplied by angular speed `Omega` in the
Campbell, frequency-response, unbalance-response, and time-response analyses.
Modal analysis can also evaluate speed-dependent stiffness at its scalar
operating `speed`, although it continues to omit damping and gyroscopic terms.
Speeds are expressed in radians per second.
