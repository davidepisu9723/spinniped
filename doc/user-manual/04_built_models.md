# Chapter 4 -- Built models

## Table of contents

- [Overview](#overview)

## Overview

`ModelBuilder.build` returns a `BuiltModel`. Deterministic and stochastic builds
use the same layout:

| Attribute | Shape or type | Meaning |
|---|---|---|
| `K`, `stiffness` | `(samples, ndof, ndof)` | Global stiffness matrices |
| `M`, `mass` | `(samples, ndof, ndof)` | Global mass matrices |
| `C`, `damping` | `(samples, ndof, ndof)` | Global viscous damping matrices |
| `G`, `gyroscopic` | `(samples, ndof, ndof)` | Global unit-speed gyroscopic matrices |
| `coordinates` | `(samples, grids, 3)` | Resolved global grid coordinates |
| `definitions` | tuple | Fully resolved deterministic definitions |
| `grid_ids` | tuple | User grid IDs in assembly order |
| `samples` | integer | Number of realizations |
| `ndof` | integer | Number of global degrees of freedom |
| `stochastic` | Boolean | Whether random sampling was used |
| `seed` | integer or `None` | Requested random seed |

Even a deterministic model has a leading sample dimension of length one. For
example, its stiffness matrix is `model.K[0]`, not `model.K`.

Every grid has six degrees of freedom in this order:

```text
x, y, z, tx, ty, tz
```

If a grid has compact assembly index `p`, its global DOFs are `6*p` through
`6*p + 5`. User grid IDs may be non-consecutive; use `model.grid_ids` to see
their assembly order.

The matrices implement the speed-dependent equation

$$
\mathbf{M}\ddot{\mathbf{q}}+
(\mathbf{C}+\Omega\mathbf{G})\dot{\mathbf{q}}+
\mathbf{K}\mathbf{q}=\mathbf{f}.
$$

`G` is a unit-speed matrix and is multiplied by angular speed `Omega` in the
Campbell, frequency-response, and time-response analyses. Speeds are expressed
in radians per second.
