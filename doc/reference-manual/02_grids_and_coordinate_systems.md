# Chapter 2 -- Grids and coordinate systems

## Table of contents

- [Coordinate-system construction](#coordinate-system-construction)
- [Nodal degrees of freedom](#nodal-degrees-of-freedom)
- [Global numbering](#global-numbering)

## Coordinate-system construction

The global frame is a right-handed Cartesian coordinate system
$(x,y,z)$. For the conventional straight-rotor configuration, $z$ is the
longitudinal shaft axis and $x$ and $y$ are transverse directions. Positive
rotations follow the right-hand rule.

Coordinate system ID `0` is reserved for the global frame. A
`CoordinateSystem` record defines a local frame using:

- `origin`: the global position of its origin;
- `x_axis`: the direction of its positive local x-axis, expressed globally;
- `xy_plane`: a second global direction lying in its local x-y plane.

The direction vectors need not have unit length, but they must be nonzero and
non-collinear. The builder orthonormalizes them into a rotation matrix. A
`Grid` stores coordinates in the frame selected by its `coordinate_system`
field, and the builder resolves them into global coordinates before assembly.

For each shaft, the element-local $z_e$ axis points from `grid_a` to `grid_b`.
The local matrices are transformed into the global frame. Thus the formulas in
the shaft theory note apply in the local element system even when the shaft is
not parallel to global $z$.

`ModelDefinition.spin_axis` supplies the global spin direction. If the supplied
nonzero vector is $\mathbf{s}$, the builder uses

$$
\widehat{\mathbf{s}}=\frac{\mathbf{s}}{\lVert\mathbf{s}\rVert},
\qquad
\boldsymbol{\Omega}=\Omega\widehat{\mathbf{s}},
$$

where $\Omega$ is angular speed in rad/s. A shaft's unit-speed gyroscopic
matrix is scaled by

$$
\gamma_e=\mathbf{e}_{z_e}\mathbin{\cdot}\widehat{\mathbf{s}}.
$$

This projection makes the assembled physical gyroscopic matrix invariant to a
reversal of the shaft's endpoint order. It also makes the contribution zero
when the shaft axis is perpendicular to the model spin axis.

## Nodal degrees of freedom

Every grid has six mechanical degrees of freedom:

$$
\mathbf{q}_i=
\begin{bmatrix}
x_i & y_i & z_i & \theta_{x_i} & \theta_{y_i} & \theta_{z_i}
\end{bmatrix}^{T}.
$$

| Local index | Symbol | Code name | Meaning |
|---:|---|---|---|
| 0 | $x_i$ | `x` | translation along x |
| 1 | $y_i$ | `y` | translation along y |
| 2 | $z_i$ | `z` | translation along z |
| 3 | $\theta_{x_i}$ | `tx` | rotation about x |
| 4 | $\theta_{y_i}$ | `ty` | rotation about y |
| 5 | $\theta_{z_i}$ | `tz` | rotation about z |

Translations and rotations in this table are expressed in the matrix's active
coordinate frame: local for a local element matrix and global after
transformation and assembly.

## Global numbering

Grid IDs are user identifiers and need not be consecutive. The builder maps
them to compact indices in `BuiltModel.grid_ids`. If $p$ is a compact grid
index, its global degrees of freedom are

$$
\begin{aligned}
x_p &= 6p, & y_p &= 6p+1, & z_p &= 6p+2,\\
\theta_{x_p} &= 6p+3, & \theta_{y_p} &= 6p+4,
& \theta_{z_p} &= 6p+5.
\end{aligned}
$$

For an element connecting compact indices $p$ and $q$, the global assembly
index vector is

$$
\mathbf{i}_e=
\begin{bmatrix}
6p & 6p+1 & 6p+2 & 6p+3 & 6p+4 & 6p+5 &
6q & 6q+1 & 6q+2 & 6q+3 & 6q+4 & 6q+5
\end{bmatrix}.
$$

Solver boundary conditions use these global integer DOF indices.
