# Chapter 5 -- Elements

## Table of contents

- [Element record types](#element-record-types)
- [Two-grid shaft ordering](#two-grid-shaft-ordering)
- [Disk and bearing matrices](#disk-and-bearing-matrices)
- [Sign and symmetry conventions](#sign-and-symmetry-conventions)

## Element record types

The `Element` union in `spinniped.records` contains `ShaftElement`,
`BearingElement`, and `DiskElement`. `LumpedMassElement` is an alias for
`DiskElement`. Element records store connectivity and property references; they
do not store numerical matrices.

The implemented matrix contributions are:

| Element | Stiffness | Mass | Damping | Gyroscopic |
|---|---|---|---|---|
| Shaft | Elastic shaft matrix | Consistent shaft matrix | Mass-proportional | Consistent lateral rotary-inertia matrix projected onto `spin_axis` |
| Bearing | Local translational coefficients, to ground or between two grids | Zero | Local translational coefficients, to ground or between two grids | Zero |
| Disk/lumped mass | Zero | Oriented concentrated mass and inertia | Mass-proportional | Oriented polar inertia projected onto `spin_axis` |

## Two-grid shaft ordering

For a shaft connecting grids $a$ and $b$, the element vector is

$$
\mathbf{q}_e=
\begin{bmatrix}
x_a & y_a & z_a & \theta_{x_a} & \theta_{y_a} & \theta_{z_a} &
x_b & y_b & z_b & \theta_{x_b} & \theta_{y_b} & \theta_{z_b}
\end{bmatrix}^{T}.
$$

In code notation:

```text
xa, ya, za, txa, tya, tza, xb, yb, zb, txb, tyb, tzb
```

The zero-based local index map is:

| Index | DOF | Index | DOF |
|---:|---|---:|---|
| 0 | $x_a$ | 6 | $x_b$ |
| 1 | $y_a$ | 7 | $y_b$ |
| 2 | $z_a$ | 8 | $z_b$ |
| 3 | $\theta_{x_a}$ | 9 | $\theta_{x_b}$ |
| 4 | $\theta_{y_a}$ | 10 | $\theta_{y_b}$ |
| 5 | $\theta_{z_a}$ | 11 | $\theta_{z_b}$ |

Many references use plane-separated orderings. Such matrices must be permuted
before use with Spinniped's local-matrix functions or assembly convention.

## Disk and bearing matrices

A disk is attached to one grid. A bearing either connects one grid to ground or
connects two grids through their relative motion. Both use the six-DOF nodal
order from Section 3.

`DiskElement.coordinate_system` and `BearingElement.coordinate_system` select
the local frame in which their property coefficients are defined. The builder
rotates their matrices into the global frame. `Grid.coordinate_system`
transforms only the grid position and does not implicitly orient attached
elements.

For a disk, $m_d$ denotes mass, $I_d$ the diametral mass moment of inertia, and
$I_p$ the polar mass moment about the local spin axis. Its lumped mass matrix is

$$
\mathbf{M}_d=
\operatorname{diag}(m_d,m_d,m_d,I_d,I_d,I_p).
$$

In its selected local frame, the disk's canonical unit-speed gyroscopic matrix
has the nonzero rotational block

$$
\mathbf{G}_d[\theta_x,\theta_y]=
\begin{bmatrix}
0 & I_p\\
-I_p & 0
\end{bmatrix}.
$$

If $\mathbf{e}_{z_d}$ is the disk coordinate system's local z-axis, the builder
scales this matrix by

$$
\gamma_d=\mathbf{e}_{z_d}\mathbin{\cdot}\widehat{\mathbf{s}}
$$

before rotating it into the global frame.

A bearing contributes translational stiffness and viscous damping:

$$
\mathbf{K}_b[0{:}3,0{:}3]=
\begin{bmatrix}
k_{xx} & k_{xy} & 0\\
k_{yx} & k_{yy} & 0\\
0 & 0 & k_{zz}
\end{bmatrix},
$$

$$
\mathbf{C}_b[0{:}3,0{:}3]=
\begin{bmatrix}
c_{xx} & c_{xy} & 0\\
c_{yx} & c_{yy} & 0\\
0 & 0 & c_{zz}
\end{bmatrix}.
$$

Cross-coupled bearing coefficients need not make these matrices symmetric.
For a bearing joining grids $a$ and $b$, either local matrix
$\mathbf{A}_b\in\{\mathbf{K}_b,\mathbf{C}_b\}$ is expanded to

$$
\begin{bmatrix}
\mathbf{A}_b & -\mathbf{A}_b\\
-\mathbf{A}_b & \mathbf{A}_b
\end{bmatrix},
$$

which acts on the relative displacement or velocity. Omitting `grid_b`
connects the first grid to ground.

The current bearing formulation has zero mass and gyroscopic matrices; the
current disk formulation has zero elastic stiffness.

## Sign and symmetry conventions

Unless a formulation states otherwise:

- displacements are positive along positive coordinate axes;
- rotations are positive by the right-hand rule;
- conservative shaft stiffness and mass matrices are symmetric;
- ideal gyroscopic matrices are skew-symmetric,
  $\mathbf{G}^{T}=-\mathbf{G}$;
- bearing stiffness and damping may be non-symmetric when cross-coupled
  coefficients differ.

The shaft bending convention is

$$
\frac{dx}{dz}=\theta_y,
\qquad
\frac{dy}{dz}=-\theta_x.
$$

It determines the signs of displacement-rotation coupling terms in the two
bending planes.
