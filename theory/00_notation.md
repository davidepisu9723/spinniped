# Notation and Conventions

This note defines the coordinate, degree-of-freedom, matrix, and realization
conventions used by Spinniped. The same conventions apply to deterministic and
stochastic models.

## 1. Declarative records and numerical models

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

## 2. Coordinate systems

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

## 3. Nodal degrees of freedom

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

## 4. Two-grid shaft ordering

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

## 5. Global degree-of-freedom numbering

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

## 6. Matrix notation and equation of motion

Element matrices carry the subscript $e$:

$$
\mathbf{M}_e,\quad \mathbf{C}_e,\quad
\mathbf{G}_e,\quad \mathbf{K}_e.
$$

Their global assembled counterparts omit the subscript:

$$
\mathbf{M},\quad \mathbf{C},\quad
\mathbf{G},\quad \mathbf{K}.
$$

For realization $s$, Spinniped uses the second-order rotor equation

$$
\mathbf{M}^{(s)}\ddot{\mathbf{q}}^{(s)}+
\left(\mathbf{C}^{(s)}+\Omega\mathbf{G}^{(s)}\right)
\dot{\mathbf{q}}^{(s)}+
\mathbf{K}^{(s)}\mathbf{q}^{(s)}=
\mathbf{f}^{(s)}.
$$

The builder stores a global unit-speed gyroscopic matrix after applying each
element's spin-axis projection. The solver supplies the scalar factor $\Omega$
in rad/s for speed-dependent analyses.

For the undamped real modal analysis, the solved problem is

$$
\mathbf{K}^{(s)}\boldsymbol{\phi}^{(s)}=
\lambda^{(s)}\mathbf{M}^{(s)}\boldsymbol{\phi}^{(s)},
\qquad
f^{(s)}=\frac{\sqrt{\lambda^{(s)}}}{2\pi}.
$$

### 6.1 Mode correspondence

Eigenvalue order alone is not a reliable physical mode label across stochastic
samples or through a Campbell speed sweep. With `track_modes=True`, Spinniped
uses the Modal Assurance Criterion between a reference vector
$\boldsymbol{\phi}_r$ and candidate $\boldsymbol{\phi}_c$:

$$
\operatorname{MAC}(\boldsymbol{\phi}_r,\boldsymbol{\phi}_c)=
\frac{\left|\boldsymbol{\phi}_r^H\boldsymbol{\phi}_c\right|^2}
{\left(\boldsymbol{\phi}_r^H\boldsymbol{\phi}_r\right)
 \left(\boldsymbol{\phi}_c^H\boldsymbol{\phi}_c\right)}.
$$

A one-to-one assignment maximizes total MAC. Modal samples are matched to
sample zero. Campbell modes are matched first between consecutive speeds, then
between each stochastic sample and sample zero at the same speed. After
assignment, each candidate is multiplied by a sign or complex phase factor so
its overlap with the reference is real and nonnegative. All comparisons use
the reduced free-DOF vectors; Campbell matching uses the displacement half of
the state eigenvector.

Mode counts must be consistent. The solver raises rather than returning ragged
object arrays when a requested count is unavailable or result shapes differ.

## 7. Disk and bearing matrices

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

## 8. Sign and symmetry conventions

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

## 9. Common symbols and units

Spinniped assumes a coherent unit system. The documentation and examples use
SI units.

| Symbol | Meaning | SI unit |
|---|---|---|
| $L$ | element length | m |
| $A$ | cross-sectional area | $\mathrm{m^2}$ |
| $E$ | Young's modulus | Pa |
| $G_s$ | shear modulus | Pa |
| $\rho$ | mass density | $\mathrm{kg/m^3}$ |
| $I$ | transverse geometric second moment | $\mathrm{m^4}$ |
| $J$ | geometric polar second moment/torsional constant | $\mathrm{m^4}$ |
| $m_d$ | disk or lumped mass | kg |
| $I_d$, $I_p$ | mass moments of inertia | $\mathrm{kg\,m^2}$ |
| $\Omega$ | spin speed | rad/s |
| $\omega$ | natural or forcing circular frequency | rad/s |
| $f$ | cyclic frequency | Hz |
| $k$ | translational stiffness | N/m |
| $c$ | translational viscous damping | N s/m |

Angles are dimensionless in the equations but represent radians. When a symbol
such as $J$ could mean either geometric or mass polar inertia, the surrounding
formulation must state which quantity is intended.

## 10. Implementation rule

Every new local matrix must declare its active coordinate frame and DOF order.
Two-grid shaft matrices must ultimately be expressed in the 12-DOF ordering in
Section 4, and one-grid matrices in the six-DOF ordering in Section 3, before
the builder assembles them.
