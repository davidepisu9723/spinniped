# Notation and Conventions

This document defines the notation used in the Spinniped theoretical notes and in the finite element implementation.

The purpose of this file is to keep symbols, coordinate directions, degree-of-freedom ordering, and matrix conventions consistent across all theory documents.

---

## 1. Coordinate system

Spinniped uses a right-handed Cartesian coordinate system:

$$
\left( x, y, z \right)
$$

where:

- $z$ is the rotor longitudinal axis, unless stated otherwise;
- $x$ and $y$ are the two transverse directions;
- positive rotations follow the right-hand rule.

For the default rotor configuration, the spin vector is aligned with the positive $z$ axis:

$$
\boldsymbol{\Omega} = \Omega \, \mathbf{e}_z
$$

where:

- $\Omega$ is the rotor angular speed;
- $\mathbf{e}_z$ is the unit vector along the $z$ axis.

If a different rotation axis is used, the element matrices must either be formulated in that local axis system or transformed consistently before global assembly.

---

## 2. Nodal degrees of freedom

Each shaft node has six mechanical degrees of freedom:

$$
\mathbf{q}_i =
\begin{bmatrix}
 x_i \\
 y_i \\
 z_i \\
 \theta_{x_i} \\
 \theta_{y_i} \\
 \theta_{z_i}
\end{bmatrix}
$$

where:

| Symbol | Code name | Meaning |
|---|---:|---|
| $x_i$ | `x` | transverse displacement of node $i$ along the global $x$ direction |
| $y_i$ | `y` | transverse displacement of node $i$ along the global $y$ direction |
| $z_i$ | `z` | transverse displacement of node $i$ along the global $z$ direction |
| $\theta_{x_i}$ | `tx` | small rotation of node $i$ about the global $x$ axis |
| $\theta_{y_i}$ | `ty` | small rotation of node $i$ about the global $y$ axis |
| $\theta_{z_i}$ | `tz` | small rotation of node $i$ about the global $z$ axis |

The nodal unknown vector is therefore ordered as:

$$
\mathbf{q}_i =
\begin{bmatrix}
 x_i & y_i & z_i & \theta_{x_i} & \theta_{y_i} & \theta_{z_i}
\end{bmatrix}^T
$$

---

## 3. Two-node element degree-of-freedom ordering

For a two-node shaft element connecting node 0 and node 1, Spinniped uses the local element vector:

$$
\mathbf{q}_e =
\begin{bmatrix}
 x_0 \\
 y_0 \\
 z_0 \\
 \theta_{x_0} \\
 \theta_{y_0} \\
 \theta_{z_0} \\
 x_1 \\
 y_1 \\
 z_1 \\
 \theta_{x_1} \\
 \theta_{y_1} \\
 \theta_{z_1}
\end{bmatrix}
$$

Equivalently:

$$
\mathbf{q}_e =
\begin{bmatrix}
 x_0 & y_0 & z_0 & \theta_{x_0} & \theta_{y_0} & \theta_{z_0} &
 x_1 & y_1 & z_1 & \theta_{x_1} & \theta_{y_1} & \theta_{z_1}
\end{bmatrix}^T
$$

In code notation, the order is:

```text
x0, y0, z0, tx0, ty0, tz0, x1, y1, z1, tx1, ty1, tz1
```

This is the reference ordering for all $12 \times 12$ shaft element matrices in Spinniped.

---

## 4. Local element index map

Using zero-based Python indexing, the local element degrees of freedom are:

| Local index | DOF |
|---:|---|
| 0 | $x_0$ |
| 1 | $y_0$ |
| 2 | $z_0$ |
| 3 | $\theta_{x_0}$ |
| 4 | $\theta_{y_0}$ |
| 5 | $\theta_{z_0}$ |
| 6 | $x_1$ |
| 7 | $y_1$ |
| 8 | $z_1$ |
| 9 | $\theta_{x_1}$ |
| 10 | $\theta_{y_1}$ |
| 11 | $\theta_{z_1}$ |


Any matrix written in a different ordering must be permuted before being used in Spinniped.

For example, a matrix written in the common plane-separated ordering

$$
\begin{bmatrix}
 x_0 & \theta_{y_0} & x_1 & \theta_{y_1} &
 y_0 & \theta_{x_0} & y_1 & \theta_{x_1}
\end{bmatrix}^T
$$

is not directly compatible with Spinniped's element convention. This is quite important: many matrices presented in textbooks and papers use this kind of dof ordering. Do not panick: you'll just have to permute the terms of the original matrix to match Spinniped's convention.

---

## 5. Global degree-of-freedom numbering

If the internal compact node index is $p$, the corresponding global degrees of freedom are:

$$
\begin{aligned}
x_{p} &= 6 p \\
y_{p} &= 6 p + 1 \\
z_{p} &= 6 p + 2 \\
\theta_{x_p} &= 6 p + 3 \\
\theta_{y_p} &= 6 p + 4 \\
\theta_{z_p} &= 6 p + 5 \\
\end{aligned}
$$

Therefore, for an element connecting compact node indices $p$ and $q$, the global assembly index vector is:

$$\mathbf{i}_e =
\begin{bmatrix}
6p & 6p+1 & 6p+2 & 6p+3 & 6p+6 & 6p+5 &
6q & 6q+1 & 6q+2 & 6q+3 & 6q+4 & 6q+5
\end{bmatrix}$$

Node IDs stored in the mesh do not have to be identical to compact node indices. If node IDs are arbitrary or non-contiguous, the assembler must use a node-to-index map.

---

## 6. Matrix notation

Element matrices are denoted with the subscript $e$:

$$
\mathbf{M}_e, \quad
\mathbf{C}_e, \quad
\mathbf{G}_e, \quad
\mathbf{K}_e
$$

where:

| Symbol | Meaning |
|---|---|
| $\mathbf{M}_e$ | element mass matrix |
| $\mathbf{C}_e$ | element damping matrix |
| $\mathbf{G}_e$ | element gyroscopic matrix |
| $\mathbf{K}_e$ | element stiffness matrix |

Global matrices are written without the element subscript:

$$
\mathbf{M}, \quad
\mathbf{C}, \quad
\mathbf{G}, \quad
\mathbf{K}
$$

The standard second-order rotor dynamic equation is written as:

$$\mathbf{M}\ddot{\mathbf{q}} + \left( \mathbf{C} + \Omega \mathbf{G} \right)\dot{\mathbf{q}} + \mathbf{K}\mathbf{q} = \mathbf{f}$$

where:

- $\mathbf{q}$ is the global displacement vector;
- $\dot{\mathbf{q}}$ is the global velocity vector;
- $\ddot{\mathbf{q}}$ is the global acceleration vector;
- $\mathbf{f}$ is the external force vector.

Depending on the formulation, the factor $\Omega$ may already be included inside the gyroscopic matrix. Each element document must state this explicitly.

---

## 7. Common scalar symbols

| Symbol | Meaning | Typical unit |
|---|---|---:|
| $L$ | element length | m |
| $A$ | cross-sectional area | m$^2$ |
| $E$ | Young's modulus | Pa |
| $G$ | shear modulus | Pa |
| $\rho$ | mass density | kg/m$^3$ |
| $I_x$ | second moment of area about the $x$ axis | m$^4$ |
| $I_y$ | second moment of area about the $y$ axis | m$^4$ |
| $I_z$ | second moment of area about the $z$ axis | m$^4$ |
| $I_p$ | polar second moment of area | m$^4$ |
| $J$ | torsional constant or polar inertia symbol, depending on context | m$^4$ or kg m$^2$ |
| $\Omega$ | rotor spin speed | rad/s |
| $\omega$ | whirl or natural circular frequency | rad/s |
| $f$ | frequency | Hz |

When $J$ is used, the document must explicitly state whether it represents a geometric torsional constant or a mass polar moment of inertia.
In many axisymmetric element types, the second moments of area are simply called $I_x = I_y = I$ and $I_z = J$. 

---

## 8. Disc and bearing notation

A concentrated disc element located at node $i$ uses the same nodal vector:

$$
\mathbf{q}_i =
\begin{bmatrix}
 x_i & y_i & z_i & \theta_{x_i} & \theta_{y_i} & \theta_{z_i}
\end{bmatrix}^T
$$

Typical disc parameters are:

| Symbol | Meaning |
|---|---|
| $m_d$ | disc mass |
| $I_{d,x}$ | disc diametral mass moment of inertia about $x$ |
| $I_{d,y}$ | disc diametral mass moment of inertia about $y$ |
| $I_{d,z}$ or $I_{d,p}$| disc polar mass moment of inertia about the spin axis, usually $z$ |

A bearing or support element located at node $i$ can be represented by a nodal stiffness and damping relation:

$$\mathbf{f}_b = - \mathbf{K}_b \mathbf{q}_i - \mathbf{C}_b \dot{\mathbf{q}}_i $$

For a simple isotropic bearing acting only on translations:

$$
\mathbf{K}_b =
\begin{bmatrix}
 k_{xx} & 0         & 0         & 0 & 0 & 0 \\
 0      & k_{yy}    & 0         & 0 & 0 & 0 \\
 0      & 0         & k_{zz}    & 0 & 0 & 0 \\
 0      & 0         & 0         & 0 & 0 & 0 \\
 0      & 0         & 0         & 0 & 0 & 0 \\
 0      & 0         & 0         & 0 & 0 & 0 \\
\end{bmatrix}
$$

with $k_{xx} = k_{yy}$ for an isotropic support. 

---

## 9. Sign and symmetry conventions

Unless stated otherwise:

- displacements are positive along the positive global axes;
- rotations are positive according to the right-hand rule;
- stiffness and mass matrices are symmetric for conservative shaft elements;
- gyroscopic matrices are generally skew-symmetric;
- damping matrices may be symmetric or non-symmetric depending on the bearing model.

For a gyroscopic matrix:

$$
\mathbf{G}^T = -\mathbf{G}
$$

for the ideal skew-symmetric case.

---

## 10. Units

Spinniped theory documents assume coherent SI units:

| Quantity | Unit |
|---|---:|
| length | m |
| mass | kg |
| force | N |
| stiffness | N/m |
| rotational stiffness | N m/rad |
| damping | N s/m |
| angular speed | rad/s |
| frequency | Hz |

Angles are treated as dimensionless in matrix equations, but their physical meaning is radians.

---

## 11. Implementation rule

Every element matrix in the theory notes must declare its degree-of-freedom order.

For Spinniped shaft elements, the required order is always:

$$
\mathbf{q}_e =
\begin{bmatrix}
 x_0 & y_0 & z_0 & \theta_{x_0} & \theta_{y_0} & \theta_{z_0} &
 x_1 & y_1 & z_1 & \theta_{x_1} & \theta_{y_1} & \theta_{z_1}
\end{bmatrix}^T
$$

In code notation:

```text
x0, y0, z0, tx0, ty0, tz0, x1, y1, z1, tx1, ty1, tz1
```

This convention takes precedence over external textbook or paper conventions.
