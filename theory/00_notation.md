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

Each shaft node has four mechanical degrees of freedom:

$$
\mathbf{q}_i =
\begin{bmatrix}
 x_i \\
 y_i \\
 \theta_{x_i} \\
 \theta_{y_i}
\end{bmatrix}
$$

where:

| Symbol | Code name | Meaning |
|---|---:|---|
| $x_i$ | `x` | transverse displacement of node $i$ along the global $x$ direction |
| $y_i$ | `y` | transverse displacement of node $i$ along the global $y$ direction |
| $\theta_{x_i}$ | `tx` | small rotation of node $i$ about the global $x$ axis |
| $\theta_{y_i}$ | `ty` | small rotation of node $i$ about the global $y$ axis |

The nodal unknown vector is therefore ordered as:

$$
\mathbf{q}_i =
\begin{bmatrix}
 x_i & y_i & \theta_{x_i} & \theta_{y_i}
\end{bmatrix}^T
$$

---

## 3. Two-node element degree-of-freedom ordering

For a two-node shaft element connecting node 1 and node 2, Spinniped uses the local element vector:

$$
\mathbf{q}_e =
\begin{bmatrix}
 x_1 \\
 y_1 \\
 \theta_{x_1} \\
 \theta_{y_1} \\
 x_2 \\
 y_2 \\
 \theta_{x_2} \\
 \theta_{y_2}
\end{bmatrix}
$$

Equivalently:

$$
\mathbf{q}_e =
\begin{bmatrix}
 x_1 & y_1 & \theta_{x_1} & \theta_{y_1} &
 x_2 & y_2 & \theta_{x_2} & \theta_{y_2}
\end{bmatrix}^T
$$

In code notation, the order is:

```text
x1, y1, tx1, ty1, x2, y2, tx2, ty2
```

This is the reference ordering for all $8 \times 8$ shaft element matrices in Spinniped.

---

## 4. Local element index map

Using zero-based Python indexing, the local element degrees of freedom are:

| Local index | DOF |
|---:|---|
| 0 | $x_1$ |
| 1 | $y_1$ |
| 2 | $\theta_{x_1}$ |
| 3 | $\theta_{y_1}$ |
| 4 | $x_2$ |
| 5 | $y_2$ |
| 6 | $\theta_{x_2}$ |
| 7 | $\theta_{y_2}$ |

Any matrix written in a different ordering must be permuted before being used in Spinniped.

For example, a matrix written in the common plane-separated ordering

$$
\begin{bmatrix}
 x_1 & \theta_{y_1} & x_2 & \theta_{y_2} &
 y_1 & \theta_{x_1} & y_2 & \theta_{x_2}
\end{bmatrix}^T
$$

is not directly compatible with Spinniped's element convention.

---

## 5. Global degree-of-freedom numbering

If the internal compact node index is $a$, the corresponding global degrees of freedom are:

$$
\begin{aligned}
dof(x_a) &= 4a \\
dof(y_a) &= 4a + 1 \\
dof(\theta_{x_a}) &= 4a + 2 \\
dof(\theta_{y_a}) &= 4a + 3
\end{aligned}
$$

Therefore, for an element connecting compact node indices $a$ and $b$, the global assembly index vector is:

$$\mathbf{i}_e =
\begin{bmatrix}
4a & 4a+1 & 4a+2 & 4a+3 &
4b & 4b+1 & 4b+2 & 4b+3
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
| $I_p$ | polar second moment of area | m$^4$ |
| $J$ | torsional constant or polar inertia symbol, depending on context | m$^4$ or kg m$^2$ |
| $\Omega$ | rotor spin speed | rad/s |
| $\omega$ | whirl or natural circular frequency | rad/s |
| $f$ | frequency | Hz |

When $J$ is used, the document must explicitly state whether it represents a geometric torsional constant or a mass polar moment of inertia.

---

## 8. Disc and bearing notation

A concentrated disc element located at node $i$ uses the same nodal vector:

$$
\mathbf{q}_i =
\begin{bmatrix}
 x_i & y_i & \theta_{x_i} & \theta_{y_i}
\end{bmatrix}^T
$$

Typical disc parameters are:

| Symbol | Meaning |
|---|---|
| $m_d$ | disc mass |
| $I_{d,x}$ | disc diametral mass moment of inertia about $x$ |
| $I_{d,y}$ | disc diametral mass moment of inertia about $y$ |
| $I_{d,p}$ | disc polar mass moment of inertia about the spin axis |

A bearing or support element located at node $i$ can be represented by a nodal stiffness and damping relation:

$$\mathbf{f}_b = - \mathbf{K}_b \mathbf{q}_i - \mathbf{C}_b \dot{\mathbf{q}}_i $$

For a simple isotropic bearing acting only on translations:

$$
\mathbf{K}_b =
\begin{bmatrix}
 k_{xx} & 0 & 0 & 0 \\
 0 & k_{yy} & 0 & 0 \\
 0 & 0 & 0 & 0 \\
 0 & 0 & 0 & 0
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
\boxed{
\mathbf{q}_e =
\begin{bmatrix}
 x_1 & y_1 & \theta_{x_1} & \theta_{y_1} &
 x_2 & y_2 & \theta_{x_2} & \theta_{y_2}
\end{bmatrix}^T
}
$$

In code notation:

```text
x1, y1, tx1, ty1, x2, y2, tx2, ty2
```

This convention takes precedence over external textbook or paper conventions.
