# Chapter 1 -- Language, coordinates, and workflow

## Table of contents

- [1.1 Reference system](#11-reference-system)
- [1.2 Symbols used throughout](#12-symbols-used-throughout)
- [1.3 The analytical workflow](#13-the-analytical-workflow)

## 1.1 Reference system

This guide follows the normative project conventions in
[Notation and conventions](../reference-manual/01_architecture_and_conventions.md):

- the undeformed shaft centerline is the global $z$-axis;
- $x$ and $y$ are transverse directions;
- translations are $x$, $y$, and $z$;
- rotations are $\theta_x$, $\theta_y$, and $\theta_z$;
- positive rotations follow the right-hand rule;
- angular speed $\Omega$ is in rad/s and cyclic frequency $f$ is in Hz;
- $\omega=2\pi f$ denotes a vibration circular frequency.

The distinction between $\Omega$ and $\omega$ matters. A rotor can spin at
$\Omega$ while vibrating at a different frequency $\omega$.

For bending, the project uses

$$
\frac{\partial x}{\partial z}=\theta_y,
\qquad
\frac{\partial y}{\partial z}=-\theta_x.
$$

Consequently, the two bending planes have equal frequencies for an isotropic
circular shaft but different displacement--rotation coupling signs.

## 1.2 Symbols used throughout

| Symbol | Meaning | SI unit |
|---|---|---|
| $z$ | position along the undeformed shaft | m |
| $t$ | time | s |
| $L$ | shaft length | m |
| $d$, $d_o$, $d_i$ | solid, outer, and inner diameters | m |
| $A$ | cross-sectional area | $\mathrm{m^2}$ |
| $I$ | second moment of area about $x$ or $y$ | $\mathrm{m^4}$ |
| $J$ | geometric polar second moment of area | $\mathrm{m^4}$ |
| $E$ | Young's modulus | Pa |
| $G$ | shear modulus | Pa |
| $\nu$ | Poisson's ratio | dimensionless |
| $\rho$ | mass density | $\mathrm{kg/m^3}$ |
| $m$ | concentrated disk mass | kg |
| $I_d$ | disk diametral mass moment of inertia | $\mathrm{kg\,m^2}$ |
| $I_p$ | disk polar mass moment of inertia | $\mathrm{kg\,m^2}$ |
| $k_b$ | stiffness of one transverse bearing | N/m |
| $\Omega$ | shaft spin speed | rad/s |
| $\omega$ | vibration circular frequency | rad/s |
| $f$ | vibration cyclic frequency, $\omega/(2\pi)$ | Hz |
| $h$ | synchronous harmonic index | dimensionless |

$J$ is a *geometric* section property, whereas $I_d$ and $I_p$ are *mass*
moments of inertia. Their similar names should not be allowed to hide their
different dimensions or physical meanings.

For a circular annulus,

$$
A=\frac{\pi(d_o^2-d_i^2)}{4},\qquad
I=\frac{\pi(d_o^4-d_i^4)}{64},\qquad
J=\frac{\pi(d_o^4-d_i^4)}{32}=2I.
$$

For an isotropic material,

$$
G=\frac{E}{2(1+\nu)}.
$$

## 1.3 The analytical workflow

Use the same sequence for every model in this guide:

1. **Idealize the machine.** Decide which masses, elasticities, constraints,
   and damping mechanisms are retained.
2. **Declare coordinates and signs.** A diagram without defined positive
   directions is incomplete.
3. **Draw the free-body diagram.** Isolate the entire body for reactions and a
   differential element for the field equation.
4. **Write equilibrium or dynamic balance.** Apply force and moment balance,
   including inertia when motion is present.
5. **Add constitutive relations.** Relate force to strain, torque to twist, or
   bending moment to curvature.
6. **Apply boundary conditions.** They select the admissible mode shapes and
   eigenvalues.
7. **Convert $\omega$ to $f$.** Use $f=\omega/(2\pi)$ only at the end.
8. **Check dimensions and limiting cases.** For example, increasing stiffness
   should raise frequency and increasing mass should lower it.
