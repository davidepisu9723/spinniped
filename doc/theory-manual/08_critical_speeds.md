# Chapter 8 -- Critical speeds and Campbell diagrams

## Table of contents

- [8.1 From modal branches to critical speeds](#81-from-modal-branches-to-critical-speeds)
- [8.2 Harmonic crossings of the extended-Jeffcott branches](#82-harmonic-crossings-of-the-extended-jeffcott-branches)
- [8.3 Existence and physical interpretation](#83-existence-and-physical-interpretation)
- [8.4 Worked extended-Jeffcott result](#84-worked-extended-jeffcott-result)

This chapter uses the extended-Jeffcott modal branches derived in
[Chapter 7](07_jeffcott_rotor.md) and finds their harmonic intersections.

## 8.1 From modal branches to critical speeds

A Campbell diagram plots each natural frequency $\omega_j(\Omega)$ vertically
against the prescribed rotor speed $\Omega$ horizontally. An excitation with
harmonic index $h$ has angular and cyclic frequencies

$$
\omega_h(\Omega)=h\Omega,
\qquad
f_h(\Omega)=\frac{h\Omega}{2\pi}.
$$

The 1x line corresponds to $h=1$. A critical speed for mode $j$ and harmonic
$h$ is a spin speed that satisfies

$$
\omega_j(\Omega_{\mathrm{crit}})
=h\Omega_{\mathrm{crit}}.
$$

Modal analysis prescribes $\Omega$ and calculates $\omega_j$; critical-speed
analysis instead searches for the intersecting value of $\Omega$.

## 8.2 Harmonic crossings of the extended-Jeffcott branches

The cylindrical lateral frequency derived in
[Section 7.4](07_jeffcott_rotor.md#74-cylindrical-lateral-translation-modes) is
constant. Its crossing with harmonic $h$ is therefore

$$
\boxed{
\Omega_{t,\mathrm{crit}}^{(h)}
=\frac{1}{h}\sqrt{\frac{k_t}{m}}
}.
$$

For the conical modes, there is no need to substitute the square-root branch
formulas again. Insert $\omega=h\Omega$ directly into the two quadratic
characteristic equations derived in
[Section 7.5](07_jeffcott_rotor.md#75-conical-modes-and-gyroscopic-coupling).
They reduce to

$$
h(hI_d+I_p)\Omega^2=k_r
\qquad\text{(backward)},
$$

$$
h(hI_d-I_p)\Omega^2=k_r
\qquad\text{(forward)}.
$$

Thus the conical critical speeds are

$$
\boxed{
\Omega_{c,-,\mathrm{crit}}^{(h)}
=\sqrt{\frac{k_r}{h(hI_d+I_p)}}
},
\qquad
\boxed{
\Omega_{c,+,\mathrm{crit}}^{(h)}
=\sqrt{\frac{k_r}{h(hI_d-I_p)}}
}.
$$

## 8.3 Existence and physical interpretation

The cylindrical lateral crossing exists for every $h>0$. The backward-conical
denominator is also positive for every $h>0$, so that crossing always exists.
The forward-conical result is real and finite only when

$$
hI_d>I_p.
$$

If $hI_d=I_p$, the harmonic line approaches the forward branch only
asymptotically. If $hI_d<I_p$, the forward branch remains above it and there is
no crossing. For the rotor used below, $I_p=2I_d$; consequently the forward
branch has no 1x crossing. This is why the worked result reports a cylindrical
lateral and a backward-conical 1x critical speed, but no forward-conical one.

A Campbell intersection identifies a possible resonance, not its amplitude.
The excitation must also have the appropriate spatial distribution and whirl
direction. Ordinary mass imbalance is a forward-rotating 1x excitation, so a
backward-branch intersection is not automatically excited; anisotropy,
coupling, or a backward-rotating force may be required. Damping and the mode
shape then determine the response magnitude.

## 8.4 Worked extended-Jeffcott result

Use

$$
E=2.0\times10^{11}\ \mathrm{Pa},\quad
d=0.02\ \mathrm m,\quad
L=1.0\ \mathrm m,
$$

$$
m=5.0\ \mathrm{kg},\quad
I_d=0.025\ \mathrm{kg\,m^2},\quad
I_p=0.05\ \mathrm{kg\,m^2},
$$

$$
k_b=1.0\times10^6\ \mathrm{N/m}.
$$

Then

| Quantity | Analytical value |
|---|---:|
| $I$ | $7.853981634\times10^{-9}\ \mathrm{m^4}$ |
| $k_t$ | $72659.042323\ \mathrm{N/m}$ |
| $k_r$ | $18164.760581\ \mathrm{N\,m/rad}$ |
| Cylindrical lateral frequency at rest | $19.185802264\ \mathrm{Hz}$ |
| Conical frequency at rest | $135.664108836\ \mathrm{Hz}$ |
| Cylindrical lateral 1x critical speed | $1151.148135859\ \mathrm{rpm}$ |
| Backward-conical 1x critical speed | $4699.542585350\ \mathrm{rpm}$ |
| Forward-conical 1x critical speed | none because $I_d<I_p$ |

The analytical mode shapes and Campbell curves are shown below:

![First stationary-shaft modes resolved into their physical degrees of freedom, beside the extended-Jeffcott Campbell diagram](images/analytical_rotordynamics.png)

*The six left panels separate the normalized stationary-shaft mode components:
axial translation $u_z$, torsional rotation $\theta_z$, transverse translations
$x$ and $y$, and their associated bending rotations. The project sign
convention gives $\theta_y=\partial x/\partial z$ and
$\theta_x=-\partial y/\partial z$. Modal amplitudes and signs are arbitrary;
the curves communicate the spatial distribution and the active DOF.*

The dashed 1x line intersects the cylindrical lateral branch and the
decreasing backward-conical branch at the two tabulated critical speeds. It
does not intersect the increasing forward-conical branch for these inertia
values.

For the direct analytical-to-finite-element comparison, run
[`02_flexible_bearing_jeffcott_comparison.py`](../../benchmark/02_flexible_bearing_jeffcott_comparison.py).
