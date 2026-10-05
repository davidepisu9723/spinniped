# Chapter 13 -- Load records

## Table of contents

- [Record types](#record-types)
- [Harmonic convention](#harmonic-convention)
- [Rotating forces](#rotating-forces)
- [Unbalance](#unbalance)
- [Builder behavior](#builder-behavior)

## Record types

The `Load` union contains three frozen, slotted records:

| Record | Amplitude data | Intended excitation |
|---|---|---|
| `HarmonicLoad` | signed local force and moment vectors | fixed-direction harmonic load |
| `RotatingLoad` | nonnegative force magnitude and rotation sense | constant-magnitude circular force |
| `Unbalance` | nonnegative mass-eccentricity product | synchronous circular force proportional to speed squared |

Every record has a positive integer `id`, references one `grid`, and selects a
`coordinate_system`. Load IDs share one namespace across the three record
types. Numeric amplitude and phase fields accept the same random-distribution
references as other model parameters.

## Harmonic convention

Harmonic quantities use

$$
\mathbf{p}(t)=\mathop{\mathrm{Re}}\left\{
\widehat{\mathbf{p}}e^{i\omega t}
\right\}.
$$

For a `HarmonicLoad`, the local phasor is

$$
\widehat{\mathbf{p}}_l=e^{i\phi}
\begin{bmatrix}
F_x & F_y & F_z & M_x & M_y & M_z
\end{bmatrix}^{T}.
$$

The force and moment entries are signed peak amplitudes and share the phase
`phase`, expressed in radians. Multiple records at one grid can represent
components having different phases.

## Rotating forces

A positive `RotatingLoad` acts in the selected frame's x-y plane and rotates
from positive x toward positive y about positive z. Its local phasor is

$$
\widehat{\mathbf{p}}_{l,+}=F e^{i\phi}
\begin{bmatrix}1 & -i & 0 & 0 & 0 & 0\end{bmatrix}^{T}.
$$

The negative sense replaces $-i$ with $+i$. Thus, a rotating force is a
particular harmonic load: two equal orthogonal components in phase quadrature.

## Unbalance

`Unbalance.mass_eccentricity` stores $m e$ in mass-length units. At rotor
angular speed $\Omega$, its synchronous force magnitude is

$$
F_u(\Omega)=m e\Omega^2,
\qquad \omega=|\Omega|.
$$

The selected frame defines phase zero through its local positive x-axis. Its
local z-axis must be parallel or antiparallel to `ModelDefinition.spin_axis`,
so its x-y plane is the physical rotor plane. The alignment sign can determine
the local circular sense when force assembly is implemented.

## Builder behavior

`ModelBuilder` validates load IDs, record types, grid and coordinate-system
references, finite vector components and phases, rotation sense, physical
bounds, and the unbalance plane. It also resolves random parameters separately
for every realization. Loads do not modify $\mathbf{K}$, $\mathbf{M}$,
$\mathbf{C}$, or $\mathbf{G}$ and are retained in each
`BuiltModel.definitions` entry.

The current response solvers still accept explicit force vectors. Conversion
of these records into complex global force vectors is intentionally deferred
to the load-enabled response implementation.
