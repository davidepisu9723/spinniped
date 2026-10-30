# Chapter 6 -- Load records

## Table of contents

- [Record types](#record-types)
- [Harmonic convention](#harmonic-convention)
- [Rotating forces](#rotating-forces)
- [Unbalance](#unbalance)
- [Builder behavior](#builder-behavior)
- [Solver assembly](#solver-assembly)

## Record types

The `Load` union contains three frozen, slotted records:

| Record | Amplitude data | Intended excitation |
|---|---|---|
| `HarmonicLoad` | signed local force and moment vectors | fixed-direction harmonic load |
| `RotatingLoad` | nonnegative force magnitude and rotation sense | circular force |
| `Unbalance` | nonnegative mass-eccentricity product | synchronous circular force proportional to speed squared |

Every record has a positive integer `id`, references one `grid`, and selects a
`coordinate_system`. Load IDs share one namespace across the three record
types. Numeric amplitude and phase fields accept the same random-distribution
references as other model parameters.

Each `HarmonicLoad.force` and `HarmonicLoad.moment` component, and
`RotatingLoad.force_magnitude`, may instead be a `SpeedDependentValue`.
`phase` and every `Unbalance` field remain speed-independent. Tabulated load
amplitudes use signed rotor speed $\Omega$, not forcing frequency $\omega$.
The outer `force` and `moment` vectors must be three-component tuples; mutable
lists are rejected at record construction.

`HarmonicLoad` construction validates its positive identity and grid reference,
nonnegative frame reference, phase parameter, vector shapes, and every
constant, random, or tabulated component. The builder retains grid/frame
existence and load-ID uniqueness checks.

`RotatingLoad` construction likewise validates its identity and reference
syntax, speed-independent phase, rotation sense, and nonnegative scalar or
tabulated magnitude. Random magnitude references are revalidated when model
resolution reconstructs the record. The builder retains only grid/frame
existence and shared load-ID uniqueness checks for this record.

`Unbalance` construction validates its identity and reference syntax, phase,
and nonnegative mass-eccentricity parameter. Random mass-eccentricity
references are revalidated after resolution. The builder retains grid/frame
existence, shared load-ID uniqueness, and the contextual requirement that the
selected frame's z-axis be parallel to the model spin axis.

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

A `RotatingLoad` acts in the selected frame's x-y plane. Let $s=+1$ denote
positive rotation from positive x toward positive y about positive z, and let
$s=-1$ denote the reverse sense. To preserve $\phi$ as the signed geometric
angle at $t=0$, its local phasor is

$$
\widehat{\mathbf{p}}_{l,s}=F e^{is\phi}
\begin{bmatrix}1 & -is & 0 & 0 & 0 & 0\end{bmatrix}^{T}.
$$

Thus, a rotating force is a particular harmonic load: two equal orthogonal
components in phase quadrature. At every frequency-response point it rotates
at the requested forcing frequency, independently of rotor speed.

## Unbalance

`Unbalance.mass_eccentricity` stores $m e$ in mass-length units. At rotor
angular speed $\Omega$, its synchronous force magnitude is

$$
F_u(\Omega)=m e\Omega^2,
\qquad \omega=|\Omega|.
$$

The selected frame defines phase zero through its local positive x-axis. Its
local z-axis must be parallel or antiparallel to `ModelDefinition.spin_axis`,
so its x-y plane is the physical rotor plane. For load-frame rotation
$\mathbf{R}$ and normalized spin direction $\widehat{\mathbf{s}}$, define

$$
a=\mathop{\mathrm{sign}}\left(
\mathbf{R}_{:,z}^{T}\widehat{\mathbf{s}}
\right),
\qquad
s=\mathop{\mathrm{sign}}(\Omega)a.
$$

At nonzero speed, unbalance uses the rotating-load phasor above with
$F=F_u(\Omega)$. At zero speed its force is exactly zero. This construction
handles reverse rotor speed and an antiparallel load-frame z-axis without
reversing the intended global rotation accidentally.

## Builder behavior

`ModelBuilder` validates collection-wide load-ID uniqueness, record types,
grid and coordinate-system existence, and the unbalance plane. It also
resolves random parameters separately for every realization; reconstruction
of each resolved record repeats its local validation. Loads do not modify
$\mathbf{K}$, $\mathbf{M}$, $\mathbf{C}$, or $\mathbf{G}$ and are retained in
each `BuiltModel.definitions` entry.

For frequency response, all selected speed-dependent load components are
interpolated at the call's scalar rotor `speed`, transformed, and summed before
the first forcing-frequency solve. Out-of-range interpolation fails before a
partial frequency sweep can be returned.

## Solver assembly

For a load-frame local-to-global rotation $\mathbf{R}$, the solver transforms
one local nodal phasor according to

$$
\widehat{\mathbf{p}}_g=
\begin{bmatrix}
\mathbf{R} & \mathbf{0}\\
\mathbf{0} & \mathbf{R}
\end{bmatrix}
\widehat{\mathbf{p}}_l.
$$

It adds all selected records at their grids in full global DOF order and only
then removes fixed DOFs. No moment-arm term is required because each record is
applied directly at its referenced grid. Assembly is repeated for every
stochastic realization so resolved amplitudes, phases, and frames remain
sample-specific.

`frequency_response` accepts either one legacy complex force vector or a
nonempty selection of `HarmonicLoad` and `RotatingLoad` IDs. It rejects
`Unbalance` records. `unbalance_response` accepts a nonempty selection of
`Unbalance` IDs and evaluates the synchronous force separately at every
supplied rotor speed. Duplicate, unknown, noninteger, and analysis-incompatible
IDs are rejected rather than ignored.
