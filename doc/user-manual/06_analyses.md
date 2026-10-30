# Chapter 6 -- Analyses

## Table of contents

- [Modal analysis](#modal-analysis)
- [Campbell analysis](#campbell-analysis)
- [Frequency response](#frequency-response)
- [Unbalance response](#unbalance-response)
- [Time response](#time-response)

`Solver.solve` is the single analysis entry point. Supply the built model to the
constructor:

```python
from spinniped import Solver

fixed = [2, 5]  # Remove axial and torsional rigid motions.
solver = Solver(model)
result = solver.solve("modal", fixed_dofs=fixed, modes=8)
```

or to an individual call:

```python
solver = Solver()
result = solver.solve(
    "modal",
    model=model,
    fixed_dofs=fixed,
    modes=8,
)
```

`fixed_dofs` contains global DOF indices. The solver removes duplicate indices,
validates their range, and solves on the remaining DOFs. Every result is a
dictionary containing `analysis`, `free_dofs`, `samples`, and `stochastic` in
addition to analysis-specific values. Solver inputs must be finite and have the
documented shapes; unknown options are rejected. Results are regular NumPy
arrays, never ragged object arrays. If an analysis cannot return a consistent
mode count, it raises an error and asks for fewer or a fixed number of modes.

Every `speed` and `speeds` input is signed rotor angular speed in rad/s. When
the model contains `SpeedDependentValue` records, the solver interpolates all
required structural operating points before the first eigensolution, linear
solve, or time integration. It rejects the complete call if any point falls
outside a table's range.

## Modal analysis

```python
result = solver.solve(
    "modal",
    fixed_dofs=fixed,
    modes=8,
    track_modes=True,
    speed=100.0,
)
```

Result values include:

- `eigenvalues`: `(samples, modes)` values of $\omega^2$ in $\mathrm{s}^{-2}$;
- `frequencies`: `(samples, modes)` frequencies in Hz;
- `eigenvectors`: `(samples, free_dofs, modes)` reduced mode shapes;
- `speed`: the scalar rotor speed used to evaluate stiffness.

Eigenvalues within `zero_tolerance` of zero are treated as rigid modes and
removed. The default tolerance is scaled from the eigenspectrum; it may be
overridden with a finite, nonnegative `zero_tolerance`. An eigenvalue below the
negative tolerance raises an instability or insufficient-constraint error
rather than being silently discarded.

The real modal route requires symmetric stiffness and mass matrices; use the
Campbell route when cross-coupled bearing coefficients make stiffness
non-symmetric. Damping and gyroscopic terms are not included in modal analysis.
The requested `modes` count is strict: requesting more positive modes than are
available raises an error.

`speed` defaults to `0.0`. It selects $\mathbf{K}(\Omega)$ when bearing
stiffness is speed-dependent, but it does not introduce damping or gyroscopic
terms into this real modal route. The selected speed must lie within every
stiffness table used by the model.

For stochastic models, `track_modes=True` (the default) matches every later
sample to sample zero by maximizing the Modal Assurance Criterion (MAC), then
aligns eigenvector signs. Set `track_modes=False` to retain independent
ascending-eigenvalue order in every sample.

## Campbell analysis

```python
result = solver.solve(
    "campbell",
    fixed_dofs=fixed,
    speeds=[0.0, 100.0, 200.0],  # rad/s
    modes=8,
    track_modes=True,
    harmonics=[1.0, 2.0],
)
```

Its analysis-specific results are:

- `speeds`: the requested angular speeds in rad/s, shape `(speeds,)`;
- `speeds_hz`: those speeds divided by $2\pi$, shape `(speeds,)`;
- `eigenvalues`: complex state eigenvalues, shape
  `(samples, speeds, modes)`;
- `frequencies`: absolute imaginary parts in Hz, shape
  `(samples, speeds, modes)`;
- `eigenvectors`: complex state vectors, shape
  `(samples, speeds, 2*free_dofs, modes)`;
- `track_modes`: whether MAC-based branch tracking was enabled;
- `harmonics`: requested positive harmonic indices $h$, shape `(harmonics,)`;
- `critical_speeds`: interpolated angular speeds in rad/s, shape
  `(samples, harmonics, modes, crossings)`;
- `critical_speeds_hz`: the same values divided by $2\pi$;
- `critical_speed_counts`: number of crossings found for each sample,
  harmonic, and mode, shape `(samples, harmonics, modes)`.

Before entering the speed loop, the solver evaluates every bearing stiffness
and damping table over the complete `speeds` input. At each loop position it
uses $\mathbf{K}(\Omega_i)$ and $\mathbf{C}(\Omega_i)$ together with the fixed
mass and unit-speed gyroscopic matrices.

With `track_modes=True`, roots after the first speed are assigned to the
previous speed by displacement-vector MAC and phase-aligned, forming continuous
modal branches through crossings. Stochastic Campbell results are then matched
at each speed to sample zero. Set `track_modes=False` to sort each speed and
sample independently by positive imaginary part. The requested `modes` count is
strict at every speed.

Tracking physically permutes the result's mode axis: array position `j` holds
the candidate assigned to reference position `j`, rather than necessarily the
$j$th frequency at that speed. The
[solver reference](../reference-manual/10_solvers.md#mode-correspondence)
documents the MAC matrix, assignment, array permutation, and phase alignment.

### Starting a symmetric stochastic sweep near zero

An isotropic rotor can have exactly repeated forward- and backward-whirl
eigenvalues at zero speed. The repeated eigenspace is physically meaningful,
but the individual eigenvectors returned by an eigensolver are not unique:
any orthogonal basis of that subspace is valid. Consequently, pairwise MAC
cannot always give those vectors consistent forward/backward labels across
stochastic realizations when zero speed is the first Campbell point.

The resulting Campbell plot can look much more uncertain than the underlying
model. Different realizations carrying the same mode index may follow opposite
whirl branches after they separate. A mean or confidence band then combines
two physically distinct branches, producing a broad fan, an intermediate mean
that is not itself a branch, or an apparent discontinuity. This is a mode-label
ambiguity, not stiffness uncertainty created by speed interpolation.

<p align="center">
  <img src="images/stochastic_campbell_zero_speed_tracking.png" width="80%" alt="The same stochastic Campbell calculation started at zero and one rpm; the zero-speed case assigns one mode label to opposite whirl branches, while the nonzero start keeps the branch correspondence narrow">
</p>

*Each thin line is the tracked fourth mode of one realization of the symmetric
speed-dependent rotor from Example 09. With a first point of exactly zero rpm,
the repeated initial eigenspace permits the same label to leave along either
whirl branch. Starting at 1 rpm breaks that exact degeneracy before labels are
established. The bearing ordinates in this comparison are perfectly correlated,
so the split cannot be attributed to independent random interpolation knots.*

For a one-direction positive-speed sweep of such a symmetric model, begin at a
small positive speed, for example `1.0` rpm, and retain zero as a knot in every
`SpeedDependentValue` table so interpolation remains valid. One rpm is not a
universal tolerance; it is simply small relative to this example's 6,000 rpm
range while large enough to distinguish the numerical whirl vectors. Verify
that reducing the starting speed does not materially change the reported
frequencies or critical speeds. Keep zero when the stationary solution itself
is required, when the modes are not repeated, or when a future subspace-based
tracking method is used.

When `track_modes=True` and the first requested speed is exactly zero, the
solver emits a `UserWarning` describing this ambiguity and suggesting a
small nonzero starting speed. The warning is advisory: the solve continues,
because zero remains valid and may be intentionally required.

When `harmonics` is supplied, critical speeds are intersections between each
tracked modal branch and

$$
f_h(\Omega)=\frac{h\Omega}{2\pi},
\qquad h>0.
$$

The solver detects exact intersections at sampled speeds and sign changes
between adjacent sampled speeds. A sign-change intersection is found by linear
interpolation of the two existing residual values. It never evaluates an
additional speed or refines the user's speed sequence. Harmonic detection
consequently requires at least two strictly increasing speeds and
`track_modes=True`.

Most branches cross a given harmonic once, but the last result axis preserves
multiple crossings in ascending speed order. Its length is the largest number
found anywhere in the result. Missing crossings are padded with `NaN`; use
`critical_speed_counts` to distinguish them explicitly. For stochastic
results, the sample axis retains corresponding tracked modes and crossing
order, allowing finite values to be summarized across realizations.

## Frequency response

For a raw force, `HarmonicLoad`, or `RotatingLoad`, input frequencies are finite
nonnegative values in Hz. Rotor `speed` is one signed operating speed in rad/s;
it scales the gyroscopic matrix but does not set those independently swept
forcing frequencies. It also selects any speed-dependent bearing coefficients
and declarative load amplitudes once for the whole forcing-frequency sweep.
`SynchronousLoad`, described below, is the order-locked exception.

The legacy `force` input may contain either all global DOFs or only the free
DOFs:

```python
force = [0.0] * model.ndof
force[6] = 1.0

result = solver.solve(
    "frequency_response",
    fixed_dofs=fixed,
    frequencies=[10.0, 20.0, 30.0],
    force=force,
    speed=100.0,  # rad/s
)
```

Alternatively, select one or more `HarmonicLoad` or `RotatingLoad` records
already stored in the model definition:

```python
result = solver.solve(
    "frequency_response",
    fixed_dofs=fixed,
    frequencies=[10.0, 20.0, 30.0],
    load_ids=[1, 2],
    speed=100.0,
)
```

Supply exactly one of `force` and `load_ids`. Selected records are transformed
from their own coordinate systems, added in global DOF order, and then reduced
by `fixed_dofs`. Stochastic load parameters are assembled from each resolved
realization rather than reused from sample zero. An `Unbalance` ID is rejected
here because its amplitude and forcing frequency depend synchronously on rotor
speed.

When a selected `HarmonicLoad` component or `RotatingLoad.force_magnitude` is
speed-dependent, it is interpolated at the scalar rotor `speed`, not at each
entry of `frequencies`. All selected declarative loads are evaluated before the
first dynamic-stiffness solve.

`result["response"]` is complex and has shape
`(samples, frequencies, free_dofs)`. The result also returns scalar `speed`
and the selected `load_ids` as a tuple; the latter is empty for the raw-vector
route.

### Order-locked synchronous loads

`SynchronousLoad` inherits the force, moment, phase, and coordinate-system
fields of `HarmonicLoad` and adds a required positive real `harmonic_order`.
For each scalar signed rotor speed $\Omega$, `frequency_response` derives

$$
\omega=\mathtt{harmonic\_order}|\Omega|,
\qquad
f=\frac{\omega}{2\pi}.
$$

Define and sweep a 2x load without supplying `frequencies`:

```python
from spinniped import SynchronousLoad

load = SynchronousLoad(
    id=4,
    grid=2,
    harmonic_order=2.0,
    force=(10.0, 0.0, 0.0),
)

result = solver.solve(
    "frequency_response",
    fixed_dofs=fixed,
    speeds=[0.0, 100.0, 200.0],
    load_ids=[4],
)
```

The returned `frequencies` array contains
$2|\mathtt{speeds}|/(2\pi)$ Hz. The result also contains signed `speeds`,
signed `speeds_hz`, and `harmonic_order=2.0`; its response shape is
`(samples, speeds, free_dofs)`. Scalar `speed` remains available for a
single-point calculation. An explicit `frequencies` argument is rejected.
Multiple selected synchronous loads must share one order and cannot be mixed
with ordinary harmonic or rotating loads. Use separate calls for 1x, 2x, 3x,
or other orders.

If a synchronous force or moment component is a `SpeedDependentValue`, its
amplitude is interpolated over the complete signed rotor `speeds` grid; the
derived forcing frequencies are not used as table coordinates. Structural and
load tables are preflighted before the first response point is solved.

## Unbalance response

Unbalance response sweeps signed rotor speeds in rad/s. At each point, the
forcing frequency and magnitude are

$$
\omega=|\Omega|,
\qquad
F_u=m e\Omega^2.
$$

```python
result = solver.solve(
    "unbalance_response",
    fixed_dofs=fixed,
    speeds=[0.0, 100.0, 200.0],
    load_ids=[3],
)
```

`load_ids` must be a nonempty sequence selecting only `Unbalance` records.
Multiple records are superimposed. Positive and negative speeds reverse the
physical rotation sense and enter the gyroscopic term with their supplied
sign. The load frame may have its local z-axis parallel or antiparallel to the
model spin axis; the solver accounts for that alignment when constructing the
rotating force.

The complete speed grid is also used to precompute speed-dependent bearing
stiffness and damping. `Unbalance.mass_eccentricity` itself does not accept a
`SpeedDependentValue`; its $\Omega^2$ amplitude law is evaluated directly.

The result contains signed `speeds`, signed `speeds_hz`, nonnegative
synchronous `frequencies`, the selected `load_ids` tuple, and complex
`response` with shape `(samples, speeds, free_dofs)`.

## Time response

Times are in seconds. `force` may be a constant vector or a callable that
accepts time and returns a full or reduced force vector:

```python
result = solver.solve(
    "time_response",
    fixed_dofs=fixed,
    times=[0.0, 0.001, 0.002],
    force=force,
    speed=100.0,
)
```

An optional `initial_state` contains the free-DOF displacements followed by
the free-DOF velocities. The returned `displacement` and `velocity` arrays have
shape `(samples, times, free_dofs)`.

The scalar `speed` selects one structural operating point for the complete
integration. Time response still accepts an explicit constant vector or
callback rather than declarative load IDs, so speed-dependent declarative load
records are not assembled by this route.
