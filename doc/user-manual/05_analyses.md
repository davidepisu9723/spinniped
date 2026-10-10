# Chapter 5 -- Analyses

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

## Modal analysis

```python
result = solver.solve(
    "modal",
    fixed_dofs=fixed,
    modes=8,
    track_modes=True,
)
```

Result values include:

- `eigenvalues`: `(samples, modes)` values of $\omega^2$ in $\mathrm{s}^{-2}$;
- `frequencies`: `(samples, modes)` frequencies in Hz;
- `eigenvectors`: `(samples, free_dofs, modes)` reduced mode shapes.

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

With `track_modes=True`, roots after the first speed are assigned to the
previous speed by displacement-vector MAC and phase-aligned, forming continuous
modal branches through crossings. Stochastic Campbell results are then matched
at each speed to sample zero. Set `track_modes=False` to sort each speed and
sample independently by positive imaginary part. The requested `modes` count is
strict at every speed.

Tracking physically permutes the result's mode axis: array position `j` holds
the candidate assigned to reference position `j`, rather than necessarily the
$j$th frequency at that speed. The
[solver reference](../reference-manual/08_solvers.md#mode-correspondence)
documents the MAC matrix, assignment, array permutation, and phase alignment.

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

Input frequencies are finite nonnegative values in Hz. Rotor `speed` is one
signed operating speed in rad/s; it scales the gyroscopic matrix but does not
set the independently swept forcing frequencies.

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

`result["response"]` is complex and has shape
`(samples, frequencies, free_dofs)`. The result also returns scalar `speed`
and the selected `load_ids` as a tuple; the latter is empty for the raw-vector
route.

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
