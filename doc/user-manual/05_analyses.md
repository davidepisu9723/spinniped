# Chapter 5 -- Analyses

## Table of contents

- [Modal analysis](#modal-analysis)
- [Campbell analysis](#campbell-analysis)
- [Frequency response](#frequency-response)
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
- `track_modes`: whether MAC-based branch tracking was enabled.

With `track_modes=True`, roots after the first speed are assigned to the
previous speed by displacement-vector MAC and phase-aligned, forming continuous
modal branches through crossings. Stochastic Campbell results are then matched
at each speed to sample zero. Set `track_modes=False` to sort each speed and
sample independently by positive imaginary part. The requested `modes` count is
strict at every speed.

## Frequency response

Input frequencies are in Hz. The force may contain either all global DOFs or
only the free DOFs:

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

`result["response"]` is complex and has shape
`(samples, frequencies, free_dofs)`.

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
