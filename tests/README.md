# Testing and Verification

Spinniped's tests serve two purposes: preventing software regressions and
checking that the finite-element equations reproduce expected mechanical
behavior. The suite exercises the declarative records, stateless local-matrix
kernels, deterministic and stochastic assembly, every solver route, and an
independent beam benchmark.

## Running the suite

From the repository root, install Spinniped and its test dependency:

```bash
python -m pip install -e ".[test]"
```

Run all tests:

```bash
python -m pytest
```

Useful focused commands include:

```bash
# Compact output
python -m pytest -q

# One test module
python -m pytest tests/test_declarative_api.py

# One test
python -m pytest tests/test_declarative_api.py::test_modal_solver_returns_real_modes_for_each_sample

# Stop at the first failure
python -m pytest -x

# Show test names and assertions verbosely
python -m pytest -vv
```

## Test summary

The table is both a table of contents and a quick map from a failure to the
part of the implementation it exercises. Counts refer to test functions;
parametrized functions run once for every listed case.

| Test file | Tests | Scope |
|---|---:|---|
| [`test_element_matrices.py`](#test_element_matricespy) | 10 | Stateless shaft, bearing, and disk matrix kernels |
| [`test_assembly.py`](#test_assemblypy) | 16 | Record validation, coordinate transforms, and global assembly |
| [`test_declarative_api.py`](#test_declarative_apipy) | 19 | Distribution resolution and every public solver route |
| [`test_plotting.py`](#test_plottingpy) | 3 | Campbell-diagram rendering and plotting-input validation |
| [`test_analytical_benchmarks.py`](#test_analytical_benchmarkspy) | 1 | Comparison with the Euler--Bernoulli beam solution |
| [`test_convergence.py`](#test_convergencepy) | 2 | Accuracy under mesh refinement |

The suite also relies on [`conftest.py`](#conftestpy-shared-test-support),
which supplies common properties, fixtures, model factories, support DOFs,
and analytical results.

## Test details

### `conftest.py`: shared test support

The shared fixtures describe a one-metre uniform steel shaft using values
that are independent of the production records:

| Quantity | Symbol | Value |
|---|---:|---:|
| Young's modulus | $E$ | $2.0\times10^{11}\ \mathrm{Pa}$ |
| Poisson's ratio | $\nu$ | $0.3$ |
| Density | $\rho$ | $7850\ \mathrm{kg/m^3}$ |
| Diameter | $d$ | $0.01\ \mathrm{m}$ |
| Total length | $L$ | $1.0\ \mathrm{m}$ |

For the solid circular section,

$$
A=\frac{\pi d^2}{4},
\qquad
I=\frac{\pi d^4}{64},
\qquad
J=2I.
$$

`shaft_properties` exposes these inputs and derived section properties;
`shaft_definition` creates one Timoshenko shaft; and `shaft_local_matrices`
evaluates its stiffness, mass, and gyroscopic kernels. `build_uniform_beam`
creates an arbitrary Euler or Timoshenko mesh through the public declarative
API. The support, analytical-frequency, and bending-pair helpers are shared by
the benchmark and convergence tests. In particular, the analytical helper is
independent of Spinniped's matrix kernels, preventing the same implementation
error from appearing on both sides of a comparison.

### `test_element_matrices.py`

This module tests stateless numerical kernels before coordinate transforms or
assembly can hide a formulation error.

#### `test_shaft_matrices_have_expected_structure`

Evaluates the reference shaft stiffness, mass, and gyroscopic matrices. Every
matrix must be finite and $12\times12$, stiffness and mass must be symmetric,
and the gyroscopic matrix must be skew-symmetric. This establishes the basic
size and reciprocity invariants expected from a two-node, six-DOF-per-node
shaft element.

#### `test_shaft_stiffness_has_six_rigid_body_modes`

Computes all eigenvalues of an unconstrained shaft stiffness matrix and counts
values smaller than $10^{-10}$ times the largest eigenvalue. Exactly six must
be zero, corresponding to three rigid translations and three rigid rotations;
the lower-bound check also rejects materially negative stiffness modes.

#### `test_shaft_mass_is_positive_definite`

Checks that the smallest eigenvalue of the consistent shaft mass matrix is
strictly positive. A failure means that a nonzero element motion could have
zero or negative kinetic energy.

#### `test_axial_and_torsional_blocks_match_closed_form`

Extracts the axial DOFs `[2, 8]` and torsional DOFs `[5, 11]` from the shaft
stiffness matrix. They are compared with $EA/L[[1,-1],[-1,1]]$ and
$GJ/L[[1,-1],[-1,1]]`, respectively, providing a direct check of magnitudes,
signs, and DOF placement against closed-form mechanics.

#### `test_density_scales_mass_and_gyroscopic_but_not_stiffness`

Re-evaluates the kernels after doubling density. Mass and gyroscopic terms
must double exactly while stiffness remains unchanged, demonstrating that
material inertia is linear in density and that density has not leaked into
the elastic formulation.

#### `test_full_shaft_gyroscopic_matrix_couples_translations_and_rotations`

For a Timoshenko shaft with rotary inertia, verifies nonzero coupling between
lateral translations and cross-section rotations at both nodes, as well as
global skew-symmetry. Repeating the calculation with rotary inertia disabled
must produce a zero gyroscopic matrix, checking both branches of the option.

#### `test_bearing_kernels_place_coefficients_in_translational_block`

Uses distinct direct and cross-coupled coefficients so swapped indices or an
accidental transpose are visible. The test checks their exact placement in
the leading $3\times3$ translational stiffness and damping blocks and confirms
that the idealized bearing contributes no mass or gyroscopic terms.

#### `test_disk_kernels_represent_lumped_mass_and_polar_gyroscopic_coupling`

Checks the disk mass diagonal against three translational masses, two
diametral inertias, and one polar inertia. It also verifies the equal and
opposite polar gyroscopic entries at `[3, 4]` and `[4, 3]`, skew-symmetry, and
the absence of disk stiffness.

#### `test_element_damping_is_mass_proportional`

Passes identity mass matrices to the shaft and disk damping kernels. The
outputs must equal the requested damping coefficient times the corresponding
mass matrix, isolating the documented mass-proportional damping rule.

#### `test_invalid_shaft_geometry_is_rejected`

Runs three nonphysical cases: zero length, equal inner and outer diameters,
and a negative inner diameter. Each must raise a `ValueError` mentioning the
positive-geometry requirement, ensuring invalid sections fail before any
matrix is formed.

### `test_assembly.py`

This module builds deliberately small `ModelDefinition` objects and inspects
the resulting global matrices and resolved coordinates.

#### `test_one_element_global_matrices_equal_local_kernels`

Builds one shaft aligned with global $z$. Because no rotation or multi-element
accumulation is required, its global stiffness and mass matrices must equal
the independently evaluated local kernels, including the leading sample axis
and the grid-ID order.

#### `test_two_elements_accumulate_at_shared_grid`

Joins two half-length shafts at a middle grid. The middle $6\times6$ block
must be the sum of the adjacent local end blocks, while the two unconnected
outer grids must have no direct coupling; together these assertions check
both accumulation and connectivity.

#### `test_nonconsecutive_grid_ids_are_mapped_by_record_order`

Uses grid IDs 10 and 70 to prove that IDs are labels rather than matrix
indices. The builder must retain record order, allocate contiguous DOF blocks,
and reproduce the same local stiffness and mass matrices.

#### `test_shaft_gyroscopic_assembly_is_invariant_to_endpoint_order`

Builds the same spinning shaft once in each endpoint order. Its gyroscopic
matrix must be nonzero and identical after global assembly, guarding against
an element-orientation sign error.

#### `test_bearing_terms_are_inserted_at_referenced_grid`

Places a grounded bearing at the middle of three grids. Distinct direct and
cross-coupled values must appear only in that grid's translational stiffness
and damping block, with the expected five nonzero entries and no mass or
gyroscopic contribution.

#### `test_two_grid_bearing_assembles_equal_and_opposite_blocks`

Connects a bearing between two grids and compares the complete global matrix
with $[[K,-K],[-K,K]]$. This sign pattern proves that the element responds to
relative displacement and applies equal-and-opposite nodal forces.

#### `test_two_grid_bearing_rejects_identical_grid_ids`

Attempts to connect both bearing ends to grid 10. The builder must reject the
self-connection with a clear `distinct grids` error instead of assembling a
degenerate contribution.

#### `test_disk_terms_are_inserted_at_referenced_grid`

Builds a one-grid disk and compares its complete global mass, damping, and
gyroscopic matrices with the local kernels. Damping must be 0.03 times mass,
and stiffness must remain zero.

#### `test_disk_matrices_follow_the_element_coordinate_system`

Defines a cyclically permuted local coordinate system whose disk axis points
along global $x$. The rotated inertia diagonal and the two gyroscopic entries
must move to the corresponding global rotational DOFs, verifying element
orientation rather than merely scalar values.

#### `test_grid_coordinates_are_transformed_to_basic_system`

Expresses a grid in a translated and rotated coordinate system and checks the
resolved global point exactly. The expected `[7, -2, 7]` result distinguishes
direction-vector axes from points offset by the coordinate-system origin.

#### `test_bearing_matrices_rotate_45_degrees_about_global_z`

Defines a bearing frame rotated 45 degrees about global $z$ and starts from
unequal local $x$ and $y$ coefficients. The test compares the transformed
stiffness and damping blocks with explicit global matrices containing the
expected symmetric cross-coupling, while confirming that rotational DOFs
remain untouched.

#### `test_rotated_grid_system_orients_shaft_along_global_x`

Maps the input frame's local $z$ axis onto global $x$ and defines both shaft
endpoints in that frame. The resolved endpoint coordinates must describe an
$x$-directed shaft, and its axial and torsional closed-form blocks must appear
on the global $x$ translation and rotation DOFs. This connects coordinate
resolution to the subsequent element-matrix orientation.

#### `test_invalid_record_references_are_rejected`

Covers duplicate grid IDs, an unknown shaft material, a shaft paired with a
bearing property, and a bearing attached to an unknown grid. Each case must
raise the appropriate type or value error with a useful relationship-specific
message.

#### `test_coincident_shaft_grids_are_rejected`

Places both shaft grids at the same coordinates. The builder must detect the
resulting zero-length element and raise a `coincident` error before attempting
a direction transform or kernel evaluation.

#### `test_builder_strictly_validates_record_types_and_finite_geometry`

Exercises six root-model validation failures: no grids, Boolean and nonpositive
IDs, a grid incorrectly stored as a material, an infinite coordinate, and a
zero spin axis. Besides rejecting bad values, the expected exception types
protect strict API semantics such as not accepting `True` as integer ID 1.

#### `test_builder_requires_a_boolean_stochastic_flag`

Passes integer, null, and string lookalikes as `stochastic`. All must raise a
`TypeError`, ensuring callers explicitly choose a Boolean mode and preventing
Python truthiness from silently changing model construction.

Distinct cross-coupled bearing values are intentional throughout this module.
For example,

$$
\mathbf{K}_b=
\begin{bmatrix}
11 & 12 & 0\\
13 & 14 & 0\\
0 & 0 & 15
\end{bmatrix}
$$

makes transposition and DOF-permutation mistakes visible.

### `test_declarative_api.py`

These end-to-end tests use a three-grid, two-shaft rotor with grounded
bearings and a central disk. Random variants replace shaft diameter, density,
bearing stiffness, and disk mass with distribution references. Output-shape
checks intentionally retain the leading realization dimension, even for a
deterministic build.

#### `test_deterministic_build_resolves_every_distribution_to_its_mean`

Builds the random-reference model with stochastic sampling disabled and
compares all four system matrices with a model containing explicit mean
values. It must produce one realization, preserve the supplied seed as
metadata, and ignore the requested stochastic sample count.

#### `test_stochastic_build_produces_reproducible_model_realizations`

Builds six realizations twice with seed 123. Corresponding matrices must match
exactly between runs, different samples must vary, every resolved definition
must contain ordinary floating-point values, and the matrix stack must have
shape `(6, 18, 18)`.

#### `test_invalid_distribution_specifications_are_rejected`

Parametrizes incomplete normal data, negative standard deviation, reversed
uniform bounds, an unsupported triangular family, nonfinite parameters, and
an unexpected parameter key. Each malformed record must fail validation with
a message identifying the violated distribution contract.

#### `test_multivariate_components_preserve_requested_correlation`

References the two components of one multivariate normal distribution from
`kxx` and `kyy`, draws 2,000 seeded samples, and computes their empirical
correlation. The result must be $0.8\pm0.04$, checking covariance construction
and component lookup without demanding an impossible exact sample statistic.

#### `test_repeated_scalar_reference_reuses_the_same_sample`

Both direct bearing stiffnesses reference the same scalar distribution. For
every realization `kxx` must equal `kyy`, proving that a distribution is drawn
once per sample and reused rather than independently redrawn at each field.

#### `test_invalid_distribution_references_are_rejected`

Checks an unknown distribution ID and an out-of-range component index. Both
must be rejected during resolution, so invalid symbolic references cannot
survive into numerical assembly.

#### `test_stochastic_sample_count_must_be_a_positive_integer`

Tries zero, a negative integer, a float, a string, and a Boolean sample count.
All must raise `ValueError`, enforcing a strictly positive integer and again
excluding Boolean values despite their integer subclass relationship.

#### `test_builder_requires_a_model_definition`

Passes a dictionary instead of the public `ModelDefinition` root record. The
builder must raise `TypeError`, making the supported declarative boundary
explicit rather than accepting loosely shaped mappings.

#### `test_modal_solver_returns_real_modes_for_each_sample`

Solves four modes for each of three stochastic realizations after two DOFs are
fixed. It checks positive real eigenvalues, eigenvector and frequency shapes,
the conversion $f=\sqrt{\lambda}/(2\pi)$, and variation between realizations.

#### `test_campbell_solver_returns_modes_at_every_speed`

Solves four state-space modes at 0, 100, and 200 rad/s. The test verifies speed
conversion to hertz, realization/speed/mode axes, doubled state-space
eigenvector size, and finite frequencies at every requested operating point.

#### `test_campbell_solver_rejects_invalid_mode_counts`

Passes zero, a negative value, a float, a Boolean, and a string as `modes`.
Each must produce the documented positive-integer error before the Campbell
eigensolution starts.

#### `test_campbell_solver_rejects_more_modes_than_are_available`

Requests 100 modes from the constrained reference model. The solver must
report that the request exceeds the available state-space modes rather than
silently truncating or returning a misleading shape.

#### `test_frequency_response_accepts_a_full_model_force`

Applies a unit force in the full 18-DOF model while constraining two DOFs. The
solver must reduce that vector internally and return finite complex responses
for two frequencies with shape `(1, 2, 16)`.

#### `test_time_response_returns_displacement_and_velocity`

Integrates a short five-point response to a full-model unit force. It checks
the returned time vector, reduced 16-DOF displacement and velocity histories,
finite values, and the expected zero initial displacement.

#### `test_time_response_reduces_a_full_model_initial_state`

Supplies concatenated displacement and velocity values for all 18 model DOFs.
At the initial time, the returned state must exactly equal those arrays with
the two constrained entries removed, proving consistent full-to-free mapping
for both halves of the state vector.

#### `test_solver_validates_entry_point_arguments`

Exercises the unified entry point without a model, with the wrong model type,
with an unknown analysis, with out-of-range or noninteger fixed DOFs, and with
a misspelled keyword. The test checks that each public misuse fails early and
with a targeted message.

#### `test_modal_rejects_cross_coupled_nonsymmetric_stiffness`

Creates a bearing for which `kxy != kyx`. Because the real generalized modal
solver assumes symmetric stiffness, it must reject this model explicitly
instead of applying an invalid symmetric eigensolver.

#### `test_modal_keeps_small_physical_modes_in_a_wide_spectrum`

Uses a one-grid mass-bearing system whose two eigenvalues are $10^{-4}$ and
$10^8$. Both must survive modal filtering, preventing a relative threshold
based on the largest eigenvalue from discarding a small but physical mode.

#### `test_response_solvers_validate_load_and_time_shapes`

Passes a force vector whose length matches neither the full nor reduced model,
then supplies only one time point to the transient solver. Both malformed
inputs must raise focused errors before numerical solution.

### `test_plotting.py`

This module uses Matplotlib's non-interactive `Agg` backend, so visualization
tests render in local and continuous-integration environments without opening
a graphical window.

#### `test_plot_campbell_draws_one_line_per_tracked_mode`

Plots the second realization from a synthetic two-sample Campbell result. It
checks that each modal column becomes one line, rad/s values are converted to
rpm, the selected sample supplies the y values, branch labels remain in mode
order, and both returned objects refer to the same figure.

#### `test_plot_campbell_draws_stochastic_statistics_and_bands`

Adds an outlying third realization so mean and median branches differ, then
checks both aggregations independently. It also verifies one confidence-band
collection per mode, exact sampled minimum and maximum curves, consistent
branch labels, and suppression of the legend when requested.

#### `test_plot_campbell_validates_plot_options`

Exercises an out-of-range sample, a noninteger sample, an unsupported speed
unit, invalid statistics and confidence levels, and non-Boolean display
options. Each case must fail with a focused exception before a figure is
created.

### `test_analytical_benchmarks.py`

#### `test_simply_supported_beam_first_four_bending_frequencies`

Builds a 16-element Euler--Bernoulli shaft and compares one member of each of
the first four degenerate bending pairs with

$$
f_m=\frac{m^2\pi}{2L^2}
\sqrt{\frac{EI}{\rho A}},
\qquad m=1,2,\ldots.
$$

For the reference values used by the tests,

$$
A=\frac{\pi(0.01)^2}{4}
=7.853981634\times10^{-5}\ \mathrm{m^2},
$$

$$
I=\frac{\pi(0.01)^4}{64}
=4.908738521\times10^{-10}\ \mathrm{m^4}.
$$

Substituting $E=2.0\times10^{11}\ \mathrm{Pa}$,
$\rho=7850\ \mathrm{kg/m^3}$, and $L=1.0\ \mathrm{m}$ gives

$$
f_m=
\frac{m^2\pi}{2(1.0)^2}
\sqrt{
\frac{(2.0\times10^{11})(4.908738521\times10^{-10})}
{(7850)(7.853981634\times10^{-5})}
}
=m^2(19.821661494\ \mathrm{Hz}).
$$

The resulting analytical references are:

| Mode $m$ | Analytical expression | Analytical frequency (Hz) |
|---:|---:|---:|
| 1 | $1^2(19.821661494)$ | 19.821661494 |
| 2 | $2^2(19.821661494)$ | 79.286645975 |
| 3 | $3^2(19.821661494)$ | 178.394953444 |
| 4 | $4^2(19.821661494)$ | 317.146583901 |

These values come only from the closed-form equation; they are not numerical
finite-element results.

Rotary inertia is disabled so the finite-element assumptions match the
analytical equation. With the six-DOF grid order `x, y, z, tx, ty, tz`, the
first grid fixes `x, y, z, tz` and the final grid fixes `x, y`; bending
rotations remain free, while axial and torsional rigid motion is removed.
Every selected numerical frequency must have relative error below $0.1\%$:

$$
\varepsilon_m=
\frac{\left|f_{m,\mathrm{num}}-f_{m,\mathrm{ana}}\right|}
{f_{m,\mathrm{ana}}}.
$$

### `test_convergence.py`

Both tests solve the same simply supported Euler--Bernoulli beam with 2, 4,
8, and 16 shaft elements. Material, section, total length, formulation, and
supports stay fixed, so changes in error measure discretization alone.

#### `test_first_frequency_is_reasonable_at_each_mesh`

For each of the four mesh sizes, selects the first physical bending frequency
from its degenerate pair and compares it with the independent analytical
value. Even the coarsest tested mesh must remain within $1\%$, catching gross
errors that could still appear to improve under refinement.

#### `test_first_frequency_converges_under_mesh_refinement`

Compares the sequence of first-mode errors. The 16-element error must be less
than the 2-element error and below $0.01\%$, demonstrating useful convergence
without requiring strictly monotonic improvement at every step, which could
make the test sensitive to harmless floating-point changes.

## Numerical tolerances

Choose a tolerance from the mathematical or physical property under test:

- use tight relative tolerances for symmetry, skew-symmetry, and exact matrix
  blocks;
- detect rigid modes relative to the largest stiffness eigenvalue rather than
  with one dimensional absolute threshold;
- account for finite-element discretization error in analytical comparisons;
- test both an improving trend and a final accuracy target in convergence
  studies;
- use seeded generators for reproducible stochastic assertions without
  asserting a particular platform-dependent sample statistic unnecessarily.

Do not loosen a tolerance solely to silence a failure. First check formulation
assumptions, units, coordinate transformations, DOF ordering, supports, mode
selection, and floating-point conditioning.

## Adding tests

When adding a record, element formulation, distribution, or analysis:

1. Test the smallest stateless numerical kernel independently.
2. Test its coordinate transformation and contribution to a minimal global
   model.
3. Test deterministic and stochastic resolution when the record has numeric
   parameters.
4. Exercise the feature through the public `ModelBuilder` and `Solver.solve`
   workflow.
5. Compare against an independent analytical or published result when one is
   available.
6. Add a convergence study when the result depends on mesh refinement.

Keep benchmark formulas independent of Spinniped's implementation. Copying a
production formula into a test can reproduce the same defect and create false
confidence.

## Verification boundaries

The suite verifies the current dense, linear rotor model and plain Monte Carlo
sampling. It does not yet constitute validation for spatial random fields,
nonlinear dynamics, Campbell branch tracking through mode crossings, sparse
assembly, or every possible bearing and coordinate-system orientation. Add
focused verification before relying on a new formulation in those regimes.
