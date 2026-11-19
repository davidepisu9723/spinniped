# Testing and verification

Spinniped's tests serve two purposes: preventing software regressions and
checking that the finite-element equations reproduce expected mechanical
behavior. This file documents verification behavior rather than deriving the
underlying mechanics; theoretical development belongs in the
[rotordynamics theory manual](../doc/theory-manual/00_theory_manual.md). The
suite exercises the declarative records, stateless local-matrix kernels,
deterministic and stochastic assembly, every solver route, and an independent
beam benchmark.

## Table of contents

- [Running the suite](#running-the-suite)
- [Test summary](#test-summary)
- [Test details](#test-details)
- [Numerical tolerances](#numerical-tolerances)
- [Adding tests](#adding-tests)
- [Verification boundaries](#verification-boundaries)

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

# Eigenvalue extraction test class
python -m pytest tests/test_eigenvalue_extraction.py::TestEigenvalueExtraction

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
| [`test_assembly.py`](#test_assemblypy) | 18 | Record validation, coordinate transforms, and global assembly |
| [`test_declarative_api.py`](#test_declarative_apipy) | 33 | Distribution resolution and core public solver routes |
| [`test_eigenvalue_extraction.py`](#test_eigenvalue_extractionpy) | 14 | Partial eigenpairs, root budgets, normalization, tracking, and error propagation |
| [`test_speed_dependence.py`](#test_speed_dependencepy) | 20 | Tabulated rotor-speed interpolation, stochastic curves, matrix deltas, and load deltas |
| [`test_load_response.py`](#test_load_responsepy) | 11 | Declarative load phasors and steady-response routes |
| [`test_plotting.py`](#test_plottingpy) | 13 | Rotor-section, frequency-response, and Campbell rendering and input validation |
| [`test_analytical_benchmarks.py`](#test_analytical_benchmarkspy) | 1 | Comparison with the Euler--Bernoulli beam solution |
| [`test_jeffcott_benchmark.py`](#test_jeffcott_benchmarkpy) | 2 | Jeffcott frequency and critical-speed comparison |
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

The section properties used by the fixtures are calculated independently from
these inputs. `shaft_properties` exposes the inputs and derived section
properties; `shaft_definition` creates one Timoshenko shaft; and
`shaft_local_matrices` evaluates its stiffness, mass, and gyroscopic kernels.
`build_uniform_beam` creates an arbitrary Euler or Timoshenko mesh through the
public declarative API. The support, analytical-frequency, and bending-pair
helpers are shared by the benchmark and convergence tests. In particular, the
analytical helper is independent of Spinniped's matrix kernels, preventing the
same implementation error from appearing on both sides of a comparison.

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

#### `test_bearing_element_rejects_invalid_local_inputs`

Parametrizes invalid bearing-element IDs, grid/property references, optional
second-grid values, and coordinate-system references. Each malformed element
must fail immediately during `BearingElement` construction.

#### `test_disk_element_rejects_invalid_local_inputs`

Parametrizes invalid disk-element IDs, grid/property references, and frame
references. Each malformed element must fail immediately during `DiskElement`
construction.

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

#### `test_coordinate_system_rejects_invalid_local_inputs`

Parametrizes invalid coordinate-system IDs, vector lengths and values, a zero
x-axis, and collinear directions. Each malformed frame must fail immediately
during `CoordinateSystem` construction.

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

#### `test_grid_rejects_invalid_local_inputs`

Parametrizes invalid grid IDs, nonfinite coordinates, and invalid coordinate-
system references. Each malformed grid must fail immediately during `Grid`
construction.

#### `test_material_rejects_invalid_local_inputs`

Parametrizes invalid material IDs, nonfinite or nonpositive elastic data, and
an out-of-range Poisson ratio. Each malformed material must fail immediately
during `Material` construction.

#### `test_shaft_property_rejects_invalid_local_inputs`

Parametrizes invalid IDs, material references, annular geometry, theory, and
rotary-inertia switches. Each malformed shaft property must fail immediately
during `ShaftProperty` construction.

#### `test_bearing_property_rejects_invalid_local_inputs`

Parametrizes invalid bearing-property IDs and malformed constant or random
coefficients. Each malformed property must fail immediately during
`BearingProperty` construction.

#### `test_disk_property_rejects_invalid_local_inputs`

Parametrizes invalid disk-property IDs and negative mass, inertia, or damping.
Each malformed property must fail immediately during `DiskProperty`
construction.

#### `test_shaft_element_rejects_invalid_local_inputs`

Parametrizes invalid shaft-element IDs and references plus identical endpoint
grids. Each malformed element must fail immediately during `ShaftElement`
construction.

#### `test_builder_strictly_validates_record_types_and_finite_geometry`

Exercises the remaining root-model validation failures: no grids, a grid
incorrectly stored as a material. The builder must retain these completeness
and mutable-collection checks.

#### `test_model_definition_requires_list_collections`

Passes a tuple to each record-collection field in turn. Every collection must
require an explicit mutable list so definitions preserve their incremental
construction contract without implicit conversion.

#### `test_model_definition_rejects_invalid_spin_axis`

Parametrizes wrong-length, nonfinite, nonnumeric, and zero spin axes. Each
malformed atomic vector must fail during `ModelDefinition` construction.

#### `test_builder_requires_a_boolean_stochastic_flag`

Passes integer, null, and string lookalikes as `stochastic`. All must raise a
`TypeError`, ensuring callers explicitly choose a Boolean mode and preventing
Python truthiness from silently changing model construction.

Distinct cross-coupled bearing values are intentional throughout this module because they make transposition and DOF-permutation mistakes visible.

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
an unexpected parameter key, plus invalid identity and field types. Each
malformed record must fail immediately during `RandomDistribution`
construction with a message identifying the violated contract.

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

#### `test_load_records_are_public_immutable_declarative_data`

Constructs harmonic, circular rotating, and unbalance records through the
package root. It checks their component, phase, sense, and coordinate-system
defaults, verifies frozen slotted storage, and confirms that separate
`ModelDefinition` instances do not share a mutable load list.

#### `test_fixed_vector_fields_require_tuples`

Passes mutable lists to coordinate-system vectors, harmonic force and moment
vectors, and the model spin axis. Every atomic vector must reject the list at
record construction and require an explicit tuple.

#### `test_builder_resolves_and_preserves_declarative_loads`

References registered distributions from harmonic-force components, rotating
force magnitude, and mass eccentricity. A deterministic build must resolve
those values into `BuiltModel.definitions` while leaving all four structural
matrix families identical to a load-free model.

#### `test_builder_rejects_invalid_load_records`

Checks the shared load-ID namespace, accepted record types, and grid and
coordinate-system existence.

#### `test_harmonic_load_rejects_invalid_local_inputs`

Parametrizes invalid harmonic-load identity, references, vector shapes,
components, and phase. Each malformed load must fail immediately during
`HarmonicLoad` construction.

#### `test_synchronous_load_inherits_the_harmonic_vector_contract`

Constructs a fractional-order synchronous load and verifies that it is an
immutable `HarmonicLoad` specialization with the same force, moment, phase,
and coordinate-system fields.

#### `test_synchronous_load_rejects_invalid_harmonic_order`

Parametrizes Boolean, textual, infinite, zero, and negative orders. Every
invalid value must fail during `SynchronousLoad` construction.

#### `test_rotating_load_rejects_invalid_local_inputs`

Parametrizes invalid rotating-load identity, references, magnitude, phase, and
rotation sense. Each malformed load must fail immediately during
`RotatingLoad` construction.

#### `test_unbalance_rejects_invalid_local_inputs`

Parametrizes invalid unbalance identity, references, mass eccentricity, and
phase. Each malformed record must fail immediately during `Unbalance`
construction.

#### `test_unbalance_plane_must_be_normal_to_the_spin_axis`

Defines a load frame whose local z-axis is transverse to the rotor spin axis.
The builder must reject it because that frame's local x-y plane is not the
physical unbalance plane.

#### `test_unbalance_plane_accepts_an_antiparallel_frame_axis`

Builds an unbalance whose load-frame z-axis opposes the model spin axis. The
frame is accepted because its x-y plane is still the rotor plane; the axis
alignment sign remains available for the later force-phasor assembly.

#### `test_modal_solver_returns_real_modes_for_each_sample`

Solves four modes for each of three stochastic realizations after two DOFs are
fixed. It checks positive real eigenvalues, eigenvector and frequency shapes,
the conversion $f=\sqrt{\lambda}/(2\pi)$, and variation between realizations.

#### `test_campbell_solver_returns_modes_at_every_speed`

Solves four state-space modes at 0, 100, and 200 rad/s. The test verifies speed
conversion to hertz, realization/speed/mode axes, doubled state-space
eigenvector size, finite frequencies at every requested operating point, and
empty critical-speed axes when no harmonics were requested.

#### `test_campbell_calculates_sampled_speed_interpolated_critical_speeds`

Solves a two-direction 10 Hz mass--bearing oscillator on a deliberately coarse
speed sequence and requests 1x and 2x harmonics. Both modal branches must cross
at 10 Hz for 1x and 5 Hz for 2x. Neither value is a supplied speed, proving
that the solver linearly interpolates existing residuals without inserting or
evaluating another speed.

#### `test_campbell_warns_only_when_tracking_starts_at_zero_speed`

Checks that tracked Campbell analysis emits the documented degeneracy warning
when its first speed is exactly zero, while a nonzero start and disabled mode
tracking remain warning-free.

#### `test_campbell_preserves_stochastic_critical_speed_correspondence`

Builds four realizations of a two-direction mass--bearing oscillator whose two
stiffness distributions occupy nonoverlapping frequency ranges. Every sample's
1x roots must equal the natural frequencies calculated from its resolved
`kxx` and `kyy`, proving that the sample and tracked-mode axes remain associated
through critical-speed detection.

#### `test_critical_speed_detection_retains_crossing_order_and_missing_values`

Passes synthetic tracked frequencies directly to the crossing detector. The
first sample crosses one harmonic three times and the second never crosses.
The result must preserve all three roots in speed order, record per-sample
counts, and fill the absent sample's crossing axis with `NaN`.

#### `test_campbell_validates_critical_speed_options`

Rejects nonpositive and duplicate harmonic indices, disabled mode tracking,
and a descending speed sequence. These conditions prevent ambiguous physical-mode
association or ordered interpolation before the state eigensolutions begin.

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

### `test_eigenvalue_extraction.py`

`TestEigenvalueExtraction` has 14 test methods and 25 parametrized cases. It
compares partial modal and damped gyroscopic roots with independent dense
solutions and checks physical equation residuals and mass normalization.

The root-budget checks verify that a six-mode modal request extracts six
eigenpairs and a six-mode Campbell request makes one ARPACK call for 14 state
roots, including conjugate partners and one extra pair. Other cases cover
rigid coordinates, repeated roots, negative stiffness, small systems, and
tracking within the extracted candidates. Singular factorizations and ARPACK
nonconvergence must propagate directly without a recovery solve.

### `test_speed_dependence.py`

This module verifies that tabulated coefficients remain declarative model data,
are resolved into transient operating-point contributions, and never add a
permanent speed dimension to the built matrices.

#### `test_speed_dependent_value_is_public_frozen_and_linearly_interpolated`

Checks the public frozen, slotted record and evaluates an unsorted sequence
containing repeated and negative speeds. Exact knots and intermediate points
must preserve request order and use linear interpolation.

#### `test_speed_dependent_value_requires_tuple_inputs`

Passes a list for each speed-table sequence in turn. Both must be rejected at
record construction because `speeds` and `values` require immutable tuples.

#### `test_speed_dependent_value_rejects_malformed_inputs`

Parametrizes insufficient and mismatched tables, repeated or descending knots,
nonfinite entries, Boolean entries, and malformed random references. Every
case must fail immediately during `SpeedDependentValue` construction.

#### `test_speed_dependent_value_rejects_extrapolation`

Accepts both table endpoints and independently rejects requests immediately
below and above the available signed rotor-speed interval.

#### `test_builder_keeps_base_matrices_three_dimensional_and_adds_deltas`

Combines constant and tabulated bearing coefficients. The permanent K, M, C,
and G arrays must remain `(samples, ndof, ndof)`, K and C must contain only the
constant terms, and a modal solve must add the interpolated stiffness exactly
once.

#### `test_rotated_two_grid_bearing_transforms_speed_dependent_deltas`

Combines a rotated bearing frame with two-grid relative motion. Interpolated K
and C deltas must act on the rotated global DOFs and retain the exact
`[[+v, -v], [-v, +v]]` block pattern.

#### `test_deterministic_build_resolves_random_speed_ordinates_to_means`

References two scalar distributions from the table ordinates. A deterministic
build must preserve the speed knots, replace both references by their means,
and interpolate the resulting deterministic curve in modal analysis.

#### `test_multivariate_speed_ordinates_are_correlated_and_reproducible`

Uses two components of one multivariate normal distribution as the endpoint
ordinates. Seeded ensembles must be identical and their sampled endpoint
correlation must reproduce the requested correlation without expanding the
built stiffness array to four dimensions.

#### `test_repeated_random_reference_reuses_one_draw_at_every_speed_knot`

Places the same scalar distribution reference at both knots. Each realization
must contain two equal resolved ordinates and consequently return the same
Campbell frequency at both operating speeds.

#### `test_campbell_uses_interpolated_stiffness_and_damping_at_every_speed`

Compares the complex eigenvalue of a damped single-DOF system at three speeds
with the closed-form roots obtained from independently interpolated K and C.

#### `test_frequency_response_uses_rotor_speed_not_forcing_frequency`

Sweeps two excitation frequencies at one rotor speed and compares both complex
responses with the same interpolated K and C. This distinguishes rotor-speed
tables from the independent harmonic-frequency grid.

#### `test_time_response_uses_one_interpolated_operating_point`

Compares a tabulated model at one scalar speed with an otherwise identical
constant-coefficient model. Their integrated displacement and velocity
histories must coincide.

#### `test_unbalance_response_uses_literal_signed_speed_for_bearing_tables`

Uses asymmetric negative and positive knots in an isotropic two-DOF rotor.
The synchronous response is checked against the signed-speed interpolation,
including reverse-force sense and the existing speed-squared load law.

#### `test_campbell_preflights_complete_speed_grid_before_first_solve`

Requests one valid point followed by an out-of-range point while replacing the
linear solver with a sentinel. Extrapolation must be rejected before even the
first otherwise-valid Campbell point is solved.

#### `test_harmonic_load_components_interpolate_transform_and_superimpose`

Interpolates force and moment components, combines them with a constant load,
applies a common complex phase, and rotates the complete nodal vector from its
declared coordinate system into global coordinates.

#### `test_frequency_response_interpolates_rotating_load_at_rotor_speed`

Checks a negative-sense circular load whose magnitude is tabulated at rotor
speed. The resulting two-component phasor and response are compared with the
closed-form dynamic-stiffness solution.

#### `test_synchronous_load_interpolates_amplitude_over_speed_sweep`

Uses a 2x synchronous load over three signed rotor speeds. Every force ordinate
must be interpolated at rotor speed while excitation frequencies are derived
independently from twice the absolute speeds.

#### `test_selected_load_is_preflighted_but_unselected_tables_are_ignored`

Shows that an out-of-range table does not affect an analysis that did not
select its load ID. Selecting it must reject the operating speed before the
first linear solve.

#### `test_rotating_load_speed_table_must_remain_nonnegative`

Constructs a rotating load with a magnitude table containing a negative
ordinate. The record must reject it immediately even when another ordinate is
valid.

#### `test_builder_rejects_mutated_speed_dependent_distribution_parameters`

Mutates a previously valid distribution parameter dictionary to contain a
speed table. The builder must revalidate mutable nested data and reject the
unsupported placement before sampling.

#### `test_model_definition_rejects_speed_dependent_spin_axis_at_declaration`

Places a speed table in the immutable model spin axis. `ModelDefinition` must
reject the unsupported component immediately.

#### `test_unbalance_rejects_speed_dependence_at_declaration`

Parametrizes the mass-eccentricity and phase fields with speed tables. Both
unsupported placements must fail immediately during `Unbalance` construction.

#### `test_grid_rejects_speed_dependent_coordinates_at_declaration`

Places a speed table in a grid coordinate. The `Grid` record must reject it
immediately because coordinates accept scalar parameters only.

#### `test_material_rejects_speed_dependent_parameters_at_declaration`

Places a speed table in material density. The `Material` record must reject it
immediately because its fields accept scalar parameters only.

#### `test_shaft_property_rejects_speed_dependence_at_declaration`

Places a speed table in shaft diameter. The `ShaftProperty` record must reject
it immediately because shaft-property fields accept scalar parameters only.

#### `test_disk_property_rejects_speed_dependence_at_declaration`

Places a speed table in disk mass. The `DiskProperty` record must reject it
immediately because disk-property fields accept scalar parameters only.

### `test_load_response.py`

This module verifies conversion of declarative load records into complex nodal
phasors and their use by both steady-response routes.

#### `test_harmonic_load_matches_an_equivalent_global_force_vector`

Rotates nonzero force and moment triples from a custom frame, applies a common
phase, and compares the declarative response with the existing raw-vector
route. The exact assembled six-DOF phasor must also match independently
transformed components.

#### `test_synchronous_load_derives_frequency_from_order_and_speed`

Solves a fractional-order load at negative rotor speed and compares the single
returned frequency and response with the closed-form dynamic-stiffness result.

#### `test_synchronous_loads_of_one_order_superimpose`

Solves two 2x records separately and together. The combined response must equal
their sum while preserving the caller's load-ID order.

#### `test_synchronous_load_sweeps_paired_speeds_and_frequencies`

Checks negative, zero, and positive speeds against direct dynamic-stiffness
solutions and verifies the paired speed, speed-Hz, frequency, and response
metadata.

#### `test_synchronous_frequency_response_rejects_ambiguous_sources`

Rejects an explicit frequency for a synchronous load, mixing synchronous and
ordinary harmonic records, combining different orders, and omitting the
frequency grid for an ordinary harmonic load.

#### `test_synchronous_frequency_response_rejects_invalid_speed_sweeps`

Parametrizes empty, scalar, and nonfinite speed inputs. Every malformed sweep
must fail before response assembly.

#### `test_rotating_load_sense_phase_and_superposition_are_preserved`

Checks positive and negative circular phasors at a nonzero initial angle, then
verifies that selecting both records produces the sum of their independent
responses while preserving caller load-ID order in the metadata.

#### `test_frequency_response_uses_each_samples_resolved_load`

Builds four realizations with random force magnitude and phase. Every response
must match its own resolved record rather than reusing sample zero's load.

#### `test_unbalance_response_matches_signed_synchronous_dynamic_stiffness`

Compares forward, zero, and reverse-speed results against direct solutions of
$K-\omega^2M+i\omega(C+\Omega G)$ with $\omega=|\Omega|$ and
$F=me\Omega^2$. A nonzero custom gyroscopic block makes signed speed part of
the numerical check.

#### `test_unbalance_uses_frame_alignment_to_follow_the_rotor`

Compares parallel and antiparallel load-frame z-axes. After transformation,
both representations of the same phase-zero unbalance must produce the same
global rotating force.

#### `test_response_load_selection_is_validated`

Rejects missing or conflicting frequency-response sources, empty, scalar,
Boolean, duplicate, unknown, and analysis-incompatible load selections,
negative forcing frequencies, and malformed unbalance speed arrays.

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

#### `test_plot_campbell_draws_deterministic_harmonics_and_critical_speeds`

Plots 1x and 2x synchronous lines for one realization and checks the expected
critical markers by harmonic and tracked mode. A mode without a 1x crossing
must not receive a marker, while marker coordinates must satisfy both the
stored critical speed and its synchronous-frequency relation.

#### `test_plot_campbell_draws_stochastic_critical_speed_statistics`

Plots three realizations through the mean-statistic route. It checks the mean
critical point, the capped diagonal one-sample-standard-deviation error bar,
the optional underlying sample scatter, and legend suppression. Missing
crossings remain excluded from the finite statistics.

#### `test_plot_frequency_response_selects_global_dofs_and_frequency_axis`

Plots an ordinary harmonic result with nonadjacent free DOFs. It checks that
the selected global indices preserve caller order, custom labels follow that
order, complex values become magnitudes, and excitation frequency is the
default horizontal coordinate.

#### `test_plot_frequency_response_uses_speed_for_synchronous_results`

Parametrizes synchronous harmonic and unbalance results over a signed speed
grid. Both routes must automatically select rotor speed and convert rad/s to
the requested display unit without changing the response ordering.

#### `test_plot_frequency_response_draws_stochastic_statistics_and_bands`

Plots the mean magnitude of a two-sample complex response with an empirical
confidence band and sampled extrema. It checks the central, minimum, and
maximum ordinates, dashed extrema, information text, and legend suppression.

#### `test_plot_frequency_response_validates_results_and_dof_selection`

Rejects a nonmapping input, an incompatible analysis, a malformed response
shape, unavailable or duplicate global DOFs, mismatched labels, and a speed
axis without one speed per response point.

#### `test_plot_frequency_response_validates_plot_options`

Exercises invalid samples, statistics, confidence levels, axis choices, speed
units, and Boolean display options. Each case must fail with a focused error
before plotting.

#### `test_plot_rotor_draws_exact_shaft_sections_disks_and_bearings`

Builds a stepped rotor containing one hollow shaft section, one solid shaft
section, a central disk, and two grounded bearings. It checks every shaft
rectangle against the resolved axial endpoints and diameters, checks the bore
against the inner diameter, and verifies that the disk and both halves of each
bearing symbol are present. Axis labels and the component legend are also
checked.

#### `test_plot_rotor_selects_resolved_sample_and_reuses_axes`

Builds three stochastic realizations of a random shaft diameter and plots the
second one into caller-owned axes. The shaft height must equal that sample's
resolved property value. The test also verifies figure reuse, preservation of
an existing title, optional grid-ID labels, and legend suppression.

#### `test_plot_rotor_validates_model_sample_and_longitudinal_geometry`

Rejects an unbuilt definition, a Boolean sample, an out-of-range sample, and a
shaft not parallel to the model spin axis. The last case ensures an oblique
element cannot silently be shortened by projection and presented as a valid
longitudinal section.

#### `test_plot_campbell_validates_plot_options`

Exercises an out-of-range sample, a noninteger sample, an unsupported speed
unit, invalid statistics and confidence levels, and non-Boolean display
options. Each case must fail with a focused exception before a figure is
created.

### `test_analytical_benchmarks.py`

This module compares a finite-element shaft with an independent analytical
reference. The derivation, boundary-condition discussion, and free-body
diagrams are kept in
[the flexural-vibration chapter](../doc/theory-manual/06_flexural_vibration.md).
The executable reference is
[`01_uniform_beam_frequencies.py`](../benchmark/01_uniform_beam_frequencies.py).

#### `test_simply_supported_beam_first_four_bending_frequencies`

Builds a 16-element Euler--Bernoulli shaft and compares one member of each of
the first four degenerate bending pairs with independently calculated,
simply-supported beam frequencies. Rotary inertia is disabled so the numerical
and analytical models use the same assumptions.

The first grid constrains `x`, `y`, `z`, and `tz`; the final grid
constrains `x` and `y`. Bending rotations remain free, while axial and
torsional rigid motion is removed.

The immutable analytical reference values are:

| Bending mode | Analytical frequency (Hz) |
|---:|---:|
| 1 | 19.821661494 |
| 2 | 79.286645975 |
| 3 | 178.394953444 |
| 4 | 317.146583901 |

These values are produced independently of Spinniped's matrix kernels. Every
selected numerical frequency must have relative error below $0.1\%$.

### `test_jeffcott_benchmark.py`

This module compares an extended flexible-bearing Jeffcott reference with an
equivalent finite-element model. The analytical derivation and diagrams are in
[the Jeffcott chapter](../doc/theory-manual/07_jeffcott_rotor.md).
The runnable comparison is
[`02_flexible_bearing_jeffcott_comparison.py`](../benchmark/02_flexible_bearing_jeffcott_comparison.py).

The fixed benchmark inputs are:

| Quantity | Value |
|---|---:|
| Young's modulus | $2.0\times10^{11}\ \mathrm{Pa}$ |
| Shaft diameter | $0.02\ \mathrm{m}$ |
| Shaft length | $1.0\ \mathrm{m}$ |
| Disk mass | $5.0\ \mathrm{kg}$ |
| Disk diametral inertia | $0.025\ \mathrm{kg\,m^2}$ |
| Disk polar inertia | $0.05\ \mathrm{kg\,m^2}$ |
| Bearing stiffness | $1.0\times10^6\ \mathrm{N/m}$ |

The resulting analytical reference values are:

| Quantity | Analytical value |
|---|---:|
| Shaft second moment | $7.853981634\times10^{-9}\ \mathrm{m^4}$ |
| Translation stiffness | $72659.042323\ \mathrm{N/m}$ |
| Rotation stiffness | $18164.760581\ \mathrm{N\,m/rad}$ |
| Cylindrical frequency at rest | $19.185802264\ \mathrm{Hz}$ |
| Conical frequency at rest | $135.664108836\ \mathrm{Hz}$ |
| Cylindrical 1x critical speed | $1151.148135859\ \mathrm{rpm}$ |
| Backward-conical 1x critical speed | $4699.542585350\ \mathrm{rpm}$ |

The comparable finite-element model uses ten Euler--Bernoulli shaft elements,
a disk at the exact midpoint, negligible positive shaft density, and grounded
flexible bearings at both ends. Shaft rotary inertia is disabled so the disk is
the only gyroscopic source. Axial and torsional DOFs are removed because they
are outside the reduced Jeffcott model.

#### `test_flexible_bearing_jeffcott_frequencies_match_finite_element_model`

Solves four real modes at rest and compares the degenerate cylindrical and
conical pairs with the tabulated analytical frequencies. An explicit zero
modal threshold retains the disk modes in the deliberately wide spectrum
created by the nearly massless shaft. Analytical and finite-element
frequencies must agree to relative tolerance $10^{-5}$.

#### `test_jeffcott_gyroscopic_shift_and_critical_speeds_match_finite_elements`

Solves at rest, at an intermediate speed, and at both tabulated analytical
critical speeds. The finite-element frequencies must reproduce the fixed
cylindrical pair and the oppositely shifting conical branches. The test also
checks that the cylindrical and backward-conical frequencies intersect the 1x
line at their respective analytical critical speeds.

Independent frequency sorting is used because the conical pair is exactly
degenerate at rest, making its initial MAC labels arbitrary. All compared
frequencies must agree to relative tolerance $10^{-5}$.

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
nonlinear dynamics, general Campbell branch identity through crossings, sparse
assembly, or every possible bearing and coordinate-system orientation. Add
focused verification before relying on a new formulation in those regimes.
