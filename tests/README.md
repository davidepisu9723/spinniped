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

## Test organization

### `conftest.py`

The shared fixtures describe a reference uniform shaft using independent
Python values and declarative records. The reference properties are:

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

`build_uniform_beam` creates `Grid`, `Material`, `ShaftProperty`, and
`ShaftElement` records, then assembles them through `ModelBuilder`. The
analytical frequency helper remains independent of the package's matrix
kernels so it cannot reproduce an implementation error on both sides of a
comparison.

### `test_element_matrices.py`

These tests call the stateless functions in `spinniped.stiffness`,
`spinniped.mass`, `spinniped.damping`, and `spinniped.gyroscopic` directly.
They verify:

- shaft matrix dimensions and finite entries;
- symmetry of shaft stiffness and mass and skew-symmetry of gyroscopic terms;
- six rigid-body modes of an unconstrained shaft stiffness matrix;
- positive definiteness of the shaft mass matrix;
- closed-form axial and torsional stiffness blocks;
- linear density scaling of mass and gyroscopic matrices;
- bearing coefficient placement and zero bearing mass/gyroscopic terms;
- disk lumped mass, polar gyroscopic coupling, and zero disk stiffness;
- mass-proportional shaft and disk damping;
- rejection of zero length and invalid annular geometry.

Local-kernel tests isolate formulation errors before coordinate transformation
or global assembly can obscure them.

### `test_assembly.py`

Assembly tests build deliberately small `ModelDefinition` instances and check
the resulting `BuiltModel`. They cover:

- equality between a one-shaft global matrix and its local kernel for an
  aligned element;
- accumulation from two shafts at a shared grid;
- absence of direct coupling between unconnected grids;
- mapping of non-consecutive grid IDs in record order;
- bearing stiffness and damping insertion at the referenced grid;
- disk mass, damping, and gyroscopic insertion;
- conversion of local grid coordinates into the global basic system;
- duplicate IDs, incompatible properties, and unknown references;
- rejection of a shaft whose two grids are coincident.

Distinct cross-coupled bearing values are intentional. For example,

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

These end-to-end tests cover the public records, sampling behavior, builder,
and unified `Solver.solve` entry point. They verify:

- deterministic resolution of normal means and uniform midpoints;
- equivalence between a distribution-at-mean model and explicit mean values;
- reproducible seeded Monte Carlo realizations;
- exact reuse of scalar distribution references;
- multivariate normal sampling with the requested component correlation;
- variation between stochastic matrix samples and storage of resolved
  definitions;
- validation of incomplete, negative-deviation, reversed-bound, and unknown
  distribution specifications;
- validation of stochastic sample counts and builder input types;
- modal output for every realization;
- Campbell output for every requested speed;
- complex frequency responses to a full-model force vector;
- time-domain displacement and velocity histories;
- rejection of missing models, unknown analysis names, invalid constrained
  DOFs, malformed forces, and insufficient time points.

Shape assertions explicitly retain the leading realization dimension, even
for deterministic builds.

### `test_analytical_benchmarks.py`

The analytical benchmark compares the first four bending frequencies of a
uniform simply supported shaft with the Euler-Bernoulli solution

$$
f_m=\frac{m^2\pi}{2L^2}
\sqrt{\frac{EI}{\rho A}},
\qquad m=1,2,\ldots.
$$

The numerical model selects `theory="euler"` and disables bending rotary
inertia to match the assumptions of this equation.

With the six-DOF grid order

```text
x, y, z, tx, ty, tz
```

the support conditions are:

| End | Fixed DOFs | Purpose |
|---|---|---|
| First grid | `x, y, z, tz` | Simple transverse support and removal of axial/torsional rigid motion |
| Final grid | `x, y` | Simple transverse support |

Bending rotations remain free. A circular shaft has equal x- and y-bending
frequencies, so each physical bending frequency occurs as a pair. The test
selects one member from each pair.

For a 16-element mesh, every one of the first four modes must have relative
error below $0.1\%$:

$$
\varepsilon_m=
\frac{\left|f_{m,\mathrm{num}}-f_{m,\mathrm{ana}}\right|}
{f_{m,\mathrm{ana}}}.
$$

### `test_convergence.py`

The convergence test solves the same simply supported model with

```text
2, 4, 8, and 16 shaft elements
```

Only the discretization changes. Material, section, total length, formulation,
and supports remain fixed. The tests require:

- each mesh's first bending frequency to be within $1\%$ of the analytical
  value;
- the 16-element mesh to be more accurate than the 2-element mesh;
- the final relative error to be below $0.01\%$.

Strict monotonic improvement at every refinement step is deliberately not
required because small floating-point changes should not make a sound
convergence test brittle.

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
