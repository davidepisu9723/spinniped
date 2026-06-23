# Testing and Verification

Spinniped uses automated tests both to detect software regressions and to
verify that its finite-element formulation reproduces expected mechanical
behaviour.

The test suite is organized in four layers:

1. element-matrix tests;
2. global assembly and boundary-condition tests;
3. analytical benchmarks;
4. mesh-convergence tests.

## Requirements

Install Spinniped's normal dependencies and `pytest`:

```bash
python -m pip install numpy scipy matplotlib pytest
```

Run the complete suite from the project root:

```bash
python -m pytest
```

Useful alternatives are:

```bash
# Compact output
python -m pytest -q

# Run one test file
python -m pytest tests/test_element_matrices.py

# Run one individual test
python -m pytest tests/test_element_matrices.py::test_mass_matrix_is_positive_definite

# Stop after the first failure
python -m pytest -x

# Display output printed by tests
python -m pytest -s
```

## Test Structure

### Shared Fixtures and Helpers

`conftest.py` contains the standard shaft properties, reusable model-building
functions, and analytical frequency equations. Pytest loads this file
automatically.

The reference shaft is a uniform, filled circular shaft with:

| Quantity | Symbol | Value |
| --- | ---: | ---: |
| Young's modulus | $E$ |  $2.0 \times 10^{11}\ \mathrm{Pa} $ |
| Poisson's ratio | $\nu$ | $0.3$ |
| Density | $\rho$ | $7850\ \mathrm{kg/m^3}$ |
| Diameter | $d$ | $0.01\ \mathrm{m}$ |
| Total length | $L$ | $1.0\ \mathrm{m}$ |

The section properties are:

$$
A = \frac{\pi d^2}{4},
\qquad
I = \frac{\pi d^4}{64}.
$$

### Element-Matrix Tests

`test_element_matrices.py` checks the local shaft matrices before global
assembly can conceal an error. It verifies:

- matrix dimensions and finite entries;
- symmetry of the stiffness and mass matrices;
- six rigid-body modes in the unconstrained stiffness matrix;
- positive definiteness of the current mass formulation;
- closed-form axial and torsional stiffness blocks;
- linear scaling of mass with density.

These tests verify mathematical structure and basic physical consistency.
They do not prove every matrix entry independently.

### Assembly Tests

`test_assembly.py` uses deliberately small models whose global matrices can be
examined directly. It verifies:

- one-element assembly;
- accumulation of two elements at a shared node;
- absence of coupling between unconnected nodes;
- mapping of non-consecutive user node IDs;
- insertion of bearing stiffness into the intended $x,y$ DOFs;
- removal and deduplication of constrained DOFs;
- rejection of elements that reference unknown nodes.

The bearing test uses distinct stiffness values:

$$
K_b =
\begin{bmatrix}
11 & 12 \\
13 & 14
\end{bmatrix}.
$$

Distinct entries make transposition and DOF-permutation errors visible.

### Analytical Benchmark

`test_analytical_benchmarks.py` compares the first four bending frequencies of
a uniform simply supported shaft with the Euler-Bernoulli analytical solution:

$$
f_m =
\frac{m^2\pi}{2L^2}
\sqrt{\frac{EI}{\rho A}},
\qquad m = 1,2,\ldots
$$

The numerical benchmark uses Euler-Bernoulli stiffness and translational mass
without bending rotary inertia, matching the assumptions of the analytical
equation.

For the six-DOF nodal order

```text
x, y, z, tx, ty, tz
```

the boundary conditions are:

| End | Constrained local DOFs | Purpose |
| --- | --- | --- |
| First node | `x, y, z, tz` | Simple transverse support and removal of axial and torsional rigid motion |
| Final node | `x, y` | Simple transverse support |

Bending rotations remain free at both supports.

Because the shaft has a circular section, its $x$- and $y$-bending
frequencies are equal and occur in pairs. The benchmark selects one frequency
from each pair before comparing it with the analytical solution.

With 16 elements, the relative error of each of the first four bending modes
must remain below $0.1\%$:

$$
\varepsilon_m =
\frac{\left|f_{m,\mathrm{num}}-f_{m,\mathrm{ana}}\right|}
{f_{m,\mathrm{ana}}}.
$$

### Mesh Convergence

`test_convergence.py` solves the same simply supported shaft using:

```text
ne = 2, 4, 8, 16
```

Only mesh density changes. Material properties, section properties, boundary
conditions, and total length remain constant.

The tests require:

- every mesh to predict the first frequency within $1\%$;
- the 16-element result to be more accurate than the 2-element result;
- the final relative error to be below $0.01\%$.

The convergence test does not demand a strictly decreasing error at every
refinement step, avoiding unnecessary sensitivity to numerical noise.

## Numerical Tolerances

Tolerances should reflect the property being tested:

- tight tolerances are used for symmetry and exact matrix-block identities;
- rigid modes are detected relative to the largest stiffness eigenvalue;
- benchmark tolerances account for finite-element discretization error;
- convergence tests check both an improving trend and a final accuracy target.

A tolerance should not be loosened merely to make a failing test pass. First
determine whether the discrepancy comes from the implementation, benchmark
assumptions, mode selection, boundary conditions, or floating-point effects.

## Adding Tests

When adding a feature:

1. test its smallest local mathematical object;
2. test how it enters the global model;
3. compare it with an independent analytical or published reference;
4. test convergence when discretization is involved;
5. document the physical assumptions and tolerance.

Keep analytical formulas independent of the Spinniped implementation. Copying
the implementation into the test can reproduce the same error on both sides
and give false confidence.

## Current Scope

The present suite verifies the existing speed-independent structural-dynamics
foundation. It does not yet verify gyroscopic matrices, speed-dependent
eigenproblems, forward and backward whirl, critical speeds, or Campbell
diagrams. Those tests should be introduced together with the corresponding
rotordynamics features.
