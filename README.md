# Spinniped

<p align="center">
  <img src="logo/spinniped-logo.png" width="240" alt="Spinniped logo">
</p>

A minimal rotordynamics finite-element toolbox for preliminary design and academic purposes.

---

## Overview

This small package provides simple classes and helpers to build rotor models and solve dynamics/rotordynamics problems.

Core modules:

- `spinniped.element` — shaft element stiffness/mass.
- `spinniped.node` — node representation with 3D coordinates.
- `spinniped.rotor` — rotor assembly, node/element management and plotting.
- `spinniped.solver` — global assembly, boundary conditions, eigen-solver.
- `spinniped.results` — plotting utilities for mode shapes.

## Requirements
A user needs only the following libraries to run Spinniped.

- Python 3.13
- NumPy
- SciPy
- Matplotlib

Install with pip:

```bash
python -m pip install numpy scipy matplotlib
```

If you are a developer, before pushing you branch you should perform some checks to make sure your contribute is valid. Therefore, install `pytest`.

```bash
python -m pip install pytest
```

## Quick Usage

You can run the provided example script to try the library:

```bash
python examples/01_simple_shaft_modal.py
```

## Conventions

- Degrees of freedom (DOFs): the project uses 6 DOFs per node by default. Element DOF ordering typically follows `[x0, y0, z0, tx0, ty0, tz0, x1, y1, z1, tx1, ty1, tz1]`.
- Stiffness and mass matrix shapes are documented in the function docstrings in the source files.

## Examples and Benchmarks

- Example scripts: `examples/01_simple_shaft_modal.py` — a runnable demo of modal analysis.
- Benchmarks: `benchmark/01_beam_freq.py` contains analytical references for verification.

## Testing and Verification

Spinniped is checked through element-level tests, assembly tests, analytical
benchmarks, and mesh-convergence studies.

See [Testing and Verification](tests/README.md) for commands, assumptions,
boundary conditions, tolerances, and guidance for adding new tests.

## Notes & To-Do

The package is currently under construction, so it supports just a small subset of features. 

### Supported element types
| Type | Description | Application |
| --- | --- | --- |
| `ShaftElement` | Axialsymmetric two-nodes beam element with filled circular section. Timoshenko and Euler-Bernoulli formulations available | Standard element for modelling mechanical shafts |
| `BallBearing` | Bearing-to-ground support element with customable stiffness and damping matrices. Optional geometric parameters for plots. | Standard element for modelling ball bearing supports |

### To-Do

1. Implement rotordynamic equation solver. Therefore:
  - Define gyroscopic matrices for each element type (in `spinniped.element.py`)
  - Implement complex eignevalue solve sequence for rotordynamics problems (in `spinniped.solver.py`)
  - Implement plotting to visualize mode shapes (in `spinniped.results.py`)
  - Implement campbell diagrams (in `spinniped.results.py`)


## Contributing

Feel free to open pull requests or issues.

My goal is to keep the dependecies at minimal level. When adding features, do not add dependecies unless they are strictly necessary. Numpy, Scipy and Matplotlib are probably enuough for any purpose for which this package was created. Pytest is also necessary to make sure your contribution does not negatively effect the behavior of the library. 

If you implement a new element type, formulation, solver, or assembly
procedure, add appropriate tests and run:

```bash
python -m pytest
```

---
