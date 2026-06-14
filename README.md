# Spinniped

A minimal rotor/shaft finite-element toolbox for modal analysis.

## Overview

This small package provides simple classes and helpers to build rotor models, assemble global matrices, solve eigenproblems, and plot mode shapes.

Core modules:

- `spinniped.element` — shaft element stiffness/mass.
- `spinniped.node` — node representation with 3D coordinates.
- `spinniped.rotor` — rotor assembly, node/element management and plotting.
- `spinniped.solver` — global assembly, boundary conditions, eigen-solver.
- `spinniped.results` — plotting utilities for mode shapes.

## Requirements

- Python 3.8+
- NumPy
- SciPy
- Matplotlib

Install with pip:

```bash
python -m pip install numpy scipy matplotlib
```

If you prefer a single requirements file, create one with the packages above.

## Quick Usage

You can run the provided example script to try the library:

```bash
python examples/01_simple_shaft_modal.py
```

## Conventions

- Degrees of freedom (DOFs): the project uses 4 DOFs per node by default. Element DOF ordering typically follows `[x1, y1, tx1, ty1, x2, y2, tx2, ty2]`.
- Stiffness and mass matrix shapes are documented in the function docstrings in the source files.

## Examples and Benchmarks

- Example scripts: `examples/01_simple_shaft_modal.py` — a runnable demo of modal analysis.
- Benchmarks: `benchmark/01_beam_freq.py` contains analytical references for verification.

## Notes & To-Do

- The file `spinniped/results.py` contains an initializer named `__init_` (single underscore) which appears to be a typo; consider renaming it to `__init__` so instances initialize correctly.

## Contributing

Feel free to open pull requests or issues. For quick testing, run the example scripts and verify plots render.

---

This README is intentionally minimal. If you want, I can:

- add a `requirements.txt` or `pyproject.toml` file,
- add a short unit test that runs the example, or
- fix the `Results.__init_` typo and verify behavior.
