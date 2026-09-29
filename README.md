# Spinniped

<p align="center">
  <img src="https://raw.githubusercontent.com/davidepisu9723/spinniped/master/logo/spinniped-logo.png" width="240" alt="Spinniped logo">
</p>

Spinniped is a small finite-element toolbox for deterministic and stochastic
rotor dynamics. It is intended for preliminary design, research, and
education.

## Table of contents

- [Installation](#installation)
- [Minimal example](#minimal-example)
- [Documentation](#documentation)
- [Tests](#tests)

## Installation

Spinniped requires Python 3.13 or newer, NumPy, and SciPy. Until a release is
available on PyPI, install a local clone in editable mode:

```bash
git clone https://github.com/davidepisu9723/spinniped.git
cd spinniped
python -m pip install -e .
```

Optional plotting and test dependencies are installed with:

```bash
python -m pip install -e ".[plot,test]"
```

## Minimal example

```python
from spinniped import Grid, Material, ModelBuilder, ModelDefinition
from spinniped import ShaftElement, ShaftProperty, Solver

definition = ModelDefinition(
    grids=[Grid(1, z=0.0), Grid(2, z=1.0)],
    materials=[Material(1, density=7850.0, young_modulus=210e9,
                        poisson_ratio=0.3)],
    properties=[ShaftProperty(1, material=1, outer_diameter=0.02)],
    elements=[ShaftElement(1, grid_a=1, grid_b=2, property=1)],
)

model = ModelBuilder().build(definition)
result = Solver(model).solve("modal", fixed_dofs=[0, 1, 2, 5], modes=4)
print(result["frequencies"][0])
```

## Documentation

- [Spinniped user manual](doc/user-manual/00_user_manual.md) explains how to
  define, build, solve, and plot models through the public API.
- [Spinniped reference manual](doc/reference-manual/00_reference_manual.md)
  documents records, validation, matrix kernels, assembly, solvers, and other
  implementation details.
- [Rotordynamics theory manual](doc/theory-manual/00_theory_manual.md) develops
  the analytical mechanics from free-body diagrams through shaft vibration and
  the extended Jeffcott rotor.
- [Examples and benchmarks](doc/user-manual/07_examples_and_benchmarks.md)
  distinguishes normal API workflows from controlled reference calculations
  and finite-element comparisons.
- [Testing and verification](tests/README.md) explains how those formulations
  are checked.

## Tests

```bash
python -m pytest
```

Spinniped is released under the MIT License.
