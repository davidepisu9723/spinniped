# Spinniped user manual

<p align="center">
  <img src="https://raw.githubusercontent.com/davidepisu9723/spinniped/master/logo/spinniped-logo.png" width="240" alt="Spinniped logo">
</p>

Spinniped is a small finite-element toolbox for deterministic and stochastic
rotor dynamics. Its input API follows a NASTRAN-like structure: grids,
coordinate systems, materials, properties, and elements are immutable Python
dataclasses whose fields contain ordinary Python values. NumPy arrays appear
only after a definition is built into a numerical model.

The project is under active development and is intended for preliminary design,
research, and education.

## Table of contents

- [Audience and scope](#audience-and-scope)
- [Chapters](#chapters)
- [Other manuals](#other-manuals)

## Audience and scope

This manual is for users who want to define, build, solve, and plot Spinniped
models through the public API. It describes software usage and result
interpretation; implementation details and analytical derivations are kept in
separate manuals.

## Chapters

1. [Installation](01_installation.md)
2. [Model-definition API](02_model_definition.md)
3. [Random parameters](03_random_parameters.md)
4. [Built models](04_built_models.md)
5. [Analyses](05_analyses.md)
6. [Plotting](06_plotting.md)
7. [Examples and benchmarks](07_examples_and_benchmarks.md)
8. [Mapping analytical models to Spinniped](08_analytical_comparisons.md)
9. [Current limitations](09_current_limitations.md)

## Other manuals

- [Spinniped reference manual](../reference-manual/00_reference_manual.md)
- [Rotordynamics theory manual](../theory-manual/00_theory_manual.md)
- [Testing and verification](../../tests/README.md)
