# Spinniped user manual

<p align="center">
  <img src="https://raw.githubusercontent.com/davidepisu9723/spinniped/master/logo/spinniped-logo.png" width="240" alt="Spinniped logo">
</p>

Spinniped is a small finite-element toolbox for deterministic and stochastic
rotor dynamics. Its input API follows a NASTRAN-like structure: grids,
coordinate systems, materials, properties, elements, and loads are immutable
Python dataclasses whose fields contain ordinary Python values. NumPy arrays
appear only after a definition is built into a numerical model.

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
3. [Speed-dependent values](03_speed_dependent_values.md)
4. [Random parameters](04_random_parameters.md)
5. [Built models](05_built_models.md)
6. [Analyses](06_analyses.md)
7. [Plotting](07_plotting.md)
8. [Examples and benchmarks](08_examples_and_benchmarks.md)
9. [Mapping analytical models to Spinniped](09_analytical_comparisons.md)
10. [Current limitations](10_current_limitations.md)

## Other manuals

- [Spinniped reference manual](../reference-manual/00_reference_manual.md)
- [Rotordynamics theory manual](../theory-manual/00_theory_manual.md)
- [Testing and verification](../../tests/README.md)
