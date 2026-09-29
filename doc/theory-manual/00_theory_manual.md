# Rotordynamics theory manual

This guide develops the analytical models that sit behind Spinniped. It is
written for a reader who may be meeting rotor dynamics for the first time. The
aim is not to memorize formulas: it is to learn a repeatable route from a
physical shaft, through a free-body diagram and governing equation, to natural
frequencies, mode shapes, and a Campbell diagram.

The guide is independent of the finite-element implementation. When it later
compares an analytical result with Spinniped, the analytical equation remains
the reference rather than being derived from a program matrix.

## Table of contents

- [Purpose](#purpose)
- [Reading order](#reading-order)
- [Chapters](#chapters)
- [Other manuals](#other-manuals)

## Purpose

This manual develops the analytical mechanics behind the models without
describing the Spinniped API or implementation.

## Reading order

The chapters add one physical effect at a time:

1. establish coordinates, forces, and static equilibrium;
2. review the simplest mass--spring oscillator;
3. derive axial vibration of a stationary shaft;
4. derive torsional vibration of a stationary shaft;
5. derive flexural vibration of a stationary shaft;
6. add a rigid disk and distinguish diametral from polar inertia;
7. add spin and obtain the gyroscopically split Jeffcott frequencies;
8. calculate critical speeds and construct a Campbell diagram;
9. identify the assumptions and limits of each analytical model.

Read the chapters in order on a first pass. Sections headed **Assumptions** are
especially important: an analytical formula is useful only when the numerical
model represents the same idealization.

## Chapters

1. [Language, coordinates, and workflow](01_notation_and_workflow.md)
2. [Stationary simply supported shaft](02_static_simply_supported_shaft.md)
3. [Single-degree-of-freedom vibration](03_single_degree_of_freedom.md)
4. [Axial vibration](04_axial_vibration.md)
5. [Torsional vibration](05_torsional_vibration.md)
6. [Flexural vibration](06_flexural_vibration.md)
7. [Extended Jeffcott rotor](07_jeffcott_rotor.md)
8. [Critical speeds and Campbell diagrams](08_critical_speeds.md)
9. [Extensions and limits](09_extensions.md)

## Other manuals

- [Spinniped user manual](../user-manual/00_user_manual.md)
- [Spinniped reference manual](../reference-manual/00_reference_manual.md)
