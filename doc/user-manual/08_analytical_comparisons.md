# Chapter 8 -- Mapping analytical models to Spinniped

## Table of contents

- [Overview](#overview)

## Overview

Analytical and finite-element results are comparable only after their
assumptions match.

| Analytical choice | Spinniped representation |
|---|---|
| Shaft along $+z$ | `Grid` coordinates with increasing `z` |
| Uniform $E$, $\rho$, $\nu$ | one `Material` |
| Uniform circular section | one `ShaftProperty` |
| Euler--Bernoulli bending | `theory="euler"` |
| Neglect bending rotary inertia | `rotary_inertia=False` |
| Simple transverse supports | constrain `x`, `y` at both end grids |
| Remove axial rigid motion | also constrain `z` at one grid |
| Remove torsional rigid motion | also constrain `tz` at one grid |
| Center rigid disk | `DiskElement` at the midpoint grid |
| Equal isotropic bearings | equal `kxx` and `kyy` in two `BearingProperty` attachments |
| Nearly massless Jeffcott shaft | use a small positive shaft density and verify disk-dominated modes |
| Constant spin | pass `speed` or `speeds` in rad/s |

With Spinniped's nodal order

```text
x, y, z, tx, ty, tz
```

the stationary simply supported benchmark constrains `x`, `y`, `z`, and `tz`
at the first grid and `x`, `y` at the final grid. The axial and torsional
constraints remove rigid motion; bending rotations remain free.

The analytical beam benchmark uses a 16-element Euler--Bernoulli mesh and
requires the first four numerical bending frequencies to lie within 0.1% of
the closed-form values. The Jeffcott benchmark uses ten elements, a midpoint
disk, flexible end bearings, disabled shaft rotary inertia, and negligible
positive shaft density. Its modal and Campbell frequencies must agree with
the reduced formulas to relative tolerance $10^{-5}$.

These are verification comparisons, not proof that every real rotor follows
the idealizations. See [Testing and verification](../../tests/README.md) for
the exact scope and acceptance criteria.
