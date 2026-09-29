# Chapter 4 -- Properties

## Table of contents

- [Property record types](#property-record-types)
- [Shaft properties](#shaft-properties)
- [Bearing properties](#bearing-properties)
- [Disk properties](#disk-properties)

## Property record types

The `Property` union in `spinniped.records` contains `ShaftProperty`,
`BearingProperty`, and `DiskProperty`. `LumpedMassProperty` is an alias for
`DiskProperty`. All are frozen, slotted dataclasses and may be shared by
multiple compatible elements.

## Shaft properties

`ShaftProperty` stores `id`, the referenced material ID, outer and inner
diameters, the beam-theory selector, the rotary-inertia switch, and a
mass-proportional damping coefficient. The implemented theory values are
`"timoshenko"` and `"euler"`. The builder requires
`outer_diameter > inner_diameter >= 0`.

## Bearing properties

`BearingProperty` stores direct and cross-coupled translational stiffness and
viscous-damping coefficients for local $x$, $y$, and $z$. Coefficients are
interpreted in the coordinate system selected by the associated bearing
element. Bearing properties contribute no mass or gyroscopic terms.

## Disk properties

`DiskProperty` stores mass, equal diametral inertia, polar inertia, and a
mass-proportional damping coefficient. Mass and inertias must be finite and
nonnegative after resolution. Disk properties contribute inertia and
gyroscopic terms but no elastic stiffness.
