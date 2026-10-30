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
`"timoshenko"` and `"euler"`. Record construction requires a positive ID and
material reference, `outer_diameter > inner_diameter >= 0`, nonnegative
damping, a supported theory, and a Boolean rotary-inertia switch. Random
parameters are checked again when their resolved record is reconstructed;
property-ID uniqueness and material lookup remain builder responsibilities.

## Bearing properties

`BearingProperty` stores direct and cross-coupled translational stiffness and
viscous-damping coefficients for local $x$, $y$, and $z$. Coefficients are
interpreted in the coordinate system selected by the associated bearing
element. Bearing properties contribute no mass or gyroscopic terms.

Each stiffness field `kxx`, `kyy`, `kzz`, `kxy`, and `kyx` and damping field
`cxx`, `cyy`, `czz`, `cxy`, and `cyx` accepts either an ordinary parameter or
a `SpeedDependentValue`. A tabulated coefficient is omitted from the constant
base matrix and represented by a transformed unit-coefficient contribution.
The solver later scales and adds that contribution at the requested signed
rotor speed. Shaft and disk property fields do not accept speed tables.

Record construction validates the positive integer bearing-property ID and
every finite constant, random-reference, or tabulated coefficient. Coefficients
remain signed because cross-coupled terms and some identified bearing models
may legitimately use either sign. Property-ID uniqueness and element/property
compatibility remain builder responsibilities.

## Disk properties

`DiskProperty` stores mass, equal diametral inertia, polar inertia, and a
mass-proportional damping coefficient. Mass and inertias must be finite and
nonnegative after resolution. Disk properties contribute inertia and
gyroscopic terms but no elastic stiffness.
