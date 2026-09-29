# Chapter 3 -- Materials

## Table of contents

- [Record implementation](#record-implementation)
- [Validation and resolution](#validation-and-resolution)

## Record implementation

`Material` is a frozen, slotted dataclass in `spinniped.records`. It stores an
integer `id` and the numeric fields `density`, `young_modulus`, and
`poisson_ratio`. Each numeric field accepts a deterministic scalar or a
distribution reference before resolution.

## Validation and resolution

The builder requires unique positive integer material IDs. Resolved density and
Young's modulus must be finite and positive; Poisson's ratio must be finite and
strictly between $-1$ and $0.5$. A resolved material contains ordinary numeric
values before any element kernel is evaluated.

Shaft properties reference materials by ID. Bearings and disks do not currently
reference a material record.
