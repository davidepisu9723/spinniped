# Chapter 7 -- Examples and benchmarks

## Table of contents

- [Examples](#examples)
- [Analytical benchmarks](#analytical-benchmarks)

## Examples

Run the included examples from the repository root:

```bash
python examples/01_simple_shaft_modal.py
python examples/02_simple_shaft_modal_with_bearings.py
python examples/03_simple_shaft_complex_eigenvalues.py
python examples/04_bearing_supported_shaft_with_disk_campbell.py
```

These examples exercise normal public-API workflows rather than controlled
analytical validation cases.

## Analytical benchmarks

Analytical reference calculations are intentionally separate from normal-use
examples:

```bash
python benchmark/01_uniform_beam_frequencies.py
python benchmark/02_flexible_bearing_jeffcott_comparison.py
```

Benchmark 01 calculates analytical bending frequencies for simply supported
and fixed--fixed Euler--Bernoulli beams. Benchmark 02 compares the
speed-invariant cylindrical modes and gyroscopically
split conical modes of a flexible-bearing Jeffcott model with an equivalent
finite-element shaft. It also marks the cylindrical and backward-conical 1x
critical speeds.

The normative conventions and implemented shaft formulation are documented in
[Reference manual](../reference-manual/00_reference_manual.md) and
[Shaft-element implementation](../reference-manual/09_shaft_element_formulation.md). For a derivation-led
introduction to shaft and rotor dynamics, see
[Rotordynamics theory manual](../theory-manual/00_theory_manual.md).
