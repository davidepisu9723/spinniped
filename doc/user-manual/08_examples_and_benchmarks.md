# Chapter 8 -- Examples and benchmarks

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
python examples/05_stochastic_campbell.py
python examples/06_unbalance_response.py
python examples/07_deterministic_speed_dependent_bearing.py
python examples/08_deterministic_speed_dependent_load.py
python examples/09_stochastic_speed_dependent_campbell.py
python examples/10_synchronous_load_response.py
```

These examples exercise normal public-API workflows rather than controlled
analytical validation cases. Example 04 uses `plot_rotor` to display the
model's longitudinal section before plotting its Campbell diagram, 1x
synchronous line, and speed-interpolated critical speeds.

Example 05 builds a reproducible 40-realization ensemble whose only uncertain
input is the direct stiffness $k_{xx}$, with the same realization used at both
supports while $k_{yy}$ remains deterministic. It plots mean tracked branches,
central empirical 90% bands, the 1x synchronous line, mean critical speeds,
critical-speed intervals, and the individual realization crossings.

Example 06 defines a declarative disk unbalance on an anisotropic Jeffcott
rotor, sweeps synchronous rotor speed, and plots the complex response magnitude
in both lateral directions.

Example 07 builds a deterministic Jeffcott model with tabulated bearing
stiffness and damping. Its Campbell sweep evaluates K and C across the complete
rotor-speed grid before solving.

Example 08 keeps the structure constant and interpolates a deterministic
harmonic-load amplitude from three representative CFD operating points. It
calls frequency response once per scalar rotor speed while keeping the forcing
frequency fixed, making the distinction between the two frequencies explicit.

Example 09 combines stochastic and speed-dependent bearing stiffness in a
symmetric Jeffcott-style rotor with ten shaft elements, one central disk, and
identical bearings at both ends. A multivariate normal distribution correlates
the stiffness ordinates across rotor speed. It plots mean Campbell branches,
empirical confidence bands, the 1x synchronous line, and the sampled critical-
speed crossings. Its sweep begins at 1 rpm because the symmetric rotor has
repeated forward/backward modes at exactly zero speed; Chapter 6 explains why
that degeneracy can make pairwise stochastic mode tracking mix the separated
whirl branches.

Example 10 defines a 2x `SynchronousLoad` using the ordinary harmonic force
vector convention. It passes the complete rotor-speed vector to one frequency-
response call without a `frequencies` argument, then plots both displacement
amplitude and the order-derived forcing frequency.

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
[Shaft-element implementation](../reference-manual/11_shaft_element_formulation.md). For a derivation-led
introduction to shaft and rotor dynamics, see
[Rotordynamics theory manual](../theory-manual/00_theory_manual.md).
