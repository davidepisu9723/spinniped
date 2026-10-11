# Chapter 10 -- Current limitations

## Table of contents

- [Overview](#overview)

## Overview

- Global matrices are dense.
- Sampling currently uses plain Monte Carlo rather than variance-reduction
  methods such as Latin hypercube sampling.
- Correlation is supported within a multivariate normal distribution, but not
  between separate distribution records.
- Speed dependence is limited to linear bearing stiffness and damping
  coefficients, harmonic force and moment components, and rotating-load
  magnitude. Geometry, material data, phases, and unbalance magnitude cannot
  be tabulated against speed.
- Speed-dependent tables use linear interpolation only and cannot extrapolate
  beyond their first and last knots.
- Analyses return numerical dictionaries rather than typed result objects.
- Rotor and Campbell plots are provided, but response plotting is not yet
  implemented.
- Declarative loads are currently limited to steady harmonic and unbalance
  response. Time response still requires an explicit vector or callback.
- Frequency response accepts one scalar rotor speed per call. A sweep of a
  speed-dependent harmonic or rotating load therefore requires one analysis
  call per rotor operating speed.
