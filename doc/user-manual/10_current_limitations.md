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
- Rotor, Campbell, and steady-response plots are provided. Dedicated
  time-history plotting is not implemented.
- Declarative loads are currently limited to steady harmonic and unbalance
  response. Time response still requires an explicit vector or callback.
- Ordinary frequency response accepts one scalar rotor speed per call. A sweep
  of a speed-dependent harmonic or rotating load therefore requires one call
  per operating speed. `SynchronousLoad` selections support a vector of paired
  rotor speeds and order-derived forcing frequencies.
- Partial Campbell extraction searches near a fixed numerical target of
  $-1\,\mathrm{s}^{-1}$ on the real axis of the complex eigenvalue plane.
  Strong damping can change its selected roots relative to sorting the entire
  spectrum by frequency. Tracking is restricted to the extracted candidates;
  there is no adaptive expansion or failure fallback.
- Time response holds rotor speed fixed and uses RK45 with the integrator's
  default tolerances. The current public route does not expose alternative
  integration methods, step controls, or tolerance settings.
