# Chapter 9 -- Current limitations

## Table of contents

- [Overview](#overview)

## Overview

- Global matrices are dense.
- Sampling currently uses plain Monte Carlo rather than variance-reduction
  methods such as Latin hypercube sampling.
- Correlation is supported within a multivariate normal distribution, but not
  between separate distribution records.
- Analyses return numerical dictionaries rather than typed result objects.
- Rotor and Campbell plots are provided, but response plotting is not yet
  implemented.
- Declarative harmonic, rotating, and unbalance load records are resolved and
  validated, but response solvers do not consume them yet; frequency and time
  response still require explicit force vectors.
