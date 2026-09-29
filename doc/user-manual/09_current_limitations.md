# Chapter 9 -- Current limitations

## Table of contents

- [Overview](#overview)

## Overview

- Global matrices are dense.
- Sampling currently uses plain Monte Carlo rather than variance-reduction
  methods such as Latin hypercube sampling.
- Correlation is supported within a multivariate normal distribution, but not
  between separate distribution records.
- Analyses return numerical dictionaries; plotting and result-object APIs are
  not currently provided.
