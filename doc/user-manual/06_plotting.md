# Chapter 6 -- Plotting

## Table of contents

- [Overview](#overview)

## Overview

Install the optional plotting dependency with:

```bash
python -m pip install -e ".[plot]"
```

`plot_campbell` accepts the result dictionary returned by the Campbell solver.
It plots one line per modal branch and defaults to rotor speed in rpm:

```python
from spinniped import plot_campbell

result = solver.solve(
    "campbell",
    fixed_dofs=fixed,
    speeds=[0.0, 100.0, 200.0],
    modes=8,
)

figure, axes = plot_campbell(result)
figure.savefig("campbell.png", dpi=150, bbox_inches="tight")
```

For a stochastic result, plot the mean or median branch together with an
empirical confidence band and the sampled minimum and maximum:

```python
figure, axes = plot_campbell(
    result,
    statistic="mean",       # alternatively "median"
    confidence=0.95,        # central 2.5th--97.5th percentile band
    show_extremes=True,
)
```

The default `statistic="sample"` plots one realization selected with
`sample=N`. Use `speed_unit="rad/s"` or `speed_unit="hz"` to change the
horizontal axis, and `ax=existing_axes` to draw into an existing layout. The
function returns the figure and axes and does not call `show()`, leaving
display and file output under caller control. Confidence bands are empirical
sample quantiles; they do not assume normally distributed frequencies.

Campbell modes are tracked by default. The solver uses displacement-vector MAC
to match each branch to the preceding speed and stores that choice in
`result["track_modes"]`. If the analysis is run with `track_modes=False`, the
plot still connects array columns, but those columns are independently sorted
at each speed and can exchange physical identity at a crossing.
