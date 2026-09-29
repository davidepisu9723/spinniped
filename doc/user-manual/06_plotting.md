# Chapter 6 -- Plotting

## Table of contents

- [Installation](#installation)
- [Rotor longitudinal section](#rotor-longitudinal-section)
- [Campbell diagrams](#campbell-diagrams)

## Installation

Install the optional plotting dependency with:

```bash
python -m pip install -e ".[plot]"
```

## Rotor longitudinal section

`plot_rotor` draws the resolved longitudinal section of a `BuiltModel`. Shaft
outer and inner diameters are plotted in data units, so steps, solid sections,
and hollow sections follow the properties assigned to each shaft element.
Disk and bearing records do not contain manufacturing dimensions: disks are
therefore shown using an equivalent radius inferred from mass and polar
inertia, while disk width and bearing housings use conventional plot symbols.

```python
from spinniped import plot_rotor

figure, axes = plot_rotor(model, show_grid_ids=True)
figure.savefig("rotor.png", dpi=150, bbox_inches="tight")
```

For a stochastic model, select the resolved realization with `sample=N`. Use
`ax=existing_axes` to add the section to an existing figure, `legend=False` to
hide the component legend, or `title=None` to preserve an existing axes title.
The function returns the figure and axes and does not call `show()`.

The longitudinal coordinate is the projection onto the model's `spin_axis`.
Every shaft must be parallel to that direction. A nonparallel shaft is rejected
because its true geometry cannot be represented by one longitudinal section.

## Campbell diagrams

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
