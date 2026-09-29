# Chapter 6 -- Plotting

## Table of contents

- [Installation](#installation)
- [Rotor longitudinal section](#rotor-longitudinal-section)
- [Campbell diagrams](#campbell-diagrams)
- [Synchronous harmonics and critical speeds](#synchronous-harmonics-and-critical-speeds)

## Installation

Matplotlib is a core Spinniped dependency and is installed automatically with
the package. Plotting functions are therefore available without an optional
installation extra.

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

## Synchronous harmonics and critical speeds

Request the desired ratios during the Campbell solution, then enable their
display in the plot:

```python
result = solver.solve(
    "campbell",
    fixed_dofs=fixed,
    speeds=speeds,
    modes=8,
    harmonics=[1.0, 2.0],
)

figure, axes = plot_campbell(
    result,
    show_harmonics=True,
)
```

For a deterministic result or `statistic="sample"`, each finite critical speed
of the selected realization is marked on its harmonic line. For a stochastic
mean plot, the marker is the mean of the sample critical speeds and a diagonal
error bar along the harmonic represents one sample standard deviation. The
error bar has perpendicular end caps. Supplying
`confidence=0.95` makes this error bar span the empirical central 95% interval.
A median plot uses the median marker and shows an error bar only when a
confidence level is requested.

The legend identifies modal branches, synchronous harmonics, and critical-speed
markers. Confidence bands and critical-speed error bars are intentionally
excluded from it. A separate information box reports whether the displayed
central values are a sample, mean, or median and gives the confidence level
when one is active.

With `show_extremes=True`, the sample minimum and maximum of each modal branch
are drawn using the same dashed linestyle and branch color. They are omitted
from the legend; the information box instead reports that the sample min--max
extremes are displayed.

Set `show_critical_samples=True` to add every finite realization as a faint
point. Samples without that crossing are omitted rather than converted to
zero. Because frequency at a critical point is fixed by the harmonic ratio,
uncertainty lies along the synchronous line; the plot therefore uses a line
error bar rather than a two-dimensional ellipse.

`plot_campbell` only visualizes the crossings stored by the solver. Neither the
solver nor the plotting function refines the speed discretization, so selecting
adequately spaced speed samples remains the user's responsibility. Set
`show_critical_speeds=False` to draw harmonic lines without their intersections.
