# Chapter 10 -- Plotting implementation

## Table of contents

- [Rotor-section plotting implementation](#rotor-section-plotting-implementation)
- [Campbell plotting implementation](#campbell-plotting-implementation)
- [Critical-speed plotting](#critical-speed-plotting)

## Rotor-section plotting implementation

`spinniped.plotting.plot_rotor` consumes a `BuiltModel`, selects one resolved
definition and coordinate array through `sample`, and projects every grid onto
the normalized model spin axis. Before drawing, each shaft direction is split
into axial and transverse components. A transverse component larger than the
floating-point geometry tolerance raises an error rather than producing a
misleading section.

Each shaft element becomes an individually outlined rectangle with the exact
resolved outer diameter. An annular shaft receives a second,
background-colored rectangle with the exact inner diameter, preserving bores
and property steps. Artists receive stable Matplotlib GIDs such as
`shaft-4`, `shaft-4-bore`, `disk-8`, and `bearing-9-upper`.

Disk records specify mass moments rather than dimensions. When mass and polar
inertia are positive, the displayed radial extent includes the equivalent
solid-cylinder radius

$$
r_{\mathrm{eq}}=\sqrt{\frac{2I_p}{m}}.
$$

The disk's axial width and both bearing housing dimensions are conventional
symbols scaled from the plotted axial span and shaft radius. They must not be
read as manufacturing dimensions. Component layers are ordered bearing,
disk, shaft, bore, and centerline so shaft sections remain visible through
nodal symbols.

Matplotlib is imported lazily. The function validates the model type, sample
index, axes object, title, and Boolean display options before creating artists.
It returns the caller-owned or newly created figure and axes without displaying
them.

## Campbell plotting implementation

`spinniped.plotting.plot_campbell` consumes the mapping returned by the
Campbell solver. It creates one line per stored modal branch and converts the
speed axis to rpm, rad/s, or Hz.

For stochastic results, sample mode selects one realization. Mean and median
modes aggregate over the sample axis; confidence bands use empirical sample
quantiles, and extrema use sampled minima and maxima. The function returns its
figure and axes and does not call `show()`.

The implementation validates result shapes, sample indices, speed units,
statistics, confidence levels, and Boolean display options before creating the
plot.

## Critical-speed plotting

With `show_harmonics=True`, `plot_campbell` reads `harmonics` and the padded
four-dimensional `critical_speeds` array from the solver result. It does not
calculate or refine intersections. Harmonic ordinates and critical-point
ordinates are reconstructed from $r\Omega/(2\pi)$ after converting only the
displayed horizontal coordinate to the selected speed unit.

The private `_sample_statistics` helper reduces any leading sample axis and
returns the selected central value, mean, median, sample variance, standard
deviation, extrema, finite count, and optional empirical quantiles. The same
helper supplies Campbell frequency bands and stochastic critical-speed
summaries, keeping missing-crossing `NaN` values out of all reductions.

Sample statistics produce one marker per harmonic, mode, and crossing order.
Mean plots use a capped diagonal error bar spanning one standard deviation
when no confidence level is requested; otherwise it spans the empirical
confidence interval. Optional faint sample points show the underlying finite
realizations. Every error bar lies along its harmonic because the two plotted
critical-point coordinates are not independent random variables.

Confidence fills and diagonal error bars use Matplotlib's `_nolegend_` label.
Critical-point markers remain in the legend with their mode and harmonic,
alongside the central modal branches and synchronous lines. Statistical mode
and confidence metadata are rendered once in an axes-relative information box
instead of being repeated for every branch.

When extrema are enabled, both the sample-minimum and sample-maximum curves use
the same dashed style and their central branch color. Both carry the
`_nolegend_` label, while the information box records that the sample min--max
envelope is present.
