# Chapter 10 -- Plotting implementation

## Table of contents

- [Rotor-section plotting implementation](#rotor-section-plotting-implementation)
- [Campbell plotting implementation](#campbell-plotting-implementation)

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
