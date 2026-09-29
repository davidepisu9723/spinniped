# Chapter 10 -- Plotting implementation

## Table of contents

- [Campbell plotting implementation](#campbell-plotting-implementation)

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
