"""Plot analysis results produced by :class:`spinniped.Solver`."""

from collections.abc import Mapping

import numpy as np


_SPEED_CONVERSIONS = {
    "rpm": (lambda speeds: speeds * 60.0 / (2.0 * np.pi), "Speed [rpm]"),
    "rad/s": (lambda speeds: speeds, "Speed [rad/s]"),
    "hz": (lambda speeds: speeds / (2.0 * np.pi), "Speed [Hz]"),
}


def plot_campbell(
    result,
    *,
    sample=0,
    statistic="sample",
    confidence=None,
    show_extremes=False,
    speed_unit="rpm",
    ax=None,
    title="Campbell diagram",
    legend=True,
):
    """Plot modal-frequency branches from a Campbell analysis.

    Parameters
    ----------
    result : mapping
        Dictionary returned by ``Solver.solve("campbell", ...)``. Its
        ``frequencies`` array must have shape ``(samples, speeds, modes)``.
    sample : int, optional
        Zero-based realization index used when ``statistic="sample"``.
        Deterministic results use zero.
    statistic : {"sample", "mean", "median"}, optional
        Central curve plotted for every modal branch. ``"sample"`` plots one
        realization; ``"mean"`` and ``"median"`` aggregate all realizations.
    confidence : float or None, optional
        Confidence level for a central empirical quantile band across
        realizations. For example, ``0.95`` shades the 2.5th--97.5th
        percentile interval. ``None`` disables confidence bands.
    show_extremes : bool, optional
        Plot the sample minimum and maximum of every modal branch.
    speed_unit : {"rpm", "rad/s", "hz"}, optional
        Unit used for the horizontal rotor-speed axis.
    ax : matplotlib.axes.Axes or None, optional
        Existing axes to draw on. When omitted, create a new figure and axes.
    title : str or None, optional
        Axes title. Use ``None`` to leave an existing title unchanged.
    legend : bool, optional
        Draw a legend containing one label per modal branch.

    Returns
    -------
    figure : matplotlib.figure.Figure
        Figure containing the diagram.
    axes : matplotlib.axes.Axes
        Axes containing the modal branches.

    Notes
    -----
    Campbell analysis tracks modal branches across speeds and stochastic
    samples by default. Statistics, extremes, and confidence bands are thus
    evaluated mode-by-mode. When the analysis was run with
    ``track_modes=False``, columns instead follow independent frequency order
    and may exchange physical identity at crossings.

    Matplotlib is imported only when this function is called, so model
    construction and solution do not require the plotting dependency.
    """
    # Restrict the input to solver-like mappings before accessing result keys.
    if not isinstance(result, Mapping):
        raise TypeError("result must be a Campbell result mapping")
    if result.get("analysis") != "campbell":
        raise ValueError("result must come from a Campbell analysis")

    # Normalize array-like values and enforce the solver's
    # (samples, speeds, modes) Campbell result convention.
    speeds = np.asarray(result.get("speeds"), dtype=float)
    frequencies = np.asarray(result.get("frequencies"), dtype=float)
    if speeds.ndim != 1 or not speeds.size or not np.isfinite(speeds).all():
        raise ValueError("result speeds must be a finite one-dimensional array")
    if (
        frequencies.ndim != 3
        or frequencies.shape[1] != speeds.size
        or not np.isfinite(frequencies).all()
    ):
        raise ValueError(
            "result frequencies must have shape (samples, speeds, modes)"
        )

    # Validate plot selections explicitly instead of relying on NumPy indexing
    # or Matplotlib to produce less-focused errors later.
    if (
        not isinstance(statistic, str)
        or statistic not in {"sample", "mean", "median"}
    ):
        raise ValueError("statistic must be 'sample', 'mean', or 'median'")
    if statistic == "sample":
        if isinstance(sample, bool) or not isinstance(
            sample, (int, np.integer)
        ):
            raise TypeError("sample must be an integer index")
        if sample < 0 or sample >= frequencies.shape[0]:
            raise IndexError("sample index is out of range")
    if confidence is not None:
        if isinstance(confidence, bool):
            raise TypeError("confidence must be a number between zero and one")
        try:
            confidence = float(confidence)
        except (TypeError, ValueError):
            raise TypeError(
                "confidence must be a number between zero and one"
            ) from None
        if not np.isfinite(confidence) or not 0.0 < confidence < 1.0:
            raise ValueError("confidence must be between zero and one")
    if not isinstance(show_extremes, bool):
        raise TypeError("show_extremes must be a boolean")
    if speed_unit not in _SPEED_CONVERSIONS:
        choices = ", ".join(repr(unit) for unit in _SPEED_CONVERSIONS)
        raise ValueError(f"speed_unit must be one of {choices}")
    if title is not None and not isinstance(title, str):
        raise TypeError("title must be a string or None")
    if not isinstance(legend, bool):
        raise TypeError("legend must be a boolean")

    # Keep Matplotlib optional for users who only build and solve models.
    try:
        import matplotlib.pyplot as plt
    except ImportError as error:
        raise ImportError(
            "plot_campbell requires Matplotlib; install Spinniped with "
            "the 'plot' extra"
        ) from error

    # Reuse caller-owned axes when composing a multi-panel figure; otherwise
    # create an independent figure suitable for direct display or saving.
    if ax is None:
        figure, axes = plt.subplots()
    else:
        if not hasattr(ax, "plot") or not hasattr(ax, "figure"):
            raise TypeError("ax must be a Matplotlib axes")
        axes = ax
        figure = axes.figure

    # Solver speeds are always radians per second. Convert only the displayed
    # horizontal coordinate, leaving the result dictionary unchanged.
    convert_speed, speed_label = _SPEED_CONVERSIONS[speed_unit]
    plotted_speeds = convert_speed(speeds)

    # Collapse only the sample axis. Speed and tracked-mode axes are preserved
    # so every output column remains one Campbell branch.
    if statistic == "sample":
        central_frequencies = frequencies[sample]
        central_label = None
    elif statistic == "mean":
        central_frequencies = np.mean(frequencies, axis=0)
        central_label = "mean"
    else:
        central_frequencies = np.median(frequencies, axis=0)
        central_label = "median"

    # Extrema are inexpensive to compute and are shared by every modal branch.
    minimum_frequencies = np.min(frequencies, axis=0)
    maximum_frequencies = np.max(frequencies, axis=0)
    if confidence is not None:
        # Split the omitted probability equally between both tails. This forms
        # an empirical central band without assuming a normal distribution.
        tail_probability = (1.0 - confidence) / 2.0
        lower_frequencies, upper_frequencies = np.quantile(
            frequencies,
            [tail_probability, 1.0 - tail_probability],
            axis=0,
        )

    # Plot every tracked mode separately. Confidence fills and extrema reuse
    # the central line's color so all representations of a mode stay grouped.
    for mode_index in range(frequencies.shape[2]):
        label = f"Mode {mode_index + 1}"
        if central_label is not None:
            label += f" {central_label}"
        line = axes.plot(
            plotted_speeds,
            central_frequencies[:, mode_index],
            label=label,
        )[0]
        color = line.get_color()

        if confidence is not None:
            confidence_percent = 100.0 * confidence
            axes.fill_between(
                plotted_speeds,
                lower_frequencies[:, mode_index],
                upper_frequencies[:, mode_index],
                color=color,
                alpha=0.2,
                label=f"Mode {mode_index + 1} {confidence_percent:g}% band",
            )

        if show_extremes:
            axes.plot(
                plotted_speeds,
                minimum_frequencies[:, mode_index],
                color=color,
                linestyle="--",
                linewidth=1.0,
                label=f"Mode {mode_index + 1} minimum",
            )
            axes.plot(
                plotted_speeds,
                maximum_frequencies[:, mode_index],
                color=color,
                linestyle=":",
                linewidth=1.0,
                label=f"Mode {mode_index + 1} maximum",
            )

    # Apply consistent engineering labels after all artists have been added.
    axes.set_xlabel(speed_label)
    axes.set_ylabel("Frequency [Hz]")
    if title is not None:
        axes.set_title(title)
    axes.grid(True)
    if legend:
        axes.legend()

    # Returning both objects lets callers further customize, save, or display
    # the plot without this utility imposing an output workflow.
    return figure, axes
