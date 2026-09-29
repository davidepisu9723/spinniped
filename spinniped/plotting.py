"""Plot Spinniped models and analysis results."""

from collections.abc import Mapping

import numpy as np

from .builder import BuiltModel
from .records import BearingElement, DiskElement, ShaftElement


_SPEED_CONVERSIONS = {
    "rpm": (lambda speeds: speeds * 60.0 / (2.0 * np.pi), "Speed [rpm]"),
    "rad/s": (lambda speeds: speeds, "Speed [rad/s]"),
    "hz": (lambda speeds: speeds / (2.0 * np.pi), "Speed [Hz]"),
}


def plot_rotor(
    model,
    *,
    sample=0,
    ax=None,
    title="Rotor longitudinal section",
    legend=True,
    show_grid_ids=False,
):
    """Plot the longitudinal section of a built rotor model.

    Parameters
    ----------
    model : BuiltModel
        Model returned by :class:`spinniped.ModelBuilder`.
    sample : int, optional
        Zero-based realization whose resolved geometry is plotted.
    ax : matplotlib.axes.Axes or None, optional
        Existing axes to draw on. When omitted, create a figure and axes.
    title : str or None, optional
        Axes title. Use ``None`` to leave an existing title unchanged.
    legend : bool, optional
        Draw one legend entry for shafts, disks, and bearings that are present.
    show_grid_ids : bool, optional
        Label every grid on the rotor centerline.

    Returns
    -------
    figure : matplotlib.figure.Figure
        Figure containing the rotor section.
    axes : matplotlib.axes.Axes
        Axes containing the section artists.

    Notes
    -----
    Shaft outer and inner diameters are drawn from the resolved
    :class:`~spinniped.ShaftProperty` values. A disk record contains inertial
    properties rather than manufacturing dimensions, so its radial size is
    represented by the equivalent solid-cylinder radius
    ``sqrt(2 * polar_inertia / mass)`` when available. Its visible axial width
    and the bearing housing dimensions are conventional plotting symbols.

    The plot is a section along ``ModelDefinition.spin_axis``. Shaft elements
    must therefore be parallel to that axis; rejecting nonparallel geometry
    avoids presenting an oblique or branched model as a misleading section.
    """
    if not isinstance(model, BuiltModel):
        raise TypeError("model must be a BuiltModel")
    if isinstance(sample, bool) or not isinstance(sample, (int, np.integer)):
        raise TypeError("sample must be an integer index")
    if sample < 0 or sample >= model.samples:
        raise IndexError("sample index is out of range")
    if title is not None and not isinstance(title, str):
        raise TypeError("title must be a string or None")
    if not isinstance(legend, bool):
        raise TypeError("legend must be a boolean")
    if not isinstance(show_grid_ids, bool):
        raise TypeError("show_grid_ids must be a boolean")

    definition = model.definitions[sample]
    coordinates = model.coordinates[sample]
    grid_index = {
        grid_id: index for index, grid_id in enumerate(model.grid_ids)
    }
    properties = {item.id: item for item in definition.properties}

    # The normalized spin direction defines the section's axial coordinate.
    spin_axis = np.asarray(definition.spin_axis, dtype=float)
    spin_axis /= np.linalg.norm(spin_axis)
    axial_positions = coordinates @ spin_axis

    shaft_sections = []
    disk_elements = []
    bearing_elements = []
    grid_radii = {grid_id: 0.0 for grid_id in model.grid_ids}

    for element in definition.elements:
        if isinstance(element, ShaftElement):
            property_record = properties[element.property]
            index_a = grid_index[element.grid_a]
            index_b = grid_index[element.grid_b]
            delta = coordinates[index_b] - coordinates[index_a]
            axial_delta = float(delta @ spin_axis)
            transverse_delta = delta - axial_delta * spin_axis
            tolerance = max(1.0e-12, 1.0e-9 * np.linalg.norm(delta))
            if np.linalg.norm(transverse_delta) > tolerance:
                raise ValueError(
                    "plot_rotor requires shaft elements parallel to spin_axis"
                )

            outer_radius = 0.5 * float(property_record.outer_diameter)
            inner_radius = 0.5 * float(property_record.inner_diameter)
            start = float(min(axial_positions[index_a], axial_positions[index_b]))
            end = float(max(axial_positions[index_a], axial_positions[index_b]))
            shaft_sections.append(
                (element, start, end, outer_radius, inner_radius)
            )
            grid_radii[element.grid_a] = max(
                grid_radii[element.grid_a], outer_radius
            )
            grid_radii[element.grid_b] = max(
                grid_radii[element.grid_b], outer_radius
            )
        elif isinstance(element, DiskElement):
            disk_elements.append(element)
        elif isinstance(element, BearingElement):
            bearing_elements.append(element)

    if not shaft_sections and not disk_elements and not bearing_elements:
        raise ValueError("model contains no shaft, disk, or bearing elements")

    # Establish data-unit sizes for symbols whose records have no geometry.
    plotted_positions = []
    for _, start, end, _, _ in shaft_sections:
        plotted_positions.extend((start, end))
    for element in disk_elements:
        plotted_positions.append(float(axial_positions[grid_index[element.grid]]))
    for element in bearing_elements:
        plotted_positions.append(float(axial_positions[grid_index[element.grid]]))
        if element.grid_b is not None:
            plotted_positions.append(
                float(axial_positions[grid_index[element.grid_b]])
            )

    axial_span = float(np.ptp(plotted_positions)) if plotted_positions else 0.0
    shaft_reference_radius = max(grid_radii.values(), default=0.0)
    if shaft_reference_radius <= 0.0:
        shaft_reference_radius = max(0.02 * axial_span, 1.0e-3)
    symbol_width = max(0.018 * axial_span, 0.35 * shaft_reference_radius)

    # Keep Matplotlib optional for users who only build and solve models.
    try:
        import matplotlib.pyplot as plt
        from matplotlib.patches import Rectangle
    except ImportError as error:
        raise ImportError(
            "plot_rotor requires Matplotlib; install Spinniped with "
            "the 'plot' extra"
        ) from error

    if ax is None:
        figure, axes = plt.subplots()
    else:
        if not hasattr(ax, "add_patch") or not hasattr(ax, "figure"):
            raise TypeError("ax must be a Matplotlib axes")
        axes = ax
        figure = axes.figure

    # Bearings sit behind the rotor. Their dimensions are schematic because
    # bearing records contain dynamic coefficients rather than housing sizes.
    bearing_label_available = True
    for element in bearing_elements:
        index_a = grid_index[element.grid]
        bearing_position = float(axial_positions[index_a])
        connected_radii = [grid_radii[element.grid]]
        if element.grid_b is not None:
            index_b = grid_index[element.grid_b]
            bearing_position = 0.5 * (
                bearing_position + float(axial_positions[index_b])
            )
            connected_radii.append(grid_radii[element.grid_b])
        local_radius = max(max(connected_radii), shaft_reference_radius * 0.5)
        inner_radius = 1.08 * local_radius
        outer_radius = max(1.8 * local_radius, 1.25 * shaft_reference_radius)
        label = "Bearing" if bearing_label_available else "_nolegend_"

        for side, lower_radius in (("upper", inner_radius), ("lower", -outer_radius)):
            bearing_patch = Rectangle(
                (bearing_position - 0.5 * symbol_width, lower_radius),
                symbol_width,
                outer_radius - inner_radius,
                facecolor="#d9a441",
                edgecolor="#6b4f16",
                hatch="///",
                linewidth=1.0,
                label=label if side == "upper" else "_nolegend_",
                zorder=1,
            )
            bearing_patch.set_gid(f"bearing-{element.id}-{side}")
            axes.add_patch(bearing_patch)
        bearing_label_available = False

    # Disk inertia defines an equivalent radial size. The axial width remains
    # a display convention because no disk-width field exists in the model.
    disk_label_available = True
    for element in disk_elements:
        property_record = properties[element.property]
        position = float(axial_positions[grid_index[element.grid]])
        local_radius = max(
            grid_radii[element.grid], shaft_reference_radius * 0.5
        )
        disk_radius = 1.75 * local_radius
        mass = float(property_record.mass)
        polar_inertia = float(property_record.polar_inertia)
        if mass > 0.0 and polar_inertia > 0.0:
            equivalent_radius = np.sqrt(2.0 * polar_inertia / mass)
            disk_radius = max(disk_radius, float(equivalent_radius))

        disk_patch = Rectangle(
            (position - 0.5 * symbol_width, -disk_radius),
            symbol_width,
            2.0 * disk_radius,
            facecolor="#6f7882",
            edgecolor="#252a2e",
            linewidth=1.2,
            label="Disk" if disk_label_available else "_nolegend_",
            zorder=2,
        )
        disk_patch.set_gid(f"disk-{element.id}")
        axes.add_patch(disk_patch)
        disk_label_available = False

    # Every shaft element gets its own outer outline, preserving steps between
    # adjacent properties. A second background-colored rectangle exposes the
    # exact bore of an annular section.
    shaft_label_available = True
    for element, start, end, outer_radius, inner_radius in shaft_sections:
        shaft_patch = Rectangle(
            (start, -outer_radius),
            end - start,
            2.0 * outer_radius,
            facecolor="#8db7d9",
            edgecolor="#24445f",
            linewidth=1.1,
            label="Shaft" if shaft_label_available else "_nolegend_",
            zorder=3,
        )
        shaft_patch.set_gid(f"shaft-{element.id}")
        axes.add_patch(shaft_patch)
        shaft_label_available = False

        if inner_radius > 0.0:
            bore_patch = Rectangle(
                (start, -inner_radius),
                end - start,
                2.0 * inner_radius,
                facecolor=axes.get_facecolor(),
                edgecolor="#24445f",
                linewidth=0.8,
                linestyle="--",
                label="_nolegend_",
                zorder=4,
            )
            bore_patch.set_gid(f"shaft-{element.id}-bore")
            axes.add_patch(bore_patch)

    axes.axhline(0.0, color="#444444", linewidth=0.7, linestyle="-.", zorder=5)
    if show_grid_ids:
        for grid_id, position in zip(
            model.grid_ids, axial_positions, strict=True
        ):
            axes.plot(
                position,
                0.0,
                marker="o",
                markersize=3.0,
                color="#222222",
                zorder=6,
            )
            axes.annotate(
                str(grid_id),
                (position, 0.0),
                xytext=(0.0, 6.0),
                textcoords="offset points",
                ha="center",
                va="bottom",
                fontsize="small",
            )

    axes.set_xlabel("Axial position along spin axis [m]")
    axes.set_ylabel("Radial position [m]")
    if title is not None:
        axes.set_title(title)
    axes.set_aspect("equal", adjustable="datalim")
    axes.autoscale_view()
    axes.margins(x=0.04, y=0.12)
    axes.grid(True, axis="x", alpha=0.25)
    if legend:
        axes.legend()

    return figure, axes


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
