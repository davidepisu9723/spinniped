"""Plot Spinniped models and analysis results."""

from collections.abc import Mapping

import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch, Rectangle
import numpy as np

from .builder import BuiltModel
from .records import BearingElement, DiskElement, ShaftElement


_SPEED_CONVERSIONS = {
    "rpm": (lambda speeds: speeds * 60.0 / (2.0 * np.pi), "Speed [rpm]"),
    "rad/s": (lambda speeds: speeds, "Speed [rad/s]"),
    "hz": (lambda speeds: speeds / (2.0 * np.pi), "Speed [Hz]"),
}


def _finite_sample_reduction(values, reducer):
    """Reduce a leading sample axis while ignoring missing realizations."""
    values = np.asarray(values, dtype=float)
    flattened = values.reshape(values.shape[0], -1)
    reduced = np.full(flattened.shape[1], np.nan)
    counts = np.zeros(flattened.shape[1], dtype=int)

    for column_index, column in enumerate(flattened.T):
        finite_values = column[np.isfinite(column)]
        counts[column_index] = len(finite_values)
        if finite_values.size:
            reduced[column_index] = reducer(finite_values)

    output_shape = values.shape[1:]
    return reduced.reshape(output_shape), counts.reshape(output_shape)


def _sample_statistics(values, *, statistic, sample, confidence):
    """Return consistent sample statistics for curves or crossing arrays."""
    mean, counts = _finite_sample_reduction(values, np.mean)
    median, _ = _finite_sample_reduction(values, np.median)
    minimum, _ = _finite_sample_reduction(values, np.min)
    maximum, _ = _finite_sample_reduction(values, np.max)
    variance, _ = _finite_sample_reduction(
        values,
        lambda data: np.var(data, ddof=1) if len(data) > 1 else 0.0,
    )

    if statistic == "sample":
        central = values[sample]
    elif statistic == "mean":
        central = mean
    else:
        central = median

    lower = upper = None
    if confidence is not None:
        tail_probability = (1.0 - confidence) / 2.0
        lower, _ = _finite_sample_reduction(
            values,
            lambda data: np.quantile(data, tail_probability),
        )
        upper, _ = _finite_sample_reduction(
            values,
            lambda data: np.quantile(data, 1.0 - tail_probability),
        )

    return {
        "central": central,
        "mean": mean,
        "median": median,
        "variance": variance,
        "standard_deviation": np.sqrt(variance),
        "minimum": minimum,
        "maximum": maximum,
        "lower": lower,
        "upper": upper,
        "count": counts,
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
    show_harmonics=False,
    show_critical_speeds=True,
    show_critical_samples=False,
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
    show_harmonics : bool, optional
        Plot the synchronous harmonic indices requested during the Campbell
        solution. For index ``h``, the plotted frequency is
        ``h * speed / (2*pi)``. The result must contain critical-speed data.
    show_critical_speeds : bool, optional
        Mark critical speeds when synchronous harmonics are shown. Sample mode
        marks the selected realization. Mean mode also draws capped diagonal
        one-standard-deviation error bars unless a confidence interval was
        requested.
    show_critical_samples : bool, optional
        Draw every finite sample critical speed as a faint point. This is most
        useful with a stochastic mean or median Campbell diagram.

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

    Matplotlib is a core Spinniped dependency and is imported with this module.
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
    if not isinstance(show_harmonics, bool):
        raise TypeError("show_harmonics must be a boolean")
    if not isinstance(show_critical_speeds, bool):
        raise TypeError("show_critical_speeds must be a boolean")
    if not isinstance(show_critical_samples, bool):
        raise TypeError("show_critical_samples must be a boolean")
    if show_critical_samples and not show_critical_speeds:
        raise ValueError(
            "show_critical_samples requires show_critical_speeds=True"
        )

    harmonic_indices = np.empty(0)
    critical_speeds = None
    if show_harmonics:
        if result.get("track_modes") is not True:
            raise ValueError(
                "critical-speed plotting requires tracked Campbell modes"
            )
        harmonic_indices = np.asarray(result.get("harmonics"), dtype=float)
        critical_speeds = np.asarray(
            result.get("critical_speeds"), dtype=float
        )
        if (
            harmonic_indices.ndim != 1
            or not harmonic_indices.size
            or not np.isfinite(harmonic_indices).all()
            or np.any(harmonic_indices <= 0.0)
        ):
            raise ValueError(
                "result must contain positive Campbell harmonic indices"
            )
        expected_prefix = (
            frequencies.shape[0],
            len(harmonic_indices),
            frequencies.shape[2],
        )
        if (
            critical_speeds.ndim != 4
            or critical_speeds.shape[:3] != expected_prefix
            or np.isinf(critical_speeds).any()
        ):
            raise ValueError(
                "result critical_speeds must have shape "
                "(samples, harmonics, modes, crossings)"
            )

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

    # One statistics helper supplies identical sample semantics to modal curves
    # and critical-speed points, including finite-only stochastic reductions.
    frequency_statistics = _sample_statistics(
        frequencies,
        statistic=statistic,
        sample=sample,
        confidence=confidence,
    )
    central_frequencies = frequency_statistics["central"]

    # Plot every tracked mode separately. Confidence fills and extrema reuse
    # the central line's color so all representations of a mode stay grouped.
    mode_colors = []
    for mode_index in range(frequencies.shape[2]):
        label = f"Mode {mode_index + 1}"
        line = axes.plot(
            plotted_speeds,
            central_frequencies[:, mode_index],
            label=label,
        )[0]
        color = line.get_color()
        mode_colors.append(color)

        if confidence is not None:
            axes.fill_between(
                plotted_speeds,
                frequency_statistics["lower"][:, mode_index],
                frequency_statistics["upper"][:, mode_index],
                color=color,
                alpha=0.2,
                label="_nolegend_",
            )

        if show_extremes:
            axes.plot(
                plotted_speeds,
                frequency_statistics["minimum"][:, mode_index],
                color=color,
                linestyle="--",
                linewidth=1.0,
                label="_nolegend_",
            )
            axes.plot(
                plotted_speeds,
                frequency_statistics["maximum"][:, mode_index],
                color=color,
                linestyle="--",
                linewidth=1.0,
                label="_nolegend_",
            )

    if show_harmonics:
        for harmonic_position, h in enumerate(harmonic_indices):
            axes.plot(
                plotted_speeds,
                h * speeds / (2.0 * np.pi),
                color="#303030",
                linestyle=(0, (5, 2 + harmonic_position % 3)),
                linewidth=1.0,
                label=f"{h:g}x synchronous",
                zorder=1,
            )

        if show_critical_speeds:
            critical_statistics = _sample_statistics(
                critical_speeds,
                statistic=statistic,
                sample=sample,
                confidence=confidence,
            )

            for harmonic_position, h in enumerate(harmonic_indices):
                for mode_index, color in enumerate(mode_colors):
                    for crossing_index in range(critical_speeds.shape[3]):
                        critical_speed = critical_statistics["central"][
                            harmonic_position, mode_index, crossing_index
                        ]
                        if not np.isfinite(critical_speed):
                            continue

                        crossing_label = (
                            f"Mode {mode_index + 1} {h:g}x critical"
                        )
                        if critical_speeds.shape[3] > 1:
                            crossing_label += f" #{crossing_index + 1}"

                        critical_artist = axes.plot(
                            convert_speed(critical_speed),
                            h * critical_speed / (2.0 * np.pi),
                            marker="o",
                            linestyle="none",
                            markersize=6.0,
                            markerfacecolor=color,
                            markeredgecolor="#202020",
                            label=crossing_label,
                            zorder=6,
                        )[0]
                        critical_artist.set_gid(
                            f"critical-h{harmonic_position}-m{mode_index}-"
                            f"c{crossing_index}"
                        )

                        lower_speed = upper_speed = None
                        if confidence is not None:
                            lower_speed = critical_statistics["lower"][
                                harmonic_position, mode_index, crossing_index
                            ]
                            upper_speed = critical_statistics["upper"][
                                harmonic_position, mode_index, crossing_index
                            ]
                        elif statistic == "mean":
                            deviation = critical_statistics[
                                "standard_deviation"
                            ][harmonic_position, mode_index, crossing_index]
                            lower_speed = critical_speed - deviation
                            upper_speed = critical_speed + deviation

                        if (
                            lower_speed is not None
                            and np.isfinite(lower_speed)
                            and np.isfinite(upper_speed)
                            and upper_speed > lower_speed
                        ):
                            spread_speeds = np.array(
                                [lower_speed, upper_speed]
                            )
                            spread_x = convert_speed(spread_speeds)
                            spread_y = (
                                h * spread_speeds / (2.0 * np.pi)
                            )
                            spread_artist = FancyArrowPatch(
                                (spread_x[0], spread_y[0]),
                                (spread_x[1], spread_y[1]),
                                arrowstyle="|-|",
                                mutation_scale=8.0,
                                color=color,
                                linewidth=1.3,
                                alpha=0.9,
                                label="_nolegend_",
                                zorder=4,
                            )
                            spread_artist.set_gid(
                                f"critical-spread-h{harmonic_position}-"
                                f"m{mode_index}-c{crossing_index}"
                            )
                            axes.add_patch(spread_artist)

                        if show_critical_samples:
                            sample_speeds = critical_speeds[
                                :, harmonic_position, mode_index, crossing_index
                            ]
                            sample_speeds = sample_speeds[
                                np.isfinite(sample_speeds)
                            ]
                            axes.scatter(
                                convert_speed(sample_speeds),
                                h * sample_speeds / (2.0 * np.pi),
                                color=color,
                                s=12.0,
                                alpha=0.3,
                                edgecolors="none",
                                label="_nolegend_",
                                zorder=5,
                            )

    # Apply consistent engineering labels after all artists have been added.
    axes.set_xlabel(speed_label)
    axes.set_ylabel("Frequency [Hz]")
    if title is not None:
        axes.set_title(title)
    axes.grid(True)

    # Keep statistical metadata out of an already component-heavy legend.
    # Deterministic sample plots need no box because they contain no ensemble
    # reduction or uncertainty interval.
    if frequencies.shape[0] > 1 or confidence is not None or show_extremes:
        if statistic == "sample":
            statistic_text = f"sample {sample}"
        else:
            statistic_text = statistic
        information = [f"Displayed statistic: {statistic_text}"]
        if confidence is not None:
            information.append(f"Confidence interval: {100.0 * confidence:g}%")
        elif show_harmonics and show_critical_speeds and statistic == "mean":
            information.append("Critical-speed error bars: ±1 std. dev.")
        if show_extremes:
            information.append("Mode extremes: sample min–max")
        axes.text(
            0.02,
            0.98,
            "\n".join(information),
            transform=axes.transAxes,
            ha="left",
            va="top",
            fontsize="small",
            bbox={
                "boxstyle": "round,pad=0.35",
                "facecolor": "white",
                "edgecolor": "#777777",
                "alpha": 0.85,
            },
            zorder=10,
        )
    if legend:
        axes.legend(loc="upper right")

    # Returning both objects lets callers further customize, save, or display
    # the plot without this utility imposing an output workflow.
    return figure, axes
