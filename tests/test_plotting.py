"""Tests for visualization of solver result dictionaries."""

import matplotlib
import numpy as np
import pytest

from spinniped import (
    BearingElement,
    BearingProperty,
    DiskElement,
    DiskProperty,
    Grid,
    Material,
    ModelBuilder,
    ModelDefinition,
    RandomDistribution,
    ShaftElement,
    ShaftProperty,
    plot_campbell,
    plot_frequency_response,
    plot_rotor,
)


# Use a non-interactive backend so plotting tests never open a window.
matplotlib.use("Agg")


def _campbell_result():
    """Return a minimal two-sample Campbell result."""
    return {
        "analysis": "campbell",
        "speeds": np.array([0.0, np.pi, 2.0 * np.pi]),
        "frequencies": np.array(
            [
                [[10.0, 20.0], [11.0, 19.0], [12.0, 18.0]],
                [[10.5, 20.5], [11.5, 19.5], [12.5, 18.5]],
            ]
        ),
        "track_modes": True,
    }


def _critical_campbell_result():
    """Return tracked stochastic branches with known critical speeds."""
    speeds_hz = np.array([0.0, 8.0, 12.0, 20.0])
    frequencies = np.array(
        [
            [[10.0, 30.0]] * 4,
            [[12.0, 30.0]] * 4,
            [[8.0, 30.0]] * 4,
        ]
    )
    critical_speeds_hz = np.array(
        [
            [[[10.0], [np.nan]], [[5.0], [15.0]]],
            [[[12.0], [np.nan]], [[6.0], [15.0]]],
            [[[8.0], [np.nan]], [[4.0], [15.0]]],
        ]
    )
    return {
        "analysis": "campbell",
        "speeds": 2.0 * np.pi * speeds_hz,
        "frequencies": frequencies,
        "track_modes": True,
        "harmonics": np.array([1.0, 2.0]),
        "critical_speeds": 2.0 * np.pi * critical_speeds_hz,
    }


def _frequency_response_result(*, analysis="frequency_response", speeds=None):
    """Return a two-sample complex response with two nonadjacent free DOFs."""
    result = {
        "analysis": analysis,
        "frequencies": np.array([10.0, 20.0, 30.0]),
        "free_dofs": np.array([0, 2]),
        "response": np.array(
            [
                [[1.0, 4.0j], [2.0, 5.0j], [3.0, 6.0j]],
                [[3.0j, 6.0], [4.0j, 7.0], [5.0j, 8.0]],
            ],
            dtype=complex,
        ),
    }
    if speeds is not None:
        result["speeds"] = np.asarray(speeds, dtype=float)
    return result


def _rotor_model():
    """Return a stepped, partly hollow rotor with disks and bearings."""
    definition = ModelDefinition(
        grids=[Grid(1, z=0.0), Grid(2, z=0.4), Grid(3, z=1.0)],
        materials=[Material(1, 7850.0, 2.0e11, 0.3)],
        properties=[
            ShaftProperty(1, 1, outer_diameter=0.10, inner_diameter=0.04),
            ShaftProperty(2, 1, outer_diameter=0.06),
            BearingProperty(3, kxx=1.0e7, kyy=1.0e7),
            DiskProperty(4, mass=10.0, diametral_inertia=0.03,
                         polar_inertia=0.05),
        ],
        elements=[
            ShaftElement(1, 1, 2, 1),
            ShaftElement(2, 2, 3, 2),
            BearingElement(3, 1, 3),
            BearingElement(4, 3, 3),
            DiskElement(5, 2, 4),
        ],
    )
    return ModelBuilder().build(definition)


def test_plot_campbell_draws_one_line_per_tracked_mode():
    """Verify Campbell branches, sample selection, and rpm conversion."""
    import matplotlib.pyplot as plt

    figure, axes = plot_campbell(_campbell_result(), sample=1)

    assert figure is axes.figure
    assert len(axes.lines) == 2
    assert np.allclose(axes.lines[0].get_xdata(), [0.0, 30.0, 60.0])
    assert np.array_equal(axes.lines[0].get_ydata(), [10.5, 11.5, 12.5])
    assert np.array_equal(axes.lines[1].get_ydata(), [20.5, 19.5, 18.5])
    assert axes.get_xlabel() == "Speed [rpm]"
    assert axes.get_ylabel() == "Frequency [Hz]"
    assert [line.get_label() for line in axes.lines] == ["Mode 1", "Mode 2"]
    plt.close(figure)


def test_plot_campbell_draws_stochastic_statistics_and_bands():
    """Verify mean branches, confidence bands, and extrema rendering."""
    import matplotlib.pyplot as plt

    result = _campbell_result()
    result["frequencies"] = np.concatenate(
        (result["frequencies"], result["frequencies"][1:2] + 20.0),
        axis=0,
    )
    figure, axes = plot_campbell(
        result,
        statistic="mean",
        confidence=0.5,
        show_extremes=True,
        speed_unit="rad/s",
        legend=False,
    )

    expected_mean = np.mean(result["frequencies"], axis=0)
    assert len(axes.lines) == 6
    assert len(axes.collections) == 2
    assert np.array_equal(axes.lines[0].get_ydata(), expected_mean[:, 0])
    assert np.array_equal(axes.lines[3].get_ydata(), expected_mean[:, 1])
    assert axes.lines[0].get_label() == "Mode 1"
    assert axes.lines[1].get_label() == "_nolegend_"
    assert axes.lines[2].get_label() == "_nolegend_"
    assert axes.lines[1].get_linestyle() == "--"
    assert axes.lines[2].get_linestyle() == "--"
    assert np.array_equal(
        axes.lines[1].get_ydata(),
        np.min(result["frequencies"], axis=0)[:, 0],
    )
    assert np.array_equal(
        axes.lines[2].get_ydata(),
        np.max(result["frequencies"], axis=0)[:, 0],
    )
    assert axes.get_legend() is None
    assert [text.get_text() for text in axes.texts] == [
        "Displayed statistic: mean\n"
        "Confidence interval: 50%\n"
        "Mode extremes: sample min–max"
    ]
    plt.close(figure)

    median_figure, median_axes = plot_campbell(
        result,
        statistic="median",
        legend=False,
    )
    assert np.array_equal(
        median_axes.lines[0].get_ydata(),
        np.median(result["frequencies"], axis=0)[:, 0],
    )
    assert median_axes.lines[0].get_label() == "Mode 1"
    assert [text.get_text() for text in median_axes.texts] == [
        "Displayed statistic: median"
    ]
    plt.close(median_figure)


def test_plot_campbell_draws_deterministic_harmonics_and_critical_speeds():
    """Verify harmonic lines and one sample's critical-speed markers."""
    import matplotlib.pyplot as plt

    result = _critical_campbell_result()
    result["frequencies"] = result["frequencies"][:1]
    result["critical_speeds"] = result["critical_speeds"][:1]
    figure, axes = plot_campbell(
        result,
        show_harmonics=True,
        speed_unit="hz",
    )

    labels = [line.get_label() for line in axes.lines]
    assert "1x synchronous" in labels
    assert "2x synchronous" in labels
    assert "Mode 1 1x critical" in labels
    assert "Mode 1 2x critical" in labels
    assert "Mode 2 2x critical" in labels
    assert "Mode 2 1x critical" not in labels

    marker = next(
        line for line in axes.lines
        if line.get_gid() == "critical-h0-m0-c0"
    )
    assert marker.get_xdata()[0] == pytest.approx(10.0)
    assert marker.get_ydata()[0] == pytest.approx(10.0)
    plt.close(figure)


def test_plot_campbell_draws_stochastic_critical_speed_statistics():
    """Verify mean, standard deviation, and sample critical-speed artists."""
    import matplotlib.pyplot as plt

    figure, axes = plot_campbell(
        _critical_campbell_result(),
        statistic="mean",
        show_harmonics=True,
        show_critical_samples=True,
        speed_unit="hz",
        legend=False,
    )

    mean_marker = next(
        line for line in axes.lines
        if line.get_gid() == "critical-h0-m0-c0"
    )
    standard_deviation = next(
        patch for patch in axes.patches
        if patch.get_gid() == "critical-spread-h0-m0-c0"
    )
    assert mean_marker.get_xdata()[0] == pytest.approx(10.0)
    assert mean_marker.get_ydata()[0] == pytest.approx(10.0)
    assert np.allclose(
        standard_deviation._posA_posB,
        [(8.0, 8.0), (12.0, 12.0)],
    )
    assert len(axes.collections) == 3
    assert axes.get_legend() is None
    assert [text.get_text() for text in axes.texts] == [
        "Displayed statistic: mean\n"
        "Critical-speed error bars: ±1 std. dev."
    ]
    plt.close(figure)


def test_plot_frequency_response_selects_global_dofs_and_frequency_axis():
    """Verify ordinary response uses frequency and preserves DOF ordering."""
    import matplotlib.pyplot as plt

    result = _frequency_response_result()
    result["response"] = result["response"][:1]
    figure, axes = plot_frequency_response(
        result,
        dofs=[2, 0],
        labels=["transverse y", "transverse x"],
    )

    assert figure is axes.figure
    assert len(axes.lines) == 2
    assert np.array_equal(axes.lines[0].get_xdata(), [10.0, 20.0, 30.0])
    assert np.array_equal(axes.lines[0].get_ydata(), [4.0, 5.0, 6.0])
    assert np.array_equal(axes.lines[1].get_ydata(), [1.0, 2.0, 3.0])
    assert [line.get_label() for line in axes.lines] == [
        "transverse y",
        "transverse x",
    ]
    assert axes.get_xlabel() == "Excitation frequency [Hz]"
    assert axes.get_ylabel() == "Response amplitude"
    plt.close(figure)


@pytest.mark.parametrize("analysis", ["frequency_response", "unbalance_response"])
def test_plot_frequency_response_uses_speed_for_synchronous_results(analysis):
    """Verify synchronous result types automatically use their speed grid."""
    import matplotlib.pyplot as plt

    speeds = np.array([-2.0 * np.pi, 0.0, 4.0 * np.pi])
    result = _frequency_response_result(analysis=analysis, speeds=speeds)
    result["response"] = result["response"][:1]
    figure, axes = plot_frequency_response(
        result,
        dofs=[0],
        speed_unit="hz",
    )

    assert np.allclose(axes.lines[0].get_xdata(), [-1.0, 0.0, 2.0])
    assert axes.get_xlabel() == "Speed [Hz]"
    plt.close(figure)


def test_plot_frequency_response_draws_stochastic_statistics_and_bands():
    """Verify magnitude statistics, confidence bands, and extrema curves."""
    import matplotlib.pyplot as plt

    result = _frequency_response_result()
    figure, axes = plot_frequency_response(
        result,
        dofs=[0],
        statistic="mean",
        confidence=0.5,
        show_extremes=True,
        legend=False,
    )

    amplitudes = np.abs(result["response"][:, :, 0])
    assert len(axes.lines) == 3
    assert len(axes.collections) == 1
    assert np.array_equal(axes.lines[0].get_ydata(), np.mean(amplitudes, axis=0))
    assert np.array_equal(axes.lines[1].get_ydata(), np.min(amplitudes, axis=0))
    assert np.array_equal(axes.lines[2].get_ydata(), np.max(amplitudes, axis=0))
    assert axes.lines[1].get_linestyle() == "--"
    assert axes.lines[2].get_linestyle() == "--"
    assert axes.get_legend() is None
    assert [text.get_text() for text in axes.texts] == [
        "Displayed statistic: mean\n"
        "Confidence interval: 50%\n"
        "Response extremes: sample min\N{EN DASH}max"
    ]
    plt.close(figure)


def test_plot_frequency_response_validates_results_and_dof_selection():
    """Verify malformed result mappings and DOF selections fail early."""
    with pytest.raises(TypeError, match="result mapping"):
        plot_frequency_response([])

    wrong_analysis = _frequency_response_result()
    wrong_analysis["analysis"] = "campbell"
    with pytest.raises(ValueError, match="must come from"):
        plot_frequency_response(wrong_analysis)

    wrong_shape = _frequency_response_result()
    wrong_shape["response"] = wrong_shape["response"][:, :-1]
    with pytest.raises(ValueError, match="shape"):
        plot_frequency_response(wrong_shape)

    result = _frequency_response_result()
    with pytest.raises(ValueError, match="not present"):
        plot_frequency_response(result, dofs=[1])
    with pytest.raises(ValueError, match="duplicate"):
        plot_frequency_response(result, dofs=[0, 0])
    with pytest.raises(ValueError, match="one entry"):
        plot_frequency_response(result, dofs=[0, 2], labels=["x"])
    with pytest.raises(ValueError, match="one finite result speed"):
        plot_frequency_response(result, x_axis="speed")


@pytest.mark.parametrize(
    ("options", "exception", "message"),
    [
        ({"sample": 2}, IndexError, "out of range"),
        ({"sample": 0.5}, TypeError, "integer index"),
        ({"statistic": "average"}, ValueError, "statistic"),
        ({"confidence": 0.0}, ValueError, "between zero and one"),
        ({"confidence": True}, TypeError, "number between zero and one"),
        ({"show_extremes": 1}, TypeError, "show_extremes must be a boolean"),
        ({"x_axis": "time"}, ValueError, "x_axis"),
        ({"speed_unit": "mph"}, ValueError, "speed_unit"),
        ({"legend": 1}, TypeError, "legend must be a boolean"),
    ],
)
def test_plot_frequency_response_validates_plot_options(
    options, exception, message
):
    """Verify invalid frequency-response plot options have focused errors."""
    with pytest.raises(exception, match=message):
        plot_frequency_response(_frequency_response_result(), **options)


def test_plot_rotor_draws_exact_shaft_sections_disks_and_bearings():
    """Verify stepped and hollow shafts plus nodal component symbols."""
    import matplotlib.pyplot as plt

    figure, axes = plot_rotor(_rotor_model())
    patches = {
        patch.get_gid(): patch
        for patch in axes.patches
        if patch.get_gid() is not None
    }

    first_shaft = patches["shaft-1"]
    first_bore = patches["shaft-1-bore"]
    second_shaft = patches["shaft-2"]
    disk = patches["disk-5"]

    assert np.allclose(first_shaft.get_xy(), (0.0, -0.05))
    assert first_shaft.get_width() == pytest.approx(0.4)
    assert first_shaft.get_height() == pytest.approx(0.10)
    assert np.allclose(first_bore.get_xy(), (0.0, -0.02))
    assert first_bore.get_height() == pytest.approx(0.04)
    assert np.allclose(second_shaft.get_xy(), (0.4, -0.03))
    assert second_shaft.get_width() == pytest.approx(0.6)
    assert second_shaft.get_height() == pytest.approx(0.06)
    assert disk.get_x() + 0.5 * disk.get_width() == pytest.approx(0.4)
    assert disk.get_height() == pytest.approx(0.2)
    assert "bearing-3-upper" in patches
    assert "bearing-3-lower" in patches
    assert "bearing-4-upper" in patches
    assert "bearing-4-lower" in patches
    assert axes.get_xlabel() == "Axial position along spin axis [m]"
    assert axes.get_ylabel() == "Radial position [m]"
    assert set(text.get_text() for text in axes.get_legend().get_texts()) == {
        "Bearing", "Disk", "Shaft"
    }
    plt.close(figure)


def test_plot_rotor_selects_resolved_sample_and_reuses_axes():
    """Verify stochastic geometry selection, labels, and axes composition."""
    import matplotlib.pyplot as plt

    definition = ModelDefinition(
        grids=[Grid(1), Grid(2, z=1.0)],
        materials=[Material(1, 7850.0, 2.0e11, 0.3)],
        properties=[
            ShaftProperty(1, 1, outer_diameter=(1,)),
        ],
        elements=[ShaftElement(1, 1, 2, 1)],
        distributions=[
            RandomDistribution(
                1,
                "shaft diameter",
                "uniform",
                {"low": 0.08, "high": 0.12},
            )
        ],
    )
    model = ModelBuilder().build(
        definition, stochastic=True, samples=3, seed=7
    )
    figure, existing_axes = plt.subplots()
    existing_axes.set_title("Existing title")

    returned_figure, returned_axes = plot_rotor(
        model,
        sample=1,
        ax=existing_axes,
        title=None,
        legend=False,
        show_grid_ids=True,
    )

    resolved_diameter = model.definitions[1].properties[0].outer_diameter
    shaft_patch = next(
        patch for patch in returned_axes.patches
        if patch.get_gid() == "shaft-1"
    )
    assert returned_figure is figure
    assert returned_axes is existing_axes
    assert shaft_patch.get_height() == pytest.approx(resolved_diameter)
    assert returned_axes.get_title() == "Existing title"
    assert returned_axes.get_legend() is None
    assert [text.get_text() for text in returned_axes.texts] == ["1", "2"]
    plt.close(figure)


def test_plot_rotor_validates_model_sample_and_longitudinal_geometry():
    """Verify invalid inputs and nonparallel shaft geometry are rejected."""
    with pytest.raises(TypeError, match="BuiltModel"):
        plot_rotor(ModelDefinition())

    model = _rotor_model()
    with pytest.raises(TypeError, match="integer index"):
        plot_rotor(model, sample=True)
    with pytest.raises(IndexError, match="out of range"):
        plot_rotor(model, sample=1)

    angled_definition = ModelDefinition(
        grids=[Grid(1), Grid(2, x=0.1, z=1.0)],
        materials=[Material(1, 7850.0, 2.0e11, 0.3)],
        properties=[ShaftProperty(1, 1, outer_diameter=0.02)],
        elements=[ShaftElement(1, 1, 2, 1)],
    )
    angled_model = ModelBuilder().build(angled_definition)
    with pytest.raises(ValueError, match="parallel to spin_axis"):
        plot_rotor(angled_model)


@pytest.mark.parametrize(
    ("options", "exception", "message"),
    [
        ({"sample": 2}, IndexError, "out of range"),
        ({"sample": 0.5}, TypeError, "integer index"),
        ({"statistic": "average"}, ValueError, "statistic"),
        ({"confidence": 0.0}, ValueError, "between zero and one"),
        ({"confidence": True}, TypeError, "number between zero and one"),
        ({"show_extremes": 1}, TypeError, "show_extremes must be a boolean"),
        ({"speed_unit": "mph"}, ValueError, "speed_unit"),
        ({"legend": 1}, TypeError, "legend must be a boolean"),
        ({"show_harmonics": 1}, TypeError, "show_harmonics"),
        (
            {"show_critical_speeds": False, "show_critical_samples": True},
            ValueError,
            "requires show_critical_speeds",
        ),
    ],
)
def test_plot_campbell_validates_plot_options(options, exception, message):
    """Verify invalid plotting selections fail with focused messages."""
    with pytest.raises(exception, match=message):
        plot_campbell(_campbell_result(), **options)
