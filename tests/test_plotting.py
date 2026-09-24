"""Tests for visualization of solver result dictionaries."""

import matplotlib
import numpy as np
import pytest

from spinniped import plot_campbell


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
    assert axes.lines[0].get_label() == "Mode 1 mean"
    assert axes.lines[1].get_label() == "Mode 1 minimum"
    assert axes.lines[2].get_label() == "Mode 1 maximum"
    assert np.array_equal(
        axes.lines[1].get_ydata(),
        np.min(result["frequencies"], axis=0)[:, 0],
    )
    assert np.array_equal(
        axes.lines[2].get_ydata(),
        np.max(result["frequencies"], axis=0)[:, 0],
    )
    assert axes.get_legend() is None
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
    assert median_axes.lines[0].get_label() == "Mode 1 median"
    plt.close(median_figure)


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
    ],
)
def test_plot_campbell_validates_plot_options(options, exception, message):
    """Verify invalid plotting selections fail with focused messages."""
    with pytest.raises(exception, match=message):
        plot_campbell(_campbell_result(), **options)
