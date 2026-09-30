"""Generate the deterministic figures embedded in the three manuals.

Run from the repository root with::

    python doc/generate_manual_figures.py

The generated shaft-DOF and whirl illustrations are non-deterministic assets
and are therefore not produced by this script.
"""

from pathlib import Path
import sys

import matplotlib

matplotlib.use("Agg")

from matplotlib import pyplot as plt
from matplotlib.lines import Line2D
import numpy as np


ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from spinniped import (  # noqa: E402
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
    Solver,
    plot_campbell,
    plot_rotor,
)


USER_IMAGES = ROOT / "doc" / "user-manual" / "images"
REFERENCE_IMAGES = ROOT / "doc" / "reference-manual" / "images"
THEORY_IMAGES = ROOT / "doc" / "theory-manual" / "images"

BLUE = "#2878B5"
ORANGE = "#E07A2D"
GREEN = "#4C956C"
PURPLE = "#8064A2"
GOLD = "#D6A23A"
CHARCOAL = "#27313A"
LIGHT = "#F4F7F9"


def _configure_style():
    plt.rcParams.update(
        {
            "figure.dpi": 120,
            "savefig.dpi": 180,
            "font.size": 10,
            "axes.titlesize": 12,
            "axes.labelsize": 10,
            "legend.fontsize": 9,
            "axes.spines.top": False,
            "axes.spines.right": False,
        }
    )


def _save(figure, path):
    path.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(path, bbox_inches="tight", facecolor="white")
    plt.close(figure)


def draw_mode_tracking():
    speed = np.linspace(0.0, 1.0, 13)
    physical_a = 16.0 + 15.0 * speed
    physical_b = 31.0 - 15.0 * speed
    before_crossing = speed <= 0.5
    after_crossing = speed >= 0.5

    figure, (sorting, tracking) = plt.subplots(1, 2, figsize=(12, 4.7), sharey=True)
    # Color identifies the independently sorted array column. Marker geometry
    # identifies the physical shape. The marker change at the crossing
    # therefore exposes the loss of physical identity directly.
    sorting.plot(
        speed[before_crossing],
        physical_a[before_crossing],
        "o-",
        color=BLUE,
        linewidth=2.2,
    )
    sorting.plot(
        speed[after_crossing],
        physical_b[after_crossing],
        "s-",
        color=BLUE,
        linewidth=2.2,
    )
    sorting.plot(
        speed[before_crossing],
        physical_b[before_crossing],
        "s-",
        color=ORANGE,
        linewidth=2.2,
    )
    sorting.plot(
        speed[after_crossing],
        physical_a[after_crossing],
        "o-",
        color=ORANGE,
        linewidth=2.2,
    )
    sorting.axvline(0.5, color="#9AA3AA", linestyle=":", linewidth=1.2)
    sorting.set_title("Sorting by frequency alone")
    sorting.set_xlabel("normalized rotor speed")
    sorting.set_ylabel("frequency [arbitrary units]")
    sorting.legend(
        handles=[
            Line2D([0], [0], color=BLUE, linewidth=2.2, label="array column 0"),
            Line2D([0], [0], color=ORANGE, linewidth=2.2, label="array column 1"),
            Line2D([0], [0], color="#737C85", marker="o", linestyle="None", label="physical shape A"),
            Line2D([0], [0], color="#737C85", marker="s", linestyle="None", label="physical shape B"),
        ],
        loc="upper center",
        ncol=2,
        frameon=False,
        columnspacing=1.0,
    )
    sorting.annotate(
        "marker shapes swap columns",
        xy=(0.5, 23.5),
        xytext=(0.5, 20.0),
        ha="center",
        color=CHARCOAL,
        arrowprops={"arrowstyle": "->", "color": CHARCOAL, "linewidth": 1.1},
    )
    sorting.text(0.13, 18.8, "A / column 0", color=BLUE, fontsize=9)
    sorting.text(0.13, 28.8, "B / column 1", color=ORANGE, fontsize=9)
    sorting.text(0.72, 18.8, "B / column 0", color=BLUE, fontsize=9)
    sorting.text(0.72, 28.8, "A / column 1", color=ORANGE, fontsize=9)

    tracking.plot(speed, physical_a, "o-", color=BLUE, label="tracked shape A")
    tracking.plot(speed, physical_b, "s-", color=ORANGE, label="tracked shape B")
    tracking.set_title("MAC-based correspondence")
    tracking.set_xlabel("normalized rotor speed")
    tracking.legend(loc="upper center", frameon=False)
    tracking.text(0.52, 22.2, "mode identities pass\nthrough the crossing", ha="center", color=CHARCOAL)

    for axes in (sorting, tracking):
        axes.grid(alpha=0.22)
        axes.set_xlim(0, 1)
    figure.suptitle("Why Campbell branches need mode-shape tracking", x=0.06, ha="left", weight="bold")
    figure.tight_layout(rect=(0, 0, 1, 0.92))
    _save(figure, REFERENCE_IMAGES / "mode_tracking_at_crossing.png")


def draw_separation_of_variables():
    z = np.linspace(0.0, 1.0, 300)
    shape = np.sin(0.5 * np.pi * z)
    time = np.linspace(0.0, 2.0, 400)
    amplitude = np.cos(2.0 * np.pi * time)

    figure, axes = plt.subplots(1, 3, figsize=(15, 4.8))
    spatial, temporal, snapshots = axes

    spatial.plot(z, shape, color=BLUE, linewidth=2.5)
    spatial.scatter([0, 1], [0, 1], color=[CHARCOAL, BLUE], zorder=3)
    spatial.set_title(r"Spatial harmonic $U(z)$")
    spatial.set_xlabel(r"position $z/L$")
    spatial.set_ylabel("relative displacement")
    spatial.text(0.03, 0.08, "fixed end", transform=spatial.transAxes, color="#52616B")
    spatial.text(0.70, 0.82, "free end", transform=spatial.transAxes, color="#52616B")

    temporal.plot(time, amplitude, color=ORANGE, linewidth=2.2)
    temporal.axhline(0, color="#A4ABB1", linewidth=0.8)
    temporal.set_title(r"Temporal harmonic $a(t)$")
    temporal.set_xlabel(r"time $t/T$")
    temporal.set_ylabel("modal amplitude")

    phases = [0.0, 0.125, 0.25, 0.5]
    styles = [(BLUE, "-"), (GREEN, "--"), (GOLD, "-."), (PURPLE, ":")]
    for phase, (color, linestyle) in zip(phases, styles):
        scale = np.cos(2.0 * np.pi * phase)
        snapshots.plot(z, scale * shape, color=color, linestyle=linestyle, linewidth=2.0, label=rf"$t/T={phase:g}$")
    snapshots.axhline(0, color="#A4ABB1", linewidth=0.8)
    snapshots.set_title(r"Product $u(z,t)=U(z)a(t)$")
    snapshots.set_xlabel(r"position $z/L$")
    snapshots.set_ylabel("instantaneous displacement")
    snapshots.legend(frameon=False, loc="lower left")

    for axes_item in axes:
        axes_item.grid(alpha=0.2)
    figure.suptitle(
        "A normal mode keeps one spatial shape while its amplitude oscillates in time",
        x=0.04,
        ha="left",
        weight="bold",
    )
    figure.tight_layout(rect=(0, 0, 1, 0.91))
    _save(figure, THEORY_IMAGES / "separation_of_variables.png")


def _rotor_definition(*, stochastic):
    number_of_elements = 10
    number_of_grids = number_of_elements + 1
    disk_grid = number_of_grids // 2 + 1
    bearing_kxx = (1,) if stochastic else 1.5e5
    distributions = []
    if stochastic:
        distributions.append(
            RandomDistribution(
                id=1,
                name="horizontal bearing stiffness",
                distribution="uniform",
                parameters={"low": 1.0e5, "high": 2.0e5},
            )
        )

    return ModelDefinition(
        grids=[Grid(id=index + 1, z=index / number_of_elements) for index in range(number_of_grids)],
        materials=[Material(id=1, density=7850.0, young_modulus=2.0e11, poisson_ratio=0.3)],
        properties=[
            ShaftProperty(
                id=1,
                material=1,
                outer_diameter=0.020,
                theory="timoshenko",
                rotary_inertia=True,
                damping=0.2,
            ),
            BearingProperty(id=2, kxx=bearing_kxx, kyy=3.0e5, cxx=500.0, cyy=500.0),
            DiskProperty(id=3, mass=5.0, diametral_inertia=0.025, polar_inertia=0.05),
        ],
        elements=[
            *[
                ShaftElement(id=index + 1, grid_a=index + 1, grid_b=index + 2, property=1)
                for index in range(number_of_elements)
            ],
            BearingElement(id=11, grid=1, property=2),
            BearingElement(id=12, grid=number_of_grids, property=2),
            DiskElement(id=13, grid=disk_grid, property=3),
        ],
        distributions=distributions,
    )


def draw_spinniped_outputs():
    deterministic = ModelBuilder().build(_rotor_definition(stochastic=False))
    rotor_figure, rotor_axes = plt.subplots(figsize=(11, 3.6))
    plot_rotor(
        deterministic,
        ax=rotor_axes,
        show_grid_ids=True,
        title="Resolved rotor geometry: bearings, shaft elements, and central disk",
    )
    rotor_figure.tight_layout()
    _save(rotor_figure, USER_IMAGES / "rotor_longitudinal_section.png")

    stochastic = ModelBuilder().build(
        _rotor_definition(stochastic=True),
        stochastic=True,
        samples=40,
        seed=42,
    )
    speeds_rpm = np.linspace(0.0, 3000.0, 25)
    result = Solver(stochastic).solve(
        "campbell",
        fixed_dofs=[2, 5],
        speeds=speeds_rpm * 2.0 * np.pi / 60.0,
        modes=4,
        harmonics=[1.0],
    )
    campbell_figure, campbell_axes = plt.subplots(figsize=(10.2, 6.1))
    plot_campbell(
        result,
        ax=campbell_axes,
        statistic="mean",
        confidence=0.90,
        show_harmonics=True,
        show_critical_samples=True,
        show_extremes=True,
        title="Stochastic Campbell diagram: uncertain horizontal bearing stiffness",
    )
    campbell_axes.set_xlim(0.0, 3000.0)
    campbell_axes.set_ylim(bottom=0.0)
    campbell_figure.tight_layout()
    _save(campbell_figure, USER_IMAGES / "stochastic_campbell.png")


def main():
    _configure_style()
    draw_mode_tracking()
    draw_separation_of_variables()
    draw_spinniped_outputs()
    print("Generated manual figures in doc/*-manual/images")


if __name__ == "__main__":
    main()
