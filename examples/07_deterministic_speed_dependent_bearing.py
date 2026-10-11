"""Plot a Campbell sweep with deterministic speed-dependent bearings."""

from matplotlib import pyplot as plt
import numpy as np

from spinniped import (
    BearingElement,
    BearingProperty,
    DiskElement,
    DiskProperty,
    Grid,
    ModelBuilder,
    ModelDefinition,
    Solver,
    SpeedDependentValue,
)


# SpeedDependentValue uses signed rotor angular speed in rad/s. Convert these
# readable source knots from rpm once when defining the tables.
speed_knots_rpm = (0.0, 3_000.0, 6_000.0)
rpm_to_rad_per_second = 2.0 * np.pi / 60.0
speed_knots = tuple(
    speed_rpm * rpm_to_rad_per_second
    for speed_rpm in speed_knots_rpm
)

# All four bearing coefficients are deterministic piecewise-linear functions
# of rotor speed. Requests outside the knot interval would be rejected.
kxx = SpeedDependentValue(
    speeds=speed_knots,
    values=(1.00e6, 1.20e6, 1.50e6),
)
kyy = SpeedDependentValue(
    speeds=speed_knots,
    values=(1.20e6, 1.35e6, 1.55e6),
)
cxx = SpeedDependentValue(
    speeds=speed_knots,
    values=(400.0, 500.0, 650.0),
)
cyy = SpeedDependentValue(
    speeds=speed_knots,
    values=(450.0, 560.0, 700.0),
)

# A one-grid Jeffcott model keeps the example focused on bearing interpolation.
# No RandomDistribution records are present, so the build has one sample.
definition = ModelDefinition(
    grids=[Grid(1)],
    properties=[
        BearingProperty(
            id=1,
            kxx=kxx,
            kyy=kyy,
            cxx=cxx,
            cyy=cyy,
        ),
        DiskProperty(
            id=2,
            mass=10.0,
            diametral_inertia=0.05,
            polar_inertia=0.10,
        ),
    ],
    elements=[
        BearingElement(id=1, grid=1, property=1),
        DiskElement(id=2, grid=1, property=2),
    ],
)

model = ModelBuilder().build(definition)

# Keep only the lateral translations. The solver validates and interpolates
# K(Omega) and C(Omega) over this complete grid before the first eigensolution.
fixed_dofs = [2, 3, 4, 5]
speeds_rpm = np.linspace(0.0, 6_000.0, 61)
speeds = speeds_rpm * rpm_to_rad_per_second
result = Solver(model).solve(
    "campbell",
    fixed_dofs=fixed_dofs,
    speeds=speeds,
    modes=2,
)

# A deterministic result still retains a leading sample dimension of one.
figure, axes = plt.subplots()
for mode, frequencies in enumerate(result["frequencies"][0].T, start=1):
    axes.plot(speeds_rpm, frequencies, label=f"mode {mode}")
for speed_rpm in speed_knots_rpm:
    axes.axvline(speed_rpm, color="black", linewidth=0.6, alpha=0.2)
axes.set_xlabel("Rotor speed [rpm]")
axes.set_ylabel("Natural frequency [Hz]")
axes.set_title("Deterministic speed-dependent bearing")
axes.grid(True, alpha=0.3)
axes.legend()
figure.tight_layout()
plt.show()
