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
    plot_campbell,
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
# of rotor speed. The stiffness tables use a representative high-speed ball-
# bearing softening trend; actual trends depend on bearing geometry, preload,
# load direction, lubrication, and temperature. Requests outside the knot
# interval would be rejected.
kxx = SpeedDependentValue(
    speeds=speed_knots,
    values=(1.50e6, 1.28e6, 1.00e6),
)
kyy =  SpeedDependentValue(
    speeds=speed_knots,
    values=(1.55e6, 1.38e6, 0.80e6),
)
cxx = SpeedDependentValue(
    speeds=speed_knots,
    values=(400.0, 500.0, 650.0),
)
cyy = cxx

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
speeds_rpm = np.linspace(0.0, 6_000.0, 100)
speeds = speeds_rpm * rpm_to_rad_per_second
result = Solver(model).solve(
    "campbell",
    fixed_dofs=fixed_dofs,
    speeds=speeds,
    modes=2,
)

# A deterministic result still retains a leading sample dimension of one.
figure, axes = plot_campbell(
    result,
    title="Deterministic speed-dependent bearing",
)
axes.set_xlim(0.0, speeds_rpm[-1])
axes.set_ylim(bottom=0.0)
figure.tight_layout()
plt.show()
