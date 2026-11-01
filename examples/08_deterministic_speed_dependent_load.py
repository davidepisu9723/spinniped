"""Plot deterministic response to a speed-dependent harmonic load."""

from matplotlib import pyplot as plt
import numpy as np

from spinniped import (
    BearingElement,
    BearingProperty,
    DiskElement,
    DiskProperty,
    Grid,
    HarmonicLoad,
    ModelBuilder,
    ModelDefinition,
    SinglePointConstraint,
    Solver,
    SpeedDependentValue,
    plot_frequency_response,
)


# Treat these amplitudes as deterministic CFD results at three rotor speeds.
# SpeedDependentValue requires tuple inputs and angular speed in rad/s,
# although rpm is retained here for the input-data labels and plot axis.
speed_knots_rpm = (0.0, 3_000.0, 6_000.0)
fluid_force_values = (50.0, 120.0, 200.0)
rpm_to_rad_per_second = 2.0 * np.pi / 60.0
speed_knots = tuple(
    speed_rpm * rpm_to_rad_per_second
    for speed_rpm in speed_knots_rpm
)
fluid_force_x = SpeedDependentValue(
    speeds=speed_knots,
    values=fluid_force_values,
)

# Keep the structure constant so the resulting change is caused only by the
# interpolated load amplitude.
definition = ModelDefinition(
    grids=[Grid(1)],
    properties=[
        BearingProperty(
            id=1,
            kxx=1.20e6,
            kyy=1.35e6,
            cxx=500.0,
            cyy=560.0,
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
    loads=[
        HarmonicLoad(
            id=1,
            grid=1,
            force=(fluid_force_x, 0.0, 0.0),
        )
    ],
    boundary_conditions=[
        SinglePointConstraint(1, dof=(2, 3, 4, 5))
    ],
)
model = ModelBuilder().build(definition)
solver = Solver(model)

# Frequency response accepts one scalar rotor speed per call. Solve one full
# forcing-frequency sweep at each CFD table knot to compare the response
# curves produced by the three tabulated load amplitudes.
forcing_frequencies = np.linspace(0.0, 80.0, 321)
figure, response_axes = plt.subplots()

for speed_rpm, speed, force_value in zip(
    speed_knots_rpm,
    speed_knots,
    fluid_force_values,
    strict=True,
):
    result = solver.solve(
        "frequency_response",
        frequencies=forcing_frequencies,
        speed=speed,
        load_ids=[1],
    )
    plot_frequency_response(
        result,
        dofs=[0],
        labels=[f"{speed_rpm:g} rpm, {force_value:g} N"],
        ax=response_axes,
        title=None,
        legend=False,
    )

response_axes.set_ylabel("x amplitude [m]")
response_axes.set_title("Speed-dependent harmonic-load response")
response_axes.legend()

figure.tight_layout()
plt.show()
