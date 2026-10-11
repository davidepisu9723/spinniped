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
    Solver,
    SpeedDependentValue,
)


# Treat these amplitudes as deterministic CFD results at three rotor speeds.
# SpeedDependentValue requires angular speed in rad/s, although rpm is retained
# here for the input-data labels and plot axis.
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
)

model = ModelBuilder().build(definition)
solver = Solver(model)

# Frequency response accepts one scalar rotor speed per call. The forcing
# frequency stays fixed at 30 Hz, demonstrating that it is independent of the
# rotor speed used to interpolate the CFD load table.
fixed_dofs = [2, 3, 4, 5]
forcing_frequency_hz = 30.0
speeds_rpm = np.linspace(0.0, 6_000.0, 61)
speeds = speeds_rpm * rpm_to_rad_per_second
x_amplitude = []
for speed in speeds:
    result = solver.solve(
        "frequency_response",
        fixed_dofs=fixed_dofs,
        frequencies=[forcing_frequency_hz],
        speed=float(speed),
        load_ids=[1],
    )
    x_amplitude.append(abs(result["response"][0, 0, 0]))

figure, response_axes = plt.subplots()
response_line = response_axes.plot(
    speeds_rpm,
    x_amplitude,
    label="x response",
)[0]
response_axes.set_xlabel("Rotor speed [rpm]")
response_axes.set_ylabel("x amplitude [m]")
response_axes.set_title(
    f"Speed-dependent harmonic load at {forcing_frequency_hz:g} Hz"
)
response_axes.grid(True, alpha=0.3)

# Show the input CFD values on a second axis. Straight segments reflect the
# interpolation law used by the solver between the tabulated points.
load_axes = response_axes.twinx()
load_line = load_axes.plot(
    speed_knots_rpm,
    fluid_force_values,
    "o--",
    color="tab:orange",
    label="tabulated force",
)[0]
load_axes.set_ylabel("Force amplitude [N]")
response_axes.legend(handles=[response_line, load_line])

figure.tight_layout()
plt.show()
