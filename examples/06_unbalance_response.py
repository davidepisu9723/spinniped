"""Plot the synchronous response of a simple unbalanced Jeffcott rotor."""

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
    Unbalance,
)


# A single disk and grounded bearing pair form the translational Jeffcott
# model. The direct coefficients differ slightly to expose both response axes.
definition = ModelDefinition(
    grids=[Grid(1)],
    properties=[
        BearingProperty(
            id=1,
            kxx=1.0e6,
            kyy=1.2e6,
            cxx=400.0,
            cyy=400.0,
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
        Unbalance(
            id=1,
            grid=1,
            mass_eccentricity=2.0e-4,
            phase=np.deg2rad(30.0),
        )
    ],
)

model = ModelBuilder().build(definition)

# Only the two lateral translations participate in this reduced Jeffcott
# calculation. Rotor speeds are supplied in rad/s; forcing is synchronous.
speeds_rpm = np.linspace(0.0, 6_000.0, 301)
speeds = speeds_rpm * 2.0 * np.pi / 60.0
result = Solver(model).solve(
    "unbalance_response",
    fixed_dofs=[2, 3, 4, 5],
    speeds=speeds,
    load_ids=[1],
)

x_amplitude = np.abs(result["response"][0, :, 0])
y_amplitude = np.abs(result["response"][0, :, 1])

figure, axes = plt.subplots()
axes.plot(speeds_rpm, x_amplitude, label="x amplitude")
axes.plot(speeds_rpm, y_amplitude, label="y amplitude")
axes.set_xlabel("Rotor speed [rpm]")
axes.set_ylabel("Displacement amplitude [m]")
axes.set_title("Synchronous unbalance response")
axes.grid(True, alpha=0.3)
axes.legend()
figure.tight_layout()
plt.show()
