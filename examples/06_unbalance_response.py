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
    SinglePointConstraint,
    Solver,
    Unbalance,
    plot_frequency_response,
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
    boundary_conditions=[
        SinglePointConstraint(1, dof=(2, 3, 4, 5))
    ],
)
model = ModelBuilder().build(definition)

# Only the two lateral translations participate in this reduced Jeffcott
# calculation. Rotor speeds are supplied in rad/s; forcing is synchronous.
speeds_rpm = np.linspace(0.0, 6_000.0, 301)
speeds = speeds_rpm * 2.0 * np.pi / 60.0
result = Solver(model).solve(
    "unbalance_response",
    speeds=speeds,
    load_ids=[1],
)

figure, axes = plot_frequency_response(
    result,
    dofs=[0, 1],
    labels=["x amplitude", "y amplitude"],
    title="Synchronous unbalance response",
)
axes.set_ylabel("Displacement amplitude [m]")
figure.tight_layout()
plt.show()
