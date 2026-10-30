"""Plot the stochastic unbalance response of a Jeffcott-style rotor."""

from matplotlib import pyplot as plt
import numpy as np

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
    Solver,
    Unbalance,
    plot_frequency_response,
)


# Twenty equal shaft elements place grid 11 exactly at the midpoint. One disk
# at that grid and identical end bearings form a symmetric Jeffcott-style
# rotor while retaining the distributed shaft mass and flexibility.
number_of_elements = 20
number_of_grids = number_of_elements + 1
shaft_length = 1.0
disk_grid = number_of_elements // 2 + 1

definition = ModelDefinition(
    grids=[
        Grid(
            id=grid_index + 1,
            z=grid_index * shaft_length / number_of_elements,
        )
        for grid_index in range(number_of_grids)
    ],
    materials=[
        Material(
            id=1,
            density=7_850.0,
            young_modulus=2.0e11,
            poisson_ratio=0.3,
        )
    ],
    properties=[
        ShaftProperty(
            id=1,
            material=1,
            outer_diameter=0.020,
            theory="timoshenko",
            rotary_inertia=True,
            damping=0.2,
        ),
        BearingProperty(
            id=2,
            # One reference is reused by both lateral coefficients and both
            # bearing elements, preserving symmetry within each realization.
            kxx=(1,),
            kyy=(1,),
            cxx=500.0,
            cyy=500.0,
        ),
        DiskProperty(
            id=3,
            mass=5.0,
            diametral_inertia=0.025,
            polar_inertia=0.050,
        ),
    ],
    elements=[
        *[
            ShaftElement(
                id=element_index + 1,
                grid_a=element_index + 1,
                grid_b=element_index + 2,
                property=1,
            )
            for element_index in range(number_of_elements)
        ],
        BearingElement(
            id=number_of_elements + 1,
            grid=1,
            property=2,
        ),
        BearingElement(
            id=number_of_elements + 2,
            grid=number_of_grids,
            property=2,
        ),
        DiskElement(
            id=number_of_elements + 3,
            grid=disk_grid,
            property=3,
        ),
    ],
    loads=[
        Unbalance(
            id=1,
            grid=disk_grid,
            mass_eccentricity=2.0e-4,
            phase=np.deg2rad(30.0),
        )
    ],
    distributions=[
        RandomDistribution(
            id=1,
            name="isotropic bearing stiffness",
            distribution="uniform",
            parameters={"low": 1.2e6, "high": 1.8e6},
        )
    ],
)

# The fixed seed makes the Monte Carlo ensemble reproducible. Bearing
# stiffness is sampled once per realization and shared by the two supports.
model = ModelBuilder().build(
    definition,
    stochastic=True,
    samples=20,
    seed=42,
)

# Unbalance response is synchronous: at every signed rotor speed, the solver
# uses abs(speed) as the excitation frequency and scales the rotating force by
# speed squared. Axial translation and torsional rotation are fixed once to
# remove the two uncoupled rigid-body coordinates.
maximum_speed_rpm = 6_000.0
speeds_rpm = np.linspace(0.0, maximum_speed_rpm, 61)
speeds = speeds_rpm * 2.0 * np.pi / 60.0
result = Solver(model).solve(
    "unbalance_response",
    fixed_dofs=[2, 5],
    speeds=speeds,
    load_ids=[1],
)

# Global x and y DOFs of the central disk are the first two coordinates of its
# grid. Plot their mean amplitudes, central empirical 90% bands, and the full
# sampled min--max envelope.
disk_x_dof = 6 * (disk_grid - 1)
disk_y_dof = disk_x_dof + 1
figure, axes = plot_frequency_response(
    result,
    dofs=[disk_x_dof, disk_y_dof],
    labels=["central-disk x", "central-disk y"],
    statistic="mean",
    confidence=0.90,
    show_extremes=True,
    title="Stochastic Jeffcott-rotor unbalance response",
)
axes.set_xlim(0.0, maximum_speed_rpm)
axes.set_ylabel("Central-disk displacement amplitude [m]")
figure.tight_layout()
plt.show()
