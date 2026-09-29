"""Plot a stochastic Campbell diagram and uncertain critical speeds."""

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
    plot_campbell,
)


# A short mesh keeps the Monte Carlo example responsive while retaining an
# exact midpoint grid for the disk. All elements referencing one property use
# the same sampled value in a given realization.
number_of_elements = 10
number_of_grids = number_of_elements + 1
shaft_length = 1.0
disk_grid = number_of_grids // 2 + 1

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
            density=7850.0,
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
            kxx=(1,),
            kyy=3.0e5,
            cxx=500.0,
            cyy=500.0,
        ),
        DiskProperty(
            id=3,
            mass=5.0,
            diametral_inertia=0.025,
            polar_inertia=0.05,
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
    distributions=[
        RandomDistribution(
            id=1,
            name="horizontal bearing stiffness",
            distribution="uniform",
            parameters={
                # The broad, positive interval makes the propagated frequency
                # bands visible at the scale of the complete Campbell diagram.
                "low": 1.0e5,
                "high": 2.0e5,
            },
        ),
    ],
)

# A fixed seed makes the demonstration reproducible. Increasing the number of
# samples improves the empirical statistics at a proportional computational
# cost; it does not alter the requested speed discretization.
model = ModelBuilder().build(
    definition,
    stochastic=True,
    samples=40,
    seed=42,
)

# The solver accepts angular speeds in rad/s. Critical speeds are detected by
# interpolation between these supplied values only; no refinement is applied.
maximum_speed_rpm = 3_000.0
speeds_rpm = np.linspace(0.0, maximum_speed_rpm, 25)
speeds = speeds_rpm * 2.0 * np.pi / 60.0
result = Solver(model).solve(
    "campbell",
    fixed_dofs=[2, 5],
    speeds=speeds,
    modes=6,
    harmonics=[1.0, 2.0],
)

# Plot mean tracked branches with a central empirical 90% band. Critical-speed
# markers are sample-wise intersections grouped by tracked mode and harmonic;
# faint points expose their Monte Carlo spread along the 1x line.
figure, axes = plot_campbell(
    result,
    statistic="mean",
    confidence=0.90,
    show_harmonics=True,
    show_critical_samples=True,
    show_extremes=True,
    title="Stochastic Campbell diagram",
)
axes.set_xlim(0.0, maximum_speed_rpm)
#axes.set_ylim(14.0, 18.0)
figure.tight_layout()
plt.show()
