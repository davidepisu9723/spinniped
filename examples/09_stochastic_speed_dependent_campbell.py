"""Plot a Campbell diagram with correlated speed-dependent bearing stiffness."""

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
    SinglePointConstraint,
    Solver,
    SpeedDependentValue,
    plot_campbell,
)


# Ten equal shaft elements place one grid exactly at the rotor midpoint. The
# central disk and identical end bearings form a symmetric Jeffcott-style rotor
# that is convenient for comparison with analytical models.
number_of_elements = 20
number_of_grids = number_of_elements + 1
shaft_length = 1.0
disk_grid = number_of_grids // 2 + 1

# The bearing stiffness softens with speed. Its three ordinates reference one
# multivariate distribution, so every realization contains a correlated curve
# rather than three independent random values. Rotor speeds are deterministic
# and expressed in rad/s, as required by SpeedDependentValue. Both table
# sequences are explicit tuples.
speed_knots_rpm = (0.0, 3_000.0, 6_000.0)
rpm_to_rad_per_second = 2.0 * np.pi / 60.0
speed_knots = tuple(
    speed_rpm * rpm_to_rad_per_second
    for speed_rpm in speed_knots_rpm
)
correlated_stiffness = SpeedDependentValue(
    speeds=speed_knots,
    values=((1, 0), (1, 1), (1, 2)),
)

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
            kxx=correlated_stiffness,
            kyy=correlated_stiffness,
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
    distributions=[
        RandomDistribution(
            id=1,
            name="correlated speed-dependent bearing stiffness",
            distribution="multivariate_normal",
            parameters={
                "mean": [1.50e6, 1.28e6, 1.00e6],
                # Use the same relative uncertainty at every speed knot.
                "stdv": [1.20e5, 1.024e5, 0.80e5],
                # Perfect correlation preserves the decreasing curve shape:
                # a realization that is stiffer at one speed is stiffer at
                # every other speed as well.
                "correlation": [
                    [1.00, 1.00, 1.00],
                    [1.00, 1.00, 1.00],
                    [1.00, 1.00, 1.00],
                ],
            },
        )
    ],
    boundary_conditions=[SinglePointConstraint(1, dof=(2, 5))],
)

# The fixed seed makes the Monte Carlo ensemble and its plotted bands
# reproducible. Both bearings and both lateral directions reuse the same
# sampled curve within a realization, preserving rotor symmetry.
model = ModelBuilder().build(
    definition,
    stochastic=True,
    samples=20,
    seed=42,
)

# Constrain the unsupported axial translation and torsional rotation. Campbell
# then evaluates each realization's interpolated bearing curve at every speed.
# Starting just above rest separates the forward and backward whirl vectors;
# at exactly zero speed their degeneracy can make stochastic mode matching
# exchange the two branches between otherwise similar realizations.
minimum_speed_rpm = 1.0
maximum_speed_rpm = speed_knots_rpm[-1]
speeds_rpm = np.linspace(minimum_speed_rpm, maximum_speed_rpm, 10)
speeds = speeds_rpm * rpm_to_rad_per_second
result = Solver(model).solve(
    "campbell",
    speeds=speeds,
    modes=10,
    harmonics=[1.0],
)

# Plot the four Jeffcott-relevant mean branches, their central empirical 90%
# bands, and every realization's critical-speed crossings along the 1x line.
figure, axes = plot_campbell(
    result,
    statistic="mean",
    confidence=0.90,
    show_harmonics=True,
    show_critical_samples=True,
    show_extremes=True,
    title="Stochastic speed-dependent bearing Campbell diagram",
)
axes.set_xlim(0.0, maximum_speed_rpm)
axes.set_ylim(bottom=0.0)
figure.tight_layout()
plt.show()
