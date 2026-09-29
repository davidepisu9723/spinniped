"""Plot the Campbell diagram of a two-bearing shaft with one disk."""

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
    ShaftElement,
    ShaftProperty,
    Solver,
    plot_campbell,
    plot_rotor,
)


# A ten-element mesh places grid 6 exactly at the shaft midpoint, where the
# rigid disk is attached.
number_of_elements = 10
number_of_grids = number_of_elements + 1
shaft_length = 1.0
disk_grid = number_of_grids // 2 + 1

# IDs following the shaft-element range keep all element records unique.
left_bearing_id = number_of_elements + 1
right_bearing_id = number_of_elements + 2
disk_element_id = number_of_elements + 3

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
            outer_diameter=0.02,
            theory="timoshenko",
            rotary_inertia=True,
            damping=0.2,
        ),
        BearingProperty(
            id=2,
            kxx=1.0e7,
            kyy=1.0e7,
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
        # Join each pair of neighboring grids with a shaft element.
        *[
            ShaftElement(
                id=element_index + 1,
                grid_a=element_index + 1,
                grid_b=element_index + 2,
                property=1,
            )
            for element_index in range(number_of_elements)
        ],
        # Grounded bearings support the two shaft ends in global x and y.
        BearingElement(
            id=left_bearing_id,
            grid=1,
            property=2,
        ),
        BearingElement(
            id=right_bearing_id,
            grid=number_of_grids,
            property=2,
        ),
        # The disk polar axis follows the default global z spin axis.
        DiskElement(
            id=disk_element_id,
            grid=disk_grid,
            property=3,
        ),
    ],
)

# Assemble one deterministic set of mass, stiffness, damping, and gyroscopic
# matrices. There are no random-distribution references in this definition.
model = ModelBuilder().build(definition)

# Plot the assembled rotor before solving it. Shaft diameters follow their
# resolved section properties; the disk and bearings use the conventional
# longitudinal-section symbols described in the plotting user guide.
rotor_figure, rotor_axes = plot_rotor(
    model,
    title="Two-bearing shaft with a central disk",
    show_grid_ids=True,
)
rotor_figure.tight_layout()

# Bearings remove lateral rigid motion. Constrain axial translation and
# torsional rotation at the first grid because those directions are unsupported.
fixed_dofs = [2, 5]

# Sweep from rest to 12,000 rpm. The solver expects angular speed in rad/s and
# tracks each mode between adjacent speed points by default.
maximum_speed_rpm = 12_000.0
speeds_rpm = np.linspace(0.0, maximum_speed_rpm, 31)
speeds = speeds_rpm * 2.0 * np.pi / 60.0
result = Solver(model).solve(
    "campbell",
    fixed_dofs=fixed_dofs,
    speeds=speeds,
    modes=8,
)

# Plot the tracked natural-frequency branches. Rotor speed defaults to rpm on
# the horizontal axis, matching the array used for the synchronous reference.
figure, axes = plot_campbell(
    result,
    title="Two-bearing shaft with a central disk",
)

# A once-per-revolution excitation has frequency speed_rpm / 60 in hertz.
# Intersections with modal branches indicate candidate 1x critical speeds.
axes.plot(
    speeds_rpm,
    speeds_rpm / 60.0,
    color="black",
    linestyle="--",
    linewidth=1.0,
    label="1x synchronous",
)
axes.set_xlim(0.0, maximum_speed_rpm)
axes.set_ylim(bottom=0.0)
axes.legend()
figure.tight_layout()
plt.show()
