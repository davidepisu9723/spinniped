"""Compare a flexible-bearing Jeffcott rotor with a finite-element model."""

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
)


# Extended Jeffcott data: a rigid central disk, a massless flexible shaft, and
# equal flexible bearings. A small positive FE shaft density approximates the
# massless assumption while keeping the assembled mass matrix nonsingular.
young_modulus = 2.0e11
shaft_diameter = 0.02
shaft_length = 1.0
shaft_density = 0.1
bearing_stiffness = 1.0e6
disk_mass = 5.0
disk_diametral_inertia = 0.025
disk_polar_inertia = 0.05
number_of_elements = 10
center_grid = number_of_elements // 2 + 1

# Combine the shaft and bearing compliances seen by center-disk translation
# and rotation. These two stiffnesses define the reduced Jeffcott reference.
second_moment = np.pi * shaft_diameter**4 / 64.0
flexural_rigidity = young_modulus * second_moment
translation_stiffness = 1.0 / (
    shaft_length**3 / (48.0 * flexural_rigidity)
    + 1.0 / (2.0 * bearing_stiffness)
)
rotation_stiffness = 1.0 / (
    shaft_length / (12.0 * flexural_rigidity)
    + 2.0 / (bearing_stiffness * shaft_length**2)
)

# Build the matching finite-element rotor. Shaft rotary inertia is disabled so
# the rigid disk supplies the model's only gyroscopic contribution.
definition = ModelDefinition(
    grids=[
        Grid(
            id=index + 1,
            z=index * shaft_length / number_of_elements,
        )
        for index in range(number_of_elements + 1)
    ],
    materials=[
        Material(1, shaft_density, young_modulus, 0.3),
    ],
    properties=[
        ShaftProperty(
            1,
            material=1,
            outer_diameter=shaft_diameter,
            theory="euler",
            rotary_inertia=False,
        ),
        BearingProperty(
            2,
            kxx=bearing_stiffness,
            kyy=bearing_stiffness,
        ),
        DiskProperty(
            3,
            mass=disk_mass,
            diametral_inertia=disk_diametral_inertia,
            polar_inertia=disk_polar_inertia,
        ),
    ],
    elements=[
        *[
            ShaftElement(index + 1, index + 1, index + 2, 1)
            for index in range(number_of_elements)
        ],
        BearingElement(number_of_elements + 1, grid=1, property=2),
        BearingElement(
            number_of_elements + 2,
            grid=number_of_elements + 1,
            property=2,
        ),
        DiskElement(number_of_elements + 3, grid=center_grid, property=3),
    ],
)
model = ModelBuilder().build(definition)

# The Jeffcott reduction describes lateral translation and bending rotation
# only, so remove axial and torsional DOFs from every grid.
fixed_dofs = [
    6 * grid_index + local_dof
    for grid_index in range(number_of_elements + 1)
    for local_dof in (2, 5)
]

# Evaluate both models through and beyond the two analytical critical speeds.
maximum_speed_rpm = 6_000.0
speeds_rpm = np.linspace(0.0, maximum_speed_rpm, 41)
speeds = speeds_rpm * 2.0 * np.pi / 60.0
result = Solver(model).solve(
    "campbell",
    fixed_dofs=fixed_dofs,
    speeds=speeds,
    modes=4,
    # Exact zero-speed degeneracy makes the initial MAC assignment arbitrary.
    track_modes=False,
)

# Evaluate the analytical cylindrical pair and gyroscopically split conical
# pair at the same speeds as the FE solution.
translation_frequency = np.sqrt(translation_stiffness / disk_mass)
analytical_frequencies = []
for speed in speeds:
    gyroscopic_term = disk_polar_inertia * speed
    discriminant = np.sqrt(
        gyroscopic_term**2
        + 4.0 * disk_diametral_inertia * rotation_stiffness
    )
    backward_conical = (
        discriminant - gyroscopic_term
    ) / (2.0 * disk_diametral_inertia)
    forward_conical = (
        discriminant + gyroscopic_term
    ) / (2.0 * disk_diametral_inertia)
    analytical_frequencies.append(
        np.sort(
            [
                translation_frequency,
                translation_frequency,
                backward_conical,
                forward_conical,
            ]
        )
        / (2.0 * np.pi)
    )
analytical_frequencies = np.asarray(analytical_frequencies)

# Solid colored curves are FE results; dashed black curves are the reduced
# Jeffcott references. Close overlap confirms the intended comparison.
figure, axes = plot_campbell(
    result,
    title="Flexible-bearing Jeffcott rotor: analytical and FE comparison",
)
for mode_index in range(analytical_frequencies.shape[1]):
    axes.plot(
        speeds_rpm,
        analytical_frequencies[:, mode_index],
        color="black",
        linestyle="--",
        linewidth=1.0,
        label="Analytical Jeffcott" if mode_index == 0 else "_nolegend_",
    )

# The 1x line identifies critical speeds where shaft speed equals a modal
# frequency. Mark the cylindrical and backward-conical analytical crossings.
axes.plot(
    speeds_rpm,
    speeds_rpm / 60.0,
    color="tab:red",
    linestyle=":",
    linewidth=1.2,
    label="1x synchronous",
)
cylindrical_critical = translation_frequency
backward_conical_critical = np.sqrt(
    rotation_stiffness / (disk_diametral_inertia + disk_polar_inertia)
)
for label, critical_speed in (
    ("Cylindrical critical", cylindrical_critical),
    ("Backward-conical critical", backward_conical_critical),
):
    critical_speed_rpm = critical_speed * 60.0 / (2.0 * np.pi)
    axes.plot(
        critical_speed_rpm,
        critical_speed / (2.0 * np.pi),
        marker="o",
        linestyle="none",
        label=label,
    )
    print(f"{label}: {critical_speed_rpm:.1f} rpm")

axes.set_xlim(0.0, maximum_speed_rpm)
axes.set_ylim(bottom=0.0)
axes.legend()
figure.tight_layout()
plt.show()
