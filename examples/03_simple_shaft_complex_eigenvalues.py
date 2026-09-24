"""Extract complex eigenvalues of a deterministic spinning shaft."""

import numpy as np

from spinniped import (
    Grid,
    Material,
    ModelBuilder,
    ModelDefinition,
    ShaftElement,
    ShaftProperty,
    Solver,
)


# A nonzero spin speed activates the shaft's gyroscopic matrix.
spin_speed_rpm = 6_000.0
spin_speed = spin_speed_rpm * 2.0 * np.pi / 60.0

# Divide a one-metre shaft into equal finite elements.
number_of_elements = 20

# Build a deterministic uniform-shaft definition. The damping coefficient
# multiplies the shaft mass matrix, giving C = damping * M.
definition = ModelDefinition(
    grids=[
        Grid(grid_index + 1, z=grid_index / number_of_elements)
        for grid_index in range(number_of_elements + 1)
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
            outer_diameter=0.01,
            theory="timoshenko",
            rotary_inertia=True,
            damping=2.0,
        )
    ],
    elements=[
        ShaftElement(
            id=element_index + 1,
            grid_a=element_index + 1,
            grid_b=element_index + 2,
            property=1,
        )
        for element_index in range(number_of_elements)
    ],
)

# Assemble one realization: no random distributions are referenced.
model = ModelBuilder().build(definition)

# Apply simple supports at both ends. The first end also removes axial and
# torsional rigid-body motion; bending rotations remain free.
last_grid_first_dof = 6 * number_of_elements
fixed_dofs = [
    0,
    1,
    2,
    5,
    last_grid_first_dof,
    last_grid_first_dof + 1,
]

# Campbell analysis forms the first-order damped gyroscopic state matrix. A
# single requested speed is sufficient when only one operating point is needed.
result = Solver(model).solve(
    "campbell",
    fixed_dofs=fixed_dofs,
    speeds=[spin_speed],
    modes=8,
)

# Arrays retain sample and speed axes. Index [0, 0] selects the only
# deterministic realization and the only requested operating speed.
complex_eigenvalues = result["eigenvalues"][0, 0]

print(f"Spin speed: {spin_speed_rpm:.0f} rpm ({spin_speed:.2f} rad/s)")
print("Mode            eigenvalue [1/s]       frequency [Hz]   damping ratio")

for mode, eigenvalue in enumerate(complex_eigenvalues, start=1):
    decay_rate = eigenvalue.real
    damped_frequency = abs(eigenvalue.imag) / (2.0 * np.pi)
    damping_ratio = -decay_rate / abs(eigenvalue)
    print(
        f"{mode:4d}  "
        f"{decay_rate:+12.5f} {eigenvalue.imag:+12.5f}j  "
        f"{damped_frequency:14.5f}  "
        f"{damping_ratio:13.6e}"
    )

