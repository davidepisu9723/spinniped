"""Compute the first modes of a simply supported uniform shaft."""

from spinniped import (
    Grid,
    Material,
    ModelBuilder,
    ModelDefinition,
    ShaftElement,
    ShaftProperty,
    SinglePointConstraint,
    Solver,
)


# Divide a one-metre shaft into equal finite elements.
number_of_elements = 21

# Create one declarative model from ordinary Python records.
definition = ModelDefinition(
    # A chain of N elements requires N + 1 grids.
    grids=[
        Grid(grid_index + 1, z=grid_index / number_of_elements)
        for grid_index in range(number_of_elements + 1)
    ],
    # All shaft elements share this steel-like material.
    materials=[
        Material(
            id=1,
            density=7850.0,
            young_modulus=2.0e11,
            poisson_ratio=0.3,
        )
    ],
    # The shaft chain also shares one circular section property.
    properties=[
        ShaftProperty(
            id=1,
            material=1,
            outer_diameter=0.01,
        )
    ],
    # Connect each grid to its immediate successor.
    elements=[
        ShaftElement(
            id=element_index + 1,
            grid_a=element_index + 1,
            grid_b=element_index + 2,
            property=1,
        )
        for element_index in range(number_of_elements)
    ],
    boundary_conditions=[
        SinglePointConstraint(1, dof=(0, 1, 2, 5)),
        SinglePointConstraint(number_of_elements + 1, dof=(0, 1)),
    ],
)

# Resolve records and assemble one deterministic set of global matrices.
model = ModelBuilder().build(definition)

# Request the first eighteen positive real modes.
result = Solver(model).solve(
    "modal",
    modes=18,
)

# A deterministic result still retains a leading sample dimension of length one.
for frequency in result["frequencies"][0]:
    print(f"Frequency: {frequency:.2f} Hz")
