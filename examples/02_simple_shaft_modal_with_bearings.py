"""Compute modal frequencies of a uniform shaft supported by bearings."""

from spinniped import (
    BearingElement,
    BearingProperty,
    Grid,
    Material,
    ModelBuilder,
    ModelDefinition,
    ShaftElement,
    ShaftProperty,
    Solver,
)


# Divide a one-metre shaft into equal finite elements.
number_of_elements = 21

# Build the NASTRAN-like input records.
definition = ModelDefinition(
    # Place all grids along the global spin axis.
    grids=[
        Grid(grid_index + 1, z=grid_index / number_of_elements)
        for grid_index in range(number_of_elements + 1)
    ],
    # Define one shared steel-like material.
    materials=[
        Material(
            id=1,
            density=7850.0,
            young_modulus=2.0e11,
            poisson_ratio=0.3,
        )
    ],
    # Define the shaft section and bearing support coefficients.
    properties=[
        ShaftProperty(
            id=1,
            material=1,
            outer_diameter=0.01,
        ),
        BearingProperty(
            id=2,
            kxx=1.0e12,
            kyy=1.0e12,
        ),
    ],
    # Create the shaft chain and one support-to-ground bearing at each end.
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
            grid=number_of_elements + 1,
            property=2,
        ),
    ],
)

# Assemble one deterministic numerical model.
model = ModelBuilder().build(definition)

# Remove only the unsupported axial and torsional rigid-body motion.
result = Solver(model).solve(
    "modal",
    fixed_dofs=[2, 5],
    modes=18,
)

# Print frequencies from the single deterministic sample.
for frequency in result["frequencies"][0]:
    print(f"Frequency: {frequency:.2f} Hz")
