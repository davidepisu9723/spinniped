"""Mesh-convergence tests for analytical beam benchmarks."""

import numpy as np
import pytest

from conftest import (
    build_uniform_beam,
    one_frequency_per_bending_pair,
    positive_frequencies,
    simply_supported_frequencies,
)


def _first_mode_error(p, ne):
    # Build the same physical shaft with a selectable number of elements.
    # Only mesh density changes; material, section, and total length do not.
    s = build_uniform_beam(p, ne, theory="euler")

    # Apply the same simply supported constraints used by the analytical
    # benchmark. The extra z and tz constraints remove unrelated rigid modes.
    s.constrain_nodes([0], constrained_dofs_per_node=[0, 1, 2, 5])
    s.constrain_nodes([ne], constrained_dofs_per_node=[0, 1])

    # Solve the finite-element model and extract the first bending frequency.
    eigenvalues, _ = s.solve_eigenproblem()
    f_num = one_frequency_per_bending_pair(
        positive_frequencies(eigenvalues), 1
    )[0]

    # Compute the exact first Euler-Bernoulli bending frequency.
    f_ana = simply_supported_frequencies(p, 1)[0]

    # Return a dimensionless error so results can be compared across meshes.
    return abs(f_num - f_ana) / f_ana


# Parametrization tells pytest to run this function as four independent test
# cases. A failure report will identify the particular value of ``ne``.
@pytest.mark.parametrize("ne", [2, 4, 8, 16])
def test_first_frequency_is_reasonable_at_each_mesh(shaft_properties, ne):
    # Even each individual mesh should produce a physically reasonable first
    # frequency. The one-percent limit also catches severe assembly or support
    # errors before the convergence trend is examined.
    error = _first_mode_error(shaft_properties, ne)
    assert error < 0.01


def test_first_frequency_converges_under_mesh_refinement(shaft_properties):
    # Use a geometric refinement sequence in which every step doubles the
    # number of elements and halves their length.
    ne_values = [2, 4, 8, 16]

    # Rebuild and solve the complete model at every mesh density.
    errors = np.array([_first_mode_error(shaft_properties, ne) for ne in ne_values])

    # The finest mesh must improve upon the coarsest mesh. This deliberately
    # avoids requiring strict monotonic decrease at every step, which can make
    # numerical tests unnecessarily brittle.
    assert errors[-1] < errors[0]

    # Besides showing an improving trend, the final mesh must achieve an
    # absolute quality target: relative error below 0.01%.
    assert errors[-1] < 1.0e-4
