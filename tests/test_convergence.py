"""Mesh-convergence tests through ModelBuilder and Solver.solve()."""

import numpy as np
import pytest

from spinniped import Solver

from conftest import (
    build_uniform_beam,
    one_frequency_per_bending_pair,
    simply_supported_fixed_dofs,
    simply_supported_frequencies,
)


def _first_mode_error(properties, elements):
    """Return relative first-mode error for one mesh density.

    Parameters
    ----------
    properties : ShaftProperties
        Uniform beam properties.
    elements : int
        Number of finite elements in the mesh.

    Returns
    -------
    float
        Absolute relative error in the first bending frequency.
    """
    # Solve the finite-element model with simply supported end conditions.
    model = build_uniform_beam(properties, elements, theory="euler")
    result = Solver(model).solve(
        "modal",
        fixed_dofs=simply_supported_fixed_dofs(elements),
        modes=2,
    )
    numerical = one_frequency_per_bending_pair(result["frequencies"][0], 1)[0]
    analytical = simply_supported_frequencies(properties, 1)[0]
    return abs(numerical - analytical) / analytical


@pytest.mark.parametrize("elements", [2, 4, 8, 16])
def test_first_frequency_is_reasonable_at_each_mesh(shaft_properties, elements):
    """Verify acceptable first-mode accuracy for every tested mesh."""
    assert _first_mode_error(shaft_properties, elements) < 0.01


def test_first_frequency_converges_under_mesh_refinement(shaft_properties):
    """Verify that first-mode error decreases as the mesh is refined."""
    # Evaluate the same physical beam at progressively finer resolutions.
    element_counts = [2, 4, 8, 16]
    errors = np.array(
        [_first_mode_error(shaft_properties, count) for count in element_counts]
    )

    assert errors[-1] < errors[0]
    assert errors[-1] < 1.0e-4
