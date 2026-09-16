"""Declarative finite-element results compared with analytical solutions."""

import numpy as np

from spinniped import Solver

from conftest import (
    build_uniform_beam,
    one_frequency_per_bending_pair,
    simply_supported_fixed_dofs,
    simply_supported_frequencies,
)


def test_simply_supported_beam_first_four_bending_frequencies(shaft_properties):
    """Compare four numerical bending modes with analytical frequencies."""
    # A fine Euler mesh isolates assembly accuracy from shear deformation.
    elements = 16
    modes = 4
    model = build_uniform_beam(shaft_properties, elements, theory="euler")
    result = Solver(model).solve(
        "modal",
        fixed_dofs=simply_supported_fixed_dofs(elements),
        modes=2 * modes,
    )

    numerical = one_frequency_per_bending_pair(result["frequencies"][0], modes)
    analytical = simply_supported_frequencies(shaft_properties, modes)
    relative_error = np.abs(numerical - analytical) / analytical

    assert np.all(relative_error < 1.0e-3), (
        "Simply supported frequency errors exceeded 0.1%: "
        f"{relative_error}"
    )
