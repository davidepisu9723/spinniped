"""Finite-element results compared with independent analytical solutions."""

import numpy as np

from conftest import (
    build_uniform_beam,
    one_frequency_per_bending_pair,
    positive_frequencies,
    simply_supported_frequencies,
)


def test_simply_supported_beam_first_four_bending_frequencies(shaft_properties):
    # ``ne`` is the number of finite elements. Sixteen elements provide enough
    # resolution to compare the first four modes with a tight tolerance.
    ne = 16

    # ``nm`` is the number of analytical bending modes to verify.
    nm = 4

    # Euler-Bernoulli theory is selected because the reference equation below
    # is derived from Euler-Bernoulli beam theory without rotary inertia.
    s = build_uniform_beam(shaft_properties, ne, theory="euler")

    # Simply supported boundary conditions require zero transverse displacement
    # at both ends while allowing bending rotations.
    #
    # At the first node, x and y translations are fixed. Axial translation z
    # is additionally fixed to remove axial rigid-body motion, and torsional
    # rotation tz is fixed to remove rigid rotation about the shaft axis.
    s.constrain_nodes([0], constrained_dofs_per_node=[0, 1, 2, 5])

    # At the final node, only transverse x and y translations are constrained.
    # The value ``ne`` is also the zero-based internal index of the final node
    # because a mesh of ne elements contains ne + 1 nodes.
    s.constrain_nodes([ne], constrained_dofs_per_node=[0, 1])

    # Solve K_ff*q = lambda*M_ff*q. Eigenvectors are not needed in this test.
    eigenvalues, _ = s.solve_eigenproblem()

    # Convert lambda to frequency and keep one value from each duplicate x/y
    # bending pair of the circular shaft.
    f_num = one_frequency_per_bending_pair(
        positive_frequencies(eigenvalues), nm
    )

    # Compute the independent closed-form Euler-Bernoulli frequencies.
    f_ana = simply_supported_frequencies(shaft_properties, nm)

    # Compare each mode through a dimensionless relative error.
    relative_error = np.abs(f_num - f_ana) / f_ana

    # Every one of the first four numerical modes must agree with the
    # analytical value within 0.1%. The custom message prints all modal errors
    # if the assertion fails.
    assert np.all(relative_error < 1.0e-3), (
        "Simply supported frequency errors exceeded 0.1%: "
        f"{relative_error}"
    )
