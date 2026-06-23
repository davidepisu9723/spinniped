"""Tests for local shaft-element mass and stiffness matrices."""

import numpy as np


def test_element_matrices_have_expected_shape_and_finite_entries(shaft_element):
    # A two-node shaft has 6 DOFs per node, therefore each local matrix must
    # contain 12 rows and 12 columns.
    assert shaft_element.K.shape == (12, 12)
    assert shaft_element.M.shape == (12, 12)

    # NaN or infinity usually indicates invalid geometry, a division by zero,
    # or another numerical error during matrix construction.
    assert np.isfinite(shaft_element.K).all()
    assert np.isfinite(shaft_element.M).all()


def test_element_matrices_are_symmetric(shaft_element):
    # For a conservative linear elastic element, Maxwell-Betti reciprocity
    # requires K = K.T. A consistent mass matrix must also be symmetric.
    assert np.allclose(shaft_element.K, shaft_element.K.T, rtol=1e-12)
    assert np.allclose(shaft_element.M, shaft_element.M.T, rtol=1e-12)


def test_stiffness_matrix_has_six_rigid_body_modes(shaft_element):
    # ``eigvalsh`` is used because K is real and symmetric. Its eigenvalues
    # represent stiffnesses along the element's independent modal directions.
    eigenvalues = np.linalg.eigvalsh(shaft_element.K)

    # Normalize the zero threshold with the largest stiffness eigenvalue.
    # Matrix entries span many orders of magnitude, so an absolute tolerance
    # would be sensitive to the chosen units and shaft dimensions.
    scale = np.max(np.abs(eigenvalues))
    number_of_zero_modes = np.count_nonzero(np.abs(eigenvalues) < 1e-10 * scale)

    # An unconstrained three-dimensional body can translate along x, y, z and
    # rotate about x, y, z without developing strain: six rigid-body modes.
    assert number_of_zero_modes == 6

    # K should be positive semidefinite. Tiny negative eigenvalues are accepted
    # because floating-point round-off can perturb an exact zero slightly.
    assert eigenvalues.min() > -1e-10 * scale


def test_mass_matrix_is_positive_definite(shaft_element):
    # Positive definiteness means every nonzero velocity field has positive
    # kinetic energy: T = 0.5 * q_dot.T * M * q_dot > 0.
    eigenvalues = np.linalg.eigvalsh(shaft_element.M)
    assert eigenvalues.min() > 0.0


def test_axial_and_torsional_stiffness_blocks_match_closed_form(shaft_element):
    # Exact two-node axial stiffness in local DOF order [z0, z1].
    expected_axial = (
        shaft_element.E
        * shaft_element.A
        / shaft_element.L
        * np.array([[1.0, -1.0], [-1.0, 1.0]])
    )

    # Exact two-node torsional stiffness in local DOF order [tz0, tz1].
    expected_torsion = (
        shaft_element.G
        * shaft_element.J
        / shaft_element.L
        * np.array([[1.0, -1.0], [-1.0, 1.0]])
    )

    # In Spinniped's 12-DOF ordering, axial displacement z occupies local
    # indices [2, 8].
    assert np.allclose(
        shaft_element.K[np.ix_([2, 8], [2, 8])], expected_axial
    )

    # Rotation about the shaft axis tz occupies local indices [5, 11].
    assert np.allclose(
        shaft_element.K[np.ix_([5, 11], [5, 11])], expected_torsion
    )


def test_density_scales_mass_but_not_stiffness(shaft_element):
    # Preserve the original matrices so the recomputed values can be compared
    # with the same element geometry and elastic properties.
    original_mass = shaft_element.M.copy()
    original_stiffness = shaft_element.K.copy()

    # Density appears linearly in the mass matrix but not in the elastic
    # stiffness matrix. Doubling rho should therefore double M and leave K
    # unchanged.
    shaft_element.rho *= 2.0
    scaled_mass = shaft_element._M()
    scaled_stiffness = shaft_element._K()

    assert np.allclose(scaled_mass, 2.0 * original_mass)
    assert np.allclose(scaled_stiffness, original_stiffness)
