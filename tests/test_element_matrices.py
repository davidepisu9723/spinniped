"""Physical and structural checks for the stateless element-matrix kernels."""

import numpy as np
import pytest

from spinniped.damping import (
    bearing_damping,
    disk_damping,
    shaft_damping,
)
from spinniped.gyroscopic import (
    bearing_gyroscopic,
    disk_gyroscopic,
    shaft_gyroscopic,
)
from spinniped.mass import bearing_mass, disk_mass, shaft_mass
from spinniped.stiffness import (
    bearing_stiffness,
    disk_stiffness,
    shaft_stiffness,
)


def test_shaft_matrices_have_expected_structure(shaft_local_matrices):
    stiffness = shaft_local_matrices["stiffness"]
    mass = shaft_local_matrices["mass"]
    gyroscopic = shaft_local_matrices["gyroscopic"]

    for matrix in (stiffness, mass, gyroscopic):
        assert matrix.shape == (12, 12)
        assert np.isfinite(matrix).all()
    assert np.allclose(stiffness, stiffness.T, rtol=1e-12)
    assert np.allclose(mass, mass.T, rtol=1e-12)
    assert np.allclose(gyroscopic, -gyroscopic.T, rtol=1e-12)


def test_shaft_stiffness_has_six_rigid_body_modes(shaft_local_matrices):
    eigenvalues = np.linalg.eigvalsh(shaft_local_matrices["stiffness"])
    scale = np.max(np.abs(eigenvalues))
    zero_modes = np.count_nonzero(np.abs(eigenvalues) < 1e-10 * scale)

    assert zero_modes == 6
    assert eigenvalues.min() > -1e-10 * scale


def test_shaft_mass_is_positive_definite(shaft_local_matrices):
    assert np.linalg.eigvalsh(shaft_local_matrices["mass"]).min() > 0.0


def test_axial_and_torsional_blocks_match_closed_form(shaft_properties):
    p = shaft_properties
    matrix = shaft_stiffness(
        length=p.length,
        outer_diameter=p.diameter,
        inner_diameter=0.0,
        young_modulus=p.young_modulus,
        poisson_ratio=p.poisson_ratio,
    )
    expected_axial = (
        p.young_modulus
        * p.area
        / p.length
        * np.array([[1.0, -1.0], [-1.0, 1.0]])
    )
    expected_torsional = (
        p.shear_modulus
        * p.polar_moment
        / p.length
        * np.array([[1.0, -1.0], [-1.0, 1.0]])
    )

    assert np.allclose(matrix[np.ix_([2, 8], [2, 8])], expected_axial)
    assert np.allclose(matrix[np.ix_([5, 11], [5, 11])], expected_torsional)


def test_density_scales_mass_and_gyroscopic_but_not_stiffness(shaft_properties):
    p = shaft_properties
    stiffness_options = {
        "length": p.length,
        "outer_diameter": p.diameter,
        "inner_diameter": 0.0,
        "young_modulus": p.young_modulus,
        "poisson_ratio": p.poisson_ratio,
    }
    mass_options = {**stiffness_options, "density": p.density}

    stiffness = shaft_stiffness(**stiffness_options)
    mass = shaft_mass(**mass_options)
    gyroscopic = shaft_gyroscopic(**mass_options)

    scaled_mass = shaft_mass(**{**mass_options, "density": 2.0 * p.density})
    scaled_gyroscopic = shaft_gyroscopic(
        **{**mass_options, "density": 2.0 * p.density}
    )
    repeated_stiffness = shaft_stiffness(**stiffness_options)

    assert np.allclose(scaled_mass, 2.0 * mass)
    assert np.allclose(scaled_gyroscopic, 2.0 * gyroscopic)
    assert np.array_equal(repeated_stiffness, stiffness)


def test_full_shaft_gyroscopic_matrix_couples_translations_and_rotations(
    shaft_properties,
):
    p = shaft_properties
    options = {
        "length": p.length,
        "outer_diameter": p.diameter,
        "inner_diameter": 0.0,
        "density": p.density,
        "young_modulus": p.young_modulus,
        "poisson_ratio": p.poisson_ratio,
        "theory": "timoshenko",
    }
    gyroscopic = shaft_gyroscopic(**options, rotary_inertia=True)
    translations = [0, 1, 6, 7]
    rotations = [3, 4, 9, 10]

    # A consistent shaft matrix contains more than nodal rotary coupling: the
    # distributed inertia also couples lateral translation and cross-section
    # rotation at both ends.
    assert np.count_nonzero(gyroscopic[np.ix_(translations, rotations)]) > 0
    assert np.count_nonzero(gyroscopic[np.ix_(rotations, translations)]) > 0
    assert np.allclose(gyroscopic, -gyroscopic.T)

    without_rotary_inertia = shaft_gyroscopic(
        **options, rotary_inertia=False
    )
    assert np.count_nonzero(without_rotary_inertia) == 0


def test_bearing_kernels_place_coefficients_in_translational_block():
    stiffness = bearing_stiffness(kxx=11.0, kxy=12.0, kyx=13.0,
                                  kyy=14.0, kzz=15.0)
    damping = bearing_damping(cxx=21.0, cxy=22.0, cyx=23.0,
                              cyy=24.0, czz=25.0)

    assert np.array_equal(
        stiffness[:3, :3],
        np.array([[11.0, 12.0, 0.0], [13.0, 14.0, 0.0], [0.0, 0.0, 15.0]]),
    )
    assert np.array_equal(
        damping[:3, :3],
        np.array([[21.0, 22.0, 0.0], [23.0, 24.0, 0.0], [0.0, 0.0, 25.0]]),
    )
    assert np.count_nonzero(bearing_mass()) == 0
    assert np.count_nonzero(bearing_gyroscopic()) == 0


def test_disk_kernels_represent_lumped_mass_and_polar_gyroscopic_coupling():
    lumped_mass = disk_mass(mass=2.0, diametral_inertia=0.1,
                            polar_inertia=0.2)
    gyroscopic = disk_gyroscopic(polar_inertia=0.2)

    assert np.array_equal(np.diag(lumped_mass), [2.0, 2.0, 2.0, 0.1, 0.1, 0.2])
    assert gyroscopic[3, 4] == pytest.approx(0.2)
    assert gyroscopic[4, 3] == pytest.approx(-0.2)
    assert np.allclose(gyroscopic, -gyroscopic.T)
    assert np.count_nonzero(disk_stiffness()) == 0


def test_element_damping_is_mass_proportional():
    shaft_mass_matrix = np.eye(12)
    disk_mass_matrix = np.eye(6)

    assert np.array_equal(
        shaft_damping(damping=0.03, mass_matrix=shaft_mass_matrix),
        0.03 * shaft_mass_matrix,
    )
    assert np.array_equal(
        disk_damping(damping=0.04, mass_matrix=disk_mass_matrix),
        0.04 * disk_mass_matrix,
    )


@pytest.mark.parametrize(
    "options",
    [
        {"length": 0.0, "outer_diameter": 0.01, "inner_diameter": 0.0},
        {"length": 1.0, "outer_diameter": 0.01, "inner_diameter": 0.01},
        {"length": 1.0, "outer_diameter": 0.01, "inner_diameter": -0.01},
    ],
)
def test_invalid_shaft_geometry_is_rejected(options, shaft_properties):
    with pytest.raises(ValueError, match="positive"):
        shaft_stiffness(
            **options,
            young_modulus=shaft_properties.young_modulus,
            poisson_ratio=shaft_properties.poisson_ratio,
        )
