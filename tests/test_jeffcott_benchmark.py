"""Finite-element comparison with the classical Jeffcott rotor."""

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
)


YOUNG_MODULUS = 2.0e11
SHAFT_DIAMETER = 0.02
SHAFT_LENGTH = 1.0
DISK_MASS = 5.0
DISK_DIAMETRAL_INERTIA = 0.025
DISK_POLAR_INERTIA = 0.05
BEARING_STIFFNESS = 1.0e6
ELEMENTS = 10


def _jeffcott_stiffnesses():
    """Return equivalent center translation and rotation stiffnesses."""
    second_moment = np.pi * SHAFT_DIAMETER**4 / 64.0
    flexural_rigidity = YOUNG_MODULUS * second_moment

    # Shaft and bearing compliances add in series for symmetric translation
    # and antisymmetric disk rotation, respectively.
    translation_stiffness = 1.0 / (
        SHAFT_LENGTH**3 / (48.0 * flexural_rigidity)
        + 1.0 / (2.0 * BEARING_STIFFNESS)
    )
    rotation_stiffness = 1.0 / (
        SHAFT_LENGTH / (12.0 * flexural_rigidity)
        + 2.0 / (BEARING_STIFFNESS * SHAFT_LENGTH**2)
    )
    return translation_stiffness, rotation_stiffness


def _jeffcott_frequencies(speeds):
    """Return analytical cylindrical and conical frequencies in hertz."""
    translation_stiffness, rotation_stiffness = _jeffcott_stiffnesses()
    translation_frequency = np.sqrt(translation_stiffness / DISK_MASS)
    frequencies = []

    for speed in np.asarray(speeds, dtype=float):
        gyroscopic_term = DISK_POLAR_INERTIA * speed
        discriminant = np.sqrt(
            gyroscopic_term**2
            + 4.0 * DISK_DIAMETRAL_INERTIA * rotation_stiffness
        )
        backward_conical = (
            discriminant - gyroscopic_term
        ) / (2.0 * DISK_DIAMETRAL_INERTIA)
        forward_conical = (
            discriminant + gyroscopic_term
        ) / (2.0 * DISK_DIAMETRAL_INERTIA)
        frequencies.append(
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

    return np.asarray(frequencies)


def _build_comparable_finite_element_model():
    """Build an FE model that closely reproduces Jeffcott assumptions."""
    grids = [
        Grid(index + 1, z=index * SHAFT_LENGTH / ELEMENTS)
        for index in range(ELEMENTS + 1)
    ]
    definition = ModelDefinition(
        grids=grids,
        materials=[
            # Kernels require positive density. This negligible value makes
            # the centered disk dominate inertia like a massless Jeffcott shaft.
            Material(1, 0.1, YOUNG_MODULUS, 0.3)
        ],
        properties=[
            ShaftProperty(
                1,
                material=1,
                outer_diameter=SHAFT_DIAMETER,
                theory="euler",
                rotary_inertia=False,
            ),
            BearingProperty(
                2,
                kxx=BEARING_STIFFNESS,
                kyy=BEARING_STIFFNESS,
            ),
            DiskProperty(
                3,
                mass=DISK_MASS,
                diametral_inertia=DISK_DIAMETRAL_INERTIA,
                polar_inertia=DISK_POLAR_INERTIA,
            ),
        ],
        elements=[
            *[
                ShaftElement(index + 1, index + 1, index + 2, 1)
                for index in range(ELEMENTS)
            ],
            BearingElement(ELEMENTS + 1, grid=1, property=2),
            BearingElement(ELEMENTS + 2, grid=ELEMENTS + 1, property=2),
            DiskElement(ELEMENTS + 3, grid=ELEMENTS // 2 + 1, property=3),
        ],
    )
    return ModelBuilder().build(definition)


def _non_lateral_dofs():
    """Remove axial and torsional motion excluded by the Jeffcott model."""
    return [
        6 * grid_index + local_dof
        for grid_index in range(ELEMENTS + 1)
        for local_dof in (2, 5)
    ]


def test_flexible_bearing_jeffcott_frequencies_match_finite_element_model():
    """Compare analytical cylindrical and conical modes with the FE model."""
    analytical_frequencies = _jeffcott_frequencies([0.0])[0]
    model = _build_comparable_finite_element_model()
    result = Solver(model).solve(
        "modal",
        fixed_dofs=_non_lateral_dofs(),
        modes=4,
        # The near-massless shaft creates very high parasitic modes; retain the
        # much smaller physical disk modes with an explicit zero threshold.
        zero_tolerance=0.0,
    )

    # At rest, isotropy produces degenerate cylindrical and conical mode pairs.
    assert np.allclose(
        np.sort(result["frequencies"][0]),
        analytical_frequencies,
        rtol=1.0e-5,
    )


def test_jeffcott_gyroscopic_shift_and_critical_speeds_match_finite_elements():
    """Compare speed-dependent splitting and 1x critical crossings."""
    translation_stiffness, rotation_stiffness = _jeffcott_stiffnesses()
    cylindrical_critical = np.sqrt(translation_stiffness / DISK_MASS)
    backward_conical_critical = np.sqrt(
        rotation_stiffness
        / (DISK_DIAMETRAL_INERTIA + DISK_POLAR_INERTIA)
    )
    speeds = np.array(
        [0.0, cylindrical_critical, 200.0, backward_conical_critical]
    )
    analytical_frequencies = _jeffcott_frequencies(speeds)
    model = _build_comparable_finite_element_model()
    result = Solver(model).solve(
        "campbell",
        fixed_dofs=_non_lateral_dofs(),
        speeds=speeds,
        modes=4,
        # The two exact conical roots are degenerate at rest. Independent
        # frequency ordering avoids an arbitrary zero-speed MAC assignment.
        track_modes=False,
    )

    # Sorting isolates frequency agreement from the labels assigned to the
    # degenerate conical pair at zero speed.
    finite_element_frequencies = np.sort(result["frequencies"][0], axis=1)
    assert np.allclose(
        finite_element_frequencies,
        analytical_frequencies,
        rtol=1.0e-5,
    )

    # The cylindrical modes remain fixed while conical gyroscopic branches
    # separate above zero speed.
    assert np.ptp(analytical_frequencies[2, 2:]) > 0.0
    assert np.ptp(finite_element_frequencies[2, 2:]) > 0.0

    # Each analytical critical speed lies on a 1x modal crossing.
    np.testing.assert_allclose(
        result["speeds_hz"][1],
        analytical_frequencies[1, 0],
        rtol=1.0e-5,
    )
    np.testing.assert_allclose(
        result["speeds_hz"][-1],
        analytical_frequencies[-1, 2],
        rtol=1.0e-5,
    )
