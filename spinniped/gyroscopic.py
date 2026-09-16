"""Unit-speed gyroscopic-matrix kernels."""

import numpy as np


def shaft_gyroscopic(
    *,
    length,
    outer_diameter,
    inner_diameter,
    density,
    young_modulus,
    poisson_ratio,
    theory="timoshenko",
    rotary_inertia=True,
):
    """Build the consistent unit-speed gyroscopic matrix of a shaft.

    Parameters
    ----------
    length : float
        Element length.
    outer_diameter : float
        Shaft outer diameter.
    inner_diameter : float
        Shaft bore diameter.
    density : float
        Material mass density.
    young_modulus : float
        Material Young's modulus.
    poisson_ratio : float
        Material Poisson's ratio.
    theory : {"timoshenko", "euler"}, optional
        Beam interpolation used by the gyroscopic coefficients.
    rotary_inertia : bool, optional
        Enable the distributed gyroscopic contribution.

    Returns
    -------
    numpy.ndarray
        Skew-symmetric local ``(12, 12)`` unit-speed matrix.

    Raises
    ------
    TypeError
        If ``rotary_inertia`` is not Boolean.
    ValueError
        If a physical parameter is invalid or the theory is unknown.
    """
    # Validate all numerical inputs before evaluating coefficient formulas.
    values = np.asarray(
        [
            length,
            outer_diameter,
            inner_diameter,
            density,
            young_modulus,
            poisson_ratio,
        ],
        dtype=float,
    )
    if not np.isfinite(values).all():
        raise ValueError("Shaft gyroscopic parameters must be finite")

    # Require an explicit Boolean switch to avoid accidental truthy strings.
    if not isinstance(rotary_inertia, (bool, np.bool_)):
        raise TypeError("rotary_inertia must be a boolean")

    # Validate geometry and material properties.
    if length <= 0 or outer_diameter <= inner_diameter or inner_diameter < 0:
        raise ValueError("Shaft length and wall dimensions must be positive")
    if density <= 0:
        raise ValueError("Density must be positive")
    if young_modulus <= 0:
        raise ValueError("Young's modulus must be positive")
    if not -1.0 < poisson_ratio < 0.5:
        raise ValueError("Poisson's ratio must be between -1 and 0.5")

    # Disabling rotary inertia also disables distributed gyroscopic inertia.
    if not rotary_inertia:
        return np.zeros((12, 12), dtype=float)

    # Compute annular-section area and diametral geometric inertia.
    area = np.pi * (outer_diameter**2 - inner_diameter**2) / 4.0
    diametral_inertia = (
        np.pi * (outer_diameter**4 - inner_diameter**4) / 64.0
    )

    # Determine the shear parameter associated with the beam interpolation.
    if theory == "timoshenko":
        shear_modulus = young_modulus / (2.0 * (1.0 + poisson_ratio))
        diameter_ratio = inner_diameter / outer_diameter
        ratio_term = (1.0 + diameter_ratio**2) ** 2
        shear_factor = (
            6.0 * (1.0 + poisson_ratio) * ratio_term
            / (
                (7.0 + 6.0 * poisson_ratio) * ratio_term
                + (20.0 + 12.0 * poisson_ratio) * diameter_ratio**2
            )
        )
        shear_parameter = (
            12.0 * young_modulus * diametral_inertia
            / (shear_factor * shear_modulus * area * length**2)
        )
    elif theory == "euler":
        # Euler-Bernoulli theory is the zero-shear-flexibility limit.
        shear_parameter = 0.0
    else:
        raise ValueError(f"Unknown shaft theory {theory!r}")

    # Define the four coefficients of the consistent lateral matrix.
    g1 = 36.0
    g2 = (3.0 - 15.0 * shear_parameter) * length
    g3 = (
        4.0 + 5.0 * shear_parameter + 10.0 * shear_parameter**2
    ) * length**2
    g4 = (
        -1.0 - 5.0 * shear_parameter + 5.0 * shear_parameter**2
    ) * length**2

    # Scale the unit-speed matrix by distributed diametral mass inertia.
    scale = (
        density
        * diametral_inertia
        / (15.0 * length * (1.0 + shear_parameter) ** 2)
    )

    # Build the skew-symmetric 8-DOF lateral gyroscopic block.
    lateral_matrix = scale * np.array(
        [
            [0.0, g1, -g2, 0.0, 0.0, -g1, -g2, 0.0],
            [-g1, 0.0, 0.0, -g2, g1, 0.0, 0.0, -g2],
            [g2, 0.0, 0.0, g3, -g2, 0.0, 0.0, g4],
            [0.0, g2, -g3, 0.0, 0.0, -g2, -g4, 0.0],
            [0.0, -g1, g2, 0.0, 0.0, g1, g2, 0.0],
            [g1, 0.0, 0.0, g2, -g1, 0.0, 0.0, g2],
            [g2, 0.0, 0.0, g4, -g2, 0.0, 0.0, g3],
            [0.0, g2, -g4, 0.0, 0.0, -g2, -g3, 0.0],
        ]
    )

    # Insert the lateral block into the package's complete 12-DOF ordering.
    matrix = np.zeros((12, 12), dtype=float)
    lateral_dofs = [0, 1, 3, 4, 6, 7, 9, 10]
    matrix[np.ix_(lateral_dofs, lateral_dofs)] = lateral_matrix

    return matrix


def bearing_gyroscopic():
    """Build the local gyroscopic matrix of a non-rotating bearing.

    Returns
    -------
    numpy.ndarray
        Zero ``(6, 6)`` matrix.
    """
    # Bearings support the rotor but do not rotate with it in this model.
    return np.zeros((6, 6), dtype=float)


def disk_gyroscopic(*, polar_inertia=0.0):
    """Build the local unit-speed gyroscopic matrix of a rigid disk.

    Parameters
    ----------
    polar_inertia : float, optional
        Disk mass moment of inertia about its local spin axis.

    Returns
    -------
    numpy.ndarray
        Skew-symmetric local ``(6, 6)`` unit-speed matrix.

    Raises
    ------
    ValueError
        If ``polar_inertia`` is non-finite or negative.
    """
    # A physically meaningful polar mass moment is finite and nonnegative.
    if not np.isfinite(polar_inertia):
        raise ValueError("Disk polar inertia must be finite")
    if polar_inertia < 0:
        raise ValueError("Disk polar inertia cannot be negative")

    # Couple transverse rotations about local x and y.
    matrix = np.zeros((6, 6), dtype=float)
    matrix[3, 4] = polar_inertia
    matrix[4, 3] = -polar_inertia

    return matrix
