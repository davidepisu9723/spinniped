"""Element stiffness-matrix kernels.

All functions are stateless and return matrices in local coordinates.
"""

import numpy as np


def shaft_stiffness(
    *,
    length,
    outer_diameter,
    inner_diameter,
    young_modulus,
    poisson_ratio,
    theory="timoshenko",
):
    """Build the local stiffness matrix of a circular shaft element.

    Parameters
    ----------
    length : float
        Element length.
    outer_diameter : float
        Shaft outer diameter.
    inner_diameter : float
        Shaft bore diameter. Use zero for a solid shaft.
    young_modulus : float
        Material Young's modulus.
    poisson_ratio : float
        Material Poisson's ratio.
    theory : {"timoshenko", "euler"}, optional
        Beam theory used for the bending terms.

    Returns
    -------
    numpy.ndarray
        Local ``(12, 12)`` stiffness matrix.

    Raises
    ------
    ValueError
        If a parameter is invalid or the beam theory is unknown.
    """
    # Collect scalar inputs so they can be validated together.
    values = np.asarray(
        [
            length,
            outer_diameter,
            inner_diameter,
            young_modulus,
            poisson_ratio,
        ],
        dtype=float,
    )

    # Matrix assembly cannot continue safely with NaN or infinite values.
    if not np.isfinite(values).all():
        raise ValueError("Shaft stiffness parameters must be finite")

    # A valid annular section requires 0 <= inner diameter < outer diameter.
    if length <= 0 or outer_diameter <= inner_diameter or inner_diameter < 0:
        raise ValueError("Shaft length and wall dimensions must be positive")

    # Linear elastic stiffness requires a positive Young's modulus.
    if young_modulus <= 0:
        raise ValueError("Young's modulus must be positive")

    # This is the admissible interval for an isotropic elastic solid.
    if not -1.0 < poisson_ratio < 0.5:
        raise ValueError("Poisson's ratio must be between -1 and 0.5")

    # Derive the isotropic shear modulus from E and nu.
    shear_modulus = young_modulus / (2.0 * (1.0 + poisson_ratio))

    # Compute the annular-section geometric properties.
    area = np.pi * (outer_diameter**2 - inner_diameter**2) / 4.0
    diametral_inertia = (
        np.pi * (outer_diameter**4 - inner_diameter**4) / 64.0
    )
    polar_inertia = 2.0 * diametral_inertia

    # Timoshenko theory modifies bending through the shear parameter phi.
    if theory == "timoshenko":
        diameter_ratio = inner_diameter / outer_diameter
        ratio_term = (1.0 + diameter_ratio**2) ** 2

        # Cowper's shear factor includes the effect of an annular bore.
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
        # Euler-Bernoulli theory neglects transverse shear deformation.
        shear_parameter = 0.0
    else:
        raise ValueError(f"Unknown shaft theory {theory!r}")

    # Allocate the full matrix in the package's 12-DOF ordering.
    matrix = np.zeros((12, 12), dtype=float)

    # Assemble the two-node axial block at local z translations.
    axial_block = (
        young_modulus
        * area
        / length
        * np.array([[1.0, -1.0], [-1.0, 1.0]])
    )
    matrix[np.ix_([2, 8], [2, 8])] = axial_block

    # Assemble the two-node torsional block at local z rotations.
    torsional_block = (
        shear_modulus
        * polar_inertia
        / length
        * np.array([[1.0, -1.0], [-1.0, 1.0]])
    )
    matrix[np.ix_([5, 11], [5, 11])] = torsional_block

    # Compute the common multiplier for both bending planes.
    bending_scale = (
        young_modulus
        * diametral_inertia
        / (length**3 * (1.0 + shear_parameter))
    )

    # Name the standard beam coefficients used in both bending blocks.
    k1 = 12.0
    k2 = 6.0 * length
    k3 = (4.0 + shear_parameter) * length**2
    k4 = (2.0 - shear_parameter) * length**2

    # Bending in the local x-z plane couples x translation and y rotation.
    bending_x = bending_scale * np.array(
        [
            [k1, k2, -k1, k2],
            [k2, k3, -k2, k4],
            [-k1, -k2, k1, -k2],
            [k2, k4, -k2, k3],
        ]
    )

    # The y-z block has opposite coupling signs because dy/dz = -theta_x.
    bending_y = bending_scale * np.array(
        [
            [k1, -k2, -k1, -k2],
            [-k2, k3, k2, k4],
            [-k1, k2, k1, k2],
            [-k2, k4, k2, k3],
        ]
    )

    # Insert both bending blocks into the complete local matrix.
    matrix[np.ix_([0, 4, 6, 10], [0, 4, 6, 10])] = bending_x
    matrix[np.ix_([1, 3, 7, 9], [1, 3, 7, 9])] = bending_y

    return matrix


def bearing_stiffness(*, kxx=0.0, kyy=0.0, kzz=0.0, kxy=0.0, kyx=0.0):
    """Build the local translational stiffness matrix of a bearing.

    Parameters
    ----------
    kxx, kyy, kzz : float, optional
        Direct translational stiffness coefficients.
    kxy, kyx : float, optional
        Cross-coupled translational stiffness coefficients.

    Returns
    -------
    numpy.ndarray
        Local ``(6, 6)`` bearing stiffness matrix.

    Raises
    ------
    ValueError
        If any coefficient is not finite.
    """
    # Validate all bearing coefficients before inserting them.
    coefficients = [kxx, kyy, kzz, kxy, kyx]
    if not np.isfinite(coefficients).all():
        raise ValueError("Bearing stiffness coefficients must be finite")

    # Start from zero because this bearing model has no rotational stiffness.
    matrix = np.zeros((6, 6), dtype=float)

    # Insert direct and cross-coupled translational coefficients.
    matrix[:3, :3] = [
        [kxx, kxy, 0.0],
        [kyx, kyy, 0.0],
        [0.0, 0.0, kzz],
    ]

    return matrix


def disk_stiffness():
    """Build the local stiffness matrix of a rigid disk.

    Returns
    -------
    numpy.ndarray
        Zero ``(6, 6)`` matrix because a rigid lumped disk stores no elastic
        strain energy.
    """
    # A lumped rigid disk contributes inertia, but no elastic stiffness.
    return np.zeros((6, 6), dtype=float)
