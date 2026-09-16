"""Element mass-matrix kernels."""

import numpy as np


def shaft_mass(
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
    """Build the consistent local mass matrix of a circular shaft.

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
        Material Young's modulus, used by the Timoshenko shape functions.
    poisson_ratio : float
        Material Poisson's ratio.
    theory : {"timoshenko", "euler"}, optional
        Beam theory used for the consistent bending terms.
    rotary_inertia : bool, optional
        Include bending rotary-inertia terms when true.

    Returns
    -------
    numpy.ndarray
        Local ``(12, 12)`` mass matrix.

    Raises
    ------
    TypeError
        If ``rotary_inertia`` is not Boolean.
    ValueError
        If a parameter is invalid or the beam theory is unknown.
    """
    # Gather scalar parameters for one finite-value check.
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
        raise ValueError("Shaft mass parameters must be finite")

    # Keep the inertia switch explicit instead of accepting truthy values.
    if not isinstance(rotary_inertia, (bool, np.bool_)):
        raise TypeError("rotary_inertia must be a boolean")

    # Validate the annular section and its material.
    if length <= 0 or outer_diameter <= inner_diameter or inner_diameter < 0:
        raise ValueError("Shaft length and wall dimensions must be positive")
    if density <= 0:
        raise ValueError("Density must be positive")
    if young_modulus <= 0:
        raise ValueError("Young's modulus must be positive")
    if not -1.0 < poisson_ratio < 0.5:
        raise ValueError("Poisson's ratio must be between -1 and 0.5")

    # Compute annular-section properties used throughout the matrix.
    area = np.pi * (outer_diameter**2 - inner_diameter**2) / 4.0
    diametral_inertia = (
        np.pi * (outer_diameter**4 - inner_diameter**4) / 64.0
    )
    polar_inertia = 2.0 * diametral_inertia

    # Evaluate the shear parameter required by Timoshenko interpolation.
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
        # Euler-Bernoulli interpolation is recovered when phi is zero.
        shear_parameter = 0.0
    else:
        raise ValueError(f"Unknown shaft theory {theory!r}")

    # Allocate the complete matrix in the standard 12-DOF ordering.
    matrix = np.zeros((12, 12), dtype=float)

    # Insert consistent axial translational inertia.
    axial_block = (
        density
        * area
        * length
        / 6.0
        * np.array([[2.0, 1.0], [1.0, 2.0]])
    )
    matrix[np.ix_([2, 8], [2, 8])] = axial_block

    # Insert consistent torsional inertia about the shaft axis.
    torsional_block = (
        density
        * polar_inertia
        * length
        / 6.0
        * np.array([[2.0, 1.0], [1.0, 2.0]])
    )
    matrix[np.ix_([5, 11], [5, 11])] = torsional_block

    # Select the consistent translational bending coefficients.
    if theory == "timoshenko":
        phi = shear_parameter
        m1 = 312.0 + 588.0 * phi + 280.0 * phi**2
        m2 = (44.0 + 77.0 * phi + 35.0 * phi**2) * length
        m3 = 108.0 + 252.0 * phi + 140.0 * phi**2
        m4 = -(26.0 + 63.0 * phi + 35.0 * phi**2) * length
        m5 = (8.0 + 14.0 * phi + 7.0 * phi**2) * length**2
        m6 = -(6.0 + 14.0 * phi + 7.0 * phi**2) * length**2
        bending_scale = (
            density * area * length / (840.0 * (1.0 + phi) ** 2)
        )
    else:
        m1 = 156.0
        m2 = 22.0 * length
        m3 = 54.0
        m4 = -13.0 * length
        m5 = 4.0 * length**2
        m6 = -3.0 * length**2
        bending_scale = density * area * length / 420.0

    # Build the translational bending block for the local x-z plane.
    bending_x = bending_scale * np.array(
        [
            [m1, m2, m3, m4],
            [m2, m5, -m4, m6],
            [m3, -m4, m1, -m2],
            [m4, m6, -m2, m5],
        ]
    )

    # Apply the package sign convention in the local y-z plane.
    bending_y = bending_scale * np.array(
        [
            [m1, -m2, m3, -m4],
            [-m2, m5, m4, m6],
            [m3, m4, m1, m2],
            [-m4, m6, m2, m5],
        ]
    )

    # Insert translational bending inertia into both lateral planes.
    x_indices = [0, 4, 6, 10]
    y_indices = [1, 3, 7, 9]
    matrix[np.ix_(x_indices, x_indices)] = bending_x
    matrix[np.ix_(y_indices, y_indices)] = bending_y

    # Stop here when bending rotary inertia has been disabled.
    if not rotary_inertia:
        return matrix

    # Compute the common scale for the rotary-inertia contribution.
    denominator = (
        (1.0 + shear_parameter) ** 2 if theory == "timoshenko" else 1.0
    )
    rotary_scale = density * diametral_inertia / (30.0 * length * denominator)

    # Define rotary coefficients for the selected beam interpolation.
    r7 = 36.0
    r8 = (3.0 - 15.0 * shear_parameter) * length
    r9 = (
        4.0 + 5.0 * shear_parameter + 10.0 * shear_parameter**2
    ) * length**2
    r10 = (
        -1.0 - 5.0 * shear_parameter + 5.0 * shear_parameter**2
    ) * length**2

    # Assemble rotary inertia in the x-z bending plane.
    rotary_x = rotary_scale * np.array(
        [
            [r7, r8, -r7, r8],
            [r8, r9, -r8, r10],
            [-r7, -r8, r7, -r8],
            [r8, r10, -r8, r9],
        ]
    )

    # Assemble rotary inertia in the y-z bending plane.
    rotary_y = rotary_scale * np.array(
        [
            [r7, -r8, -r7, -r8],
            [-r8, r9, r8, r10],
            [-r7, r8, r7, r8],
            [-r8, r10, r8, r9],
        ]
    )

    # Add rotary terms to the translational bending terms already present.
    matrix[np.ix_(x_indices, x_indices)] += rotary_x
    matrix[np.ix_(y_indices, y_indices)] += rotary_y

    return matrix


def bearing_mass():
    """Build the local mass matrix of a massless bearing.

    Returns
    -------
    numpy.ndarray
        Zero ``(6, 6)`` mass matrix.
    """
    # The current bearing model represents support coefficients only.
    return np.zeros((6, 6), dtype=float)


def disk_mass(*, mass, diametral_inertia=0.0, polar_inertia=0.0):
    """Build the local lumped mass matrix of a rigid axisymmetric disk.

    Parameters
    ----------
    mass : float
        Disk translational mass.
    diametral_inertia : float, optional
        Mass moment of inertia about either local transverse axis.
    polar_inertia : float, optional
        Mass moment of inertia about the local spin axis.

    Returns
    -------
    numpy.ndarray
        Diagonal ``(6, 6)`` lumped mass matrix.

    Raises
    ------
    ValueError
        If any value is non-finite or negative.
    """
    # Validate the three physical inertia parameters.
    values = [mass, diametral_inertia, polar_inertia]
    if not np.isfinite(values).all():
        raise ValueError("Disk mass and inertias must be finite")
    if mass < 0 or diametral_inertia < 0 or polar_inertia < 0:
        raise ValueError("Disk mass and inertias cannot be negative")

    # Translational mass repeats on x, y, z; diametral inertia repeats on rx, ry.
    diagonal = [
        mass,
        mass,
        mass,
        diametral_inertia,
        diametral_inertia,
        polar_inertia,
    ]
    return np.diag(diagonal).astype(float)
