"""Element viscous-damping-matrix kernels."""

import numpy as np


def shaft_damping(*, damping=0.0, mass_matrix=None):
    """Build a mass-proportional shaft damping matrix.

    Parameters
    ----------
    damping : float, optional
        Mass-proportional damping coefficient.
    mass_matrix : array_like, optional
        Local ``(12, 12)`` shaft mass matrix. It is required when ``damping``
        is nonzero.

    Returns
    -------
    numpy.ndarray
        Local ``(12, 12)`` damping matrix.

    Raises
    ------
    ValueError
        If the coefficient or supplied mass matrix is invalid.
    """
    # Reject invalid coefficients before multiplying a physical matrix.
    if not np.isfinite(damping):
        raise ValueError("Shaft damping must be finite")

    # A zero coefficient can produce a zero matrix without a mass matrix.
    if mass_matrix is None:
        if damping:
            raise ValueError("mass_matrix is required when shaft damping is nonzero")
        return np.zeros((12, 12), dtype=float)

    # Convert and validate the caller's matrix explicitly.
    mass_matrix = np.asarray(mass_matrix, dtype=float)
    if mass_matrix.shape != (12, 12) or not np.isfinite(mass_matrix).all():
        raise ValueError("shaft mass_matrix must be a finite 12 by 12 matrix")

    # Apply simple mass-proportional damping.
    return damping * mass_matrix


def bearing_damping(*, cxx=0.0, cyy=0.0, czz=0.0, cxy=0.0, cyx=0.0):
    """Build the local translational damping matrix of a bearing.

    Parameters
    ----------
    cxx, cyy, czz : float, optional
        Direct viscous damping coefficients.
    cxy, cyx : float, optional
        Cross-coupled viscous damping coefficients.

    Returns
    -------
    numpy.ndarray
        Local ``(6, 6)`` bearing damping matrix.

    Raises
    ------
    ValueError
        If any coefficient is not finite.
    """
    # Validate all coefficients before inserting them into the matrix.
    coefficients = [cxx, cyy, czz, cxy, cyx]
    if not np.isfinite(coefficients).all():
        raise ValueError("Bearing damping coefficients must be finite")

    # Bearings currently damp translations but not nodal rotations.
    matrix = np.zeros((6, 6), dtype=float)
    matrix[:3, :3] = [
        [cxx, cxy, 0.0],
        [cyx, cyy, 0.0],
        [0.0, 0.0, czz],
    ]

    return matrix


def disk_damping(*, damping=0.0, mass_matrix=None):
    """Build a mass-proportional disk damping matrix.

    Parameters
    ----------
    damping : float, optional
        Mass-proportional damping coefficient.
    mass_matrix : array_like, optional
        Local ``(6, 6)`` disk mass matrix. It is required when ``damping`` is
        nonzero.

    Returns
    -------
    numpy.ndarray
        Local ``(6, 6)`` damping matrix.

    Raises
    ------
    ValueError
        If the coefficient or supplied mass matrix is invalid.
    """
    # Validate the scalar damping coefficient first.
    if not np.isfinite(damping):
        raise ValueError("Disk damping must be finite")

    # A zero coefficient can produce a zero matrix without a mass matrix.
    if mass_matrix is None:
        if damping:
            raise ValueError("mass_matrix is required when disk damping is nonzero")
        return np.zeros((6, 6), dtype=float)

    # Convert and validate the supplied lumped mass matrix.
    mass_matrix = np.asarray(mass_matrix, dtype=float)
    if mass_matrix.shape != (6, 6) or not np.isfinite(mass_matrix).all():
        raise ValueError("disk mass_matrix must be a finite 6 by 6 matrix")

    # Apply simple mass-proportional damping.
    return damping * mass_matrix
