"""
Benchmarking the natural frequencies of a simple shaft against analytical solutions for a cantilever beam. 
"""
import numpy as np

def simply_supported_beam(E, I, L, rho, A, m):
    """Analytical natural frequency for an m-th mode simply-supported beam.

    Parameters
    ----------
    E : float
        Young's modulus (Pa).
    I : float
        Second moment of area (m^4).
    L : float
        Beam length (m).
    rho : float
        Material density (kg/m^3).
    A : float
        Cross-sectional area (m^2).
    m : int
        Mode number (1-based).

    Returns
    -------
    f : float
        Natural frequency in Hz for the specified mode.
    """
    k = {1: 9.87, 2: 39.48, 3: 88.82, 4: 156.96}  # First four mode constants for simply supported beam
    f = k[m]/(2*np.pi) * np.sqrt(E*I/(rho*A*L**4))
    return f

def fix_ends_beam(E, I, L, rho, A, m):
    """Analytical natural frequency for an m-th mode fixed-ends beam.

    Parameters
    ----------
    E, I, L, rho, A, m : see ``simply_supported_beam``.

    Returns
    -------
    f : float
        Natural frequency in Hz for the specified mode.
    """
    k = {1: 22.4, 2: 61.7, 3: 120.9, 4: 200.1}  # First four mode constants for fixed ends beam
    f = k[m]/(2*np.pi) * np.sqrt(E*I/(rho*A*L**4))
    return f


for m in range(1,5):

    E = 2.0e11
    nu = 0.3
    rho = 7850
    d = 0.01
    L = 1.0

    A = np.pi / 4.0 * d**2
    I = np.pi / 64.0 * d**4

    f_simply_supported = simply_supported_beam(E, I, L, rho, A, m)
    f_fixed_ends = fix_ends_beam(E, I, L, rho, A, m)

    print(f"Mode {m}: Simply Supported Frequency = {f_simply_supported:.2f} Hz, Fixed Ends Frequency = {f_fixed_ends:.2f} Hz")