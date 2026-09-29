"""Print analytical bending frequencies for a uniform reference beam."""

import numpy as np


def simply_supported_beam(
    young_modulus,
    second_moment,
    length,
    density,
    area,
    mode,
):
    """Return a simply supported Euler-Bernoulli bending frequency.

    Parameters
    ----------
    young_modulus : float
        Young's modulus in pascals.
    second_moment : float
        Section second moment of area in metres to the fourth power.
    length : float
        Beam length in metres.
    density : float
        Material density in kilograms per cubic metre.
    area : float
        Cross-sectional area in square metres.
    mode : int
        One-based bending mode number from one through four.

    Returns
    -------
    float
        Natural frequency in hertz.
    """
    # A simply supported beam has the exact characteristic value beta*L=n*pi.
    characteristic_value = mode * np.pi

    # Apply the closed-form Euler-Bernoulli frequency equation.
    return (
        characteristic_value**2
        / (2.0 * np.pi)
        * np.sqrt(
            young_modulus
            * second_moment
            / (density * area * length**4)
        )
    )


def fixed_ends_beam(
    young_modulus,
    second_moment,
    length,
    density,
    area,
    mode,
):
    """Return a fixed-fixed Euler-Bernoulli bending frequency.

    Parameters
    ----------
    young_modulus : float
        Young's modulus in pascals.
    second_moment : float
        Section second moment of area in metres to the fourth power.
    length : float
        Beam length in metres.
    density : float
        Material density in kilograms per cubic metre.
    area : float
        Cross-sectional area in square metres.
    mode : int
        One-based bending mode number from one through four.

    Returns
    -------
    float
        Natural frequency in hertz.
    """
    # These are the first roots of cos(beta*L) cosh(beta*L) = 1.
    characteristic_roots = {
        1: 4.730040744862704,
        2: 7.853204624095838,
        3: 10.99560783800167,
        4: 14.137165491257,
    }

    # Apply the corresponding closed-form beam frequency equation.
    return (
        characteristic_roots[mode] ** 2
        / (2.0 * np.pi)
        * np.sqrt(
            young_modulus
            * second_moment
            / (density * area * length**4)
        )
    )


# Define one uniform circular steel-like reference beam.
YOUNG_MODULUS = 2.0e11
DENSITY = 7850.0
DIAMETER = 0.01
LENGTH = 1.0

# Compute its circular-section geometric properties.
AREA = np.pi * DIAMETER**2 / 4.0
SECOND_MOMENT = np.pi * DIAMETER**4 / 64.0

# Print the first four frequencies for both support conditions.
for mode_number in range(1, 5):
    simply_supported_frequency = simply_supported_beam(
        YOUNG_MODULUS,
        SECOND_MOMENT,
        LENGTH,
        DENSITY,
        AREA,
        mode_number,
    )
    fixed_ends_frequency = fixed_ends_beam(
        YOUNG_MODULUS,
        SECOND_MOMENT,
        LENGTH,
        DENSITY,
        AREA,
        mode_number,
    )

    print(
        f"Mode {mode_number}: "
        f"Simply Supported Frequency = {simply_supported_frequency:.2f} Hz, "
        f"Fixed Ends Frequency = {fixed_ends_frequency:.2f} Hz"
    )
