"""Shared declarative-model fixtures and independent analytical helpers."""

from dataclasses import dataclass

import numpy as np
import pytest

from spinniped import (
    Grid,
    Material,
    ModelBuilder,
    ModelDefinition,
    ShaftElement,
    ShaftProperty,
)
from spinniped.gyroscopic import shaft_gyroscopic
from spinniped.mass import shaft_mass
from spinniped.stiffness import shaft_stiffness


@dataclass(frozen=True)
class ShaftProperties:
    """Reference SI properties used throughout the beam tests."""

    young_modulus: float = 2.0e11
    poisson_ratio: float = 0.3
    density: float = 7850.0
    diameter: float = 0.01
    length: float = 1.0

    @property
    def area(self):
        """Return the circular cross-sectional area.

        Returns
        -------
        float
            Cross-sectional area.
        """
        return np.pi * self.diameter**2 / 4.0

    @property
    def second_moment(self):
        """Return the transverse second moment of area.

        Returns
        -------
        float
            Second moment about either transverse centroidal axis.
        """
        return np.pi * self.diameter**4 / 64.0

    @property
    def polar_moment(self):
        """Return the polar second moment of area.

        Returns
        -------
        float
            Polar second moment about the shaft axis.
        """
        return 2.0 * self.second_moment

    @property
    def shear_modulus(self):
        """Return the isotropic shear modulus.

        Returns
        -------
        float
            Shear modulus derived from Young's modulus and Poisson's ratio.
        """
        return self.young_modulus / (2.0 * (1.0 + self.poisson_ratio))


@pytest.fixture
def shaft_properties():
    """Return the shared reference shaft properties fixture.

    Returns
    -------
    ShaftProperties
        Reference SI properties used by beam and shaft tests.
    """
    return ShaftProperties()


@pytest.fixture
def shaft_definition(shaft_properties):
    """One declarative shaft record joining two grids."""
    p = shaft_properties
    return ModelDefinition(
        grids=[Grid(id=1, z=0.0), Grid(id=2, z=p.length)],
        materials=[
            Material(
                id=1,
                density=p.density,
                young_modulus=p.young_modulus,
                poisson_ratio=p.poisson_ratio,
            )
        ],
        properties=[
            ShaftProperty(
                id=1,
                material=1,
                outer_diameter=p.diameter,
                theory="timoshenko",
                rotary_inertia=True,
            )
        ],
        elements=[ShaftElement(id=1, grid_a=1, grid_b=2, property=1)],
    )


@pytest.fixture
def shaft_local_matrices(shaft_properties):
    """Local numerical kernels evaluated from the reference Python values."""
    p = shaft_properties
    common = {
        "length": p.length,
        "outer_diameter": p.diameter,
        "inner_diameter": 0.0,
        "young_modulus": p.young_modulus,
        "poisson_ratio": p.poisson_ratio,
        "theory": "timoshenko",
    }
    mass_options = {
        **common,
        "density": p.density,
        "rotary_inertia": True,
    }
    return {
        "stiffness": shaft_stiffness(**common),
        "mass": shaft_mass(**mass_options),
        "gyroscopic": shaft_gyroscopic(**mass_options),
    }


def build_uniform_beam(properties, elements, theory="euler"):
    """Build a straight uniform beam along global z using declarative records."""
    if theory not in {"euler", "timoshenko"}:
        raise ValueError(f"Unknown beam theory: {theory}")
    if not isinstance(elements, int) or elements < 1:
        raise ValueError("elements must be a positive integer")

    p = properties
    grids = [
        Grid(id=index + 1, z=index * p.length / elements)
        for index in range(elements + 1)
    ]
    definition = ModelDefinition(
        grids=grids,
        materials=[
            Material(
                id=1,
                density=p.density,
                young_modulus=p.young_modulus,
                poisson_ratio=p.poisson_ratio,
            )
        ],
        properties=[
            ShaftProperty(
                id=1,
                material=1,
                outer_diameter=p.diameter,
                theory=theory,
                # The analytical reference neglects bending rotary inertia.
                rotary_inertia=False,
            )
        ],
        elements=[
            ShaftElement(
                id=index + 1,
                grid_a=index + 1,
                grid_b=index + 2,
                property=1,
            )
            for index in range(elements)
        ],
    )
    return ModelBuilder().build(definition)


def simply_supported_fixed_dofs(elements):
    """Return global support DOFs for the six-DOF grid ordering."""
    final_grid_offset = 6 * elements
    return [0, 1, 2, 5, final_grid_offset, final_grid_offset + 1]


def simply_supported_frequencies(properties, modes):
    """Independent Euler--Bernoulli frequencies of a uniform simple beam."""
    indices = np.arange(1, modes + 1, dtype=float)
    return indices**2 * np.pi / (2.0 * properties.length**2) * np.sqrt(
        properties.young_modulus
        * properties.second_moment
        / (properties.density * properties.area)
    )


def one_frequency_per_bending_pair(frequencies, modes):
    """Select one member of each degenerate x/y bending-mode pair."""
    return np.asarray(frequencies)[0 : 2 * modes : 2]
