"""Shared fixtures and independent analytical helpers for Spinniped tests."""

from dataclasses import dataclass

import numpy as np
import pytest

from spinniped import element, node, rotor, solver


# This immutable container collects the material and geometric data used by
# several tests. Keeping one standard set of data prevents small differences
# between tests from obscuring the physical quantity that is being verified.
@dataclass(frozen=True)
class ShaftProperties:
    # Young's modulus [Pa].
    E: float = 2.0e11

    # Poisson's ratio [-].
    nu: float = 0.3

    # Material density [kg/m^3].
    rho: float = 7850.0

    # Diameter of the filled circular shaft section [m].
    d: float = 0.01

    # Total shaft length used by the benchmark model [m].
    L: float = 1.0

    @property
    def A(self):
        # Cross-sectional area of a filled circular section [m^2].
        return np.pi * self.d**2 / 4.0

    @property
    def I(self):
        # Second moment of area about either transverse centroidal axis [m^4].
        # For a circular section, Ix and Iy have the same value.
        return np.pi * self.d**4 / 64.0


# A pytest fixture is created automatically whenever a test function declares
# an argument named ``shaft_properties``. By default, pytest calls this
# function separately for every test, so tests remain independent.
@pytest.fixture
def shaft_properties():
    return ShaftProperties()


# This fixture supplies one isolated two-node shaft element. It is mainly used
# by the local-matrix and one-element assembly tests.
@pytest.fixture
def shaft_element(shaft_properties):
    # ``p`` is a short alias for the standard shaft properties.
    p = shaft_properties

    # Node identifiers 1 and 2 are sufficient here because this fixture tests
    # the element itself; their physical coordinates are not needed until the
    # element is inserted into a mesh.
    return element.ShaftElement(
        eid=1,
        E=p.E,
        nu=p.nu,
        rho=p.rho,
        d=p.d,
        L=p.L,
        n1=1,
        n2=2,
    )


def build_uniform_beam(p, ne, theory="euler"):
    """Build and assemble a straight, uniform beam along the global z axis."""
    # A chain of ``ne`` two-node elements requires one more node than elements.
    nn = ne + 1

    # Nodes are equally spaced from z = 0 to z = L. The x and y coordinates
    # retain their default value of zero, making a straight axial mesh.
    meshgrid = node.Grid()
    for i in range(nn):
        meshgrid.add_node(z=i * p.L / ne)

    # The Rotor object groups the mesh and the finite elements before they are
    # passed to the solver.
    r = rotor.Rotor()

    # All elements have the same length because the mesh is uniform.
    Le = p.L / ne

    # Connect each node to its immediate neighbour, producing the topology:
    #
    #     node 1 -- element 1 -- node 2 -- element 2 -- ... -- node nn
    for eid in range(ne):
        e = element.ShaftElement(
            eid=eid + 1,
            E=p.E,
            nu=p.nu,
            rho=p.rho,
            d=p.d,
            L=Le,
            n1=int(meshgrid.Nodes["id"][eid]),
            n2=int(meshgrid.Nodes["id"][eid + 1]),
        )

        if theory == "euler":
            # The public constructor currently has no theory argument.
            # Therefore, replace the constructor's default Timoshenko matrices
            # with Euler-Bernoulli matrices explicitly.
            #
            # Rotary inertia is disabled so that the numerical model matches
            # the classical Euler-Bernoulli analytical frequency formula used
            # by the benchmark tests.
            e.K = e._K(theory=1)
            e.M = e._M(theory=1, rotary_inertia=False)
        elif theory != "timoshenko":
            # Reject misspelled or unsupported theory names early; otherwise a
            # benchmark could silently test a different formulation.
            raise ValueError(f"Unknown beam theory: {theory}")

        # Add the completed element to the rotor assembly.
        r.add_element(e)

    # Attach the mesh only after its nodes and the element connectivity have
    # been fully defined.
    r.add_meshgrid(meshgrid)

    # The Solver copies the rotor's elements and mesh, then assembles the local
    # element matrices into full global K and M matrices.
    s = solver.Solver()
    s.add_rotor(r)
    s.assemble_global_matrices()

    # Boundary conditions are intentionally not applied here. Each benchmark
    # applies the supports appropriate to the physical problem it tests.
    return s


def simply_supported_frequencies(p, nm):
    """Euler-Bernoulli frequencies for a uniform simply supported beam."""
    # ``m`` contains the one-based bending mode numbers: 1, 2, ..., nm.
    m = np.arange(1, nm + 1, dtype=float)

    # Analytical solution:
    #
    # f_m = (m^2*pi)/(2*L^2) * sqrt(E*I/(rho*A))
    #
    # This helper is independent of Spinniped's element matrices. That
    # independence is important: a benchmark should not reproduce the same
    # implementation that it is intended to verify.
    return m**2 * np.pi / (2.0 * p.L**2) * np.sqrt(
        p.E * p.I / (p.rho * p.A)
    )


def positive_frequencies(eigenvalues, rtol=1.0e-10):
    """Convert positive eigenvalues to frequencies while ignoring rigid modes."""
    # Ensure slicing and vectorized arithmetic work even when the caller
    # provides a Python list.
    eigenvalues = np.asarray(eigenvalues)

    # The tolerance is relative to the largest eigenvalue. This is more robust
    # than comparing with one fixed absolute number because eigenvalue units
    # and magnitudes vary strongly with geometry and material properties.
    scale = max(float(np.max(np.abs(eigenvalues))), 1.0)

    # Remove zero or numerically tiny eigenvalues associated with rigid-body
    # motion. Small negative values caused by round-off are removed as well.
    positive_eigenvalues = eigenvalues[eigenvalues > rtol * scale]

    # For K*q = lambda*M*q, lambda = omega^2. Convert angular frequency
    # omega [rad/s] into cyclic frequency f [Hz].
    return np.sqrt(positive_eigenvalues) / (2.0 * np.pi)


def one_frequency_per_bending_pair(frequencies, nm):
    """Select one frequency from each x/y pair of an axisymmetric shaft."""
    # A circular shaft has equal bending stiffness in x and y. Consequently,
    # every bending frequency normally appears twice in the eigensolution.
    # Taking indices 0, 2, 4, ... keeps one representative from each pair.
    return np.asarray(frequencies)[0 : 2 * nm : 2]
