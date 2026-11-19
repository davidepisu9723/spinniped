"""Partial extraction accuracy, root budgets, and mode correspondence."""

from dataclasses import replace

import numpy as np
import pytest
from scipy.linalg import eig, eigh

from spinniped import ModelBuilder, Solver
from conftest import build_uniform_beam


@pytest.fixture
def eigenvalue_model(shaft_properties):
    """A supported shaft with nonzero gyroscopic and damping matrices."""
    definition = build_uniform_beam(shaft_properties, 10).definitions[0]
    definition = replace(
        definition,
        properties=[replace(
            definition.properties[0], theory="timoshenko",
            rotary_inertia=True, damping=0.2,
        )],
    )
    return ModelBuilder().build(definition)


def dense_state_roots(stiffness, mass, damping):
    """Independent full state solution used only as a numerical reference."""
    size = len(mass)
    state = np.block([
        [np.zeros_like(mass), np.eye(size)],
        [-np.linalg.solve(mass, stiffness), -np.linalg.solve(mass, damping)],
    ])
    values, _ = eig(state)
    values = values[values.imag > 1.0e-8]
    return values[np.argsort(values.imag)]


class TestEigenvalueExtraction:
    """Verify root budgets, accuracy, and direct failure reporting."""

    @pytest.mark.parametrize("modes", [1, 6, None])
    def test_modal_matches_full_solution(self, eigenvalue_model, modes):
        model, free, _ = Solver._apply_boundary_conditions(eigenvalue_model)
        stiffness = model.K[0][np.ix_(free, free)]
        mass = model.M[0][np.ix_(free, free)]
        expected, _ = eigh(stiffness, mass)
        result = Solver(eigenvalue_model).solve("modal", modes=modes)
        count = len(free) if modes is None else modes
        np.testing.assert_allclose(
            result["eigenvalues"][0], expected[:count], rtol=2.0e-7,
        )
        vectors = result["eigenvectors"][0]
        np.testing.assert_allclose(
            vectors.T @ mass @ vectors, np.eye(count), atol=1.0e-10,
        )
        residual = stiffness @ vectors - (
            mass @ vectors * result["eigenvalues"][0]
        )
        assert np.linalg.norm(residual) / (
            np.linalg.norm(stiffness) * np.linalg.norm(vectors)
        ) < 1.0e-12

    def test_modal_only_requests_a_subset(self, eigenvalue_model, monkeypatch):
        import scipy.linalg

        extracted = []
        original = scipy.linalg.eigh

        def record_subset(*args, **kwargs):
            values, vectors = original(*args, **kwargs)
            extracted.append(len(values))
            return values, vectors

        monkeypatch.setattr(scipy.linalg, "eigh", record_subset)
        Solver(eigenvalue_model).solve("modal", modes=6)
        assert extracted == [6]

    @pytest.mark.parametrize("speed", [0.0, 0.1, 200.0, 600.0])
    def test_campbell_matches_full_solution(self, eigenvalue_model, speed):
        model, free, _ = Solver._apply_boundary_conditions(eigenvalue_model)
        stiffness, mass, damping, gyro = (
            matrix[0][np.ix_(free, free)]
            for matrix in (model.K, model.M, model.C, model.G)
        )
        effective_damping = damping + speed * gyro
        expected = dense_state_roots(stiffness, mass, effective_damping)[:6]
        result = Solver(eigenvalue_model).solve(
            "campbell", speeds=[speed], modes=6,
        )
        values = result["eigenvalues"][0, 0]
        np.testing.assert_allclose(values, expected, rtol=2.0e-7, atol=1.0e-6)
        vectors = result["eigenvectors"][0, 0]
        displacement, velocity = np.split(vectors, 2)
        np.testing.assert_allclose(
            velocity, displacement * values, rtol=1.0e-7, atol=1.0e-9,
        )
        residual = (
            stiffness @ displacement
            + (effective_damping @ displacement) * values
            + (mass @ displacement) * values**2
        )
        scale = (
            np.linalg.norm(stiffness)
            + np.max(np.abs(values)) * np.linalg.norm(effective_damping)
            + np.max(np.abs(values))**2 * np.linalg.norm(mass)
        ) * np.linalg.norm(displacement)
        assert np.linalg.norm(residual) / scale < 1.0e-10

    def test_campbell_avoids_dense_state_solve(self, eigenvalue_model, monkeypatch):
        import scipy.linalg
        import scipy.sparse.linalg

        roots = []
        original = scipy.sparse.linalg.eigs

        def record_roots(operator, **kwargs):
            roots.append((operator.shape[0], kwargs["k"]))
            return original(operator, **kwargs)

        def unexpected_dense(*args, **kwargs):
            raise AssertionError("requested six modes triggered a full solve")

        monkeypatch.setattr(scipy.sparse.linalg, "eigs", record_roots)
        monkeypatch.setattr(scipy.linalg, "eig", unexpected_dense)
        Solver(eigenvalue_model).solve("campbell", speeds=[200.0], modes=6)
        assert roots == [(120, 14)]

    @pytest.mark.parametrize("analysis", ["modal", "campbell"])
    def test_tracking_reorders_the_extracted_modes(self, analysis):
        mass = np.eye(24)
        stiffness = np.diag(np.arange(1.0, 25.0)**2)
        reference = np.eye(24)[:, [2, 0, 1]]
        if analysis == "modal":
            values, vectors = Solver._modal_eigenpairs(
                stiffness, mass, 3,
            )
            expected = [9.0, 1.0, 4.0]
        else:
            values, vectors = Solver._campbell_eigenpairs(
                stiffness, mass, np.zeros_like(mass), 3,
            )
            vectors = vectors[:24]
            expected = [3.0j, 1.0j, 2.0j]
        order = Solver._mode_order(reference, vectors)
        np.testing.assert_allclose(values[order], expected, rtol=1.0e-10)

    def test_modal_skips_rigid_modes_and_honors_tolerance(self):
        mass = np.eye(24)
        stiffness = np.diag([0.0] * 10 + list(np.arange(1.0, 15.0)**2))
        values, _ = Solver._modal_eigenpairs(stiffness, mass, 3)
        np.testing.assert_allclose(values[:3], [1.0, 4.0, 9.0])
        values, _ = Solver._modal_eigenpairs(
            stiffness, mass, 3, zero_tolerance=4.0,
        )
        np.testing.assert_allclose(values[:3], [9.0, 16.0, 25.0])

    def test_modal_does_not_hide_a_negative_root(self):
        stiffness = np.diag([-1000.0] + list(np.arange(1.0, 24.0)))
        with pytest.raises(ValueError, match="negative eigenvalues"):
            Solver._modal_eigenpairs(stiffness, np.eye(24), 2)

    @pytest.mark.parametrize("modes", [None, 1, 2])
    def test_campbell_small_systems(self, modes):
        values, _ = Solver._campbell_eigenpairs(
            np.diag([4.0, 9.0]), np.eye(2), np.zeros((2, 2)), modes,
        )
        expected = [2.0j] if modes == 1 else [2.0j, 3.0j]
        np.testing.assert_allclose(values, expected)

    def test_campbell_all_modes_preserve_small_roots_in_a_wide_spectrum(self):
        frequencies = np.geomspace(1.0e-3, 1.0e6, 24)
        values, _ = Solver._campbell_eigenpairs(
            np.diag(frequencies**2), np.eye(24), np.zeros((24, 24)), None,
        )
        np.testing.assert_allclose(values, 1j * frequencies, rtol=1.0e-10)

    def test_campbell_singular_shift_raises(self):
        from scipy.linalg import LinAlgWarning

        stiffness = np.diag(np.arange(1.0, 25.0)**2)
        damping = np.diag([2.0] + [0.0] * 23)
        with pytest.raises(LinAlgWarning, match="[Ss]ingular"):
            Solver._campbell_eigenpairs(stiffness, np.eye(24), damping, 3)

    @pytest.mark.parametrize("damping", [0.0, 0.2])
    @pytest.mark.parametrize("free", [1, 2])
    def test_campbell_filters_zero_roots_without_expanding(self, damping, free):
        stiffness = np.diag([0.0] * free + list(np.arange(1.0, 25.0 - free)**2))
        mass = np.eye(24)
        values, _ = Solver._campbell_eigenpairs(
            stiffness, mass, damping * mass, 3,
        )
        # Eight state roots include two roots per free coordinate. The solver
        # returns the oscillating candidates without retrying for more roots.
        expected = -damping / 2 + 1j * np.sqrt(
            np.arange(1.0, 5.0 - free)**2 - (damping / 2)**2
        )
        np.testing.assert_allclose(values, expected, atol=1.0e-10)

    def test_campbell_retains_repeated_root_multiplicity(self):
        frequencies = np.repeat(np.arange(1.0, 9.0), 3)
        mass = np.eye(len(frequencies))
        values, _ = Solver._campbell_eigenpairs(
            np.diag(frequencies**2), mass, 0.2 * mass, 6,
        )
        expected = -0.1 + 1j * np.sqrt(frequencies[:6]**2 - 0.01)
        np.testing.assert_allclose(values[:6], expected, atol=1.0e-10)

    def test_campbell_propagates_nonconvergence(self, monkeypatch):
        import scipy.linalg
        import scipy.sparse.linalg

        def nonconvergent(*args, **kwargs):
            raise scipy.sparse.linalg.ArpackNoConvergence(
                "test", np.empty(0), np.empty((48, 0)),
            )

        def unexpected_dense(*args, **kwargs):
            raise AssertionError("iteration failure triggered a dense solve")

        monkeypatch.setattr(scipy.sparse.linalg, "eigs", nonconvergent)
        monkeypatch.setattr(scipy.linalg, "eig", unexpected_dense)
        with pytest.raises(scipy.sparse.linalg.ArpackNoConvergence):
            Solver._campbell_eigenpairs(
                np.diag(np.arange(1.0, 25.0)**2), np.eye(24), np.zeros((24, 24)), 3,
            )

    def test_campbell_rejects_unavailable_oscillating_modes(self, eigenvalue_model):
        with pytest.raises(ValueError, match="Requested 100 modes"):
            Solver(eigenvalue_model).solve(
                "campbell", speeds=[100.0], modes=100,
            )
