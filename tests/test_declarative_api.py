"""End-to-end tests for records, random sampling, and the unified solver."""

import numpy as np
import pytest

from spinniped import (
    BearingElement,
    BearingProperty,
    DiskElement,
    DiskProperty,
    Grid,
    Material,
    ModelBuilder,
    ModelDefinition,
    RandomDistribution,
    ShaftElement,
    ShaftProperty,
    Solver,
    build_model,
)


def _rotor_definition(*, random=False):
    """Return a small rotor definition with optional random references.

    Parameters
    ----------
    random : bool, optional
        Insert registered distribution references when true.

    Returns
    -------
    ModelDefinition
        Three-grid shaft, bearing, and disk model.
    """
    diameter = (1,) if random else 0.01
    density = (2,) if random else 7850.0
    bearing_stiffness = (3,) if random else 1.0e8
    disk_mass = (4,) if random else 2.0
    distributions = (
        [
            RandomDistribution(
                1, "shaft diameter", "normal", {"mean": 0.01, "stdv": 0.0001}
            ),
            RandomDistribution(
                2, "steel density", "uniform", {"low": 7800.0, "high": 7900.0}
            ),
            RandomDistribution(
                3,
                "bearing stiffness",
                "normal",
                {"mean": 1.0e8, "stdv": 1.0e6},
            ),
            RandomDistribution(
                4, "disk mass", "normal", {"mean": 2.0, "stdv": 0.02}
            ),
        ]
        if random
        else []
    )
    return ModelDefinition(
        grids=[Grid(1, z=0.0), Grid(2, z=0.5), Grid(3, z=1.0)],
        materials=[
            Material(
                1,
                density=density,
                young_modulus=2.0e11,
                poisson_ratio=0.3,
            )
        ],
        properties=[
            ShaftProperty(
                1,
                material=1,
                outer_diameter=diameter,
                theory="euler",
                rotary_inertia=False,
            ),
            BearingProperty(
                2,
                kxx=bearing_stiffness,
                kyy=bearing_stiffness,
                cxx=100.0,
                cyy=100.0,
            ),
            DiskProperty(
                3,
                mass=disk_mass,
                diametral_inertia=0.01,
                polar_inertia=0.02,
            ),
        ],
        elements=[
            ShaftElement(1, 1, 2, 1),
            ShaftElement(2, 2, 3, 1),
            BearingElement(3, 1, 2),
            BearingElement(4, 3, 2),
            DiskElement(5, 2, 3),
        ],
        distributions=distributions,
    )


def test_deterministic_build_resolves_every_distribution_to_its_mean():
    """Verify deterministic builds use distribution expected values."""
    resolved = ModelBuilder().build(
        _rotor_definition(random=True), stochastic=False, samples=20, seed=7
    )
    explicit_mean = build_model(_rotor_definition(random=False))

    assert resolved.samples == 1
    assert not resolved.stochastic
    assert resolved.seed == 7
    assert resolved.K.shape == (1, 18, 18)
    assert np.allclose(resolved.K, explicit_mean.K)
    assert np.allclose(resolved.M, explicit_mean.M)
    assert np.allclose(resolved.C, explicit_mean.C)
    assert np.allclose(resolved.G, explicit_mean.G)


def test_stochastic_build_produces_reproducible_model_realizations():
    """Verify seeded stochastic builds reproduce all realizations."""
    first = ModelBuilder().build(
        _rotor_definition(random=True),
        stochastic=True,
        samples=6,
        seed=123,
    )
    repeated = ModelBuilder().build(
        _rotor_definition(random=True),
        stochastic=True,
        samples=6,
        seed=123,
    )

    assert first.samples == 6
    assert first.stochastic
    assert first.seed == 123
    assert len(first.definitions) == 6
    assert first.K.shape == first.M.shape == first.C.shape == first.G.shape
    assert first.K.shape == (6, 18, 18)
    assert np.array_equal(first.K, repeated.K)
    assert np.array_equal(first.M, repeated.M)
    assert not np.array_equal(first.K[0], first.K[1])
    assert not np.array_equal(first.M[0], first.M[1])

    sampled_diameters = [
        definition.properties[0].outer_diameter
        for definition in first.definitions
    ]
    sampled_densities = [
        definition.materials[0].density for definition in first.definitions
    ]
    sampled_disk_masses = [
        definition.properties[2].mass for definition in first.definitions
    ]
    assert all(isinstance(value, float) for value in sampled_diameters)
    assert len(set(sampled_diameters)) > 1
    assert len(set(sampled_densities)) > 1
    assert len(set(sampled_disk_masses)) > 1


@pytest.mark.parametrize(
    ("distribution", "message"),
    [
        (RandomDistribution(1, "bad", "normal", {"mean": 1.0}), "requires only"),
        (
            RandomDistribution(
                1, "bad", "normal", {"mean": 1.0, "stdv": -0.1}
            ),
            "cannot be negative",
        ),
        (
            RandomDistribution(
                1, "bad", "uniform", {"low": 2.0, "high": 1.0}
            ),
            "high >= low",
        ),
        (
            RandomDistribution(
                1, "bad", "triangular", {"low": 0.0, "high": 1.0}
            ),
            "Unsupported distribution",
        ),
        (
            RandomDistribution(
                1, "bad", "normal", {"mean": np.nan, "stdv": 0.1}
            ),
            "finite values",
        ),
        (
            RandomDistribution(
                1,
                "bad",
                "normal",
                {"mean": 1.0, "stdv": 0.1, "units": "m"},
            ),
            "requires only",
        ),
    ],
)
def test_invalid_distribution_specifications_are_rejected(distribution, message):
    """Verify malformed distribution records are rejected."""
    definition = _rotor_definition()
    definition.properties[0] = ShaftProperty(
        1, material=1, outer_diameter=(1,)
    )
    definition.distributions.append(distribution)

    with pytest.raises(ValueError, match=message):
        ModelBuilder().build(definition)


def test_multivariate_components_preserve_requested_correlation():
    """Verify sampled components reproduce the requested correlation."""
    definition = _rotor_definition()
    definition.properties[1] = BearingProperty(2, kxx=(17, 0), kyy=(17, 1))
    definition.distributions.append(
        RandomDistribution(
            17,
            "correlated bearing stiffness",
            "multivariate_normal",
            {
                "mean": [1.0e8, 1.2e8],
                "stdv": [1.0e7, 2.0e7],
                "correlation": [[1.0, 0.8], [0.8, 1.0]],
            },
        )
    )

    model = ModelBuilder().build(
        definition, stochastic=True, samples=2_000, seed=42
    )
    kxx = [sample.properties[1].kxx for sample in model.definitions]
    kyy = [sample.properties[1].kyy for sample in model.definitions]

    assert np.corrcoef(kxx, kyy)[0, 1] == pytest.approx(0.8, abs=0.04)


def test_repeated_scalar_reference_reuses_the_same_sample():
    """Verify repeated references retrieve the same scalar draw."""
    model = ModelBuilder().build(
        _rotor_definition(random=True), stochastic=True, samples=10, seed=4
    )

    for definition in model.definitions:
        bearing = definition.properties[1]
        assert bearing.kxx == bearing.kyy


@pytest.mark.parametrize(
    ("reference", "message"),
    [
        ((99,), "unknown distribution ID"),
        ((1, 1), "out of range"),
    ],
)
def test_invalid_distribution_references_are_rejected(reference, message):
    """Verify missing distributions and invalid components are rejected."""
    definition = _rotor_definition(random=True)
    definition.properties[0] = ShaftProperty(1, 1, reference)

    with pytest.raises(ValueError, match=message):
        ModelBuilder().build(definition)


@pytest.mark.parametrize("samples", [0, -1, 1.5, "5", True])
def test_stochastic_sample_count_must_be_a_positive_integer(samples):
    """Verify stochastic sample counts are positive integers."""
    with pytest.raises(ValueError, match="positive integer"):
        ModelBuilder().build(
            _rotor_definition(random=True),
            stochastic=True,
            samples=samples,
        )


def test_builder_requires_a_model_definition():
    """Verify the builder rejects non-model root objects."""
    with pytest.raises(TypeError, match="ModelDefinition"):
        ModelBuilder().build({"grids": []})


def test_modal_solver_returns_real_modes_for_each_sample():
    """Verify modal analysis returns real modes for every realization."""
    model = ModelBuilder().build(
        _rotor_definition(random=True),
        stochastic=True,
        samples=3,
        seed=12,
    )
    result = Solver(model).solve("modal", fixed_dofs=[2, 5], modes=4)

    assert result["analysis"] == "modal"
    assert result["samples"] == 3
    assert result["stochastic"]
    assert result["eigenvalues"].shape == (3, 4)
    assert result["eigenvectors"].shape == (3, 16, 4)
    assert result["frequencies"].shape == (3, 4)
    assert np.isrealobj(result["eigenvalues"])
    assert np.all(result["eigenvalues"] > 0.0)
    assert np.allclose(
        result["frequencies"],
        np.sqrt(result["eigenvalues"]) / (2.0 * np.pi),
    )
    assert not np.array_equal(result["frequencies"][0], result["frequencies"][1])


def test_campbell_solver_returns_modes_at_every_speed():
    """Verify Campbell analysis returns branches at every speed."""
    model = ModelBuilder().build(_rotor_definition())
    result = Solver().solve(
        "campbell",
        model=model,
        fixed_dofs=[2, 5],
        speeds=[0.0, 100.0, 200.0],
        modes=4,
    )

    assert np.array_equal(result["speeds"], [0.0, 100.0, 200.0])
    assert np.allclose(result["speeds_hz"], result["speeds"] / (2.0 * np.pi))
    assert result["speeds_hz"].shape == (3,)
    assert result["eigenvalues"].shape == (1, 3, 4)
    assert result["frequencies"].shape == (1, 3, 4)
    assert result["eigenvectors"].shape == (1, 3, 32, 4)
    assert result["track_modes"] is True
    assert np.isfinite(result["frequencies"]).all()
    assert result["harmonics"].shape == (0,)
    assert result["critical_speeds"].shape == (1, 0, 4, 0)


def test_campbell_calculates_sampled_speed_interpolated_critical_speeds():
    """Verify tracked branches intersect requested harmonic indices."""
    natural_frequency = 10.0
    stiffness = (2.0 * np.pi * natural_frequency) ** 2
    definition = ModelDefinition(
        grids=[Grid(1)],
        properties=[
            BearingProperty(1, kxx=stiffness, kyy=stiffness),
            DiskProperty(2, mass=1.0, diametral_inertia=1.0),
        ],
        elements=[BearingElement(1, 1, 1), DiskElement(2, 1, 2)],
    )
    model = ModelBuilder().build(definition)
    result = Solver(model).solve(
        "campbell",
        fixed_dofs=[2, 3, 4, 5],
        speeds=2.0 * np.pi * np.array([0.0, 12.0, 20.0]),
        modes=2,
        harmonics=[1.0, 2.0],
    )

    assert np.array_equal(result["harmonics"], [1.0, 2.0])
    assert result["critical_speeds"].shape == (1, 2, 2, 1)
    assert np.all(result["critical_speed_counts"] == 1)
    assert np.allclose(
        result["critical_speeds_hz"][0, 0, :, 0],
        natural_frequency,
    )
    assert np.allclose(
        result["critical_speeds_hz"][0, 1, :, 0],
        natural_frequency / 2.0,
    )
    assert not np.any(
        np.isclose(
            result["critical_speeds_hz"][0, 0, 0, 0],
            [0.0, 12.0, 20.0],
        )
    )


def test_campbell_preserves_stochastic_critical_speed_correspondence():
    """Verify sample roots remain associated with tracked physical modes."""
    definition = ModelDefinition(
        grids=[Grid(1)],
        properties=[
            BearingProperty(1, kxx=(1,), kyy=(2,)),
            DiskProperty(2, mass=1.0, diametral_inertia=1.0),
        ],
        elements=[BearingElement(1, 1, 1), DiskElement(2, 1, 2)],
        distributions=[
            RandomDistribution(
                1,
                "x stiffness",
                "uniform",
                {
                    "low": (2.0 * np.pi * 8.0) ** 2,
                    "high": (2.0 * np.pi * 10.0) ** 2,
                },
            ),
            RandomDistribution(
                2,
                "y stiffness",
                "uniform",
                {
                    "low": (2.0 * np.pi * 12.0) ** 2,
                    "high": (2.0 * np.pi * 14.0) ** 2,
                },
            ),
        ],
    )
    model = ModelBuilder().build(
        definition, stochastic=True, samples=4, seed=9
    )
    result = Solver(model).solve(
        "campbell",
        fixed_dofs=[2, 3, 4, 5],
        speeds=2.0 * np.pi * np.array([0.0, 15.0]),
        modes=2,
        harmonics=[1.0],
    )

    expected = np.array(
        [
            [
                np.sqrt(sample.properties[0].kxx) / (2.0 * np.pi),
                np.sqrt(sample.properties[0].kyy) / (2.0 * np.pi),
            ]
            for sample in model.definitions
        ]
    )
    assert np.allclose(result["critical_speeds_hz"][:, 0, :, 0], expected)


def test_critical_speed_detection_retains_crossing_order_and_missing_values():
    """Verify multiple roots are ordered and absent roots use NaN padding."""
    speeds = 2.0 * np.pi * np.array([0.0, 1.0, 2.0, 3.0])
    synchronous = speeds / (2.0 * np.pi)
    frequencies = np.empty((2, 4, 1))
    frequencies[0, :, 0] = synchronous + [0.5, -0.5, 0.5, -0.5]
    frequencies[1, :, 0] = synchronous + 1.0

    critical, counts = Solver._critical_speeds(
        speeds, frequencies, np.array([1.0])
    )

    assert critical.shape == (2, 1, 1, 3)
    assert np.array_equal(counts[:, 0, 0], [3, 0])
    assert np.allclose(
        critical[0, 0, 0] / (2.0 * np.pi), [0.5, 1.5, 2.5]
    )
    assert np.isnan(critical[1, 0, 0]).all()


@pytest.mark.parametrize(
    ("options", "exception", "message"),
    [
        ({"harmonics": [0.0]}, ValueError, "positive finite"),
        ({"harmonics": [1.0, 1.0]}, ValueError, "duplicate"),
        (
            {"harmonics": [1.0], "track_modes": False},
            ValueError,
            "track_modes=True",
        ),
        (
            {"harmonics": [1.0], "speeds": [1.0, 0.0]},
            ValueError,
            "strictly increasing",
        ),
    ],
)
def test_campbell_validates_critical_speed_options(options, exception, message):
    """Verify critical-speed calculation requires coherent grid options."""
    model = ModelBuilder().build(_rotor_definition())
    solve_options = {"speeds": [0.0, 1.0], "modes": 2, **options}

    with pytest.raises(exception, match=message):
        Solver(model).solve(
            "campbell", fixed_dofs=[2, 5], **solve_options
        )


@pytest.mark.parametrize("modes", [0, -1, 1.5, True, "4"])
def test_campbell_solver_rejects_invalid_mode_counts(modes):
    """Verify Campbell analysis validates requested mode counts."""
    model = ModelBuilder().build(_rotor_definition())

    with pytest.raises(ValueError, match="modes must be a positive integer"):
        Solver(model).solve(
            "campbell", fixed_dofs=[2, 5], speeds=[0.0], modes=modes
        )


def test_campbell_solver_rejects_more_modes_than_are_available():
    """Verify Campbell analysis rejects unavailable mode counts."""
    model = ModelBuilder().build(_rotor_definition())

    with pytest.raises(ValueError, match="Requested 100 modes"):
        Solver(model).solve(
            "campbell", fixed_dofs=[2, 5], speeds=[0.0], modes=100
        )


def test_frequency_response_accepts_a_full_model_force():
    """Verify frequency response reduces a full-model force vector."""
    model = ModelBuilder().build(_rotor_definition())
    force = np.zeros(model.ndof)
    force[6] = 1.0
    result = Solver(model).solve(
        "frequency_response",
        fixed_dofs=[2, 5],
        frequencies=[10.0, 20.0],
        force=force,
        speed=50.0,
    )

    assert np.array_equal(result["frequencies"], [10.0, 20.0])
    assert result["response"].shape == (1, 2, 16)
    assert np.iscomplexobj(result["response"])
    assert np.isfinite(result["response"]).all()


def test_time_response_returns_displacement_and_velocity():
    """Verify time response returns displacement and velocity histories."""
    model = ModelBuilder().build(_rotor_definition())
    force = np.zeros(model.ndof)
    force[6] = 1.0
    times = np.linspace(0.0, 0.002, 5)
    result = Solver(model).solve(
        "time_response",
        fixed_dofs=[2, 5],
        times=times,
        force=force,
    )

    assert np.array_equal(result["times"], times)
    assert result["displacement"].shape == (1, 5, 16)
    assert result["velocity"].shape == (1, 5, 16)
    assert np.isfinite(result["displacement"]).all()
    assert np.isfinite(result["velocity"]).all()
    assert np.count_nonzero(result["displacement"][:, 0]) == 0


def test_time_response_reduces_a_full_model_initial_state():
    """Verify time response reduces full-model initial conditions."""
    model = ModelBuilder().build(_rotor_definition())
    fixed = [2, 5]
    free = np.setdiff1d(np.arange(model.ndof), fixed)
    full_displacement = np.linspace(1.0e-7, 1.8e-6, model.ndof)
    full_velocity = np.linspace(-1.8e-4, -1.0e-5, model.ndof)
    full_initial_state = np.concatenate((full_displacement, full_velocity))

    result = Solver(model).solve(
        "time_response",
        fixed_dofs=fixed,
        times=[0.0, 1.0e-8],
        force=np.zeros(model.ndof),
        initial_state=full_initial_state,
    )

    assert np.array_equal(result["displacement"][0, 0], full_displacement[free])
    assert np.array_equal(result["velocity"][0, 0], full_velocity[free])


def test_solver_validates_entry_point_arguments():
    """Verify the unified solver validates its public arguments."""
    model = ModelBuilder().build(_rotor_definition())

    with pytest.raises(TypeError, match="ModelBuilder"):
        Solver().solve("modal")
    with pytest.raises(TypeError, match="ModelBuilder"):
        Solver().solve("modal", model=object())
    with pytest.raises(ValueError, match="Unknown analysis"):
        Solver(model).solve("static")
    with pytest.raises(ValueError, match="out-of-range"):
        Solver(model).solve("modal", fixed_dofs=[model.ndof])
    with pytest.raises(TypeError, match="only integers"):
        Solver(model).solve("modal", fixed_dofs=[True])
    with pytest.raises(TypeError, match="only integers"):
        Solver(model).solve("modal", fixed_dofs=[1.5])
    with pytest.raises(TypeError, match="unexpected keyword"):
        Solver(model).solve("modal", fixed_dofs=[2, 5], modse=3)


def test_modal_rejects_cross_coupled_nonsymmetric_stiffness():
    """Verify real modal analysis rejects nonsymmetric stiffness."""
    definition = _rotor_definition()
    definition.properties[1] = BearingProperty(
        2, kxx=1.0e8, kyy=1.0e8, kxy=2.0e6, kyx=1.0e6
    )
    model = ModelBuilder().build(definition)

    with pytest.raises(ValueError, match="requires symmetric stiffness"):
        Solver(model).solve("modal", fixed_dofs=[2, 5])


def test_modal_keeps_small_physical_modes_in_a_wide_spectrum():
    """Verify modal filtering retains small positive physical modes."""
    definition = ModelDefinition(
        grids=[Grid(1)],
        properties=[
            BearingProperty(1, kxx=1.0e-4, kyy=1.0e8),
            DiskProperty(
                2, mass=1.0, diametral_inertia=1.0, polar_inertia=1.0
            ),
        ],
        elements=[BearingElement(1, 1, 1), DiskElement(2, 1, 2)],
    )
    model = ModelBuilder().build(definition)
    result = Solver(model).solve("modal", fixed_dofs=[2, 3, 4, 5])

    assert result["eigenvalues"].shape == (1, 2)
    assert np.allclose(result["eigenvalues"][0], [1.0e-4, 1.0e8])


def test_response_solvers_validate_load_and_time_shapes():
    """Verify response solvers reject malformed load and time arrays."""
    model = ModelBuilder().build(_rotor_definition())
    solver = Solver(model)

    with pytest.raises(ValueError, match="force must"):
        solver.solve(
            "frequency_response",
            fixed_dofs=[2, 5],
            frequencies=[10.0],
            force=[1.0, 2.0],
        )
    with pytest.raises(ValueError, match="times must"):
        solver.solve(
            "time_response",
            fixed_dofs=[2, 5],
            times=[0.0],
            force=np.zeros(model.ndof),
        )
