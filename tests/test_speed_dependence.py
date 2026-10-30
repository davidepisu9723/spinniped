"""Tests for tabulated rotor-speed-dependent coefficients and loads."""

from dataclasses import FrozenInstanceError

import numpy as np
import pytest

from spinniped import (
    BearingElement,
    BearingProperty,
    CoordinateSystem,
    DiskElement,
    DiskProperty,
    Grid,
    HarmonicLoad,
    Material,
    ModelBuilder,
    ModelDefinition,
    RandomDistribution,
    RotatingLoad,
    ShaftProperty,
    Solver,
    SpeedDependentValue,
    Unbalance,
)
from spinniped.builder import _interpolate_speed_dependent_value


FIXED_TO_X = [1, 2, 3, 4, 5]
FIXED_TO_XY = [2, 3, 4, 5]


def _one_dof_definition(
    *,
    stiffness,
    damping=0.0,
    mass=2.0,
    loads=(),
    distributions=(),
):
    """Return a one-grid model whose only retained DOF can be global x."""
    return ModelDefinition(
        grids=[Grid(1)],
        properties=[
            BearingProperty(1, kxx=stiffness, cxx=damping),
            DiskProperty(
                2,
                mass=mass,
                diametral_inertia=1.0,
                polar_inertia=1.0,
            ),
        ],
        elements=[BearingElement(1, 1, 1), DiskElement(2, 1, 2)],
        distributions=list(distributions),
        loads=list(loads),
    )


def _two_dof_unbalance_model(stiffness, damping):
    """Return an isotropic x-y model for signed synchronous response."""
    definition = ModelDefinition(
        grids=[Grid(1)],
        properties=[
            BearingProperty(
                1,
                kxx=stiffness,
                kyy=stiffness,
                cxx=damping,
                cyy=damping,
            ),
            DiskProperty(
                2,
                mass=2.0,
                diametral_inertia=1.0,
                polar_inertia=0.0,
            ),
        ],
        elements=[BearingElement(1, 1, 1), DiskElement(2, 1, 2)],
        loads=[Unbalance(1, 1, mass_eccentricity=2.0e-3, phase=0.2)],
    )
    return ModelBuilder().build(definition)


def test_speed_dependent_value_is_public_frozen_and_linearly_interpolated():
    """Verify public record behavior, signed knots, order, and repetition."""
    value = SpeedDependentValue(
        speeds=(-100, 0, 100),
        values=(10, 20, 50),
    )

    assert not hasattr(value, "__dict__")
    assert value.speeds == (-100, 0, 100)
    assert value.values == (10, 20, 50)
    with pytest.raises(FrozenInstanceError):
        value.values = (1.0, 2.0, 3.0)

    evaluated = _interpolate_speed_dependent_value(
        value,
        np.array([100.0, -50.0, 0.0, -50.0]),
        "test value",
    )
    assert np.array_equal(evaluated, [50.0, 15.0, 20.0, 15.0])


@pytest.mark.parametrize(
    ("speeds", "values", "message"),
    [
        ([-100, 0, 100], (10, 20, 50), "speeds must be a tuple"),
        ((-100, 0, 100), [10, 20, 50], "values must be a tuple"),
    ],
)
def test_speed_dependent_value_requires_tuple_inputs(
    speeds,
    values,
    message,
):
    """Verify speed-table sequences are explicit immutable tuples."""
    with pytest.raises(TypeError, match=message):
        SpeedDependentValue(speeds=speeds, values=values)


@pytest.mark.parametrize(
    ("speeds", "values", "exception", "message"),
    [
        (
            (0.0,),
            (1.0,),
            ValueError,
            "at least two speeds",
        ),
        (
            (0.0, 1.0),
            (1.0,),
            ValueError,
            "same length",
        ),
        (
            (0.0, 0.0),
            (1.0, 2.0),
            ValueError,
            "strictly increasing",
        ),
        (
            (1.0, 0.0),
            (1.0, 2.0),
            ValueError,
            "strictly increasing",
        ),
        (
            (0.0, np.inf),
            (1.0, 2.0),
            ValueError,
            "finite",
        ),
        (
            (0.0, 1.0),
            (1.0, np.nan),
            ValueError,
            "finite",
        ),
        (
            (False, 1.0),
            (1.0, 2.0),
            TypeError,
            "speeds must contain only real numeric values",
        ),
        (
            (0.0, 1.0),
            (1.0, True),
            TypeError,
            "real number or random reference",
        ),
        (
            (0.0, 1.0),
            ((1, 0, 2), 1.0),
            TypeError,
            "random reference",
        ),
        (
            (0.0, 1.0),
            ((0,), 1.0),
            ValueError,
            "positive distribution ID",
        ),
        (
            (0.0, 1.0),
            ((1, -1), 1.0),
            ValueError,
            "nonnegative component index",
        ),
    ],
)
def test_speed_dependent_value_rejects_malformed_inputs(
    speeds,
    values,
    exception,
    message,
):
    """Verify malformed speed tables fail during record construction."""
    with pytest.raises(exception, match=message):
        SpeedDependentValue(speeds=speeds, values=values)


@pytest.mark.parametrize("speed", [-0.01, 100.01])
def test_speed_dependent_value_rejects_extrapolation(speed):
    """Verify interpolation accepts endpoints but rejects either outer side."""
    value = SpeedDependentValue((0, 100), (400, 900))
    assert np.array_equal(
        _interpolate_speed_dependent_value(
            value, np.array([0.0, 100.0]), "bearing"
        ),
        [400.0, 900.0],
    )

    with pytest.raises(ValueError, match="outside its speed range"):
        _interpolate_speed_dependent_value(
            value,
            np.array([speed]),
            "bearing",
        )


def test_builder_keeps_base_matrices_three_dimensional_and_adds_deltas():
    """Verify tabulated coefficients stay out of permanent K and C arrays."""
    stiffness = SpeedDependentValue((0, 100), (400, 900))
    damping = SpeedDependentValue((0, 100), (4, 14))
    definition = _one_dof_definition(
        stiffness=stiffness,
        damping=damping,
    )
    definition.properties.insert(1, BearingProperty(3, kxx=25.0, cxx=2.0))
    definition.elements.insert(1, BearingElement(3, 1, 3))

    model = ModelBuilder().build(definition)

    assert model.K.shape == model.M.shape == model.C.shape == model.G.shape
    assert model.K.shape == (1, 6, 6)
    assert model.K[0, 0, 0] == 25.0
    assert model.C[0, 0, 0] == 2.0
    assert model.definitions[0].properties[0].kxx == stiffness

    result = Solver(model).solve(
        "modal",
        fixed_dofs=FIXED_TO_X,
        speed=50.0,
    )
    assert result["eigenvalues"][0, 0] == pytest.approx((25.0 + 650.0) / 2.0)
    assert model.K[0, 0, 0] == 25.0
    assert model.C[0, 0, 0] == 2.0


def test_rotated_two_grid_bearing_transforms_speed_dependent_deltas():
    """Verify transformed relative-bearing bases retain every block sign."""
    frame = CoordinateSystem(
        7,
        x_axis=(0.0, 1.0, 0.0),
        xy_plane=(-1.0, 0.0, 0.0),
    )
    definition = ModelDefinition(
        grids=[Grid(1), Grid(2)],
        coordinate_systems=[frame],
        properties=[
            BearingProperty(
                1,
                kxx=SpeedDependentValue((0, 100), (10.0, 30.0)),
                cxx=SpeedDependentValue((0, 100), (2.0, 6.0)),
            )
        ],
        elements=[
            BearingElement(
                1,
                grid=1,
                grid_b=2,
                property=1,
                coordinate_system=7,
            )
        ],
    )
    model = ModelBuilder().build(definition)

    stiffness_delta = Solver._speed_matrix_deltas(
        model,
        model._stiffness_speed_contributions,
        np.array([50.0]),
    )[0, 0]
    damping_delta = Solver._speed_matrix_deltas(
        model,
        model._damping_speed_contributions,
        np.array([50.0]),
    )[0, 0]
    y_dofs = [1, 7]

    assert np.allclose(
        stiffness_delta[np.ix_(y_dofs, y_dofs)],
        [[20.0, -20.0], [-20.0, 20.0]],
    )
    assert np.allclose(
        damping_delta[np.ix_(y_dofs, y_dofs)],
        [[4.0, -4.0], [-4.0, 4.0]],
    )
    assert np.count_nonzero(stiffness_delta) == 4
    assert np.count_nonzero(damping_delta) == 4


def test_deterministic_build_resolves_random_speed_ordinates_to_means():
    """Verify scalar references at knots resolve before interpolation."""
    value = SpeedDependentValue((0, 100), ((1,), (2,)))
    distributions = [
        RandomDistribution(1, "low-speed stiffness", "normal", {"mean": 400.0, "stdv": 20.0}),
        RandomDistribution(2, "high-speed stiffness", "normal", {"mean": 900.0, "stdv": 30.0}),
    ]
    model = ModelBuilder().build(
        _one_dof_definition(
            stiffness=value,
            distributions=distributions,
        )
    )

    resolved = model.definitions[0].properties[0].kxx
    assert resolved.speeds == (0, 100)
    assert resolved.values == (400.0, 900.0)

    result = Solver(model).solve(
        "modal",
        fixed_dofs=FIXED_TO_X,
        speed=50.0,
    )
    assert result["eigenvalues"][0, 0] == pytest.approx(650.0 / 2.0)


def test_multivariate_speed_ordinates_are_correlated_and_reproducible():
    """Verify a joint draw becomes one deterministic curve per realization."""
    value = SpeedDependentValue((0, 100), ((17, 0), (17, 1)))
    distribution = RandomDistribution(
        17,
        "correlated stiffness curve",
        "multivariate_normal",
        {
            "mean": [400.0, 900.0],
            "stdv": [40.0, 90.0],
            "correlation": [[1.0, 0.75], [0.75, 1.0]],
        },
    )
    definition = _one_dof_definition(
        stiffness=value,
        distributions=[distribution],
    )

    first = ModelBuilder().build(
        definition,
        stochastic=True,
        samples=1_000,
        seed=42,
    )
    repeated = ModelBuilder().build(
        definition,
        stochastic=True,
        samples=1_000,
        seed=42,
    )
    first_values = np.asarray(
        [sample.properties[0].kxx.values for sample in first.definitions]
    )
    repeated_values = np.asarray(
        [sample.properties[0].kxx.values for sample in repeated.definitions]
    )

    assert np.array_equal(first_values, repeated_values)
    assert np.corrcoef(first_values.T)[0, 1] == pytest.approx(0.75, abs=0.05)
    assert first.K.shape == (1_000, 6, 6)


def test_repeated_random_reference_reuses_one_draw_at_every_speed_knot():
    """Verify a repeated scalar reference is not resampled across speed."""
    value = SpeedDependentValue((0, 100), ((8,), (8,)))
    definition = _one_dof_definition(
        stiffness=value,
        distributions=[
            RandomDistribution(
                8,
                "common stiffness",
                "uniform",
                {"low": 400.0, "high": 900.0},
            )
        ],
    )
    model = ModelBuilder().build(
        definition,
        stochastic=True,
        samples=8,
        seed=9,
    )

    for sample in model.definitions:
        assert sample.properties[0].kxx.values[0] == sample.properties[0].kxx.values[1]

    result = Solver(model).solve(
        "campbell",
        fixed_dofs=FIXED_TO_X,
        speeds=[0.0, 100.0],
        modes=1,
    )
    assert np.allclose(result["frequencies"][:, 0], result["frequencies"][:, 1])


def test_campbell_uses_interpolated_stiffness_and_damping_at_every_speed():
    """Verify K(Omega) and C(Omega) enter each state-space eigensolution."""
    stiffness = SpeedDependentValue((0, 100), (800.0, 1800.0))
    damping = SpeedDependentValue((0, 100), (4.0, 8.0))
    model = ModelBuilder().build(
        _one_dof_definition(stiffness=stiffness, damping=damping)
    )
    speeds = np.array([0.0, 50.0, 100.0])

    result = Solver(model).solve(
        "campbell",
        fixed_dofs=FIXED_TO_X,
        speeds=speeds,
        modes=1,
    )

    interpolated_k = np.array([800.0, 1300.0, 1800.0])
    interpolated_c = np.array([4.0, 6.0, 8.0])
    expected = -interpolated_c / 4.0 + 1j * np.sqrt(
        interpolated_k / 2.0 - (interpolated_c / 4.0) ** 2
    )
    assert np.allclose(result["eigenvalues"][0, :, 0], expected)


def test_frequency_response_uses_rotor_speed_not_forcing_frequency():
    """Verify one operating point is reused across the excitation sweep."""
    stiffness = SpeedDependentValue((0, 100), (1_000.0, 3_000.0))
    damping = SpeedDependentValue((0, 100), (10.0, 30.0))
    model = ModelBuilder().build(
        _one_dof_definition(stiffness=stiffness, damping=damping)
    )
    frequencies = np.array([1.0, 2.0])
    force = np.zeros(6)
    force[0] = 1.0

    result = Solver(model).solve(
        "frequency_response",
        fixed_dofs=FIXED_TO_X,
        frequencies=frequencies,
        force=force,
        speed=50.0,
    )

    omega = 2.0 * np.pi * frequencies
    expected = 1.0 / (2_000.0 - 2.0 * omega**2 + 1j * omega * 20.0)
    assert np.allclose(result["response"][0, :, 0], expected)


def test_time_response_uses_one_interpolated_operating_point():
    """Verify scalar-speed integration matches an equivalent constant model."""
    speed_model = ModelBuilder().build(
        _one_dof_definition(
            stiffness=SpeedDependentValue((0, 100), (4.0, 16.0)),
            damping=SpeedDependentValue((0, 100), (0.2, 0.6)),
            mass=1.0,
        )
    )
    constant_model = ModelBuilder().build(
        _one_dof_definition(stiffness=10.0, damping=0.4, mass=1.0)
    )
    options = {
        "fixed_dofs": FIXED_TO_X,
        "times": np.linspace(0.0, 0.5, 8),
        "force": np.zeros(6),
        "initial_state": [1.0, 0.0],
        "speed": 50.0,
    }

    speed_result = Solver(speed_model).solve("time_response", **options)
    constant_result = Solver(constant_model).solve("time_response", **options)

    assert np.allclose(
        speed_result["displacement"],
        constant_result["displacement"],
    )
    assert np.allclose(speed_result["velocity"], constant_result["velocity"])


def test_unbalance_response_uses_literal_signed_speed_for_bearing_tables():
    """Verify negative speed is interpolated at negative, not absolute, knots."""
    stiffness = SpeedDependentValue(
        (-100, 0, 100),
        (1_000.0, 2_000.0, 5_000.0),
    )
    damping = SpeedDependentValue((-100, 0, 100), (4.0, 8.0, 20.0))
    model = _two_dof_unbalance_model(stiffness, damping)
    speeds = np.array([-50.0, 0.0, 100.0])

    result = Solver(model).solve(
        "unbalance_response",
        fixed_dofs=FIXED_TO_XY,
        speeds=speeds,
        load_ids=[1],
    )

    expected = []
    for speed, k_value, c_value in zip(
        speeds,
        [1_500.0, 2_000.0, 5_000.0],
        [6.0, 8.0, 20.0],
        strict=True,
    ):
        if speed == 0.0:
            expected.append(np.zeros(2, dtype=complex))
            continue
        sense = np.sign(speed)
        force_x = 2.0e-3 * speed**2 * np.exp(1j * sense * 0.2)
        force = np.array([force_x, -1j * sense * force_x])
        dynamic_stiffness = (
            k_value - 2.0 * abs(speed) ** 2 + 1j * abs(speed) * c_value
        )
        expected.append(force / dynamic_stiffness)

    assert np.allclose(result["response"][0], expected)


def test_campbell_preflights_complete_speed_grid_before_first_solve(monkeypatch):
    """Verify a late out-of-range point prevents every numerical solve."""
    model = ModelBuilder().build(
        _one_dof_definition(
            stiffness=SpeedDependentValue((0, 100), (400.0, 900.0))
        )
    )
    solve_called = False

    def unexpected_solve(*args, **kwargs):
        nonlocal solve_called
        solve_called = True
        raise AssertionError("linear solve started before speed preflight")

    monkeypatch.setattr(np.linalg, "solve", unexpected_solve)

    with pytest.raises(ValueError, match="outside its speed range"):
        Solver(model).solve(
            "campbell",
            fixed_dofs=FIXED_TO_X,
            speeds=[50.0, 101.0],
            modes=1,
        )
    assert not solve_called


def test_harmonic_load_components_interpolate_transform_and_superimpose():
    """Verify F(Omega) for force and moment components in a rotated frame."""
    phase = 0.25
    frame = CoordinateSystem(
        7,
        x_axis=(0.0, 1.0, 0.0),
        xy_plane=(-1.0, 0.0, 0.0),
    )
    varying = HarmonicLoad(
        1,
        grid=1,
        force=(SpeedDependentValue((0, 100), (2.0, 6.0)), 1.0, 0.0),
        moment=(0.0, 0.0, SpeedDependentValue((0, 100), (4.0, 8.0))),
        phase=phase,
        coordinate_system=7,
    )
    constant = HarmonicLoad(
        2,
        grid=1,
        force=(3.0, 0.0, 0.0),
        phase=phase,
        coordinate_system=7,
    )
    definition = _one_dof_definition(
        stiffness=1_000.0,
        loads=[varying, constant],
    )
    definition.coordinate_systems.append(frame)
    model = ModelBuilder().build(definition)

    assembled = Solver._assemble_harmonic_loads(
        model,
        0,
        (1, 2),
        speed=50.0,
    )
    expected_local = np.array([7.0, 1.0, 0.0, 0.0, 0.0, 6.0])
    expected_global = np.array([-1.0, 7.0, 0.0, 0.0, 0.0, 6.0])
    assert np.allclose(assembled, expected_global * np.exp(1j * phase))
    assert not np.array_equal(expected_local, expected_global)


def test_frequency_response_interpolates_rotating_load_at_rotor_speed():
    """Verify a tabulated circular-force magnitude drives the full response."""
    load = RotatingLoad(
        4,
        grid=1,
        force_magnitude=SpeedDependentValue((0, 100), (2.0, 6.0)),
        phase=0.3,
        sense="negative",
    )
    definition = ModelDefinition(
        grids=[Grid(1)],
        properties=[
            BearingProperty(1, kxx=2_000.0, kyy=2_000.0),
            DiskProperty(2, mass=2.0, diametral_inertia=1.0),
        ],
        elements=[BearingElement(1, 1, 1), DiskElement(2, 1, 2)],
        loads=[load],
    )
    model = ModelBuilder().build(definition)
    frequency = 2.0

    result = Solver(model).solve(
        "frequency_response",
        fixed_dofs=FIXED_TO_XY,
        frequencies=[frequency],
        speed=50.0,
        load_ids=[4],
    )

    radial = 4.0 * np.exp(-1j * 0.3)
    force = np.array([radial, 1j * radial])
    omega = 2.0 * np.pi * frequency
    expected = force / (2_000.0 - 2.0 * omega**2)
    assert np.allclose(result["response"][0, 0], expected)


def test_selected_load_is_preflighted_but_unselected_tables_are_ignored(
    monkeypatch,
):
    """Verify load range checks occur before solve and are selection scoped."""
    out_of_range = HarmonicLoad(
        1,
        grid=1,
        force=(SpeedDependentValue((0, 10), (1.0, 2.0)), 0.0, 0.0),
    )
    static = HarmonicLoad(2, grid=1, force=(1.0, 0.0, 0.0))
    model = ModelBuilder().build(
        _one_dof_definition(
            stiffness=1_000.0,
            loads=[out_of_range, static],
        )
    )

    # An unselected table must not constrain an otherwise valid analysis.
    valid = Solver(model).solve(
        "frequency_response",
        fixed_dofs=FIXED_TO_X,
        frequencies=[1.0],
        speed=50.0,
        load_ids=[2],
    )
    assert np.isfinite(valid["response"]).all()

    solve_called = False

    def unexpected_solve(*args, **kwargs):
        nonlocal solve_called
        solve_called = True
        raise AssertionError("linear solve started before load preflight")

    monkeypatch.setattr(np.linalg, "solve", unexpected_solve)
    with pytest.raises(ValueError, match="outside its speed range"):
        Solver(model).solve(
            "frequency_response",
            fixed_dofs=FIXED_TO_X,
            frequencies=[1.0],
            speed=50.0,
            load_ids=[1],
        )
    assert not solve_called


def test_rotating_load_speed_table_must_remain_nonnegative():
    """Verify every source ordinate satisfies the magnitude constraint."""
    with pytest.raises(ValueError, match="values must be nonnegative"):
        RotatingLoad(
            1,
            grid=1,
            force_magnitude=SpeedDependentValue((0, 100), (1.0, -1.0)),
        )


def test_builder_rejects_mutated_speed_dependent_distribution_parameters():
    """Verify mutable distribution dictionaries remain guarded at build time."""
    distribution = RandomDistribution(
        1,
        "mutable",
        "normal",
        {"mean": 0.0, "stdv": 1.0},
    )
    distribution.parameters["mean"] = SpeedDependentValue((0, 1), (1.0, 2.0))
    definition = ModelDefinition(grids=[Grid(1)], distributions=[distribution])

    with pytest.raises(TypeError, match="does not accept SpeedDependentValue"):
        ModelBuilder().build(definition)


def test_model_definition_rejects_speed_dependent_spin_axis_at_declaration():
    """Verify the immutable spin axis cannot contain a rotor-speed table."""
    table = SpeedDependentValue((0, 1), (1.0, 2.0))

    with pytest.raises(TypeError, match="real number or random reference"):
        ModelDefinition(spin_axis=(0.0, 0.0, table))


@pytest.mark.parametrize("field", ["mass_eccentricity", "phase"])
def test_unbalance_rejects_speed_dependence_at_declaration(field):
    """Verify unbalance parameters cannot contain rotor-speed tables."""
    arguments = {"mass_eccentricity": 1.0e-4, "phase": 0.0}
    arguments[field] = SpeedDependentValue((0, 1), (1.0e-4, 2.0e-4))

    with pytest.raises(TypeError, match="real number or random reference"):
        Unbalance(1, grid=1, **arguments)


def test_harmonic_load_rejects_speed_dependent_phase_at_declaration():
    """Verify harmonic phase cannot contain a rotor-speed table."""
    table = SpeedDependentValue((0, 1), (0.0, 1.0))

    with pytest.raises(TypeError, match="real number or random reference"):
        HarmonicLoad(1, grid=1, phase=table)


def test_grid_rejects_speed_dependent_coordinates_at_declaration():
    """Verify grid coordinates accept scalar parameters, not speed tables."""
    table = SpeedDependentValue((0, 1), (0, 1))

    with pytest.raises(TypeError, match="real number or random reference"):
        Grid(1, x=table)


def test_material_rejects_speed_dependent_parameters_at_declaration():
    """Verify material parameters cannot contain speed tables."""
    table = SpeedDependentValue((0, 1), (1, 2))

    with pytest.raises(TypeError, match="real number or random reference"):
        Material(1, density=table, young_modulus=1.0, poisson_ratio=0.3)


def test_shaft_property_rejects_speed_dependence_at_declaration():
    """Verify shaft-property parameters cannot contain speed tables."""
    table = SpeedDependentValue((0, 1), (0.1, 0.2))

    with pytest.raises(TypeError, match="real number or random reference"):
        ShaftProperty(1, material=1, outer_diameter=table)


def test_disk_property_rejects_speed_dependence_at_declaration():
    """Verify disk-property parameters cannot contain speed tables."""
    table = SpeedDependentValue((0, 1), (1.0, 2.0))

    with pytest.raises(TypeError, match="real number or random reference"):
        DiskProperty(1, mass=table)
