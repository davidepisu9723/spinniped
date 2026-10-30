"""Tests for declarative harmonic, rotating, and unbalance excitation."""

import numpy as np
import pytest

from spinniped import (
    BearingElement,
    BearingProperty,
    BuiltModel,
    CoordinateSystem,
    DiskElement,
    DiskProperty,
    Grid,
    HarmonicLoad,
    ModelBuilder,
    ModelDefinition,
    RandomDistribution,
    RotatingLoad,
    Solver,
    SynchronousLoad,
    Unbalance,
)


def _one_grid_model(
    loads,
    *,
    coordinate_systems=None,
    distributions=None,
    stochastic=False,
    samples=1,
    seed=None,
):
    """Build one fully inertial grid with anisotropic grounded bearings."""
    definition = ModelDefinition(
        grids=[Grid(10)],
        coordinate_systems=list(coordinate_systems or []),
        properties=[
            BearingProperty(
                1,
                kxx=4_000.0,
                kyy=6_000.0,
                kzz=8_000.0,
                cxx=4.0,
                cyy=6.0,
                czz=8.0,
            ),
            DiskProperty(
                2,
                mass=2.0,
                diametral_inertia=0.5,
                polar_inertia=0.8,
            ),
        ],
        elements=[
            BearingElement(1, grid=10, property=1),
            DiskElement(2, grid=10, property=2),
        ],
        distributions=list(distributions or []),
        loads=list(loads),
    )
    return ModelBuilder().build(
        definition,
        stochastic=stochastic,
        samples=samples,
        seed=seed,
    )


def test_harmonic_load_matches_an_equivalent_global_force_vector():
    """Verify force, moment, phase, and frame rotation in a full response."""
    phase = 0.3
    frame = CoordinateSystem(
        7,
        x_axis=(0.0, 1.0, 0.0),
        xy_plane=(-1.0, 0.0, 0.0),
    )
    load = HarmonicLoad(
        17,
        grid=10,
        force=(1.0, 2.0, 3.0),
        moment=(4.0, 5.0, 6.0),
        phase=phase,
        coordinate_system=7,
    )
    model = _one_grid_model([load], coordinate_systems=[frame])
    expected_force = (
        np.array([-2.0, 1.0, 3.0, -5.0, 4.0, 6.0])
        * np.exp(1j * phase)
    )

    assembled_force = Solver._assemble_harmonic_loads(model, 0, (17,))
    declarative = Solver(model).solve(
        "frequency_response",
        frequencies=[7.0],
        speed=12.0,
        load_ids=[17],
    )
    explicit = Solver(model).solve(
        "frequency_response",
        frequencies=[7.0],
        speed=12.0,
        force=expected_force,
    )

    assert np.allclose(assembled_force, expected_force)
    assert np.allclose(declarative["response"], explicit["response"])
    assert declarative["load_ids"] == (17,)
    assert declarative["speed"] == 12.0


def test_synchronous_load_derives_frequency_from_order_and_speed():
    """Verify frequency is order-locked and the harmonic vector is preserved."""
    speed = -20.0
    harmonic_order = 2.5
    phase = 0.3
    force_amplitude = 10.0
    model = _one_grid_model(
        [
            SynchronousLoad(
                21,
                grid=10,
                harmonic_order=harmonic_order,
                force=(force_amplitude, 0.0, 0.0),
                phase=phase,
            )
        ]
    )

    result = Solver(model).solve(
        "frequency_response",
        fixed_dofs=[1, 2, 3, 4, 5],
        speed=speed,
        load_ids=[21],
    )

    angular_frequency = harmonic_order * abs(speed)
    dynamic_stiffness = (
        4_000.0
        - 2.0 * angular_frequency**2
        + 1j * angular_frequency * 4.0
    )
    expected = (
        force_amplitude * np.exp(1j * phase) / dynamic_stiffness
    )

    assert result["frequencies"] == pytest.approx(
        [angular_frequency / (2.0 * np.pi)]
    )
    assert result["harmonic_order"] == harmonic_order
    assert result["speed"] == speed
    assert result["load_ids"] == (21,)
    assert result["response"].shape == (1, 1, 1)
    assert result["response"][0, 0, 0] == pytest.approx(expected)


def test_synchronous_loads_of_one_order_superimpose():
    """Verify same-order records share one derived-frequency solve."""
    loads = [
        SynchronousLoad(
            1,
            grid=10,
            harmonic_order=2.0,
            force=(3.0, 0.0, 0.0),
        ),
        SynchronousLoad(
            2,
            grid=10,
            harmonic_order=2.0,
            force=(4.0, 0.0, 0.0),
            phase=0.2,
        ),
    ]
    solver = Solver(_one_grid_model(loads))
    options = {
        "fixed_dofs": [1, 2, 3, 4, 5],
        "speed": 30.0,
    }

    first = solver.solve("frequency_response", load_ids=[1], **options)
    second = solver.solve("frequency_response", load_ids=[2], **options)
    combined = solver.solve("frequency_response", load_ids=[2, 1], **options)

    assert np.allclose(
        combined["response"],
        first["response"] + second["response"],
    )
    assert combined["load_ids"] == (2, 1)


def test_synchronous_frequency_response_rejects_ambiguous_sources():
    """Verify order-derived response cannot accept contradictory frequencies."""
    model = _one_grid_model(
        [
            HarmonicLoad(1, grid=10, force=(1.0, 0.0, 0.0)),
            SynchronousLoad(
                2,
                grid=10,
                harmonic_order=1.0,
                force=(1.0, 0.0, 0.0),
            ),
            SynchronousLoad(
                3,
                grid=10,
                harmonic_order=2.0,
                force=(1.0, 0.0, 0.0),
            ),
        ]
    )
    solver = Solver(model)
    options = {"fixed_dofs": [1, 2, 3, 4, 5], "speed": 20.0}

    with pytest.raises(ValueError, match="frequencies must be omitted"):
        solver.solve(
            "frequency_response",
            frequencies=[1.0],
            load_ids=[2],
            **options,
        )
    with pytest.raises(ValueError, match="cannot mix"):
        solver.solve("frequency_response", load_ids=[1, 2], **options)
    with pytest.raises(ValueError, match="share one harmonic_order"):
        solver.solve("frequency_response", load_ids=[2, 3], **options)
    with pytest.raises(ValueError, match="requires frequencies"):
        solver.solve("frequency_response", load_ids=[1], **options)


def test_rotating_load_sense_phase_and_superposition_are_preserved():
    """Verify exact circular phasors and additive selected-load response."""
    magnitude = 9.0
    phase = 0.4
    positive = RotatingLoad(
        11,
        grid=10,
        force_magnitude=magnitude,
        phase=phase,
        sense="positive",
    )
    negative = RotatingLoad(
        42,
        grid=10,
        force_magnitude=magnitude,
        phase=phase,
        sense="negative",
    )
    model = _one_grid_model([positive, negative])

    positive_force = Solver._assemble_harmonic_loads(model, 0, (11,))
    negative_force = Solver._assemble_harmonic_loads(model, 0, (42,))
    expected_positive_x = magnitude * np.exp(1j * phase)
    expected_negative_x = magnitude * np.exp(-1j * phase)
    assert np.allclose(
        positive_force[:2],
        [expected_positive_x, -1j * expected_positive_x],
    )
    assert np.allclose(
        negative_force[:2],
        [expected_negative_x, 1j * expected_negative_x],
    )

    options = {
        "fixed_dofs": [2, 3, 4, 5],
        "frequencies": [6.0],
        "speed": 20.0,
    }
    positive_result = Solver(model).solve(
        "frequency_response", load_ids=[11], **options
    )
    negative_result = Solver(model).solve(
        "frequency_response", load_ids=[42], **options
    )
    combined_result = Solver(model).solve(
        "frequency_response", load_ids=[42, 11], **options
    )

    assert np.allclose(
        combined_result["response"],
        positive_result["response"] + negative_result["response"],
    )
    assert combined_result["load_ids"] == (42, 11)


def test_frequency_response_uses_each_samples_resolved_load():
    """Verify stochastic amplitudes and phases are not reused from sample zero."""
    distributions = [
        RandomDistribution(
            1,
            "force amplitude",
            "uniform",
            {"low": 8.0, "high": 12.0},
        ),
        RandomDistribution(
            2,
            "force phase",
            "uniform",
            {"low": 0.1, "high": 0.5},
        ),
    ]
    model = _one_grid_model(
        [HarmonicLoad(1, grid=10, force=((1,), 0.0, 0.0), phase=(2,))],
        distributions=distributions,
        stochastic=True,
        samples=4,
        seed=8,
    )
    frequency = 5.0
    result = Solver(model).solve(
        "frequency_response",
        fixed_dofs=[1, 2, 3, 4, 5],
        frequencies=[frequency],
        load_ids=[1],
    )

    omega = 2.0 * np.pi * frequency
    dynamic_stiffness = 4_000.0 - 2.0 * omega**2 + 1j * omega * 4.0
    expected = np.array(
        [
            sample.loads[0].force[0] * np.exp(1j * sample.loads[0].phase)
            / dynamic_stiffness
            for sample in model.definitions
        ]
    )

    assert result["response"].shape == (4, 1, 1)
    assert np.allclose(result["response"][:, 0, 0], expected)
    assert len(np.unique(expected)) == 4


def test_unbalance_response_matches_signed_synchronous_dynamic_stiffness():
    """Verify speed-squared forcing, reverse rotation, and gyroscopic scaling."""
    stiffness = np.diag([1_600.0, 2_100.0, 1.0, 1.0, 1.0, 1.0])
    mass = np.diag([2.0, 3.0, 1.0, 1.0, 1.0, 1.0])
    damping = np.diag([5.0, 7.0, 0.0, 0.0, 0.0, 0.0])
    gyroscopic = np.zeros((6, 6))
    gyroscopic[0, 1] = 0.2
    gyroscopic[1, 0] = -0.2
    mass_eccentricity = 2.5e-3
    phase = 0.35
    definition = ModelDefinition(
        grids=[Grid(10)],
        loads=[Unbalance(8, 10, mass_eccentricity, phase=phase)],
    )
    model = BuiltModel(
        stiffness=stiffness[None, :, :],
        mass=mass[None, :, :],
        damping=damping[None, :, :],
        gyroscopic=gyroscopic[None, :, :],
        grid_ids=(10,),
        coordinates=np.zeros((1, 1, 3)),
        definitions=(definition,),
        stochastic=False,
        seed=None,
    )
    speeds = np.array([-18.0, 0.0, 25.0])

    result = Solver(model).solve(
        "unbalance_response",
        fixed_dofs=[2, 3, 4, 5],
        speeds=speeds,
        load_ids=[8],
    )

    expected = []
    for speed in speeds:
        omega = abs(speed)
        dynamic_stiffness = (
            stiffness[:2, :2]
            - omega**2 * mass[:2, :2]
            + 1j
            * omega
            * (damping[:2, :2] + speed * gyroscopic[:2, :2])
        )
        if speed == 0.0:
            force = np.zeros(2, dtype=complex)
        else:
            sense = np.sign(speed)
            force_x = (
                mass_eccentricity
                * speed**2
                * np.exp(1j * sense * phase)
            )
            force = np.array([force_x, -1j * sense * force_x])
        expected.append(np.linalg.solve(dynamic_stiffness, force))

    assert result["response"].shape == (1, 3, 2)
    assert np.allclose(result["response"][0], expected)
    assert np.allclose(result["speeds_hz"], speeds / (2.0 * np.pi))
    assert np.allclose(result["frequencies"], np.abs(speeds) / (2.0 * np.pi))
    assert result["load_ids"] == (8,)


def test_unbalance_uses_frame_alignment_to_follow_the_rotor():
    """Verify antiparallel load axes produce the same global rotating force."""
    opposite_frame = CoordinateSystem(
        9,
        x_axis=(1.0, 0.0, 0.0),
        xy_plane=(0.0, -1.0, 0.0),
    )
    model = _one_grid_model(
        [
            Unbalance(1, grid=10, mass_eccentricity=1.0e-3),
            Unbalance(
                2,
                grid=10,
                mass_eccentricity=1.0e-3,
                coordinate_system=9,
            ),
        ],
        coordinate_systems=[opposite_frame],
    )

    parallel = Solver._assemble_unbalances(model, 0, (1,), 30.0)
    antiparallel = Solver._assemble_unbalances(model, 0, (2,), 30.0)

    assert np.allclose(parallel, antiparallel)
    assert np.allclose(parallel[:2], [0.9, -0.9j])


def test_response_load_selection_is_validated():
    """Verify source exclusivity, ID validity, and analysis compatibility."""
    model = _one_grid_model(
        [
            HarmonicLoad(1, grid=10, force=(1.0, 0.0, 0.0)),
            RotatingLoad(2, grid=10, force_magnitude=1.0),
            Unbalance(3, grid=10, mass_eccentricity=1.0e-3),
        ]
    )
    solver = Solver(model)
    frequency_options = {
        "fixed_dofs": [2, 3, 4, 5],
        "frequencies": [5.0],
    }

    with pytest.raises(ValueError, match="exactly one"):
        solver.solve("frequency_response", **frequency_options)
    with pytest.raises(ValueError, match="exactly one"):
        solver.solve(
            "frequency_response",
            force=np.zeros(6),
            load_ids=[1],
            **frequency_options,
        )
    with pytest.raises(ValueError, match="at least one"):
        solver.solve(
            "frequency_response", load_ids=[], **frequency_options
        )
    with pytest.raises(TypeError, match="sequence"):
        solver.solve(
            "frequency_response", load_ids=1, **frequency_options
        )
    with pytest.raises(TypeError, match="only integers"):
        solver.solve(
            "frequency_response", load_ids=[True], **frequency_options
        )
    with pytest.raises(ValueError, match="duplicate"):
        solver.solve(
            "frequency_response", load_ids=[1, 1], **frequency_options
        )
    with pytest.raises(ValueError, match="unknown load ID"):
        solver.solve(
            "frequency_response", load_ids=[99], **frequency_options
        )
    with pytest.raises(ValueError, match="not Unbalance"):
        solver.solve(
            "frequency_response", load_ids=[3], **frequency_options
        )
    with pytest.raises(ValueError, match="not HarmonicLoad"):
        solver.solve(
            "unbalance_response",
            fixed_dofs=[2, 3, 4, 5],
            speeds=[10.0],
            load_ids=[1],
        )
    with pytest.raises(ValueError, match="nonnegative"):
        solver.solve(
            "frequency_response",
            frequencies=[-1.0],
            load_ids=[1],
            fixed_dofs=[2, 3, 4, 5],
        )
    with pytest.raises(ValueError, match="speeds must"):
        solver.solve(
            "unbalance_response",
            fixed_dofs=[2, 3, 4, 5],
            speeds=np.nan,
            load_ids=[3],
        )
