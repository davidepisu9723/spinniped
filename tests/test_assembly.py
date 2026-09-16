"""Tests for declarative validation, coordinate handling, and global assembly."""

from dataclasses import replace

import numpy as np
import pytest

from spinniped import (
    BearingElement,
    BearingProperty,
    CoordinateSystem,
    DiskElement,
    DiskProperty,
    Grid,
    Material,
    ModelBuilder,
    ModelDefinition,
    ShaftElement,
    ShaftProperty,
)
from spinniped.gyroscopic import disk_gyroscopic
from spinniped.mass import disk_mass, shaft_mass
from spinniped.stiffness import bearing_stiffness, shaft_stiffness


def _shaft_records(properties, *, grid_ids=(1, 2), length=None):
    p = properties
    length = p.length if length is None else length
    definition = ModelDefinition(
        grids=[Grid(grid_ids[0], z=0.0), Grid(grid_ids[1], z=length)],
        materials=[Material(7, p.density, p.young_modulus, p.poisson_ratio)],
        properties=[
            ShaftProperty(
                id=9,
                material=7,
                outer_diameter=p.diameter,
                theory="euler",
                rotary_inertia=False,
            )
        ],
        elements=[ShaftElement(11, grid_ids[0], grid_ids[1], 9)],
    )
    return definition


def _local_shaft_matrices(properties, length):
    p = properties
    common = {
        "length": length,
        "outer_diameter": p.diameter,
        "inner_diameter": 0.0,
        "young_modulus": p.young_modulus,
        "poisson_ratio": p.poisson_ratio,
        "theory": "euler",
    }
    return (
        shaft_stiffness(**common),
        shaft_mass(
            **common,
            density=p.density,
            rotary_inertia=False,
        ),
    )


def test_one_element_global_matrices_equal_local_kernels(shaft_properties):
    definition = _shaft_records(shaft_properties)
    model = ModelBuilder().build(definition)
    local_stiffness, local_mass = _local_shaft_matrices(
        shaft_properties, shaft_properties.length
    )

    assert model.grid_ids == (1, 2)
    assert model.K.shape == (1, 12, 12)
    assert np.allclose(model.K[0], local_stiffness)
    assert np.allclose(model.M[0], local_mass)


def test_two_elements_accumulate_at_shared_grid(shaft_properties):
    p = shaft_properties
    element_length = p.length / 2.0
    definition = _shaft_records(p, length=element_length)
    definition = replace(
        definition,
        grids=[Grid(10, z=0.0), Grid(30, z=element_length), Grid(70, z=p.length)],
        elements=[
            ShaftElement(1, 10, 30, 9),
            ShaftElement(2, 30, 70, 9),
        ],
    )
    model = ModelBuilder().build(definition)
    local_stiffness, _ = _local_shaft_matrices(p, element_length)

    assert model.grid_ids == (10, 30, 70)
    assert model.K.shape == (1, 18, 18)
    assert np.allclose(
        model.K[0, 6:12, 6:12],
        local_stiffness[6:12, 6:12] + local_stiffness[:6, :6],
    )
    assert np.count_nonzero(model.K[0, :6, 12:18]) == 0


def test_nonconsecutive_grid_ids_are_mapped_by_record_order(shaft_properties):
    definition = _shaft_records(shaft_properties, grid_ids=(10, 70))
    model = ModelBuilder().build(definition)
    local_stiffness, local_mass = _local_shaft_matrices(
        shaft_properties, shaft_properties.length
    )

    assert model.grid_ids == (10, 70)
    assert np.allclose(model.K[0], local_stiffness)
    assert np.allclose(model.M[0], local_mass)


def test_shaft_gyroscopic_assembly_is_invariant_to_endpoint_order(
    shaft_properties,
):
    forward_definition = _shaft_records(shaft_properties)
    forward_definition = replace(
        forward_definition,
        properties=[
            replace(forward_definition.properties[0], rotary_inertia=True)
        ],
        spin_axis=(0.0, 0.0, 3.0),
    )
    reverse_definition = replace(
        forward_definition,
        elements=[ShaftElement(11, grid_a=2, grid_b=1, property=9)],
    )

    forward = ModelBuilder().build(forward_definition)
    reverse = ModelBuilder().build(reverse_definition)

    assert np.count_nonzero(forward.G) > 0
    assert np.allclose(reverse.G, forward.G)


def test_bearing_terms_are_inserted_at_referenced_grid():
    definition = ModelDefinition(
        grids=[Grid(10), Grid(30), Grid(70)],
        properties=[
            BearingProperty(
                id=4,
                kxx=11.0,
                kxy=12.0,
                kyx=13.0,
                kyy=14.0,
                kzz=15.0,
                cxx=21.0,
                cxy=22.0,
                cyx=23.0,
                cyy=24.0,
                czz=25.0,
            )
        ],
        elements=[BearingElement(id=8, grid=30, property=4)],
    )
    model = ModelBuilder().build(definition)

    expected_stiffness = np.array(
        [[11.0, 12.0, 0.0], [13.0, 14.0, 0.0], [0.0, 0.0, 15.0]]
    )
    expected_damping = np.array(
        [[21.0, 22.0, 0.0], [23.0, 24.0, 0.0], [0.0, 0.0, 25.0]]
    )
    assert np.array_equal(model.K[0, 6:9, 6:9], expected_stiffness)
    assert np.array_equal(model.C[0, 6:9, 6:9], expected_damping)
    assert np.count_nonzero(model.K) == 5
    assert np.count_nonzero(model.M) == 0
    assert np.count_nonzero(model.G) == 0


def test_two_grid_bearing_assembles_equal_and_opposite_blocks():
    prop = BearingProperty(
        id=4,
        kxx=11.0,
        kxy=12.0,
        kyx=13.0,
        kyy=14.0,
        kzz=15.0,
    )
    definition = ModelDefinition(
        grids=[Grid(10), Grid(30)],
        properties=[prop],
        elements=[BearingElement(id=8, grid=10, grid_b=30, property=4)],
    )
    model = ModelBuilder().build(definition)
    local = bearing_stiffness(
        kxx=prop.kxx,
        kxy=prop.kxy,
        kyx=prop.kyx,
        kyy=prop.kyy,
        kzz=prop.kzz,
    )
    expected = np.block([[local, -local], [-local, local]])

    assert np.array_equal(model.K[0], expected)


def test_two_grid_bearing_rejects_identical_grid_ids():
    definition = ModelDefinition(
        grids=[Grid(10)],
        properties=[BearingProperty(id=4, kxx=1.0)],
        elements=[BearingElement(id=8, grid=10, grid_b=10, property=4)],
    )

    with pytest.raises(ValueError, match="distinct grids"):
        ModelBuilder().build(definition)


def test_disk_terms_are_inserted_at_referenced_grid():
    definition = ModelDefinition(
        grids=[Grid(42)],
        properties=[
            DiskProperty(
                id=5,
                mass=2.0,
                diametral_inertia=0.1,
                polar_inertia=0.2,
                damping=0.03,
            )
        ],
        elements=[DiskElement(id=6, grid=42, property=5)],
    )
    model = ModelBuilder().build(definition)
    expected_mass = disk_mass(
        mass=2.0, diametral_inertia=0.1, polar_inertia=0.2
    )

    assert np.array_equal(model.M[0], expected_mass)
    assert np.array_equal(model.C[0], 0.03 * expected_mass)
    assert np.array_equal(model.G[0], disk_gyroscopic(polar_inertia=0.2))
    assert np.count_nonzero(model.K) == 0


def test_disk_matrices_follow_the_element_coordinate_system():
    definition = ModelDefinition(
        grids=[Grid(42)],
        coordinate_systems=[
            # Local axes (x, y, z) point along global (y, z, x), so the
            # disk's polar axis is global x.
            CoordinateSystem(
                id=2,
                origin=(10.0, -4.0, 3.0),
                x_axis=(0.0, 1.0, 0.0),
                xy_plane=(0.0, 0.0, 1.0),
            )
        ],
        properties=[
            DiskProperty(
                id=5,
                mass=2.0,
                diametral_inertia=0.1,
                polar_inertia=0.2,
            )
        ],
        elements=[DiskElement(id=6, grid=42, property=5, coordinate_system=2)],
        spin_axis=(1.0, 0.0, 0.0),
    )
    model = ModelBuilder().build(definition)

    assert np.allclose(np.diag(model.M[0]), [2.0, 2.0, 2.0, 0.2, 0.1, 0.1])
    assert model.G[0, 4, 5] == pytest.approx(0.2)
    assert model.G[0, 5, 4] == pytest.approx(-0.2)
    assert np.count_nonzero(model.G[0]) == 2


def test_grid_coordinates_are_transformed_to_basic_system():
    definition = ModelDefinition(
        coordinate_systems=[
            CoordinateSystem(
                id=2,
                origin=(10.0, -4.0, 3.0),
                # Axis fields are direction vectors, not points measured
                # from the origin. Local (x, y, z) maps to global (y, -x, z).
                x_axis=(0.0, 1.0, 0.0),
                xy_plane=(-1.0, 0.0, 0.0),
            )
        ],
        grids=[Grid(id=10, x=2.0, y=3.0, z=4.0, coordinate_system=2)],
    )
    model = ModelBuilder().build(definition)

    assert np.array_equal(model.coordinates, [[[7.0, -2.0, 7.0]]])


@pytest.mark.parametrize(
    ("definition", "message"),
    [
        (
            ModelDefinition(grids=[Grid(1), Grid(1)]),
            "Duplicate grid ID 1",
        ),
        (
            ModelDefinition(
                grids=[Grid(1), Grid(2)],
                properties=[ShaftProperty(1, material=99, outer_diameter=0.01)],
                elements=[ShaftElement(1, 1, 2, 1)],
            ),
            "unknown material 99",
        ),
        (
            ModelDefinition(
                grids=[Grid(1)],
                properties=[BearingProperty(1)],
                elements=[ShaftElement(1, 1, 99, 1)],
            ),
            "requires a ShaftProperty",
        ),
        (
            ModelDefinition(
                grids=[Grid(1)],
                properties=[BearingProperty(1)],
                elements=[BearingElement(1, grid=99, property=1)],
            ),
            "unknown grid 99",
        ),
    ],
)
def test_invalid_record_references_are_rejected(definition, message):
    with pytest.raises((TypeError, ValueError), match=message):
        ModelBuilder().build(definition)


def test_coincident_shaft_grids_are_rejected(shaft_properties):
    definition = _shaft_records(shaft_properties)
    definition = replace(definition, grids=[Grid(1), Grid(2)])

    with pytest.raises(ValueError, match="coincident"):
        ModelBuilder().build(definition)


@pytest.mark.parametrize(
    ("definition", "exception", "message"),
    [
        (ModelDefinition(), ValueError, "at least one grid"),
        (ModelDefinition(grids=[Grid(True)]), TypeError, "IDs must be integers"),
        (ModelDefinition(grids=[Grid(0)]), ValueError, "IDs must be positive"),
        (
            ModelDefinition(grids=[Grid(1)], materials=[Grid(2)]),
            TypeError,
            "Every material record",
        ),
        (
            ModelDefinition(grids=[Grid(1, x=np.inf)]),
            ValueError,
            "coordinates must be finite",
        ),
        (
            ModelDefinition(grids=[Grid(1)], spin_axis=(0.0, 0.0, 0.0)),
            ValueError,
            "spin_axis cannot be zero",
        ),
    ],
)
def test_builder_strictly_validates_record_types_and_finite_geometry(
    definition, exception, message
):
    with pytest.raises(exception, match=message):
        ModelBuilder().build(definition)


@pytest.mark.parametrize("stochastic", [0, 1, None, "yes"])
def test_builder_requires_a_boolean_stochastic_flag(stochastic):
    with pytest.raises(TypeError, match="must be a boolean"):
        ModelBuilder().build(
            ModelDefinition(grids=[Grid(1)]), stochastic=stochastic
        )
