# Chapter 2 -- Model-definition API

## Table of contents

- [How the records fit together](#how-the-records-fit-together)
- [Record reference](#record-reference)
- [Coordinate systems](#coordinate-systems)
- [Spin-axis convention](#spin-axis-convention)

A model definition is data, not a collection of stateful finite-element
objects. Records refer to one another by integer IDs:

- `Grid` defines a point and its input coordinate system.
- `CoordinateSystem` defines a local Cartesian frame.
- `Material` defines density, Young's modulus, and Poisson's ratio.
- `ShaftProperty`, `BearingProperty`, and `DiskProperty` define reusable
  physical data.
- `ShaftElement`, `BearingElement`, and `DiskElement` define connectivity and
  reference a compatible property.
- `SpeedDependentValue` tabulates selected linear coefficients against rotor
  angular speed.
- `HarmonicLoad`, `RotatingLoad`, and `Unbalance` define nodal excitation
  data.
- `ModelDefinition` holds the record lists.

`LumpedMassElement` and `LumpedMassProperty` are aliases for the corresponding
disk records.

## How the records fit together

The left side of this diagram shows the ID references used by one typical
rotor definition. The right side shows when declarative records become
numerical arrays and, finally, analysis results.

```mermaid
flowchart LR
    subgraph records["ID-linked definition records"]
        direction TB
        material["Material 1"] --> shaft_property["ShaftProperty 1"]
        shaft_property --> shaft_element["ShaftElement"]
        distribution["RandomDistribution 17"] --> bearing_property["BearingProperty: kxx = (17,)"]
        speed_table["SpeedDependentValue"] --> bearing_property
        bearing_property --> bearing_element["BearingElement"]
        disk_property["DiskProperty 3"] --> disk_element["DiskElement"]
        grids["Grid IDs"] --> shaft_element
        grids --> bearing_element
        grids --> disk_element
        grids --> loads["HarmonicLoad / RotatingLoad / Unbalance"]
    end

    shaft_element --> definition["ModelDefinition"]
    bearing_element --> definition
    disk_element --> definition
    loads --> definition
    definition --> builder["ModelBuilder<br/>resolve, validate, assemble"]
    builder --> built["BuiltModel<br/>K, M, C, G"]
    built --> solver["Solver"]
    solver --> result["result arrays"]
    result --> plots["plots"]
```

The arrows among definition records are integer references, not copied
numerical data. `ModelBuilder` is the boundary between that declarative graph
and the assembled numerical model.

This example defines a three-grid shaft with a disk, two bearings, and two
loads:

```python
from spinniped import (
    BearingElement,
    BearingProperty,
    DiskElement,
    DiskProperty,
    Grid,
    HarmonicLoad,
    Material,
    ModelDefinition,
    RandomDistribution,
    ShaftElement,
    ShaftProperty,
    SpeedDependentValue,
    Unbalance,
)

bearing_stiffness = SpeedDependentValue(
    speeds=(0.0, 100.0, 200.0),
    values=(1.00e8, 0.95e8, 0.85e8),
)

definition = ModelDefinition(
    grids=[
        Grid(1, z=0.0),
        Grid(2, z=0.5),
        Grid(3, z=1.0),
    ],
    materials=[
        Material(
            id=1,
            density=7850.0,
            young_modulus=(1,),
            poisson_ratio=0.3,
        ),
    ],
    properties=[
        ShaftProperty(
            id=1,
            material=1,
            outer_diameter=(2,),
            theory="timoshenko",
        ),
        BearingProperty(
            id=2,
            kxx=bearing_stiffness,
            kyy=bearing_stiffness,
            cxx=100.0,
            cyy=100.0,
        ),
        DiskProperty(
            id=3,
            mass=2.0,
            diametral_inertia=0.01,
            polar_inertia=0.02,
        ),
    ],
    elements=[
        ShaftElement(1, grid_a=1, grid_b=2, property=1),
        ShaftElement(2, grid_a=2, grid_b=3, property=1),
        BearingElement(3, grid=1, property=2),
        BearingElement(4, grid=3, property=2),
        DiskElement(5, grid=2, property=3),
    ],
    loads=[
        HarmonicLoad(
            id=1,
            grid=2,
            force=(10.0, 0.0, 0.0),
            phase=0.0,
        ),
        Unbalance(
            id=2,
            grid=2,
            mass_eccentricity=2.0e-4,
            phase=0.0,
        ),
    ],
    distributions=[
        RandomDistribution(
            id=1,
            name="steel Young modulus",
            distribution="normal",
            parameters={"mean": 210e9, "stdv": 5e9},
        ),
        RandomDistribution(
            id=2,
            name="shaft diameter",
            distribution="uniform",
            parameters={"low": 0.019, "high": 0.021},
        ),
    ],
    spin_axis=(0.0, 0.0, 1.0),
)
```

## Record reference

All dimensions and coefficients must use one coherent unit system; examples
and theory notes use SI.

| Record | Fields |
|---|---|
| `Grid` | `id`; coordinates `x`, `y`, `z` (default `0.0`); `coordinate_system` (default global ID `0`) |
| `CoordinateSystem` | `id`; `origin`; direction `x_axis`; direction `xy_plane` |
| `Material` | `id`; `density`; `young_modulus`; `poisson_ratio` |
| `ShaftProperty` | `id`; material reference `material`; `outer_diameter`; `inner_diameter` (default `0.0`); `theory`; `rotary_inertia`; `damping` |
| `BearingProperty` | `id`; stiffnesses `kxx`, `kyy`, `kzz`, `kxy`, `kyx`; damping coefficients `cxx`, `cyy`, `czz`, `cxy`, `cyx` |
| `DiskProperty` | `id`; `mass`; `diametral_inertia`; `polar_inertia`; `damping` |
| `ShaftElement` | `id`; end references `grid_a`, `grid_b`; shaft `property` reference |
| `BearingElement` | `id`; first `grid`; bearing `property`; optional `grid_b`; property `coordinate_system` |
| `DiskElement` | `id`; nodal `grid`; disk `property`; property `coordinate_system` |
| `HarmonicLoad` | `id`; nodal `grid`; three-component `force` and `moment`; common `phase`; `coordinate_system` |
| `RotatingLoad` | `id`; nodal `grid`; `force_magnitude`; `phase`; `sense`; `coordinate_system` |
| `Unbalance` | `id`; nodal `grid`; `mass_eccentricity`; `phase`; `coordinate_system` |
| `SpeedDependentValue` | strictly increasing rotor `speeds`; corresponding scalar `values` |
| `RandomDistribution` | unique `id`; unique `name`; `distribution`; family-specific `parameters` |

`ShaftProperty.theory` accepts `"timoshenko"` (the default) and `"euler"`.
Bending rotary inertia is enabled by default. Shaft and disk `damping` values
are coefficients in the mass-proportional relation
$\mathbf{C}=\alpha\mathbf{M}$. Bearing coefficients default to zero and occupy
the translational x-y-z block.

All bearing stiffness and damping coefficients accept a
`SpeedDependentValue`. So do the individual force and moment components of a
`HarmonicLoad` and `RotatingLoad.force_magnitude`. Other fields remain scalar
or random parameters. [Chapter 3](03_speed_dependent_values.md) documents the
supported placements, interpolation rules, and operating-speed convention.

`ModelDefinition` groups the `grids`, `coordinate_systems`, `materials`,
`properties`, `elements`, `distributions`, and `loads` lists; each defaults to
an empty list. Its `spin_axis` field is the global spin direction and defaults
to `(0.0, 0.0, 1.0)`.

Lists and tuples have distinct roles in the declarative API. The fields above
are lists because they are record collections that may be assembled
incrementally. Fixed value sequences must be tuples: `SpeedDependentValue`
`speeds` and `values`, the three `CoordinateSystem` vectors,
`HarmonicLoad.force` and `moment`, and `ModelDefinition.spin_axis`. Mutable
lists are rejected for these atomic value fields rather than converted
implicitly.

IDs must be unique positive integers within each record category. A shaft
element must reference a `ShaftProperty`, a bearing element a
`BearingProperty`, and a disk element a `DiskProperty`. A shaft property also
references a material ID. Harmonic, rotating, and unbalance records share one
load-ID namespace. The builder validates these relationships.

`Grid` validates its own positive integer ID, coordinate parameter structure,
and coordinate-system reference type when it is constructed. Model-wide
checks such as duplicate grid IDs and unknown coordinate-system references
remain builder responsibilities.

## Coordinate systems

`Grid.x`, `Grid.y`, and `Grid.z` are interpreted in the coordinate system
selected by `Grid.coordinate_system`. System `0` is reserved for the global
Cartesian system and must not be redefined. A custom system is defined by a
global origin, a direction for its positive x-axis, and another direction in
its x-y plane. The two direction vectors need not be normalized, but they must
be nonzero and non-collinear:

```python
from spinniped import CoordinateSystem, Grid

frame = CoordinateSystem(
    id=10,
    origin=(1.0, 0.0, 0.0),
    x_axis=(0.0, 1.0, 0.0),
    xy_plane=(0.0, 0.0, 1.0),
)

point = Grid(20, x=0.5, coordinate_system=10)
```

Shaft matrices are formulated in an element-local frame whose z-axis runs
from `grid_a` to `grid_b`, then transformed into the global frame during
assembly. Shaft elements therefore need not be parallel to the global z-axis.

`BearingElement.coordinate_system` rotates bearing coefficients from the
selected local frame into the global frame. With `grid_b=None` (the default),
the bearing connects `grid` to ground. Supplying `grid_b` creates a relative
two-grid bearing with the block pattern

$$
\begin{bmatrix}
\mathbf{A}_b & -\mathbf{A}_b\\
-\mathbf{A}_b & \mathbf{A}_b
\end{bmatrix},
$$

for both bearing stiffness and damping. `DiskElement.coordinate_system`
similarly orients disk inertia and its local spin axis. These element
coordinate-system fields default to global system `0`. A grid's own
`coordinate_system` locates the grid; it does not implicitly orient an attached
bearing or disk.

Load coordinate systems follow the same separation. A grid's coordinate
system locates the application point, while each load's `coordinate_system`
orients its components or rotating plane. `HarmonicLoad.force` and `moment`
are signed peak amplitudes in that frame and share one phase angle in radians.
For `RotatingLoad`, `sense="positive"` means rotation from local positive x
toward local positive y about local positive z; `"negative"` reverses it.
An `Unbalance` frame may point parallel or antiparallel to `spin_axis`, but its
local x-y plane must be the rotor plane.

Loads do not contribute to the structural matrices. Steady response analyses
select records through their shared load IDs and assemble the corresponding
sample-specific complex force vectors.

## Spin-axis convention

`ModelDefinition.spin_axis` is a nonzero direction expressed in global
coordinates. The builder normalizes it, so it defines orientation rather than
speed. For each shaft, the unit-speed gyroscopic matrix is multiplied by the
projection of this direction onto the shaft's local axis. Reversing
`grid_a`/`grid_b` therefore does not reverse the physical global gyroscopic
effect. A shaft perpendicular to the spin axis has zero gyroscopic
contribution in the current formulation.

For a disk, the same projection uses the local z-axis selected by
`DiskElement.coordinate_system`. The scalar `speed` or each Campbell `speeds`
value later supplies the angular-speed magnitude in rad/s.
