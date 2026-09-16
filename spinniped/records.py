"""Declarative, NASTRAN-like records used to define a Spinniped model.

All record fields contain only Python scalars, tuples, lists, dictionaries, or
``None``.  A numeric field may be replaced by a distribution dictionary such
as ``{"distribution": "normal", "mean": 1.0, "std": 0.1}``.
"""

from dataclasses import dataclass, field
from typing import TypeAlias


Parameter: TypeAlias = int | float | dict[str, object]


@dataclass(frozen=True, slots=True)
class Grid:
    """Define one finite-element grid point.

    Parameters
    ----------
    id : int
        Unique grid identifier.
    x, y, z : int, float, or dict, optional
        Coordinates in the selected input coordinate system.
    coordinate_system : int, optional
        ID of the coordinate system in which coordinates are expressed.
    """

    # Grid IDs are referenced by element connectivity records.
    id: int

    # Coordinates may be deterministic scalars or distribution dictionaries.
    x: Parameter = 0.0
    y: Parameter = 0.0
    z: Parameter = 0.0

    # Coordinate system zero denotes the implicit global frame.
    coordinate_system: int = 0


@dataclass(frozen=True, slots=True)
class CoordinateSystem:
    """Define a right-handed Cartesian coordinate system.

    Parameters
    ----------
    id : int
        Unique coordinate-system identifier. Zero is reserved globally.
    origin : tuple, optional
        Global coordinates of the local origin.
    x_axis : tuple, optional
        Direction of the local x axis in global components.
    xy_plane : tuple, optional
        Direction lying in the desired local xy plane.
    """

    # The builder converts these pure-Python vectors into an orthonormal basis.
    id: int
    origin: tuple[Parameter, Parameter, Parameter] = (0.0, 0.0, 0.0)
    x_axis: tuple[Parameter, Parameter, Parameter] = (1.0, 0.0, 0.0)
    xy_plane: tuple[Parameter, Parameter, Parameter] = (0.0, 1.0, 0.0)


@dataclass(frozen=True, slots=True)
class Material:
    """Define isotropic material properties.

    Parameters
    ----------
    id : int
        Unique material identifier.
    density : int, float, or dict
        Mass density.
    young_modulus : int, float, or dict
        Young's modulus.
    poisson_ratio : int, float, or dict
        Poisson's ratio.
    """

    # Each numerical field may be deterministic or stochastic.
    id: int
    density: Parameter
    young_modulus: Parameter
    poisson_ratio: Parameter


@dataclass(frozen=True, slots=True)
class ShaftProperty:
    """Define reusable circular-shaft section data.

    Parameters
    ----------
    id : int
        Unique property identifier.
    material : int
        Referenced :class:`Material` ID.
    outer_diameter, inner_diameter : int, float, or dict
        Annular section diameters.
    theory : {"timoshenko", "euler"}, optional
        Beam theory used by shaft kernels.
    rotary_inertia : bool, optional
        Include bending rotary inertia and shaft gyroscopic terms.
    damping : int, float, or dict, optional
        Mass-proportional damping coefficient.
    """

    # Properties are shared by every shaft element that references their ID.
    id: int
    material: int
    outer_diameter: Parameter
    inner_diameter: Parameter = 0.0
    theory: str = "timoshenko"
    rotary_inertia: bool = True
    damping: Parameter = 0.0


@dataclass(frozen=True, slots=True)
class BearingProperty:
    """Define translational bearing stiffness and damping coefficients.

    Parameters
    ----------
    id : int
        Unique property identifier.
    kxx, kyy, kzz, kxy, kyx : int, float, or dict, optional
        Direct and cross-coupled stiffness coefficients.
    cxx, cyy, czz, cxy, cyx : int, float, or dict, optional
        Direct and cross-coupled viscous damping coefficients.
    """

    # Coefficients are expressed in the bearing element's coordinate system.
    id: int
    kxx: Parameter = 0.0
    kyy: Parameter = 0.0
    kzz: Parameter = 0.0
    kxy: Parameter = 0.0
    kyx: Parameter = 0.0
    cxx: Parameter = 0.0
    cyy: Parameter = 0.0
    czz: Parameter = 0.0
    cxy: Parameter = 0.0
    cyx: Parameter = 0.0


@dataclass(frozen=True, slots=True)
class DiskProperty:
    """Define rigid-disk or lumped-mass inertia data.

    Parameters
    ----------
    id : int
        Unique property identifier.
    mass : int, float, or dict
        Translational mass.
    diametral_inertia : int, float, or dict, optional
        Mass moment about either transverse local axis.
    polar_inertia : int, float, or dict, optional
        Mass moment about the local spin axis.
    damping : int, float, or dict, optional
        Mass-proportional damping coefficient.
    """

    # A disk property contains inertia data but no stiffness.
    id: int
    mass: Parameter
    diametral_inertia: Parameter = 0.0
    polar_inertia: Parameter = 0.0
    damping: Parameter = 0.0


# A lumped mass uses the same numerical formulation as a rigid disk.
LumpedMassProperty = DiskProperty


@dataclass(frozen=True, slots=True)
class ShaftElement:
    """Connect two grids using one shaft property.

    Parameters
    ----------
    id : int
        Unique element identifier.
    grid_a, grid_b : int
        Endpoint grid IDs.
    property : int
        Referenced :class:`ShaftProperty` ID.
    """

    # Connectivity order defines the element's temporary local positive axis.
    id: int
    grid_a: int
    grid_b: int
    property: int


@dataclass(frozen=True, slots=True)
class BearingElement:
    """Connect a bearing to ground or between two grids.

    Parameters
    ----------
    id : int
        Unique element identifier.
    grid : int
        First connected grid ID.
    property : int
        Referenced :class:`BearingProperty` ID.
    grid_b : int or None, optional
        Second grid ID. ``None`` represents a bearing connected to ground.
    coordinate_system : int, optional
        Frame in which bearing coefficients are defined.
    """

    # A missing second grid produces a conventional support-to-ground element.
    id: int
    grid: int
    property: int
    grid_b: int | None = None
    coordinate_system: int = 0


@dataclass(frozen=True, slots=True)
class DiskElement:
    """Attach a rigid disk or lumped mass to one grid.

    Parameters
    ----------
    id : int
        Unique element identifier.
    grid : int
        Connected grid ID.
    property : int
        Referenced :class:`DiskProperty` ID.
    coordinate_system : int, optional
        Frame whose local z direction is the disk polar axis.
    """

    # Disk matrices are rotated from this input frame during assembly.
    id: int
    grid: int
    property: int
    coordinate_system: int = 0


# A lumped-mass element is an alias with identical storage and behavior.
LumpedMassElement = DiskElement

# These unions document valid contents of model property and element lists.
Property: TypeAlias = ShaftProperty | BearingProperty | DiskProperty
Element: TypeAlias = ShaftElement | BearingElement | DiskElement


@dataclass(frozen=True, slots=True)
class ModelDefinition:
    """Group every declarative record required to build a model.

    Parameters
    ----------
    grids : list of Grid, optional
        Grid records in desired global matrix order.
    coordinate_systems : list of CoordinateSystem, optional
        User-defined input frames.
    materials : list of Material, optional
        Isotropic material records.
    properties : list of property records, optional
        Shaft, bearing, and disk property records.
    elements : list of element records, optional
        Shaft, bearing, and disk connectivity records.
    spin_axis : tuple, optional
        Global rotor spin direction. The builder normalizes this vector.
    """

    # Default factories prevent model instances from sharing mutable lists.
    grids: list[Grid] = field(default_factory=list)
    coordinate_systems: list[CoordinateSystem] = field(default_factory=list)
    materials: list[Material] = field(default_factory=list)
    properties: list[Property] = field(default_factory=list)
    elements: list[Element] = field(default_factory=list)
    spin_axis: tuple[Parameter, Parameter, Parameter] = (0.0, 0.0, 1.0)
