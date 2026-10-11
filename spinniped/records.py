"""Declarative, NASTRAN-like records used to define a Spinniped model.

All record fields contain only Python scalars, tuples, lists, dictionaries,
dataclasses, or ``None``. A numeric field may reference a registered random
distribution with ``(distribution_id,)`` or one multivariate component with
``(distribution_id, component)``. Selected linear coefficients may instead be
tabulated against rotor speed with :class:`SpeedDependentValue`.
"""

from dataclasses import dataclass, field
from typing import TypeAlias


RandomReference: TypeAlias = tuple[int] | tuple[int, int]
Parameter: TypeAlias = int | float | RandomReference
Vector3: TypeAlias = tuple[Parameter, Parameter, Parameter]


@dataclass(frozen=True, slots=True)
class SpeedDependentValue:
    """Define one scalar value at tabulated rotor angular speeds.

    Parameters
    ----------
    speeds : tuple of int or float
        Strictly increasing rotor angular speeds in radians per second.
    values : tuple of int, float, or random reference
        Value at each corresponding speed. Random references are resolved once
        per model realization before interpolation.

    Notes
    -----
    Interpolation is linear and extrapolation is not permitted. This wrapper is
    accepted only by explicitly documented linear bearing and load fields.
    """

    speeds: tuple[int | float, ...]
    values: tuple[Parameter, ...]

    def __post_init__(self):
        """Snapshot ordinary input sequences to preserve value semantics.

        Invalid non-sequence inputs are left in place so the builder can
        report them together with the field path. Strings are likewise left
        for the dedicated numeric-table validator to reject.
        """
        for name in ("speeds", "values"):
            sequence = getattr(self, name)
            if isinstance(sequence, (str, bytes)):
                continue
            try:
                normalized = tuple(sequence)
            except TypeError:
                continue
            object.__setattr__(self, name, normalized)


SpeedDependentParameter: TypeAlias = Parameter | SpeedDependentValue
SpeedDependentVector3: TypeAlias = tuple[
    SpeedDependentParameter,
    SpeedDependentParameter,
    SpeedDependentParameter,
]


@dataclass(frozen=True, slots=True)
class RandomDistribution:
    """Define a named scalar or multivariate random distribution.

    Parameters
    ----------
    id : int
        Unique integer used by parameter references.
    name : str
        Unique, human-readable distribution name.
    distribution : {"normal", "uniform", "multivariate_normal"}
        Probability distribution family.
    parameters : dict
        Distribution parameters. Normal variables require ``mean`` and
        ``stdv``; uniform variables require ``low`` and ``high``;
        multivariate normal variables require ``mean``, ``stdv``, and
        ``correlation``.
    """

    # IDs make references compact, while names keep reports understandable.
    id: int
    name: str
    distribution: str
    parameters: dict[str, object]


@dataclass(frozen=True, slots=True)
class Grid:
    """Define one finite-element grid point.

    Parameters
    ----------
    id : int
        Unique grid identifier.
    x, y, z : int, float, or tuple, optional
        Coordinates in the selected input coordinate system.
    coordinate_system : int, optional
        ID of the coordinate system in which coordinates are expressed.
    """

    # Grid IDs are referenced by element connectivity records.
    id: int

    # Coordinates may be deterministic scalars or distribution references.
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
    density : int, float, or tuple
        Mass density.
    young_modulus : int, float, or tuple
        Young's modulus.
    poisson_ratio : int, float, or tuple
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
    outer_diameter, inner_diameter : int, float, or tuple
        Annular section diameters.
    theory : {"timoshenko", "euler"}, optional
        Beam theory used by shaft kernels.
    rotary_inertia : bool, optional
        Include bending rotary inertia and shaft gyroscopic terms.
    damping : int, float, or tuple, optional
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
    kxx, kyy, kzz, kxy, kyx : int, float, tuple, or SpeedDependentValue, optional
        Direct and cross-coupled stiffness coefficients.
    cxx, cyy, czz, cxy, cyx : int, float, tuple, or SpeedDependentValue, optional
        Direct and cross-coupled viscous damping coefficients.
    """

    # Coefficients are expressed in the bearing element's coordinate system.
    id: int
    kxx: SpeedDependentParameter = 0.0
    kyy: SpeedDependentParameter = 0.0
    kzz: SpeedDependentParameter = 0.0
    kxy: SpeedDependentParameter = 0.0
    kyx: SpeedDependentParameter = 0.0
    cxx: SpeedDependentParameter = 0.0
    cyy: SpeedDependentParameter = 0.0
    czz: SpeedDependentParameter = 0.0
    cxy: SpeedDependentParameter = 0.0
    cyx: SpeedDependentParameter = 0.0


@dataclass(frozen=True, slots=True)
class DiskProperty:
    """Define rigid-disk or lumped-mass inertia data.

    Parameters
    ----------
    id : int
        Unique property identifier.
    mass : int, float, or tuple
        Translational mass.
    diametral_inertia : int, float, or tuple, optional
        Mass moment about either transverse local axis.
    polar_inertia : int, float, or tuple, optional
        Mass moment about the local spin axis.
    damping : int, float, or tuple, optional
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


# Loads are independent records even where their physical behavior is related.
# Keeping them as peers makes their different amplitude laws explicit: a
# harmonic load has signed components, a rotating load has a circular-force
# magnitude, and unbalance scales with rotor speed squared. The first two may
# select their amplitudes from a rotor-speed table.
@dataclass(frozen=True, slots=True)
class HarmonicLoad:
    """Define a harmonic nodal load.

    Parameters
    ----------
    id : int
        Unique load identifier.
    grid : int
        ID of the grid where the load is applied.
    force : tuple, optional
        Signed peak force amplitudes ``(fx, fy, fz)`` in the selected
        coordinate system, in force units. Each component may be a
        :class:`SpeedDependentValue`.
    moment : tuple, optional
        Signed peak moment amplitudes ``(mx, my, mz)`` in the selected
        coordinate system, in force-length units. Each component may be a
        :class:`SpeedDependentValue`.
    phase : int, float, or tuple, optional
        Common phase angle in radians.
    coordinate_system : int, optional
        ID of the frame in which force and moment components are expressed.

    Notes
    -----
    Harmonic-response assembly uses the convention
    ``p(t) = Re(p_hat * exp(i * omega * t))``. All six components share
    ``phase``; multiple records can represent components with different
    phases.
    """

    id: int
    grid: int
    force: SpeedDependentVector3 = (0.0, 0.0, 0.0)
    moment: SpeedDependentVector3 = (0.0, 0.0, 0.0)
    phase: Parameter = 0.0
    coordinate_system: int = 0


@dataclass(frozen=True, slots=True)
class RotatingLoad:
    """Define a constant-magnitude circular rotating force.

    Parameters
    ----------
    id : int
        Unique load identifier.
    grid : int
        ID of the grid where the force is applied.
    force_magnitude : int, float, tuple, or SpeedDependentValue
        Nonnegative peak radial force magnitude.
    phase : int, float, or tuple, optional
        Signed initial geometric angle in radians from the selected frame's
        positive x-axis toward its positive y-axis.
    sense : {"positive", "negative"}, optional
        Rotation sense about the selected frame's positive z-axis. Positive
        rotation proceeds from positive x toward positive y.
    coordinate_system : int, optional
        ID of the frame whose x-y plane contains the rotating force.
    """

    id: int
    grid: int
    force_magnitude: SpeedDependentParameter
    phase: Parameter = 0.0
    sense: str = "positive"
    coordinate_system: int = 0


@dataclass(frozen=True, slots=True)
class Unbalance:
    """Define a synchronous mass-unbalance excitation.

    Parameters
    ----------
    id : int
        Unique load identifier.
    grid : int
        ID of the grid carrying the unbalance.
    mass_eccentricity : int, float, or tuple
        Nonnegative unbalance magnitude ``m * e`` in mass-length units.
    phase : int, float, or tuple, optional
        Signed initial eccentricity angle in radians from the selected frame's
        positive x-axis toward its positive y-axis.
    coordinate_system : int, optional
        ID of the frame whose x-y plane is the unbalance plane.

    Notes
    -----
    Unbalance is a synchronous rotating force with magnitude
    ``mass_eccentricity * speed**2``. Its rotation follows the rotor rather
    than defining an independent sense.
    """

    id: int
    grid: int
    mass_eccentricity: Parameter
    phase: Parameter = 0.0
    coordinate_system: int = 0


# These unions document valid contents of model record lists.
Property: TypeAlias = ShaftProperty | BearingProperty | DiskProperty
Element: TypeAlias = ShaftElement | BearingElement | DiskElement
Load: TypeAlias = HarmonicLoad | RotatingLoad | Unbalance


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
    distributions : list of RandomDistribution, optional
        Registered probability distributions referenced by numeric fields.
    spin_axis : tuple, optional
        Global rotor spin direction. The builder normalizes this vector.
    loads : list of load records, optional
        Harmonic, rotating, and unbalance excitations associated with grids.
    """

    # Default factories prevent model instances from sharing mutable lists.
    grids: list[Grid] = field(default_factory=list)
    coordinate_systems: list[CoordinateSystem] = field(default_factory=list)
    materials: list[Material] = field(default_factory=list)
    properties: list[Property] = field(default_factory=list)
    elements: list[Element] = field(default_factory=list)
    distributions: list[RandomDistribution] = field(default_factory=list)
    spin_axis: tuple[Parameter, Parameter, Parameter] = (0.0, 0.0, 1.0)
    # Appended to preserve the positional order of the pre-existing fields.
    loads: list[Load] = field(default_factory=list)
