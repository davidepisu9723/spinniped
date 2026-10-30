"""Declarative, NASTRAN-like records used to define a Spinniped model.

All record fields contain only Python scalars, tuples, lists, dictionaries,
dataclasses, or ``None``. A numeric field may reference a registered random
distribution with ``(distribution_id,)`` or one multivariate component with
``(distribution_id, component)``. Selected linear coefficients may instead be
tabulated against rotor speed with :class:`SpeedDependentValue`.
"""

from dataclasses import dataclass, field
from math import isfinite
from numbers import Real
from typing import TypeAlias

import numpy as np


RandomReference: TypeAlias = tuple[int] | tuple[int, int]
Parameter: TypeAlias = int | float | RandomReference
Vector3: TypeAlias = tuple[int | float, int | float, int | float]


def _validate_parameter(value, label):
    """Validate one finite scalar or structurally valid random reference."""
    if isinstance(value, tuple):
        if (
            len(value) not in (1, 2)
            or any(
                isinstance(item, bool) or not isinstance(item, int)
                for item in value
            )
        ):
            raise TypeError(f"{label} must be a real number or random reference")
        if value[0] <= 0 or (len(value) == 2 and value[1] < 0):
            raise ValueError(
                f"{label} random reference requires a positive distribution "
                "ID and a nonnegative component index"
            )
    elif isinstance(value, bool) or not isinstance(value, Real):
        raise TypeError(f"{label} must be a real number or random reference")
    elif not isfinite(value):
        raise ValueError(f"{label} must be finite")


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
        """Validate the complete local speed-table contract."""
        # Tables are immutable paired sequences, so their structure can be
        # rejected immediately without any surrounding model information.
        if not isinstance(self.speeds, tuple):
            raise TypeError("speeds must be a tuple")
        if not isinstance(self.values, tuple):
            raise TypeError("values must be a tuple")
        if len(self.speeds) < 2:
            raise ValueError("SpeedDependentValue requires at least two speeds")
        if len(self.values) != len(self.speeds):
            raise ValueError("speeds and values must have the same length")

        # Knots are deterministic. Validate their numerical domain and order
        # before the table can be placed in another record.
        if any(
            isinstance(speed, bool) or not isinstance(speed, Real)
            for speed in self.speeds
        ):
            raise TypeError("speeds must contain only real numeric values")
        if not all(isfinite(speed) for speed in self.speeds):
            raise ValueError("speeds must contain only finite values")
        if any(
            current <= previous
            for previous, current in zip(
                self.speeds,
                self.speeds[1:],
            )
        ):
            raise ValueError("speeds must be strictly increasing")

        # Values may be finite numbers or unresolved random references.
        # Reference existence and upper component bounds require ModelBuilder.
        for index, value in enumerate(self.values):
            _validate_parameter(value, f"values[{index}]")


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

    def __post_init__(self):
        """Validate the complete local distribution specification."""
        # Identity and naming are intrinsic to the record. Registry-wide
        # uniqueness remains the responsibility of ModelBuilder.
        if isinstance(self.id, bool) or not isinstance(self.id, int):
            raise TypeError("Distribution IDs must be integers")
        if self.id <= 0:
            raise ValueError("Distribution IDs must be positive")
        if not isinstance(self.name, str) or not self.name.strip():
            raise ValueError("Distribution names must be non-empty strings")
        if not isinstance(self.distribution, str):
            raise TypeError("Distribution type must be a string")
        if not isinstance(self.parameters, dict):
            raise TypeError("Distribution parameters must be a dictionary")

        def finite_array(value, label, dimensions):
            """Return one finite numeric parameter array."""
            try:
                array = np.asarray(value, dtype=float)
            except (TypeError, ValueError) as error:
                raise TypeError(
                    f"{label} must contain only numeric values"
                ) from error
            if array.ndim not in dimensions:
                raise ValueError(f"{label} has an invalid shape")
            if not np.isfinite(array).all():
                raise ValueError(f"{label} must contain only finite values")
            return array

        # Family schemas and parameter relationships require no model context,
        # so reject malformed specifications before registration or sampling.
        family = self.distribution.lower()
        parameters = self.parameters
        if family == "normal":
            if set(parameters) != {"mean", "stdv"}:
                raise ValueError(
                    "A normal distribution requires only 'mean' and 'stdv'"
                )
            finite_array(parameters["mean"], "mean", (0,))
            deviation = finite_array(parameters["stdv"], "stdv", (0,))
            if deviation < 0:
                raise ValueError(
                    "Distribution standard deviation cannot be negative"
                )
        elif family == "uniform":
            if set(parameters) != {"low", "high"}:
                raise ValueError(
                    "A uniform distribution requires only 'low' and 'high'"
                )
            lower = finite_array(parameters["low"], "low", (0,))
            upper = finite_array(parameters["high"], "high", (0,))
            if upper < lower:
                raise ValueError("Uniform distribution requires high >= low")
        elif family == "multivariate_normal":
            required = {"mean", "stdv", "correlation"}
            if set(parameters) != required:
                raise ValueError(
                    "A multivariate normal distribution requires only "
                    "'mean', 'stdv', and 'correlation'"
                )
            mean = finite_array(parameters["mean"], "mean", (1,))
            deviation = finite_array(parameters["stdv"], "stdv", (1,))
            correlation = finite_array(
                parameters["correlation"], "correlation", (2,)
            )
            component_count = mean.size
            if component_count == 0:
                raise ValueError(
                    "A multivariate distribution needs at least one component"
                )
            if deviation.shape != mean.shape:
                raise ValueError("mean and stdv must have the same length")
            if correlation.shape != (component_count, component_count):
                raise ValueError(
                    "correlation shape must match the number of components"
                )
            if np.any(deviation < 0):
                raise ValueError(
                    "Distribution standard deviations cannot be negative"
                )
            if not np.allclose(correlation, correlation.T):
                raise ValueError("correlation must be symmetric")
            if not np.allclose(np.diag(correlation), 1.0):
                raise ValueError("correlation diagonal entries must equal one")
            if np.any(np.abs(correlation) > 1.0):
                raise ValueError(
                    "correlation coefficients must lie between -1 and 1"
                )
            if np.linalg.eigvalsh(correlation).min() < -1e-12:
                raise ValueError("correlation must be positive semidefinite")
        else:
            raise ValueError(f"Unsupported distribution {family!r}")


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

    def __post_init__(self):
        """Validate the complete local grid specification."""
        # Grid identity is intrinsic; collection-wide uniqueness is checked
        # when ModelBuilder registers the record.
        if isinstance(self.id, bool) or not isinstance(self.id, int):
            raise TypeError("Grid IDs must be integers")
        if self.id <= 0:
            raise ValueError("Grid IDs must be positive")

        # Coordinates may remain symbolic random references until build time,
        # but their scalar/reference structure is locally decidable.
        for name in ("x", "y", "z"):
            _validate_parameter(getattr(self, name), f"Grid {self.id} {name}")

        # Frame zero is the global system; user-defined frame IDs are positive.
        if isinstance(self.coordinate_system, bool) or not isinstance(
            self.coordinate_system, int
        ):
            raise TypeError("Grid coordinate_system must be an integer ID")
        if self.coordinate_system < 0:
            raise ValueError("Grid coordinate_system cannot be negative")


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
    origin: Vector3 = (0.0, 0.0, 0.0)
    x_axis: Vector3 = (1.0, 0.0, 0.0)
    xy_plane: Vector3 = (0.0, 1.0, 0.0)

    def __post_init__(self):
        """Validate the complete local coordinate-system specification."""
        if isinstance(self.id, bool) or not isinstance(self.id, int):
            raise TypeError("Coordinate system IDs must be integers")
        if self.id <= 0:
            raise ValueError("Coordinate system IDs must be positive")

        # Reference frames are deterministic immutable triples. Stochastic
        # orientation would implicitly randomize every attached model record.
        if not isinstance(self.origin, tuple):
            raise TypeError("origin must be a tuple")
        if not isinstance(self.x_axis, tuple):
            raise TypeError("x_axis must be a tuple")
        if not isinstance(self.xy_plane, tuple):
            raise TypeError("xy_plane must be a tuple")
        for name in ("origin", "x_axis", "xy_plane"):
            vector = getattr(self, name)
            if len(vector) != 3:
                raise ValueError(f"{name} must contain three components")
            for index, component in enumerate(vector):
                label = f"{name}[{index}]"
                if isinstance(component, bool) or not isinstance(
                    component, Real
                ):
                    raise TypeError(f"{label} must be a real number")
                if not isfinite(component):
                    raise ValueError(f"{label} must be finite")

        # With deterministic components, frame geometry is fully decidable.
        x_direction = np.asarray(self.x_axis, dtype=float)
        xy_direction = np.asarray(self.xy_plane, dtype=float)
        x_norm = np.linalg.norm(x_direction)
        if x_norm < 1.0e-14:
            raise ValueError("Coordinate system has a zero-length x axis")
        cross_norm = np.linalg.norm(
            np.cross(x_direction / x_norm, xy_direction)
        )
        if cross_norm < 1.0e-14:
            raise ValueError("Coordinate system has collinear axes")


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

    def __post_init__(self):
        """Validate the complete local isotropic-material specification."""
        if isinstance(self.id, bool) or not isinstance(self.id, int):
            raise TypeError("Material IDs must be integers")
        if self.id <= 0:
            raise ValueError("Material IDs must be positive")

        # Parameters may remain random references. Reconstructed sampled
        # records pass through the numerical bounds below automatically.
        for name in ("density", "young_modulus", "poisson_ratio"):
            _validate_parameter(getattr(self, name), f"Material {self.id} {name}")
        if isinstance(self.density, Real) and self.density <= 0.0:
            raise ValueError("Material density must be positive")
        if isinstance(self.young_modulus, Real) and self.young_modulus <= 0.0:
            raise ValueError("Material young_modulus must be positive")
        if isinstance(self.poisson_ratio, Real) and not (
            -1.0 < self.poisson_ratio < 0.5
        ):
            raise ValueError(
                "Material poisson_ratio must lie strictly between -1 and 0.5"
            )


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

    def __post_init__(self):
        """Validate the complete local shaft-property specification."""
        if isinstance(self.id, bool) or not isinstance(self.id, int):
            raise TypeError("Shaft property IDs must be integers")
        if self.id <= 0:
            raise ValueError("Shaft property IDs must be positive")
        if isinstance(self.material, bool) or not isinstance(self.material, int):
            raise TypeError("Shaft property material must be an integer ID")
        if self.material <= 0:
            raise ValueError("Shaft property material ID must be positive")

        # Sampled properties are reconstructed and rechecked against the same
        # annular-section and damping bounds.
        for name in ("outer_diameter", "inner_diameter", "damping"):
            _validate_parameter(getattr(self, name), f"Shaft property {self.id} {name}")
        if isinstance(self.outer_diameter, Real) and self.outer_diameter <= 0.0:
            raise ValueError("Shaft outer_diameter must be positive")
        if isinstance(self.inner_diameter, Real) and self.inner_diameter < 0.0:
            raise ValueError("Shaft inner_diameter cannot be negative")
        if isinstance(self.damping, Real) and self.damping < 0.0:
            raise ValueError("Shaft damping cannot be negative")
        if isinstance(self.outer_diameter, Real) and isinstance(
            self.inner_diameter, Real
        ) and self.inner_diameter >= self.outer_diameter:
            raise ValueError("Shaft outer_diameter must exceed inner_diameter")
        if not isinstance(self.theory, str):
            raise TypeError("Shaft theory must be a string")
        if self.theory not in {"timoshenko", "euler"}:
            raise ValueError("Shaft theory must be 'timoshenko' or 'euler'")
        if not isinstance(self.rotary_inertia, bool):
            raise TypeError("Shaft rotary_inertia must be a boolean")


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

    def __post_init__(self):
        """Validate that force and moment are tuples."""
        if not isinstance(self.force, tuple):
            raise TypeError("force must be a tuple")
        if not isinstance(self.moment, tuple):
            raise TypeError("moment must be a tuple")


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

    def __post_init__(self):
        """Validate that the spin axis is a tuple."""
        if not isinstance(self.spin_axis, tuple):
            raise TypeError("spin_axis must be a tuple")
