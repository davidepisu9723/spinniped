"""Model validation, random sampling, coordinate handling, and assembly."""

from dataclasses import dataclass, field, fields, is_dataclass, replace

import numpy as np

from . import damping, gyroscopic, mass, stiffness
from .records import (
    BearingElement,
    BearingProperty,
    CoordinateSystem,
    DiskElement,
    DiskProperty,
    Grid,
    HarmonicLoad,
    Material,
    ModelDefinition,
    RandomDistribution,
    RotatingLoad,
    ShaftElement,
    ShaftProperty,
    SpeedDependentValue,
    Unbalance,
)


# Every grid owns three translations and three rotations.
DOFS_PER_GRID = 6


@dataclass(frozen=True, slots=True)
class _SpeedDependentMatrixContribution:
    """Store one transformed unit-coefficient matrix contribution."""

    value: SpeedDependentValue
    matrix: np.ndarray
    indices: tuple[int, ...]
    path: str


@dataclass(frozen=True, slots=True)
class BuiltModel:
    """Numerical model or Monte Carlo ensemble produced by ``ModelBuilder``.

    Parameters
    ----------
    stiffness, mass, damping, gyroscopic : numpy.ndarray
        Global matrices with shape ``(samples, ndof, ndof)``. Stiffness and
        damping contain the speed-independent base terms; gyroscopic is the
        unit-speed matrix.
    grid_ids : tuple of int
        Grid IDs in global matrix order.
    coordinates : numpy.ndarray
        Resolved global coordinates with shape ``(samples, grids, 3)``.
    definitions : tuple of ModelDefinition
        Fully resolved deterministic definition for every sample.
    stochastic : bool
        Whether the model represents a sampled stochastic ensemble.
    seed : int or None
        Random seed used during construction.
    _stiffness_speed_contributions, _damping_speed_contributions : tuple
        Internal per-sample unit-matrix terms used to form transient
        speed-dependent matrix deltas. The public equation matrices remain
        three-dimensional.
    """

    stiffness: np.ndarray
    mass: np.ndarray
    damping: np.ndarray
    gyroscopic: np.ndarray
    grid_ids: tuple[int, ...]
    coordinates: np.ndarray
    definitions: tuple[ModelDefinition, ...]
    stochastic: bool
    seed: int | None
    _stiffness_speed_contributions: tuple[
        tuple[_SpeedDependentMatrixContribution, ...], ...
    ] = field(default_factory=tuple, repr=False)
    _damping_speed_contributions: tuple[
        tuple[_SpeedDependentMatrixContribution, ...], ...
    ] = field(default_factory=tuple, repr=False)

    def __post_init__(self):
        """Validate cross-field dimensions and finite numerical contents.

        Raises
        ------
        TypeError
            If matrix, coordinate, definition, or flag types are invalid.
        ValueError
            If shapes disagree or numerical values are non-finite.
        """
        # Group the four equation matrices for uniform validation.
        matrices = (
            self.stiffness,
            self.mass,
            self.damping,
            self.gyroscopic,
        )

        # Built numerical models expose NumPy arrays, never nested lists.
        if any(not isinstance(matrix, np.ndarray) for matrix in matrices):
            raise TypeError("BuiltModel matrices must be NumPy arrays")

        # The leading dimension always represents Monte Carlo samples.
        if any(matrix.ndim != 3 for matrix in matrices):
            raise ValueError(
                "BuiltModel matrices must have shape (samples, ndof, ndof)"
            )

        # K, M, C, and G must describe exactly the same systems.
        if len({matrix.shape for matrix in matrices}) != 1:
            raise ValueError("BuiltModel matrices must all have the same shape")

        # Unpack the shared dimensions for the remaining checks.
        sample_count, row_count, column_count = self.stiffness.shape

        # Global matrices are square and contain six DOFs per grid.
        if (
            sample_count < 1
            or row_count != column_count
            or row_count != DOFS_PER_GRID * len(self.grid_ids)
        ):
            raise ValueError(
                "BuiltModel matrix dimensions are inconsistent with its grids"
            )

        # Coordinates must be stored once for each sample and grid.
        if not isinstance(self.coordinates, np.ndarray):
            raise TypeError("BuiltModel coordinates must be a NumPy array")
        expected_coordinate_shape = (sample_count, len(self.grid_ids), 3)
        if self.coordinates.shape != expected_coordinate_shape:
            raise ValueError("BuiltModel coordinates have an inconsistent shape")

        # Grid IDs define one unique global matrix ordering.
        if (
            not isinstance(self.grid_ids, tuple)
            or len(set(self.grid_ids)) != len(self.grid_ids)
        ):
            raise ValueError("BuiltModel grid_ids must be a tuple of unique IDs")

        # Each numerical realization retains its resolved input records.
        if not isinstance(self.definitions, tuple):
            raise TypeError("BuiltModel definitions must be a tuple")
        if len(self.definitions) != sample_count:
            raise ValueError(
                "BuiltModel requires one resolved definition per sample"
            )

        # Speed-dependent terms are stored per realization while the base
        # matrices retain their established sample-first three-dimensional
        # shape. An empty outer tuple remains valid for manually constructed
        # static BuiltModel instances.
        contribution_groups = (
            self._stiffness_speed_contributions,
            self._damping_speed_contributions,
        )
        for groups in contribution_groups:
            if not isinstance(groups, tuple):
                raise TypeError(
                    "BuiltModel speed-dependent contributions must be tuples"
                )
            if groups and len(groups) != sample_count:
                raise ValueError(
                    "BuiltModel requires one speed-contribution group per sample"
                )
            for sample_terms in groups:
                if not isinstance(sample_terms, tuple):
                    raise TypeError(
                        "BuiltModel sample speed contributions must be tuples"
                    )
                for term in sample_terms:
                    if not isinstance(term, _SpeedDependentMatrixContribution):
                        raise TypeError(
                            "BuiltModel contains an invalid speed contribution"
                        )
                    if (
                        not isinstance(term.matrix, np.ndarray)
                        or term.matrix.ndim != 2
                        or term.matrix.shape[0] != term.matrix.shape[1]
                        or term.matrix.shape[0] != len(term.indices)
                    ):
                        raise ValueError(
                            "Speed contribution matrices must be square and "
                            "match their DOF indices"
                        )
                    if not np.isfinite(term.matrix).all():
                        raise ValueError(
                            "Speed contribution matrices must contain finite values"
                        )
                    if (
                        len(set(term.indices)) != len(term.indices)
                        or any(
                            index < 0 or index >= row_count
                            for index in term.indices
                        )
                    ):
                        raise ValueError(
                            "Speed contribution DOF indices are inconsistent "
                            "with the model"
                        )

        # A manually instantiated BuiltModel must not silently discard tables
        # present on bearing properties. Only ModelBuilder has the geometry
        # and coordinate information needed to create these private bases.
        stiffness_names = ("kxx", "kyy", "kzz", "kxy", "kyx")
        damping_names = ("cxx", "cyy", "czz", "cxy", "cyx")
        for sample, definition in enumerate(self.definitions):
            if (
                not isinstance(definition, ModelDefinition)
                or not isinstance(definition.properties, (list, tuple))
                or not isinstance(definition.elements, (list, tuple))
            ):
                continue
            bearing_properties = {
                item.id: item
                for item in definition.properties
                if isinstance(item, BearingProperty)
            }
            expected_counts = [0, 0]
            for element in definition.elements:
                if not isinstance(element, BearingElement):
                    continue
                property_record = bearing_properties.get(element.property)
                if property_record is None:
                    continue
                expected_counts[0] += sum(
                    isinstance(
                        getattr(property_record, name),
                        SpeedDependentValue,
                    )
                    for name in stiffness_names
                )
                expected_counts[1] += sum(
                    isinstance(
                        getattr(property_record, name),
                        SpeedDependentValue,
                    )
                    for name in damping_names
                )
            actual_counts = [
                len(groups[sample]) if groups else 0
                for groups in contribution_groups
            ]
            if actual_counts != expected_counts:
                raise ValueError(
                    "Speed-dependent bearing definitions must be assembled "
                    "by ModelBuilder"
                )

        # Keep the stochastic flag unambiguous.
        if not isinstance(self.stochastic, bool):
            raise TypeError("BuiltModel stochastic must be a boolean")

        # Solvers require finite matrices and coordinates.
        if any(not np.isfinite(matrix).all() for matrix in matrices):
            raise ValueError("BuiltModel matrices must contain only finite values")
        if not np.isfinite(self.coordinates).all():
            raise ValueError("BuiltModel coordinates must contain only finite values")

    @property
    def samples(self):
        """Return the number of deterministic realizations.

        Returns
        -------
        int
            Size of the leading sample dimension.
        """
        return self.stiffness.shape[0]

    @property
    def ndof(self):
        """Return the number of global degrees of freedom.

        Returns
        -------
        int
            Row count of each global matrix.
        """
        return self.stiffness.shape[1]

    @property
    def K(self):
        """Return the speed-independent global base stiffness matrices.

        Returns
        -------
        numpy.ndarray
            Array with shape ``(samples, ndof, ndof)``.
        """
        return self.stiffness

    @property
    def M(self):
        """Return the global mass matrices.

        Returns
        -------
        numpy.ndarray
            Array with shape ``(samples, ndof, ndof)``.
        """
        return self.mass

    @property
    def C(self):
        """Return the speed-independent global base damping matrices.

        Returns
        -------
        numpy.ndarray
            Array with shape ``(samples, ndof, ndof)``.
        """
        return self.damping

    @property
    def G(self):
        """Return the global unit-speed gyroscopic matrices.

        Returns
        -------
        numpy.ndarray
            Array with shape ``(samples, ndof, ndof)``.
        """
        return self.gyroscopic


def _speed_dependent_arrays(value, path):
    """Convert one already validated, resolved table to numerical arrays."""
    if not isinstance(value, SpeedDependentValue):
        raise TypeError(f"{path} must be a SpeedDependentValue")
    # SpeedDependentValue owns structural and numeric validation. At this
    # stage random ordinates have also been resolved, so only conversion is
    # needed for interpolation and field-specific physical checks.
    return (
        np.asarray(value.speeds, dtype=float),
        np.asarray(value.values, dtype=float),
    )


def _interpolate_speed_dependent_value(value, speeds, path):
    """Linearly evaluate one resolved speed table without extrapolation.

    All requested speeds are checked before interpolation, allowing solvers to
    preflight a complete operating-speed grid before the first numerical solve.
    """
    knots, values = _speed_dependent_arrays(value, path)
    try:
        requested = np.asarray(speeds, dtype=float)
    except (TypeError, ValueError) as error:
        raise TypeError(
            f"{path} evaluation speeds must be real numeric values"
        ) from error
    if requested.ndim != 1 or not requested.size:
        raise ValueError(
            f"{path} evaluation speeds must be a non-empty one-dimensional sequence"
        )
    if not np.isfinite(requested).all():
        raise ValueError(f"{path} evaluation speeds must be finite")
    if np.any(requested < knots[0]) or np.any(requested > knots[-1]):
        raise ValueError(
            f"{path} cannot be evaluated outside its speed range "
            f"[{knots[0]}, {knots[-1]}] rad/s"
        )
    return np.interp(requested, knots, values)


def _contains_speed_dependent_value(value):
    """Return whether an arbitrarily nested record value contains a table."""
    if isinstance(value, SpeedDependentValue):
        return True
    if isinstance(value, (tuple, list)):
        return any(_contains_speed_dependent_value(item) for item in value)
    if isinstance(value, dict):
        return any(_contains_speed_dependent_value(item) for item in value.values())
    return False


def _reject_speed_dependent_value(value, path):
    """Reject a speed table in a field whose formulation is not linear."""
    if _contains_speed_dependent_value(value):
        raise TypeError(f"{path} does not accept SpeedDependentValue")


def _validate_speed_dependent_placements(definition):
    """Restrict speed tables to supported linear bearing and load fields."""
    for distribution in definition.distributions:
        if not isinstance(distribution, RandomDistribution):
            continue
        _reject_speed_dependent_value(
            distribution.parameters,
            f"Distribution {distribution.id} parameters",
        )


def _sample_distribution(specification, rng, stochastic):
    """Generate one realization from a validated distribution record.

    Parameters
    ----------
    specification : RandomDistribution
        Distribution record to validate and sample.
    rng : numpy.random.Generator
        Random generator used for stochastic construction.
    stochastic : bool
        Draw random values when true; otherwise return expected values.

    Returns
    -------
    tuple of float
        One value for a scalar distribution or one value per component.

    """
    # RandomDistribution owns schema validation. Sampling only converts its
    # already valid parameters and selects stochastic draws or expected values.
    distribution_name = specification.distribution.lower()
    parameters = specification.parameters

    if distribution_name == "normal":
        mean = float(parameters["mean"])
        deviation = float(parameters["stdv"])
        draw = rng.normal(mean, deviation) if stochastic else mean
    elif distribution_name == "uniform":
        lower = float(parameters["low"])
        upper = float(parameters["high"])
        draw = (
            rng.uniform(lower, upper)
            if stochastic
            else (lower + upper) / 2.0
        )
    else:
        mean = np.asarray(parameters["mean"], dtype=float)
        deviation = np.asarray(parameters["stdv"], dtype=float)
        correlation = np.asarray(parameters["correlation"], dtype=float)
        covariance = np.outer(deviation, deviation) * correlation
        draw = rng.multivariate_normal(mean, covariance) if stochastic else mean

    # Internally, scalar and vector distributions share one tuple representation.
    values = np.atleast_1d(draw).astype(float)
    if not np.isfinite(values).all():
        raise ValueError("Distribution produced a non-finite value")
    return tuple(values.tolist())


def _sample_distributions(distributions, rng, stochastic):
    """Validate and sample every registered distribution exactly once.

    Parameters
    ----------
    distributions : list of RandomDistribution
        Model-level distribution registry.
    rng : numpy.random.Generator
        Random generator used for stochastic construction.
    stochastic : bool
        Draw random values when true; otherwise return expected values.

    Returns
    -------
    dict
        Mapping from distribution ID to its sampled component tuple.
    """
    sampled = {}
    names = set()
    for specification in distributions:
        # Registry entries must be explicit public distribution records.
        if not isinstance(specification, RandomDistribution):
            raise TypeError("distributions must contain RandomDistribution records")
        if specification.id in sampled:
            raise ValueError(f"Duplicate distribution ID {specification.id}")
        if specification.name in names:
            raise ValueError(f"Duplicate distribution name {specification.name!r}")
        # Store one joint draw; every reference will read from this same tuple.
        sampled[specification.id] = _sample_distribution(
            specification, rng, stochastic
        )
        names.add(specification.name)
    return sampled


def _is_distribution_reference(value):
    """Return whether a value follows the random-reference tuple syntax.

    Parameters
    ----------
    value : object
        Candidate value from a declarative model record.

    Returns
    -------
    bool
        ``True`` for ``(distribution_id,)`` and
        ``(distribution_id, component)`` tuples.
    """
    # Boolean values are integers in Python, but are not valid IDs or indices.
    return (
        isinstance(value, tuple)
        and len(value) in (1, 2)
        and all(
            isinstance(item, int) and not isinstance(item, bool)
            for item in value
        )
    )


def _resolve_distribution_reference(reference, sampled_distributions, path):
    """Retrieve one sampled value through a distribution reference.

    Parameters
    ----------
    reference : tuple of int
        Scalar ``(distribution_id,)`` or component
        ``(distribution_id, component)`` reference.
    sampled_distributions : dict
        Sampled component tuples indexed by distribution ID.
    path : str
        Human-readable model location used in validation errors.

    Returns
    -------
    float
        Referenced scalar value for the current realization.

    Raises
    ------
    ValueError
        If the ID is unknown, a multivariate component is omitted, or the
        requested component is outside the sampled vector.
    """
    # The first tuple item always identifies the registered distribution.
    distribution_id = reference[0]
    if distribution_id not in sampled_distributions:
        raise ValueError(f"{path}: unknown distribution ID {distribution_id}")

    # Scalar references implicitly select the only available component.
    sampled_values = sampled_distributions[distribution_id]
    component = 0 if len(reference) == 1 else reference[1]

    # Multivariate variables require an explicit component to avoid ambiguity.
    if len(reference) == 1 and len(sampled_values) != 1:
        raise ValueError(
            f"{path}: multivariate distribution {distribution_id} "
            "requires a component index"
        )

    # Reject negative and oversized indices before indexing the sampled tuple.
    if component < 0 or component >= len(sampled_values):
        raise ValueError(
            f"{path}: component {component} is out of range for "
            f"distribution {distribution_id}"
        )

    return sampled_values[component]


def _resolve_dataclass(value, sampled_distributions, path):
    """Reconstruct a dataclass after recursively resolving all fields.

    Parameters
    ----------
    value : dataclass instance
        Declarative record to resolve.
    sampled_distributions : dict
        Sampled component tuples indexed by distribution ID.
    path : str
        Human-readable model location used in validation errors.

    Returns
    -------
    object
        Dataclass of the same type containing resolved field values.
    """
    # ``replace`` supports the frozen and slotted records exposed by the API.
    resolved_fields = {
        field.name: _resolve(
            getattr(value, field.name),
            sampled_distributions,
            f"{path}.{field.name}",
        )
        for field in fields(value)
    }
    return replace(value, **resolved_fields)


def _resolve(value, sampled_distributions, path="model"):
    """Recursively traverse model data and resolve random references.

    Parameters
    ----------
    value : object
        Dataclass, list, tuple, distribution reference, or scalar.
    sampled_distributions : dict
        Values indexed by registered distribution ID.
    path : str, optional
        Human-readable location used in validation errors.

    Returns
    -------
    object
        Value with the same structure and no distribution references.
    """
    # Distribution records define the registry and must remain unchanged.
    if isinstance(value, RandomDistribution):
        return value

    # A speed table owns two structural tuples. Resolve only its ordinates so
    # a two-integer speed grid such as ``(0, 100)`` cannot be mistaken for a
    # random-distribution reference by the generic tuple syntax.
    if isinstance(value, SpeedDependentValue):
        resolved_values = tuple(
            _resolve(
                item,
                sampled_distributions,
                f"{path}.values[{index}]",
            )
            for index, item in enumerate(value.values)
        )
        return replace(value, values=resolved_values)

    # Reference validation is isolated from the recursive traversal logic.
    if _is_distribution_reference(value):
        return _resolve_distribution_reference(
            value,
            sampled_distributions,
            path,
        )

    # Preserve structural tuples such as coordinate-system vectors.
    if isinstance(value, tuple):
        return tuple(
            _resolve(item, sampled_distributions, f"{path}[{index}]")
            for index, item in enumerate(value)
        )

    # Preserve lists while resolving every record or value they contain.
    if isinstance(value, list):
        return [
            _resolve(item, sampled_distributions, f"{path}[{index}]")
            for index, item in enumerate(value)
        ]

    # Delegate dataclass reconstruction to keep this dispatcher concise.
    if is_dataclass(value):
        return _resolve_dataclass(value, sampled_distributions, path)

    # Ordinary deterministic scalars require no conversion here.
    return value


def _unique(
    items,
    label,
    expected_type=None,
    *,
    validate_id=True,
    validated_id_types=(),
):
    """Index records by ID while validating type and uniqueness.

    Parameters
    ----------
    items : iterable
        Records containing an integer ``id`` field.
    label : str
        Record-family name used in error messages.
    expected_type : type or tuple of type, optional
        Accepted record class or classes.
    validate_id : bool, optional
        Validate IDs for record classes that do not yet own that invariant.
    validated_id_types : tuple of type, optional
        Record classes whose constructors already validate their IDs.
    validated_id_types : tuple of type, optional
        Record classes whose constructors already validate their IDs.

    Returns
    -------
    dict
        Insertion-ordered mapping from record ID to record.
    """
    # Build an ordinary dictionary, which preserves input order in Python.
    records_by_id = {}

    for item in items:
        # Reject records placed in the wrong ModelDefinition collection.
        if expected_type is not None and not isinstance(item, expected_type):
            if isinstance(expected_type, tuple):
                expected_names = ", ".join(
                    record_type.__name__ for record_type in expected_type
                )
            else:
                expected_names = expected_type.__name__
            raise TypeError(
                f"Every {label} record must be one of: {expected_names}"
            )

        # All indexed values must be dataclass records carrying an ID.
        if not is_dataclass(item) or not hasattr(item, "id"):
            raise TypeError(
                f"Every {label} record must be a dataclass with an ID"
            )

        if validate_id and not isinstance(item, validated_id_types):
            # Boolean values are excluded even though bool subclasses int.
            if isinstance(item.id, bool) or not isinstance(item.id, int):
                raise TypeError(f"{label.capitalize()} IDs must be integers")

            # ID zero is reserved for the implicit global coordinate system.
            if item.id < 0 or (item.id == 0 and label != "coordinate system"):
                raise ValueError(f"{label.capitalize()} IDs must be positive")

        # Duplicate IDs would make reference resolution ambiguous.
        if item.id in records_by_id:
            raise ValueError(f"Duplicate {label} ID {item.id}")

        records_by_id[item.id] = item

    return records_by_id


def _coordinate_system_matrix(coordinate_system):
    """Convert a coordinate-system record to an origin and rotation matrix.

    Parameters
    ----------
    coordinate_system : CoordinateSystem
        Frame defined by an origin, local x direction, and local xy-plane
        direction.

    Returns
    -------
    origin : numpy.ndarray
        Global origin with shape ``(3,)``.
    rotation : numpy.ndarray
        Matrix whose columns are local basis vectors in global coordinates.
    """
    # Convert the record's pure-Python vectors into numerical arrays.
    origin = np.asarray(coordinate_system.origin, dtype=float)
    x_direction = np.asarray(coordinate_system.x_axis, dtype=float)
    xy_direction = np.asarray(coordinate_system.xy_plane, dtype=float)

    # Normalize the first local basis direction.
    x_norm = np.linalg.norm(x_direction)
    local_x = x_direction / x_norm

    # The cross product constructs a right-handed local z direction.
    local_z = np.cross(local_x, xy_direction)
    z_norm = np.linalg.norm(local_z)
    local_z /= z_norm

    # Recompute local y so the final basis is exactly orthonormal.
    local_y = np.cross(local_z, local_x)

    # Columns map local vector components into the global coordinate system.
    rotation = np.column_stack((local_x, local_y, local_z))
    return origin, rotation


def _grid_coordinates(definition):
    """Resolve every grid coordinate into the global frame.

    Parameters
    ----------
    definition : ModelDefinition
        Resolved deterministic model definition.

    Returns
    -------
    numpy.ndarray
        Global grid coordinates with shape ``(grids, 3)``.
    """
    # Coordinate system zero is the implicit global Cartesian frame.
    systems = {0: (np.zeros(3), np.eye(3))}

    # Validate user-defined frames and reserve ID zero.
    coordinate_systems = _unique(
        definition.coordinate_systems,
        "coordinate system",
        CoordinateSystem,
        validate_id=False,
    )

    # Convert each user-defined record into numerical frame data.
    for coordinate_system in coordinate_systems.values():
        systems[coordinate_system.id] = _coordinate_system_matrix(
            coordinate_system
        )

    # Transform grid-local coordinates into global coordinates.
    global_coordinates = []
    for grid in definition.grids:
        # Every referenced input frame must exist.
        if grid.coordinate_system not in systems:
            raise ValueError(
                f"Grid {grid.id} references unknown coordinate system "
                f"{grid.coordinate_system}"
            )

        # Retrieve the selected frame and transform the validated coordinate.
        origin, rotation = systems[grid.coordinate_system]
        local_coordinate = np.array([grid.x, grid.y, grid.z], dtype=float)

        # Apply x_global = origin + R * x_local.
        global_coordinates.append(origin + rotation @ local_coordinate)

    return np.asarray(global_coordinates)


def _nodal_transform(rotation):
    """Build a six-DOF global-to-local nodal transformation.

    Parameters
    ----------
    rotation : numpy.ndarray
        ``(3, 3)`` local-to-global basis matrix.

    Returns
    -------
    numpy.ndarray
        ``(6, 6)`` transformation for translations and rotations.
    """
    # Rows of R.T convert global vector components into local components.
    transform = np.zeros((6, 6), dtype=float)
    transform[:3, :3] = rotation.T
    transform[3:, 3:] = rotation.T
    return transform


def _shaft_transform(point_a, point_b):
    """Build the length and transformation of a two-grid shaft.

    Parameters
    ----------
    point_a, point_b : array_like
        Global endpoint coordinates.

    Returns
    -------
    length : float
        Distance between endpoints.
    transform : numpy.ndarray
        ``(12, 12)`` global-to-local element transformation.
    """
    # The local z axis follows the element from grid A to grid B.
    local_z = np.asarray(point_b, dtype=float) - np.asarray(point_a, dtype=float)
    length = np.linalg.norm(local_z)
    if length <= 0:
        raise ValueError("A shaft element cannot have coincident grids")
    local_z /= length

    # Select a stable trial vector that is not almost parallel to local z.
    if abs(local_z[2]) < 0.9:
        trial_vector = np.array([0.0, 0.0, 1.0])
    else:
        trial_vector = np.array([1.0, 0.0, 0.0])

    # Complete a right-handed orthonormal local frame.
    local_x = np.cross(trial_vector, local_z)
    local_x /= np.linalg.norm(local_x)
    local_y = np.cross(local_z, local_x)

    # Rows convert global vector components into local components.
    rotation = np.vstack((local_x, local_y, local_z))

    # Apply the same rotation to both translations and rotations at both grids.
    transform = np.zeros((12, 12), dtype=float)
    for offset in (0, 3, 6, 9):
        transform[offset : offset + 3, offset : offset + 3] = rotation

    return length, transform


def _assemble(definition):
    """Assemble one resolved definition into global K0, M, C0, and G matrices.

    Parameters
    ----------
    definition : ModelDefinition
        Deterministic definition containing no distribution references.

    Returns
    -------
    stiffness_matrix, mass_matrix, damping_matrix, gyroscopic_matrix : numpy.ndarray
        Global ``(ndof, ndof)`` base matrices. Tabulated bearing terms are
        omitted from stiffness and damping and returned separately.
    grid_ids : tuple of int
        Grid IDs in assembly order.
    coordinates : numpy.ndarray
        Resolved global grid coordinates.
    stiffness_speed_contributions, damping_speed_contributions : tuple
        Transformed unit-coefficient matrix terms evaluated by solvers at
        requested rotor speeds.
    """
    # Validate and index every record collection.
    grids = _unique(definition.grids, "grid", Grid, validate_id=False)
    if not grids:
        raise ValueError("A model requires at least one grid")

    materials = _unique(
        definition.materials,
        "material",
        Material,
        validate_id=False,
    )
    properties = _unique(
        definition.properties,
        "property",
        (ShaftProperty, BearingProperty, DiskProperty),
        validated_id_types=(ShaftProperty, BearingProperty, DiskProperty),
    )
    _unique(
        definition.elements,
        "element",
        (ShaftElement, BearingElement, DiskElement),
        validated_id_types=(ShaftElement, BearingElement, DiskElement),
    )
    loads = _unique(
        definition.loads,
        "load",
        (HarmonicLoad, RotatingLoad, Unbalance),
        validated_id_types=(HarmonicLoad, RotatingLoad, Unbalance),
    )

    # Record insertion order determines global grid and DOF order.
    grid_ids = tuple(grids)
    grid_index = {grid_id: index for index, grid_id in enumerate(grid_ids)}

    # Resolve all grid positions before computing element geometry.
    coordinates = _grid_coordinates(definition)

    # Build numerical coordinate frames for nodal elements.
    coordinate_systems = {0: (np.zeros(3), np.eye(3))}
    for coordinate_system in definition.coordinate_systems:
        coordinate_systems[coordinate_system.id] = _coordinate_system_matrix(
            coordinate_system
        )

    # Normalize the model-wide physical spin direction.
    spin_axis = np.asarray(definition.spin_axis, dtype=float)
    spin_norm = np.linalg.norm(spin_axis)
    spin_axis /= spin_norm

    # Allocate dense global matrices for this realization.
    number_of_dofs = DOFS_PER_GRID * len(grid_ids)
    global_stiffness = np.zeros((number_of_dofs, number_of_dofs))
    global_mass = np.zeros((number_of_dofs, number_of_dofs))
    global_damping = np.zeros((number_of_dofs, number_of_dofs))
    global_gyroscopic = np.zeros((number_of_dofs, number_of_dofs))
    stiffness_speed_contributions = []
    damping_speed_contributions = []

    def add(global_matrix, local_matrix, indices):
        """Accumulate one local matrix into selected global indices.

        Parameters
        ----------
        global_matrix : numpy.ndarray
            Matrix being assembled in place.
        local_matrix : numpy.ndarray
            Element matrix to accumulate.
        indices : sequence of int
            Global DOF indices corresponding to local rows and columns.
        """
        global_matrix[np.ix_(indices, indices)] += local_matrix

    def integer_reference(value, label):
        """Validate one record reference as a pure integer ID.

        Parameters
        ----------
        value : object
            Reference value to validate.
        label : str
            Human-readable field name for an error message.

        Returns
        -------
        int
            Validated reference value.
        """
        if isinstance(value, bool) or not isinstance(value, int):
            raise TypeError(f"{label} must be an integer ID")
        return value

    # Validate declarative excitations without adding them to structural
    # matrices. Response solvers assemble sample-specific load vectors.
    for load in loads.values():
        if load.grid not in grid_index:
            raise ValueError(
                f"Load {load.id} references unknown grid {load.grid}"
            )
        if load.coordinate_system not in coordinate_systems:
            raise ValueError(
                f"Load {load.id} references unknown coordinate system "
                f"{load.coordinate_system}"
            )

        if isinstance(load, Unbalance):
            # The load frame supplies phase zero, but its transverse plane
            # must be the rotor plane. Parallel and antiparallel z axes are
            # physically equivalent and their sign can set the later phasor.
            load_axis = coordinate_systems[load.coordinate_system][1][:, 2]
            alignment = abs(float(load_axis @ spin_axis))
            if not np.isclose(alignment, 1.0, rtol=1.0e-10, atol=1.0e-12):
                raise ValueError(
                    f"Unbalance {load.id} coordinate-system z-axis must be "
                    "parallel to spin_axis"
                )

    # Assemble each declarative element record.
    for element in definition.elements:
        # Resolve the element's property reference first.
        if element.property not in properties:
            raise ValueError(
                f"Element {element.id} references unknown property "
                f"{element.property}"
            )
        element_property = properties[element.property]

        # Assemble a two-grid shaft element.
        if isinstance(element, ShaftElement):
            if not isinstance(element_property, ShaftProperty):
                raise TypeError(
                    f"Shaft element {element.id} requires a ShaftProperty"
                )

            # Validate and resolve shaft grid and material references.
            integer_reference(
                element_property.material,
                f"Shaft property {element_property.id} material",
            )
            try:
                index_a = grid_index[element.grid_a]
                index_b = grid_index[element.grid_b]
            except KeyError as error:
                raise ValueError(
                    f"Element {element.id} references unknown grid {error.args[0]}"
                ) from None

            if element_property.material not in materials:
                raise ValueError(
                    f"Shaft property {element_property.id} references unknown "
                    f"material {element_property.material}"
                )
            material = materials[element_property.material]

            # Derive element geometry and its global-to-local transformation.
            length, transform = _shaft_transform(
                coordinates[index_a], coordinates[index_b]
            )

            # Project the global spin vector onto the element's positive axis.
            local_spin = float(transform[2, :3] @ spin_axis)

            # Collect common kernel arguments in one readable mapping.
            shaft_arguments = {
                "length": length,
                "outer_diameter": element_property.outer_diameter,
                "inner_diameter": element_property.inner_diameter,
                "young_modulus": material.young_modulus,
                "poisson_ratio": material.poisson_ratio,
                "density": material.density,
                "theory": element_property.theory,
                "rotary_inertia": element_property.rotary_inertia,
            }

            # Stiffness does not use density or the inertia switch.
            stiffness_arguments = {
                key: value
                for key, value in shaft_arguments.items()
                if key not in ("density", "rotary_inertia")
            }

            # Evaluate all four local shaft matrices.
            local_stiffness = stiffness.shaft_stiffness(**stiffness_arguments)
            local_mass = mass.shaft_mass(**shaft_arguments)
            local_damping = damping.shaft_damping(
                damping=element_property.damping,
                mass_matrix=local_mass,
            )
            local_gyroscopic = (
                local_spin * gyroscopic.shaft_gyroscopic(**shaft_arguments)
            )

            # Concatenate the six global DOFs of both endpoint grids.
            shaft_indices = list(
                range(DOFS_PER_GRID * index_a, DOFS_PER_GRID * index_a + 6)
            ) + list(
                range(DOFS_PER_GRID * index_b, DOFS_PER_GRID * index_b + 6)
            )

            # Rotate every local matrix and add it to its global positions.
            matrix_pairs = (
                (global_stiffness, local_stiffness),
                (global_mass, local_mass),
                (global_damping, local_damping),
                (global_gyroscopic, local_gyroscopic),
            )
            for global_matrix, local_matrix in matrix_pairs:
                rotated_matrix = transform.T @ local_matrix @ transform
                add(global_matrix, rotated_matrix, shaft_indices)

        # Assemble a bearing connected to ground or to a second grid.
        elif isinstance(element, BearingElement):
            if not isinstance(element_property, BearingProperty):
                raise TypeError(
                    f"Bearing element {element.id} requires a BearingProperty"
                )

            if element.grid not in grid_index:
                raise ValueError(
                    f"Element {element.id} references unknown grid {element.grid}"
                )
            if element.coordinate_system not in coordinate_systems:
                raise ValueError(
                    f"Element {element.id} references unknown coordinate system "
                    f"{element.coordinate_system}"
                )

            # The base matrices contain only constant coefficients. Each
            # tabulated coefficient is replaced by zero and represented by a
            # transformed unit-coefficient basis stored below.
            stiffness_names = ("kxx", "kyy", "kzz", "kxy", "kyx")
            damping_names = ("cxx", "cyy", "czz", "cxy", "cyx")
            stiffness_arguments = {
                name: (
                    0.0
                    if isinstance(getattr(element_property, name), SpeedDependentValue)
                    else getattr(element_property, name)
                )
                for name in stiffness_names
            }
            damping_arguments = {
                name: (
                    0.0
                    if isinstance(getattr(element_property, name), SpeedDependentValue)
                    else getattr(element_property, name)
                )
                for name in damping_names
            }
            local_stiffness = stiffness.bearing_stiffness(**stiffness_arguments)
            local_damping = damping.bearing_damping(**damping_arguments)
            local_mass = mass.bearing_mass()
            local_gyroscopic = gyroscopic.bearing_gyroscopic()

            # Rotate bearing coefficients from their input frame to global DOFs.
            rotation = coordinate_systems[element.coordinate_system][1]
            transform = _nodal_transform(rotation)
            bearing_matrices = [
                transform.T @ local_stiffness @ transform,
                transform.T @ local_mass @ transform,
                transform.T @ local_damping @ transform,
                transform.T @ local_gyroscopic @ transform,
            ]
            speed_terms = []
            for name in stiffness_names:
                value = getattr(element_property, name)
                if isinstance(value, SpeedDependentValue):
                    unit_local = stiffness.bearing_stiffness(**{name: 1.0})
                    speed_terms.append(
                        (
                            "stiffness",
                            value,
                            transform.T @ unit_local @ transform,
                            f"Bearing property {element_property.id} {name}",
                        )
                    )
            for name in damping_names:
                value = getattr(element_property, name)
                if isinstance(value, SpeedDependentValue):
                    unit_local = damping.bearing_damping(**{name: 1.0})
                    speed_terms.append(
                        (
                            "damping",
                            value,
                            transform.T @ unit_local @ transform,
                            f"Bearing property {element_property.id} {name}",
                        )
                    )

            # Start with the six DOFs of the first grid.
            index_a = grid_index[element.grid]
            bearing_indices = list(
                range(DOFS_PER_GRID * index_a, DOFS_PER_GRID * index_a + 6)
            )

            # A second grid converts support-to-ground into relative coupling.
            if element.grid_b is not None:
                if element.grid_b not in grid_index:
                    raise ValueError(
                        f"Element {element.id} references unknown grid "
                        f"{element.grid_b}"
                    )
                index_b = grid_index[element.grid_b]
                bearing_indices += list(
                    range(
                        DOFS_PER_GRID * index_b,
                        DOFS_PER_GRID * index_b + 6,
                    )
                )

                # q_relative = q_a - q_b produces [[B, -B], [-B, B]].
                relative_transform = np.block([[np.eye(6), -np.eye(6)]])
                bearing_matrices = [
                    relative_transform.T @ matrix @ relative_transform
                    for matrix in bearing_matrices
                ]
                speed_terms = [
                    (
                        target,
                        value,
                        relative_transform.T @ matrix @ relative_transform,
                        path,
                    )
                    for target, value, matrix, path in speed_terms
                ]

            # Add bearing K, M, C, and G in the same order as their globals.
            global_matrices = (
                global_stiffness,
                global_mass,
                global_damping,
                global_gyroscopic,
            )
            for global_matrix, local_matrix in zip(
                global_matrices, bearing_matrices, strict=True
            ):
                add(global_matrix, local_matrix, bearing_indices)

            # Store compact transformed local bases. Solvers expand them into
            # global delta matrices only for the requested operating speeds.
            for target, value, matrix, path in speed_terms:
                contribution = _SpeedDependentMatrixContribution(
                    value=value,
                    matrix=matrix,
                    indices=tuple(bearing_indices),
                    path=path,
                )
                if target == "stiffness":
                    stiffness_speed_contributions.append(contribution)
                else:
                    damping_speed_contributions.append(contribution)

        # Assemble a one-grid rigid disk or lumped mass.
        elif isinstance(element, DiskElement):
            if not isinstance(element_property, DiskProperty):
                raise TypeError(
                    f"Disk element {element.id} requires a DiskProperty"
                )

            if element.grid not in grid_index:
                raise ValueError(
                    f"Element {element.id} references unknown grid {element.grid}"
                )
            if element.coordinate_system not in coordinate_systems:
                raise ValueError(
                    f"Element {element.id} references unknown coordinate system "
                    f"{element.coordinate_system}"
                )

            # Build the disk's local orientation and spin projection.
            rotation = coordinate_systems[element.coordinate_system][1]
            transform = _nodal_transform(rotation)
            local_spin = float(rotation[:, 2] @ spin_axis)

            # Evaluate the local disk mass once because damping reuses it.
            local_mass = mass.disk_mass(
                mass=element_property.mass,
                diametral_inertia=element_property.diametral_inertia,
                polar_inertia=element_property.polar_inertia,
            )

            # Pair each global matrix with its local disk contribution.
            disk_matrix_pairs = (
                (global_stiffness, stiffness.disk_stiffness()),
                (global_mass, local_mass),
                (
                    global_damping,
                    damping.disk_damping(
                        damping=element_property.damping,
                        mass_matrix=local_mass,
                    ),
                ),
                (
                    global_gyroscopic,
                    local_spin
                    * gyroscopic.disk_gyroscopic(
                        polar_inertia=element_property.polar_inertia
                    ),
                ),
            )

            # The disk acts on the six DOFs of its referenced grid.
            grid_position = grid_index[element.grid]
            disk_indices = list(
                range(
                    DOFS_PER_GRID * grid_position,
                    DOFS_PER_GRID * grid_position + 6,
                )
            )

            # Rotate and assemble every local disk matrix.
            for global_matrix, local_matrix in disk_matrix_pairs:
                rotated_matrix = transform.T @ local_matrix @ transform
                add(global_matrix, rotated_matrix, disk_indices)

        # This branch protects assembly if a new record bypasses validation.
        else:
            raise TypeError(
                f"Unsupported element record {type(element).__name__}"
            )

    return (
        global_stiffness,
        global_mass,
        global_damping,
        global_gyroscopic,
        grid_ids,
        coordinates,
        tuple(stiffness_speed_contributions),
        tuple(damping_speed_contributions),
    )


class ModelBuilder:
    """Build deterministic numerical models or Monte Carlo ensembles."""

    def build(
        self,
        definition,
        *,
        stochastic=False,
        samples=1,
        seed=None,
    ):
        """Resolve and assemble a declarative model definition.

        Parameters
        ----------
        definition : ModelDefinition
            Model records to validate, sample, and assemble.
        stochastic : bool, optional
            Draw random parameters and build a Monte Carlo ensemble. The
            default uses distribution means and builds one realization.
        samples : int, optional
            Number of Monte Carlo realizations when ``stochastic=True``.
        seed : int or None, optional
            Seed passed to NumPy's random generator.

        Returns
        -------
        BuiltModel
            Global matrices and resolved definitions for all realizations.
        """
        # Only the declarative root record can be assembled.
        if not isinstance(definition, ModelDefinition):
            raise TypeError("definition must be a ModelDefinition")

        # Avoid ambiguous truthy values such as the string ``"false"``.
        if not isinstance(stochastic, bool):
            raise TypeError("stochastic must be a boolean")

        # Stochastic construction requires an explicit positive sample count.
        if stochastic and (
            isinstance(samples, bool)
            or not isinstance(samples, int)
            or samples < 1
        ):
            raise ValueError("samples must be a positive integer")
        if not stochastic:
            samples = 1

        # Placement is independent of random realizations, so validate it
        # once while preserving the contextual diagnostics of the guarded
        # per-sample build path below.
        try:
            _validate_speed_dependent_placements(definition)
        except (TypeError, ValueError) as error:
            sample_label = "sample 0" if stochastic else "deterministic model"
            raise type(error)(
                f"Failed to build {sample_label}: {error}"
            ) from error

        # One generator supplies reproducible scalar and multivariate draws.
        random_generator = np.random.default_rng(seed)

        # Accumulate resolved records, matrices, and coordinates by sample.
        resolved_definitions = []
        assembled_matrices = []
        sample_coordinates = []
        sample_stiffness_speed_contributions = []
        sample_damping_speed_contributions = []

        for sample_index in range(samples):
            try:
                # Generate each registered distribution once for this sample.
                sampled_distributions = _sample_distributions(
                    definition.distributions,
                    random_generator,
                    stochastic,
                )

                # Replace every tuple reference with its generated component.
                resolved_definition = _resolve(
                    definition,
                    sampled_distributions,
                )

                # Assemble one complete deterministic numerical model.
                (
                    stiffness_matrix,
                    mass_matrix,
                    damping_matrix,
                    gyroscopic_matrix,
                    grid_ids,
                    coordinates,
                    stiffness_speed_contributions,
                    damping_speed_contributions,
                ) = _assemble(resolved_definition)
            except (TypeError, ValueError) as error:
                # Attach the failing realization without hiding the root cause.
                sample_label = (
                    f"sample {sample_index}"
                    if stochastic
                    else "deterministic model"
                )
                raise type(error)(
                    f"Failed to build {sample_label}: {error}"
                ) from error

            # Preserve this realization and its four global matrices.
            resolved_definitions.append(resolved_definition)
            assembled_matrices.append(
                (
                    stiffness_matrix,
                    mass_matrix,
                    damping_matrix,
                    gyroscopic_matrix,
                )
            )
            sample_coordinates.append(coordinates)
            sample_stiffness_speed_contributions.append(
                stiffness_speed_contributions
            )
            sample_damping_speed_contributions.append(
                damping_speed_contributions
            )

        # Stack each matrix family along the leading sample dimension.
        stiffness_matrices = np.stack(
            [matrices[0] for matrices in assembled_matrices]
        )
        mass_matrices = np.stack(
            [matrices[1] for matrices in assembled_matrices]
        )
        damping_matrices = np.stack(
            [matrices[2] for matrices in assembled_matrices]
        )
        gyroscopic_matrices = np.stack(
            [matrices[3] for matrices in assembled_matrices]
        )

        # Return one immutable numerical container for solver consumption.
        return BuiltModel(
            stiffness=stiffness_matrices,
            mass=mass_matrices,
            damping=damping_matrices,
            gyroscopic=gyroscopic_matrices,
            grid_ids=grid_ids,
            coordinates=np.stack(sample_coordinates),
            definitions=tuple(resolved_definitions),
            stochastic=stochastic,
            seed=seed,
            _stiffness_speed_contributions=tuple(
                sample_stiffness_speed_contributions
            ),
            _damping_speed_contributions=tuple(
                sample_damping_speed_contributions
            ),
        )


def build_model(definition, **options):
    """Build a model using a temporary ``ModelBuilder`` instance.

    Parameters
    ----------
    definition : ModelDefinition
        Declarative input records.
    **options
        Keyword arguments forwarded to :meth:`ModelBuilder.build`.

    Returns
    -------
    BuiltModel
        Assembled deterministic model or stochastic ensemble.
    """
    # This convenience function keeps simple scripts concise.
    return ModelBuilder().build(definition, **options)
