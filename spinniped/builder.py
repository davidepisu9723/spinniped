"""Model validation, random sampling, coordinate handling, and assembly."""

from dataclasses import dataclass, fields, is_dataclass, replace

import numpy as np

from . import damping, gyroscopic, mass, stiffness
from .records import (
    BearingElement,
    BearingProperty,
    CoordinateSystem,
    DiskElement,
    DiskProperty,
    Grid,
    Material,
    ModelDefinition,
    RandomDistribution,
    ShaftElement,
    ShaftProperty,
)


# Every grid owns three translations and three rotations.
DOFS_PER_GRID = 6


@dataclass(frozen=True, slots=True)
class BuiltModel:
    """Numerical model or Monte Carlo ensemble produced by ``ModelBuilder``.

    Parameters
    ----------
    stiffness, mass, damping, gyroscopic : numpy.ndarray
        Global matrices with shape ``(samples, ndof, ndof)``.
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
        """Return the global stiffness matrices.

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
        """Return the global damping matrices.

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


def _finite_array(value, label, *, dimensions):
    """Convert distribution parameters to a finite floating-point array.

    Parameters
    ----------
    value : object
        Scalar or sequence to convert.
    label : str
        Parameter name used in validation messages.
    dimensions : tuple of int
        Permitted numbers of array dimensions.

    Returns
    -------
    numpy.ndarray
        Validated floating-point array.

    Raises
    ------
    TypeError
        If the parameter cannot be converted to floating point.
    ValueError
        If its shape or contents are invalid.
    """
    # Normalize scalars and sequences before checking their common properties.
    try:
        array = np.asarray(value, dtype=float)
    except (TypeError, ValueError) as error:
        raise TypeError(f"{label} must contain only numeric values") from error
    if array.ndim not in dimensions:
        raise ValueError(f"{label} has an invalid shape")
    if not np.isfinite(array).all():
        raise ValueError(f"{label} must contain only finite values")
    return array


def _sample_distribution(specification, rng, stochastic):
    """Generate one scalar or vector realization from a distribution.

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

    Raises
    ------
    TypeError
        If the distribution record or parameters have invalid types.
    ValueError
        If the distribution definition is incomplete or inconsistent.
    """
    if not isinstance(specification, RandomDistribution):
        raise TypeError("distributions must contain RandomDistribution records")
    if not isinstance(specification.name, str) or not specification.name.strip():
        raise ValueError("Distribution names must be non-empty strings")
    if not isinstance(specification.distribution, str):
        raise TypeError("Distribution type must be a string")
    if not isinstance(specification.parameters, dict):
        raise TypeError("Distribution parameters must be a dictionary")

    # Family names are case-insensitive, while parameter names remain strict.
    distribution_name = specification.distribution.lower()
    parameters = specification.parameters

    # A scalar normal distribution produces exactly one component.
    if distribution_name == "normal":
        if set(parameters) != {"mean", "stdv"}:
            raise ValueError("A normal distribution requires only 'mean' and 'stdv'")
        mean = _finite_array(parameters["mean"], "mean", dimensions=(0,))
        deviation = _finite_array(parameters["stdv"], "stdv", dimensions=(0,))
        if deviation < 0:
            raise ValueError("Distribution standard deviation cannot be negative")
        draw = rng.normal(mean, deviation) if stochastic else mean

    # A scalar uniform distribution is represented by its two bounds.
    elif distribution_name == "uniform":
        if set(parameters) != {"low", "high"}:
            raise ValueError("A uniform distribution requires only 'low' and 'high'")
        lower = _finite_array(parameters["low"], "low", dimensions=(0,))
        upper = _finite_array(parameters["high"], "high", dimensions=(0,))
        if upper < lower:
            raise ValueError("Uniform distribution requires high >= low")
        draw = (
            rng.uniform(lower, upper)
            if stochastic
            else (lower + upper) / 2.0
        )

    # A multivariate normal jointly generates all correlated components.
    elif distribution_name == "multivariate_normal":
        required = {"mean", "stdv", "correlation"}
        if set(parameters) != required:
            raise ValueError(
                "A multivariate normal distribution requires only "
                "'mean', 'stdv', and 'correlation'"
            )
        mean = _finite_array(parameters["mean"], "mean", dimensions=(1,))
        deviation = _finite_array(parameters["stdv"], "stdv", dimensions=(1,))
        correlation = _finite_array(
            parameters["correlation"], "correlation", dimensions=(2,)
        )
        component_count = mean.size
        if component_count == 0:
            raise ValueError("A multivariate distribution needs at least one component")
        if deviation.shape != mean.shape:
            raise ValueError("mean and stdv must have the same length")
        if correlation.shape != (component_count, component_count):
            raise ValueError("correlation shape must match the number of components")
        if np.any(deviation < 0):
            raise ValueError("Distribution standard deviations cannot be negative")
        if not np.allclose(correlation, correlation.T):
            raise ValueError("correlation must be symmetric")
        if not np.allclose(np.diag(correlation), 1.0):
            raise ValueError("correlation diagonal entries must equal one")
        if np.any(np.abs(correlation) > 1.0):
            raise ValueError("correlation coefficients must lie between -1 and 1")
        if np.linalg.eigvalsh(correlation).min() < -1e-12:
            raise ValueError("correlation must be positive semidefinite")
        # Convert standard deviations and correlation into covariance units.
        covariance = np.outer(deviation, deviation) * correlation
        draw = rng.multivariate_normal(mean, covariance) if stochastic else mean

    else:
        raise ValueError(f"Unsupported distribution {distribution_name!r}")

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
        if isinstance(specification.id, bool) or not isinstance(specification.id, int):
            raise TypeError("Distribution IDs must be integers")
        if specification.id <= 0:
            raise ValueError("Distribution IDs must be positive")
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


def _resolve_reference(reference, sampled_distributions, path):
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

    # Reference validation is isolated from the recursive traversal logic.
    if _is_distribution_reference(value):
        return _resolve_reference(
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


def _unique(items, label, expected_type=None):
    """Index records by ID while validating type and uniqueness.

    Parameters
    ----------
    items : iterable
        Records containing an integer ``id`` field.
    label : str
        Record-family name used in error messages.
    expected_type : type or tuple of type, optional
        Accepted record class or classes.

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

        # Boolean values are intentionally excluded even though bool subclasses int.
        if isinstance(item.id, bool) or not isinstance(item.id, int):
            raise TypeError(f"{label.capitalize()} IDs must be integers")

        # ID zero is reserved only for the implicit global coordinate system.
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

    # Every coordinate-system vector is three-dimensional.
    if (
        origin.shape != (3,)
        or x_direction.shape != (3,)
        or xy_direction.shape != (3,)
    ):
        raise ValueError(
            f"Coordinate system {coordinate_system.id} vectors must have "
            "three entries"
        )

    # Reject invalid numerical frame definitions.
    all_values = np.concatenate((origin, x_direction, xy_direction))
    if not np.isfinite(all_values).all():
        raise ValueError(
            f"Coordinate system {coordinate_system.id} must contain finite values"
        )

    # Normalize the first local basis direction.
    x_norm = np.linalg.norm(x_direction)
    if x_norm < 1.0e-14:
        raise ValueError(
            f"Coordinate system {coordinate_system.id} has a zero-length x axis"
        )
    local_x = x_direction / x_norm

    # The cross product constructs a right-handed local z direction.
    local_z = np.cross(local_x, xy_direction)
    z_norm = np.linalg.norm(local_z)
    if z_norm < 1.0e-14:
        raise ValueError(
            f"Coordinate system {coordinate_system.id} has collinear axes"
        )
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
    )
    if 0 in coordinate_systems:
        raise ValueError(
            "Coordinate system ID 0 is reserved for the global system"
        )

    # Convert each user-defined record into numerical frame data.
    for coordinate_system in coordinate_systems.values():
        systems[coordinate_system.id] = _coordinate_system_matrix(
            coordinate_system
        )

    # Transform grid-local coordinates into global coordinates.
    global_coordinates = []
    for grid in definition.grids:
        # Coordinate-system references are integer IDs.
        if (
            isinstance(grid.coordinate_system, bool)
            or not isinstance(grid.coordinate_system, int)
        ):
            raise TypeError(
                f"Grid {grid.id} coordinate_system must be an integer ID"
            )

        # Every referenced input frame must exist.
        if grid.coordinate_system not in systems:
            raise ValueError(
                f"Grid {grid.id} references unknown coordinate system "
                f"{grid.coordinate_system}"
            )

        # Retrieve the selected frame and validate the local coordinate.
        origin, rotation = systems[grid.coordinate_system]
        local_coordinate = np.array([grid.x, grid.y, grid.z], dtype=float)
        if not np.isfinite(local_coordinate).all():
            raise ValueError(f"Grid {grid.id} coordinates must be finite")

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
    """Assemble one resolved definition into global K, M, C, and G matrices.

    Parameters
    ----------
    definition : ModelDefinition
        Deterministic definition containing no distribution references.

    Returns
    -------
    stiffness_matrix, mass_matrix, damping_matrix, gyroscopic_matrix : numpy.ndarray
        Global ``(ndof, ndof)`` matrices.
    grid_ids : tuple of int
        Grid IDs in assembly order.
    coordinates : numpy.ndarray
        Resolved global grid coordinates.
    """
    # Validate and index every record collection.
    grids = _unique(definition.grids, "grid", Grid)
    if not grids:
        raise ValueError("A model requires at least one grid")

    materials = _unique(definition.materials, "material", Material)
    properties = _unique(
        definition.properties,
        "property",
        (ShaftProperty, BearingProperty, DiskProperty),
    )
    _unique(
        definition.elements,
        "element",
        (ShaftElement, BearingElement, DiskElement),
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
    if spin_axis.shape != (3,) or not np.isfinite(spin_axis).all():
        raise ValueError("spin_axis must contain three finite values")
    spin_norm = np.linalg.norm(spin_axis)
    if spin_norm < 1.0e-14:
        raise ValueError("spin_axis cannot be zero")
    spin_axis /= spin_norm

    # Allocate dense global matrices for this realization.
    number_of_dofs = DOFS_PER_GRID * len(grid_ids)
    global_stiffness = np.zeros((number_of_dofs, number_of_dofs))
    global_mass = np.zeros((number_of_dofs, number_of_dofs))
    global_damping = np.zeros((number_of_dofs, number_of_dofs))
    global_gyroscopic = np.zeros((number_of_dofs, number_of_dofs))

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

    # Assemble each declarative element record.
    for element in definition.elements:
        # Resolve the element's property reference first.
        integer_reference(element.property, f"Element {element.id} property")
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
            integer_reference(element.grid_a, f"Element {element.id} grid_a")
            integer_reference(element.grid_b, f"Element {element.id} grid_b")
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

            # Validate all bearing references.
            integer_reference(element.grid, f"Element {element.id} grid")
            integer_reference(
                element.coordinate_system,
                f"Element {element.id} coordinate_system",
            )
            if element.grid_b is not None:
                integer_reference(element.grid_b, f"Element {element.id} grid_b")

            if element.grid not in grid_index:
                raise ValueError(
                    f"Element {element.id} references unknown grid {element.grid}"
                )
            if element.coordinate_system not in coordinate_systems:
                raise ValueError(
                    f"Element {element.id} references unknown coordinate system "
                    f"{element.coordinate_system}"
                )

            # Build the local bearing matrices from property coefficients.
            local_stiffness = stiffness.bearing_stiffness(
                kxx=element_property.kxx,
                kyy=element_property.kyy,
                kzz=element_property.kzz,
                kxy=element_property.kxy,
                kyx=element_property.kyx,
            )
            local_damping = damping.bearing_damping(
                cxx=element_property.cxx,
                cyy=element_property.cyy,
                czz=element_property.czz,
                cxy=element_property.cxy,
                cyx=element_property.cyx,
            )
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
                if element.grid_b == element.grid:
                    raise ValueError(
                        f"Bearing element {element.id} must connect distinct grids"
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

        # Assemble a one-grid rigid disk or lumped mass.
        elif isinstance(element, DiskElement):
            if not isinstance(element_property, DiskProperty):
                raise TypeError(
                    f"Disk element {element.id} requires a DiskProperty"
                )

            # Validate disk grid and orientation references.
            integer_reference(element.grid, f"Element {element.id} grid")
            integer_reference(
                element.coordinate_system,
                f"Element {element.id} coordinate_system",
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

        # One generator supplies reproducible scalar and multivariate draws.
        random_generator = np.random.default_rng(seed)

        # Accumulate resolved records, matrices, and coordinates by sample.
        resolved_definitions = []
        assembled_matrices = []
        sample_coordinates = []

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
