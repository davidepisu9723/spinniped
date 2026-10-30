"""Unified analyses for deterministic models and stochastic ensembles."""

import warnings

import numpy as np

from .builder import (
    BuiltModel,
    _coordinate_system_matrix,
    _interpolate_speed_dependent_value,
)
from .records import (
    HarmonicLoad,
    RotatingLoad,
    SpeedDependentValue,
    SynchronousLoad,
    Unbalance,
)


class Solver:
    """Run every supported analysis through a single ``solve`` entry point.

    Parameters
    ----------
    model : BuiltModel or None, optional
        Numerical model stored for subsequent calls to :meth:`solve`.

    Notes
    -----
    Every analysis loops over the model's leading sample dimension. Therefore,
    deterministic and stochastic models follow the same numerical path.
    """

    # These names form the stable public analysis vocabulary.
    _ANALYSES = (
        "modal",
        "campbell",
        "frequency_response",
        "unbalance_response",
        "time_response",
    )

    def __init__(self, model: BuiltModel | None = None):
        """Initialize a solver with an optional built model.

        Parameters
        ----------
        model : BuiltModel or None, optional
            Model used when :meth:`solve` does not receive an override.
        """
        # Store the model without performing an analysis during construction.
        self.model = model

    def solve(self, analysis="modal", *, model=None, fixed_dofs=None, **options):
        """Run one selected analysis for every model realization.

        Parameters
        ----------
        analysis : str, optional
            One of ``modal``, ``campbell``, ``frequency_response``,
            ``unbalance_response``, or ``time_response``.
        model : BuiltModel or None, optional
            Per-call model override. By default, use the constructor model.
        fixed_dofs : sequence of int or None, optional
            Global degrees of freedom removed before solution.
        **options
            Analysis-specific keyword arguments.

        Returns
        -------
        dict
            Analysis arrays plus common metadata.

        Raises
        ------
        TypeError
            If no valid built model is available.
        ValueError
            If the analysis name or boundary conditions are invalid.
        """
        # A per-call model takes precedence over the constructor value.
        numerical_model = self.model if model is None else model
        if not isinstance(numerical_model, BuiltModel):
            raise TypeError("model must be a model returned by ModelBuilder")

        # Convert fixed global DOFs into a sorted free-DOF array.
        free_dofs = self._free_dofs(numerical_model.ndof, fixed_dofs)

        # Normalize the public name before selecting an analysis strategy.
        if not isinstance(analysis, str):
            raise TypeError("analysis must be a string")

        normalized_analysis = analysis.lower()
        
        # Select the requested strategy explicitly so the dispatch remains
        # visible when reading this public entry point.
        if normalized_analysis == "modal":
            result = self._solve_modal(numerical_model, free_dofs, **options)

        elif normalized_analysis == "campbell":
            result = self._solve_campbell(numerical_model, free_dofs, **options)

        elif normalized_analysis == "frequency_response":
            result = self._solve_frequency_response(
                numerical_model,
                free_dofs,
                **options,
            )

        elif normalized_analysis == "unbalance_response":
            result = self._solve_unbalance_response(
                numerical_model,
                free_dofs,
                **options,
            )

        elif normalized_analysis == "time_response":
            result = self._solve_time_response(
                numerical_model,
                free_dofs,
                **options,
            )

        else:
            raise ValueError(
                f"Unknown analysis {analysis!r}; choose from {self._ANALYSES}"
            )

        # Add metadata common to all analysis result dictionaries.
        result.update(
            analysis=normalized_analysis,
            free_dofs=free_dofs,
            samples=numerical_model.samples,
            stochastic=numerical_model.stochastic,
        )
        return result

    @staticmethod
    def _free_dofs(number_of_dofs, fixed_dofs):
        """Validate constraints and return complementary free DOFs.

        Parameters
        ----------
        number_of_dofs : int
            Total number of global degrees of freedom.
        fixed_dofs : sequence of int or None
            Global indices to remove.

        Returns
        -------
        numpy.ndarray
            Sorted unique free-DOF indices.
        """
        # Treat a missing constraint sequence as an empty sequence.
        if fixed_dofs is None:
            fixed_values = []
        else:
            try:
                fixed_values = list(fixed_dofs)
            except TypeError:
                raise TypeError(
                    "fixed_dofs must be a sequence of integers"
                ) from None

        # Reject floats and Booleans instead of silently converting them.
        invalid_type = any(
            isinstance(value, (bool, np.bool_))
            or not isinstance(value, (int, np.integer))
            for value in fixed_values
        )
        if invalid_type:
            raise TypeError("fixed_dofs must contain only integers")

        # Convert validated values to an integer NumPy array.
        fixed_array = np.asarray(fixed_values, dtype=int)
        if fixed_array.ndim != 1:
            raise ValueError("fixed_dofs must be a one-dimensional sequence")

        # All indices must address existing global DOFs.
        if fixed_array.size:
            if fixed_array.min() < 0 or fixed_array.max() >= number_of_dofs:
                raise ValueError("fixed_dofs contains an out-of-range index")

        # Remove duplicates before computing the set difference.
        unique_fixed = np.unique(fixed_array)
        all_dofs = np.arange(number_of_dofs)
        free_dofs = np.setdiff1d(all_dofs, unique_fixed)

        # A fully constrained model has no unknowns to solve.
        if not free_dofs.size:
            raise ValueError("At least one degree of freedom must remain free")

        return free_dofs

    @staticmethod
    def _reduced(matrix, sample, free_dofs):
        """Extract one sample's free-free matrix block.

        Parameters
        ----------
        matrix : numpy.ndarray
            Matrix ensemble with shape ``(samples, ndof, ndof)``.
        sample : int
            Realization index.
        free_dofs : numpy.ndarray
            Retained global DOF indices.

        Returns
        -------
        numpy.ndarray
            Reduced square matrix.
        """
        # np.ix_ forms the Cartesian product of selected rows and columns.
        return matrix[sample][np.ix_(free_dofs, free_dofs)]

    @staticmethod
    def _scalar_speed(value):
        """Return one finite signed rotor angular speed."""
        if isinstance(value, (bool, np.bool_)):
            raise TypeError("speed must be a real scalar")
        try:
            speed = float(value)
        except (TypeError, ValueError):
            raise TypeError("speed must be a real scalar") from None
        if not np.isfinite(speed):
            raise ValueError("speed must be finite")
        return speed

    @staticmethod
    def _speed_matrix_deltas(model, contributions, speeds):
        """Resolve all tabulated coefficients before constructing deltas.

        The returned transient array has shape
        ``(samples, speeds, ndof, ndof)``. It belongs to one analysis call and
        never changes the permanent three-dimensional BuiltModel matrices.
        """
        if not contributions or not any(contributions):
            return None

        # First interpolate every table across the complete requested grid.
        # No eigensolve or linear solve can begin until this pass succeeds.
        resolved = []
        for sample, sample_terms in enumerate(contributions):
            for term in sample_terms:
                coefficients = _interpolate_speed_dependent_value(
                    term.value,
                    speeds,
                    term.path,
                )
                resolved.append((sample, term, coefficients))

        deltas = np.zeros(
            (model.samples, len(speeds), model.ndof, model.ndof),
            dtype=float,
        )
        for sample, term, coefficients in resolved:
            matrix_indices = np.ix_(term.indices, term.indices)
            for speed_index, coefficient in enumerate(coefficients):
                deltas[sample, speed_index][matrix_indices] += (
                    coefficient * term.matrix
                )
        return deltas

    @classmethod
    def _reduced_at_speed(
        cls,
        base,
        deltas,
        sample,
        speed_index,
        free_dofs,
    ):
        """Return one reduced base matrix plus its optional speed delta."""
        matrix = cls._reduced(base, sample, free_dofs)
        if deltas is not None:
            matrix = matrix + deltas[sample, speed_index][
                np.ix_(free_dofs, free_dofs)
            ]
        return matrix

    def _solve_modal(
        self,
        model,
        free_dofs,
        modes=None,
        zero_tolerance=None,
        track_modes=True,
        speed=0.0,
    ):
        """Solve real symmetric generalized eigenproblems.

        Parameters
        ----------
        model : BuiltModel
            Numerical model or stochastic ensemble.
        free_dofs : numpy.ndarray
            Retained global DOF indices.
        modes : int or None, optional
            Number of positive modes to return.
        zero_tolerance : float or None, optional
            Absolute threshold used to remove rigid-body eigenvalues.
        track_modes : bool, optional
            Match stochastic modes to the first sample using the MAC.
        speed : float, optional
            Rotor angular speed in radians per second used to evaluate
            speed-dependent stiffness. Modal analysis still omits damping and
            gyroscopic terms.

        Returns
        -------
        dict
            Eigenvalues, frequencies in hertz, and reduced eigenvectors.
        """
        # Import SciPy lazily so model-definition workflows stay lightweight.
        from scipy.linalg import eigh

        # Validate options shared by every realization.
        requested_modes = self._validate_modes(modes)
        if not isinstance(track_modes, bool):
            raise TypeError("track_modes must be a boolean")
        speed = self._scalar_speed(speed)

        # Resolve every sample's stiffness at the operating point before the
        # first generalized eigensolution is attempted.
        stiffness_deltas = self._speed_matrix_deltas(
            model,
            model._stiffness_speed_contributions,
            np.array([speed]),
        )

        # Store one eigensolution per realization.
        sample_eigenvalues = []
        sample_eigenvectors = []

        for sample in range(model.samples):

            print(f"Analyzing Sample {sample}")

            # Extract the free-free stiffness and mass matrices.
            stiffness_matrix = self._reduced_at_speed(
                model.K,
                stiffness_deltas,
                sample,
                0,
                free_dofs,
            )
            mass_matrix = self._reduced(model.M, sample, free_dofs)

            # The Hermitian solver is valid only for symmetric matrices.
            if not np.allclose(
                stiffness_matrix,
                stiffness_matrix.T,
                rtol=1.0e-10,
                atol=1.0e-12,
            ):
                raise ValueError(
                    "modal analysis requires symmetric stiffness; use "
                    "campbell for cross-coupled bearings"
                )
            if not np.allclose(
                mass_matrix,
                mass_matrix.T,
                rtol=1.0e-10,
                atol=1.0e-12,
            ):
                raise ValueError("modal analysis requires symmetric mass")

            # Solve K q = lambda M q with mass-normalized eigenvectors.
            eigenvalues, eigenvectors = eigh(stiffness_matrix, mass_matrix)

            # Derive a round-off-aware default rigid-mode threshold.
            if zero_tolerance is None:
                eigenvalue_scale = max(
                    float(np.max(np.abs(eigenvalues))),
                    1.0,
                )
                tolerance = (
                    10.0
                    * np.finfo(float).eps
                    * max(len(eigenvalues), 1)
                    * eigenvalue_scale
                )
            else:
                tolerance = float(zero_tolerance)
                if not np.isfinite(tolerance) or tolerance < 0:
                    raise ValueError(
                        "zero_tolerance must be finite and nonnegative"
                    )

            # Significant negative eigenvalues indicate an unstable model.
            if np.any(eigenvalues < -tolerance):
                raise ValueError(
                    "modal analysis found negative eigenvalues; the model "
                    "is unstable or insufficiently constrained"
                )

            # Remove zero and round-off-scale rigid-body eigenpairs.
            positive_modes = eigenvalues > tolerance
            eigenvalues = eigenvalues[positive_modes]
            eigenvectors = eigenvectors[:, positive_modes]

            # Later samples are matched against the first sample before slicing.
            if track_modes and sample_eigenvectors:
                mode_order = self._mode_order(
                    sample_eigenvectors[0],
                    eigenvectors,
                )
                eigenvalues = eigenvalues[mode_order]
                eigenvectors = eigenvectors[:, mode_order]
                eigenvectors = self._phase_align(
                    sample_eigenvectors[0],
                    eigenvectors,
                )
            elif requested_modes is not None:
                # The reference sample selects the branches that will be tracked.
                if len(eigenvalues) < requested_modes:
                    raise ValueError(
                        f"Requested {requested_modes} modes but only "
                        f"{len(eigenvalues)} positive modes are available"
                    )
                eigenvalues = eigenvalues[:requested_modes]
                eigenvectors = eigenvectors[:, :requested_modes]

            # Preserve the tracked solution for stacking and later matching.
            sample_eigenvalues.append(eigenvalues)
            sample_eigenvectors.append(eigenvectors)

        # Convert angular eigenvalues lambda=omega^2 into frequency in hertz.
        sample_frequencies = [
            np.sqrt(eigenvalues) / (2.0 * np.pi)
            for eigenvalues in sample_eigenvalues
        ]

        return {
            "speed": speed,
            "eigenvalues": self._stack(sample_eigenvalues),
            "frequencies": self._stack(sample_frequencies),
            "eigenvectors": self._stack(sample_eigenvectors),
        }

    def _solve_campbell(
        self,
        model,
        free_dofs,
        *,
        speeds,
        modes=None,
        track_modes=True,
        harmonics=None,
    ):
        """Solve damped gyroscopic state problems over angular speeds.

        Parameters
        ----------
        model : BuiltModel
            Numerical model or stochastic ensemble.
        free_dofs : numpy.ndarray
            Retained global DOF indices.
        speeds : array_like
            Rotor angular speeds in radians per second.
        modes : int or None, optional
            Number of positive-imaginary branches to return.
        track_modes : bool, optional
            Track branches across speeds and stochastic samples using the MAC.
        harmonics : array_like or None, optional
            Positive harmonic indices for critical-speed detection. For index
            ``h``, the synchronous line is ``h * speed / (2*pi)``. Crossings
            are linearly interpolated only between the supplied speed points;
            the solver never refines the speed grid.

        Returns
        -------
        dict
            Speeds, complex state eigenvalues, frequencies, and eigenvectors.
        """
        # Import the general complex eigensolver only for Campbell analysis.
        from scipy.linalg import eig

        # Convert and validate the requested speed grid.
        angular_speeds = np.asarray(speeds, dtype=float)
        if (
            angular_speeds.ndim != 1
            or not angular_speeds.size
            or not np.isfinite(angular_speeds).all()
        ):
            raise ValueError(
                "speeds must be a non-empty one-dimensional sequence"
            )

        # Validate branch-selection and tracking options.
        requested_modes = self._validate_modes(modes)
        if not isinstance(track_modes, bool):
            raise TypeError("track_modes must be a boolean")
        harmonic_indices = self._validate_harmonics(harmonics)

        # At zero spin, an axisymmetric rotor can have repeated forward- and
        # backward-whirl modes whose individual eigenvectors are not unique.
        # Pairwise MAC may then initialize inconsistent branch labels.
        if track_modes and angular_speeds[0] == 0.0:
            warnings.warn(
                "Campbell mode tracking starts at exactly zero rotor speed. "
                "In axisymmetric models, degenerate forward/backward modes "
                "can make pairwise MAC correspondence ambiguous and produce "
                "irregular branches or misleading stochastic bands. If the "
                "diagram looks unexpected and the stationary point is not "
                "required, start at a small nonzero speed, for example "
                "1 rpm (approximately 0.10472 rad/s).",
                UserWarning,
                stacklevel=3,
            )

        if harmonic_indices.size:
            if not track_modes:
                raise ValueError(
                    "critical-speed calculation requires track_modes=True"
                )
            if len(angular_speeds) < 2:
                raise ValueError(
                    "critical-speed calculation requires at least two speeds"
                )
            if np.any(np.diff(angular_speeds) <= 0.0):
                raise ValueError(
                    "critical-speed calculation requires strictly increasing "
                    "speeds"
                )

        # Resolve every requested operating point before the first state-space
        # eigensolution. This also rejects extrapolation without returning a
        # partially solved Campbell sweep.
        stiffness_deltas = self._speed_matrix_deltas(
            model,
            model._stiffness_speed_contributions,
            angular_speeds,
        )
        damping_deltas = self._speed_matrix_deltas(
            model,
            model._damping_speed_contributions,
            angular_speeds,
        )

        # Collect speed-dependent solutions for every stochastic sample.
        all_sample_values = []
        all_sample_vectors = []

        for sample in range(model.samples):
            
            print(f"Analyzing Sample {sample}")

            # Extract this realization's reduced equation matrices.
            mass_matrix = self._reduced(model.M, sample, free_dofs)
            gyroscopic_matrix = self._reduced(model.G, sample, free_dofs)

            # Prepare constant state-space blocks.
            reduced_size = len(free_dofs)
            zero_matrix = np.zeros_like(mass_matrix)
            identity_matrix = np.eye(reduced_size)

            # Accumulate one eigensolution for every speed.
            speed_values = []
            speed_vectors = []

            for speed_index, angular_speed in enumerate(angular_speeds):
                print(f"    Speed: {angular_speed:5.0f} rad/s")
                stiffness_matrix = self._reduced_at_speed(
                    model.K,
                    stiffness_deltas,
                    sample,
                    speed_index,
                    free_dofs,
                )
                damping_matrix = self._reduced_at_speed(
                    model.C,
                    damping_deltas,
                    sample,
                    speed_index,
                    free_dofs,
                )
                # Form the first-order state matrix without explicitly inverting M.
                lower_left = -np.linalg.solve(
                    mass_matrix,
                    stiffness_matrix,
                )
                lower_right = -np.linalg.solve(
                    mass_matrix,
                    damping_matrix + angular_speed * gyroscopic_matrix,
                )
                state_matrix = np.block(
                    [
                        [zero_matrix, identity_matrix],
                        [lower_left, lower_right],
                    ]
                )

                # Solve for complex state eigenvalues and right eigenvectors.
                eigenvalues, eigenvectors = eig(
                    state_matrix,
                    check_finite=False,
                    right=True,
                )

                # Retain one eigenvalue from each complex-conjugate pair.
                positive_imaginary = np.imag(eigenvalues) > 1.0e-8
                eigenvalues = eigenvalues[positive_imaginary]
                eigenvectors = eigenvectors[:, positive_imaginary]

                # Frequency ordering defines the initial reference branches.
                frequency_order = np.argsort(np.imag(eigenvalues))
                eigenvalues = eigenvalues[frequency_order]
                eigenvectors = eigenvectors[:, frequency_order]

                # Track later speeds against displacement components at prior speed.
                if track_modes and speed_vectors:
                    mode_order = self._mode_order(
                        speed_vectors[-1][:reduced_size],
                        eigenvectors[:reduced_size],
                    )
                    eigenvalues = eigenvalues[mode_order]
                    eigenvectors = eigenvectors[:, mode_order]
                    eigenvectors = self._phase_align(
                        speed_vectors[-1],
                        eigenvectors,
                    )
                elif requested_modes is not None:
                    # Select reference branches only after all candidates exist.
                    if len(eigenvalues) < requested_modes:
                        raise ValueError(
                            f"Requested {requested_modes} modes but only "
                            f"{len(eigenvalues)} are available at speed "
                            f"{angular_speed}"
                        )
                    eigenvalues = eigenvalues[:requested_modes]
                    eigenvectors = eigenvectors[:, :requested_modes]

                speed_values.append(eigenvalues)
                speed_vectors.append(eigenvectors)

            # A sample must produce a rectangular speed-by-mode result.
            all_sample_values.append(self._stack(speed_values))
            all_sample_vectors.append(self._stack(speed_vectors))

        # Match stochastic realizations to corresponding first-sample branches.
        if track_modes and len(all_sample_values) > 1:
            reference_vectors = all_sample_vectors[0]

            for sample in range(1, len(all_sample_values)):
                for speed_index in range(len(angular_speeds)):
                    mode_order = self._mode_order(
                        reference_vectors[speed_index, :reduced_size],
                        all_sample_vectors[sample][
                            speed_index, :reduced_size
                        ],
                    )

                    # Apply the same permutation to eigenvalues and eigenvectors.
                    all_sample_values[sample][speed_index] = (
                        all_sample_values[sample][speed_index, mode_order]
                    )
                    all_sample_vectors[sample][speed_index] = (
                        all_sample_vectors[sample][speed_index][:, mode_order]
                    )

                    # Remove arbitrary complex phase relative to the reference.
                    all_sample_vectors[sample][speed_index] = self._phase_align(
                        reference_vectors[speed_index],
                        all_sample_vectors[sample][speed_index],
                    )

        # Stack the sample dimension after all branch matching is complete.
        eigenvalues = self._stack(all_sample_values)
        eigenvectors = self._stack(all_sample_vectors)
        frequencies = np.abs(np.imag(eigenvalues)) / (2.0 * np.pi)

        # Intersect each tracked modal branch with every requested synchronous
        # harmonic. The result retains crossing order and uses NaN padding so
        # samples with missing crossings remain rectangular numerical arrays.
        critical_speeds, critical_speed_counts = self._critical_speeds(
            angular_speeds,
            frequencies,
            harmonic_indices,
        )

        # Report speed in both angular and cyclic units for clear comparisons.
        return {
            "speeds": angular_speeds,
            "speeds_hz": angular_speeds / (2.0 * np.pi),
            "eigenvalues": eigenvalues,
            "frequencies": frequencies,
            "eigenvectors": eigenvectors,
            "track_modes": track_modes,
            "harmonics": harmonic_indices,
            "critical_speeds": critical_speeds,
            "critical_speeds_hz": critical_speeds / (2.0 * np.pi),
            "critical_speed_counts": critical_speed_counts,
        }

    @staticmethod
    def _validate_harmonics(harmonics):
        """Return validated positive harmonic indices."""
        if harmonics is None:
            return np.empty(0, dtype=float)
        if isinstance(harmonics, (bool, np.bool_)):
            raise TypeError(
                "harmonics must contain only positive finite numbers"
            )

        try:
            raw_harmonics = np.asarray(harmonics, dtype=object)
            if any(
                isinstance(value, (bool, np.bool_))
                for value in raw_harmonics.flat
            ):
                raise TypeError
            harmonic_indices = np.asarray(harmonics, dtype=float)
        except (TypeError, ValueError):
            raise TypeError(
                "harmonics must contain only positive finite numbers"
            ) from None
        if harmonic_indices.ndim == 0:
            harmonic_indices = harmonic_indices.reshape(1)
        if harmonic_indices.ndim != 1:
            raise ValueError("harmonics must be a one-dimensional sequence")
        if not harmonic_indices.size:
            return harmonic_indices
        if not np.isfinite(harmonic_indices).all() or np.any(
            harmonic_indices <= 0.0
        ):
            raise ValueError("harmonics must contain positive finite values")
        if len(np.unique(harmonic_indices)) != len(harmonic_indices):
            raise ValueError("harmonics must not contain duplicate indices")
        return harmonic_indices

    @staticmethod
    def _critical_speeds(speeds, frequencies, harmonics):
        """Intersect tracked Campbell branches with harmonic lines.

        Detection uses only exact grid-point matches and sign changes between
        adjacent supplied speeds. Each sign-change root is obtained by linear
        interpolation; no model evaluation or grid refinement is performed.
        """
        sample_count, _, mode_count = frequencies.shape
        harmonic_count = len(harmonics)
        crossing_lists = [
            [
                [[] for _ in range(mode_count)]
                for _ in range(harmonic_count)
            ]
            for _ in range(sample_count)
        ]
        counts = np.zeros(
            (sample_count, harmonic_count, mode_count), dtype=int
        )

        for sample in range(sample_count):
            for harmonic_position, h in enumerate(harmonics):
                synchronous = h * speeds / (2.0 * np.pi)
                for mode in range(mode_count):
                    residual = frequencies[sample, :, mode] - synchronous
                    scale = max(float(np.max(np.abs(residual))), 1.0)
                    tolerance = (
                        10.0 * np.finfo(float).eps * len(speeds) * scale
                    )
                    is_zero = np.abs(residual) <= tolerance
                    roots = []

                    for index in range(len(speeds) - 1):
                        if is_zero[index]:
                            roots.append(float(speeds[index]))
                        if is_zero[index] or is_zero[index + 1]:
                            continue
                        if np.signbit(residual[index]) == np.signbit(
                            residual[index + 1]
                        ):
                            continue

                        # Linear interpolation of the residual on the user's
                        # existing interval; this is detection, not refinement.
                        fraction = -residual[index] / (
                            residual[index + 1] - residual[index]
                        )
                        roots.append(
                            float(
                                speeds[index]
                                + fraction
                                * (speeds[index + 1] - speeds[index])
                            )
                        )

                    if is_zero[-1]:
                        roots.append(float(speeds[-1]))

                    # Consecutive exact points can represent the same sampled
                    # intersection. Preserve order while removing duplicates.
                    unique_roots = []
                    for root in roots:
                        if not unique_roots or not np.isclose(
                            root,
                            unique_roots[-1],
                            rtol=1.0e-12,
                            atol=1.0e-12,
                        ):
                            unique_roots.append(root)

                    crossing_lists[sample][harmonic_position][mode] = (
                        unique_roots
                    )
                    counts[sample, harmonic_position, mode] = len(unique_roots)

        maximum_crossings = int(counts.max(initial=0))
        critical_speeds = np.full(
            (
                sample_count,
                harmonic_count,
                mode_count,
                maximum_crossings,
            ),
            np.nan,
        )
        for sample in range(sample_count):
            for harmonic_position in range(harmonic_count):
                for mode in range(mode_count):
                    roots = crossing_lists[sample][harmonic_position][mode]
                    critical_speeds[
                        sample,
                        harmonic_position,
                        mode,
                        : len(roots),
                    ] = roots

        return critical_speeds, counts

    def _solve_frequency_response(
        self,
        model,
        free_dofs,
        *,
        frequencies=None,
        force=None,
        load_ids=None,
        speed=0.0,
    ):
        """Solve steady harmonic response at requested frequencies.

        Parameters
        ----------
        model : BuiltModel
            Numerical model or stochastic ensemble.
        free_dofs : numpy.ndarray
            Retained global DOF indices.
        frequencies : array_like or None, optional
            Excitation frequencies in hertz. Required for an explicit force,
            :class:`HarmonicLoad`, or :class:`RotatingLoad`; omitted for
            :class:`SynchronousLoad`, whose frequency is order-derived.
        force : array_like or None, optional
            Complex force vector in full or reduced DOF ordering. Supply
            exactly one of ``force`` and ``load_ids``.
        load_ids : sequence of int or None, optional
            IDs of compatible harmonic load records to assemble. Selected
            synchronous loads must all share one harmonic order and cannot be
            mixed with ordinary harmonic or rotating loads.
        speed : float, optional
            Rotor angular speed in radians per second.

        Returns
        -------
        dict
            Frequencies and complex displacement responses.
        """
        # Gyroscopic damping and every SpeedDependentValue use the same signed
        # rotor angular speed, independently of the forcing-frequency grid.
        speed = self._scalar_speed(speed)

        # Preserve the raw-vector route while making declarative selection
        # explicit. Combining both sources would obscure accidental double
        # application, so callers must choose exactly one.
        if (force is None) == (load_ids is None):
            raise ValueError(
                "frequency_response requires exactly one of force or load_ids"
            )

        selected_load_ids = ()
        reduced_force = None
        harmonic_order = None
        if load_ids is None:
            if frequencies is None:
                raise ValueError(
                    "frequency_response requires frequencies for an explicit "
                    "force"
                )
            # One caller-supplied vector is shared by every realization.
            reduced_force = self._free_vector(
                force,
                model.ndof,
                free_dofs,
                "force",
                dtype=complex,
            )
        else:
            selected_load_ids = self._validate_load_ids(
                model,
                load_ids,
                (HarmonicLoad, SynchronousLoad, RotatingLoad),
                "frequency_response",
            )
            records = {load.id: load for load in model.definitions[0].loads}
            selected_records = [records[load_id] for load_id in selected_load_ids]
            synchronous_records = [
                load for load in selected_records
                if isinstance(load, SynchronousLoad)
            ]
            if synchronous_records and len(synchronous_records) != len(
                selected_records
            ):
                raise ValueError(
                    "frequency_response cannot mix SynchronousLoad with "
                    "HarmonicLoad or RotatingLoad records"
                )
            if synchronous_records:
                if frequencies is not None:
                    raise ValueError(
                        "frequencies must be omitted for SynchronousLoad records"
                    )
                harmonic_orders = {
                    float(load.harmonic_order) for load in synchronous_records
                }
                if len(harmonic_orders) != 1:
                    raise ValueError(
                        "selected SynchronousLoad records must share one "
                        "harmonic_order"
                    )
                harmonic_order = harmonic_orders.pop()
                frequencies = np.array(
                    [harmonic_order * abs(speed) / (2.0 * np.pi)]
                )
            elif frequencies is None:
                raise ValueError(
                    "frequency_response requires frequencies for HarmonicLoad "
                    "or RotatingLoad records"
                )

        # Validate either the caller-supplied grid or the single frequency
        # derived above from a synchronous load's harmonic order.
        frequencies = np.asarray(frequencies, dtype=float)
        if (
            frequencies.ndim != 1
            or not frequencies.size
            or not np.isfinite(frequencies).all()
            or np.any(frequencies < 0.0)
        ):
            raise ValueError(
                "frequencies must be a non-empty one-dimensional sequence "
                "of finite nonnegative values"
            )

        # Precompute every structural and declarative-load operating-point
        # value before the first dynamic-stiffness solve.
        evaluation_speeds = np.array([speed])
        stiffness_deltas = self._speed_matrix_deltas(
            model,
            model._stiffness_speed_contributions,
            evaluation_speeds,
        )
        damping_deltas = self._speed_matrix_deltas(
            model,
            model._damping_speed_contributions,
            evaluation_speeds,
        )
        declarative_forces = None
        if selected_load_ids:
            declarative_forces = np.stack(
                [
                    self._assemble_harmonic_loads(
                        model,
                        sample,
                        selected_load_ids,
                        speed=speed,
                    )[free_dofs]
                    for sample in range(model.samples)
                ]
            )

        # Solve the complex dynamic stiffness for each sample and frequency.
        all_responses = []
        angular_frequencies = 2.0 * np.pi * frequencies

        for sample in range(model.samples):

            print(f"Analyzing Sample {sample}")
            mass_matrix = self._reduced(model.M, sample, free_dofs)
            stiffness_matrix = self._reduced_at_speed(
                model.K,
                stiffness_deltas,
                sample,
                0,
                free_dofs,
            )
            damping_matrix = self._reduced_at_speed(
                model.C,
                damping_deltas,
                sample,
                0,
                free_dofs,
            )
            gyroscopic_matrix = self._reduced(model.G, sample, free_dofs)

            # Declarative parameters are resolved independently per sample.
            if selected_load_ids:
                sample_force = declarative_forces[sample]
            else:
                sample_force = reduced_force

            sample_responses = []
            for angular_frequency in angular_frequencies:

                print(f"    Angular Frequency: {angular_frequency:5.0f} rad/s")

                dynamic_stiffness = self._dynamic_stiffness(
                    stiffness_matrix,
                    mass_matrix,
                    damping_matrix,
                    gyroscopic_matrix,
                    speed,
                    angular_frequency,
                )
                sample_responses.append(
                    np.linalg.solve(dynamic_stiffness, sample_force)
                )

            all_responses.append(np.stack(sample_responses))

        result = {
            "frequencies": frequencies,
            "speed": speed,
            "load_ids": selected_load_ids,
            "response": np.stack(all_responses),
        }
        if harmonic_order is not None:
            result["harmonic_order"] = harmonic_order
        return result

    def _solve_unbalance_response(
        self,
        model,
        free_dofs,
        *,
        speeds,
        load_ids,
    ):
        """Solve synchronous response to selected mass-unbalance records.

        Parameters
        ----------
        model : BuiltModel
            Numerical model or stochastic ensemble.
        free_dofs : numpy.ndarray
            Retained global DOF indices.
        speeds : array_like
            Signed rotor angular speeds in radians per second.
        load_ids : sequence of int
            IDs of :class:`Unbalance` records to assemble and superimpose.

        Returns
        -------
        dict
            Speeds, synchronous frequencies, selected IDs, and complex
            displacement responses.
        """
        # Speed sign controls rotation sense and the gyroscopic term, while
        # its absolute value is the synchronous excitation frequency.
        speeds = np.asarray(speeds, dtype=float)
        if (
            speeds.ndim != 1
            or not speeds.size
            or not np.isfinite(speeds).all()
        ):
            raise ValueError(
                "speeds must be a non-empty one-dimensional sequence of "
                "finite values"
            )

        selected_load_ids = self._validate_load_ids(
            model,
            load_ids,
            (Unbalance,),
            "unbalance_response",
        )

        # Bearing tables are evaluated over the complete signed speed grid
        # before any synchronous response solve begins.
        stiffness_deltas = self._speed_matrix_deltas(
            model,
            model._stiffness_speed_contributions,
            speeds,
        )
        damping_deltas = self._speed_matrix_deltas(
            model,
            model._damping_speed_contributions,
            speeds,
        )

        all_responses = []
        angular_frequencies = np.abs(speeds)

        for sample in range(model.samples):

            print(f"Analyzing Sample {sample}")

            mass_matrix = self._reduced(model.M, sample, free_dofs)
            gyroscopic_matrix = self._reduced(model.G, sample, free_dofs)

            # Geometry, phase, and sampled mass eccentricity are constant
            # within one realization. Assemble one coefficient for each spin
            # sense, then apply Omega**2 at every requested point.
            positive_force_coefficient = self._assemble_unbalances(
                model,
                sample,
                selected_load_ids,
                1.0,
            )[free_dofs]
            negative_force_coefficient = self._assemble_unbalances(
                model,
                sample,
                selected_load_ids,
                -1.0,
            )[free_dofs]
            zero_force = np.zeros(len(free_dofs), dtype=complex)

            sample_responses = []
            for speed_index, (speed, angular_frequency) in enumerate(
                zip(speeds, angular_frequencies, strict=True)
            ):
                print(f"    Speed: {speed} rad/s")
                stiffness_matrix = self._reduced_at_speed(
                    model.K,
                    stiffness_deltas,
                    sample,
                    speed_index,
                    free_dofs,
                )
                damping_matrix = self._reduced_at_speed(
                    model.C,
                    damping_deltas,
                    sample,
                    speed_index,
                    free_dofs,
                )
                if speed > 0.0:
                    sample_force = speed**2 * positive_force_coefficient
                elif speed < 0.0:
                    sample_force = speed**2 * negative_force_coefficient
                else:
                    sample_force = zero_force
                dynamic_stiffness = self._dynamic_stiffness(
                    stiffness_matrix,
                    mass_matrix,
                    damping_matrix,
                    gyroscopic_matrix,
                    speed,
                    angular_frequency,
                )
                sample_responses.append(
                    np.linalg.solve(
                        dynamic_stiffness,
                        sample_force,
                    )
                )

            all_responses.append(np.stack(sample_responses))

        return {
            "speeds": speeds,
            "speeds_hz": speeds / (2.0 * np.pi),
            "frequencies": angular_frequencies / (2.0 * np.pi),
            "load_ids": selected_load_ids,
            "response": np.stack(all_responses),
        }

    @staticmethod
    def _validate_load_ids(model, load_ids, accepted_types, analysis):
        """Validate and preserve one explicit load-record selection."""
        try:
            requested_ids = list(load_ids)
        except TypeError:
            raise TypeError("load_ids must be a sequence of integers") from None

        if not requested_ids:
            raise ValueError("load_ids must contain at least one ID")
        if any(
            isinstance(load_id, (bool, np.bool_))
            or not isinstance(load_id, (int, np.integer))
            for load_id in requested_ids
        ):
            raise TypeError("load_ids must contain only integers")

        normalized_ids = tuple(int(load_id) for load_id in requested_ids)
        if len(set(normalized_ids)) != len(normalized_ids):
            raise ValueError("load_ids must not contain duplicate IDs")

        available_loads = {
            load.id: load for load in model.definitions[0].loads
        }
        unknown_ids = [
            load_id for load_id in normalized_ids
            if load_id not in available_loads
        ]
        if unknown_ids:
            raise ValueError(f"load_ids contains unknown load ID {unknown_ids[0]}")

        invalid_loads = [
            available_loads[load_id]
            for load_id in normalized_ids
            if not isinstance(available_loads[load_id], accepted_types)
        ]
        if invalid_loads:
            accepted_names = " or ".join(
                accepted_type.__name__ for accepted_type in accepted_types
            )
            invalid_name = type(invalid_loads[0]).__name__
            raise ValueError(
                f"{analysis} accepts only {accepted_names} records, not "
                f"{invalid_name}"
            )

        return normalized_ids

    @staticmethod
    def _load_coordinate_systems(definition):
        """Return local-to-global rotations indexed by coordinate-system ID."""
        rotations = {0: np.eye(3)}
        for coordinate_system in definition.coordinate_systems:
            rotations[coordinate_system.id] = _coordinate_system_matrix(
                coordinate_system
            )[1]
        return rotations

    @staticmethod
    def _add_nodal_load(global_load, grid_position, rotation, local_load):
        """Rotate and accumulate one six-component nodal load in place."""
        offset = 6 * grid_position
        global_load[offset : offset + 3] += rotation @ local_load[:3]
        global_load[offset + 3 : offset + 6] += rotation @ local_load[3:]

    @staticmethod
    def _load_value_at_speed(value, speed, path):
        """Return one scalar load value at the requested rotor speed."""
        if isinstance(value, SpeedDependentValue):
            if speed is None:
                raise ValueError(f"{path} requires a rotor speed")
            return float(
                _interpolate_speed_dependent_value(
                    value,
                    np.array([speed]),
                    path,
                )[0]
            )
        return float(value)

    @classmethod
    def _assemble_harmonic_loads(
        cls,
        model,
        sample,
        load_ids,
        *,
        speed=None,
    ):
        """Assemble component-vector or rotating loads into global DOFs."""
        definition = model.definitions[sample]
        records = {load.id: load for load in definition.loads}
        rotations = cls._load_coordinate_systems(definition)
        grid_positions = {
            grid_id: position
            for position, grid_id in enumerate(model.grid_ids)
        }
        global_load = np.zeros(model.ndof, dtype=complex)

        for load_id in load_ids:
            load = records[load_id]
            rotation = rotations[load.coordinate_system]

            if isinstance(load, (HarmonicLoad, SynchronousLoad)):
                components = tuple(load.force) + tuple(load.moment)
                local_load = np.asarray(
                    [
                        cls._load_value_at_speed(
                            value,
                            speed,
                            (
                                f"Harmonic load {load.id} "
                                f"{'force' if component < 3 else 'moment'}"
                                f"[{component if component < 3 else component - 3}]"
                            ),
                        )
                        for component, value in enumerate(components)
                    ],
                    dtype=complex,
                )
                local_load *= np.exp(1j * float(load.phase))
            else:
                sense = 1.0 if load.sense == "positive" else -1.0
                magnitude = cls._load_value_at_speed(
                    load.force_magnitude,
                    speed,
                    f"Rotating load {load.id} force_magnitude",
                )
                radial_phasor = magnitude * np.exp(
                    1j * sense * float(load.phase)
                )
                local_load = np.zeros(6, dtype=complex)
                local_load[0] = radial_phasor
                local_load[1] = -1j * sense * radial_phasor

            cls._add_nodal_load(
                global_load,
                grid_positions[load.grid],
                rotation,
                local_load,
            )

        return global_load

    @classmethod
    def _assemble_unbalances(cls, model, sample, load_ids, speed):
        """Assemble synchronous unbalance forces at one signed rotor speed."""
        definition = model.definitions[sample]
        records = {load.id: load for load in definition.loads}
        rotations = cls._load_coordinate_systems(definition)
        grid_positions = {
            grid_id: position
            for position, grid_id in enumerate(model.grid_ids)
        }
        spin_axis = np.asarray(definition.spin_axis, dtype=float)
        spin_axis /= np.linalg.norm(spin_axis)
        global_load = np.zeros(model.ndof, dtype=complex)

        # At zero speed every unbalance force is exactly zero, independently
        # of phase or frame orientation.
        if speed == 0.0:
            return global_load

        for load_id in load_ids:
            load = records[load_id]
            rotation = rotations[load.coordinate_system]
            axis_alignment = float(rotation[:, 2] @ spin_axis)
            sense = 1.0 if speed * axis_alignment > 0.0 else -1.0
            force_magnitude = float(load.mass_eccentricity) * speed**2
            radial_phasor = force_magnitude * np.exp(
                1j * sense * float(load.phase)
            )
            local_load = np.zeros(6, dtype=complex)
            local_load[0] = radial_phasor
            local_load[1] = -1j * sense * radial_phasor

            cls._add_nodal_load(
                global_load,
                grid_positions[load.grid],
                rotation,
                local_load,
            )

        return global_load

    @staticmethod
    def _dynamic_stiffness(
        stiffness,
        mass,
        damping,
        gyroscopic,
        speed,
        angular_frequency,
    ):
        """Return ``K - omega**2 M + i omega (C + Omega G)``."""
        return (
            stiffness
            - angular_frequency**2 * mass
            + 1j * angular_frequency * (damping + speed * gyroscopic)
        )

    def _solve_time_response(
        self,
        model,
        free_dofs,
        *,
        times,
        force,
        speed=0.0,
        initial_state=None,
    ):
        """Integrate the second-order rotor equation in time.

        Parameters
        ----------
        model : BuiltModel
            Numerical model or stochastic ensemble.
        free_dofs : numpy.ndarray
            Retained global DOF indices.
        times : array_like
            Strictly increasing output times.
        force : array_like or callable
            Constant force vector or ``force(time)`` callback.
        speed : float, optional
            Rotor angular speed in radians per second.
        initial_state : array_like or None, optional
            Concatenated displacement and velocity in full or reduced order.

        Returns
        -------
        dict
            Times, displacements, and velocities for every sample.
        """
        # Import the integrator only when time response is requested.
        from scipy.integrate import solve_ivp

        # Validate the requested output-time grid.
        times = np.asarray(times, dtype=float)
        if (
            times.ndim != 1
            or len(times) < 2
            or not np.isfinite(times).all()
            or np.any(np.diff(times) <= 0)
        ):
            raise ValueError(
                "times must contain at least two strictly increasing values"
            )

        # The constant rotor speed selects one structural operating point for
        # the complete integration.
        speed = self._scalar_speed(speed)

        # The reduced state concatenates q and q_dot.
        reduced_size = len(free_dofs)

        if initial_state is None:
            # Default to zero displacement and zero velocity.
            reduced_initial_state = np.zeros(2 * reduced_size)
        else:
            reduced_initial_state = np.asarray(initial_state, dtype=float)

            # Reduce a full [q, q_dot] vector when one is supplied.
            if reduced_initial_state.shape == (2 * model.ndof,):
                full_displacement = reduced_initial_state[: model.ndof]
                full_velocity = reduced_initial_state[model.ndof :]
                reduced_initial_state = np.concatenate(
                    (
                        full_displacement[free_dofs],
                        full_velocity[free_dofs],
                    )
                )

            # Otherwise the state must already use reduced ordering.
            if reduced_initial_state.shape != (2 * reduced_size,):
                raise ValueError(
                    "initial_state must contain displacement and velocity "
                    "at every free DOF"
                )
            if not np.isfinite(reduced_initial_state).all():
                raise ValueError(
                    "initial_state must contain only finite values"
                )

        # Validate a constant force once; callbacks are validated per call.
        constant_force = None
        if not callable(force):
            constant_force = self._free_vector(
                force,
                model.ndof,
                free_dofs,
                "force",
                dtype=float,
            )

        # Reject any out-of-range structural interpolation before invoking the
        # first time integrator.
        evaluation_speeds = np.array([speed])
        stiffness_deltas = self._speed_matrix_deltas(
            model,
            model._stiffness_speed_contributions,
            evaluation_speeds,
        )
        damping_deltas = self._speed_matrix_deltas(
            model,
            model._damping_speed_contributions,
            evaluation_speeds,
        )

        # Integrate every deterministic realization independently.
        sample_states = []

        for sample in range(model.samples):
            
            print(f"Analyzing Sample {sample}")

            mass_matrix = self._reduced(model.M, sample, free_dofs)
            stiffness_matrix = self._reduced_at_speed(
                model.K,
                stiffness_deltas,
                sample,
                0,
                free_dofs,
            )
            effective_damping = (
                self._reduced_at_speed(
                    model.C,
                    damping_deltas,
                    sample,
                    0,
                    free_dofs,
                )
                + speed * self._reduced(model.G, sample, free_dofs)
            )

            def right_hand_side(time, state):
                """Evaluate the first-order state derivative.

                Parameters
                ----------
                time : float
                    Current integration time.
                state : numpy.ndarray
                    Concatenated reduced displacement and velocity.

                Returns
                -------
                numpy.ndarray
                    Concatenated velocity and acceleration.
                """
                # Evaluate or reuse the external load vector.
                if callable(force):
                    load = self._free_vector(
                        force(time),
                        model.ndof,
                        free_dofs,
                        "force",
                        dtype=float,
                    )
                else:
                    load = constant_force

                # Split the state without allocating named copies.
                displacement = state[:reduced_size]
                velocity = state[reduced_size:]

                # Rearrange M q_ddot = f - C_eff q_dot - K q.
                residual_force = (
                    load
                    - effective_damping @ velocity
                    - stiffness_matrix @ displacement
                )
                acceleration = np.linalg.solve(mass_matrix, residual_force)

                return np.concatenate((velocity, acceleration))

            # Integrate only across the requested time interval.
            solution = solve_ivp(
                right_hand_side,
                (times[0], times[-1]),
                reduced_initial_state,
                t_eval=times,
            )

            # Propagate integrator failure as a clear runtime error.
            if not solution.success:
                raise RuntimeError(solution.message)

            # Transpose to the conventional (time, state) result ordering.
            sample_states.append(solution.y.T)

        # Add the leading sample dimension after every integration succeeds.
        states = np.stack(sample_states)

        return {
            "times": times,
            "displacement": states[:, :, :reduced_size],
            "velocity": states[:, :, reduced_size:],
        }

    @staticmethod
    def _free_vector(value, number_of_dofs, free_dofs, name, *, dtype):
        """Validate and reduce a full or free-DOF vector.

        Parameters
        ----------
        value : array_like
            Full or already reduced vector.
        number_of_dofs : int
            Full global vector length.
        free_dofs : numpy.ndarray
            Retained global indices.
        name : str
            Field name used in validation messages.
        dtype : numpy dtype
            Required output data type.

        Returns
        -------
        numpy.ndarray
            Finite vector in reduced DOF order.
        """
        # Convert the caller's vector to the required real or complex type.
        vector = np.asarray(value, dtype=dtype)

        # Reduce a full global vector when necessary.
        if vector.shape == (number_of_dofs,):
            vector = vector[free_dofs]

        # Any other length is ambiguous and therefore rejected.
        if vector.shape != (len(free_dofs),):
            raise ValueError(
                f"{name} must contain either ndof or number-of-free-dofs entries"
            )

        # Linear solvers cannot consume NaN or infinite loads safely.
        if not np.isfinite(vector).all():
            raise ValueError(f"{name} must contain only finite values")

        return vector

    @staticmethod
    def _validate_modes(modes):
        """Validate an optional requested mode count.

        Parameters
        ----------
        modes : int or None
            Requested number of modes.

        Returns
        -------
        int or None
            Normalized Python integer or ``None``.
        """
        # None means that all available positive modes should be retained.
        if modes is None:
            return None

        # Reject Booleans, non-integral values, zero, and negative counts.
        if (
            isinstance(modes, (bool, np.bool_))
            or not isinstance(modes, (int, np.integer))
            or modes < 1
        ):
            raise ValueError("modes must be a positive integer")

        return int(modes)

    @staticmethod
    def _mode_order(reference, candidate):
        """Match candidate modes to reference modes using the MAC.

        Parameters
        ----------
        reference : numpy.ndarray
            Reference mode shapes stored by column.
        candidate : numpy.ndarray
            Candidate mode shapes stored by column.

        Returns
        -------
        numpy.ndarray
            Candidate-column indices aligned with reference columns.
        """
        # Import assignment machinery only when tracking is enabled.
        from scipy.optimize import linear_sum_assignment

        # Every tracked branch needs at least one candidate branch.
        if reference.shape[1] > candidate.shape[1]:
            raise ValueError(
                "Fewer candidate modes are available than tracked branches; "
                "request fewer modes or improve constraints"
            )

        # Compute the normalized squared modal inner products.
        numerator = np.abs(reference.conj().T @ candidate) ** 2
        reference_norm = np.sum(np.abs(reference) ** 2, axis=0)
        candidate_norm = np.sum(np.abs(candidate) ** 2, axis=0)
        denominator = np.outer(reference_norm, candidate_norm)

        # Zero-norm columns receive a zero MAC instead of producing NaN.
        modal_assurance = np.divide(
            numerator,
            denominator,
            out=np.zeros_like(numerator, dtype=float),
            where=denominator > 0,
        )

        # Maximize total MAC by minimizing its negative assignment cost.
        reference_rows, candidate_columns = linear_sum_assignment(
            -modal_assurance
        )

        # Return candidates in original reference-branch order.
        return candidate_columns[np.argsort(reference_rows)]

    @staticmethod
    def _phase_align(reference, candidate):
        """Remove arbitrary eigenvector signs or complex phases.

        Parameters
        ----------
        reference : numpy.ndarray
            Reference eigenvectors stored by column.
        candidate : numpy.ndarray
            Matched candidate eigenvectors stored by column.

        Returns
        -------
        numpy.ndarray
            Candidate vectors phase-aligned with the reference.
        """
        # Column-wise complex overlap contains the relative phase.
        overlap = np.sum(reference.conj() * candidate, axis=0)

        # Default to no phase change for numerically orthogonal vectors.
        factors = np.ones_like(overlap)
        nonzero_overlap = np.abs(overlap) > np.finfo(float).eps

        # Rotate each candidate so its reference overlap becomes real-positive.
        factors[nonzero_overlap] = (
            np.conj(overlap[nonzero_overlap])
            / np.abs(overlap[nonzero_overlap])
        )

        return candidate * factors

    @staticmethod
    def _stack(values):
        """Stack equal-shaped analysis arrays and reject ragged results.

        Parameters
        ----------
        values : sequence of numpy.ndarray
            Arrays expected to have identical shapes.

        Returns
        -------
        numpy.ndarray
            Arrays stacked along a new leading dimension.

        Raises
        ------
        ValueError
            If mode counts or other result dimensions are inconsistent.
        """
        try:
            return np.stack(values)
        except ValueError:
            raise ValueError(
                "Analysis produced inconsistent result shapes; request a "
                "fixed number of modes"
            ) from None
