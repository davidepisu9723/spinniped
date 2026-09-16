"""Unified analyses for deterministic models and stochastic ensembles."""

import numpy as np

from .builder import BuiltModel


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
            One of ``modal``, ``campbell``, ``frequency_response``, or
            ``time_response``.
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

        # Resolve the analysis name to its private implementation.
        try:
            normalized_analysis = analysis.lower()
            handler = getattr(self, f"_solve_{normalized_analysis}")
        except (AttributeError, TypeError):
            raise ValueError(
                f"Unknown analysis {analysis!r}; choose from {self._ANALYSES}"
            ) from None

        # Convert fixed global DOFs into a sorted free-DOF array.
        free_dofs = self._free_dofs(numerical_model.ndof, fixed_dofs)

        # Delegate only the options understood by the selected handler.
        result = handler(numerical_model, free_dofs, **options)

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

    def _solve_modal(
        self,
        model,
        free_dofs,
        modes=None,
        zero_tolerance=None,
        track_modes=True,
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

        # Store one eigensolution per realization.
        sample_eigenvalues = []
        sample_eigenvectors = []

        for sample in range(model.samples):
            # Extract the free-free stiffness and mass matrices.
            stiffness_matrix = self._reduced(model.K, sample, free_dofs)
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

        # Collect speed-dependent solutions for every stochastic sample.
        all_sample_values = []
        all_sample_vectors = []

        for sample in range(model.samples):
            # Extract this realization's reduced equation matrices.
            mass_matrix = self._reduced(model.M, sample, free_dofs)
            stiffness_matrix = self._reduced(model.K, sample, free_dofs)
            damping_matrix = self._reduced(model.C, sample, free_dofs)
            gyroscopic_matrix = self._reduced(model.G, sample, free_dofs)

            # Prepare constant state-space blocks.
            reduced_size = len(free_dofs)
            zero_matrix = np.zeros_like(mass_matrix)
            identity_matrix = np.eye(reduced_size)

            # Accumulate one eigensolution for every speed.
            speed_values = []
            speed_vectors = []

            for angular_speed in angular_speeds:
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

        # Report speed in both angular and cyclic units for clear comparisons.
        return {
            "speeds": angular_speeds,
            "speeds_hz": angular_speeds / (2.0 * np.pi),
            "eigenvalues": eigenvalues,
            "frequencies": np.abs(np.imag(eigenvalues)) / (2.0 * np.pi),
            "eigenvectors": eigenvectors,
        }

    def _solve_frequency_response(
        self,
        model,
        free_dofs,
        *,
        frequencies,
        force,
        speed=0.0,
    ):
        """Solve steady harmonic response at requested frequencies.

        Parameters
        ----------
        model : BuiltModel
            Numerical model or stochastic ensemble.
        free_dofs : numpy.ndarray
            Retained global DOF indices.
        frequencies : array_like
            Excitation frequencies in hertz.
        force : array_like
            Complex force vector in full or reduced DOF ordering.
        speed : float, optional
            Rotor angular speed in radians per second.

        Returns
        -------
        dict
            Frequencies and complex displacement responses.
        """
        # Validate the cyclic-frequency grid.
        frequencies = np.asarray(frequencies, dtype=float)
        if (
            frequencies.ndim != 1
            or not frequencies.size
            or not np.isfinite(frequencies).all()
        ):
            raise ValueError(
                "frequencies must be a non-empty one-dimensional sequence"
            )

        # Gyroscopic damping is scaled by a finite angular speed.
        if not np.isfinite(speed):
            raise ValueError("speed must be finite")

        # Accept either a full global force or an already reduced force.
        reduced_force = self._free_vector(
            force,
            model.ndof,
            free_dofs,
            "force",
            dtype=complex,
        )

        # Solve the complex dynamic stiffness for each sample and frequency.
        all_responses = []
        angular_frequencies = 2.0 * np.pi * frequencies

        for sample in range(model.samples):
            mass_matrix = self._reduced(model.M, sample, free_dofs)
            stiffness_matrix = self._reduced(model.K, sample, free_dofs)
            effective_damping = (
                self._reduced(model.C, sample, free_dofs)
                + speed * self._reduced(model.G, sample, free_dofs)
            )

            sample_responses = []
            for angular_frequency in angular_frequencies:
                # Z(omega) = K - omega^2 M + i omega (C + Omega G).
                dynamic_stiffness = (
                    stiffness_matrix
                    - angular_frequency**2 * mass_matrix
                    + 1j * angular_frequency * effective_damping
                )
                sample_responses.append(
                    np.linalg.solve(dynamic_stiffness, reduced_force)
                )

            all_responses.append(np.stack(sample_responses))

        return {
            "frequencies": frequencies,
            "response": np.stack(all_responses),
        }

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

        # Gyroscopic damping requires a finite angular speed.
        if not np.isfinite(speed):
            raise ValueError("speed must be finite")

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

        # Integrate every deterministic realization independently.
        sample_states = []

        for sample in range(model.samples):
            mass_matrix = self._reduced(model.M, sample, free_dofs)
            stiffness_matrix = self._reduced(model.K, sample, free_dofs)
            effective_damping = (
                self._reduced(model.C, sample, free_dofs)
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
