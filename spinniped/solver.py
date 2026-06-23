import numpy as np

# Import custom structured dtypes used by the project.
from spinniped import dtypes
from spinniped import node


class Solver:
    """
    Solver class for rotor Finite Element Analysis (FEA).

    The solver stores rotor elements, assembles the global stiffness and mass
    matrices, applies boundary conditions, and solves the modal eigenvalue
    problem.

    Notes
    -----
    The current implementation assumes 4 Degree Of Freedom (DOF) per node.

    The expected DOF ordering at each node is:

        0 -> transverse displacement in x
        1 -> rotation related to x-plane bending
        2 -> transverse displacement in y
        3 -> rotation related to y-plane bending

    Therefore, the global DOF indices for node `n` (internal sequential id) are:

        4*n, 4*n + 1, 4*n + 2, 4*n + 3
    """

    def __init__(self):
        """
        Initialize an empty solver.

        This method creates the main containers used by the solver.

        Attributes
        ----------
        K_global : ndarray or None
            Full global stiffness matrix. It is created during assembly.

        M_global : ndarray or None
            Full global mass matrix. It is created during assembly.

        ListOfElements : list
            List of finite elements to be assembled.

        meshgrid : ndarray or Grid-like object
            Mesh grid containing the nodes used by the solver.

        fixed_dofs : list of int
            Global DOF indices constrained by boundary conditions.

        free_dofs : list of int
            Global DOF indices that remain active after applying constraints.

        K_ff : ndarray or None
            Reduced stiffness matrix containing only free DOFs.

        M_ff : ndarray or None
            Reduced mass matrix containing only free DOFs.
        """

        # Global stiffness matrix.
        # It is initialized later by `assemble_global_matrices()`.
        self.K_global = None

        # Global mass matrix.
        # It is initialized later by `assemble_global_matrices()`.
        self.M_global = None

        # List containing all elements added to the solver.
        self.ListOfElements = []

        # Mesh grid used by the solver.
        self.meshgrid = None

        # List of constrained global DOF indices.
        self.fixed_dofs = []

        # List of unconstrained global DOF indices.
        self.free_dofs = []

        # Reduced stiffness matrix after boundary conditions.
        self.K_ff = None

        # Reduced mass matrix after boundary conditions.
        self.M_ff = None

    def _get_internal_id(self, user_id):
        """
        Return the internal node index associated with a user-defined node ID.

        Parameters
        ----------
        user_id : int
            Node ID provided by the user.

        Returns
        -------
        int
            Internal zero-based node index used by the solver.

        Notes
        -----
        The solver uses internal node indices to assemble matrices.

        User IDs may not be zero-based or consecutive, so they must be mapped
        to internal indices before computing global DOF positions.
        """

        # Find the position of the user node ID in the sorted mesh.
        idx = np.searchsorted(self.meshgrid.Nodes['id'], user_id)

        # Check that the node ID really exists.
        if idx >= self.meshgrid.Nodes.size or self.meshgrid.Nodes['id'][idx] != user_id:
            raise ValueError(f"Node ID {user_id} not found in mesh.")

        # Return the internal node index.
        return self.meshgrid.Nodes['internal_id'][idx]

    def assemble_global_matrices(self):
        """
        Assemble the global stiffness and mass matrices.

        This method loops over all elements stored in `ListOfElements`.
        Each element contributes its local matrices to the correct position
        inside the global matrices.

        Supported element types
        -----------------------
        shaft
            Adds both stiffness and mass contributions.

        ball-bearing
            Adds only stiffness contribution.

        Notes
        -----
        Shaft elements are assumed to connect two nodes.

        Ball-bearing elements are currently assumed to act on one node and
        only on the first two DOFs of that node.
        """

        # Count the number of elements stored in the solver.
        ne = len(self.ListOfElements)

        # Count the number of nodes in the mesh.
        nn = self.meshgrid.Nodes.shape[0]

        # Allocate the global stiffness matrix.
        # Matrix size is 6 DOFs per node.
        self.K_global = np.zeros((6 * nn, 6 * nn))

        # Allocate the global mass matrix.
        # Matrix size is 4 DOFs per node.
        self.M_global = np.zeros((6 * nn, 6 * nn))

        # Loop over all elements to assemble their contributions.
        for e, elem in enumerate(self.ListOfElements):

            # Assemble a shaft element.
            if elem.etype == "shaft":

                # Local stiffness matrix of the shaft element.
                K_e = elem.K

                # Local mass matrix of the shaft element.
                M_e = elem.M

                # Convert the first user node ID to an internal node index.
                p = self._get_internal_id(elem.n1)

                # Convert the second user node ID to an internal node index.
                q = self._get_internal_id(elem.n2)

                # Build the list of global DOF indices for the two-node element.
                # Each node has 4 DOFs, so a shaft element has 8 DOFs.
                dof_indices = [
                    6 * p,
                    6 * p + 1,
                    6 * p + 2,
                    6 * p + 3,
                    6 * p + 4,
                    6 * p + 5,
                    6 * q,
                    6 * q + 1,
                    6 * q + 2,
                    6 * q + 3,
                    6 * q + 4,
                    6 * q + 5
                ]

                # Add the local stiffness matrix into the global stiffness matrix.
                self.K_global[np.ix_(dof_indices, dof_indices)] += K_e

                # Add the local mass matrix into the global mass matrix.
                self.M_global[np.ix_(dof_indices, dof_indices)] += M_e

            # Assemble a ball-bearing element.
            elif elem.etype == "ball-bearing":

                # Local stiffness matrix of the bearing element.
                K_e = elem.K

                # Convert the bearing node ID to an internal node index.
                p = self._get_internal_id(elem.n1)

                # Bearings currently act only on the first two DOFs of the node.
                dof_indices = [
                    6 * p,
                    6 * p + 1,
                ]

                # Add the bearing stiffness matrix into the global stiffness matrix.
                self.K_global[np.ix_(dof_indices, dof_indices)] += K_e

    def apply_boundary_conditions(self, fixed_dofs: list | None = None):
        """
        Apply boundary conditions to the global matrices.

        Parameters
        ----------
        fixed_dofs : list of int, optional
            Global DOF indices to constrain.

            If `fixed_dofs` is `None`, the method keeps the already stored
            constraints and only rebuilds the reduced matrices.

        Effects
        -------
        This method updates:

        fixed_dofs
            Sorted list of constrained DOFs.

        free_dofs
            Sorted list of unconstrained DOFs.

        K_ff
            Stiffness matrix reduced to free DOFs.

        M_ff
            Mass matrix reduced to free DOFs.

        Notes
        -----
        Boundary conditions are applied by removing constrained rows and
        columns from the global matrices.
        """

        # Create a set containing all global DOF indices.
        all_dofs = set(range(self.K_global.shape[0]))

        # If no new constraints are provided, reuse the existing constraints.
        if fixed_dofs is None:

            # Remove duplicates and sort the existing fixed DOFs.
            self.fixed_dofs = sorted(set(self.fixed_dofs))

        # If new constraints are provided, merge them with the existing ones.
        else:

            # Add the new fixed DOFs, remove duplicates, and sort the result.
            self.fixed_dofs = sorted(set(self.fixed_dofs) | set(fixed_dofs))

        # Free DOFs are all DOFs that are not fixed.
        self.free_dofs = sorted(all_dofs - set(self.fixed_dofs))

        # Extract the reduced stiffness matrix using only free DOFs.
        self.K_ff = self.K_global[np.ix_(self.free_dofs, self.free_dofs)]

        # Extract the reduced mass matrix using only free DOFs.
        self.M_ff = self.M_global[np.ix_(self.free_dofs, self.free_dofs)]

    def add_rotor(self, rotor):
        """
        Add the elements and mesh of a rotor to the solver.

        Parameters
        ----------
        rotor : Rotor
            Rotor object containing elements and mesh data.

        Effects
        -------
        The rotor elements are appended to `ListOfElements`.

        The rotor mesh is inserted into or assigned to `self.meshgrid`.

        Notes
        -----
        This method does not assemble matrices immediately.
        Call `assemble_global_matrices()` after adding the rotor.
        """

        # Add all rotor elements.
        self.ListOfElements.extend(rotor.ListOfElements)

        # If no mesh exists yet, copy the rotor mesh and exit.
        if self.meshgrid is None:
            self.meshgrid = rotor.meshgrid
            return

        # Check for duplicated user node IDs.
        duplicated_ids = np.intersect1d(self.meshgrid.Nodes['id'], rotor.meshgrid.Nodes['id'])

        # Stop if duplicated node IDs are found.
        if duplicated_ids.size > 0:
            raise ValueError(f"Duplicated node IDs found: {duplicated_ids.tolist()}")

        # Merge the two structured arrays.
        self.meshgrid.Nodes = np.concatenate((self.meshgrid.Nodes, rotor.meshgrid))

            # Sort the merged mesh by node ID.
        self.meshgrid.Nodes.sort(order='id')

    def constrain_nodes(self, constrained_nodes, constrained_dofs_per_node=None):
        """
        Constrain selected DOFs on selected nodes.

        Parameters
        ----------
        constrained_nodes : sequence of int
            Internal node indices to constrain.

        constrained_dofs_per_node : sequence of int, optional
            Local DOF indices to constrain at each selected node.

            If not provided, all 4 DOFs of each node are constrained.

        Notes
        -----
        This method converts node-based constraints into global DOF indices.

        For a node `n`, the global DOF index is:

            4*n + local_dof
        """

        # If no specific DOFs are provided, constrain all DOFs at each node.
        if constrained_dofs_per_node is None:

            # Default: fix all 4 DOFs of each selected node.
            constrained_dofs_per_node = [0, 1, 2, 3, 4, 5]

        # List that will store the global DOF indices to constrain.
        constrained_dofs = []

        # Loop over all selected nodes.
        for node in constrained_nodes:

            # Convert local DOF indices into global DOF indices.
            constrained_dofs.extend([
                6 * node + dof
                for dof in constrained_dofs_per_node
            ])

        # Apply the computed global constraints.
        self.apply_boundary_conditions(constrained_dofs)

    def solve_eigenproblem(self):
        """
        Solve the generalized eigenvalue problem.

        The solved problem is:

            K_ff x = lambda M_ff x

        where:

            K_ff
                Reduced stiffness matrix.

            M_ff
                Reduced mass matrix.

            lambda
                Eigenvalue.

            x
                Eigenvector.

        Returns
        -------
        eigenvalues : ndarray
            Eigenvalues of the reduced system.

        eigenvectors : ndarray
            Eigenvectors of the reduced system.

        Notes
        -----
        This method calls `apply_boundary_conditions()` before solving.
        This ensures that reduced matrices exist even when the user did not
        explicitly constrain any node.
        """

        # Rebuild reduced matrices before solving.
        # This is useful when the model has no fixed nodes but is supported
        # only by bearing elements.
        self.apply_boundary_conditions()

        # Import the symmetric generalized eigenvalue solver from SciPy.
        from scipy.linalg import eigh

        # Solve the generalized eigenvalue problem.
        eigenvalues, eigenvectors = eigh(self.K_ff, self.M_ff)

        # Return eigenvalues and eigenvectors to the caller.
        return eigenvalues, eigenvectors
