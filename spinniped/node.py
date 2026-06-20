import numpy as np
from warnings import warn
from spinniped import dtypes


class Grid:
    """
    Mesh grid container based on a NumPy structured array.

    The grid stores nodes as records inside `self.Nodes`.

    Each node has:

        id
            User-defined node identifier.

        internal_id
            Zero-based internal node index assigned by the program.

        coordinates
            Node coordinates in 3D space, stored as [x, y, z].

    Notes
    -----
    User node IDs are used by the user when defining elements.

    Internal node IDs are used by the solver to assemble global matrices.
    They are assigned in insertion order and are not necessarily equal to
    the user node IDs.
    """

    def __init__(self):
        """
        Initialize an empty grid.

        Attributes
        ----------
        Nodes : ndarray
            Structured array containing all grid nodes.

        internal_id : int
            Counter used to assign a unique internal ID to each new node.
        """

        # Structured array used to store all grid nodes.
        # The fields are defined in `dtypes.grid_dtype`.
        self.Nodes = np.empty(0, dtype=dtypes.grid_dtype)

        # Counter used to assign internal node IDs.
        # This ID is used internally by the solver.
        self.internal_id = 0

    def add_node(self, id: int = 0, x: float = 0, y: float = 0, z: float = 0):
        """
        Add a node to the grid.

        Parameters
        ----------
        id : int, optional
            User-defined node ID.

            If `id <= 0`, a new node ID is automatically assigned.

            If the requested ID already exists, the node is automatically
            renumbered to the next available ID.

        x : float, optional
            Node x-coordinate.

        y : float, optional
            Node y-coordinate.

        z : float, optional
            Node z-coordinate.

        Returns
        -------
        new_row : ndarray
            Structured array record containing the inserted node.

        Notes
        -----
        Nodes are inserted so that `self.Nodes` remains sorted by user node ID.

        The internal node ID is assigned from `self.internal_id` before the
        counter is increased.
        """

        # Convert the user-provided ID to NumPy int64.
        nid = np.int64(id)

        # If no valid positive ID is provided, generate one automatically.
        if nid <= 0:

            # If the grid already contains nodes, use the next ID after the maximum.
            if self.Nodes.shape[0] > 0:

                # Assign the next available user node ID.
                nid = np.max(self.Nodes['id']) + 1

            # If the grid is empty, start numbering from 1.
            else:

                # Assign the first user node ID.
                nid = 1

        # Check whether the requested user node ID already exists.
        if nid in self.Nodes['id']:

            # If the ID already exists, assign the next ID after the current maximum.
            nid = np.max(self.Nodes['id']) + 1

            # Warn the user that the requested ID was already present.
            warn(
                f"Grid point with id {id} is already present in the mesh. "
                f"Renumbered to {nid}"
            )

        # Find the insertion index that keeps nodes sorted by user node ID.
        idx = np.searchsorted(self.Nodes['id'], nid)

        # Create the new structured-array row.
        # The field order must match `dtypes.grid_dtype`.
        new_row = np.array(
            (nid, self.internal_id, np.array([x, y, z])),
            dtype=dtypes.grid_dtype,
        )

        # Insert the new node in the sorted position.
        self.Nodes = np.insert(self.Nodes, idx, new_row)

        # Increase the internal ID counter for the next inserted node.
        self.internal_id += 1

        # Return the inserted node record.
        return new_row
