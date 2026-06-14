import numpy as np
class Node:

    def __init__(self, id, x, y, z):
        """Create a node with 3D coordinates.

        Parameters
        ----------
        id : int
            Unique node identifier.
        x, y, z : float
            Coordinates of the node in 3D space.
        """
        self.id = id
        self.coordinates = np.array([x, y, z])  # Initialize coordinates for each node (x, y, z)