import numpy as np
class Node:

    def __init__(self, id, x, y, z):
        self.id = id
        self.coordinates = np.array([x, y, z])  # Initialize coordinates for each node (x, y, z)