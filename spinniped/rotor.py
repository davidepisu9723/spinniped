
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

class Rotor:
    def __init__(self):

        self.ListOfElements = []
        self.ListOfNodes = []  # Initialize nodes attribute
        self.node_map = {}

    def add_element(self, element):
        # Add element
        self.ListOfElements.append(element)
    
    def add_node(self, node):
        self.ListOfNodes.append(node)
        self.node_map[node.id] = node

    def _get_node(self, node_id):
        if node_id in self.node_map:
            return self.node_map[node_id]
        if 0 <= node_id < len(self.ListOfNodes):
            return self.ListOfNodes[node_id]
        raise IndexError(f"Node id {node_id} not found in rotor")

    def _create_cylinder(self, p1, p2, radius, num_points=30):
        """
        Create cylinder surfaces (side and end caps) between two points.
        
        Args:
            p1: Starting point (x1, y1, z1)
            p2: Ending point (x2, y2, z2)
            radius: Cylinder radius
            num_points: Number of points around the circumference
            
        Returns:
            List of (X, Y, Z) arrays for surface plotting
        """
        # Direction vector
        d = p2 - p1
        length = np.linalg.norm(d)
        
        if length == 0:
            return []
            
        # Normalize direction
        d = d / length
        
        # Create two perpendicular vectors
        if abs(d[2]) < 0.9:
            perp1 = np.array([-d[1], d[0], 0])
        else:
            perp1 = np.array([1, 0, 0])
        perp1 = perp1 / np.linalg.norm(perp1)
        perp2 = np.cross(d, perp1)
        perp2 = perp2 / np.linalg.norm(perp2)
        
        surfaces = []
        
        # Create cylinder side surface
        theta = np.linspace(0, 2*np.pi, num_points)
        z_line = np.linspace(0, length, num_points)
        
        theta_grid, z_grid = np.meshgrid(theta, z_line)
        
        # Cylinder in standard position
        X = radius * np.cos(theta_grid)
        Y = radius * np.sin(theta_grid)
        Z = z_grid
        
        # Transform to actual position
        X_new = np.zeros_like(X)
        Y_new = np.zeros_like(Y)
        Z_new = np.zeros_like(Z)
        
        for i in range(X.shape[0]):
            for j in range(X.shape[1]):
                point = p1 + X[i, j] * perp1 + Y[i, j] * perp2 + Z[i, j] * d
                X_new[i, j] = point[0]
                Y_new[i, j] = point[1]
                Z_new[i, j] = point[2]
        
        surfaces.append((X_new, Y_new, Z_new))
        
        # Create end cap at p1
        theta_cap = np.linspace(0, 2*np.pi, num_points)
        r_cap = np.linspace(0, radius, int(num_points/2))
        theta_cap_grid, r_cap_grid = np.meshgrid(theta_cap, r_cap)
        
        X_cap1 = r_cap_grid * np.cos(theta_cap_grid)
        Y_cap1 = r_cap_grid * np.sin(theta_cap_grid)
        Z_cap1 = np.zeros_like(X_cap1)
        
        X_cap1_new = np.zeros_like(X_cap1)
        Y_cap1_new = np.zeros_like(Y_cap1)
        Z_cap1_new = np.zeros_like(Z_cap1)
        
        for i in range(X_cap1.shape[0]):
            for j in range(X_cap1.shape[1]):
                point = p1 + X_cap1[i, j] * perp1 + Y_cap1[i, j] * perp2
                X_cap1_new[i, j] = point[0]
                Y_cap1_new[i, j] = point[1]
                Z_cap1_new[i, j] = point[2]
        
        surfaces.append((X_cap1_new, Y_cap1_new, Z_cap1_new))
        
        # Create end cap at p2
        X_cap2_new = np.zeros_like(X_cap1)
        Y_cap2_new = np.zeros_like(Y_cap1)
        Z_cap2_new = np.zeros_like(Z_cap1)
        
        for i in range(X_cap1.shape[0]):
            for j in range(X_cap1.shape[1]):
                point = p2 + X_cap1[i, j] * perp1 + Y_cap1[i, j] * perp2
                X_cap2_new[i, j] = point[0]
                Y_cap2_new[i, j] = point[1]
                Z_cap2_new[i, j] = point[2]
        
        surfaces.append((X_cap2_new, Y_cap2_new, Z_cap2_new))
        
        return surfaces

    def plot(self):
        """
        3D plot of the rotor geometry with solid cylinders representing each shaft element.
        """
        fig = plt.figure(figsize=(12, 8))
        ax = fig.add_subplot(111, projection='3d')
        ax.view_init(vertical_axis='x')
        
        # Plot each element as a solid cylinder
        for e in self.ListOfElements:
            # Each shaft element is represented by a cylinder between its two nodes
            radius = e.d / 2.0  # radius
            
            # Get coordinates of node 1 and node 2
            n1 = self._get_node(e.n1)
            n2 = self._get_node(e.n2)
            p1 = n1.coordinates
            p2 = n2.coordinates
            
            # Create cylinder surfaces (side and end caps)
            surfaces = self._create_cylinder(p1, p2, radius, num_points=20)
            
            # Plot all surfaces (side and end caps)
            for X, Y, Z in surfaces:
                ax.plot_surface(X, Y, Z, alpha=0.8, color='steelblue', edgecolor='none')
        
        # Plot nodes as points
        for node in self.ListOfNodes:
            ax.scatter(*node.coordinates, color='red', s=20)
        
        # Set labels and title
        ax.set_xlabel('X')
        ax.set_ylabel('Y')
        ax.set_zlabel('Z')
        ax.set_title('3D Shaft Visualization')
        
        # Make axes equal aspect ratio
        max_range = 0
        for node in self.ListOfNodes:
            max_range = max(max_range, np.max(np.abs(node.coordinates)))
        
        #if max_range > 0:
        #    ax.set_xlim([-max_range*1.1, max_range*1.1])
        #    ax.set_ylim([-max_range*1.1, max_range*1.1])
        #    ax.set_zlim([-max_range*1.1, max_range*1.1])
        
        plt.axis('equal')
        plt.tight_layout()
        plt.show()
