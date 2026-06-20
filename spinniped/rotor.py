import numpy as np
import matplotlib.pyplot as plt

# Import 3D plotting support for Matplotlib.
# This import may be required by some Matplotlib versions to enable 3D axes.
from mpl_toolkits.mplot3d import Axes3D

import sys

# Add the current project root to the Python import path.
# This allows importing local modules such as `spinniped`.
sys.path.append('./')

# Import custom structured dtypes used by the project.
from spinniped import dtypes


class Rotor:
    """
    Rotor assembly container.

    This class stores the rotor mesh, the rotor elements, and basic rotor-level
    properties such as total length, total mass, and rotation axis.

    It also provides simple 3D plotting utilities for visual inspection of the
    rotor geometry.

    Notes
    -----
    The rotor geometry is stored through a `Grid` object assigned to
    `self.meshgrid`.

    The mesh nodes are expected to be stored in:

        self.meshgrid.Nodes

    Each node is expected to contain at least:

        id
            User-defined node identifier.

        coordinates
            Node coordinates stored as [x, y, z].

    The plotting functions currently support:

        shaft
            Plotted as a solid cylinder between two nodes.

        ball-bearing
            Plotted as a hollow cylinder centered on one node.
    """

    def __init__(self, rotation_axis: list[float] | None = None):
        """
        Initialize an empty rotor.

        Parameters
        ----------
        rotation_axis : list of float, optional
            Direction of the rotor rotation axis.

            If not provided, the global Z axis is used:

                [0.0, 0.0, 1.0]

        Attributes
        ----------
        ListOfElements : list
            List containing all elements added to the rotor.

        meshgrid : Grid or None
            Mesh grid associated with the rotor.

        L : float
            Rotor length, computed later from the mesh.

        m : float
            Rotor mass. Currently initialized but not computed here.

        rotation_axis : ndarray
            Unit vector defining the rotor rotation axis.
        """

        # List containing shaft, bearing, disk, or other rotor elements.
        self.ListOfElements = []

        # Mesh grid associated with this rotor.
        # It is assigned later using `add_meshgrid()`.
        self.meshgrid = None

        # Rotor length.
        # It is computed later by `_get_length()`.
        self.L = 0.0

        # Rotor mass.
        # It is initialized here but not computed in this class.
        self.m = 0.0

        # If no rotation axis is provided, use the global Z axis.
        if rotation_axis is None:
            rotation_axis = np.array([0.0, 0.0, 1.0])

        # Convert the rotation axis to a NumPy array.
        self.rotation_axis = np.asarray(rotation_axis, dtype=float)

        # Normalize the rotation axis to unit length.
        self.rotation_axis = self.rotation_axis / np.linalg.norm(self.rotation_axis)

    def add_element(self, element):
        """
        Add one element to the rotor.

        Parameters
        ----------
        element : object
            Rotor element to append.

            The element is expected to define an `etype` attribute.

            Supported element types in the current plotting method are:

                "shaft"
                "ball-bearing"
        """

        # Append the element to the rotor element list.
        self.ListOfElements.append(element)

    def add_meshgrid(self, meshgrid):
        """
        Assign a mesh grid to the rotor.

        Parameters
        ----------
        meshgrid : Grid
            Grid object containing the rotor nodes.

        Notes
        -----
        The mesh grid is expected to expose the node structured array through:

            meshgrid.Nodes
        """

        # Store the mesh grid reference inside the rotor.
        self.meshgrid = meshgrid

    def _get_node(self, node_id):
        """
        Return a mesh node from its user-defined node ID.

        Parameters
        ----------
        node_id : int
            User-defined node ID.

        Returns
        -------
        node : ndarray record
            Structured-array row corresponding to the requested node.

        Raises
        ------
        IndexError
            Raised if the node ID is not found in the rotor mesh.

        Notes
        -----
        This method assumes that `self.meshgrid.Nodes['id']` is sorted.

        The search is performed with `np.searchsorted()`, so the node ID field
        must remain sorted for this method to work correctly.
        """

        # Check whether the requested node ID exists in the mesh.
        if np.isin(node_id, self.meshgrid.Nodes['id']):

            # Find the position of the node ID in the sorted node array.
            idx = np.searchsorted(self.meshgrid.Nodes['id'], node_id)

            # Return the structured-array row of the requested node.
            return self.meshgrid.Nodes[idx]

        # Stop if the requested node ID is not available.
        else:
            raise IndexError(f"Node id {node_id} not found in rotor")

    def _get_length(self):
        """
        Compute the rotor length from the mesh coordinates.

        The current implementation computes the length along the global Z axis.

        Effects
        -------
        Updates
            self.L

        Notes
        -----
        This method does not currently project the rotor on `rotation_axis`.
        It simply computes:

            max(z) - min(z)

        Therefore, it assumes that the rotor axial direction is aligned with
        the global Z axis.
        """

        # Initialize the minimum coordinate with very large values.
        coord_min = np.finfo(np.float64).max + np.zeros(3)

        # Initialize the maximum coordinate with very small values.
        coord_max = -np.finfo(np.float64).max + np.zeros(3)

        # Loop over all nodes in the mesh.
        for node in self.meshgrid.Nodes:

            # Update the minimum x, y, z coordinates.
            coord_min = np.minimum(coord_min, node['coordinates'])

            # Update the maximum x, y, z coordinates.
            coord_max = np.maximum(coord_max, node['coordinates'])

        # Compute rotor length along the global Z direction.
        L = np.abs(coord_max[2] - coord_min[2])

        # Store the computed rotor length.
        self.L = L

    def _create_shaft(self, p1, p2, radius, num_points=30):
        """
        Create plotting surfaces for a solid cylindrical shaft element.

        Parameters
        ----------
        p1 : array_like
            Start point of the shaft cylinder.

        p2 : array_like
            End point of the shaft cylinder.

        radius : float
            Shaft radius.

        num_points : int, optional
            Number of discretization points used around the circumference.

        Returns
        -------
        surfaces : list of tuple
            List of surfaces to plot.

            Each surface is returned as:

                (X, Y, Z)

            where `X`, `Y`, and `Z` are 2D arrays compatible with
            `ax.plot_surface()`.

        Notes
        -----
        The returned surfaces are:

            1. Cylinder lateral surface
            2. First circular end cap
            3. Second circular end cap
        """

        # Compute the shaft direction vector.
        d = p2 - p1

        # Compute the shaft length.
        length = np.linalg.norm(d)

        # If both points are equal, no cylinder can be created.
        if length == 0:
            return []

        # Convert the shaft direction to a unit vector.
        d = d / length

        # Build the first vector perpendicular to the shaft direction.
        if abs(d[2]) < 0.9:

            # Use a vector in the XY plane when the shaft is not almost vertical.
            perp1 = np.array([-d[1], d[0], 0])

        else:

            # Use the global X direction when the shaft is almost aligned with Z.
            perp1 = np.array([1, 0, 0])

        # Normalize the first perpendicular vector.
        perp1 = perp1 / np.linalg.norm(perp1)

        # Build the second perpendicular vector using a cross product.
        perp2 = np.cross(d, perp1)

        # Normalize the second perpendicular vector.
        perp2 = perp2 / np.linalg.norm(perp2)

        # List that will contain all generated surfaces.
        surfaces = []

        # Angular coordinate around the cylinder.
        theta = np.linspace(0, 2 * np.pi, num_points)

        # Axial coordinate along the cylinder length.
        z_line = np.linspace(0, length, num_points)

        # Create 2D parameter grids for angle and axial position.
        theta_grid, z_grid = np.meshgrid(theta, z_line)

        # Local X coordinate of a standard cylinder.
        X = radius * np.cos(theta_grid)

        # Local Y coordinate of a standard cylinder.
        Y = radius * np.sin(theta_grid)

        # Local Z coordinate of a standard cylinder.
        Z = z_grid

        # Allocate transformed global X coordinates.
        X_new = np.zeros_like(X)

        # Allocate transformed global Y coordinates.
        Y_new = np.zeros_like(Y)

        # Allocate transformed global Z coordinates.
        Z_new = np.zeros_like(Z)

        # Transform each local cylinder point into the global coordinate system.
        for i in range(X.shape[0]):
            for j in range(X.shape[1]):

                # Map the local cylinder point to the actual shaft position.
                point = p1 + X[i, j] * perp1 + Y[i, j] * perp2 + Z[i, j] * d

                # Store global X coordinate.
                X_new[i, j] = point[0]

                # Store global Y coordinate.
                Y_new[i, j] = point[1]

                # Store global Z coordinate.
                Z_new[i, j] = point[2]

        # Add the cylinder lateral surface.
        surfaces.append((X_new, Y_new, Z_new))

        # Angular coordinate for the circular end caps.
        theta_cap = np.linspace(0, 2 * np.pi, num_points)

        # Radial coordinate for the circular end caps.
        r_cap = np.linspace(0, radius, int(num_points / 2))

        # Create the polar parameter grid for the end caps.
        theta_cap_grid, r_cap_grid = np.meshgrid(theta_cap, r_cap)

        # Local X coordinate of the first cap.
        X_cap1 = r_cap_grid * np.cos(theta_cap_grid)

        # Local Y coordinate of the first cap.
        Y_cap1 = r_cap_grid * np.sin(theta_cap_grid)

        # Local Z coordinate of the first cap.
        Z_cap1 = np.zeros_like(X_cap1)

        # Allocate transformed global X coordinates for the first cap.
        X_cap1_new = np.zeros_like(X_cap1)

        # Allocate transformed global Y coordinates for the first cap.
        Y_cap1_new = np.zeros_like(Y_cap1)

        # Allocate transformed global Z coordinates for the first cap.
        Z_cap1_new = np.zeros_like(Z_cap1)

        # Transform the first end cap to point `p1`.
        for i in range(X_cap1.shape[0]):
            for j in range(X_cap1.shape[1]):

                # Map the local cap point to the global coordinate system.
                point = p1 + X_cap1[i, j] * perp1 + Y_cap1[i, j] * perp2

                # Store global X coordinate.
                X_cap1_new[i, j] = point[0]

                # Store global Y coordinate.
                Y_cap1_new[i, j] = point[1]

                # Store global Z coordinate.
                Z_cap1_new[i, j] = point[2]

        # Add the first end cap surface.
        surfaces.append((X_cap1_new, Y_cap1_new, Z_cap1_new))

        # Allocate transformed global X coordinates for the second cap.
        X_cap2_new = np.zeros_like(X_cap1)

        # Allocate transformed global Y coordinates for the second cap.
        Y_cap2_new = np.zeros_like(Y_cap1)

        # Allocate transformed global Z coordinates for the second cap.
        Z_cap2_new = np.zeros_like(Z_cap1)

        # Transform the second end cap to point `p2`.
        for i in range(X_cap1.shape[0]):
            for j in range(X_cap1.shape[1]):

                # Map the local cap point to the global coordinate system.
                point = p2 + X_cap1[i, j] * perp1 + Y_cap1[i, j] * perp2

                # Store global X coordinate.
                X_cap2_new[i, j] = point[0]

                # Store global Y coordinate.
                Y_cap2_new[i, j] = point[1]

                # Store global Z coordinate.
                Z_cap2_new[i, j] = point[2]

        # Add the second end cap surface.
        surfaces.append((X_cap2_new, Y_cap2_new, Z_cap2_new))

        # Return all shaft plotting surfaces.
        return surfaces

    def _create_bearing(
        self,
        center,
        inner_diameter,
        outer_diameter,
        length,
        axis=None,
        num_points=30,
    ):
        """
        Create plotting surfaces for a hollow cylindrical bearing.

        The bearing is represented as an annular cylinder centered on one
        rotor node.

        Parameters
        ----------
        center : array_like
            Center point of the bearing.

        inner_diameter : float
            Inner diameter of the bearing.

        outer_diameter : float
            Outer diameter of the bearing.

        length : float
            Bearing axial length.

        axis : array_like, optional
            Direction of the bearing axis.

            If not provided, `self.rotation_axis` is used.

        num_points : int, optional
            Number of discretization points around the circumference.

        Returns
        -------
        surfaces : list of tuple
            List of surfaces to plot.

            Each surface is returned as:

                (X, Y, Z)

            where `X`, `Y`, and `Z` are 2D arrays compatible with
            `ax.plot_surface()`.

        Notes
        -----
        The returned surfaces are:

            1. Outer cylindrical surface
            2. Inner cylindrical surface
            3. First annular end face
            4. Second annular end face
        """

        # Convert the bearing center to a floating-point NumPy array.
        center = np.asarray(center, dtype=float)

        # Convert inner diameter to inner radius.
        inner_radius = float(inner_diameter) / 2

        # Convert outer diameter to outer radius.
        outer_radius = float(outer_diameter) / 2

        # Convert bearing length to float.
        length = float(length)

        # If the inner radius is invalid, assign a default value.
        if inner_radius <= 0.0:
            inner_radius = self.L / 20

        # If the outer radius is invalid, assign a default value.
        if outer_radius <= inner_radius:
            outer_radius = inner_radius * 1.1

        # If the bearing length is invalid, assign a default value.
        if length <= 0.0:
            length = self.L / 50

        # If no bearing axis is provided, use the rotor rotation axis.
        if axis is None:
            axis = self.rotation_axis

        # If a bearing axis is provided, convert it to a NumPy array.
        else:
            axis = np.asarray(axis, dtype=float)

        # Compute the norm of the bearing axis.
        axis_norm = np.linalg.norm(axis)

        # Stop if the axis vector has zero length.
        if axis_norm == 0.0:
            raise ValueError("axis vector must be non-zero")

        # Normalize the bearing axis.
        axis = axis / axis_norm

        # Build the first unit vector perpendicular to the bearing axis.
        if abs(axis[2]) < 0.9:

            # Use a vector in the XY plane when the axis is not almost vertical.
            perp1 = np.array([-axis[1], axis[0], 0.0])

        else:

            # Use the global X direction when the axis is almost aligned with Z.
            perp1 = np.array([1.0, 0.0, 0.0])

        # Normalize the first perpendicular vector.
        perp1 = perp1 / np.linalg.norm(perp1)

        # Build the second perpendicular vector.
        perp2 = np.cross(axis, perp1)

        # Normalize the second perpendicular vector.
        perp2 = perp2 / np.linalg.norm(perp2)

        # List that will contain all generated bearing surfaces.
        surfaces = []

        # Half of the bearing axial length.
        half_length = 0.5 * length

        # Number of points along the bearing axis.
        axial_points = max(2, int(num_points / 2))

        # Number of points along the radial direction of the annular faces.
        radial_points = max(2, int(num_points / 4))

        def transform(radius_grid, theta_grid, axial_grid):
            """
            Transform local cylindrical coordinates into global coordinates.

            Parameters
            ----------
            radius_grid : ndarray
                Radius values.

            theta_grid : ndarray
                Angular values.

            axial_grid : ndarray
                Axial values measured along the bearing axis.

            Returns
            -------
            X, Y, Z : ndarray
                Global coordinate arrays.
            """

            # Local coordinate along the first perpendicular direction.
            local_1 = radius_grid * np.cos(theta_grid)

            # Local coordinate along the second perpendicular direction.
            local_2 = radius_grid * np.sin(theta_grid)

            # Convert local cylindrical coordinates to global Cartesian points.
            points = (
                center
                + local_1[..., None] * perp1
                + local_2[..., None] * perp2
                + axial_grid[..., None] * axis
            )

            # Return global X, Y, and Z coordinate arrays.
            return points[..., 0], points[..., 1], points[..., 2]

        # Angular coordinate around the bearing.
        theta = np.linspace(0.0, 2.0 * np.pi, num_points)

        # Axial coordinate from one bearing end to the other.
        axial = np.linspace(-half_length, half_length, axial_points)

        # Create angular and axial grids for the cylindrical surfaces.
        theta_grid, axial_grid = np.meshgrid(theta, axial)

        # Create a constant-radius grid for the outer surface.
        radius_grid = np.full_like(theta_grid, outer_radius)

        # Add the outer cylindrical surface.
        surfaces.append(transform(radius_grid, theta_grid, axial_grid))

        # Create a constant-radius grid for the inner surface.
        radius_grid = np.full_like(theta_grid, inner_radius)

        # Add the inner cylindrical surface.
        surfaces.append(transform(radius_grid, theta_grid, axial_grid))

        # Radial coordinate for annular end faces.
        radii = np.linspace(inner_radius, outer_radius, radial_points)

        # Create angular and radial grids for the annular faces.
        theta_cap_grid, radius_cap_grid = np.meshgrid(theta, radii)

        # Create the first end-face axial grid.
        axial_cap_grid = np.full_like(theta_cap_grid, -half_length)

        # Add the first annular end face.
        surfaces.append(transform(radius_cap_grid, theta_cap_grid, axial_cap_grid))

        # Create the second end-face axial grid.
        axial_cap_grid = np.full_like(theta_cap_grid, half_length)

        # Add the second annular end face.
        surfaces.append(transform(radius_cap_grid, theta_cap_grid, axial_cap_grid))

        # Return all bearing plotting surfaces.
        return surfaces

    def plot(self):
        """
        Plot the rotor assembly in 3D.

        Shaft elements are drawn as solid cylinders.
        Ball bearings are drawn as hollow cylinders.
        Nodes are drawn as points.

        Raises
        ------
        ValueError
            Raised if the computed rotor length is zero or negative.

        Notes
        -----
        The plotting method currently uses Matplotlib 3D plotting.

        The rotor length check depends on `_get_length()`, which currently
        computes length along the global Z axis.
        """

        # Create a new Matplotlib figure.
        fig = plt.figure(figsize=(12, 8))

        # Create a 3D axis.
        ax = fig.add_subplot(111, projection='3d')

        # Set the vertical axis used by the 3D view.
        ax.view_init(vertical_axis='x')

        # Compute the rotor length before plotting.
        self._get_length()

        # Stop if the rotor has zero or negative axial length.
        if self.L <= 0:
            raise ValueError("Rotor lenght can't be zero")

        # Loop over all rotor elements.
        for elem in self.ListOfElements:

            # Plot shaft elements.
            if elem.etype == "shaft":

                # Convert shaft diameter to radius.
                radius = elem.d / 2.0

                # Get the first node of the shaft element.
                n1 = self._get_node(elem.n1)

                # Get the second node of the shaft element.
                n2 = self._get_node(elem.n2)

                # Extract the coordinates of the first node.
                p1 = n1['coordinates']

                # Extract the coordinates of the second node.
                p2 = n2['coordinates']

                # Create the shaft cylinder surfaces.
                surfaces = self._create_shaft(p1, p2, radius, num_points=20)

                # Plot all shaft surfaces.
                for X, Y, Z in surfaces:
                    ax.plot_surface(
                        X,
                        Y,
                        Z,
                        alpha=0.8,
                        color='steelblue',
                        edgecolor='none',
                    )

            # Plot ball-bearing elements.
            elif elem.etype == "ball-bearing":

                # Get the bearing node.
                n1 = self._get_node(elem.n1)

                # Create the bearing hollow-cylinder surfaces.
                surfaces = self._create_bearing(
                    center=n1['coordinates'],
                    inner_diameter=elem.di,
                    outer_diameter=elem.do,
                    length=elem.L,
                )

                # Plot all bearing surfaces.
                for X, Y, Z in surfaces:
                    ax.plot_surface(
                        X,
                        Y,
                        Z,
                        alpha=0.8,
                        color='red',
                        edgecolor='none',
                    )

        # Plot all mesh nodes as red points.
        for node in self.meshgrid.Nodes:
            ax.scatter(*node['coordinates'], color='red', s=20)

        # Set the X-axis label.
        ax.set_xlabel('X')

        # Set the Y-axis label.
        ax.set_ylabel('Y')

        # Set the Z-axis label.
        ax.set_zlabel('Z')

        # Set the plot title.
        ax.set_title('3D Shaft Visualization')

        # Initialize the maximum coordinate range.
        max_range = 0

        # Find the largest absolute coordinate value among all nodes.
        for node in self.meshgrid.Nodes:
            max_range = max(max_range, np.max(np.abs(node['coordinates'])))

        # Optional manual axis limits.
        # These lines are currently disabled.
        # if max_range > 0:
        #     ax.set_xlim([-max_range * 1.1, max_range * 1.1])
        #     ax.set_ylim([-max_range * 1.1, max_range * 1.1])
        #     ax.set_zlim([-max_range * 1.1, max_range * 1.1])

        # Ask Matplotlib to use equal axis scaling.
        plt.axis('equal')

        # Improve plot spacing.
        plt.tight_layout()

        # Display the figure.
        plt.show()