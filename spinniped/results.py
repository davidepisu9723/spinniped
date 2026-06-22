import numpy as np
import matplotlib.pyplot as plt


class Results:
    def __init__(self):
        """Initialize results container.

        Notes
        -----
        This method appears to be a misspelled initializer (`__init_`)
        and may not be called when instantiating `Results`. It sets
        placeholders for eigen-data used by plotting routines.
        """
        self.eigenvalues = None
        self.eigenvectors = None
        self.free_dofs = None
        self.nn = None

    def plot_mode_shapes(self, eigenvectors, eigenvalues, free_dofs, nn, list_of_modes=None):
        """Plot mode shapes extracted from eigenvectors.

        Parameters
        ----------
        eigenvectors : ndarray
            Matrix of eigenvectors (columns are modes).
        eigenvalues : ndarray
            Array of eigenvalues corresponding to modes.
        free_dofs : array_like
            Indices of free degrees of freedom in the global system.
        nn : int
            Number of nodes in the model.
        list_of_modes : sequence, optional
            Specific mode indices to plot. If ``None``, the first
            up to five modes are plotted.

        Notes
        -----
        Only the x-direction displacements (every 4th DOF) are
        extracted and plotted for each selected mode.
        """
        # Global dof array
        u = np.zeros(6*nn)

        plt.figure(figsize=(10, 6))
        if list_of_modes is None:
            modes_to_plot = range(min(5, eigenvectors.shape[1])) # Plot first 5 modes by default)
        else:
            modes_to_plot = list_of_modes

        for m in modes_to_plot:  # Plot x dir modes only
            # Get only x displacements for plotting
            u[free_dofs] = eigenvectors[:, m]
            mode_shape = u[::6]  # Extract x displacements (every 4th entry corresponds to x displacement)
            plt.plot(mode_shape, label=f'Mode {m+1} ({np.sqrt(eigenvalues[m])/(2*np.pi):.2f} Hz)')
        plt.xlabel('Node')
        plt.ylabel('Displacement (m)')
        plt.title('Mode Shapes of the Shaft')
        plt.legend()
        plt.grid()
        plt.show()