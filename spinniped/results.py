import numpy as np
import matplotlib.pyplot as plt


class Results:
    def __init_(self):
        self.eigenvalues = None
        self.eigenvectors = None
        self.free_dofs = None
        self.nn = None

    def plot_mode_shapes(self, eigenvectors, eigenvalues, free_dofs, nn, list_of_modes=None):
        # Global dof array
        u = np.zeros(4*nn)

        plt.figure(figsize=(10, 6))
        if list_of_modes is None:
            modes_to_plot = range(min(5, eigenvectors.shape[1])) # Plot first 5 modes by default)
        else:
            modes_to_plot = list_of_modes

        for m in modes_to_plot:  # Plot x dir modes only
            # Get only x displacements for plotting
            u[free_dofs] = eigenvectors[:, m]
            mode_shape = u[::4]  # Extract x displacements (every 4th entry corresponds to x displacement)
            plt.plot(mode_shape, label=f'Mode {m+1} ({np.sqrt(eigenvalues[m])/(2*np.pi):.2f} Hz)')
        plt.xlabel('Node')
        plt.ylabel('Displacement (m)')
        plt.title('Mode Shapes of the Shaft')
        plt.legend()
        plt.grid()
        plt.show()