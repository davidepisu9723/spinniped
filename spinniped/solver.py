import numpy as np 

class Solver:
    
    def __init__(self):
        """Initialize a solver instance for rotor finite element analysis.

        Attributes
        ----------
        K_global, M_global : ndarray
            Global stiffness and mass matrices (assembled later).
        ListOfElements : list
            Elements to be assembled into global matrices.
        fixed_dofs, free_dofs : list
            Lists of constrained and free global DOF indices.
        K_ff, M_ff : ndarray
            Reduced matrices after applying boundary conditions.
        """
        self.K_global = None
        self.M_global = None
        self.ListOfElements = []
        self.ListOfNodes = []


        self.fixed_dofs = []
        self.free_dofs = []

        self.K_ff = None
        self.M_ff = None
    
    def assemble_global_matrices(self):
        """Assemble global stiffness and mass matrices from elements.

        Notes
        -----
        Assumes each element provides `.K` and `.M` attributes for
        its element matrices and that element nodes are numbered
        consecutively. Global DOF ordering uses 4 DOFs per node.
        """
        # Now assemble the global stiffness and mass matrices
        ne = len(self.ListOfElements)
        nn = len(self.ListOfNodes)
        self.K_global = np.zeros((4*nn, 4*nn))
        self.M_global = np.zeros((4*nn, 4*nn))

        for e,elem in enumerate(self.ListOfElements):

            if elem.etype == "shaft":   
                K_e = elem.K
                M_e = elem.M
                n1 = elem.n1
                n2 = elem.n2
                dof_indices = [4*n1, 4*n1+1, 4*n1+2, 4*n1+3, 4*n2, 4*n2+1, 4*n2+2, 4*n2+3]
                self.K_global[np.ix_(dof_indices, dof_indices)] += K_e
                self.M_global[np.ix_(dof_indices, dof_indices)] += M_e
            elif elem.etype == "ball-bearing":
                K_e = elem.K
                n1 = elem.n1
                dof_indices = [4*n1, 4*n1+1]
                self.K_global[np.ix_(dof_indices, dof_indices)] += K_e



    def apply_boundary_conditions(self, fixed_dofs: list | None = None):
        """Apply boundary conditions by specifying fixed DOF indices.

        Parameters
        ----------
        fixed_dofs : sequence of int
            Global DOF indices to be constrained.

        Effects
        -------
        Updates `fixed_dofs` and `free_dofs`, and computes reduced
        matrices `K_ff` and `M_ff` for the free DOFs.
        """

        # Initialize all DOFs
        all_dofs = set(range(self.K_global.shape[0]))

        # Apply boundary conditions (fix specified dofs)
        # Maintain unique, sorted lists
        if fixed_dofs is None:
            self.fixed_dofs = sorted(set(self.fixed_dofs))
        else:
            self.fixed_dofs = sorted(set(self.fixed_dofs) | set(fixed_dofs))
        
        self.free_dofs = sorted(all_dofs - set(self.fixed_dofs))

        self.K_ff = self.K_global[np.ix_(self.free_dofs, self.free_dofs)]   
        self.M_ff = self.M_global[np.ix_(self.free_dofs, self.free_dofs)]

    def add_rotor(self, rotor):
        """Ingest elements from a `Rotor` instance into the solver.

        Parameters
        ----------
        rotor : Rotor
            Rotor instance whose elements will be added to the solver.
        """
        # Add rotor's elements to the solver's list of elements
        self.ListOfElements.extend(rotor.ListOfElements)
        self.ListOfNodes.extend(rotor.ListOfNodes)

        

    def constrain_nodes(self, constrained_nodes, constrained_dofs_per_node=None):
        """Convenience to constrain whole nodes by specifying node ids.

        Parameters
        ----------
        constrained_nodes : sequence of int
            Node indices to constrain.
        constrained_dofs_per_node : sequence of int, optional
            DOF indices per node to fix (default all 4 DOFs: [0,1,2,3]).
        """
        # Constrain specified nodes by fixing their degrees of freedom
        if constrained_dofs_per_node is None:
            constrained_dofs_per_node = [0, 1, 2, 3]  # Fix all dofs by default

        constrained_dofs = []
        for node in constrained_nodes:
            constrained_dofs.extend([4*node + dof for dof in constrained_dofs_per_node])
            
        self.apply_boundary_conditions(constrained_dofs)

    def solve_eigenproblem(self):
        """Solve the generalized eigenvalue problem K_ff x = lambda M_ff x.

        Returns
        -------
        eigenvalues : ndarray
            Array of eigenvalues.
        eigenvectors : ndarray
            Matrix whose columns are eigenvectors corresponding to
            `eigenvalues`.
        """
        # Before calling the solver, evalutate K_ff and M_ff in case the user defined no constrained node (i.e. he's just using bearings)
        self.apply_boundary_conditions()

        # Solve the eigenvalue problem  
        from scipy.linalg import eigh
        eigenvalues, eigenvectors = eigh(self.K_ff, self.M_ff)
        return eigenvalues, eigenvectors
    
