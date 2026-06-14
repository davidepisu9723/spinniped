import numpy as np 

class Solver:
    
    def __init__(self):
        self.K_global = None
        self.M_global = None
        self.ListOfElements = []

        self.fixed_dofs = []
        self.free_dofs = []

        self.K_ff = None
        self.M_ff = None
    
    def assemble_global_matrices(self):
        # Now assemble the global stiffness and mass matrices
        ne = len(self.ListOfElements)
        nn = ne + 1
        self.K_global = np.zeros((4*nn, 4*nn))
        self.M_global = np.zeros((4*nn, 4*nn))

        for e in range(ne):
            K_e = self.ListOfElements[e].K
            M_e = self.ListOfElements[e].M
            n1 = self.ListOfElements[e].n1
            n2 = self.ListOfElements[e].n2
            dof_indices = [4*n1, 4*n1+1, 4*n1+2, 4*n1+3, 4*n2, 4*n2+1, 4*n2+2, 4*n2+3]
            self.K_global[np.ix_(dof_indices, dof_indices)] += K_e
            self.M_global[np.ix_(dof_indices, dof_indices)] += M_e


    def apply_boundary_conditions(self, fixed_dofs):

        # Initialize all DOFs
        all_dofs = set(range(self.K_global.shape[0]))

        # Apply boundary conditions (fix specified dofs)
        # Maintain unique, sorted lists
        self.fixed_dofs = sorted(set(self.fixed_dofs) | set(fixed_dofs))
        self.free_dofs = sorted(all_dofs - set(self.fixed_dofs))
        
        self.K_ff = self.K_global[np.ix_(self.free_dofs, self.free_dofs)]   
        self.M_ff = self.M_global[np.ix_(self.free_dofs, self.free_dofs)]

    def add_rotor(self, rotor):
        # Add rotor's elements to the solver's list of elements
        self.ListOfElements.extend(rotor.ListOfElements)
        

    def constrain_nodes(self, constrained_nodes, constrained_dofs_per_node=None):
        # Constrain specified nodes by fixing their degrees of freedom
        if constrained_dofs_per_node is None:
            constrained_dofs_per_node = [0, 1, 2, 3]  # Fix all dofs by default

        constrained_dofs = []
        for node in constrained_nodes:
            constrained_dofs.extend([4*node + dof for dof in constrained_dofs_per_node])
            
        self.apply_boundary_conditions(constrained_dofs)

    def solve_eigenproblem(self):
        # Solve the eigenvalue problem  
        from scipy.linalg import eigh
        eigenvalues, eigenvectors = eigh(self.K_ff, self.M_ff)
        return eigenvalues, eigenvectors
    
