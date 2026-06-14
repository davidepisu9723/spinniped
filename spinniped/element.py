import numpy as np

class ShaftElement:
    
    def __init__(self, 
                 n: int, 
                 E: float,
                 nu: float,
                 rho: float,
                 d: float,
                 L: float,
                 n1: int,
                 n2: int):
        
        self.n = n
        self.E = E
        self.nu = nu
        self.rho = rho
        self.d = d
        self.L = L
        self.slenderness = L/d

        if self.slenderness < 0.:
            raise ValueError(f"Warning: Element {n} has low slenderness ratio ({self.slenderness:.2f}). Consider using a different element type for better accuracy.")

        self.n1 = n1
        self.n2 = n2    

        self.A = np.pi / 4.0 * d**2
        self.I = np.pi / 64.0 * d**4
        self.G = E / (2.0 * (1.0 + nu))
        self.kappa = 6*(1+nu) / (7 + 6*nu)
        self.phi = 12.0 * E * self.I / (self.kappa * self.G * self.A * self.L**2)

        self.M = self.M()
        self.K = self.K()   
    
    def K_timosh(self):
        """Returns the stiffness matrix of the shaft element."""


        # Compute stiffness matrix for timoshenko beam theory
        # For now just use 8 degrees of freedom (4 per node) [x_1, y_1, tx_1, ty_1, x_2, y_2, tx_2, ty_2]
        K = np.zeros((8, 8))

        k1 = 12.0
        k2 = 6.0 * self.L
        k3 = (4.0 + self.phi) * self.L**2
        k4 = (2.0 - self.phi) * self.L**2

        K = self.E * self.I / (self.L**3 * (1.0 + self.phi)) * np.array(
            [
                [ k1, 0.0,  k2, 0.0, -k1, 0.0,  k2, 0.0],
                [0.0,  k1, 0.0,  k2, 0.0, -k1, 0.0,  k2],
                [ k2, 0.0,  k3, 0.0, -k2, 0.0,  k4, 0.0],
                [0.0,  k2, 0.0,  k3, 0.0, -k2, 0.0,  k4],
                [-k1, 0.0, -k2, 0.0,  k1, 0.0, -k2, 0.0],
                [0.0, -k1, 0.0, -k2, 0.0,  k1, 0.0, -k2],
                [ k2, 0.0,  k4, 0.0, -k2, 0.0,  k3, 0.0],
                [0.0,  k2, 0.0,  k4, 0.0, -k2, 0.0,  k3],
            ],
            dtype=float,
        )
        return K
    
    def K_euler(self):
        """Returns the stiffness matrix of the shaft element."""


        # Compute stiffness matrix for timoshenko beam theory
        # For now just use 8 degrees of freedom (4 per node) [x_1, y_1, tx_1, ty_1, x_2, y_2, tx_2, ty_2]
        K = np.zeros((8, 8))

        k1 = 12.0
        k2 = 6.0 * self.L
        k3 = 4.0 * self.L**2
        k4 = 2.0 * self.L**2

        K = self.E * self.I / self.L**3 * np.array(
            [
                [ k1, 0.0,  k2, 0.0, -k1, 0.0,  k2, 0.0],
                [0.0,  k1, 0.0,  k2, 0.0, -k1, 0.0,  k2],
                [ k2, 0.0,  k3, 0.0, -k2, 0.0,  k4, 0.0],
                [0.0,  k2, 0.0,  k3, 0.0, -k2, 0.0,  k4],
                [-k1, 0.0, -k2, 0.0,  k1, 0.0, -k2, 0.0],
                [0.0, -k1, 0.0, -k2, 0.0,  k1, 0.0, -k2],
                [ k2, 0.0,  k4, 0.0, -k2, 0.0,  k3, 0.0],
                [0.0,  k2, 0.0,  k4, 0.0, -k2, 0.0,  k3],
            ],
            dtype=float,
        )
        return K
    
    def K(self):
        """Returns the stiffness matrix of the shaft element."""
        #The order of the dofs is [x_1, y_1, tx_1, ty_1, x_2, y_2, tx_2, ty_2]
        return self.K_timosh()
    
    def printK(self):
        """Prints the stiffness matrix of the shaft element."""
        K = self.K
        print("Stiffness matrix K:")
        for i in range(K.shape[0]):
            for j in range(K.shape[1]):
                print(f"{K[i, j]:.3e}", end=" ")
            print()

    def M_trans(self):
        """Returns the mass matrix of the shaft element."""
        #The order of the dofs is [x_1, y_1, tx_1, ty_1, x_2, y_2, tx_2, ty_2]

        m1 = 312.0 + 588.0 * self.phi + 280.0 * self.phi**2
        m2 = (44.0 + 77.0 * self.phi + 35.0 * self.phi**2) * self.L
        m3 = 108.0 + 252.0 * self.phi + 140.0 * self.phi**2
        m4 = -(26.0 + 63.0 * self.phi + 35.0 * self.phi**2) * self.L
        m5 = (8.0 + 14.0 * self.phi + 7.0 * self.phi**2) * self.L**2
        m6 = -(6.0 + 14.0 * self.phi + 7.0 * self.phi**2) * self.L**2

        M_trans = self.rho * self.A * self.L / (840.0 * (1.0 + self.phi)**2) * np.array(
            [
                [ m1, 0.0,  m2, 0.0,  m3, 0.0,  m4, 0.0],
                [0.0,  m1, 0.0,  m2, 0.0,  m3, 0.0,  m4],
                [ m2, 0.0,  m5, 0.0, -m4, 0.0,  m6, 0.0],
                [0.0,  m2, 0.0,  m5, 0.0, -m4, 0.0,  m6],
                [ m3, 0.0, -m4, 0.0,  m1, 0.0, -m2, 0.0],
                [0.0,  m3, 0.0, -m4, 0.0,  m1, 0.0, -m2],
                [ m4, 0.0,  m6, 0.0, -m2, 0.0,  m5, 0.0],
                [0.0,  m4, 0.0,  m6, 0.0, -m2, 0.0,  m5],
            ],
            dtype=float,
        )

        return M_trans

    def M_rot(self):

        """Returns the mass matrix of the shaft element."""
        m7 = 36.0
        m8 = (3.0 - 15.0 * self.phi * self.phi) * self.L
        m9 = (4.0 + 5.0 * self.phi + 10.0 * self.phi**2) * self.L**2
        m10 = (-1.0 - 5.0 * self.phi + 5.0 * self.phi**2) * self.L**2

        M_rot = self.rho * self.I / (30.0 * self.L * (1.0 + self.phi)**2) * np.array(
        [
            [ m7, 0.0,  m8, 0.0, -m7, 0.0,  m8, 0.0],
            [0.0,  m7, 0.0,  m8, 0.0, -m7, 0.0,  m8],
            [ m8, 0.0,  m9, 0.0, -m8, 0.0,  m10, 0.0],
            [0.0,  m8, 0.0,  m9, 0.0, -m8, 0.0,  m10],
            [-m7, 0.0, -m8, 0.0,  m7, 0.0, -m8, 0.0],
            [0.0, -m7, 0.0, -m8, 0.0,  m7, 0.0, -m8],
            [ m8, 0.0,  m10, 0.0, -m8, 0.0,  m9, 0.0],
            [0.0,  m8, 0.0,  m10, 0.0, -m8, 0.0,  m9],
        ],
        dtype=float,
        )

        return M_rot

    def M(self):
        """Returns the mass matrix of the shaft element."""
        return self.M_trans()# + self.M_rot()    
    
    def printM(self):
        """Prints the mass matrix of the shaft element."""
        M = self.M
        print("Mass matrix M:")
        for i in range(M.shape[0]):
            for j in range(M.shape[1]):
                print(f"{M[i, j]:.3e}", end=" ")
            print()