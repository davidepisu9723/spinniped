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
        """Create a shaft element and precompute properties.

        Parameters
        ----------
        n : int
            Element identifier.
        E : float
            Young's modulus (Pa).
        nu : float
            Poisson's ratio.
        rho : float
            Material density (kg/m^3).
        d : float
            Outer diameter (m).
        L : float
            Element length (m).
        n1 : int
            Node index at element start.
        n2 : int
            Node index at element end.

        Notes
        -----
        Computes cross-sectional area, second moment of area, shear
        correction factor and element mass/stiffness matrices.
        """
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
        """Stiffness matrix using Timoshenko beam theory.

        Returns
        -------
        K : ndarray, shape (8, 8)
            Element stiffness matrix ordered as
            [x1, y1, tx1, ty1, x2, y2, tx2, ty2].
        """


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
        """Stiffness matrix using Euler-Bernoulli beam theory.

        Returns
        -------
        K : ndarray, shape (8, 8)
            Element stiffness matrix ordered as
            [x1, y1, tx1, ty1, x2, y2, tx2, ty2].
        """


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
        """Select and return the element stiffness matrix.

        Returns
        -------
        K : ndarray, shape (8, 8)
            Element stiffness matrix (Timoshenko by default).

        Notes
        -----
        The ordering of DOFs is [x_1, y_1, tx_1, ty_1, x_2, y_2, tx_2, ty_2].
        """
        return self.K_timosh()
    
    def printK(self):
        """Print the element stiffness matrix to standard output.

        Notes
        -----
        Prints entries in scientific notation with three decimal
        places for each matrix entry.
        """
        K = self.K
        print("Stiffness matrix K:")
        for i in range(K.shape[0]):
            for j in range(K.shape[1]):
                print(f"{K[i, j]:.3e}", end=" ")
            print()

    def M_trans(self):
        """Translational consistent mass matrix for the element.

        Returns
        -------
        M_trans : ndarray, shape (8, 8)
            Consistent translational mass matrix ordered as
            [x1, y1, tx1, ty1, x2, y2, tx2, ty2].
        """
        # The order of the dofs is [x_1, y_1, tx_1, ty_1, x_2, y_2, tx_2, ty_2]

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
        """Rotational inertia contribution to the element mass matrix.

        Returns
        -------
        M_rot : ndarray, shape (8, 8)
            Rotational mass terms used in rotordynamics. Currently
            this contribution is computed but not added by default.
        """
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
        """Return the element mass matrix used in analyses.

        Returns
        -------
        M : ndarray, shape (8, 8)
            Element mass matrix. Currently only the translational
            contribution is returned (rotational inertia omitted).
        """
        return self.M_trans()# + self.M_rot()    
    
    def printM(self):
        """Print the element mass matrix to standard output.

        Notes
        -----
        Prints entries in scientific notation with three decimal
        places for each matrix entry.
        """
        M = self.M
        print("Mass matrix M:")
        for i in range(M.shape[0]):
            for j in range(M.shape[1]):
                print(f"{M[i, j]:.3e}", end=" ")
            print()