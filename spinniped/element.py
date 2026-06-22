import numpy as np
from spinniped import dtypes

class ShaftElement:
    
    def __init__(self, 
                 eid: int, 
                 E: float,
                 nu: float,
                 rho: float,
                 d: float,
                 L: float,
                 n1: int,
                 n2: int,
                 etype: str = "shaft"):
        """Create a shaft element and precompute properties.

        Parameters
        ----------
        eid : int
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
        self.eid = eid
        self.E = E
        self.nu = nu
        self.rho = rho
        self.d = d
        self.L = L
        self.slenderness = L/d
        self.etype = etype

        if self.slenderness < 0.:
            raise ValueError(f"Element {eid} has low slenderness ratio ({self.slenderness:.2f}). Consider using a different element type for better accuracy.")

        self.n1 = n1
        self.n2 = n2    

        self.A = np.pi / 4.0 * d**2
        self.I = np.pi / 64.0 * d**4
        self.J = np.pi / 32.0 * d**4

        self.G = E / (2.0 * (1.0 + nu))
        self.kappa = 6*(1+nu) / (7 + 6*nu)
        self.phi = 12.0 * E * self.I / (self.kappa * self.G * self.A * self.L**2)

        self.M = self._M()
        self.K = self._K()   
    
    def _K_timosh(self):
        """Stiffness matrix using Timoshenko beam theory.

        Returns
        -------
        K : ndarray, shape (12, 12)
            Element stiffness matrix ordered as
            [x0, y0, z0, tx0, ty0, tz0, x1, y1, z1, tx1, ty1, tz1].

        """
        E = self.E
        G = self.G
        A = self.A
        L = self.L
        phi = self.phi
        I = self.I
        J = self.J

        # Compute stiffness matrix for timoshenko beam theory
        K = np.zeros((12, 12), dtype=np.float64)

        # Axial consistent stiffness matrix
        # DOF order: [z0, z1]
        K_ax = E*A /L * np.array([
            [1.0, -1.0],
            [-1.0, 1.0]
        ], dtype=np.float64)

        # Torsional consistent stiffness matrix
        # DOF order: [tz0, tz1]
        K_tors =  (G * J / L) * np.array([
        [1.0, -1.0],
        [-1.0, 1.0]
        ], dtype=np.float64)

        # Terms of bending stiffness matrices
        k0 = E * I / (L**3 * (1.0 + phi))
        k1 = 12.0
        k2 = 6.0 * L
        k3 = (4.0 + phi) * L**2
        k4 = (2.0 - phi) * L**2

        # Bending stiffness matrix in x-z plane
        # DOF order: [x0, ty0, x1, ty1]
        K_bend_x =  k0 * np.array([
        [ k1,    k2,    -k1,    k2],
        [ k2,    k3,    -k2,    k4],
        [-k1,   -k2,     k1,   -k2],
        [ k2,    k4,    -k2,    k3],
        ], dtype=np.float64) 

        # Bending stiffness matrix in y-z plane
        # DOF order: [y0, tx0, y1, tx1]
        #
        # Sign convention is already corrected because:
        #
        #     dy/dz = -tx
        #
        K_bend_y = k0 * np.array([
        [ k1,   -k2,    -k1,    -k2],
        [-k2,    k3,     k2,     k4],
        [-k1,    k2,     k1,     k2],
        [-k2,    k4,     k2,     k3],
        ], dtype=np.float64) 

        # Insert matrices permuting dofs according to Spinniped's convention 
        K[np.ix_([2,8],[2,8])] = K_ax
        K[np.ix_([5,11],[5,11])] = K_tors
        K[np.ix_([0,4,6,10],[0,4,6,10])] = K_bend_x
        K[np.ix_([1,3,7,9],[1,3,7,9])] = K_bend_y

        return K
    
    
    def _K_euler(self):
        """Stiffness matrix using Euler-Bernoulli beam theory.

        Returns
        -------
        K : ndarray, shape (12, 12)
            Element stiffness matrix ordered as
            [x0, y0, z0, tx0, ty0, tz0, x1, y1, z1, tx1, ty1, tz1].

        """

        # Compute stiffness matrix for Euler-Bernoulli beam theory
        K = np.zeros((12, 12))

        E = self.E
        G = self.G
        A = self.A
        L = self.L
        I = self.I
        J = self.J

        # Compute stiffness matrix for euler-bernoulli beam theory
        K = np.zeros((12, 12), dtype=np.float64)

        # Axial consistent stiffness matrix
        # DOF order: [z0, z1]
        K_ax = E*A /L * np.array([
            [1.0, -1.0],
            [-1.0, 1.0]
        ], dtype=np.float64)

        # Torsional consistent stiffness matrix
        # DOF order: [tz0, tz1]        
        K_tors =  (G * J / L) * np.array([
        [1.0, -1.0],
        [-1.0, 1.0]
        ], dtype=np.float64)

        # Terms of bending stiffness matrices
        k0 = E * I / (L**3)
        k1 = 12.0
        k2 = 6.0 * L
        k3 = 4.0 * L**2
        k4 = 2.0 * L**2

        # Bending stiffness matrix in x-z plane
        # DOF order: [x0, ty0, x1, ty1]
        K_bend_x =  k0 * np.array([
        [ k1,    k2,    -k1,    k2],
        [ k2,    k3,    -k2,    k4],
        [-k1,   -k2,     k1,   -k2],
        [ k2,    k4,    -k2,    k3],
        ], dtype=np.float64) 

        # Bending stiffness matrix in y-z plane
        # DOF order: [y0, tx0, y1, tx1]
        #
        # Sign convention is already corrected because:
        #
        #     dy/dz = -tx
        #
        K_bend_y = k0 * np.array([
        [ k1,   -k2,    -k1,    -k2],
        [-k2,    k3,     k2,     k4],
        [-k1,    k2,     k1,     k2],
        [-k2,    k4,     k2,     k3],
        ], dtype=np.float64) 
    
        # Insert matrices permuting dofs according to Spinniped's convention 
        K[np.ix_([2,8],[2,8])] = K_ax
        K[np.ix_([5,11],[5,11])] = K_tors
        K[np.ix_([0,4,6,10],[0,4,6,10])] = K_bend_x
        K[np.ix_([1,3,7,9],[1,3,7,9])] = K_bend_y

        return K
    
    def _K(self, theory: int = 0):
        """
        Select and return the element stiffness matrix.

        Parameters
        ----------
        theory : int
            Flag that selects the assiomatic theroy for beam element:
            0 -> Timoshenko beam theory, includes shear deformation effects.
            1 -> Euler-Bernoulli beam theory, does NOT include shear deformation effects.

        Returns
        -------
        K : ndarray, shape (12, 12)
            Element stiffness matrix (Timoshenko by default).

        Notes
        -----
        The ordering of DOFs is 
        [x0, y0, z0, tx0, ty0, tz0, x1, y1, z1, tx1, ty1, tz1].

        """
        if theory == 0:
            return self._K_timosh()
        elif theory in [1,2]:
            return self._K_euler()
        else:
            raise NotImplementedError(f"Method {theory} is not implemented. Use 0 for Timoshenko, 1 for Euler-Bernoulli, 2 for Rayleigh")
    
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

    def _M_trans_timosh(self):
        """Translational mass matrix for the element.

        Returns
        -------
        M_trans : ndarray, shape (12, 12)
            Consistent translational mass matrix ordered as
            [x0, y0, z0, tx0, ty0, tz0, x1, y1, z1, tx1, ty1, tz1].
        """

        rho = self.rho 
        A = self.A 
        L = self.L
        J = self.J
        phi = self.phi

        M_trans = np.zeros((12, 12), dtype=np.float64)

        # Axial consistent mass matrix
        # DOF order: [z0, z1]
        M_ax = rho * A * L / 6 * np.array([
            [2,     1],
            [1,     2]
        ], dtype = np.float64)

        # Torsional consistent inertia matrix
        # DOF order: [tz0, tz1]
        #
        # This is kept inside M_trans because torsional modes need this
        # inertia even when bending rotary inertia is disabled.
        M_tors = rho * J * L / 6 * np.array([
            [2,     1],
            [1,     2]
        ], dtype = np.float64)

        # Terms of bending translational matrices
        m1 = 312.0 + 588.0 * phi + 280.0 * phi**2
        m2 = (44.0 + 77.0 * phi + 35.0 * phi**2) * L
        m3 = 108.0 + 252.0 * phi + 140.0 * phi**2
        m4 = -(26.0 + 63.0 * phi + 35.0 * phi**2) * L
        m5 = (8.0 + 14.0 * phi + 7.0 * phi**2) * L**2
        m6 = -(6.0 + 14.0 * phi + 7.0 * phi**2) * L**2

        # Bending translational inertia in x-z plane
        # DOF order: [y0, tx0, y1, tx1]
        M_bendx = rho * A * L / (840.0 * (1.0 + phi)**2) * np.array(
            [
                [m1,  m2,  m3,  m4],
                [m2,  m5, -m4,  m6],
                [m3, -m4,  m1, -m2],
                [m4,  m6, -m2,  m5],
            ],
            dtype=np.float64,
        )

        # Bending translational inertia in y-z plane
        # DOF order: [y0, tx0, y1, tx1]
        #
        # Sign convention is already corrected because:
        #
        #     dy/dz = -tx
        #
        M_bendy = rho * A * L / (840.0 * (1.0 + phi)**2) * np.array(
            [
                [ m1, -m2,  m3, -m4],
                [-m2,  m5,  m4,  m6],
                [ m3,  m4,  m1,  m2],
                [-m4,  m6,  m2,  m5],
            ],
            dtype=np.float64,
        )


        # Insert matrices permuting dofs according to Spinniped's convention
        M_trans[np.ix_([2,8],[2,8])] = M_ax        
        M_trans[np.ix_([5,11],[5,11])] = M_tors
        M_trans[np.ix_([0,4,6,10],[0,4,6,10])] = M_bendx
        M_trans[np.ix_([1,3,7,9],[1,3,7,9])] = M_bendy

        return M_trans

    def _M_trans_euler(self):
        """Translational mass matrix using Euler-Bernoulli beam theory.

        Returns
        -------
        M_trans : ndarray, shape (12, 12)
            Consistent translational mass matrix ordered as
            [x0, y0, z0, tx0, ty0, tz0, x1, y1, z1, tx1, ty1, tz1].

        Notes
        -----
        This matrix contains:
        - axial inertia;
        - torsional inertia;
        - bending translational inertia in the x-z plane;
        - bending translational inertia in the y-z plane.

        It does not contain bending rotary inertia. That contribution is
        computed separately by M_rot_euler().
        """

        rho = self.rho
        A = self.A
        L = self.L
        J = self.J

        M_trans = np.zeros((12, 12), dtype=np.float64)

        # Axial consistent mass matrix
        # DOF order: [z0, z1]
        M_ax = rho * A * L / 6.0 * np.array([
            [2.0, 1.0],
            [1.0, 2.0],
        ], dtype=np.float64)

        # Torsional consistent inertia matrix
        # DOF order: [tz0, tz1]
        #
        # This is kept inside M_trans because torsional modes need this
        # inertia even when bending rotary inertia is disabled.
        M_tors = rho * J * L / 6.0 * np.array([
            [2.0, 1.0],
            [1.0, 2.0],
        ], dtype=np.float64)

        # Euler-Bernoulli bending translational inertia
        m0 = rho * A * L / 420.0

        m1 = 156.0
        m2 = 22.0 * L
        m3 = 54.0
        m4 = -13.0 * L
        m5 = 4.0 * L**2
        m6 = -3.0 * L**2

        # Bending translational inertia in x-z plane
        # DOF order: [x0, ty0, x1, ty1]
        M_bendx = m0 * np.array([
            [m1,  m2,  m3,  m4],
            [m2,  m5, -m4,  m6],
            [m3, -m4,  m1, -m2],
            [m4,  m6, -m2,  m5],
        ], dtype=np.float64)

        # Bending translational inertia in y-z plane
        # DOF order: [y0, tx0, y1, tx1]
        #
        # Sign convention is already corrected because:
        #
        #     dy/dz = -tx
        #
        M_bendy = m0 * np.array([
            [ m1, -m2,  m3, -m4],
            [-m2,  m5,  m4,  m6],
            [ m3,  m4,  m1,  m2],
            [-m4,  m6,  m2,  m5],
        ], dtype=np.float64)

        # Insert matrices permuting dofs according to Spinniped's convention
        M_trans[np.ix_([2, 8], [2, 8])] = M_ax
        M_trans[np.ix_([5, 11], [5, 11])] = M_tors
        M_trans[np.ix_([0, 4, 6, 10], [0, 4, 6, 10])] = M_bendx
        M_trans[np.ix_([1, 3, 7, 9], [1, 3, 7, 9])] = M_bendy

        return M_trans

    def _M_rot_timosh(self):
        """Rotational inertia contribution to the element mass matrix.

        Returns
        -------
        M_rot : ndarray, shape (12, 12)
            Rotational mass terms used in rotordynamics. Currently
            this contribution is computed but not added by default.
        """

        M_rot = np.zeros((12,12),dtype = np.float64)

        rho = self.rho
        L = self.L
        I = self.I
        phi = self.phi

        m0 = rho * I / (30.0 * L * (1.0 + phi)**2) 
        m7 = 36.0
        m8 = (3.0 - 15.0 * phi) * L
        m9 = (4.0 + 5.0 * phi + 10.0 * phi**2) * L**2
        m10 = (-1.0 - 5.0 * phi + 5.0 * phi**2) * L**2

        # Rotary inertia contribution for x-z bending
        # DOF order: [x0, ty0, x1, ty1]
        M_rot_x = m0 * np.array(
            [
                [ m7,  m8, -m7,  m8],
                [ m8,  m9, -m8,  m10],
                [-m7, -m8,  m7, -m8],
                [ m8,  m10, -m8,  m9],
            ],
            dtype=np.float64,
        )

        # Rotary inertia contribution for y-z bending
        # DOF order: [y0, tx0, y1, tx1]
        #
        # Sign convention is already corrected because:
        #
        #     dy/dz = -tx
        #
        M_rot_y = m0 * np.array(
            [
                [ m7, -m8, -m7, -m8],
                [-m8,  m9,  m8,  m10],
                [-m7,  m8,  m7,  m8],
                [-m8,  m10, m8,  m9],
            ],
            dtype=np.float64,
        )

        # Insert matrices permuting dofs according to Spinniped's convention
        M_rot[np.ix_([0, 4, 6, 10],[0, 4, 6, 10])] = M_rot_x
        M_rot[np.ix_([1, 3, 7, 9],[1, 3, 7, 9])] = M_rot_y

        return M_rot
    
    def _M_rot_euler(self):
        """Rotary inertia contribution using Euler-Bernoulli kinematics.

        Returns
        -------
        M_rot : ndarray, shape (12, 12)
            Bending rotary inertia contribution ordered as
            [x0, y0, z0, tx0, ty0, tz0, x1, y1, z1, tx1, ty1, tz1].

        Notes
        -----
        Strict Euler-Bernoulli theory neglects shear deformation. If this
        matrix is added, the mass model includes rotary inertia and becomes
        closer to a Rayleigh-beam-type inertia model.

        This matrix does not contain torsional inertia. Torsional inertia is
        already included in M_trans_euler().
        """

        M_rot = np.zeros((12, 12), dtype=np.float64)

        rho = self.rho
        L = self.L
        I = self.I

        m0 = rho * I / (30.0 * L)

        m7 = 36.0
        m8 = 3.0 * L
        m9 = 4.0 * L**2
        m10 = -1.0 * L**2

        # Rotary inertia contribution for x-z bending
        # DOF order: [x0, ty0, x1, ty1]
        M_rot_x = m0 * np.array([
            [ m7,  m8, -m7,  m8],
            [ m8,  m9, -m8,  m10],
            [-m7, -m8,  m7, -m8],
            [ m8,  m10, -m8,  m9],
        ], dtype=np.float64)

        # Rotary inertia contribution for y-z bending
        # DOF order: [y0, tx0, y1, tx1]
        #
        # Sign convention is already corrected because:
        #
        #     dy/dz = -tx
        #
        M_rot_y = m0 * np.array([
            [ m7, -m8, -m7, -m8],
            [-m8,  m9,  m8,  m10],
            [-m7,  m8,  m7,  m8],
            [-m8,  m10, m8,  m9],
        ], dtype=np.float64)

        M_rot[np.ix_([0, 4, 6, 10], [0, 4, 6, 10])] = M_rot_x
        M_rot[np.ix_([1, 3, 7, 9], [1, 3, 7, 9])] = M_rot_y

        return M_rot

    def _M(self, theory: int = 0, rotary_inertia: bool = True):
        """Return the element mass matrix used in analyses.

        Parameters
        ----------
        theory : int
            Flag that selects the assiomatic theroy for beam element:
            0 -> Timoshenko beam theory, includes rotary inertia effects.
            1 -> Euler-Bernoulli beam theory, usually does NOT include rotary inertia effects, unless specified by user
        rotary_intertia : bool
            Flag that activates (or deactivates) rotary inertia effects (relevant for high-frequency bending modes).

        Returns
        -------
        M : ndarray, shape (12, 12)
            Element mass matrix. 
        """
        if theory == 0:
            M_trans = self._M_trans_timosh()
        elif theory == 1:
            M_trans = self._M_trans_euler()
        else:
            raise NotImplementedError(f"Method {theory} is not implemented. Use 0 for Timoshenko, 1 for Euler-Bernoulli")
        
        # If no inertial effect is required, return the current translational mass matrix
        if not rotary_inertia:
            return M_trans
        
        # Otherwise, compute the rotary inertia mass matrix and sum it to the M_trans       
        if theory == 0:
            M_rot = self._M_rot_timosh()
        elif theory == 1:
            M_rot = self._M_rot_euler()
        else:
            raise NotImplementedError(f"Method {theory} is not implemented. Use 0 for Timoshenko, 1 for Euler-Bernoulli")

        return M_trans + M_rot    
    
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



class BallBearingElement:

    def __init__(self, 
                eid: int, 
                n1: int,
                n2: int | None = None,
                kxx: float = 0.0,
                kyy: float = 0.0,
                kxy: float = 0.0,
                kyx: float = 0.0,
                cxx: float = 0.0,
                cyy: float = 0.0,
                cxy: float = 0.0,
                cyx: float = 0.0,
                di: float = 0.0,
                do: float = 0.0,
                L: float = 0.0,
                etype: str = "ball-bearing"
                ):
    
        self.eid = eid
        self.n1 = n1
        self.K = np.array([[kxx, kxy],[kyx, kyy]])
        self.C = np.array([[cxx, cxy],[cyx, cyy]])
        self.di = di
        self.do = do
        self.L = L
        self.etype = etype
    