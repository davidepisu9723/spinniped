import numpy as np

grid_dtype = np.dtype([('id',np.int64),('internal_id',np.int64),('coordinates',np.float64,(3,))])


DOF_X_BENDING = [0, 4, 6, 10]   # [x0, ty0, x1, ty1]
DOF_Y_BENDING = [1, 3, 7, 9]    # [y0, tx0, y1, tx1]
DOF_AXIAL     = [2, 8]          # [z0, z1]
DOF_TORSION   = [5, 11]         # [tz0, tz1]
