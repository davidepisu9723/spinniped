
import numpy as np
import matplotlib.pyplot as plt
import sys

sys.path.append('./')  # Add parent directory to path to import spinniped
from spinniped import element
from spinniped import node

ne = 21
nn = ne + 1
L = 1.0
ele = []
meshgrid = node.Grid()

# Create nodes once and connect each element sequentially
for i in range(nn):
    meshgrid.add_node(z = i*L/ne)

for e in range(ne):
    ele.append(element.ShaftElement(
        eid=e,
        E=2.0e11,
        nu=0.3,
        rho=7850,
        d=0.01,
        L=L/ne,
        n1=meshgrid.Nodes['id'][e],
        n2=meshgrid.Nodes['id'][e+1]
    ))

# Create a rotor and add elements to it
from spinniped import rotor
r = rotor.Rotor()
for e in ele:
    r.add_element(e)    

r.add_meshgrid(meshgrid)

r.plot()

# Create a solver and add the rotor to it
from spinniped import solver
s = solver.Solver()
s.add_rotor(r)

s.assemble_global_matrices()
s.constrain_nodes([0, nn-1],constrained_dofs_per_node=[0,1,2,5])  # Constrain the first node
#s.constrain_nodes([nn-1],constrained_dofs_per_node=[0,1])  # Constrain last node

eigenvalues, eigenvectors = s.solve_eigenproblem()

for m in range(min(len(eigenvalues),18)):
    print(f"Frequency: {np.sqrt(eigenvalues[m])/(2*np.pi):.2f}")

from spinniped import results
res = results.Results()
res.plot_mode_shapes(eigenvectors, eigenvalues, s.free_dofs, nn)