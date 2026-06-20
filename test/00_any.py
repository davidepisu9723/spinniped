
import numpy as np
import sys

sys.path.append('./')  # Add parent directory to path to import spinniped
from spinniped import node


nn = 10
L = 1.0
meshgrid = node.Grid()

for i in range(nn):
    meshgrid.add_node(z = i*L/(nn-1))

print(meshgrid.Nodes)