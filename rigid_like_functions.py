import numpy as np
from rigid_link import Rigid_Link
r_b = np.array([0,0,1]).T
orientation = np.array([56,2,3,4]).T

rigid_link = Rigid_Link.rigid_link(2)

rigid_link.orientation[:,0] = orientation
print(rigid_link.orientation)
print(rigid_link.compute_Q_I(0))

print


#def compute_L_matrix(rigid_link)