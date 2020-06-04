import numpy as np
import matplotlib.pyplot as plt
import h5py
import json
import os
import mayavi
from mayavi import mlab

nx = 20
ny = 20
nz = 20
x = np.arange(20)
y = np.arange(20)
z = np.arange(20)
X,Y,Z = np.meshgrid(x,y,z)
phi = np.zeros([nx,ny,nz])

phi[5:10,5:10,5:10] = 1.0
phi[5:10,5:10,5:10] = 1.0

# mlab.contour3d(phi,contours=8,opacity=.2 )
mlab.pipeline.volume(mlab.pipeline.scalar_field(phi))

mlab.show()
