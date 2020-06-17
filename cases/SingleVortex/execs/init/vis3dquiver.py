import numpy as np
import h5py
from mayavi import mlab
f = h5py.File('visual.h5','r')
u = np.array(f['visual']['u'])
v = np.array(f['visual']['v'])
w = np.array(f['visual']['w'])
mlab.quiver3d(u, v, w)
nx, ny, nz = u.shape
f.close()
mlab.outline(extent=[0,nx,0,ny,0,nz])
mlab.show()
