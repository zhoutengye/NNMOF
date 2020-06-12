import numpy as np
import h5py
import matplotlib.pyplot as plt
f = h5py.File('visual.h5','r')
for key in f['visual']:
    vis = np.array(f['visual'][key][10])
    plt.contourf(vis)

plt.show()
# nx, ny, nz = vis.shape
f.close()
