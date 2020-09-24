import numpy as np                                                                                                                                           
import h5py
from mayavi import mlab

# Plot
def vis3d(vol):
    mlab.contour3d(vol,contours=8,opacity=.2)
    nx, ny, nz = vol.shape
    mlab.outline(extent=[0,nx,0,ny,0,nz])

    cam = mlab.gcf().scene.camera
    cam.position = [118.04698401457188, 154.5798768424432, 241.63946391675876]
    cam.focal_point = [50.0, 50.0, 25.0]
    cam.view_angle = 30.0
    cam.view_up = [-0.12850042543437323, 0.9082950116114877, -0.3981052782178007]
    cam.clipping_range = [135.68325605947476, 394.5042430005132]
    cam.compute_view_plane_normal()

    mlab.show()

f = h5py.File('mof.h5','r')
exact = np.array(f['visual']['vis01'])
exact = np.transpose(exact, (1, 2, 0))
f.close()
vis3d(exact)

f = h5py.File('mof.h5','r')
exact = np.array(f['visual']['vis02'])
exact = np.transpose(exact, (1, 2, 0))
f.close()
vis3d(exact)

f = h5py.File('dtmof.h5','r')
exact = np.array(f['visual']['vis02'])
exact = np.transpose(exact, (1, 2, 0))
f.close()
vis3d(exact)

f = h5py.File('elvira.h5','r')
exact = np.array(f['visual']['vis02'])
exact = np.transpose(exact, (1, 2, 0))
f.close()
vis3d(exact)
