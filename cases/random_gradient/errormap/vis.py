from scipy.interpolate import griddata
import numpy as np
from mayavi import mlab

x,y,z = np.mgrid[0:1:0.02,0:1:0.02,0:1:0.02]

f=np.loadtxt('error.dat')
points = f[:,0:3]
values = f[:,3]
e = griddata(points,values,(x,y,z))

mlab.contour3d(x,y,z,e,opacity=0.2)
mlab.show()
