from scipy.interpolate import griddata  
import numpy as np  
from mayavi import mlab  

x,y,z = np.mgrid[0:0.5:0.01,0:0.5:0.01,0:0.5:0.01]
x2,y2,z2 = np.mgrid[0:1.0:0.01,0:1.0:0.01,0:1.0:0.01]


f=np.loadtxt('error_1.dat')  
points = f[:,0:3]  
values = f[:,3]  
e = griddata(points,values,(x,y,z))  

el = e[::-1,:,:]
e1 = np.concatenate((e,el),axis=0)
eu = e1[:,::-1,:]
e2 = np.concatenate((e1,eu),axis=1)
ef = e2[:,:,::-1]
e3 = np.concatenate((e2,ef),axis=2)

src = mlab.pipeline.scalar_field(x2,y2,z2,e3)
# mlab.pipeline.iso_surface(src, contours=[0.1,0.5,0.9,1.3], opacity=0.2)
mlab.pipeline.iso_surface(src, contours=[0.1,0.3,0.5], opacity=0.2)
mlab.pipeline.image_plane_widget(src,
                            _orientation='x_axes',
                            slice_index=50,
                        )
mlab.outline(extent=[0,1,0,1,0,1])
mlab.colorbar()
mlab.show()
