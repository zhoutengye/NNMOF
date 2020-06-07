import numpy as np

V = 1.0/6.0

xmin = V
xmax = 1.0/3.0
nn   = 100.0

# [x,y] = np.mgrid[xmin:xmax:(xmax-xmin)/nn,xmin:xmax:(xmax-xmin)/nn]
# z = V/9.0/x/y

# for i in range(100):
#     for j in range(100):
#         if z[i,j] > xmax:
#             z[i,j] = 0

dd = 0.05

x = np.linspace(0.0,1.0/3.0,dd)

for x1 in x:
    y = linsapce()


# View it.
from mayavi import mlab
import PyQt5
mlab.points3d(x, y, z)
mlab.show()

