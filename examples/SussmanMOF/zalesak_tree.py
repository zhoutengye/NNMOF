import numpy as np
import math
import matplotlib as mpl
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt
import matplotlib.colors as colors

class distance():

    def __init__(self,types,point,para):
        self.x    = point[0]
        self.y    = point[1]
        self.para = para
        if types == 'zalesak':
            self.status = self.zalesak()
        if types == 'circle':
            self.status = self.circle()

    def zalesak(self):
        x_center = self.para[0]
        y_center = self.para[1]
        radius   = self.para[2]
        notch_width  = self.para[3]
        notch_height = self.para[4]
        d1 = math.sqrt( ( self.x - self.para[0] )**2.0 + ( self.y - self.para[1] )**2.0 ) - self.para[2]
        d2 = - max( 0.5 -self.para[3]/2-self.x , self.x-0.5 - self.para[3]/2 )
        d3 = self.para[4] - self.y 
        return(- max(d1,min(d2,d3)))

    def circle(self):
        x_center = self.para[0]
        y_center = self.para[1]
        radius   = self.para[2]
        return(math.sqrt( ( self.x - self.para[1] )**2.0 + ( self.y - self.para[2] )**2.0 ) - self.para[3])


xc = np.arange(0+0.005,1-0.005,0.01)
yc = np.arange(0+0.005,1-0.005,0.01)

f = np.zeros((99,99))

global shape_type
global shape_para
shape_type = 'zalesak'
shape_para = [0.5,0.75,0.15,0.05,0.85]

Y, X = np.meshgrid(xc,yc)

print( distance(shape_type,[xc[10],yc[10]],shape_para).status)

# for i in range(0,99):
#     for j in range(0,99):
#         # f[i,j] = object_distance.zalesak(xc[i],yc[j])
#         a = object_distance('zalesak',xc[i],yc[j])
#         if (a.status>0):
#         # if (get_shape(object_distance.zalesak,xc[i],yc[j])>0):
#             f[i,j] = 1
#         else:
#             f[i,j] = 0


# fig = plt.figure()
# plt.contour(X, Y, f, 0)
# plt.show()