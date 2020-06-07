import numpy as np
from numpy import pi, sin, cos, mgrid
import math
import matplotlib.pyplot as plt
import mayavi
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm

nn = 100

x = np.linspace(0,1/3,nn+1)
y = np.linspace(0,1/3,nn+1)

[X,Y] = np.meshgrid(x,y)

V = np.ones((nn+1,nn+1))/2

for i in range(nn+1):
	for j in range(nn+1):
		V[    i,    j] = x[i]*y[j]*9/2


fig = plt.figure()

plt.contour(X,Y,V,np.logspace(0,8,21,base =2)/pow(2,8),colors='k')

# ax = fig.gca(projection='3d')
# surf = ax.plot_surface(X, Y, V, cmap=cm.coolwarm,
#                        linewidth=0, antialiased=False,alpha=0.5)
# ax.contour(X,Y,V,np.logspace(0,8,21,base =2)/pow(2,8),colors='k')
# ax.contour(X, Y, V, np.logspace(0,8,21,base =2)/pow(2,8), lw=3, cmap="autumn_r", linestyles="solid", offset=-1)

plt.show()