import numpy as np
from numpy import pi, sin, cos, mgrid
import math
import matplotlib.pyplot as plt
import mayavi
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm

nn = 300

x = np.linspace(0,1,nn+1)
y = np.linspace(0,1,nn+1)

[X,Y] = np.meshgrid(x,y)

V = np.ones((nn+1,nn+1))/2

for i in range(int(nn/3)+1):
	for j in range(int(nn/3)+1):
		V[    i,    j] = x[i]*y[j]*9/2
		V[i+int(nn/3*2),j+int(nn/3*2)] = (1/3-x[i])*(1/3-y[j])*9/2
		V[i+int(nn/3*2),       j] = (1/3-x[i])*y[j]*9/2
		V[    i,j+int(nn/3*2)] = x[i]*(1/3-y[j])*9/2

		V[i+int(nn/3),    j] = y[j] / (6*pow((0.5 - x[i+int(nn/3)]),2)+0.5)
		V[i,    j+int(nn/3)] = x[i] / (6*pow((0.5 - y[j+int(nn/3)]),2)+0.5)
		V[i+int(nn/3),j+int(nn/3*2)] = (1/3-y[j]) / (6*pow((0.5 - x[i+int(nn/3)]),2)+0.5)
		V[i+int(nn/3*2),j+int(nn/3)] = (1/3-x[i]) / (6*pow((0.5 - y[j+int(nn/3)]),2)+0.5)


fig = plt.figure()

ax = fig.gca(projection='3d')
surf = ax.plot_surface(X, Y, V, cmap=cm.coolwarm,
                       linewidth=0, antialiased=False,alpha=0.5)
ax.contour(X,Y,V,np.logspace(0,8,21,base =2)/pow(2,8),colors='k')
ax.contour(X, Y, V, np.logspace(0,8,21,base =2)/pow(2,8), lw=3, cmap="autumn_r", linestyles="solid", offset=-1)

plt.show()