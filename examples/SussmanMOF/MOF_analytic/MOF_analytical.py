import numpy as np
from numpy import pi, sin, cos, mgrid
import math
import matplotlib.pyplot as plt
# import mayavi

V  = 1.0/16.0
ny = math.sqrt(2.0)/2.0
nx = math.sqrt(2.0)/2.0


gx1 = []
gy1 = []

for V in np.arange(0,0.5,0.05):
	gx = []
	gy = []
	for nx in np.arange(0,1,0.001):
		ny    = math.sqrt(1.00-pow(nx,2.0))
		alpha = math.sqrt(2.0*V*ny/nx)
		beta  = math.sqrt(2.0*V*nx/ny)
		if (alpha<=1 and beta<=1):
			gx.append(1.0/3.0*alpha)
			gy.append(1.0/3.0*beta)
	plt.plot(gx,gy,'k')

	gx = []
	gy = []
	for nx in np.arange(0,1,0.001):
		ny    = math.sqrt(1.0-pow(nx,2.0))
		alpha = math.sqrt(2.0*V*ny/nx)
		beta  = math.sqrt(2.0*V*nx/ny)
		if (alpha<=1 and beta<=1):
			gx.append(1.0- 1.0/3.0*alpha)
			gy.append(1.0- 1.0/3.0*beta)
	plt.plot(gx,gy,'k')

	gx = []
	gy = []
	for nx in np.arange(0.0,1,0.001):
		ny    = math.sqrt(1.0-pow(nx,2.0))
		alpha = math.sqrt(2.0*V*ny/nx)
		beta  = math.sqrt(2.0*V*nx/ny)
		if (alpha<=1 and beta<=1):
			gx.append(1.0/3.0*alpha)
			gy.append(1.0- 1.0/3.0*beta)
	plt.plot(gx,gy,'k')

	gx = []
	gy = []
	for nx in np.arange(0,1,0.001):
		ny    = math.sqrt(1.0-pow(nx,2.0))
		alpha = math.sqrt(2.0*V*ny/nx)
		beta  = math.sqrt(2.0*V*nx/ny)
		if (alpha<=1 and beta<=1):
			gx.append(1.0-1.0/3.0*alpha)
			gy.append(1.0/3.0*beta)
	plt.plot(gx,gy,'k')

	gx1 = []
	gy1 = []
	for nx in np.arange(0,1,0.001):
		ny    = math.sqrt(1.0-pow(nx,2.0))
		gamma = V + nx/2.0/ny
		sigma = V - nx/2.0/ny
		if ((0.5-1.0/12.0/V*nx/ny>=1.0/3.0) and (0.5-1.0/12.0/V*nx/ny<=2.0/3.0)): 
			gx1.append(0.5-1.0/12.0/V*nx/ny)
			gy1.append(V/2+1.0/24.0/V*pow(nx/ny,2.0))
	plt.plot(gx1,gy1,'r')

	gx1 = []
	gy1 = []
	for nx in np.arange(0,1,0.001):
		ny    = math.sqrt(1.0-pow(nx,2.0))
		gamma = V + nx/2/ny
		sigma = V - nx/2/ny
		if ((0.5-1.0/12.0/V*nx/ny>=1.0/3.0) and (0.5-1.0/12.0/V*nx/ny<=2.0/3.0)): 
			gx1.append(0.5+1.0/12.0/V*nx/ny)
			gy1.append(V/2.0+1.0/24.0/V*pow(nx/ny,2.0))
	plt.plot(gx1,gy1,'r')

	gx1 = []
	gy1 = []
	for nx in np.arange(0,1,0.001):
		ny    = math.sqrt(1.0-pow(nx,2.0))
		gamma = V + nx/2.0/ny
		sigma = V - nx/2.0/ny
		if ((0.5-1.0/12.0/V*nx/ny>=1.0/3.0) and (0.5-1.0/12.0/V*nx/ny<=2.0/3.0)): 
			gy1.append(0.5-1.0/12.0/V*nx/ny)
			gx1.append(V/2.0+1.0/24.0/V*pow(nx/ny,2.0))
	plt.plot(gx1,gy1,'r')

	gx1 = []
	gy1 = []
	for nx in np.arange(0,1,0.001):
		ny    = math.sqrt(1-pow(nx,2.0))
		gamma = V + nx/2.0/ny
		sigma = V - nx/2.0/ny
		if ((0.5-1.0/12.0/V*nx/ny>=1.0/3.0) and (0.5-1.0/12.0/V*nx/ny<=2.0/3.0)): 
			gy1.append(0.5+1.0/12.0/V*nx/ny)
			gx1.append(V/2.0+1.0/24.0/V*pow(nx/ny,2.0))
	plt.plot(gx1,gy1,'r')

	gx1 = []
	gy1 = []
	for nx in np.arange(0,1,0.001):
		ny    = math.sqrt(1.0-pow(nx,2.0))
		gamma = V + nx/2.0/ny
		sigma = V - nx/2.0/ny
		if ((0.5-1.0/12.0/V*nx/ny>=1.0/3.0) and (0.5-1.0/12.0/V*nx/ny<=2.0/3.0)): 
			gy1.append(0.5-1.0/12.0/V*nx/ny)
			gx1.append(1.0-(V/2.0+1.0/24.0/V*pow(nx/ny,2.0)))
	plt.plot(gx1,gy1,'r')

	gx1 = []
	gy1 = []
	for nx in np.arange(0,1,0.001):
		ny    = math.sqrt(1.0-pow(nx,2.0))
		gamma = V + nx/2.0/ny
		sigma = V - nx/2.0/ny
		if ((0.5-1.0/12.0/V*nx/ny>=1.0/3.0) and (0.5-1.0/12.0/V*nx/ny<=2.0/3.0)): 
			gy1.append(0.5+1.0/12.0/V*nx/ny)
			gx1.append(1-(V/2.0+1.0/24.0/V*pow(nx/ny,2.0)))
	plt.plot(gx1,gy1,'r')

	gx1 = []
	gy1 = []
	for nx in np.arange(0,1,0.001):
		ny    = math.sqrt(1.0-pow(nx,2.0))
		gamma = V + nx/2.0/ny
		sigma = V - nx/2.0/ny
		if ((0.5-1.0/12.0/V*nx/ny>=1.0/3.0) and (0.5-1.0/12.0/V*nx/ny<=2.0/3.0)): 
			gx1.append(0.5-1.0/12.0/V*nx/ny)
			gy1.append(1.0-(V/2.0+1.0/24.0/V*pow(nx/ny,2.0)))
	plt.plot(gx1,gy1,'r')

	gx1 = []
	gy1 = []
	for nx in np.arange(0,1,0.001):
		ny    = math.sqrt(1.0-pow(nx,2.0))
		gamma = V + nx/2.0/ny
		sigma = V - nx/2.0/ny
		if ((0.5-1.0/12.0/V*nx/ny>=1.0/3.0) and (0.5-1.0/12.0/V*nx/ny<=2.0/3.0)): 
			gx1.append(0.5+1.0/12.0/V*nx/ny)
			gy1.append(1.0-(V/2.0+1.0/24.0/V*pow(nx/ny,2)))
	plt.plot(gx1,gy1,'r')

plt.savefig("sketch.eps")