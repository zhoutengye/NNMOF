import numpy as np
from numpy import pi, sin, cos, mgrid
import math
import matplotlib.pyplot as plt
import mayavi

V  = 1/16
ny = math.sqrt(2)/2
nx = math.sqrt(2)/2


gx1 = []
gy1 = []

for V in np.arange(0,0.5,0.05):
	gx = []
	gy = []
	for nx in np.arange(0,1,0.001):
		ny    = math.sqrt(1-pow(nx,2))
		alpha = math.sqrt(2*V*ny/nx)
		beta  = math.sqrt(2*V*nx/ny)
		if (alpha<=1 and beta<=1):
			gx.append(1/3*alpha)
			gy.append(1/3*beta)
	plt.plot(gx,gy,'k')

	gx = []
	gy = []
	for nx in np.arange(0,1,0.001):
		ny    = math.sqrt(1-pow(nx,2))
		alpha = math.sqrt(2*V*ny/nx)
		beta  = math.sqrt(2*V*nx/ny)
		if (alpha<=1 and beta<=1):
			gx.append(1- 1/3*alpha)
			gy.append(1- 1/3*beta)
	plt.plot(gx,gy,'k')

	gx = []
	gy = []
	for nx in np.arange(0,1,0.001):
		ny    = math.sqrt(1-pow(nx,2))
		alpha = math.sqrt(2*V*ny/nx)
		beta  = math.sqrt(2*V*nx/ny)
		if (alpha<=1 and beta<=1):
			gx.append(1/3*alpha)
			gy.append(1- 1/3*beta)
	plt.plot(gx,gy,'k')

	gx = []
	gy = []
	for nx in np.arange(0,1,0.001):
		ny    = math.sqrt(1-pow(nx,2))
		alpha = math.sqrt(2*V*ny/nx)
		beta  = math.sqrt(2*V*nx/ny)
		if (alpha<=1 and beta<=1):
			gx.append(1-1/3*alpha)
			gy.append(1/3*beta)
	plt.plot(gx,gy,'k')

	gx1 = []
	gy1 = []
	for nx in np.arange(0,1,0.001):
		ny    = math.sqrt(1-pow(nx,2))
		gamma = V + nx/2/ny
		sigma = V - nx/2/ny
		if ((0.5-1/12/V*nx/ny>=1/3) and (0.5-1/12/V*nx/ny<=2/3)): 
			gx1.append(0.5-1/12/V*nx/ny)
			gy1.append(V/2+1/24/V*pow(nx/ny,2))
	plt.plot(gx1,gy1,'r')

	gx1 = []
	gy1 = []
	for nx in np.arange(0,1,0.001):
		ny    = math.sqrt(1-pow(nx,2))
		gamma = V + nx/2/ny
		sigma = V - nx/2/ny
		if ((0.5-1/12/V*nx/ny>=1/3) and (0.5-1/12/V*nx/ny<=2/3)): 
			gx1.append(0.5+1/12/V*nx/ny)
			gy1.append(V/2+1/24/V*pow(nx/ny,2))
	plt.plot(gx1,gy1,'r')

	gx1 = []
	gy1 = []
	for nx in np.arange(0,1,0.001):
		ny    = math.sqrt(1-pow(nx,2))
		gamma = V + nx/2/ny
		sigma = V - nx/2/ny
		if ((0.5-1/12/V*nx/ny>=1/3) and (0.5-1/12/V*nx/ny<=2/3)): 
			gy1.append(0.5-1/12/V*nx/ny)
			gx1.append(V/2+1/24/V*pow(nx/ny,2))
	plt.plot(gx1,gy1,'r')

	gx1 = []
	gy1 = []
	for nx in np.arange(0,1,0.001):
		ny    = math.sqrt(1-pow(nx,2))
		gamma = V + nx/2/ny
		sigma = V - nx/2/ny
		if ((0.5-1/12/V*nx/ny>=1/3) and (0.5-1/12/V*nx/ny<=2/3)): 
			gy1.append(0.5+1/12/V*nx/ny)
			gx1.append(V/2+1/24/V*pow(nx/ny,2))
	plt.plot(gx1,gy1,'r')

	gx1 = []
	gy1 = []
	for nx in np.arange(0,1,0.001):
		ny    = math.sqrt(1-pow(nx,2))
		gamma = V + nx/2/ny
		sigma = V - nx/2/ny
		if ((0.5-1/12/V*nx/ny>=1/3) and (0.5-1/12/V*nx/ny<=2/3)): 
			gy1.append(0.5-1/12/V*nx/ny)
			gx1.append(1-(V/2+1/24/V*pow(nx/ny,2)))
	plt.plot(gx1,gy1,'r')

	gx1 = []
	gy1 = []
	for nx in np.arange(0,1,0.001):
		ny    = math.sqrt(1-pow(nx,2))
		gamma = V + nx/2/ny
		sigma = V - nx/2/ny
		if ((0.5-1/12/V*nx/ny>=1/3) and (0.5-1/12/V*nx/ny<=2/3)): 
			gy1.append(0.5+1/12/V*nx/ny)
			gx1.append(1-(V/2+1/24/V*pow(nx/ny,2)))
	plt.plot(gx1,gy1,'r')

	gx1 = []
	gy1 = []
	for nx in np.arange(0,1,0.001):
		ny    = math.sqrt(1-pow(nx,2))
		gamma = V + nx/2/ny
		sigma = V - nx/2/ny
		if ((0.5-1/12/V*nx/ny>=1/3) and (0.5-1/12/V*nx/ny<=2/3)): 
			gx1.append(0.5-1/12/V*nx/ny)
			gy1.append(1-(V/2+1/24/V*pow(nx/ny,2)))
	plt.plot(gx1,gy1,'r')

	gx1 = []
	gy1 = []
	for nx in np.arange(0,1,0.001):
		ny    = math.sqrt(1-pow(nx,2))
		gamma = V + nx/2/ny
		sigma = V - nx/2/ny
		if ((0.5-1/12/V*nx/ny>=1/3) and (0.5-1/12/V*nx/ny<=2/3)): 
			gx1.append(0.5+1/12/V*nx/ny)
			gy1.append(1-(V/2+1/24/V*pow(nx/ny,2)))
	plt.plot(gx1,gy1,'r')

plt.savefig("sketch.eps")