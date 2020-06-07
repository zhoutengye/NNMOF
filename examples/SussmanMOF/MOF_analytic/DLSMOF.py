import numpy as np
from numpy import pi, sin, cos, mgrid
import math
import matplotlib.pyplot as plt

# The input Volume and Centroid
V= 1/16.0
px = 1/6.0
py = 1/8.0

# Grid information
lx = 1.0
ly = 1.0
n = 100
dx = lx/float(n)
dy = ly/float(n) 
x  = np.arange(0,1.0/3.0+dx,dx)
y  = np.arange(0,1.0/3.0+dy,dy)
xc = np.arange(dx/2.0,1.0/3.0+dx/2.0,dy)
yc = np.arange(dx/2.0,1.0/3.0+dx/2.0,dy)
f  = np.zeros([n+1,n+1])
fx  = np.zeros([n+1,n+1])
fy  = np.zeros([n+1,n+1])
ffx  = np.zeros([n+1,n+1])
ffy  = np.zeros([n+1,n+1])
for i in range(len(xc)):
	for j in range(len(yc)):
		f[i,j] = x[i] * y[j] * 9.0 / 2.0;
		ffx[i,j] = yc[j] * 9.0 / 2.0;
		ffy[i,j] = xc[i] * 9.0 / 2.0;

fx[1:n-1,1:n-1]= (f[2:n,1:n-1] - f[0:n-2,1:n-1])/dx/2.0
fy[1:n-1,1:n-1]= (f[1:n-1,2:n] - f[1:n-1,0:n-2])/dy/2.0

# Find the initial cell location
for i in range(len(x)):
	if px < x[i]:
		iic = i-1
		x_sub = x[i-1]
		x_sup = x[i]
		break

for j in range(len(y)):
	if py < y[j]:
		jjc = j-1
		y_sub = y[j-1]
		y_sup = y[j]
		break
ffc = f[iic,jjc]
ff = f[iic-1:iic+2,jjc-1:jjc+2]


# Determine the marching direction
if ffc > V:
	sign = -1.0
else:
	sign = 1.0

xk = px
yk = py
xk1 = -1.0
xk2 = -1.0
iteraction = 0

## Find the closest point by improved Newton iteration
while (iteraction < 100 ):

	#  Bilinear interpolation of f, fx, fy
	iic = int(xk / dx)
	jjc = int(yk / dy)
	pxx = xk - x[iic]
	pyy = yk - y[jjc]

	f00 = f[iic  ,jjc]
	f01 = f[iic+1,jjc]
	f10 = f[iic  ,jjc+1]
	f11 = f[iic+1,jjc+1]
	A   = f00 + (f10-f00)*pxx/dx
	B   = f01 + (f11-f01)*pxx/dx
	fi  = A + (B-A)*pyy/dy

	fx00 = fx[iic  ,jjc]
	fx01 = fx[iic+1,jjc]
	fx10 = fx[iic  ,jjc+1]
	fx11 = fx[iic+1,jjc+1]
	Ax   = fx00 + (fx10-fx00)*pxx/dx
	Bx   = fx01 + (fx11-fx01)*pxx/dx
	fix  = Ax + (Bx-Ax)*pyy/dy

	fy00 = fy[iic  ,jjc]
	fy01 = fy[iic+1,jjc]
	fy10 = fy[iic  ,jjc+1]
	fy11 = fy[iic+1,jjc+1]
	Ay   = fy00 + (fy10-fy00)*pxx/dx
	By   = fy01 + (fy11-fy01)*pxx/dx
	fiy  = Ay + (By-Ay)*pyy/dy

	# fi  = xk * yk * 9.0 / 2.0;
	# fix = yk * 9.0 / 2.0;
	# fiy = xk * 9.0 / 2.0;
	
	# Newton iteration of D. Chopp(1999)
	dp2  = fix*fix + fiy*fiy
	dxdp = (px-xk)*fix  + (py-yk)*fiy

	sigma11 = -(fi-V)*fix/dp2
	sigma12 = -(fi-V)*fiy/dp2
	xk12 = xk + sigma11 
	yk12 = yk + sigma12
	sigma21 = px - xk - dxdp*fix/dp2
	sigma22 = py - yk - dxdp*fiy/dp2
	xk1 = xk12 + sigma21 
	yk1 = yk12 + sigma22

	xk = xk1
	yk = yk1

	print(xk1,yk1)
	iteraction = iteraction +1


# fiy = 

# fi = f00*(dx-pxx)*(dy-pyy) + f10*pxx*(dy-pyy) + f10*(dx-pxx)*pyy + f11*pxx*pyy
print(px,py)
# print(px / dx, py / dy)

# # Determine whether close enough 
# if (f[iic-1:iic+2,jjc-1:jjc+2] - V >0).all() == True:
# 	condition = 1
# else:
# 	if (f[iic-1:iic+2,jjc-1:jjc+2] - V <0).all() == True:
# 		condition = -1
# 	else:
# 		condition = 0



# # Marching until it fall with criteria
# if sign == 1.0:
# 	while (condition != 0):
# 		t = min(abs(dd / fx[iic,jjc]),abs(dd / fy[iic,jjc]))
# 		print(px,py)
# 		px = px + sign * t * fx[iic,jjc]
# 		py = py + sign * t * fy[iic,jjc]
# 		print(px,py)
# 		print(iic)
# 		if px > x[iic+1]:
# 			iic = iic+1
# 		if py > y[jjc+1]:
# 			jjc = jjc+1
# 		print(iic)

# 		if (f[iic-1:iic+2,jjc-1:jjc+2] - V >0).all() == True:
# 			condition = 1
# 		else:
# 			if (f[iic-1:iic+2,jjc-1:jjc+2] - V <0).all() == True:
# 				condition = -1
# 			else:
# 				condition = 0
# else:
# 	while (condition != 0):
# 		t = min(abs(dd / fx[iic,jjc]),abs(dd / fy[iic,jjc]))/10.0
# 		print(px,py)
# 		px = px + sign * t * fx[iic,jjc]
# 		py = py + sign * t * fy[iic,jjc]
# 		print(px,py)
# 		print(iic)
# 		print(jjc)
# 		if px < x[iic]:
# 			iic = iic-1
# 		if py < y[jjc]:
# 			jjc = jjc-1
# 		print(iic)

# 		if (f[iic-1:iic+2,jjc-1:jjc+2] - V >0).all() == True:
# 			condition = 1
# 		else:
# 			if (f[iic-1:iic+2,jjc-1:jjc+2] - V <0).all() == True:
# 				condition = -1
# 			else:
# 				condition = 0

# print(xc[iic],yc[jjc],f[iic,jjc])


# Find 