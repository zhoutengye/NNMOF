import numpy as np
from numpy import pi, sin, cos, mgrid
import math
import matplotlib.pyplot as plt
# import mayavi

V  = 1/10.0
ny = math.sqrt(2.0)/2.0
nx = math.sqrt(2.0)/2.0


gx1 = []
gy1 = []

sx = []
sy = []

# for V in np.arange(0,0.5,0.05):
# 	gx = []
# 	gy = []
# 	for nx in np.arange(0,1,0.001):
# 		ny    = math.sqrt(1-pow(nx,2))
# 		alpha = math.sqrt(2*V*ny/nx)
# 		beta  = math.sqrt(2*V*nx/ny)
# 		if (alpha<=1 and beta<=1):
# 			gx.append(1/3*alpha)
# 			gy.append(1/3*beta)
# 	plt.plot(gx,gy,'k')

V= 1.0/16.0

px = 1.0/6.0
py = 1.0/6.0

c1 = 3.0/2.0/V
c2 = 0
c3 = 0
c4 = -(py + 1.0/3.0*px*px/py)
c5 = 2.0/9.0*V

# c1 = 1
# c2 = -10 
# c3 = 35
# c4 = - 50
# c5 = 24

coeff = [c1,c2,c3,c4,c5]

roots = np.roots(coeff)
print(roots)

print(roots[np.isreal(roots)])

x1 = np.sqrt(2.0/9.0*V)

x2 = 0

ii  = 0

while np.abs(x1-x2)>0.0000001:
	x2 = - ( c1 * pow(x1,4.0) + c5 ) / c4
	x3 = x1
	x1 = x2
	x2 = x3
	ii = ii+1

print(ii,x2)

# plt.show()