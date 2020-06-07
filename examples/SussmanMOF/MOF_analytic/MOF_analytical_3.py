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

sx = []
sy = []

V= 1/4

px = 1/9
py = 1/4

# x1 = np.sqrt(1/28)
# y1 = V*2/9/x1

x1 = np.sqrt(V*2/9)
y1 = V*2/9/x1

x2 = 0
y2 = 0

ii = 0

print(x1,y1)

while np.abs(x1-x2)>0.00000001 or np.abs(y1-y2)>0.00000001:

	delta11 = - ( x1*y1 - V*2/9 ) * y1  / (x1*x1+y1*y1)
	delta12 = - ( x1*y1 - V*2/9 ) * x1  / (x1*x1+y1*y1)
	x11 = x1 + delta11
	y11 = y1 + delta12

	delta21 = px - x1 - ( (px-x1)* x1 + (py-y1)*y1 ) / (x1*x1+y1*y1) * y1
	delta22 = py - y1 - ( (px-x1)* x1 + (py-y1)*y1 ) / (x1*x1+y1*y1) * x1
	x2 = x11 + delta21
	y2 = y11 + delta22

	print(ii,delta11,delta12)
	
	x3 = x1
	x1 = x2
	x2 = x3
	y3 = y1
	y1 = y2
	y3 = y3
	ii = ii+1

	

print(ii,x1,y1)

# plt.show()