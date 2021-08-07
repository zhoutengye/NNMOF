import numpy as np
from pyevtk.hl import gridToVTK, pointsToVTK

nn=100

pts = np.loadtxt('pts1.dat')
dp2 = np.ascontiguousarray(pts[:,3])
dt2 = np.ascontiguousarray(pts[:,4])
x2 = np.ascontiguousarray(pts[:,0])
y2 = np.ascontiguousarray(pts[:,1])
z2 = np.ascontiguousarray(pts[:,2])
f2 = np.ascontiguousarray(pts[:,5])

X = np.ascontiguousarray(np.reshape(x2,[nn,nn,nn]))
Y = np.ascontiguousarray(np.reshape(y2,[nn,nn,nn]))
Z = np.ascontiguousarray(np.reshape(z2,[nn,nn,nn]))
P = np.ascontiguousarray(np.reshape(dp2,[nn,nn,nn]))
T = np.ascontiguousarray(np.reshape(dt2,[nn,nn,nn]))
F = np.ascontiguousarray(np.reshape(f2,[nn,nn,nn]))

X2, Y2, Z2=np.mgrid[-0.5:0.5:100j, -0.5:0.5:100j, -0.5:0.5:100j]

# gridToVTK("pt2", X2, Y2, Z2, pointData={"f": F})
# gridToVTK("pt2", X2, Y2, Z2, pointData={"dp": P,"dt":T})

gridToVTK("pt2", X2, Y2, Z2, pointData={"f": F})
gridToVTK("pt2", X2, Y2, Z2, pointData={"dp": P,"dt":T})
