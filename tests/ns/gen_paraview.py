from pyevtk.vtk import VtkFile, VtkStructuredGrid
import numpy as np
from pyevtk.vtk import VtkGroup
import h5py
import os
import sys
f = h5py.File('output.h5','r')
vars = list(f.keys())
scalar_vars = []
scalar_vars_string = []
for item in vars:
    if (item != 'u' and item != 'v' and item != 'w'):
        scalar_vars.append(item)
        scalar_vars_string.append(str(item))
times = list(f['phi'].keys())
nz = f.attrs['nx'][0]
ny = f.attrs['ny'][0]
nx = f.attrs['nz'][0]
dz = f.attrs['dx'][0]
dy = f.attrs['dy'][0]
dx = f.attrs['dz'][0]
x = np.arange(nx) * dx + dx * 0.5
y = np.arange(ny) * dy + dy * 0.5
z = np.arange(nz) * dz + dz * 0.5
X, Y, Z = np.mgrid[dz*0.5:nz*dx+dz*0.5:dz,
                   dy*0.5:ny*dy+dy*0.5:dy,
                   dx*0.5:nx*dz+dx*0.5:dx
                   ]
start, end = (1,1,1), (nx, ny, nz) #Modify 0->1
if (len(sys.argv) ==2):
    new_dir = sys.argv[1]
else:
    new_dir = 'paraview'
new_data = new_dir + '/data'
os.system('mkdir -p ' + new_dir)
os.system('mkdir -p ' + new_data)
for step in times:
    filename=new_data + '/'+ str(step)
    w = VtkFile(filename, VtkStructuredGrid) #evtk_test0
    w.openGrid(start = start, end = end)
    w.openPiece( start = start, end = end)
    w.openData('Point', scalars = scalar_vars_string, vectors = 'Velocity')
    for key in scalar_vars:
        w.addData(str(key),np.array(f[key][step]))
    vx = np.array(f['u'][step])
    vy = np.array(f['v'][step])
    vz = np.array(f['w'][step])
    w.addData('Velocity', (vx,vy,vz))
    w.closeData('Point')
    w.openElement('Points')
    w.addData('points', (X, Y, Z))
    w.closeElement('Points')
    w.closePiece()
    w.closeGrid()
    for key in scalar_vars:
        w.appendData(data = np.array(f[key][step]))
    w.appendData(data = (vx,vy,vz))
    w.appendData((X, Y, Z))
    w.save()
    print('file: '+filename+' added')
g = VtkGroup(new_dir+'/group')
for step in times:
    g.addFile(filepath = new_data + '/' + str(step)+'.vts', sim_time = float(step))
g.save()
print('group file: group.pvd added')
