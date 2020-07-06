from pyevtk.vtk import VtkFile, VtkRectilinearGrid
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

nx, ny, nz = f['u'][times[0]].shape

#Read data
x = np.arange(nx) / 20.0
y = np.arange(ny) / 20.0
z = np.arange(nz) / 20.0

#Define some parameters
start, end = (1,1,1), (nx, ny, nz) #Modify 0->1

if (len(sys.argv) ==2):
    new_dir = sys.argv[1]
else:
    new_dir = 'paraview'
new_data = new_dir + '/data'
os.system('mkdir -p ' + new_dir)
os.system('mkdir -p ' + new_data)

#Add to write file.
for step in times:
    filename=new_data + '/'+ str(step)
    w = VtkFile(filename, VtkRectilinearGrid) #evtk_test0

    w.openGrid(start = start, end = end)
    w.openPiece( start = start, end = end)

    w.openData("Point", scalars = scalar_vars_string, vectors = "Velocity")
    for key in scalar_vars:
        w.addData(str(key),np.array(f[key][step]))
    vx = np.array(f['u'][step])
    vy = np.array(f['v'][step])
    vz = np.array(f['w'][step])
    w.addData("Velocity", (vx,vy,vz))
    w.closeData("Point")

    # Coordinates of cell vertices
    w.openElement("Coordinates")
    w.addData("x_coordinates", x);
    w.addData("y_coordinates", y);
    w.addData("z_coordinates", z);
    w.closeElement("Coordinates");

    w.closePiece()
    w.closeGrid()

    # Need to modify parameters
    for key in scalar_vars:
        w.appendData(data = np.array(f[key][step]))
    w.appendData(data = (vx,vy,vz))
    w.appendData(x)
    w.appendData(y)
    w.appendData(z)
    w.save()
    print("file: "+filename+" added")

g = VtkGroup(new_dir+"/group")
for step in times:
    g.addFile(filepath = new_data + '/' + str(step)+'.vtr', sim_time = float(step))
g.save()

print("group file: group.pvd added")
