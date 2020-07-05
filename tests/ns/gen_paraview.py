from pyevtk.vtk import VtkFile, VtkRectilinearGrid
import numpy as np
from pyevtk.vtk import VtkGroup
import h5py

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
x = np.arange(nx)
y = np.arange(ny)
z = np.arange(nz)

#Define some parameters
start, end = (1,1,1), (nx, ny, nz) #Modify 0->1


#Add to write file. 
for step in times:
    filename=str(step)+'.vtr'
    w = VtkFile(str(step), VtkRectilinearGrid) #evtk_test0

    # w = VtkFile("./filename", VtkRectilinearGrid) #evtk_test0
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
    w.appendData(x).appendData(y).appendData(z)
    w.save()

g = VtkGroup("./group")
for step in times:
    g.addFile(filepath = str(step)+'.vtr', sim_time = float(step))
g.save()


