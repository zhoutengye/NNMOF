from pyevtk.vtk import VtkFile, VtkRectilinearGrid
import numpy as np
from pyevtk.vtk import VtkGroup
import h5py
import os
import sys
from mayavi import mlab
import subprocess

preplot_path = '/usr/local/tecplot/360ex_2018r1/bin/'

f = h5py.File('output.h5','r')
vars = list(f.keys())
scalar_vars = []
scalar_vars_string = []
for item in vars:
    scalar_vars.append(item)
    scalar_vars_string.append(str(item))
times = list(f['phi'].keys())

time_series1 = []
time_series2 = []
for item in times:
    time_series1.append(float(item))
time_series1 = np.sort(np.array(time_series1))
time_series2 = ["%.6f" % x for x in time_series1]
times = time_series2

nx, ny, nz = f['u'][times[0]].shape
nz = f.attrs['nz'][0]
ny = f.attrs['ny'][0]
nx = f.attrs['nx'][0]
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

num_xyz = nx * ny * nz

if (len(sys.argv) ==2):
    new_dir = sys.argv[1]
else:
    new_dir = 'tecplot'
os.system('mkdir -p ' + new_dir)

Var_Line = "VARIABLES = " + "\"X\" " + ",\"Y\" " + ",\"Z\" "
for key in scalar_vars_string:
    Var_Line += ",\"" + key + "\" "

os.system('mkdir -p ' + new_dir)

f2 = open(new_dir + '/tecplot.dat','w')
f2.close()
# Write file
nn = 1
for step in times:
    Zone_Line = "ZONE T=VEL"  + ",I = " + str(nx) + " ,J = " + str(ny) + " ,K = " + str(nz) + ',f=point'
    f2 = open(new_dir + '/tecplot.dat','a')
    f2.write(Var_Line+'\n')
    f2.write(Zone_Line+'\n')
    f2.close()
    len_col = len(scalar_vars) + 3
    data_frame = np.zeros([num_xyz, len_col])
    data_frame[:,0] = Z.flatten()
    data_frame[:,1] = Y.flatten()
    data_frame[:,2] = X.flatten()
    i = 3
    for key in scalar_vars:
        data = np.array(f[key][step])
        data_frame[:,i] = data.flatten()
        i += 1
    with open(new_dir + '/tecplot.dat', 'ab') as f_handle:
        np.savetxt(f_handle, data_frame)

os.system(preplot_path + 'preplot ' + new_dir + '/tecplot.dat ' + new_dir + '/tecplot.plt && rm ' + new_dir + '/tecplot.dat')
