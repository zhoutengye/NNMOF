import numpy as np

exact_f = np.loadtxt('exact_f.dat')
exact_angle = np.loadtxt('exact_angle.dat')
exact_centroid = np.loadtxt('exact_centroid.dat')
delta_angle = np.loadtxt('delta_angle.dat')
initial_angle = np.loadtxt('initial_angle.dat')

np.save('exact_f.npy',exact_f)
np.save('exact_angle.npy',exact_angle)
np.save('exact_centroid.npy',exact_centroid)
np.save('delta_angle.npy',delta_angle)
np.save('initial_angle.npy',initial_angle)