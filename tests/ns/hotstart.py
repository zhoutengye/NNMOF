import h5py
import numpy as np

dataset = 'modify here'
input_name = 'modify here'
input_name += '.h5'

f=h5py.File('output.h5','r')
f2 = h5py.File(input_name,'w')

for key in f:
    f2.create_group(key)
    grp = f2[key]
    grp.create_dataset("init", data=np.array(f[key][dataset]))

f.close()
f2.close()

