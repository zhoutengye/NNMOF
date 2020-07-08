import h5py
import numpy as np
import sys
import os

hottime = .true.
input_name = 'input.h5'
f=h5py.File('output.h5','r')
var_list = f.keys()
times = f['u'].keys()
time_series1 = []
for item in times:
    time_series1.append(float(item))
time_series1 = np.sort(np.array(time_series1))
time_series2 = ['%.6f' % x for x in time_series1]
if hottime:
    dataset = time_series2[-1]
    hot_data = np.array(f[key][dataset])
    f.close()
else:
    dataset = ''
    hot_data = np.array(f[key][dataset])
    os.system('cp output.h5 output_old.h5')
    f.close()
    f=h5py.File('output.h5','r+')
    d_frame = np.where(time_series1 == 0.11)[0][0]
    for item in time_series2[fd:]:
        if item[0] == '0':
            item = item[1:]
        for key in var_list:
            del f[key][item]
    f.close
if (dataset[0] == '0'):
    dataset = dataset[1:]
    hot_data = np.array(f[key][dataset])
f2 = h5py.File(input_name,'w')
for key in f:
    f2.create_group(key)
    grp = f2[key]
    grp.create_dataset('init', data=hot_data)
f2.close()
