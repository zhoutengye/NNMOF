import h5py
import numpy as np
import sys


if (len(sys.argv) ==2):
    input_name = sys.argv[1]+'.h5'
else:
    input_name = 'input.h5'


f=h5py.File('output.h5','r')


times = f['u'].keys()
time_series1 = []
for item in times:
    time_series1.append(float(item))
time_series1 = np.sort(np.array(time_series1))
time_series2 = ['%.6f' % x for x in time_series1]
dataset = time_series2[-1]
if (dataset[0] == '0'):
	dataset = dataset[1:]

## By defatult, it will sxtract the last frame. 
## To specify the frame, uncomment and type the name of the data set below
# dataset = ''


f2 = h5py.File(input_name,'w')

for key in f:
    f2.create_group(key)
    grp = f2[key]
    grp.create_dataset('init', data=np.array(f[key][dataset]))

f.close()
f2.close()

