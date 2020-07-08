import h5py
import numpy as np
import sys
import os
os.system('cp output.h5 output_old.h5')
f=h5py.File('output_old.h5','r')
list_vars=f.keys()
f3=h5py.File('output.h5','r+')
times = f['u'].keys()
time_series1 = []
for item in times:
    time_series1.append(float(item))
time_series1 = np.sort(np.array(time_series1))
time_series2 = ['%.6f' % x for x in time_series1]
if (len(sys.argv) ==3):
    input_name = sys.argv[1]+'.h5'
    dataset = sys.argv[2]
    if (dataset[0] == '0'):
        dataset = dataset[1:]
    d_frame = np.where(time_series1 == float(dataset))
    if (d_frame[0][0] < len(time_series2)-1):
        for item in time_series2[d_frame[0][0]+1:]:
            if(item[0]=='0'):
                item = item[1:]
            for key in list_vars:
                del f3[key][item]
else:
    input_name = 'input.h5'
    dataset = time_series2[-1]
    if (dataset[0] == '0'):
        dataset = dataset[1:]
f2 = h5py.File(input_name,'w')
for key in f:
    f2.create_group(key)
    grp = f2[key]
    grp.create_dataset('init', data=np.array(f[key][dataset]))
f.close()
f2.close()
f3.close()
f4=open('starttime.txt','w')
f4.write(dataset)
f4.close()
print('New input file ready')
