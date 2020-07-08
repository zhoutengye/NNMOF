## Two command line input arguments
# first one is the name for the input file for hotstart
# second one is the name start data set.
# If no command line input argument, it uses default name pick the last dataframe.
import h5py
import numpy as np
import sys
# make a copy of the old output file
os.system('cp output.h5 output_old.h5')
# open files
f=h5py.File('output_old.h5','r')
f3=h5py.File('output.h5','r+')
# get informations
list_vars = f.keys()
times = f['u'].keys()
time_series1 = []
# sort time series
for item in times:
    time_series1.append(float(item))
time_series1 = np.sort(np.array(time_series1))
time_series2 = ['%.6f' % x for x in time_series1]

# determine the hotstart frame.
# If not the last frame, delete the superfulours data
if (len(sys.argv) ==3):
    input_name = sys.argv[1]+'.h5'
    dataset = sys.argv[2]
    if (dataset[0] == '0'):
	    dataset = dataset[1:]
    d_frame = np.where[float(dataset)]
    if (d_frame < len(time_series2)-1):
        for item in time_series2(d_frame+1:):
            if (dataset[0] == '0'):
            dataset = dataset[1:
            for key in list_vars:
                del f3[list_vars][item]
else:
    input_name = 'input.h5'
    dataset = time_series2[-1]

# write data to the new initial file
f2 = h5py.File(input_name,'w')
for key in f:
    f2.create_group(key)
    grp = f2[key]
    grp.create_dataset('init', data=np.array(f[key][dataset]))

# close h5files
f.close()
f2.close()
f3.close()


