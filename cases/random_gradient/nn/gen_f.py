import numpy as np
import matplotlib.pyplot as plt
num_sampling = 1000000
num_sampling2 = int(num_sampling / 2)
num_sampling4 = int(num_sampling / 4)
f = np.random.normal(loc=0.0, scale=0.1, size=num_sampling2) 
f1 = (f + 0.5) / 2
f = np.random.normal(loc=0.0, scale=0.1, size=num_sampling2) 
f2 = (f + 0.5) / 2+ 0.5
f=np.concatenate([f1,f2])
x = np.arange(0,1,0.1)
np.savetxt('fdata.dat',f)
