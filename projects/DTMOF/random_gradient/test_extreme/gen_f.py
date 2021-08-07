import numpy as np
import matplotlib.pyplot as plt
num_sampling = 1000000
num_sampling2 = int(num_sampling / 2)
num_sampling4 = int(num_sampling / 4)
f = np.random.normal(loc=0.0, scale=0.1, size=num_sampling2) 
f3 = 0.5-f[f>0]
f4 = -0.5-f[f<=0]
f1 = (np.concatenate([f3,f4]) + 0.5) / 2
f = np.random.normal(loc=0.0, scale=0.1, size=num_sampling2) 
f3 = 0.5-f[f>0]
f4 = -0.5-f[f<=0]
f2 = (np.concatenate([f3,f4]) + 0.5) / 2+ 0.5
f=np.concatenate([f1,f2])
x = np.arange(0,1,0.1)
m=plt.hist(f,bins=100)
np.savetxt('fdata.dat',f)