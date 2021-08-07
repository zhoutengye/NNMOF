import numpy as np
import matplotlib.pyplot as plt

num_sampling = 10000000
# num_sampling2 = int(num_sampling / 2)
# num_sampling4 = int(num_sampling / 4)
f = np.random.uniform(low=0.0,high=1.0,size=num_sampling)
m=plt.hist(f,bins=100)
np.savetxt('fdata.dat',f)
