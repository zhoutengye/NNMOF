import numpy as np
import matplotlib.pyplot as plt
import os

nodes = np.zeros(29)
e1 = np.zeros(29)
e2 = np.zeros(29)
for i in range(29):
  de = int(i+1)

  f1 = open('batch/dt_'+str(de)+'.dat')
  ll = 0
  for line in f1:
    if ll == 1:
      nodes[i] = int(line)
      f1.close()
      break
    ll = ll+1

  f2 = open('batch/dt_'+str(de)+'score.dat')
  for line in f2:
    st1 = line.split('.')[1]
    st2 = line.split('.')[2]
    e1[i] = str('0.'+st1)
    e2[i] = str('0.'+st2)

# print(e1)
# plt.plot(e1)
# plt.plot(e2)
plt.plot(nodes)
plt.show()