import os

f=open('batch_fortran.sh','w')

for i in range(30):
  de = str(i+1)
  f.write('cp batch/dt_' + str(de) + '.dat' + ' dt_coef.dat\n')
  f.write('./gen_data.exec\n')
  f.write('cp dt_time.dat batch/dt_time' + str(de) + '.dat\n')
f.close()