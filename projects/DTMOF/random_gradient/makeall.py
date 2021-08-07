import numpy as np
import matplotlib.pyplot as plt

dir_list = ['extreme','uniform','exponential']
f_make = open('makeall.sh','w')
for di in dir_list:
    f_make.write("cd "+di+'\n')
    for i in range(3):
        f_make.write("make clean\n")
        f_make.write("sed -i '1d' Makefile\n")
        f_make.write("echo 'TEST_NUM = "+str(i+1)+"' | cat - Makefile > temp && mv temp Makefile\n")
        f_make.write("make\n")
    f_make.write("make clean\n")
    f_make.write("cd ..\n")
f_make.close()

f_run = open('runeall.sh','w')
for di in dir_list:
    f_run.write("cd "+di+'\n')
    f_run.write("python gen_f.py\n")
    f_run.write("./TEST1.exec\n")
    f_run.write("cd ..\n")
f_run.write("\n")
for di in dir_list:
    f_run.write("cd "+di+'\n')
    f_run.write("./TEST2.exec\n")
    f_run.write("./TEST3.exec\n")
    f_run.write("cd ..\n")
f_run.close()
