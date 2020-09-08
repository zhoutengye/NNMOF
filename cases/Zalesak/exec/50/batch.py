import numpy as np
import json


with open('setup.json') as f:
    data = json.load(f)

input_path = data["inputfile_path"]
ngrid = data["grid"]
case_list = data["exec_list"]
input_namelist_file = input_path+"input"+str(ngrid)+".namelist"
input_data_file = input_path+"input"+str(ngrid)+".h5"

f = open('batch.sh','w')
for case in case_list:
    f.write("cd "+ case + "\n")
    f.write("make clean\n")
    f.write("make \n")
    f.write("cp " + input_namelist_file + " . \n")
    f.write("cp " + input_data_file + " . \n")
    f.write("./Zalesak* input" + str(ngrid) + "\n")
    f.write("cd ..\n")
