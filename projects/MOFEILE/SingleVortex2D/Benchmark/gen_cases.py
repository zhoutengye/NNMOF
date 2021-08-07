import numpy as np
import json
import sys
import os

with open('setup.json') as f:
    data = json.load(f)

reso = '16'
cfl = '0.3'

input_path = data["inputfile_path"]
case_list = data["exec_list"]
bench_dir = data["BenchDir"]
exec_path = data["exec_path"]
script_dict = data["CaseParams"]
c_case_list = script_dict.keys()


input_namelist_file = input_path+"/cfl_" + cfl + "/input" + reso + ".namelist"
input_data_file = input_path+"input" + reso + ".h5"
input_ml_file = input_path+"dt_coef.dat"

exec_line = "./*.exec input" + reso + " \n"

for case in case_list:
    if case not in c_case_list:
        print("*** Fatal Error ***")
        print("Case " + case + " not in the parameter list, please check \n")
        sys.exit("Error message")
    os.system("mkdir -p " + case )
    print("cp " + exec_path + case + ".exec " + case + '/')
    os.system("cp " + exec_path + case + ".exec " + case + '/')
    os.system("cp " + input_namelist_file + ' '+ case + '/')
    os.system("cp " + input_data_file + " " + case + '/')
    os.system("cp " + input_ml_file + " " + case + '/')
    main_file = case + "/main.f90"
    new_f = []

f = open('run_cases.sh','w')
for case in case_list:
    f.write("echo 'start case "+ case + "'\n")
    f.write("cd "+ case + "\n")
    f.write(exec_line)
    f.write("cd ..\n")
    f.write("echo 'end case "+ case + "'\n")
f.close()

f = open('clean_cases.sh','w')
for case in case_list:
    f.write("rm -r " + case + "\n")
f.close()
