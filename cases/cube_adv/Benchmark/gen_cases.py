import numpy as np
import json
import sys
import os

with open('setup.json') as f:
    data = json.load(f)

input_path = data["inputfile_path"]
case_list = data["exec_list"]
bench_dir = data["BenchDir"]
script_dict = data["CaseParams"]
c_case_list = script_dict.keys()

input_namelist_file = input_path+"input50.namelist"
input_data_file = input_path+"input50.h5"
input_ml_file = input_path+"dt_coef.dat"

for case in case_list:
    if case not in c_case_list:
        print("*** Fatal Error ***")
        print("Case " + case + " not in the parameter list, please check \n")
        sys.exit("Error message")
    os.system("cp -r " + bench_dir + " " + case )
    os.system("cp " + input_namelist_file + ' '+ case + '/')
    os.system("cp " + input_data_file + " " + case + '/')
    os.system("cp " + input_ml_file + " " + case + '/')
    os.system("cp " + input_namelist_file + ' '+ case + '/')
    main_file = case + "/main.f90"
    new_f = []
    f2 = open(case + '/main.f90','r')
    replace_flag = 0
    for line in f2:
        if "!>1" in line:
            replace_flag = 1
            for n_line in script_dict[case]["Reconstruction"]:
                new_f.append(n_line+"\n")
        if "!<1" in line:
            replace_flag = 0
        if "!>2" in line:
            replace_flag = 1
            if script_dict[case]["Type"] == "MOF":
                new_f.append("Call MOFCIAM(Phi, u, v, w, nl, dl, dt,rank, cx, cy, cz)\n")
            elif script_dict[case]["Type"] == "VOF":
                new_f.append("Call VOFCIAM(Phi, u, v, w, nl, dl, dt,rank)\n")
            elif script_dict[case]["Type"] == "THINC":
                new_f.append("Call VOFTHINC(Phi, u, v, w, nl, dl, dt,rank)\n")
            else:
                print("*** Fatal Error ***")
                print("*** VOF type not supported, should be MOF, VOF or THINC ***")
        if "!<2" in line:
            replace_flag = 0
        if replace_flag == 0:
                new_f.append(line)
    f2.close()
    f3 = open(case + '/main.f90','w')
    for line in new_f:
        f3.write(line)
        # print(line)
    f3.close()

f = open('run_cases.sh','w')
for case in case_list:
    f.write("echo 'start case "+ case + "'\n")
    f.write("cd "+ case + "\n")
    f.write("rm Linear* \n")
    f.write("make clean\n")
    f.write("make \n")
    f.write("./Linear* input\n")
    f.write("cd ..\n")
    f.write("echo 'end case "+ case + "'\n")
f.close()

f = open('clean_cases.sh','w')
for case in case_list:
    f.write("rm -r " + case + "\n")
f.close()
