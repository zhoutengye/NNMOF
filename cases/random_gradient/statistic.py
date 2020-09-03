import numpy as np
import pandas as pd

# Directories
ds_list = [
    "exponential",
    "uniform",
    "extreme"
]

# Cases
case_list = [
    "Lemoine, BFGS, analytic_info.dat",
    "Lemoine, BFGS, analytic_old_info.dat",
    "Lemoine, BFGS, numerical_info.dat",
    "Lemoine, BFGS, numerical_old_info.dat",
    "Lemoine, Gauss-Newton_info.dat",
    "Lemoine, Gauss-Newton_old_info.dat"
]

# Dictionary that stores everything
stat = {}
for ds in ds_list:
    da = {}
    for case in case_list:
        f=open(ds+"/"+case,'r')
        cas = {}
        for line in f:
            if "time" in line:
                time = float(line.split('=')[1])
                cas.update({"time":time})
            if "num_iter1" in line:
                iter = float(line.split('=')[1])
                cas.update({"iter":iter})
            if "num_iter2" in line:
                ls = float(line.split('=')[1])
                cas.update({"ls":ls})
            if "centroid_error" in line:
                error = float(line.split('=')[1])
                cas.update({"error":error})
        da.update({case:cas})
    stat.update({ds:da})

# Export error
errors = {}
for di in stat:
    da = {}
    for cas in stat[di]:
        da.update({cas:stat[di][cas]["error"]})
    errors.update({di:da})

pd_error = pd.DataFrame.from_dict(errors)
pd_error.to_csv("error.csv")

# Export time
times = {}
for cas in case_list:
    iter1 = 0.0
    iter2 = 0.0
    ti = 0.0
    ca = {}
    for di in ds_list:
        iter1 += stat[di][cas]["iter"]
        iter2 += stat[di][cas]["ls"]
        ti += stat[di][cas]["time"]
    ca.update({"iter":iter1/3.0})
    ca.update({"ls":iter2/3.0})
    ca.update({"time":ti/3.0})
    times.update({cas:ca})
for cas in case_list:
    ti2 = stat[di][cas]["time"] / stat[di]["Lemoine, Gauss-Newton_info.dat"]["time"]
    times[cas].update({"time2":ti2})
pd_error = pd.DataFrame.from_dict(times,orient='index')
pd_error.to_csv("time.csv") 