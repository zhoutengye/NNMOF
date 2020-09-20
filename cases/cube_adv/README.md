# workflow

- in **init**, make and run with single processor to generate the initial condition. A file named **output.h5** (which is actually input file) will be generated.

- Copy the file **output.h5** to the **inputs** and rename to coressponding namelist file

- copy the **dt_coef** file to the input 

- copy the Benchmark directory and rename, modify necessary fortran, python, json, then 

- run **gen_case.py** to generate cases, **run_cases** to bath run case, **clean_cases** to clean case

- Post processing using Jupyter-notebook

