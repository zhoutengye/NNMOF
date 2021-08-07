# Random gradient
## Basic ideas
Randomly generate data corresponds with the linear reconstruction.
Three data-sets in three directories
- extreme 
- uniform
- normal

Each data set is evaluated with 6 mathods

| Method       | gradient  | initial guess |
| ------------ | --------- | ------------- |
| Gauss-Newton | analytic  | original      |
| BFGS         | analytic  | original      |
| BFGS         | numerical | original      |
| Gauss-Newton | analytic  | improved      |
| BFGS         | analytic  | improved      |
| BFGS         | numerical | improved      |

In each eirectory, the executable will be compiled with **Makefile**,

## Work flow
In each directory,
- Compile with **Makefile** to get three executables:
    
    **TEST1.exec**: generate test data
    **TEST2.exec**: executa MOF with improved initial guess
    **TEST3.exec**: executa MOF with original initial guess

- Generate a random list of volume fraction with **gen_f.py**.
- Run **TEST_1**
- Run **TEST_2** and **TEST_3** and output **.dat** data.
- Post-processing

## Implementation by command

- step 1: generate batch-run script

**python makeall.py**

- step 2: compile all executables

**bash makeall.sh**

- step 3: generate all test data set and evaluate MOF

**bash runeall.sh**

- step 4: post processing, gather all output data

**python statistics.py**

