# Workflow

1. copy the **namelist** and **h5** file from ```Zalesk/exec/init``` directory

2. Go to ```MOFBFGS``` and ```MOFLemoineGN``` directory, make and then rename the executable with it marked as **improved**.

3. Go to ```src/vof_func.f90```, uncomment the last line in 
   subroutine ```Initial_Guess``` 

4. Go to the two cases, make again and rename th executable../s
   as **original**.

5. Revert the file in step 3.

- If other cases is to be done, it would be similar with Zalesak

- There are not a lot of data to deal with, so the cpu time are simpled obtaiend from the print-screen.