**f2py** is used to compile the Fortran code to a dynamic library so that python can call it.

To compile the dynamic library, execute

f2py -c -m mof mof.f90
