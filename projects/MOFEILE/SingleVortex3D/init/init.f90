#include 'param.h'
Module InitSH
  use ModGlobal, only : sp
  use ModTools
Contains
  Real(8) Function ShapeRV(x,y,z)
    Implicit None
    Real(sp) :: x, y, z
    Real(sp) :: lv1, lv2, lv3
    ShapeRV = 0.15_sp - sqrt(x**2+y**2+z**2)
  End Function ShapeRV

End Module InitSH

