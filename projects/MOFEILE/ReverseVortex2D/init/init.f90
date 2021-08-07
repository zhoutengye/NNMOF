#include 'param.h'
Module InitSH
  use ModGlobal, only : sp
  use ModTools
Contains
  Real(8) Function ShapeSingleVortex(x,y,z)
    Implicit None
    Real(sp) :: x, y, z
    Real(sp) :: lv1, lv2, lv3
    ShapeSingleVortex = 0.15_sp - sqrt(x**2+y**2)
  End Function ShapeSingleVortex

End Module InitSH

