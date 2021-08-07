#include 'param.h'
Module InitSH
  use ModGlobal, only : sp
  use ModTools
Contains
  Real(8) Function ShapeZalesak(x,y,z)
    Implicit None
    Real(sp) :: x, y, z
    Real(sp) :: lv1, lv2, lv3
    lv1 = 0.2_sp - sqrt(x**2+y**2)
    lv2 = min(0.04_sp-y, 0.04_sp+y)
    lv3 = 0.1_sp-x
    ShapeZalesak = min( -min(lv3, lv2), lv1)
  End Function ShapeZalesak

End Module InitSH

