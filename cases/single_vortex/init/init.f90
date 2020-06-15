Module InitSH
  use ModTools
Contains
  Real(8) Function ShapeRiderKothe(x,y,z)
    Implicit None
    Real(sp) :: x, y, z
    ShapeRiderKothe = 0.15_sp - sqrt(x**2+y**2+z**2)
  End Function ShapeRiderKothe


End Module InitSH

