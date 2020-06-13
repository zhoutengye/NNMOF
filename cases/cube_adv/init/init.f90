Module InitSH
  use ModTools
Contains
  Real(8) Function ShapeTest(x,y,z)
    Implicit None
    Real(8) :: x, y, z
    ShapeTest = 1.0_sp-x-y-z
  End Function ShapeTest

  Real(8) Function ShapeSphere(x,y,z)
    Implicit None
    Real(8) :: x, y, z
    ShapeSphere = 1.0_sp - sqrt(x**2+y**2+z**2)
  End Function ShapeSphere

End Module InitSH

