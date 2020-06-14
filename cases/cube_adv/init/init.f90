Module InitSH
  use ModTools
Contains
  Real(8) Function ShapeTest(x,y,z)
    Implicit None
    Real(8) :: x, y, z
    ShapeTest = 1.0_sp-x-y-z
    ShapeTest = x+y+z - 1.0_sp
  End Function ShapeTest

  Real(8) Function ShapeCube1(x,y,z)
    Implicit None
    Real(8) :: x, y, z
    ShapeCube1 = min(x, min(y, z))
  End Function ShapeCube1

  Real(8) Function ShapeCube2(x,y,z)
    Implicit None
    Real(8) :: x, y, z
    ShapeCube2 = min(0.5d0-x, min(0.5d0-y, 0.5d0-z))
  End Function ShapeCube2

  Real(8) Function ShapeSphere(x,y,z)
    Implicit None
    Real(8) :: x, y, z
    Real(8) :: sph1, sph2
    sph1 = 1.0_sp - sqrt(x**2+y**2+z**2)
    sph2 = sqrt(x**2+y**2+z**2) - 0.5_sp
    ShapeSphere = min(sph1, sph2)
  End Function ShapeSphere

  Real(8) Function ShapeTCube1(x,y,z)
    Implicit None
    Real(8) :: x, y, z
    ShapeTCube1 = 1.0_sp - (x+y+z)
  End Function ShapeTCube1

  Real(8) Function ShapeTCube2(x,y,z)
    Implicit None
    Real(8) :: x, y, z
    ShapeTCube2 = 0.5_sp - (x+y+z)
  End Function ShapeTCube2

  Real(8) Function ShapeLA1(x,y,z)
    Implicit None
    Real(8) :: x, y, z
    Real(8) :: la1, la2
    la1 = 1.0 - sqrt(y**2+z**2) - x / 2.0_sp
    la2 =  sqrt(y**2+z**2) + x / 2.0_sp - 0.5_sp
    ShapeLA1 = min(la1,la2)
  End Function ShapeLA1

  Real(8) Function GShapeCube(x,y,z)
    Implicit None
    Real(8) :: x, y, z
    Real(8) :: xm1, ym1, zm1
    Real(8) :: xm2, ym2, zm2
    Real(8) :: sp1, sp2
    xm1 = min(x-5.0_sp, 45.0_sp-x)
    ym1 = min(y-5.0_sp, 45.0_sp-y)
    zm1 = min(z-5.0_sp, 45.0_sp-z)
    xm2 = - min(x-15.0_sp, 35.0_sp-x)
    ym2 = - min(y-15.0_sp, 35.0_sp-y)
    zm2 = - min(z-15.0_sp, 35.0_sp-z)
    sp1 = min(xm1, min(ym1,zm1))
    sp2 = max(xm2, max(ym2,zm2))
    GShapeCube1 = min(sp1,sp2)
  End Function GShapeCube


  Real(8) Function GShapeSphere(x,y,z)
    Implicit None
    Real(8) :: x, y, z
    Real(8) :: sph1, sph2
    sph1 = 1.0_sp - sqrt(x**2+y**2+z**2)
    sph2 = sqrt(x**2+y**2+z**2) - 0.5_sp
    ShapeSphere = min(sph1, sph2)
  End Function GShapeSphere

  ! Real(8) Function GShapeTCube1(x,y,z)
  !   Implicit None
  !   Real(8) :: x, y, z
  !   ShapeTCube1 = 1.0_sp - (x+y+z)
  ! End Function GShapeTCube1

  ! Real(8) Function GShapeTCube2(x,y,z)
  !   Implicit None
  !   Real(8) :: x, y, z
  !   ShapeTCube2 = 0.5_sp - (x+y+z)
  ! End Function GShapeTCube2

  ! Real(8) Function GShapeLA1(x,y,z)
  !   Implicit None
  !   Real(8) :: x, y, z
  !   Real(8) :: la1, la2
  !   la1 = 1.0 - sqrt(y**2+z**2) - x / 2.0_sp
  !   la2 =  sqrt(y**2+z**2) + x / 2.0_sp - 0.5_sp
  !   ShapeLA1 = min(la1,la2)
  ! End Function GShapeLA1

End Module InitSH

