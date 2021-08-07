Module InitSH
  Use MOdGlobal, only: sp
  use ModTools
Contains
  Real(8) Function ShapeTest(x,y,z)
    Implicit None
    Real(8) :: x, y, z
    ShapeTest = 4.0_sp-x-y-z
    ShapeTest = x+y+z - 4.0_sp
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
    xm1 = min(x-0.05_sp, 0.45_sp-x)
    ym1 = min(y-0.05_sp, 0.45_sp-y)
    zm1 = min(z-0.05_sp, 0.45_sp-z)
    xm2 = - min(x-0.15_sp, 0.35_sp-x)
    ym2 = - min(y-0.15_sp, 0.35_sp-y)
    zm2 = - min(z-0.15_sp, 0.35_sp-z)
    sp1 = min(xm1, min(ym1,zm1))
    sp2 = max(xm2, max(ym2,zm2))
    GShapeCube = min(sp1,sp2)
  End Function GShapeCube


  Real(8) Function GShapeSphere(x,y,z)
    Implicit None
    Real(8) :: x, y, z
    Real(8) :: sph1, sph2
    x = x - 0.25_sp
    y = y - 0.75_sp
    z = z - 0.25_sp
    sph1 = 0.2_sp - sqrt(x**2+y**2+z**2)
    sph2 = sqrt(x**2+y**2+z**2) - 0.1_sp
    GShapeSphere = min(sph1, sph2)
  End Function GShapeSphere

  Real(8) Function GShapeTCube(x,y,z)
    Implicit None
    Real(8) :: x, y, z
    Real(8) :: d1, d2, d3, d4, d5, d6, d7, d8
    Real(8) :: l1, l2
    x = x - 0.75_sp
    y = y - 0.25_sp
    z = z - 0.25_sp
    d1 = 0.20_sp + ( x+y+z)
    d2 = 0.20_sp + (-x+y+z)
    d3 = 0.20_sp + ( x-y+z)
    d4 = 0.20_sp + (-x-y+z)
    d5 = 0.20_sp + ( x+y-z)
    d6 = 0.20_sp + (-x+y-z)
    d7 = 0.20_sp + ( x-y-z)
    d8 = 0.20_sp + (-x-y-z)
    d1 = min(d1,d2)
    d3 = min(d3,d4)
    d5 = min(d5,d6)
    d7 = min(d7,d8)
    l1 = min(min(d1,d3),min(d5,d7))
    d1 = - 0.1_sp + ( x+y+z)
    d2 = - 0.1_sp + (-x+y+z)
    d3 = - 0.1_sp + ( x-y+z)
    d4 = - 0.1_sp + (-x-y+z)
    d5 = - 0.1_sp + ( x+y-z)
    d6 = - 0.1_sp + (-x+y-z)
    d7 = - 0.1_sp + ( x-y-z)
    d8 = - 0.1_sp + (-x-y-z)
    d1 = max(d1,d2)
    d3 = max(d3,d4)
    d5 = max(d5,d6)
    d7 = max(d7,d8)
    l2 = max(max(d1,d3),max(d5,d7))
    GShapeTCube = min(l1,l2)
  End Function GShapeTCube

  Real(8) Function GShapeLA(x,y,z)
    Implicit None
    Real(8) :: x, y, z
    Real(8) :: d1, d2, la1, la2
    x = x-0.75_sp
    y = y-0.75_sp
    z = z-0.25_sp
    d1 = 0.1_sp - sqrt(z**2+y**2) - x / 2.0_sp
    d2 = x+0.2_sp
    la1 = min(d1,d2)
    d1 = sqrt(z**2+y**2) + x / 2.0_sp - 0.02_sp
    d2 = -x-0.08_sp
    la2 = max(d1,d2)
    la1 = min(la1,la2)
    d1 = sqrt(z**2+y**2) + x / 2.0_sp - 0.02_sp
    d2 = -x-0.12_sp
    d2 =  x+0.15_sp
    la2 = max(d1,d2)
    la2 = min(la1,la2)
    GShapeLA = min(la1,la2)
  End Function GShapeLA

End Module InitSH

