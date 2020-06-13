#include "param.h"
# if defined TEST1
# define test_case test1
# elif defined TEST2
# define test_case test2
# elif defined TEST3
# define test_case test3
# elif defined TEST4
# define test_case test4
# endif
program test
  use ModGlobal
  Implicit None

  call test_case

end program test

!===============
! test1: MISC:
!    (1) Norm2Angle and Angle2Norm
!    (2) Flood_BackwardC
!    (3) FindCentroid
!    (4) MOFZY
!===============
Subroutine test1
  ! Use ModGlobal
  Use ModTools
  Use InitSH
  Implicit None

  ! Call Init(inputfield=.false.)

  Type(Octree) :: gr

  gr%level = 1
  gr%maxlevel = 5
  gr%xc(1) = 0.5_sp
  gr%xc(2) = 0.5_sp
  gr%xc(3) = 0.5_sp
  gr%dx(1) = 1.0_sp
  gr%dx(2) = 1.0_sp
  gr%dx(3) = 1.0_sp
  ShapeLevelSet => ShapeTest

  Call VolumeCentroidQuadTree(gr)

  print *, gr%vof
  print *, gr%centroid
  ! print *, 'tt'

  gr%level = 1
  gr%maxlevel = 5
  gr%xc(1) = 0.5_sp
  gr%xc(2) = 0.5_sp
  gr%xc(3) = 0.5_sp
  gr%dx(1) = 1.0_sp
  gr%dx(2) = 1.0_sp
  gr%dx(3) = 1.0_sp
  ShapeLevelSet => ShapeTest
  Call VolumeQuadTree(gr)
  print *, gr%vof

  ! Call MPI_FINALIZE(ierr)

end Subroutine test1

Subroutine test2
  Use ModTools
  Use InitSH
  Implicit None
  Real(sp) :: sph(20,20,20)
  Real(sp) :: xc(20)
  Type(Octree) :: gr
  Integer :: i,j,k
  Integer :: flag
  Real(sp) :: dx

  Call Init(inputfield=.false.)

  Do k = 1, 20
    xc(k) = dble(k) / dble(20) - 0.025_sp
  End Do

  dx = 0.05_sp

  ShapeLevelSet => ShapeSphere

  Do k = 1, 20
    Do j = 1, 20
      Do i = 1, 20
        flag = Inout(xc(i), xc(j), xc(k), dx, dx, dx)
        If (flag .eq. 1) Then
          sph(i,j,k) = 1.0
        Elseif (flag .eq. 3) Then
          sph(i,j,k) = 0.0
        Elseif (flag .eq. 2) Then
          gr%level = 1
          gr%maxlevel = 10
          gr%xc(1) = xc(i)
          gr%xc(2) = xc(j)
          gr%xc(3) = xc(k)
          gr%dx(1) = dx
          gr%dx(2) = dx
          gr%dx(3) = dx
          Call VolumeQuadTree(gr)
          sph(i,j,k) = gr%vof
        EndIf
        phi(i,j,k) = sph(i,j,k)
      End Do
    End Do
  End Do

  Call Visual3DContour(f1=phi)

End Subroutine test2


