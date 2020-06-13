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

  Call VolumeCentroidOctree(gr)

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
  Call VolumeOctree(gr)
  print *, gr%vof

  ! Call MPI_FINALIZE(ierr)

end Subroutine test1

Subroutine test2
  Use ModTools
  Use InitSH
  Implicit None
  Real(sp) :: sph(20,20,20)
  Real(sp) :: cube1(20,20,20)
  Real(sp) :: cube2(20,20,20)
  Real(sp) :: tcube1(20,20,20)
  Real(sp) :: tcube2(20,20,20)
  Real(sp) :: LA1(40,20,20)
  Real(sp), Allocatable :: xc(:)
  Type(Octree) :: gr
  Integer :: i,j,k
  Integer :: icc,jcc,kcc
  Integer :: flag
  Real(sp) :: dx

  Call Init(inputfield=.false.)


  ! ---------------------------------
  ! Hollow cube
  Allocate(xc(20))
  Do k = 1, 20
    xc(k) = dble(k) / dble(20) - 0.025_sp
  End Do
  dx = 0.05_sp
  ShapeLevelSet => ShapeCube1
  Do k = 1, 20
    Do j = 1, 20
      Do i = 1, 20
        flag = Inout(xc(i), xc(j), xc(k), dx/2.0_sp, dx/2.0_sp, dx/2.0_sp)
        If (flag .eq. 1) Then
          cube1(i,j,k) = 1.0
        Elseif (flag .eq. 3) Then
          cube1(i,j,k) = 0.0
        Elseif (flag .eq. 2) Then
          gr%level = 1
          gr%maxlevel = 5
          gr%xc(1) = xc(i)
          gr%xc(2) = xc(j)
          gr%xc(3) = xc(k)
          gr%dx(1) = dx
          gr%dx(2) = dx
          gr%dx(3) = dx
          Call VolumeOctree(gr)
          cube1(i,j,k) = gr%vof
        EndIf
      End Do
    End Do
  End Do
  ShapeLevelSet => ShapeCube2
  Do k = 1, 20
    Do j = 1, 20
      Do i = 1, 20
        flag = Inout(xc(i), xc(j), xc(k), dx/2.0_sp, dx/2.0_sp, dx/2.0_sp)
        If (flag .eq. 1) Then
          cube2(i,j,k) = 1.0
        Elseif (flag .eq. 3) Then
          cube2(i,j,k) = 0.0
        Elseif (flag .eq. 2) Then
          gr%level = 1
          gr%maxlevel = 5
          gr%xc(1) = xc(i)
          gr%xc(2) = xc(j)
          gr%xc(3) = xc(k)
          gr%dx(1) = dx
          gr%dx(2) = dx
          gr%dx(3) = dx
          Call VolumeOctree(gr)
          cube2(i,j,k) = gr%vof
        EndIf
      End Do
    End Do
  End Do

  icc = 25
  jcc = 25
  kcc = 25
  Do k = 1, 20
    Do j = 1, 20
      Do i = 1, 20
        Phi(icc+i  , jcc+j  , kcc+k  ) = cube1(i,j,k) - cube2(i,j,k)
        Phi(icc+1-i, jcc+j  , kcc+k  ) = cube1(i,j,k) - cube2(i,j,k)
        Phi(icc+i  , jcc+1-j, kcc+k  ) = cube1(i,j,k) - cube2(i,j,k)
        Phi(icc+1-i, jcc+1-j, kcc+k  ) = cube1(i,j,k) - cube2(i,j,k)
        Phi(icc+i  , jcc+j  , kcc+1-k) = cube1(i,j,k) - cube2(i,j,k)
        Phi(icc+1-i, jcc+j  , kcc+1-k) = cube1(i,j,k) - cube2(i,j,k)
        Phi(icc+i  , jcc+1-j, kcc+1-k) = cube1(i,j,k) - cube2(i,j,k)
        Phi(icc+1-i, jcc+1-j, kcc+1-k) = cube1(i,j,k) - cube2(i,j,k)
      End Do
    End Do
  End Do
  ! End Hollow cube
  ! ---------------------------------
  
  !----------------------------------
  ! Hollow Sphere
  Do k = 1, 20
    xc(k) = dble(k) / dble(20) - 0.025_sp
  End Do
  dx = 0.05_sp
  ShapeLevelSet => Shapesphere
  icc = 25
  jcc = 75
  kcc = 25
  Do k = 1, 20
    Do j = 1, 20
      Do i = 1, 20
        flag = Inout(xc(i), xc(j), xc(k), dx/2.0_sp, dx/2.0_sp, dx/2.0_sp)
        If (flag .eq. 1) Then
          sph(i,j,k) = 1.0
        Elseif (flag .eq. 3) Then
          sph(i,j,k) = 0.0
        Elseif (flag .eq. 2) Then
          gr%level = 1
          gr%maxlevel = 5
          gr%xc(1) = xc(i)
          gr%xc(2) = xc(j)
          gr%xc(3) = xc(k)
          gr%dx(1) = dx
          gr%dx(2) = dx
          gr%dx(3) = dx
          Call VolumeOctree(gr)
          sph(i,j,k) = gr%vof
        EndIf
        Phi(icc+i  , jcc+j  , kcc+k  ) = sph(i,j,k)
        Phi(icc+1-i, jcc+j  , kcc+k  ) = sph(i,j,k)
        Phi(icc+i  , jcc+1-j, kcc+k  ) = sph(i,j,k)
        Phi(icc+1-i, jcc+1-j, kcc+k  ) = sph(i,j,k)
        Phi(icc+i  , jcc+j  , kcc+1-k) = sph(i,j,k)
        Phi(icc+1-i, jcc+j  , kcc+1-k) = sph(i,j,k)
        Phi(icc+i  , jcc+1-j, kcc+1-k) = sph(i,j,k)
        Phi(icc+1-i, jcc+1-j, kcc+1-k) = sph(i,j,k)
      End Do
    End Do
  End Do
  ! End hollow Sphere
  !----------------------------------


  !----------------------------------
  ! Tilt hollow cube
  Do k = 1, 20
    xc(k) = dble(k) / dble(20) - 0.025_sp
  End Do
  dx = 0.05_sp
  ShapeLevelSet => ShapeTCube1
  Do k = 1, 20
    Do j = 1, 20
      Do i = 1, 20
        flag = Inout(xc(i), xc(j), xc(k), dx/2.0_sp, dx/2.0_sp, dx/2.0_sp)
        If (flag .eq. 1) Then
          tcube1(i,j,k) = 1.0
        Elseif (flag .eq. 3) Then
          tcube1(i,j,k) = 0.0
        Elseif (flag .eq. 2) Then
          gr%level = 1
          gr%maxlevel = 5
          gr%xc(1) = xc(i)
          gr%xc(2) = xc(j)
          gr%xc(3) = xc(k)
          gr%dx(1) = dx
          gr%dx(2) = dx
          gr%dx(3) = dx
          Call VolumeOctree(gr)
          tcube1(i,j,k) = gr%vof
        EndIf
      End Do
    End Do
  End Do
  ShapeLevelSet => ShapeTCube2
  Do k = 1, 20
    Do j = 1, 20
      Do i = 1, 20
        flag = Inout(xc(i), xc(j), xc(k), dx/2.0_sp, dx/2.0_sp, dx/2.0_sp)
        If (flag .eq. 1) Then
          tcube2(i,j,k) = 1.0
        Elseif (flag .eq. 3) Then
          tcube2(i,j,k) = 0.0
        Elseif (flag .eq. 2) Then
          gr%level = 1
          gr%maxlevel = 5
          gr%xc(1) = xc(i)
          gr%xc(2) = xc(j)
          gr%xc(3) = xc(k)
          gr%dx(1) = dx
          gr%dx(2) = dx
          gr%dx(3) = dx
          Call VolumeOctree(gr)
          tcube2(i,j,k) = gr%vof
        EndIf
      End Do
    End Do
  End Do

  icc = 75
  jcc = 25
  kcc = 25
  Do k = 1, 20
    Do j = 1, 20
      Do i = 1, 20
        Phi(icc+i  , jcc+j  , kcc+k  ) = tcube1(i,j,k) - tcube2(i,j,k)
        Phi(icc+1-i, jcc+j  , kcc+k  ) = tcube1(i,j,k) - tcube2(i,j,k)
        Phi(icc+i  , jcc+1-j, kcc+k  ) = tcube1(i,j,k) - tcube2(i,j,k)
        Phi(icc+1-i, jcc+1-j, kcc+k  ) = tcube1(i,j,k) - tcube2(i,j,k)
        Phi(icc+i  , jcc+j  , kcc+1-k) = tcube1(i,j,k) - tcube2(i,j,k)
        Phi(icc+1-i, jcc+j  , kcc+1-k) = tcube1(i,j,k) - tcube2(i,j,k)
        Phi(icc+i  , jcc+1-j, kcc+1-k) = tcube1(i,j,k) - tcube2(i,j,k)
        Phi(icc+1-i, jcc+1-j, kcc+1-k) = tcube1(i,j,k) - tcube2(i,j,k)
      End Do
    End Do
  End Do
  ! End tilt hollow cube
  !----------------------------------


  !----------------------------------
  ! Tilt hollow cube
  Deallocate(xc)
  Allocate(xc(40))
  Do k = 1, 40
    xc(k) = dble(k) / dble(20) - 0.025_sp
  End Do
  dx = 0.05_sp
  ShapeLevelSet => ShapeLA1
  Do k = 1, 20
    Do j = 1, 20
      Do i = 1, 40
        flag = Inout(xc(i), xc(j), xc(k), dx/2.0_sp, dx/2.0_sp, dx/2.0_sp)
        If (flag .eq. 1) Then
          LA1(i,j,k) = 1.0
        Elseif (flag .eq. 3) Then
          LA1(i,j,k) = 0.0
        Elseif (flag .eq. 2) Then
          gr%level = 1
          gr%maxlevel = 5
          gr%xc(1) = xc(i)
          gr%xc(2) = xc(j)
          gr%xc(3) = xc(k)
          gr%dx(1) = dx
          gr%dx(2) = dx
          gr%dx(3) = dx
          Call VolumeOctree(gr)
          LA1(i,j,k) = gr%vof
        EndIf
      End Do
    End Do
  End Do

  icc = 55
  jcc = 75
  kcc = 25
  Do k = 1, 20
    Do j = 1, 20
      Do i = 1, 40
        Phi(icc+i, jcc+j  , kcc+k  ) = LA1(i,j,k)
        Phi(icc+i, jcc+1-j, kcc+k  ) = LA1(i,j,k)
        Phi(icc+i, jcc+j  , kcc+1-k) = LA1(i,j,k)
        Phi(icc+i, jcc+1-j, kcc+1-k) = LA1(i,j,k)
      End Do
    End Do
  End Do
  ! End Lettter A
  !----------------------------------

  Block 
    Integer :: nn
    Character(80) :: data_name
    data_name = 'test'
    print *, n_vars
    do nn = 1, n_vars
      Select Case(Trim(h5_output_field(nn)%groupname))
      Case('phi')
        Call HDF5WriteData(h5_output_field(nn), phi,data_name)
      Case('u')
        Call HDF5WriteData(h5_output_field(nn), u,data_name)
      Case('v')
        Call HDF5WriteData(h5_output_field(nn), v,data_name)
      Case('w')
        Call HDF5WriteData(h5_output_field(nn), w,data_name)
      Case('cx')
        Call HDF5WriteData(h5_output_field(nn), cx,data_name)
      Case('cy')
        Call HDF5WriteData(h5_output_field(nn), cy,data_name)
      Case('cz')
        Call HDF5WriteData(h5_output_field(nn), cz,data_name)
      End Select
    end do
  End Block

  Call Visual3DContour(f1=phi)


  End Subroutine test2

