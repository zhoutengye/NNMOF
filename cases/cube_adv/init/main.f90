#include "param.h"
# if defined TEST1
# define test_case test1
# elif defined TEST2
# define test_case InitVolume
# elif defined TEST3
# define test_case InitVolumecentroid1
# elif defined TEST4
# define test_case InitVolumecentroid2
# endif
# define epsc 1.0e-12
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
  Integer :: i

  ! Call Init(inputfield=.false.)

  Type(Octree) :: gr

  do i = 1,15
  gr%level = 1
  gr%maxlevel = i
  gr%xc(1) = 0.5_sp
  gr%xc(2) = 0.5_sp
  gr%xc(3) = 0.5_sp
  gr%dx(1) = 1.0_sp
  gr%dx(2) = 1.0_sp
  gr%dx(3) = 1.0_sp
  ShapeLevelSet => ShapeTest

  Call VolumeCentroidOctree(gr)

  print *, i, gr%vof, gr%centroid
  End do
  ! print *, 'tt'

  ! gr%level = 1
  ! gr%maxlevel = 5
  ! gr%xc(1) = 0.5_sp
  ! gr%xc(2) = 0.5_sp
  ! gr%xc(3) = 0.5_sp
  ! gr%dx(1) = 1.0_sp
  ! gr%dx(2) = 1.0_sp
  ! gr%dx(3) = 1.0_sp
  ! ShapeLevelSet => ShapeTest
  ! Call VolumeOctree(gr)
  ! print *, gr%vof

  ! Call MPI_FINALIZE(ierr)

end Subroutine test1

Subroutine InitVolume
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
  Integer :: Octreelevel = 5

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
        flag = Inout(xc(i), xc(j), xc(k), dx, dx, dx)
        If (flag .eq. 1) Then
          cube1(i,j,k) = 1.0
        Elseif (flag .eq. 2) Then
          cube1(i,j,k) = 0.0
        Elseif (flag .eq. 3) Then
          gr%level = 1
          gr%maxlevel = Octreelevel
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
        flag = Inout(xc(i), xc(j), xc(k), dx, dx, dx)
        If (flag .eq. 1) Then
          cube2(i,j,k) = 1.0
        Elseif (flag .eq. 2) Then
          cube2(i,j,k) = 0.0
        Elseif (flag .eq. 3) Then
          gr%level = 1
          gr%maxlevel = Octreelevel
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
        flag = Inout(xc(i), xc(j), xc(k), dx, dx, dx)
        If (flag .eq. 1) Then
          sph(i,j,k) = 1.0
        Elseif (flag .eq. 2) Then
          sph(i,j,k) = 0.0
        Elseif (flag .eq. 3) Then
          gr%level = 1
          gr%maxlevel = OctreeLevel
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
        flag = Inout(xc(i), xc(j), xc(k), dx, dx, dx)
        If (flag .eq. 1) Then
          tcube1(i,j,k) = 1.0
        Elseif (flag .eq. 2) Then
          tcube1(i,j,k) = 0.0
        Elseif (flag .eq. 3) Then
          gr%level = 1
          gr%maxlevel = OctreeLevel
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
        Elseif (flag .eq. 2) Then
          tcube2(i,j,k) = 0.0
        Elseif (flag .eq. 3) Then
          gr%level = 1
          gr%maxlevel = OctreeLevel
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
        flag = Inout(xc(i), xc(j), xc(k), dx, dx, dx)
        If (flag .eq. 1) Then
          LA1(i,j,k) = 1.0
        Elseif (flag .eq. 2) Then
          LA1(i,j,k) = 0.0
        Elseif (flag .eq. 3) Then
          gr%level = 1
          gr%maxlevel = OctreeLevel
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
    data_name = 'init'
    print *, n_vars
    u = 1.0_sp
    v = 1.0_sp
    w = 1.0_sp
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

  Call MPI_FINALIZE(ierr)

  End Subroutine InitVolume

Subroutine InitVolumeCentroid1
  Use ModTools
  Use InitSH
  Implicit None
  Real(sp) :: sph(20,20,20)
  Real(sp) :: cube1(20,20,20)
  Real(sp) :: cube2(20,20,20)
  Real(sp) :: tcube1(20,20,20)
  Real(sp) :: tcube2(20,20,20)
  Real(sp) :: LA1(40,20,20)
  Real(sp), Dimension(20,20,20) :: cxsph, cysph, czsph
  Real(sp), Dimension(20,20,20) :: cxcube, cycube, czcube
  Real(sp), Dimension(20,20,20) :: cxcube1, cycube1, czcube1
  Real(sp), Dimension(20,20,20) :: cxcube2, cycube2, czcube2
  Real(sp), Dimension(20,20,20) :: cxtcube, cytcube, cztcube
  Real(sp), Dimension(20,20,20) :: cxtcube1, cytcube1, cztcube1
  Real(sp), Dimension(20,20,20) :: cxtcube2, cytcube2, cztcube2
  Real(sp), Dimension(40,20,20) :: cxLA1, cyLA1, czLA1
  Real(sp), Allocatable :: xc(:)
  Type(Octree) :: gr
  Integer :: i,j,k
  Integer :: icc,jcc,kcc
  Integer :: flag
  Real(sp) :: dx
  Integer :: Octreelevel = 5

  Call Init(inputfield=.false.)

  cx = 0.0_sp
  cy = 0.0_sp
  cz = 0.0_sp

  ! ! ---------------------------------
  ! ! Hollow cube
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
          cxcube1(i,j,k) = 0.5_sp
          cycube1(i,j,k) = 0.5_sp
          czcube1(i,j,k) = 0.5_sp
        Elseif (flag .eq. 2) Then
          cube1(i,j,k) = 0.0
        Elseif (flag .eq. 3) Then
          gr%level = 1
          gr%maxlevel = Octreelevel
          gr%xc(1) = xc(i)
          gr%xc(2) = xc(j)
          gr%xc(3) = xc(k)
          gr%dx(1) = dx
          gr%dx(2) = dx
          gr%dx(3) = dx
          Call VolumeCentroidOctree(gr)
          cube1(i,j,k) = gr%vof
          cxcube1(i,j,k) = gr%centroid(1) / dx
          cycube1(i,j,k) = gr%centroid(2) / dx
          czcube1(i,j,k) = gr%centroid(3) / dx
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
          cxcube2(i,j,k) = 0.5_sp
          cycube2(i,j,k) = 0.5_sp
          czcube2(i,j,k) = 0.5_sp
        Elseif (flag .eq. 2) Then
          cube2(i,j,k) = 0.0
        Elseif (flag .eq. 3) Then
          gr%level = 1
          gr%maxlevel = Octreelevel
          gr%xc(1) = xc(i)
          gr%xc(2) = xc(j)
          gr%xc(3) = xc(k)
          gr%dx(1) = dx
          gr%dx(2) = dx
          gr%dx(3) = dx
          Call VolumeCentroidOctree(gr)
          cube2(i,j,k) = gr%vof
          cxcube2(i,j,k) = gr%centroid(1) / dx
          cycube2(i,j,k) = gr%centroid(2) / dx
          czcube2(i,j,k) = gr%centroid(3) / dx
        EndIf
      End Do
    End Do
  End Do

  Do k = 1, 20
    Do j = 1, 20
      Do i = 1, 20
        if ( (cube1(i,j,k) - cube2(i,j,k)) .lt. epsc) Then
          cxcube(i,j,k) = 0.0_sp
          cycube(i,j,k) = 0.0_sp
          czcube(i,j,k) = 0.0_sp
        elseif ( (cube1(i,j,k) - cube2(i,j,k)) .gt. 1.0_sp-epsc) Then
          cxcube(i,j,k) = 0.5_sp
          cycube(i,j,k) = 0.5_sp
          czcube(i,j,k) = 0.5_sp
        else
          ! print *, 'yes', i,j,k, (cube1(i,j,k) - cube2(i,j,k))
          cxcube(i,j,k) = ( cube1(i,j,k)*cxcube1(i,j,k) - cube2(i,j,k)*cxcube2(i,j,k) ) / (cube1(i,j,k) - cube2(i,j,k))
          cycube(i,j,k) = cube1(i,j,k)*cycube1(i,j,k) - cube2(i,j,k)*cycube2(i,j,k)
          czcube(i,j,k) = cube1(i,j,k)*czcube1(i,j,k) - cube2(i,j,k)*czcube2(i,j,k)
          if ( cxcube(i,j,k) .lt. 0.0_sp .or. cxcube(i,j,k) .gt. 0.5_sp) print *, 'Hollow cube, cx', i, j, k
          if ( cycube(i,j,k) .lt. 0.0_sp .or. cycube(i,j,k) .gt. 0.5_sp) print *, 'Hollow cube, cx', i, j, k
          if ( czcube(i,j,k) .lt. 0.0_sp .or. czcube(i,j,k) .gt. 0.5_sp) print *, 'Hollow cube, cx', i, j, k
          cxcube(i,j,k) = min(0.0_sp, max(0.5,cxcube(i,j,k)))
          cycube(i,j,k) = min(0.0_sp, max(0.5,cycube(i,j,k)))
          czcube(i,j,k) = min(0.0_sp, max(0.5,czcube(i,j,k)))
        endif
      End Do
    End Do
  End Do

  icc = 25
  jcc = 25
  kcc = 25
  Do k = 1, 20
    Do j = 1, 20
      Do i = 1, 20
        ! phi
        Phi(icc+i  , jcc+j  , kcc+k  ) = cube1(i,j,k) - cube2(i,j,k)
        Phi(icc+1-i, jcc+j  , kcc+k  ) = cube1(i,j,k) - cube2(i,j,k)
        Phi(icc+i  , jcc+1-j, kcc+k  ) = cube1(i,j,k) - cube2(i,j,k)
        Phi(icc+1-i, jcc+1-j, kcc+k  ) = cube1(i,j,k) - cube2(i,j,k)
        Phi(icc+i  , jcc+j  , kcc+1-k) = cube1(i,j,k) - cube2(i,j,k)
        Phi(icc+1-i, jcc+j  , kcc+1-k) = cube1(i,j,k) - cube2(i,j,k)
        Phi(icc+i  , jcc+1-j, kcc+1-k) = cube1(i,j,k) - cube2(i,j,k)
        Phi(icc+1-i, jcc+1-j, kcc+1-k) = cube1(i,j,k) - cube2(i,j,k)
        if(cxcube(i,j,k) .gt. epsc) then
        ! cx
          cx(icc+i  , jcc+j  , kcc+k  ) = cxcube(i,j,k)
          cx(icc+1-i, jcc+j  , kcc+k  ) = 1.0_sp - cxcube(i,j,k)
          cx(icc+i  , jcc+1-j, kcc+k  ) = cxcube(i,j,k)
          cx(icc+1-i, jcc+1-j, kcc+k  ) = 1.0_sp - cxcube(i,j,k)
          cx(icc+i  , jcc+j  , kcc+1-k) = cxcube(i,j,k)
          cx(icc+1-i, jcc+j  , kcc+1-k) = 1.0_sp - cxcube(i,j,k)
          cx(icc+i  , jcc+1-j, kcc+1-k) = cxcube(i,j,k)
          cx(icc+1-i, jcc+1-j, kcc+1-k) = 1.0_sp - cxcube(i,j,k)
        endif
        ! cy
        if(cycube(i,j,k) .gt. epsc) then
          cy(icc+i  , jcc+j  , kcc+k  ) = cycube(i,j,k)
          cy(icc+1-i, jcc+j  , kcc+k  ) = cycube(i,j,k)
          cy(icc+i  , jcc+1-j, kcc+k  ) = 1.0_sp - cycube(i,j,k)
          cy(icc+1-i, jcc+1-j, kcc+k  ) = 1.0_sp - cycube(i,j,k)
          cy(icc+i  , jcc+j  , kcc+1-k) = cycube(i,j,k)
          cy(icc+1-i, jcc+j  , kcc+1-k) = cycube(i,j,k)
          cy(icc+i  , jcc+1-j, kcc+1-k) = 1.0_sp - cycube(i,j,k)
          cy(icc+1-i, jcc+1-j, kcc+1-k) = 1.0_sp - cycube(i,j,k)
        endif
        ! cz
        if(czcube(i,j,k) .gt. epsc) then
          cz(icc+i  , jcc+j  , kcc+k  ) = czcube(i,j,k)
          cz(icc+1-i, jcc+j  , kcc+k  ) = czcube(i,j,k)
          cz(icc+i  , jcc+1-j, kcc+k  ) = czcube(i,j,k)
          cz(icc+1-i, jcc+1-j, kcc+k  ) = czcube(i,j,k)
          cz(icc+i  , jcc+j  , kcc+1-k) = 1.0_sp - czcube(i,j,k)
          cz(icc+1-i, jcc+j  , kcc+1-k) = 1.0_sp - czcube(i,j,k)
          cz(icc+i  , jcc+1-j, kcc+1-k) = 1.0_sp - czcube(i,j,k)
          cz(icc+1-i, jcc+1-j, kcc+1-k) = 1.0_sp - czcube(i,j,k)
        endif
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
  cxsph = 0.0_sp
  cysph = 0.0_sp
  czsph = 0.0_sp
  Do k = 1, 20
    Do j = 1, 20
      Do i = 1, 20
        flag = Inout(xc(i), xc(j), xc(k), dx, dx, dx)
        If (flag .eq. 1) Then
          sph(i,j,k) = 1.0_sp
          cxsph(i,j,k) = 0.5_sp
          cysph(i,j,k) = 0.5_sp
          czsph(i,j,k) = 0.5_sp
        Elseif (flag .eq. 2) Then
          sph(i,j,k) = 0.0
        Elseif (flag .eq. 3) Then
          gr%level = 1
          gr%maxlevel = OctreeLevel
          gr%xc(1) = xc(i)
          gr%xc(2) = xc(j)
          gr%xc(3) = xc(k)
          gr%dx(1) = dx
          gr%dx(2) = dx
          gr%dx(3) = dx
          Call VolumeCentroidOctree(gr)
          sph(i,j,k) = gr%vof
          cxsph(i,j,k) = gr%centroid(1)
          cysph(i,j,k) = gr%centroid(2)
          czsph(i,j,k) = gr%centroid(3)
        EndIf
      End Do
    End Do
  End Do

  Do k = 1, 20
    Do j = 1, 20
      Do i = 1, 20
        ! phi
        Phi(icc+i  , jcc+j  , kcc+k  ) = sph(i,j,k)
        Phi(icc+1-i, jcc+j  , kcc+k  ) = sph(i,j,k)
        Phi(icc+i  , jcc+1-j, kcc+k  ) = sph(i,j,k)
        Phi(icc+1-i, jcc+1-j, kcc+k  ) = sph(i,j,k)
        Phi(icc+i  , jcc+j  , kcc+1-k) = sph(i,j,k)
        Phi(icc+1-i, jcc+j  , kcc+1-k) = sph(i,j,k)
        Phi(icc+i  , jcc+1-j, kcc+1-k) = sph(i,j,k)
        Phi(icc+1-i, jcc+1-j, kcc+1-k) = sph(i,j,k)
        if(sph(i,j,k) .gt. epsc) then
          cx(icc+i  , jcc+j  , kcc+k  ) = cxsph(i,j,k)
          cx(icc+1-i, jcc+j  , kcc+k  ) = 1.0_sp - cxsph(i,j,k)
          cx(icc+i  , jcc+1-j, kcc+k  ) = cxsph(i,j,k)
          cx(icc+1-i, jcc+1-j, kcc+k  ) = 1.0_sp - cxsph(i,j,k)
          cx(icc+i  , jcc+j  , kcc+1-k) = cxsph(i,j,k)
          cx(icc+1-i, jcc+j  , kcc+1-k) = 1.0_sp - cxsph(i,j,k)
          cx(icc+i  , jcc+1-j, kcc+1-k) = cxsph(i,j,k)
          cx(icc+1-i, jcc+1-j, kcc+1-k) = 1.0_sp - cxsph(i,j,k)
          cy(icc+i  , jcc+j  , kcc+k  ) = cysph(i,j,k)
          cy(icc+1-i, jcc+j  , kcc+k  ) = cysph(i,j,k)
          cy(icc+i  , jcc+1-j, kcc+k  ) = 1.0_sp - cysph(i,j,k)
          cy(icc+1-i, jcc+1-j, kcc+k  ) = 1.0_sp - cysph(i,j,k)
          cy(icc+i  , jcc+j  , kcc+1-k) = cysph(i,j,k)
          cy(icc+1-i, jcc+j  , kcc+1-k) = cysph(i,j,k)
          cy(icc+i  , jcc+1-j, kcc+1-k) = 1.0_sp - cysph(i,j,k)
          cy(icc+1-i, jcc+1-j, kcc+1-k) = 1.0_sp - cysph(i,j,k)
          cz(icc+i  , jcc+j  , kcc+k  ) = czsph(i,j,k)
          cz(icc+1-i, jcc+j  , kcc+k  ) = czsph(i,j,k)
          cz(icc+i  , jcc+1-j, kcc+k  ) = czsph(i,j,k)
          cz(icc+1-i, jcc+1-j, kcc+k  ) = czsph(i,j,k)
          cz(icc+i  , jcc+j  , kcc+1-k) = 1.0_sp - czsph(i,j,k)
          cz(icc+1-i, jcc+j  , kcc+1-k) = 1.0_sp - czsph(i,j,k)
          cz(icc+i  , jcc+1-j, kcc+1-k) = 1.0_sp - czsph(i,j,k)
          cz(icc+1-i, jcc+1-j, kcc+1-k) = 1.0_sp - czsph(i,j,k)
        endif
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
          tcube1(i,j,k) = 1.0_sp
          cxtcube1(i,j,k) = 0.5_sp
          cytcube1(i,j,k) = 0.5_sp
          cztcube1(i,j,k) = 0.5_sp
        Elseif (flag .eq. 2) Then
          tcube1(i,j,k) = 0.0
        Elseif (flag .eq. 3) Then
          gr%level = 1
          gr%maxlevel = OctreeLevel
          gr%xc(1) = xc(i)
          gr%xc(2) = xc(j)
          gr%xc(3) = xc(k)
          gr%dx(1) = dx
          gr%dx(2) = dx
          gr%dx(3) = dx
          cxtcube1(i,j,k) = gr%centroid(1)
          cytcube1(i,j,k) = gr%centroid(2)
          cztcube1(i,j,k) = gr%centroid(3)
          if ( cxtcube1(i,j,k) .gt.0.5_sp) print *,cxtcube1(i,j,k)
          Call VolumeCentroidOctree(gr)
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
          tcube2(i,j,k) = 1.0_sp
          cxtcube2(i,j,k) = 0.5_sp
          cytcube2(i,j,k) = 0.5_sp
          cztcube2(i,j,k) = 0.5_sp
        Elseif (flag .eq. 2) Then
          tcube2(i,j,k) = 0.0
        Elseif (flag .eq. 3) Then
          gr%level = 1
          gr%maxlevel = OctreeLevel
          gr%xc(1) = xc(i)
          gr%xc(2) = xc(j)
          gr%xc(3) = xc(k)
          gr%dx(1) = dx
          gr%dx(2) = dx
          gr%dx(3) = dx
          cxtcube2(i,j,k) = gr%centroid(1)
          cytcube2(i,j,k) = gr%centroid(2)
          cztcube2(i,j,k) = gr%centroid(3)
          Call VolumeCentroidOctree(gr)
          tcube2(i,j,k) = gr%vof
        EndIf
      End Do
    End Do
  End Do


  ! Do k = 1, 20
  !   Do j = 1, 20
  !     Do i = 1, 20
  !       if ( (tcube1(i,j,k) - tcube2(i,j,k)) .lt. epsc) Then
  !         cxtcube(i,j,k) = 0.0_sp
  !         cytcube(i,j,k) = 0.0_sp
  !         cztcube(i,j,k) = 0.0_sp
  !       elseif ( (tcube1(i,j,k) - tcube2(i,j,k)) .gt. 1.0_sp-epsc) Then
  !         cxtcube(i,j,k) = 0.5_sp
  !         cytcube(i,j,k) = 0.5_sp
  !         cztcube(i,j,k) = 0.5_sp
  !       else
  !         cxtcube(i,j,k) = ( tcube1(i,j,k)*cxtcube1(i,j,k) - tcube2(i,j,k)*cxtcube2(i,j,k) ) / (tcube1(i,j,k) - tcube2(i,j,k))
  !         cytcube(i,j,k) = ( tcube1(i,j,k)*cytcube1(i,j,k) - tcube2(i,j,k)*cytcube2(i,j,k) ) / (tcube1(i,j,k) - tcube2(i,j,k))
  !         cztcube(i,j,k) = ( tcube1(i,j,k)*cztcube1(i,j,k) - tcube2(i,j,k)*cztcube2(i,j,k) ) / (tcube1(i,j,k) - tcube2(i,j,k))
  !         cxtcube(i,j,k) = min(0.0_sp, max(1.0,cxtcube(i,j,k)))
  !         cytcube(i,j,k) = min(0.0_sp, max(1.0,cytcube(i,j,k)))
  !         cztcube(i,j,k) = min(0.0_sp, max(1.0,cztcube(i,j,k)))
  !       endif
  !     End Do
  !   End Do
  ! End Do

  icc = 75
  jcc = 25
  kcc = 25
  Do k = 1, 20
    Do j = 1, 20
      Do i = 1, 20
        ! phi
        Phi(icc+i  , jcc+j  , kcc+k  ) = tcube1(i,j,k) - tcube2(i,j,k)
        Phi(icc+1-i, jcc+j  , kcc+k  ) = tcube1(i,j,k) - tcube2(i,j,k)
        Phi(icc+i  , jcc+1-j, kcc+k  ) = tcube1(i,j,k) - tcube2(i,j,k)
        Phi(icc+1-i, jcc+1-j, kcc+k  ) = tcube1(i,j,k) - tcube2(i,j,k)
        Phi(icc+i  , jcc+j  , kcc+1-k) = tcube1(i,j,k) - tcube2(i,j,k)
        Phi(icc+1-i, jcc+j  , kcc+1-k) = tcube1(i,j,k) - tcube2(i,j,k)
        Phi(icc+i  , jcc+1-j, kcc+1-k) = tcube1(i,j,k) - tcube2(i,j,k)
        Phi(icc+1-i, jcc+1-j, kcc+1-k) = tcube1(i,j,k) - tcube2(i,j,k)
        cx(icc+i  , jcc+j  , kcc+k  ) = cxtcube2(i,j,k)
        cy(icc+i  , jcc+j  , kcc+k  ) = cytcube2(i,j,k)
        cz(icc+i  , jcc+j  , kcc+k  ) = cxtcube2(i,j,k)
        if(tcube1(i,j,k) - tcube2(i,j,k) .gt. epsc) then
        ! cx
          ! cx(icc+i  , jcc+j  , kcc+k  ) = cxtcube2(i,j,k)
          ! cx(icc+1-i, jcc+j  , kcc+k  ) = 1.0_sp - cxtcube(i,j,k)
          ! cx(icc+i  , jcc+1-j, kcc+k  ) = cxtcube(i,j,k)
          ! cx(icc+1-i, jcc+1-j, kcc+k  ) = 1.0_sp - cxtcube(i,j,k)
          ! cx(icc+i  , jcc+j  , kcc+1-k) = cxtcube(i,j,k)
          ! cx(icc+1-i, jcc+j  , kcc+1-k) = 1.0_sp - cxtcube(i,j,k)
          ! cx(icc+i  , jcc+1-j, kcc+1-k) = cxtcube(i,j,k)
          ! cx(icc+1-i, jcc+1-j, kcc+1-k) = 1.0_sp - cxtcube(i,j,k)
        ! endif
        ! ! cy
        ! if(cytcube(i,j,k) .gt. epsc) then
          ! cy(icc+i  , jcc+j  , kcc+k  ) = cytcube(i,j,k)
          ! cy(icc+1-i, jcc+j  , kcc+k  ) = cytcube(i,j,k)
          ! cy(icc+i  , jcc+1-j, kcc+k  ) = 1.0_sp - cytcube(i,j,k)
          ! cy(icc+1-i, jcc+1-j, kcc+k  ) = 1.0_sp - cytcube(i,j,k)
          ! cy(icc+i  , jcc+j  , kcc+1-k) = cytcube(i,j,k)
          ! cy(icc+1-i, jcc+j  , kcc+1-k) = cytcube(i,j,k)
          ! cy(icc+i  , jcc+1-j, kcc+1-k) = 1.0_sp - cytcube(i,j,k)
          ! cy(icc+1-i, jcc+1-j, kcc+1-k) = 1.0_sp - cytcube(i,j,k)
        ! endif
        ! ! cz
        ! if(cztcube(i,j,k) .gt. epsc) then
          ! cz(icc+i  , jcc+j  , kcc+k  ) = cztcube(i,j,k)
          ! cz(icc+1-i, jcc+j  , kcc+k  ) = cztcube(i,j,k)
          ! cz(icc+i  , jcc+1-j, kcc+k  ) = cztcube(i,j,k)
          ! cz(icc+1-i, jcc+1-j, kcc+k  ) = cztcube(i,j,k)
          ! cz(icc+i  , jcc+j  , kcc+1-k) = 1.0_sp - cztcube(i,j,k)
          ! cz(icc+1-i, jcc+j  , kcc+1-k) = 1.0_sp - cztcube(i,j,k)
          ! cz(icc+i  , jcc+1-j, kcc+1-k) = 1.0_sp - cztcube(i,j,k)
          ! cz(icc+1-i, jcc+1-j, kcc+1-k) = 1.0_sp - cztcube(i,j,k)
        endif
      End Do
    End Do
  End Do
  Call Visual3DContour(f1=cx)
  Call Visual3DContour(f1=cy)
  Call Visual3DContour(f1=cz)


  ! Call Visual3DContour(f1=cx)
  ! End tilt hollow cube
  !----------------------------------


  ! !----------------------------------
  ! ! Tilt hollow cube
  ! Deallocate(xc)
  ! Allocate(xc(40))
  ! Do k = 1, 40
  !   xc(k) = dble(k) / dble(20) - 0.025_sp
  ! End Do
  ! dx = 0.05_sp
  ! ShapeLevelSet => ShapeLA1
  ! Do k = 1, 20
  !   Do j = 1, 20
  !     Do i = 1, 40
  !       flag = Inout(xc(i), xc(j), xc(k), dx/2.0_sp, dx/2.0_sp, dx/2.0_sp)
  !       If (flag .eq. 1) Then
  !         LA1(i,j,k) = 1.0
  !       Elseif (flag .eq. 3) Then
  !         LA1(i,j,k) = 0.0
  !       Elseif (flag .eq. 2) Then
  !         gr%level = 1
  !         gr%maxlevel = OctreeLevel
  !         gr%xc(1) = xc(i)
  !         gr%xc(2) = xc(j)
  !         gr%xc(3) = xc(k)
  !         gr%dx(1) = dx
  !         gr%dx(2) = dx
  !         gr%dx(3) = dx
  !         Call VolumeCentroidOctree(gr)
  !         LA1(i,j,k) = gr%vof
  !       EndIf
  !     End Do
  !   End Do
  ! End Do

  ! icc = 55
  ! jcc = 75
  ! kcc = 25
  ! Do k = 1, 20
  !   Do j = 1, 20
  !     Do i = 1, 40
  !       Phi(icc+i, jcc+j  , kcc+k  ) = LA1(i,j,k)
  !       Phi(icc+i, jcc+1-j, kcc+k  ) = LA1(i,j,k)
  !       Phi(icc+i, jcc+j  , kcc+1-k) = LA1(i,j,k)
  !       Phi(icc+i, jcc+1-j, kcc+1-k) = LA1(i,j,k)
  !     End Do
  !   End Do
  ! End Do
  ! ! End Lettter A
  ! !----------------------------------

  Block 
    Integer :: nn
    Character(80) :: data_name
    data_name = 'init'
    print *, n_vars
    u = 1.0_sp
    v = 1.0_sp
    w = 1.0_sp
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

  ! Call Visual3DContour(f1=cx)

  Call MPI_FINALIZE(ierr)

  End Subroutine InitVolumeCentroid1

Subroutine InitVolumeCentroid2
  Use ModTools
  Use InitSH
  Implicit None
  Real(sp), Dimension(0:101,0:101,0:51) :: LS
  Real(sp), Dimension(0:101,0:101,0:51) :: LSsph
  Real(sp), Dimension(0:101,0:101,0:51) :: LScube1
  Real(sp), Dimension(0:101,0:101,0:51) :: LScube2
  Real(sp), Dimension(0:101,0:101,0:51) :: LStcube1
  Real(sp), Dimension(0:101,0:101,0:51) :: LStcube2
  Real(sp), Dimension(0:101,0:101,0:51) :: LSla
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
  Integer :: Octreelevel = 5

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
        flag = Inout(xc(i), xc(j), xc(k), dx, dx, dx)
        If (flag .eq. 1) Then
          cube1(i,j,k) = 1.0
        Elseif (flag .eq. 2) Then
          cube1(i,j,k) = 0.0
        Elseif (flag .eq. 3) Then
          gr%level = 1
          gr%maxlevel = Octreelevel
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
        flag = Inout(xc(i), xc(j), xc(k), dx, dx, dx)
        If (flag .eq. 1) Then
          cube2(i,j,k) = 1.0
        Elseif (flag .eq. 2) Then
          cube2(i,j,k) = 0.0
        Elseif (flag .eq. 3) Then
          gr%level = 1
          gr%maxlevel = Octreelevel
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
        flag = Inout(xc(i), xc(j), xc(k), dx, dx, dx)
        If (flag .eq. 1) Then
          sph(i,j,k) = 1.0
        Elseif (flag .eq. 2) Then
          sph(i,j,k) = 0.0
        Elseif (flag .eq. 3) Then
          gr%level = 1
          gr%maxlevel = OctreeLevel
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
        flag = Inout(xc(i), xc(j), xc(k), dx, dx, dx)
        If (flag .eq. 1) Then
          tcube1(i,j,k) = 1.0
        Elseif (flag .eq. 2) Then
          tcube1(i,j,k) = 0.0
        Elseif (flag .eq. 3) Then
          gr%level = 1
          gr%maxlevel = OctreeLevel
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
        Elseif (flag .eq. 2) Then
          tcube2(i,j,k) = 0.0
        Elseif (flag .eq. 3) Then
          gr%level = 1
          gr%maxlevel = OctreeLevel
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
        flag = Inout(xc(i), xc(j), xc(k), dx, dx, dx)
        If (flag .eq. 1) Then
          LA1(i,j,k) = 1.0
        Elseif (flag .eq. 2) Then
          LA1(i,j,k) = 0.0
        Elseif (flag .eq. 3) Then
          gr%level = 1
          gr%maxlevel = OctreeLevel
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

  !----------------------------------
  ! Hollow cube ls
  deallocate(xc)
  allocate(xc(100))
  Do k = 1, 100
    xc(k) = dble(k) - 0.5_sp
  End Do
  Do k = 1, 50
    Do j = 1, 100
      Do i = 1, 100
        LScube1(i,j,k) = GShapeCube1(xc(i),xc(j),xc(k))
        Phi(i,j,k) = heaviside(LScube1(i,j,k),1.0_sp)
        ! LScube2(i,j,k) = GShapeCube2(xc(i),xc(j),xc(k))
        ! Phi(i,j,k) = heaviside(min(LScube1(i,j,k), 1.0_sp),LScube2(i,j,k))
      End Do
    End Do
  End Do
  Call Visual3DContour(f1=phi)
  ! End hollow cube ls
  !----------------------------------

  Block 
    Integer :: nn
    Character(80) :: data_name
    data_name = 'init'
    print *, n_vars
    u = 1.0_sp
    v = 1.0_sp
    w = 1.0_sp
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


  Call MPI_FINALIZE(ierr)

  End Subroutine InitVolumeCentroid2
