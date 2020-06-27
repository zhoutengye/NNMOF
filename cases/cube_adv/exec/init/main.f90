#include "param.h"
# if defined TEST1
# define test_case test1
# elif defined TEST2
# define test_case InitVolume
# elif defined TEST3
# define test_case InitVolumecentroid
# elif defined TEST4
# define test_case InitVolumecentroid
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

  do i = 1,10
    gr%level = 1
    gr%maxlevel = i
    gr%xc(1) = 1.5_sp
    gr%xc(2) = 1.5_sp
    gr%xc(3) = 1.5_sp
    gr%dx(1) = 1.0_sp
    gr%dx(2) = 1.0_sp
    gr%dx(3) = 1.0_sp
    ShapeLevelSet => Shapetest

    ! Call VolumeCentroidOctree(gr)
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
  Real(sp), Allocatable :: xc(:)
  Type(Octree) :: gr
  Integer :: i,j,k
  Integer :: flag
  Real(sp) :: dx
  Integer :: Octreelevel = 3

  Call Init(inputfield=.false.)


  ! ---------------------------------
  ! Hollow cube
  Allocate(xc(100))
  Do k = 1, 100
    xc(k) = dble(k) - 0.5_sp
  End Do
  dx = 1.0_sp
  ShapeLevelSet => GShapeCube
  Do k = 1, 50
    Do j = 1, 50
      Do i = 1, 50
        flag = Inout(xc(i), xc(j), xc(k), dx, dx, dx)
        If (flag .eq. 1) Then
          phi(i,j,k) = 1.0
        Elseif (flag .eq. 2) Then
          phi(i,j,k) = 0.0
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
          phi(i,j,k) = gr%vof
        EndIf
      End Do
    End Do
  End Do
  ! End Hollow cube
  ! ---------------------------------

  !----------------------------------
  ! Hollow Sphere
  ShapeLevelSet => GShapesphere
  Do k = 1, 50
    Do j = 51, 100
      Do i = 1, 50
        flag = Inout(xc(i), xc(j), xc(k), dx, dx, dx)
        If (flag .eq. 1) Then
          phi(i,j,k) = 1.0
        Elseif (flag .eq. 2) Then
          phi(i,j,k) = 0.0
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
          phi(i,j,k) = gr%vof
        EndIf
      End Do
    End Do
  End Do
  ! End hollow Sphere
  !----------------------------------


  !----------------------------------
  ! Tilt hollow cube
  ShapeLevelSet => GShapeTCube
  Do k = 1, 50
    Do j = 1, 50
      Do i = 51, 100
        flag = Inout(xc(i), xc(j), xc(k), dx, dx, dx)
        If (flag .eq. 1) Then
          phi(i,j,k) = 1.0
        Elseif (flag .eq. 2) Then
          phi(i,j,k) = 0.0
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
          phi(i,j,k) = gr%vof
        EndIf
      End Do
    End Do
  End Do
  ! End tilt hollow cube
  !----------------------------------


  !----------------------------------
  ! Tilt letter A
  ShapeLevelSet => GShapeLA
  Do k = 1, 50
    Do j = 51, 100
      Do i = 51, 100
        flag = Inout(xc(i), xc(j), xc(k), dx, dx, dx)
        If (flag .eq. 1) Then
          phi(i,j,k) = 1.0
        Elseif (flag .eq. 2) Then
          phi(i,j,k) = 0.0
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
          phi(i,j,k) = gr%vof
        EndIf
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


Subroutine InitVolumeCentroid
  Use ModTools
  Use InitSH
  Implicit None
  Real(sp), Allocatable :: xc(:)
  Type(Octree) :: gr
  Integer :: i,j,k
  Integer :: flag
  Real(sp) :: dx
  Integer :: Octreelevel = 5

  Call Init(inputfield=.false.)


  ! ---------------------------------
  ! Hollow cube
  Allocate(xc(100))
  Do k = 1, 100
    xc(k) = dble(k) - 0.5_sp
  End Do
  dx = 1.0_sp
  ShapeLevelSet => GShapeCube
  Do k = 1, 50
    Do j = 1, 50
      Do i = 1, 50
        flag = Inout(xc(i), xc(j), xc(k), dx, dx, dx)
        If (flag .eq. 1) Then
          phi(i,j,k) = 1.0_sp
          cx(i,j,k)  = 0.5_sp
          cy(i,j,k)  = 0.5_sp
          cz(i,j,k)  = 0.5_sp
        Elseif (flag .eq. 2) Then
          phi(i,j,k) = 0.0_sp
          cx(i,j,k)  = 0.0_sp
          cy(i,j,k)  = 0.0_sp
          cz(i,j,k)  = 0.0_sp
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
          phi(i,j,k) = gr%vof
          if ( phi(i,j,k) .gt. epsc ) cx(i,j,k)  = gr%centroid(1) + 0.5_sp
          if ( phi(i,j,k) .gt. epsc ) cy(i,j,k)  = gr%centroid(2) + 0.5_sp
          if ( phi(i,j,k) .gt. epsc ) cz(i,j,k)  = gr%centroid(3) + 0.5_sp
        EndIf
      End Do
    End Do
  End Do
  ! End Hollow cube
  ! ---------------------------------

  !----------------------------------
  ! Hollow Sphere
  ShapeLevelSet => GShapesphere
  Do k = 1, 50
    Do j = 51, 100
      Do i = 1, 50
        flag = Inout(xc(i), xc(j), xc(k), dx, dx, dx)
        If (flag .eq. 1) Then
          phi(i,j,k) = 1.0_sp
          cx(i,j,k)  = 0.5_sp
          cy(i,j,k)  = 0.5_sp
          cz(i,j,k)  = 0.5_sp
        Elseif (flag .eq. 2) Then
          phi(i,j,k) = 0.0_sp
          cx(i,j,k)  = 0.0_sp
          cy(i,j,k)  = 0.0_sp
          cz(i,j,k)  = 0.0_sp
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
          phi(i,j,k) = gr%vof
          if ( phi(i,j,k) .gt. epsc ) cx(i,j,k)  = gr%centroid(1) + 0.5_sp
          if ( phi(i,j,k) .gt. epsc ) cy(i,j,k)  = gr%centroid(2) + 0.5_sp
          if ( phi(i,j,k) .gt. epsc ) cz(i,j,k)  = gr%centroid(3) + 0.5_sp
        EndIf
      End Do
    End Do
  End Do
  ! End hollow Sphere
  !----------------------------------


  !----------------------------------
  ! Tilt hollow cube
  ShapeLevelSet => GShapeTCube
  Do k = 1, 50
    Do j = 1, 50
      Do i = 51, 100
        flag = Inout(xc(i), xc(j), xc(k), dx, dx, dx)
        If (flag .eq. 1) Then
          phi(i,j,k) = 1.0_sp
          cx(i,j,k)  = 0.0_sp
          cy(i,j,k)  = 0.0_sp
          cz(i,j,k)  = 0.0_sp
        Elseif (flag .eq. 2) Then
          phi(i,j,k) = 0.0_sp
          cx(i,j,k)  = 0.0_sp
          cy(i,j,k)  = 0.0_sp
          cz(i,j,k)  = 0.0_sp
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
          phi(i,j,k) = gr%vof
          cx(i,j,k)  = gr%centroid(1) + 0.5_sp
          cy(i,j,k)  = gr%centroid(2) + 0.5_sp
          cz(i,j,k)  = gr%centroid(3) + 0.5_sp
        EndIf
      End Do
    End Do
  End Do
  ! End tilt hollow cube
  !----------------------------------


  !----------------------------------
  ! Tilt letter A
  ShapeLevelSet => GShapeLA
  Do k = 1, 50
    Do j = 51, 100
      Do i = 51, 100
        flag = Inout(xc(i), xc(j), xc(k), dx, dx, dx)
        If (flag .eq. 1) Then
          phi(i,j,k) = 1.0_sp
          cx(i,j,k)  = 0.5_sp
          cy(i,j,k)  = 0.5_sp
          cz(i,j,k)  = 0.5_sp
        Elseif (flag .eq. 2) Then
          phi(i,j,k) = 0.0_sp
          cx(i,j,k)  = 0.0_sp
          cy(i,j,k)  = 0.0_sp
          cz(i,j,k)  = 0.0_sp
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
          phi(i,j,k) = gr%vof
           cx(i,j,k)  = gr%centroid(1) + 0.5_sp
           cy(i,j,k)  = gr%centroid(2) + 0.5_sp
           cz(i,j,k)  = gr%centroid(3) + 0.5_sp
        EndIf
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
  Call Visual3DContour(f1=cx)
  Call Visual3DContour(f1=cy)
  Call Visual3DContour(f1=cz)

  Call MPI_FINALIZE(ierr)

End Subroutine InitVolumeCentroid
