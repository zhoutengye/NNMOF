# define epsc 1.0e-12
program RiderKotheInit
  use ModGlobal
  Implicit None
  Integer :: ng

  ! numer of grids
  ng = 100

  call InitVolume(ng)

end program RiderKotheInit


Subroutine InitVolume(ng)
  Use ModTools
  Use InitSH
  Implicit None
  Integer :: ng
  Real(sp), Allocatable :: xc(:)
  Type(Octree) :: gr
  Integer :: i,j,k
  Integer :: flag
  Real(sp) :: dx
  Integer :: Octreelevel = 3
  Real(sp) :: xx, yy, zz

  Call Init(inputfield=.false.)


  Allocate(xc(ng))
  Do k = 1, 100
    xc(k) = dble(k)/ dble(ng) - 0.5_sp
  End Do
  dx = 1.0_sp / dble(ng)
  ShapeLevelSet => ShapeRiderKothe
  Do k = 1, ng
    Do j = 1, ng
      Do i = 1, ng
        xx = xc(i) - 0.28_sp
        yy = xc(j) - 0.28_sp
        zz = xc(k) - 0.28_sp
        flag = Inout(xx, yy, zz, dx, dx, dx)
        If (flag .eq. 1) Then
          phi(i,j,k) = 1.0
          cx(i,j,k)  = 0.5_sp
          cy(i,j,k)  = 0.5_sp
          cz(i,j,k)  = 0.5_sp
        Elseif (flag .eq. 2) Then
          phi(i,j,k) = 0.0
          cx(i,j,k)  = 0.0_sp
          cy(i,j,k)  = 0.0_sp
          cz(i,j,k)  = 0.0_sp
        Elseif (flag .eq. 3) Then
          gr%level = 1
          gr%maxlevel = Octreelevel
          gr%xc(1) = xx
          gr%xc(2) = yy
          gr%xc(3) = zz
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
