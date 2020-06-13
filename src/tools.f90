# include "param.h"
Module ModTools
  Use ModGlobal
  Implicit None

  Type Octree
    Real(sp) :: xc(3)
    Real(sp) :: dx(3)
    Real(sp) :: centroids(2,2,2,3)
    Real(sp) :: vofs(2,2,2)
    Real(sp) :: centroid(3)
    Real(sp) :: vof
    Integer  :: level
    Integer  :: maxlevel
  End Type Octree
  PROCEDURE(InterfaceLS), POINTER :: ShapeLevelSet
  Interface
    Real(sp) Function InterfaceLS(x,y,z)
      import
      Real(sp) :: x, y, z
    End Function  InterfaceLS
  End Interface

Contains
  Subroutine Visual3DContour(f1, f2, f3, f4, f5, slice_dir, slice_coord)
    Implicit None
    Real(sp), Intent(In) :: f1(0:nl(1)+1,0:nl(2)+1,0:nl(3)+1)
    Real(sp), Intent(In), optional :: f2(0:nl(1)+1,0:nl(2)+1,0:nl(3)+1)
    Real(sp), Intent(In), optional :: f3(0:nl(1)+1,0:nl(2)+1,0:nl(3)+1)
    Real(sp), Intent(In), optional :: f4(0:nl(1)+1,0:nl(2)+1,0:nl(3)+1)
    Real(sp), Intent(In), optional :: f5(0:nl(1)+1,0:nl(2)+1,0:nl(3)+1)
    Character,  Intent(In), optional :: slice_dir
    Integer,  Intent(In), optional :: slice_coord
    Character(5) :: slice_num
    Type(HDF5File)  :: h5_visual_file
    Type(HDF5Group) :: h5_visual_group
    Character(80) :: data_name

    INTEGER(HID_T) :: plist_id      ! Property list identifier 
    Integer :: h5error

    !    Remove the h5 file if existts
    if (myid .eq. 0) Call system('rm -f visual.h5')
    !    Create input file
    Call h5open_f(h5error)
    Call h5pcreate_f(H5P_FILE_ACCESS_F, plist_id, h5error)
    Call h5pset_fapl_mpio_f(plist_id, MPI_COMM_WORLD, MPI_INFO_NULL, h5error)
    ! Create the file collectively.
    h5_visual_file%filename = "visual.h5"
    Call h5fcreate_f(h5_visual_file%filename, H5F_ACC_TRUNC_F, h5_visual_file%file_id, h5error, access_prp = plist_id)
    Call h5pclose_f(plist_id, h5error)
    ! Create group
    h5_visual_group%groupname = "visual"
    Call HDF5CreateGroup(h5_visual_file, h5_visual_group)
    ! Write dataset
    data_name = 'vis01'
    Call HDF5WriteData(h5_visual_group, f1, data_name)
    if(present(f2))then
      data_name = 'vis02'
      Call HDF5WriteData(h5_visual_group, f2, data_name)
    end if
    if(present(f3))then
      data_name = 'vis03'
      Call HDF5WriteData(h5_visual_group, f3, data_name)
    end if
    if(present(f4))then
      data_name = 'vis04'
      Call HDF5WriteData(h5_visual_group, f4, data_name)
    end if
    if(present(f5))then
      data_name = 'vis05'
      Call HDF5WriteData(h5_visual_group, f5, data_name)
    End If

    Call h5gclose_f(h5_visual_group%group_id, h5error)
    Call h5fclose_f(h5_visual_file%file_id, h5error)

    write(slice_num , '(i5)') slice_coord
    print *, slice_num

    If ( myid .eq. 0 ) Then
      open(10,file='vis3dcontour.py',status='unknown')
      Write(10,'(a)') "import numpy as np"
      Write(10,'(a)') "import h5py"
      Write(10,'(a)') "from mayavi import mlab"
      Write(10,'(a)') "f = h5py.File('visual.h5','r')"
      Write(10,'(a)') "for key in f['visual']:"
      Write(10,'(a)') "    vis = np.array(f['visual'][key])"
      Write(10,'(a)') "    mlab.contour3d(vis,contours=8,opacity=.2 )"
      Write(10,'(a)') "nx, ny, nz = vis.shape"
      Write(10,'(a)') "f.close()"
      if(present(slice_dir))then
        Write(10,'(a)') "mlab.pipeline.image_plane_widget(mlab.pipeline.scalar_field(vis),plane_orientation='"&
            //slice_dir//"_axes',slice_index="//slice_num//")"
      endif
      Write(10,'(a)') "mlab.outline(extent=[0,nx,0,ny,0,nz])"
      Write(10,'(a)') "mlab.show()"
      close(10)
      Call system('python vis3dcontour.py')
      ! Call system('rm vis3dcontour.py')
    End If

  End Subroutine Visual3DContour

 Subroutine Visual2DContour(f1, f2, f3, f4, f5, slice_dir, slice_coord)
    Implicit None
    Real(sp), Intent(In) :: f1(0:nl(1)+1,0:nl(2)+1,0:nl(3)+1)
    Real(sp), Intent(In), optional :: f2(0:nl(1)+1,0:nl(2)+1,0:nl(3)+1)
    Real(sp), Intent(In), optional :: f3(0:nl(1)+1,0:nl(2)+1,0:nl(3)+1)
    Real(sp), Intent(In), optional :: f4(0:nl(1)+1,0:nl(2)+1,0:nl(3)+1)
    Real(sp), Intent(In), optional :: f5(0:nl(1)+1,0:nl(2)+1,0:nl(3)+1)
    Integer, optional :: slice_dir
    Integer, optional :: slice_coord
    Character(5) :: slice_num
    Type(HDF5File)  :: h5_visual_file
    Type(HDF5Group) :: h5_visual_group
    Character(80) :: data_name

    INTEGER(HID_T) :: plist_id      ! Property list identifier 
    Integer :: h5error

    !    Remove the h5 file if existts
    if (myid .eq. 0) Call system('rm -f visual.h5')
    !    Create input file
    Call h5open_f(h5error)
    Call h5pcreate_f(H5P_FILE_ACCESS_F, plist_id, h5error)
    Call h5pset_fapl_mpio_f(plist_id, MPI_COMM_WORLD, MPI_INFO_NULL, h5error)
    ! Create the file collectively.
    h5_visual_file%filename = "visual.h5"
    Call h5fcreate_f(h5_visual_file%filename, H5F_ACC_TRUNC_F, h5_visual_file%file_id, h5error, access_prp = plist_id)
    Call h5pclose_f(plist_id, h5error)
    ! Create group
    h5_visual_group%groupname = "visual"
    Call HDF5CreateGroup(h5_visual_file, h5_visual_group)
    ! Write dataset
    data_name = 'vis01'
    Call HDF5WriteData(h5_visual_group, f1, data_name)
    if(present(f2))then
      data_name = 'vis02'
      Call HDF5WriteData(h5_visual_group, f2, data_name)
    end if
    if(present(f3))then
      data_name = 'vis03'
      Call HDF5WriteData(h5_visual_group, f3, data_name)
    end if
    if(present(f4))then
      data_name = 'vis04'
      Call HDF5WriteData(h5_visual_group, f4, data_name)
    end if
    if(present(f5))then
      data_name = 'vis05'
      Call HDF5WriteData(h5_visual_group, f5, data_name)
    End If

    Call h5gclose_f(h5_visual_group%group_id, h5error)
    Call h5fclose_f(h5_visual_file%file_id, h5error)

    write(slice_num , '(i5)') slice_coord-1
    print *, slice_num

    If ( myid .eq. 0 ) Then
      open(10,file='vis2dcontour.py',status='unknown')
      Write(10,'(a)') "import numpy as np"
      Write(10,'(a)') "import h5py"
      Write(10,'(a)') "import matplotlib.pyplot as plt"
      Write(10,'(a)') "f = h5py.File('visual.h5','r')"
      Write(10,'(a)') "for key in f['visual']:"
      If ( slice_dir .eq. 3) Then
        Write(10,'(a)') "    vis = np.array(f['visual'][key]["//slice_num//",:,:])"
      Else If ( slice_dir .eq. 2) Then
        Write(10,'(a)') "    vis = np.array(f['visual'][key][:,"//slice_num//",:])"
      Else If ( slice_dir .eq. 1) Then
        Write(10,'(a)') "    vis = np.array(f['visual'][key][:,:,"//slice_num//"])"
      End If
      Write(10,'(a)') "    plt.contourf(vis)"
      Write(10,'(a)') "f.close()"
      ! Write(10,'(a)') "plt.colorbar()"
      Write(10,'(a)') "plt.show()"
      close(10)
      Call system('python vis2dcontour.py')
      ! Call system('rm vis2dcontour.py')
    End If

  End Subroutine Visual2DContour

  Recursive Subroutine VolumeCentroidOctree(Tree)
    Implicit None
    Type(Octree) :: Tree
    Type(Octree) :: TreeChild
    Real(sp) :: xc, yc, zc
    Real(sp) :: halfdx, halfdy, halfdz
    Integer  :: flag
    Integer  :: level
    Integer  :: dir
    Integer  :: i, j, k
    Real(sp) :: moment(3)

    level = Tree%level + 1
    halfdx = Tree%dx(1) / 2.0_sp
    halfdy = Tree%dx(2) / 2.0_sp
    halfdz = Tree%dx(3) / 2.0_sp

    Do k = 1, 2
      Do j = 1, 2
        Do i = 1, 2
          ! actually pass xc +|- 1/4dx
          xc = Tree%xc(1) + halfdx * (dble(i) - 1.5_sp)
          yc = Tree%xc(2) + halfdy * (dble(j) - 1.5_sp)
          zc = Tree%xc(3) + halfdz * (dble(k) - 1.5_sp)
          flag = InOut(xc, yc, zc, halfdx, halfdy, halfdz)
          ! Dark cell
          If ( flag .eq. 1) Then
            Tree%vofs(i,j,k) = 1.0_sp
            Tree%centroids(i,j,k,1)  = xc
            Tree%centroids(i,j,k,2)  = yc
            Tree%centroids(i,j,k,3)  = zc
          ! light cell
          ElseIf ( flag .eq. 2) Then
            Tree%vofs(i,j,k) = 0.0_sp
            Tree%centroids(i,j,k,1)  = 0.0_sp
            Tree%centroids(i,j,k,2)  = 0.0_sp
            Tree%centroids(i,j,k,3)  = 0.0_sp
          ! interface cell
          Else If ( flag .eq. 3) Then
            ! stop at maximum level
            If (level >= tree%maxlevel) Then
              Tree%vofs(i,j,k) = 1.0_sp
              Tree%centroids(i,j,k,1)  = 0.5_sp
              Tree%centroids(i,j,k,2)  = 0.5_sp
              Tree%centroids(i,j,k,3)  = 0.5_sp
            ! Recursive octree
            Else
              TreeChild%level = level
              TreeChild%maxlevel = Tree%maxlevel
              TreeChild%xc(1) = xc
              TreeChild%xc(2) = yc
              TreeChild%xc(3) = zc
              TreeChild%dx(1) = halfdx
              TreeChild%dx(2) = halfdy
              TreeChild%dx(3) = halfdz
              Call VolumeCentroidOctree(TreeChild)
              Tree%vofs(i,j,k)  = TreeChild%vof
              Tree%centroids(i,j,k,1)  = TreeChild%centroid(1)
              Tree%centroids(i,j,k,2)  = TreeChild%centroid(2)
              Tree%centroids(i,j,k,3)  = TreeChild%centroid(3)
            End If
          EndIf
        End Do
      End Do
    End Do

    Tree%vof = 0.0_sp
    moment = 0
    Do k = 1, 2
      Do j = 1, 2
        Do i = 1, 2
          Tree%vof = Tree%vof + Tree%vofs(i,j,k)
          Do dir = 1,3
            moment(dir) = moment(dir) + Tree%vofs(i,j,k) * Tree%Centroids(i,j,k,dir)
          End Do
        End Do
      End Do
    End Do

    Tree%vof = Tree%vof / 8.0_sp

    Do dir = 1,3
      Tree%centroid(dir) = moment(dir) / (Tree%vof + 1e-12) / 8.0_sp
    End Do

    If (Tree%vof <= 1e-12) Then
      Tree%centroid(1:3) = 0.0_sp
      Tree%vof = 0.0_sp
    ElseIf (Tree%vof >= 1.0_sp - 1e-12) Then
      Tree%centroid(1:3) = 0.5_sp
      Tree%vof = 1.0_sp
    End If

  End Subroutine VolumeCentroidOctree

 Recursive Subroutine VolumeOctree(Tree)
    Implicit None
    Type(Octree) :: Tree
    Type(Octree) :: TreeChild
    Real(sp) :: xc, yc, zc
    Real(sp) :: halfdx, halfdy, halfdz
    Integer  :: flag
    Integer  :: level
    Integer  :: i, j, k

    level = Tree%level + 1
    halfdx = Tree%dx(1) / 2.0_sp
    halfdy = Tree%dx(2) / 2.0_sp
    halfdz = Tree%dx(3) / 2.0_sp

    ! print *, level, Tree%xc, Tree%dx

    Do k = 1, 2
      Do j = 1, 2
        Do i = 1, 2
          ! actually pass xc +|- 1/4dx
          xc = Tree%xc(1) + halfdx * (dble(i) - 1.5_sp)
          yc = Tree%xc(2) + halfdy * (dble(j) - 1.5_sp)
          zc = Tree%xc(3) + halfdz * (dble(k) - 1.5_sp)
          flag = InOut(xc, yc, zc, halfdx, halfdy, halfdz)
          ! Dark cell
          If ( flag .eq. 1) Then
            Tree%vofs(i,j,k) = 1.0_sp
          ! light cell
            ! print  *, '111', xc, yc, zc
          ElseIf ( flag .eq. 2) Then
            Tree%vofs(i,j,k) = 0.0_sp
            ! print  *, '000', xc, yc, zc , sqrt(xc**2+yc**2+zc**2)
          ! interface cell
          Else If ( flag .eq. 3) Then
            ! stop at maximum level
            If (level >= tree%maxlevel) Then
              ! Tree%vofs(i,j,k) = 1.0_sp
              Tree%vofs(i,j,k) = heaviside(ShapeLevelSet(xc,yc,zc),halfdx)
            ! Recursive octree
            Else
              TreeChild%level = level
              TreeChild%maxlevel = Tree%maxlevel
              TreeChild%xc(1) = xc
              TreeChild%xc(2) = yc
              TreeChild%xc(3) = zc
              TreeChild%dx(1) = halfdx
              TreeChild%dx(2) = halfdy
              TreeChild%dx(3) = halfdz
              Call VolumeOctree(TreeChild)
              Tree%vofs(i,j,k)  = TreeChild%vof
            End If
          EndIf
        End Do
      End Do
    End Do

    Tree%vof = 0.0_sp
    Do k = 1, 2
      Do j = 1, 2
        Do i = 1, 2
          Tree%vof = Tree%vof + Tree%vofs(i,j,k)
        End Do
      End Do
    End Do

    Tree%vof = Tree%vof / 8.0_sp

    If (Tree%vof <= 1e-12) Then
      Tree%vof = 0.0_sp
    ElseIf (Tree%vof >= 1.0_sp - 1e-12) Then
      Tree%vof = 1.0_sp
    End If


  End Subroutine VolumeOctree

  Integer Function InOut(xc, yc, zc, dx, dy, dz)
    Implicit None
    Real(sp) :: xc, yc, zc
    Real(sp) :: dx, dy, dz

    Real(sp) :: xx, yy, zz
    Real(sp) :: lsoctree(2,2,2)
    Logical  :: flag(2,2,2)
    Integer  :: i, j, k

    Do k = 1, 2
      Do j = 1, 2
        Do i = 1, 2
          xx = xc + dx * (dble(i) - 1.5_sp)
          yy = yc + dy * (dble(j) - 1.5_sp)
          zz = zc + dz * (dble(k) - 1.5_sp)
          lsoctree(i,j,k) = ShapeLevelSet(xx, yy, zz)
          If ( lsoctree(i,j,k) .gt. 1.0e-10) Then
            flag(i,j,k) = .true.
          Else
            flag(i,j,k) = .false.
          End If
        End Do
      End Do
    End Do

    If ( all(flag) ) Then
      Inout = 1
    ElseIf ( any(flag) ) Then
      Inout = 3
    Else
      Inout = 2
    EndIf

 End Function InOut

 Real(sp) Function ShapeLS(x,y,z)
   Implicit None
   Real(sp) :: x, y, z
   Print *, "x=", x, 'y=', y, 'z=', z
   Print *, "Please specify the shape function"
   ShapeLs = 8888.888888888_sp
 End Function ShapeLS

 Real(sp) Function Heaviside(dis,h)
   Implicit None
   Real(sp) :: dis
   Real(sp) :: h
   If ( dis .ge. h) Then
     heaviside =  1.0_sp
   Else If ( dis .le. -h) Then
     heaviside = 0.0_sp
   Else
     heaviside = 0.5_sp * ( 1.0_sp + dis / h + 1.0_sp / Pi * sin(Pi * dis / h) )
   End If
 End Function Heaviside

End Module ModTools
