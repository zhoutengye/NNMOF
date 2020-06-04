Module ModTools
  Use ModGlobal
  Implicit None
Contains
  Subroutine Visual3DContour(f1, f2, f3, f4, f5)
    Implicit None
    Real(sp), Intent(In) :: f1(0:nl(1)+1,0:nl(2)+1,0:nl(3)+1)
    Real(sp), Intent(In), optional :: f2(0:nl(1)+1,0:nl(2)+1,0:nl(3)+1)
    Real(sp), Intent(In), optional :: f3(0:nl(1)+1,0:nl(2)+1,0:nl(3)+1)
    Real(sp), Intent(In), optional :: f4(0:nl(1)+1,0:nl(2)+1,0:nl(3)+1)
    Real(sp), Intent(In), optional :: f5(0:nl(1)+1,0:nl(2)+1,0:nl(3)+1)
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

    If ( myid .eq. 0 ) Then
      open(10,file='vis.py',status='unknown')
      Write(10,'(a)') "import numpy as np"
      Write(10,'(a)') "import h5py"
      Write(10,'(a)') "from mayavi import mlab"
      Write(10,'(a)') "f = h5py.File('visual.h5')"
      Write(10,'(a)') "for key in f['visual']:"
      Write(10,'(a)') "    vis = np.array(f['visual'][key])"
      Write(10,'(a)') "    mlab.contour3d(vis,contours=8,opacity=.2 )"
      Write(10,'(a)') "nx, ny, nz = vis.shape"
      Write(10,'(a)') "f.close()"
      Write(10,'(a)') "mlab.outline(extent=[0,nx,0,ny,0,nz])"
      Write(10,'(a)') "mlab.show()"
      close(10)
      Call system('python vis.py')
      Call system('rm vis.py')
    End If

  End Subroutine Visual3DContour

End Module ModTools
