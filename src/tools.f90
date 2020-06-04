Module ModTools
  Use ModGlobal
  Implicit None
Contains
  Subroutine Visual3D(f)
    Implicit None
    Real(sp), Intent(In) :: f(0:nl(1)+1,0:nl(2)+1,0:nl(3)+1)
    Type(HDF5File)  :: h5_visual_file
    Type(HDF5Group) :: h5_visual_group
    Character(80) :: data_name

    INTEGER(HID_T) :: plist_id      ! Property list identifier 
    logical :: inputfield
    Integer :: h5error

    !    Remove the h5 file if existts
    if (myid .eq. 0) Call system('rm -f visual.h5')
    !    Create input file
    Call h5open_f(h5error)
    Call h5pcreate_f(H5P_FILE_ACCESS_F, plist_id, h5error)
    Call h5pset_fapl_mpio_f(plist_id, MPI_COMM_WORLD, MPI_INFO_NULL, h5error)
    ! Create the file collectively.
    h5_visual_file%filename = "visial.h5"
    Call h5fcreate_f(h5_visual_file%filename, H5F_ACC_TRUNC_F, h5_visual_file%file_id, h5error, access_prp = plist_id)
    Call h5pclose_f(plist_id, h5error)
    ! Create group
    h5_visual_group%groupname = "group.h5"
    Call HDF5CreateGroup(h5_visual_file, h5_visual_group)
    ! Write dataset
    data_name = 'visual'
    ! Call HDF5WriteData(h5_visual_group, f, data_name)
  End Subroutine Visual3D
End Module ModTools
