  end subroutine read_input

#if defined (HDF5IO)
Subroutine H5Init(h5file)
  use hdf5
  use global
  Implicit None
  INTEGER(HID_T) :: plist_id      ! Property list identifier 
  Type(hdf5file) :: h5file

  Integer :: h5error

  h5_total_dims(1) = Mglob
  h5_total_dims(2) = Nglob
  h5_block_dims(1) = Mloc - 2 * Nghost
  h5_block_dims(2) = Nloc - 2 * Nghost


  h5_offset(1) = npx*int(mglob/px) + min(npx, mod( mglob, px ))
  h5_offset(2) = npy*int(nglob/py) + min(npy, mod( nglob, py ))


  if (myid .eq. 0 ) Call system('rm -f'//trim(h5file%filename))

  Call MPI_barrier(MPI_COMM_WORLD, h5error)

  Call h5open_f(h5error)
  CALL h5pcreate_f(H5P_FILE_ACCESS_F, plist_id, h5error)
  CALL h5pset_fapl_mpio_f(plist_id, MPI_COMM_WORLD, MPI_INFO_NULL, h5error)

  ! Create the file collectively.
  CALL h5fcreate_f(h5file%filename, H5F_ACC_TRUNC_F, h5file%file_id, h5error, access_prp = plist_id)
  CALL h5pclose_f(plist_id, h5error)


  IF(OUT_DEPTH.OR.BREAKWATER)THEN
    h5group_dep%groupname   = 'dep'
    Call HDF5CreateGroup(h5file, h5group_dep)
    h5group_breakwater%groupname   = 'breakwater'
    Call HDF5CreateGroup(h5file, h5group_breakwater)
  ENDIF

  IF(OUT_ETA)THEN
    h5group_eta%groupname = 'eta'
    Call HDF5CreateGroup(h5file, h5group_eta)
  end IF

  IF(OUT_Hmax)THEN
    h5group_hmax%groupname = 'Hmax'
    Call HDF5CreateGroup(h5file, h5group_hmax)
  ENDIF

  IF(OUT_Hmin)THEN
    h5group_hmin%groupname = 'Hmin'
    Call HDF5CreateGroup(h5file, h5group_hmin)
  ENDIF

  IF(OUT_Umax)THEN
    h5group_Umax%groupname = 'Umax'
    call HDF5CreateGroup(h5file, h5group_umax)
  ENDIF

  IF(OUT_MFmax)THEN
    h5group_MFmax%groupname = 'MFmax'
    call HDF5CreateGroup(h5file, h5group_MFmax)
  ENDIF

  IF(OUT_VORmax)THEN
    h5group_VORmax%groupname = 'VORmax'
    call HDF5CreateGroup(h5file, h5group_VORmax)
  ENDIF

  IF(OUT_MASK)THEN
    h5group_mask%groupname = 'mask'
    call HDF5CreateGroup(h5file, h5group_mask)
  ENDIF

  IF(OUT_MASK9)THEN
    h5group_mask9%groupname = 'mask9'
    call HDF5CreateGroup(h5file, h5group_mask9)
  ENDIF

  IF(OUT_U)THEN
    h5group_u%groupname   = 'U'
    Call HDF5CreateGroup(h5file, h5group_u)
  EndIf
  IF(OUT_V)THEN
    h5group_v%groupname   = 'V'
    Call HDF5CreateGroup(h5file, h5group_v)
  end IF
  IF(OUT_P)THEN
    h5group_P%groupname   = 'P'
    Call HDF5CreateGroup(h5file, h5group_P)
  end IF
  IF(OUT_Q)THEN
    h5group_Q%groupname   = 'Q'
    Call HDF5CreateGroup(h5file, h5group_Q)
  end IF

  
  IF(OUT_AGE)THEN
    IF(SHOW_BREAKING)THEN
      h5group_age%groupname   = 'age'
      call HDF5CreateGroup(h5file, h5group_age)
    ENDIF
  ENDIF

  IF(OUT_ROLLER)THEN
    h5group_roller%groupname   = 'roller'
    call HDF5CreateGroup(h5file, h5group_roller)
  ENDIF

  IF(OUT_UNDERTOW)THEN
    h5group_undertow_u%groupname   = 'undertow_u'
    h5group_undertow_v%groupname   = 'undertow_v'
    call HDF5CreateGroup(h5file, h5group_undertow_u)
    call HDF5CreateGroup(h5file, h5group_undertow_v)
  ENDIF

  IF(VISCOSITY_BREAKING)THEN
    IF(OUT_NU)THEN
      h5group_nubrk%groupname   = 'nu_break'
      call HDF5CreateGroup(h5file, h5group_nubrk)
    ENDIF
  ENDIF

# if defined (METEO)
   IF(OUT_METEO)THEN

    IF(WindHollandModel)THEN
      h5group_Pstorm%groupname   = 'P_storm'
      h5group_Ustorm%groupname   = 'U_storm'
      h5group_Vstorm%groupname   = 'V_storm'
      call HDF5CreateGroup(h5file, h5group_Pstorm)
      call HDF5CreateGroup(h5file, h5group_Ustorm)
      call HDF5CreateGroup(h5file, h5group_Vstorm)
    ENDIF

    IF(MeteoGausian)THEN
      h5group_Pstorm%groupname   = 'P_storm'
      call HDF5CreateGroup(h5file, h5group_Pstorm)
    ENDIF

    IF(SlideModel)THEN
      h5group_Pstorm%groupname   = 'P_storm'
      call HDF5CreateGroup(h5file, h5group_Pstorm)
    ENDIF

  ENDIF
# endif
     
# if defined (VESSEL)
     IF(OUT_VESSEL)THEN
# if defined (VESSEL_PANEL_SOURCE)
       h5group_Fves%groupname   = 'Fves'
       call HDF5CreateGroup(h5file, h5group_Fves)
# else
       h5group_Pves%groupname   = 'Pves'
       call HDF5CreateGroup(h5file, h5group_Pves)
# endif
   ! end panel source

     ENDIF
# endif

! sediment
# if defined (SEDIMENT)
     h5group_c%groupname   = 'c'
     h5group_pick%groupname   = 'pick'
     h5group_dego%groupname   = 'depo'
     h5group_Dchgs%groupname   = 'Dchgs'
     h5group_DchgB%groupname   = 'DchgB'
     h5group_BedFx%groupname   = 'Bed_Fx'
     h5group_BedFy%groupname   = 'Bed_Fy'
     call HDF5CreateGroup(h5file, h5group_c)
     call HDF5CreateGroup(h5file, h5group_pick)
     call HDF5CreateGroup(h5file, h5group_depo)
     call HDF5CreateGroup(h5file, h5group_Dchgs)
     call HDF5CreateGroup(h5file, h5group_DchgB)
     call HDF5CreateGroup(h5file, h5group_BedFx)
     call HDF5CreateGroup(h5file, h5group_BedFy)
# endif


End Subroutine H5Init

Subroutine H5Finalize(h5file)
  use hdf5
  use global
  Implicit None
  Integer :: error
  Type(hdf5file) :: h5file

  Call HDF5_WRITE_PARAMETERS(h5file)

  IF(OUT_DEPTH.OR.BREAKWATER)THEN
    h5group_dep%groupname   = 'dep'
    Call HDF5CloseGroup( h5group_dep)
    h5group_breakwater%groupname   = 'breakwater'
    Call HDF5CloseGroup( h5group_breakwater)
  ENDIF

  IF(OUT_ETA)THEN
    h5group_eta%groupname = 'eta'
    Call HDF5CloseGroup( h5group_eta)
  end IF

  IF(OUT_Hmax)THEN
    h5group_hmax%groupname = 'Hmax'
    Call HDF5CloseGroup( h5group_hmax)
  ENDIF

  IF(OUT_Hmin)THEN
    h5group_hmin%groupname = 'Hmin'
    Call HDF5CloseGroup( h5group_hmin)
  ENDIF

  IF(OUT_Umax)THEN
    call HDF5CloseGroup( h5group_umax)
  ENDIF

  IF(OUT_MFmax)THEN
    call HDF5CloseGroup( h5group_MFmax)
  ENDIF

  IF(OUT_VORmax)THEN
    call HDF5CloseGroup( h5group_VORmax)
  ENDIF

  IF(OUT_MASK)THEN
    call HDF5CloseGroup( h5group_mask)
  ENDIF

  IF(OUT_MASK9)THEN
    call HDF5CloseGroup( h5group_mask9)
  ENDIF

  IF(OUT_U)THEN
    h5group_u%groupname   = 'U'
    Call HDF5CloseGroup( h5group_u)
  EndIf
  IF(OUT_V)THEN
    h5group_v%groupname   = 'V'
    Call HDF5CloseGroup( h5group_v)
  end IF
  IF(OUT_P)THEN
    h5group_P%groupname   = 'P'
    Call HDF5CloseGroup( h5group_P)
  end IF
  IF(OUT_Q)THEN
    h5group_Q%groupname   = 'Q'
    Call HDF5CloseGroup( h5group_Q)
  end IF

  
  IF(OUT_AGE)THEN
    IF(SHOW_BREAKING)THEN
      h5group_age%groupname   = 'age'
      call HDF5CloseGroup( h5group_age)
    ENDIF
  ENDIF

  IF(OUT_ROLLER)THEN
    h5group_roller%groupname   = 'roller'
    call HDF5CloseGroup( h5group_roller)
  ENDIF

  IF(OUT_UNDERTOW)THEN
    h5group_undertow_u%groupname   = 'undertow_u'
    h5group_undertow_v%groupname   = 'undertow_v'
    call HDF5CloseGroup( h5group_undertow_u)
    call HDF5CloseGroup( h5group_undertow_v)
  ENDIF

  IF(VISCOSITY_BREAKING)THEN
    IF(OUT_NU)THEN
      h5group_nubrk%groupname   = 'nu_break'
      call HDF5CloseGroup( h5group_nubrk)
    ENDIF
  ENDIF

# if defined (METEO)
   IF(OUT_METEO)THEN

    IF(WindHollandModel)THEN
      h5group_Pstorm%groupname   = 'P_storm'
      h5group_Ustorm%groupname   = 'U_storm'
      h5group_Vstorm%groupname   = 'V_storm'
      call HDF5CloseGroup( h5group_Pstorm)
      call HDF5CloseGroup( h5group_Ustorm)
      call HDF5CloseGroup( h5group_Vstorm)
    ENDIF

    IF(MeteoGausian)THEN
      h5group_Pstorm%groupname   = 'P_storm'
      call HDF5CloseGroup( h5group_Pstorm)
    ENDIF

    IF(SlideModel)THEN
      h5group_Pstorm%groupname   = 'P_storm'
      call HDF5CloseGroup( h5group_Pstorm)
    ENDIF

  ENDIF
# endif
     
# if defined (VESSEL)
     IF(OUT_VESSEL)THEN
# if defined (VESSEL_PANEL_SOURCE)
       h5group_Fves%groupname   = 'Fves'
       call HDF5CloseGroup( h5group_Fves)
# else
       h5group_Pves%groupname   = 'Pves'
       call HDF5CloseGroup( h5group_Pves)
# endif
   ! end panel source

     ENDIF
# endif

! sediment
# if defined (SEDIMENT)
     h5group_c%groupname   = 'c'
     h5group_pick%groupname   = 'pick'
     h5group_dego%groupname   = 'depo'
     h5group_Dchgs%groupname   = 'Dchgs'
     h5group_DchgB%groupname   = 'DchgB'
     h5group_BedFx%groupname   = 'Bed_Fx'
     h5group_BedFy%groupname   = 'Bed_Fy'
     call HDF5CloseGroup( h5group_c)
     call HDF5CloseGroup( h5group_pick)
     call HDF5CloseGroup( h5group_depo)
     call HDF5CloseGroup( h5group_Dchgs)
     call HDF5CloseGroup( h5group_DchgB)
     call HDF5CloseGroup( h5group_BedFx)
     call HDF5CloseGroup( h5group_BedFy)
# endif

     

  Call h5fclose_f(h5file%file_id, error)

#if !defined (FASTHDF5IO)
  call Easy_HDF5
#endif

End Subroutine H5Finalize

Subroutine HDF5CreateGroup(h5file, h5group)
  use global
  Implicit None

  Type(hdf5file) :: h5file
  Type(hdf5group) :: h5group

  Integer :: h5error

  CALL h5gcreate_f(h5file%file_id, h5group%groupname, h5group%group_id, h5error)

End Subroutine HDF5CreateGroup

Subroutine HDF5CloseGroup(h5group)
  use global
  Implicit None

  Type(hdf5group) :: h5group

  Integer :: h5error

  CALL h5gclose_f(h5group%group_id, h5error)

End Subroutine HDF5CloseGroup

SUBROUTINE HDF5PutFile(h5group, data, data_name)
  USE GLOBAL
  IMPLICIT NONE

  Type(hdf5group) :: h5group
  Real(SP), Intent(In) :: data(Mloc,Nloc)
  Character(80) :: data_name
  Character(80) :: data_name1

  INTEGER(HID_T) :: filespace     ! Dataspace identifier in file 
  INTEGER(HID_T) :: memspace      ! Dataspace identifier in memory
  INTEGER(HID_T) :: plist_id      ! Property list identifier 
  INTEGER(HID_T) :: dset_id       ! Data set identifier 

  Integer :: mblock, nblock

  integer :: rank = 2

  Integer :: h5error

  data_name1 = trim(h5group%groupname) // '_' // trim(data_name)

  ! Create the data space for the  dataset. 
  CALL h5screate_simple_f(rank, h5_total_dims, filespace, h5error)

  ! Create the dataset with default properties.
# if defined DOUBLE_PRECISION
  CALL h5dcreate_f(h5group%group_id, data_name1, H5T_NATIVE_DOUBLE, filespace, &
      dset_id, h5error)
#else
  CALL h5dcreate_f(h5group%group_id, data_name1, H5T_NATIVE_REAL, filespace, &
      dset_id, h5error)
#endif
  CALL h5sclose_f(filespace, h5error)

  ! Create the data space for the  dataset. 
  CALL h5screate_simple_f(rank, h5_block_dims, memspace, h5error)


  ! Select hyperslab in the file.
    !
    CALL h5dget_space_f(dset_id, filespace, h5error)
    CALL h5sselect_hyperslab_f (filespace, H5S_SELECT_SET_F, h5_offset, h5_block_dims, h5error)

    CALL h5pcreate_f(H5P_DATASET_XFER_F, plist_id, h5error) 
    CALL h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, h5error)


# if defined DOUBLE_PRECISION
    CALL h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, data(ibeg:iend,jbeg:jend), h5_total_dims ,error, &
        file_space_id = filespace, mem_space_id = memspace, xfer_prp=plist_id)
# else
    CALL h5dwrite_f(dset_id, H5T_NATIVE_REAL, data(ibeg:iend,jbeg:jend), h5_total_dims ,h5error, &
          file_space_id = filespace, mem_space_id = memspace, xfer_prp=plist_id)
# endif


    ! Close dataspaces.
    CALL h5sclose_f(filespace, h5error)
    CALL h5sclose_f(memspace, h5error)
    
    ! Close the dataset.
    CALL h5dclose_f(dset_id, h5error)

    Call h5pclose_f(plist_id, h5error)


end SUBROUTINE HDF5PutFile

SUBROUTINE HDF5_WRITE_PARAMETERS(h5file)
  USE GLOBAL
  Implicit None
  Type(hdf5file) :: h5file
  Character(80) :: tmp_name 


  ! write the total grid information
  tmp_name = 'nx'
  Call add_Parameter_integer(h5file%file_id, tmp_name, mglob, 1)
  tmp_name = 'ny'
  Call add_Parameter_integer(h5file%file_id, tmp_name, nglob, 1)
  ! write the grid resolution
  tmp_name = 'dx'
  Call add_Parameter_float(h5file%file_id, tmp_name, dx, 1)
  tmp_name = 'dy'
  Call add_Parameter_float(h5file%file_id, tmp_name, dy, 1)
  ! write the count of data
  tmp_name = 'out_interval'
  Call add_Parameter_float(h5file%file_id, tmp_name, PLOT_INTV,1)
  tmp_name = 'out_count'
  Call add_Parameter_integer(h5file%file_id, tmp_name, icount, 1)

END SUBROUTINE HDF5_WRITE_PARAMETERS

Subroutine add_Parameter_character(id, p_name, p_value)
  use HDF5
  Implicit None
  Integer(HID_T) :: id
  Character(80) :: p_name
  Character(80) :: p_value
  INTEGER(HID_T) :: aspace_id     ! Attribute Dataspace identifier
  INTEGER(HID_T) :: attr_id       ! Attribute Dataspace identifier
  INTEGER(HID_T) :: atype_id      ! Attribute Dataspace identifier
  INTEGER(HSIZE_T), DIMENSION(1) :: adims = (/1/) ! Attribute dimension
  Integer(Size_T) :: attrlen=80
  INTEGER(HSIZE_T), DIMENSION(1) :: data_dims = (/1/)
  INTEGER     ::   arank = 1                      ! Attribure rank

  Integer :: error

  CALL h5screate_simple_f(arank, adims, aspace_id, error)
  CALL h5tcopy_f(H5T_NATIVE_CHARACTER, atype_id, error)
  CALL h5tset_size_f(atype_id, attrlen, error)
  CALL h5acreate_f(id, p_name, atype_id, aspace_id, attr_id, error)
  CALL h5awrite_f(attr_id, atype_id, p_value, data_dims, error)
  CALL h5aclose_f(attr_id, error)
  CALL h5tclose_f(atype_id, error)
  CALL h5sclose_f(aspace_id, error)

End Subroutine add_Parameter_character

Subroutine add_Parameter_Integer(id, p_name, p_value, dim)
  use HDF5
  Implicit None
  Integer(HID_T) :: id
  Integer :: dim
  Character(80) :: p_name
  integer :: p_value
  INTEGER(HID_T) :: aspace_id     ! Attribute Dataspace identifier
  INTEGER(HID_T) :: attr_id       ! Attribute Dataspace identifier
  INTEGER(HID_T) :: atype_id      ! Attribute Dataspace identifier
  INTEGER(HSIZE_T), DIMENSION(1) :: adims = (/1/) ! Attribute dimension
  INTEGER(HSIZE_T), DIMENSION(1) :: data_dims = (/1/)
  INTEGER     ::   arank = 1                      ! Attribure rank

  Integer :: error

  adims = dim
  data_dims = dim

  CALL h5screate_simple_f(arank, adims, aspace_id, error)
  CALL h5tcopy_f(H5T_NATIVE_INTEGER, atype_id, error)
  CALL h5acreate_f(id, p_name, atype_id, aspace_id, attr_id, error)
  CALL h5awrite_f(attr_id, atype_id, p_value, data_dims, error)
  CALL h5aclose_f(attr_id, error)
  CALL h5tclose_f(atype_id, error)
  CALL h5sclose_f(aspace_id, error)

End Subroutine add_Parameter_Integer

Subroutine add_Parameter_Float(id, p_name, p_value, dim)
  use param
  use HDF5
  Implicit None
  Integer(HID_T) :: id
  Integer :: dim
  Character(80) :: p_name
  Real(SP) :: p_value
  INTEGER(HID_T) :: aspace_id     ! Attribute Dataspace identifier
  INTEGER(HID_T) :: attr_id       ! Attribute Dataspace identifier
  INTEGER(HID_T) :: atype_id      ! Attribute Dataspace identifier
  INTEGER(HSIZE_T), DIMENSION(1) :: adims = (/1/) ! Attribute dimension
  INTEGER(HSIZE_T), DIMENSION(1) :: data_dims = (/1/)
  INTEGER     ::   arank = 1                      ! Attribure rank

  Integer :: error

  adims = dim
  data_dims = dim

  CALL h5screate_simple_f(arank, adims, aspace_id, error)
# if defined (DOUBLE_PRECISION)
  CALL h5tcopy_f(H5T_NATIVE_DOUBLE, atype_id, error)
# else
  CALL h5tcopy_f(H5T_NATIVE_REAL, atype_id, error)
# endif
  CALL h5acreate_f(id, p_name, atype_id, aspace_id, attr_id, error)
  CALL h5awrite_f(attr_id, atype_id, p_value, data_dims, error)
  CALL h5aclose_f(attr_id, error)
  CALL h5tclose_f(atype_id, error)
  CALL h5sclose_f(aspace_id, error)

End Subroutine add_Parameter_Float

Subroutine Easy_HDF5
Implicit None
character(80) :: py_file
py_file = trim(RESULT_FOLDER) // 'easyhdf5.py'
if (myid .eq. 0) Then
  open(4533, file=py_file, status='unknown')
  Write(4533,'(A)') "import numpy as np"
  Write(4533,'(A)') "import h5py"
  Write(4533,'(A)') "f=h5py.File('field.h5')"
  Write(4533,'(A)') "vars = f.keys()"
  Write(4533,'(A)') "import os"
  Write(4533,'(A)') "from os import listdir"
  Write(4533,'(A)') "from os.path import isfile, join"
  Write(4533,'(A)') "mypath = '.'"
  Write(4533,'(A)') "files = [f for f in listdir(mypath) if isfile(join(mypath, f))]"
  Write(4533,'(A)') "for name in files:"
  Write(4533,'(A)') "    for item in vars:"
  Write(4533,'(A)') "        if item.lower() in name.lower():"
  Write(4533,'(A)') "            f2 = np.loadtxt(name)"
  Write(4533,'(A)') "            new_name = name"
  Write(4533,'(A)') "            f[item].update({new_name:f2})"
  Write(4533,'(A)') "for name in files:"
  Write(4533,'(A)') "    for item in vars:"
  Write(4533,'(A)') "        if item.lower() in name.lower():"
  Write(4533,'(A)') "            os.system('rm -f '+ name )"
  Write(4533,'(A)') "print('Successfully collect all data to field.h5')"
  close(4533)
  call system("cd " // trim(RESULT_FOLDER) // " && python easyhdf5.py")
End If
Call mpi_barrier(MPI_COMM_WORLD, ier)
End Subroutine Easy_HDF5

#endif

END MODULE PARALLEL_FIELD_IO
