!================================
! Basic module. Includes:
!     (1) Global variables
!     (2) Initialization
!     (3) Read parameters from namelist
!     (4) HDF5 for field IO
!     (5) MPI initialization and exchange
!     (6) Boundary conditions
!
!   - Note:
!   - (2) Calls (3, 4, 5, 6)
!   -  Seach the above key value to reach the start and end code blocks
!
!================================
# include "param.h"

Module ModGlobal
  Use mpi
  Use hdf5
  Implicit None
  save
  ! Precisoin controls single or double
  Integer, Parameter :: sp = setsp

  !! -----(1) Global variables------------
  !! MPI variables
  Integer  :: myid
  Integer  :: left,right,front,back
  Integer, dimension(2) :: coord
  Integer  :: comm_cart,ierr
  Integer  :: xhalo,yhalo
  Integer  :: status(MPI_STATUS_SIZE)
  Integer  :: blockx, blocky
  Integer  :: nproc
  Integer  :: dims(2)

  !! Computational parameters
  Integer  :: n(3)
  Integer  :: nl(3)
  Real(sp) :: dl(3)
  Real(sp) :: scale = 1
  Real(sp) :: dt
  Real(sp) :: tstart
  Real(sp) :: tend
  Real(sp) :: cfl

  !! Control output
  Real(SP) :: output_inteval
  Integer  :: startframe

  !! HDF5 parameters
  INTEGER(HSIZE_T), DIMENSION(3) :: h5_total_dims
  INTEGER(HSIZE_T), DIMENSION(3) :: h5_block_dims
  INTEGER(HSSIZE_T), DIMENSION(3) :: h5_offset
  Type HDF5File
    Integer(HID_T) :: file_id
    Character(100) :: filename
  End Type HDF5File
  Type HDF5Group
    Integer(HID_T) :: group_id
    Character(100) :: groupname
    Logical        :: Input
    Logical        :: Output
  End Type HDF5Group
  Type(HDF5File)  :: h5_input
  Type(HDF5File)  :: h5_output
  Integer :: n_vars
  Type(HDF5Group), allocatable :: h5_input_field(:)
  Type(HDF5Group), allocatable :: h5_output_field(:)
  Character(80) :: output_name, output_path

  !! Field data
  Type Field
    Character(100) :: name
    Integer :: bound_type(3,2)
    Integer :: bound_value(3,2)
    Integer :: lohi(3,2)
  Contains
    procedure :: SetBC
    procedure :: SetBCS
  End Type Field
  Real(SP), Allocatable :: phi(:,:,:)
  Type(Field)           :: phi_bc
  Real(SP), Allocatable :: u(:,:,:)
  Type(Field)           :: u_bc
  Real(SP), Allocatable :: v(:,:,:)
  Type(Field)           :: v_bc
  Real(SP), Allocatable :: w(:,:,:)
  Type(Field)           :: w_bc
  Real(SP), Allocatable :: cx(:,:,:)
  Type(Field)           :: cx_bc
  Real(SP), Allocatable :: cy(:,:,:)
  Type(Field)           :: cy_bc
  Real(SP), Allocatable :: cz(:,:,:)
  Type(Field)           :: cz_bc

Contains

  !! -----(2) Initialization------------
  Subroutine Init(inputfield)
    Implicit None
    Logical, intent(in) :: inputfield
    ! Initializa MPI
    Call MPI_INIT(ierr)
    Call MPI_COMM_RANK(MPI_COMM_WORLD,myid,ierr)
    Call MPI_COMM_SIZE(MPI_COMM_WORLD,nproc,ierr)
    ! Read parameters from namelist
    Call Init_param
    ! Initialize mpi blocks
    Call InitMPI(n)
    ! Initialize the shape of the block variables
    Call InitShape
    ! Initialize mpi blocks
    Call H5Init(inputfield)
    ! Load the initial field
    Call H5LoadInit
    ! Initialize Boundary conditions
    Call InitBound

  End Subroutine Init
  !! -----End Initialization------------

  !! -----(3) Read parameters from namelist------------
  !======================
  ! Initialize parameters
  !======================
  Subroutine Init_param
    Implicit None
    Integer  :: nx, ny, nz
    Real(sp) :: dx, dy, dz
    Real(sp) :: px, py
    Logical  :: io_phi, io_u, io_v, io_w, io_cx, io_cy, io_cz
    Integer  :: nn
    Character(80) :: input_name
    Character(80) :: file_name

    namelist /mpivar/ px, py
    namelist /gridvar/ nx,ny,nz,dx,dy,dz
    namelist /compvar/ dt, tstart, tend
    namelist /outvar/ output_inteval, startframe
    namelist /iofield/ n_vars, io_phi, io_u, io_v, io_w, io_cx, io_cy, io_cz, output_path, output_name

    Call getarg(1,input_name)
    if (INPUT_NAME .eq. '') Then
      file_name = 'input.namelist'
    Else
      file_name = trim(input_name)//'.namelist'
    endif

    ! Read at processor 0
    if (myid .eq. 0) then
      Open(10, file=file_name)
      Read(10, nml = mpivar)
      Read(10, nml = gridvar)
      Read(10, nml = compvar)
      Read(10, nml = outvar)
      Read(10, nml = iofield)
      dims(1) = px
      dims(2) = py
      n(1) = nx
      n(2) = ny
      n(3) = nz
      dl(1) = dx
      dl(2) = dy
      dl(3) = dz
    end if

    ! Broad values to all processors
    Call MPI_BCAST(dims, 2, MPI_REAL_SP, 0, MPI_COMM_WORLD, ierr)
    Call MPI_barrier(MPI_COMM_WORLD, ierr)

    Call MPI_BCAST(n, 3, MPI_INT, 0, MPI_COMM_WORLD, ierr)
    Call MPI_barrier(MPI_COMM_WORLD, ierr)
    Call MPI_BCAST(dl, 3, MPI_REAL_SP, 0, MPI_COMM_WORLD, ierr)
    Call MPI_barrier(MPI_COMM_WORLD, ierr)

    Call MPI_BCAST(dt, 1, MPI_REAL_SP, 0, MPI_COMM_WORLD, ierr)
    Call MPI_barrier(MPI_COMM_WORLD, ierr)
    Call MPI_BCAST(tstart, 1, MPI_REAL_SP, 0, MPI_COMM_WORLD, ierr)
    Call MPI_barrier(MPI_COMM_WORLD, ierr)
    Call MPI_BCAST(tend, 1, MPI_REAL_SP, 0, MPI_COMM_WORLD, ierr)

    Call MPI_barrier(MPI_COMM_WORLD, ierr)
    Call MPI_BCAST(output_inteval, 1, MPI_REAL_SP, 0, MPI_COMM_WORLD, ierr)
    Call MPI_barrier(MPI_COMM_WORLD, ierr)
    Call MPI_BCAST(startframe, 1, MPI_INT, 0, MPI_COMM_WORLD, ierr)
    Call MPI_barrier(MPI_COMM_WORLD, ierr)

    Call MPI_BCAST(output_name, 80, MPI_CHARACTER, 0, MPI_COMM_WORLD, ierr)
    Call MPI_barrier(MPI_COMM_WORLD, ierr)
    Call MPI_BCAST(output_path, 80, MPI_CHARACTER, 0, MPI_COMM_WORLD, ierr)
    Call MPI_barrier(MPI_COMM_WORLD, ierr)
    Call MPI_BCAST(n_vars, 1, MPI_INT, 0, MPI_COMM_WORLD, ierr)
    Call MPI_barrier(MPI_COMM_WORLD, ierr)
    Call MPI_BCAST(io_phi, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
    Call MPI_barrier(MPI_COMM_WORLD, ierr)
    Call MPI_BCAST(io_u, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
    Call MPI_barrier(MPI_COMM_WORLD, ierr)
    Call MPI_BCAST(io_v, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
    Call MPI_barrier(MPI_COMM_WORLD, ierr)
    Call MPI_BCAST(io_w, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
    Call MPI_barrier(MPI_COMM_WORLD, ierr)
    Call MPI_BCAST(io_cx, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
    Call MPI_barrier(MPI_COMM_WORLD, ierr)
    Call MPI_BCAST(io_cy, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
    Call MPI_barrier(MPI_COMM_WORLD, ierr)
    Call MPI_BCAST(io_cz, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
    Call MPI_barrier(MPI_COMM_WORLD, ierr)

    ! Determine input variables for HDF5
    Allocate(h5_input_field(n_vars))
    Allocate(h5_output_field(n_vars))
    nn = 1
    If (io_phi) Then
      h5_input_field(nn)%groupname = 'phi'
      h5_output_field(nn)%groupname = 'phi'
      nn = nn + 1
    End If
    If (io_u) Then
      h5_input_field(nn)%groupname = 'u'
      h5_output_field(nn)%groupname = 'u'
      nn = nn + 1
    End If
    If (io_v) Then
      h5_input_field(nn)%groupname = 'v'
      h5_output_field(nn)%groupname = 'v'
      nn = nn + 1
    End If
    If (io_w) Then
      h5_input_field(nn)%groupname = 'w'
      h5_output_field(nn)%groupname = 'w'
      nn = nn + 1
    End If
    If (io_cx) Then
      h5_input_field(nn)%groupname = 'cx'
      h5_output_field(nn)%groupname = 'cx'
      nn = nn + 1
    End If
    If (io_cy) Then
      h5_input_field(nn)%groupname = 'cy'
      h5_output_field(nn)%groupname = 'cy'
      nn = nn + 1
    End If
    If (io_cz) Then
      h5_input_field(nn)%groupname = 'cz'
      h5_output_field(nn)%groupname = 'cz'
      nn = nn + 1
    End If

  End Subroutine Init_param
  !! -----End Read parameters from namelist------------


  !! -----(4) HDF5 for field IO------------
  !================
  ! Initialization of HDF5 files, includes
  !   1. Read initial condition from input file
  !   2. Create output file
  !================
  Subroutine H5Init(inputfield)
    use hdf5
    Implicit None
    INTEGER(HID_T) :: plist_id      ! Property list identifier 
    logical :: inputfield

    Integer :: h5error

    h5_total_dims(1) = n(1)
    h5_total_dims(2) = n(2)
    h5_total_dims(3) = n(3)
    h5_block_dims(1) = nl(1)
    h5_block_dims(2) = nl(2)
    h5_block_dims(3) = nl(3)


    h5_offset(1) = int(nl(1)) * coord(1)
    h5_offset(2) = int(nl(2)) * coord(2)
    h5_offset(3) = 0

    Call MPI_barrier(MPI_COMM_WORLD, h5error)

    if ( inputfield ) then 
    !   1. Create input file
    Call h5open_f(h5error)
    Call h5pcreate_f(H5P_FILE_ACCESS_F, plist_id, h5error)
    Call h5pset_fapl_mpio_f(plist_id, MPI_COMM_WORLD, MPI_INFO_NULL, h5error)

    ! Create the file collectively.
    h5_input%filename = "input.h5"
    Call h5fopen_f(h5_input%filename, H5F_ACC_RDONLY_F, h5_input%file_id, h5error, access_prp = plist_id)
    Call h5pclose_f(plist_id, h5error)
    end if

    !   2. Create output file
    Call h5open_f(h5error)
    Call h5pcreate_f(H5P_FILE_ACCESS_F, plist_id, h5error)
    Call h5pset_fapl_mpio_f(plist_id, MPI_COMM_WORLD, MPI_INFO_NULL, h5error)

    ! Create the file collectively.
    h5_output%filename = "output.h5"
    if (myid .eq. 0 ) Call system('rm -f '//trim(h5_output%filename))
    Call h5fcreate_f(h5_output%filename, H5F_ACC_TRUNC_F, h5_output%file_id, h5error, access_prp = plist_id)
    Call h5pclose_f(plist_id, h5error)

  End Subroutine H5Init

  !====================
  ! Load initial conditions from hdf5 file
  ! from the input h5 file with 'init' dataset in each group
  !====================
  Subroutine H5LoadInit
    Implicit None
    Integer :: i
    Character(80) :: data_name
    data_name = 'init'
    do i = 1, n_vars
      if ( Trim(h5_input_field(i)%groupname) .eq. 'phi' ) Then
        Call HDF5OpenGroup(h5_input, h5_input_field(i))
        Call HDF5ReadData(h5_input_field(i), phi, data_name)
      elseif ( Trim(h5_input_field(i)%groupname) .eq. 'u' ) Then
        Call HDF5OpenGroup(h5_input, h5_input_field(i))
        Call HDF5ReadData(h5_input_field(i), u, data_name)
      elseif ( Trim(h5_input_field(i)%groupname) .eq. 'v' ) Then
        Call HDF5OpenGroup(h5_input, h5_input_field(i))
        Call HDF5ReadData(h5_input_field(i), v, data_name)
      elseif ( Trim(h5_input_field(i)%groupname) .eq. 'w' ) Then
        Call HDF5OpenGroup(h5_input, h5_input_field(i))
        Call HDF5ReadData(h5_input_field(i), v, data_name)
      elseif ( Trim(h5_input_field(i)%groupname) .eq. 'cx' ) Then
        Call HDF5OpenGroup(h5_input, h5_input_field(i))
        Call HDF5ReadData(h5_input_field(i), cx, data_name)
      elseif ( Trim(h5_input_field(i)%groupname) .eq. 'cy' ) Then
        Call HDF5OpenGroup(h5_input, h5_input_field(i))
        Call HDF5ReadData(h5_input_field(i), cy, data_name)
      elseif ( Trim(h5_input_field(i)%groupname) .eq. 'cz' ) Then
        Call HDF5OpenGroup(h5_input, h5_input_field(i))
        Call HDF5ReadData(h5_input_field(i), cy, data_name)
      endif
    end do

  End Subroutine H5LoadInit

  Subroutine HDF5OpenGroup(h5file, h5group)
    Implicit None
    Type(hdf5file) :: h5file
    Type(hdf5group) :: h5group
    Integer :: h5error

    CALL h5gopen_f(h5file%file_id, h5group%groupname, h5group%group_id, h5error)

  End Subroutine HDF5OpenGroup

  Subroutine HDF5CreateGroup(h5file, h5group)
    Implicit None
    Type(hdf5file) :: h5file
    Type(hdf5group) :: h5group
    Integer :: h5error

    CALL h5gcreate_f(h5file%file_id, h5group%groupname, h5group%group_id, h5error)

  End Subroutine HDF5CreateGroup

  Subroutine HDF5WriteData(h5group, data, data_name)
    IMPLICIT NONE

    Type(hdf5group) :: h5group
    Real(sp), Intent(In) :: data(:,:,:)
    Real(sp) :: data_out(nl(1),nl(2),nl(3))
    Character(80) :: data_name

    INTEGER(HID_T) :: filespace     ! Dataspace identifier in file 
    INTEGER(HID_T) :: memspace      ! Dataspace identifier in memory
    INTEGER(HID_T) :: plist_id      ! Property list identifier 
    INTEGER(HID_T) :: dset_id       ! Data set identifier 

    integer :: rank = 3

    Integer :: h5error

    ! Create the data space for the  dataset. 
    CALL h5screate_simple_f(rank, h5_total_dims, filespace, h5error)

    ! Create the dataset with default properties.
    CALL h5dcreate_f(h5group%group_id, data_name,HDF5_REAL_SP, filespace, &
        dset_id, h5error)

    CALL h5sclose_f(filespace, h5error)

    ! Create the data space for the  dataset. 
    CALL h5screate_simple_f(rank, h5_block_dims, memspace, h5error)


    ! Select hyperslab in the file.
    !
    CALL h5dget_space_f(dset_id, filespace, h5error)
    CALL h5sselect_hyperslab_f (filespace, H5S_SELECT_SET_F, h5_offset, h5_block_dims, h5error)

    CALL h5pcreate_f(H5P_DATASET_XFER_F, plist_id, h5error)
    CALL h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, h5error)

    data_out(:,:,:) = data(1:nl(1),1:nl(2),1:nl(3))
    CALL h5dwrite_f(dset_id, HDF5_REAL_SP, data_out, h5_total_dims ,h5error, &
        file_space_id = filespace, mem_space_id = memspace, xfer_prp=plist_id)


    ! Close dataspaces.
    CALL h5sclose_f(filespace, h5error)
    CALL h5sclose_f(memspace, h5error)

    ! Close the dataset.
    CALL h5dclose_f(dset_id, h5error)

    Call h5pclose_f(plist_id, h5error)


  end SUBROUTINE HDF5WriteData

  Subroutine HDF5ReadData(h5group, data, data_name)
    IMPLICIT NONE

    Type(hdf5group) :: h5group
    Real(sp), Intent(InOut) :: data(0:nl(1)+1,0:nl(2)+1,0:nl(3)+1)
    Real(sp) :: data_out(nl(1),nl(2),nl(3)) 
    Character(80) :: data_name

    INTEGER(HID_T) :: filespace     ! Dataspace identifier in file 
    INTEGER(HID_T) :: memspace      ! Dataspace identifier in memory
    INTEGER(HID_T) :: plist_id      ! Property list identifier 
    INTEGER(HID_T) :: dset_id       ! Data set identifier 

    integer :: rank = 3

    Integer :: h5error


    data_out = 0.0_sp

    ! Create the data space for the  dataset. 
    CALL h5screate_simple_f(rank, h5_total_dims, filespace, h5error)

    ! Create the dataset with default properties.
    CALL h5dopen_f(h5group%group_id, data_name, dset_id, ierr)
    

    CALL h5sclose_f(filespace, h5error)

    ! Create the data space for the  dataset. 
    CALL h5screate_simple_f(rank, h5_block_dims, memspace, h5error)

    ! Select hyperslab in the file.
    !
    CALL h5dget_space_f(dset_id, filespace, h5error)
    CALL h5sselect_hyperslab_f (filespace, H5S_SELECT_SET_F, h5_offset, h5_block_dims, h5error)

    CALL h5pcreate_f(H5P_DATASET_XFER_F, plist_id, h5error)
    CALL h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, h5error)

    CALL h5dread_f(dset_id, HDF5_REAL_SP, data_out, h5_total_dims ,h5error, &
        file_space_id = filespace, mem_space_id = memspace, xfer_prp=plist_id)
    data(1:nl(1),1:nl(2),1:nl(3)) = data_out(1:nl(1),1:nl(2),1:nl(3))

    ! Close dataspaces.
    CALL h5sclose_f(filespace, h5error)
    CALL h5sclose_f(memspace, h5error)

    ! Close the dataset.
    CALL h5dclose_f(dset_id, h5error)

    Call h5pclose_f(plist_id, h5error)


  end SUBROUTINE HDF5ReadData
  !! -----End HDF5 for field IO------------

  !!------(5) MPI initialization and exchange------
  !======================
  ! Initialize mpi with 2D pencil like blocks
  ! Adopted from CaNS
  !======================
  Subroutine Initmpi(n)
    implicit none
    integer, intent(in), dimension(3) :: n
    integer :: ntx,nty,ntz
    logical, dimension(3) :: periods
    logical :: reorder = .true.
    !
    periods(:) = .false.
    call MPI_CART_CREATE( MPI_COMM_WORLD, 2, dims, &
        periods, reorder, comm_cart, ierr)
    call MPI_CART_COORDS( comm_cart, myid, 2, coord, ierr)
    !
    call MPI_CART_SHIFT(comm_cart,0,1,left,right,ierr)
    call MPI_CART_SHIFT(comm_cart,1,1,front,back,ierr)
    !
    nl(1) = n(1) / dims(1)
    nl(2) = n(2) / dims(2)
    nl(3) = n(3)
    ntx = nl(1)+2
    nty = nl(2)+2
    ntz = nl(3)+2
    !
    ! Definitions of datatypes for velocity and pressure b.c.'s
    ! Note: array(i,j,k) is basically a 1-dim array;
    !       k is the outer and i is the inner loop counter =>
    !         * for fixed i, (j1+1)*(k1+1) blocks of 1 element,
    !           with (i1+1) elements between start and end
    !         * for fixed j, (k1+1) blocks of (i1+1) elements,
    !           with (i1+1)*(j1+1) elements between start and end
    !
    call MPI_TYPE_VECTOR(nty*ntz,1  ,ntx    ,MPI_REAL_SP,xhalo,ierr)
    call MPI_TYPE_VECTOR(ntz    ,ntx,ntx*nty,MPI_REAL_SP,yhalo,ierr)
    call MPI_TYPE_COMMIT(xhalo,ierr)
    call MPI_TYPE_COMMIT(yhalo,ierr)
    return
  End Subroutine Initmpi

  Subroutine InitShape
    Implicit None
    ! Allocate variables
    Allocate(phi(0:nl(1)+1,0:nl(2)+1,0:nl(3)+1))
    Allocate(  u(0:nl(1)+1,0:nl(2)+1,0:nl(3)+1))
    Allocate(  v(0:nl(1)+1,0:nl(2)+1,0:nl(3)+1))
    Allocate(  w(0:nl(1)+1,0:nl(2)+1,0:nl(3)+1))
    Allocate( cx(0:nl(1)+1,0:nl(2)+1,0:nl(3)+1))
    Allocate( cy(0:nl(1)+1,0:nl(2)+1,0:nl(3)+1))
    Allocate( cz(0:nl(1)+1,0:nl(2)+1,0:nl(3)+1))
    phi = 0.0_sp
    u   = 0.0_sp
    v   = 0.0_sp
    w   = 0.0_sp
    cx  = 0.0_sp
    cy  = 0.0_sp
    cz  = 0.0_sp
    End Subroutine InitShape

  !============================
  ! mpi exchange for one direction
  !============================
  subroutine updthalo(n,idir,f)
    implicit none
    integer , dimension(2), intent(in) :: n
    integer , intent(in) :: idir
    Real(sp), Intent(InOut) :: f(0:,0:,0:)
    !integer :: requests(4), statuses(MPI_STATUS_SIZE,4)
    !
    !  this subroutine updates the halos that store info
    !  from the neighboring computational sub-domain
    !
    select case(idir)
    case(1) ! x direction
      call MPI_SENDRECV(f(1     ,0,0),1,xhalo,left ,0, &
                        f(n(1)+1,0,0),1,xhalo,right,0, &
                        comm_cart,status,ierr)
      call MPI_SENDRECV(f(n(1),0,0),1,xhalo,right,0, &
                        f(0   ,0,0),1,xhalo,left ,0, &
                        comm_cart,status,ierr)
         !call MPI_IRECV(p(0     ,0,0),1,xhalo,left ,1, &
         !               comm_cart,requests(2),error)
         !call MPI_IRECV(p(n(1)+1,0,0),1,xhalo,right,0, &
         !               comm_cart,requests(1),error)
         !call MPI_ISSEND(p(n(1),0,0),1,xhalo,right,1, &
         !               comm_cart,requests(4),error)
         !call MPI_ISSEND(p(1   ,0,0),1,xhalo,left ,0, &
         !               comm_cart,requests(3),error)
         !call MPI_WAITALL(4, requests, statuses, error)
    case(2) ! y direction
      call MPI_SENDRECV(f(0,1     ,0),1,yhalo,front,0, &
                        f(0,n(2)+1,0),1,yhalo,back ,0, &
                        comm_cart,status,ierr)
      call MPI_SENDRECV(f(0,n(2),0),1,yhalo,back ,0, &
                        f(0,0   ,0),1,yhalo,front,0, &
                        comm_cart,status,ierr)
         !call MPI_IRECV(p(0,n(2)+1,0),1,yhalo,back ,0, &
         !               comm_cart,requests(1),error)
         !call MPI_IRECV(p(0,0     ,0),1,yhalo,front,1, &
         !               comm_cart,requests(2),error)
         !call MPI_ISSEND(p(0,1   ,0),1,yhalo,front,0, &
         !               comm_cart,requests(3),error)
         !call MPI_ISSEND(p(0,n(2),0),1,yhalo,back ,1, &
         !               comm_cart,requests(4),error)
         !call MPI_WAITALL(4, requests, statuses, error)
    end select
    return
  end subroutine updthalo
  !!------End MPI initialization and exchange------
  !!------Boundary conditions---------
  Subroutine InitBound
    Implicit None
    phi_bc%name = 'phi'
    phi_bc%lohi(:,1) = 0
    phi_bc%lohi(:,2) = nl(:)+1

    cx_bc%name = 'cx'
    cx_bc%lohi(:,1) = 0
    cx_bc%lohi(:,2) = nl(:)+1

    cy_bc%name = 'cy'
    cy_bc%lohi(:,1) = 0
    cy_bc%lohi(:,2) = nl(:)+1

    cz_bc%name = 'cz'
    cz_bc%lohi(:,1) = 0
    cz_bc%lohi(:,2) = nl(:)+1

    u_bc%name = 'u'
    u_bc%lohi(:,1) = 0
    u_bc%lohi(:,2) = nl(:)+1
    u_bc%lohi(1,2) = nl(1)

    v_bc%name = 'v'
    v_bc%lohi = 0
    v_bc%lohi(:,2) = nl(:)+1
    v_bc%lohi(2,2) = nl(2)

    w_bc%name = 'w'
    w_bc%lohi = 0
    w_bc%lohi(:,2) = nl(:)+1
    w_bc%lohi(3,2) = nl(3)
    ! At present, the values are not imported from file, but assigned directly here.
    ! For VOF problem, always set to 0 Nuemann
    phi_bc%bound_type(:,:) = 2
    phi_bc%bound_value(:,:) = 0

    cx_bc%bound_type(:,:) = 2
    cx_bc%bound_value(:,:) = 0

    cy_bc%bound_type(:,:) = 2
    cy_bc%bound_value(:,:) = 0

    u_bc%bound_type(:,:) = 2
    u_bc%bound_value(:,:) = 0

    v_bc%bound_type(:,:) = 2
    v_bc%bound_value(:,:) = 0

    w_bc%bound_type(:,:) = 2
    w_bc%bound_value(:,:) = 0

  End Subroutine InitBound

  !=======================
  ! Call all MPI and physical boundary condition
  !=======================
  Subroutine SetBCS(self,f)
    Implicit None
    Class(Field) :: self
    Integer :: nexch(2)
    Real(sp), Intent(InOut) :: f(0:nl(1)+1,0:nl(2)+1,0:nl(3)+1)
    nexch = nl(1:2)
    call updthalo(nexch,1,f)
    call updthalo(nexch,2,f)
    if(left .eq.MPI_PROC_NULL)  Call self%setBC(1, -1, f)
    if(right .eq.MPI_PROC_NULL) Call self%setBC(1, 1, f)
    if(front .eq.MPI_PROC_NULL) Call self%setBC(2, -1, f)
    if(back .eq.MPI_PROC_NULL)  Call self%setBC(2, 1, f)
    Call self%setBC(3, -1, f)
    Call self%setBC(3, 1, f)
  End Subroutine SetBCS

  !===============================
  ! Set bounday values.
  ! Only fixed Dirichelt and Nuemann is supported
  ! Need further test if appleid to problems other than VOF
  !     bound_type:
  !        1: Dilichlet
  !        2: Nuemann
  !     idir:
  !        1,2,3: x,y,z
  !     isign:
  !        -1: means lower boundary
  !         1: means upper boundary
  !===============================
  Subroutine SetBC(self, idir, isign, f)
    Implicit None
    Class(Field) :: self
    Integer, Intent(in) :: idir
    Integer, Intent(in) :: isign
    Real(sp), Intent(InOut) :: f(0:nl(1)+1,0:nl(2)+1,0:nl(3)+1)
    Integer :: ind
    Integer :: hl
    if (isign .eq. -1 ) then
      ind = self%lohi(idir,1)
      hl = 1
    else  if (isign .eq. 1 ) then
      ind = self%lohi(idir,2)
      hl = 2
    else
      if (myid .eq. 0) print *,'Wrong bc isign value, should be -1 or 1'
      Call MPI_FINALIZE(ierr)
    end if
    Select Case(idir)
    Case(1)
      Select Case(self%bound_type(idir,hl))
      Case(1)
        f(ind,:,:) = self%bound_value(idir,hl)
      Case(2)
        f(ind,:,:) = f(ind-isign,:,:) + &
            isign * f(ind-isign,:,:) * self%bound_value(idir,hl) / dl(1)
      End Select
    Case(2)
      Select Case(self%bound_type(idir,hl))
      Case(1)
        f(:,ind,:) = self%bound_value(idir,hl)
      Case(2)
        f(:,ind,:) = f(:,ind-isign,:) + &
            isign * f(:,ind-isign,:) * self%bound_value(idir,hl) / dl(2)
      End Select
    Case(3)
      Select Case(self%bound_type(idir,hl))
      Case(1)
        f(:,:,ind) = self%bound_value(idir,hl)
      Case(2)
        f(:,:,ind) = f(:,:,ind-isign) + &
            isign * f(:,:,ind-isign) * self%bound_value(idir,hl) / dl(3)
      End Select
    End Select
  End Subroutine SetBC
End Module ModGlobal
