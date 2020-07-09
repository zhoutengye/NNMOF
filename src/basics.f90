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
  Integer  :: left,right,front,back,top,bottom
  Integer, dimension(3) :: coord
  logical, dimension(3) :: periods
  Integer  :: comm_cart,ierr
  Integer  :: xhalo,yhalo,zhalo
  Integer  :: status(MPI_STATUS_SIZE)
  Integer  :: blockx, blocky
  Integer  :: nproc
  Integer  :: dims(3)
  Integer  :: nglevel = 1 ! Now only support one layer of ghost
  Integer  :: nlo(3) = 0
  Integer  :: nhi(3) = 0

  !! Computational parameters
  Integer  :: n(3)
  Integer  :: nl(3)
  Real(sp) :: dl(3)
  Real(sp) :: scale = 1
  Real(sp) :: dt
  Real(sp) :: dt0
  Real(sp) :: tstart
  Real(sp) :: time
  Real(sp) :: tend
  Real(sp) :: cfl

  !! Control output
  Real(SP) :: output_inteval
  Real(sp) :: time_out
  Real(sp) :: startoutputtime
  Integer  :: output_type = 0
  Logical  :: hotstart = .false.
  Integer  :: hotstart_type
  Real(sp)  :: hotstart_time

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
    Real(sp) :: bound_value(3,2)
    Integer :: lohi(3,2)
    Logical :: ghost_flag(3,2)
  Contains
    procedure :: SetBC
    procedure :: SetBCS
  End Type Field
  Real(SP), Allocatable :: phi(:,:,:)
  Real(SP), Allocatable :: p(:,:,:)
  Real(SP), Allocatable :: u(:,:,:)
  Real(SP), Allocatable :: v(:,:,:)
  Real(SP), Allocatable :: w(:,:,:)
  Real(SP), Allocatable :: cx(:,:,:)
  Real(SP), Allocatable :: cy(:,:,:)
  Real(SP), Allocatable :: cz(:,:,:)

  Type(Field)           :: phi_bc
  Type(Field)           :: p_bc
  Type(Field)           :: u_bc
  Type(Field)           :: v_bc
  Type(Field)           :: w_bc

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
    ! Initialize HDF5 file
    Call H5Init(inputfield)

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
    Integer  :: px, py, pz
    Logical  :: io_phi = .false.
    Logical  :: io_p   = .false.
    Logical  :: io_u   = .false.
    Logical  :: io_v   = .false.
    Logical  :: io_w   = .false.
    Logical  :: io_cx  = .false.
    Logical  :: io_cy  = .false.
    Logical  :: io_cz  = .false.
    Logical  :: periodx = .false.
    Logical  :: periody = .false.
    Logical  :: periodz = .false.
    Integer  :: nn, i
    logical :: file_exists
    Character(80) :: input_name
    Character(80) :: file_name
    Character(10) :: output_format

    namelist /mpivar/ px, py, pz
    namelist /gridvar/ nx,ny,nz,dx,dy,dz, periodx, periody, periodz
    namelist /compvar/ dt, tstart, tend, hotstart, hotstart_type, hotstart_time
    namelist /outvar/ output_inteval, startoutputtime, output_format
    namelist /iofield/ n_vars, io_phi, io_p, io_u, io_v, io_w, io_cx, io_cy, io_cz, output_path, output_name

    Call getarg(1,input_name)
    if (INPUT_NAME .eq. '') Then
      file_name = 'input.namelist'
      h5_input%filename = "input.h5"
    Else
      file_name = trim(input_name)//'.namelist'
      h5_input%filename = trim(input_name)//'.h5'
    endif

    ! Read at processor 0
    If (myid .eq. 0) Then
      INQUIRE(FILE=file_name, EXIST=file_exists)
      If ( file_exists ) Then
        Open(10, file=file_name)
        write(6,*) file_exists
        Read(10, nml = mpivar)
        Read(10, nml = gridvar)
        Read(10, nml = compvar)
        Read(10, nml = outvar)
        Read(10, nml = iofield)
        dims(1) = max(1,px)
        dims(2) = max(1,py)
        dims(3) = max(1,pz)
        n(1) = nx
        n(2) = ny
        n(3) = nz
        dl(1) = dx
        dl(2) = dy
        dl(3) = dz
        periods(1) = periodx
        periods(2) = periody
        periods(3) = periodz
        close(10)
        if (trim(output_format) .eq. 'paraview') then
          output_type = 1
        else if (trim(output_format) .eq. 'tecplot') then
          output_type = 2
        endif
        If (dims(1)*dims(2)*dims(3) .ne. nproc) Then
          print *, "======Fatal Error=============================="
          print *, "nx * ny * nz != np, please check"
          print *, "==============================================="
          Call MPI_Finalize(ierr)
          stop
        End If
      Else
        print *, "======Fatal Error=============================="
        print *, Trim(file_name), " does not exist, please check"
        print *, "==============================================="
        Call MPI_Finalize(ierr)
        stop
      End If
    End If

    ! Broad values to all processors
    Call MPI_BCAST(dims, 3, MPI_REAL_SP, 0, MPI_COMM_WORLD, ierr)
    Call MPI_barrier(MPI_COMM_WORLD, ierr)
    Call MPI_BCAST(periods, 3, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
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
    Call MPI_BCAST(hotstart, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
    Call MPI_barrier(MPI_COMM_WORLD, ierr)
    Call MPI_BCAST(hotstart_type, 1, MPI_INT, 0, MPI_COMM_WORLD, ierr)
    Call MPI_barrier(MPI_COMM_WORLD, ierr)
    Call MPI_BCAST(hotstart_time, 1, MPI_REAL_SP, 0, MPI_COMM_WORLD, ierr)
    Call MPI_barrier(MPI_COMM_WORLD, ierr)

    Call MPI_barrier(MPI_COMM_WORLD, ierr)
    Call MPI_BCAST(output_inteval, 1, MPI_REAL_SP, 0, MPI_COMM_WORLD, ierr)
    Call MPI_barrier(MPI_COMM_WORLD, ierr)
    Call MPI_BCAST(startoutputtime, 1, MPI_REAL_SP, 0, MPI_COMM_WORLD, ierr)
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
    Call MPI_BCAST(io_p, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
    Call MPI_barrier(MPI_COMM_WORLD, ierr)

    dt0 = dt
    time_out = startoutputtime

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
    If (io_p) Then
      h5_input_field(nn)%groupname = 'p'
      h5_output_field(nn)%groupname = 'p'
      nn = nn + 1
    End If

    If (nn-1 .ne. n_vars) then
      If (myid .eq. 0) then
        print *, '==========================================='
        print *, 'Wrong n_vars = ', n_vars, ' in namelist, '
        print *, 'the following variables are in io field'
        do i = 1, nn-1
          print *, i, h5_input_field(i)%groupname
        end do
        print *, '==========================================='
      End If
      Call MPI_Finalize(ierr)
      Stop
    End If

  End Subroutine Init_param

  Subroutine Setup_HotStart
    Implicit None
    Character(80) :: input_name
    Character(80) :: h_time
    Integer :: i
    write(h_time,'(F0.6)') hotstart_time

    Call getarg(1,input_name)
    if (INPUT_NAME .eq. '') Then
      input_name = 'input'
    endif
    If ( myid .eq. 0 ) Then
      Open(10,file='hotstart.py')
      Write(10,'(a)') "import h5py"
      Write(10,'(a)') "import numpy as np"
      Write(10,'(a)') "import sys"
      Write(10,'(a)') "import os"
      Write(10,'(a)') "os.system('cp output.h5 output_old.h5')"
      Write(10,'(a)') "f=h5py.File('output_old.h5','r')"
      Write(10,'(a)') "list_vars=f.keys()"
      Write(10,'(a)') "f3=h5py.File('output.h5','r+')"
      Write(10,'(a)') "times = f['u'].keys()"
      Write(10,'(a)') "time_series1 = []"
      Write(10,'(a)') "for item in times:"
      Write(10,'(a)') "    time_series1.append(float(item))"
      Write(10,'(a)') "time_series1 = np.sort(np.array(time_series1))"
      Write(10,'(a)') "time_series2 = ['%.6f' % x for x in time_series1]"
      Write(10,'(a)') "if (len(sys.argv) ==3):"
      Write(10,'(a)') "    input_name = sys.argv[1]+'.h5'"
      Write(10,'(a)') "    dataset = sys.argv[2]"
      Write(10,'(a)') "    if (dataset[0] == '0'):"
	    Write(10,'(a)') "        dataset = dataset[1:]"
      Write(10,'(a)') "    d_frame = np.where(time_series1 == float(dataset))"
      Write(10,'(a)') "    if (d_frame[0][0] < len(time_series2)-1):"
      Write(10,'(a)') "        for item in time_series2[d_frame[0][0]+1:]:"
      Write(10,'(a)') "            if(item[0]=='0'):"
      Write(10,'(a)') "                item = item[1:]"
      Write(10,'(a)') "            for key in list_vars:"
      Write(10,'(a)') "                del f3[key][item]"
      Write(10,'(a)') "else:"
      Write(10,'(a)') "    input_name = '"//trim(input_name)//".h5'"
      Write(10,'(a)') "    dataset = time_series2[-1]"
      Write(10,'(a)') "    if (dataset[0] == '0'):"
      Write(10,'(a)') "        dataset = dataset[1:]"
      Write(10,'(a)') "f2 = h5py.File(input_name,'w')"
      Write(10,'(a)') "for key in f:"
      Write(10,'(a)') "    f2.create_group(key)"
      Write(10,'(a)') "    grp = f2[key]"
      Write(10,'(a)') "    grp.create_dataset('init', data=np.array(f[key][dataset]))"
      Write(10,'(a)') "f.close()"
      Write(10,'(a)') "f2.close()"
      Write(10,'(a)') "f3.close()"
      Write(10,'(a)') "f4=open('starttime.txt','w')"
      Write(10,'(a)') "f4.write(dataset)"
      Write(10,'(a)') "f4.close()"
      Write(10,'(a)') "print('New input file ready')"
      close(10)
      If (hotstart_type == 1) Then
        Call system('python hotstart.py '//trim(input_name)//' '//trim(h_time))
      Else
        Call system('python hotstart.py')
      End If
    End If
    Call MPI_Barrier(MPI_COMM_WORLD, ierr)
    Do i = 1, nproc
      If (myid .eq. 0) Then
        open(10,file='starttime.txt',status='old')
      End If
      Read(10,*) time
      Close(10)
      Call MPI_Barrier(MPI_COMM_WORLD, ierr)
    End Do
    time_out = max(time+output_inteval,startoutputtime)

  End Subroutine Setup_HotStart
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
    h5_offset(3) = int(nl(3)) * coord(3)

    Call MPI_barrier(MPI_COMM_WORLD, h5error)

    If ( hotstart ) Then
      Call Setup_HotStart
    EndIf

    if ( inputfield ) then
      !   1. Open input file
      Call h5open_f(h5error)
      Call h5pcreate_f(H5P_FILE_ACCESS_F, plist_id, h5error)
      Call h5pset_fapl_mpio_f(plist_id, MPI_COMM_WORLD, MPI_INFO_NULL, h5error)
      ! Open the file collectively.
      Call h5fopen_f(h5_input%filename, H5F_ACC_RDONLY_F, h5_input%file_id, h5error, access_prp = plist_id)
      Call h5pclose_f(plist_id, h5error)
    end if

    !   2. Create output file
    Call h5open_f(h5error)
    Call h5pcreate_f(H5P_FILE_ACCESS_F, plist_id, h5error)
    Call h5pset_fapl_mpio_f(plist_id, MPI_COMM_WORLD, MPI_INFO_NULL, h5error)

    h5_output%filename = "output.h5"
    If ( hotstart ) Then
    ! Create the file collectively.
      Call h5fopen_f(h5_output%filename, H5F_ACC_RDONLY_F, h5_output%file_id, h5error, access_prp = plist_id)
      Call h5pclose_f(plist_id, h5error)
    Else
      if (myid .eq. 0 ) Call system('rm -f '//trim(h5_output%filename))
      Call h5fcreate_f(h5_output%filename, H5F_ACC_TRUNC_F, h5_output%file_id, h5error, access_prp = plist_id)
      Call h5pclose_f(plist_id, h5error)
    EndIf

    ! Load the initial field
    Call H5LoadInit(Inputfield)
    If ( hotstart .eqv. .false.) Then
      ! Write attributes to the output file
      Call HDF5_WRITE_PARAMETERS(h5_output)
      Call h5fclose_f(h5_output%file_id, h5error)
    EndIf


  End Subroutine H5Init

  !====================
  ! Load initial conditions from hdf5 file
  ! from the input h5 file with 'init' dataset in each group
  !====================
  Subroutine H5LoadInit(inputfield)
    Implicit None
    logical :: inputfield
    Integer :: i, h5error
    Character(80) :: data_name
    data_name = 'init'
    If (inputfield) Then
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
          Call HDF5ReadData(h5_input_field(i), w, data_name)
        elseif ( Trim(h5_input_field(i)%groupname) .eq. 'cx' ) Then
          Call HDF5OpenGroup(h5_input, h5_input_field(i))
          Call HDF5ReadData(h5_input_field(i), cx, data_name)
        elseif ( Trim(h5_input_field(i)%groupname) .eq. 'cy' ) Then
          Call HDF5OpenGroup(h5_input, h5_input_field(i))
          Call HDF5ReadData(h5_input_field(i), cy, data_name)
        elseif ( Trim(h5_input_field(i)%groupname) .eq. 'cz' ) Then
          Call HDF5OpenGroup(h5_input, h5_input_field(i))
          Call HDF5ReadData(h5_input_field(i), cz, data_name)
        elseif ( Trim(h5_input_field(i)%groupname) .eq. 'p' ) Then
          Call HDF5OpenGroup(h5_input, h5_input_field(i))
          Call HDF5ReadData(h5_input_field(i), p, data_name)
        endif
      end do
    End If

    If ( hotstart ) Then
      do i = 1, n_vars
        if ( Trim(h5_input_field(i)%groupname) .eq. 'phi' ) Then
          Call HDF5OpenGroup(h5_output, h5_output_field(i))
        elseif ( Trim(h5_input_field(i)%groupname) .eq. 'u' ) Then
          Call HDF5OpenGroup(h5_output, h5_output_field(i))
        elseif ( Trim(h5_input_field(i)%groupname) .eq. 'v' ) Then
          Call HDF5OpenGroup(h5_output, h5_output_field(i))
        elseif ( Trim(h5_input_field(i)%groupname) .eq. 'w' ) Then
          Call HDF5OpenGroup(h5_output, h5_output_field(i))
        elseif ( Trim(h5_input_field(i)%groupname) .eq. 'cx' ) Then
          Call HDF5OpenGroup(h5_output, h5_output_field(i))
        elseif ( Trim(h5_input_field(i)%groupname) .eq. 'cy' ) Then
          Call HDF5OpenGroup(h5_output, h5_output_field(i))
        elseif ( Trim(h5_input_field(i)%groupname) .eq. 'cz' ) Then
          Call HDF5OpenGroup(h5_output, h5_output_field(i))
        elseif ( Trim(h5_input_field(i)%groupname) .eq. 'p' ) Then
          Call HDF5OpenGroup(h5_output, h5_output_field(i))
        endif
        Call H5gclose_f(H5_output_field(i)%group_id, h5error)
      end do
    Else
      do i = 1, n_vars
        if ( Trim(h5_input_field(i)%groupname) .eq. 'phi' ) Then
          Call HDF5CreateGroup(h5_output, h5_output_field(i))
        elseif ( Trim(h5_input_field(i)%groupname) .eq. 'u' ) Then
          Call HDF5CreateGroup(h5_output, h5_output_field(i))
        elseif ( Trim(h5_input_field(i)%groupname) .eq. 'v' ) Then
          Call HDF5CreateGroup(h5_output, h5_output_field(i))
        elseif ( Trim(h5_input_field(i)%groupname) .eq. 'w' ) Then
          Call HDF5CreateGroup(h5_output, h5_output_field(i))
        elseif ( Trim(h5_input_field(i)%groupname) .eq. 'cx' ) Then
          Call HDF5CreateGroup(h5_output, h5_output_field(i))
        elseif ( Trim(h5_input_field(i)%groupname) .eq. 'cy' ) Then
          Call HDF5CreateGroup(h5_output, h5_output_field(i))
        elseif ( Trim(h5_input_field(i)%groupname) .eq. 'cz' ) Then
          Call HDF5CreateGroup(h5_output, h5_output_field(i))
        elseif ( Trim(h5_input_field(i)%groupname) .eq. 'p' ) Then
          Call HDF5CreateGroup(h5_output, h5_output_field(i))
        endif
        Call H5gclose_f(H5_output_field(i)%group_id, h5error)
      end do
    End If

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

  Subroutine HDF5WriteData(h5file, h5group, data, data_name)
    IMPLICIT NONE

    Type(hdf5file) :: h5file
    Type(hdf5group) :: h5group
    Real(sp), Intent(In) :: data(0:nl(1)+1,0:nl(2)+1,0:nl(3)+1)
    Real(sp) :: data_out(nl(1),nl(2),nl(3))
    Character(80) :: data_name

    INTEGER(HID_T) :: filespace     ! Dataspace identifier in file 
    INTEGER(HID_T) :: memspace      ! Dataspace identifier in memory
    INTEGER(HID_T) :: plist_id      ! Property list identifier 
    INTEGER(HID_T) :: dset_id       ! Data set identifier 

    integer :: rank = 3

    Integer :: h5error

    ! Open the file collectively.
    Call h5open_f(h5error)
    Call h5pcreate_f(H5P_FILE_ACCESS_F, plist_id, h5error)
    Call h5pset_fapl_mpio_f(plist_id, MPI_COMM_WORLD, MPI_INFO_NULL, h5error)
    Call h5fopen_f(h5file%filename, H5F_ACC_RDWR_F, h5file%file_id, h5error, access_prp = plist_id)
    Call h5pclose_f(plist_id, h5error)

    ! Open group
    Call h5gopen_f(h5file%file_id, h5group%groupname, h5group%group_id, h5error)

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

    ! Close the dataset, space, group and file
    CALL h5dclose_f(dset_id, h5error)
    Call h5pclose_f(plist_id, h5error)
    Call h5gclose_f(h5group%group_id, h5error)
    Call h5fclose_f(h5file%file_id, h5error)

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

  Subroutine HDF5WriteFrame(data_name)
    IMPLICIT NONE
    Character(80) :: data_name
    Integer :: nn

    do nn = 1, n_vars
      Select Case(Trim(h5_output_field(nn)%groupname))
      Case('phi')
        Call HDF5WriteData(h5_output, h5_output_field(nn), phi,data_name)
      Case('p')
        Call HDF5WriteData(h5_output, h5_output_field(nn), p,data_name)
      Case('u')
        Call HDF5WriteData(h5_output, h5_output_field(nn), u,data_name)
      Case('v')
        Call HDF5WriteData(h5_output, h5_output_field(nn), v,data_name)
      Case('w')
        Call HDF5WriteData(h5_output, h5_output_field(nn), w,data_name)
      Case('cx')
        Call HDF5WriteData(h5_output, h5_output_field(nn), cx,data_name)
      Case('cy')
        Call HDF5WriteData(h5_output, h5_output_field(nn), cy,data_name)
      Case('cz')
        Call HDF5WriteData(h5_output, h5_output_field(nn), cz,data_name)
      End Select
    end do

  end SUBROUTINE HDF5WriteFrame

  Subroutine WriteFieldData
    Implicit None
    Character(80) :: frame_name
    If (time >= time_out ) Then
      write(frame_name,'(F0.6)') time_out
      frame_name = trim(frame_name)
      Call HDF5WriteFrame(frame_name)
      If (myid .eq. 0 ) Print *, "Write data ",time_out," to output.h5"
      time_out = time_out + output_inteval
    End If
  End Subroutine WriteFieldData

  SUBROUTINE HDF5_WRITE_PARAMETERS(h5file)
    Implicit None
    Type(hdf5file) :: h5file
    Character(80) :: tmp_name

    ! write the total grid information
    tmp_name = 'nx'
    Call add_Parameter_integer(h5file%file_id, tmp_name, n(1), 1)
    tmp_name = 'ny'
    Call add_Parameter_integer(h5file%file_id, tmp_name, n(2), 1)
    tmp_name = 'nz'
    Call add_Parameter_integer(h5file%file_id, tmp_name, n(3), 1)
    ! write the grid resolution
    tmp_name = 'dx'
    Call add_Parameter_float(h5file%file_id, tmp_name, dl(1), 1)
    tmp_name = 'dy'
    Call add_Parameter_float(h5file%file_id, tmp_name, dl(2), 1)
    tmp_name = 'dz'
    Call add_Parameter_float(h5file%file_id, tmp_name, dl(3), 1)
    tmp_name = 'out_inteval'
    Call add_Parameter_float(h5file%file_id, tmp_name, dl(3), 1)

  END SUBROUTINE HDF5_WRITE_PARAMETERS

  Subroutine add_Parameter_character(id, p_name, p_value)
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

  Subroutine H5TOTecplot
    Implicit None
    If (myid .eq. 0) Then
      open(10, file='gen_tecplot.py',status='unknown')
      Write(10,'(a)') "import numpy as np"
      Write(10,'(a)') "import h5py"
      Write(10,'(a)') "import os"
      Write(10,'(a)') "import sys"
      Write(10,'(a)') "f = h5py.File('output.h5','r')"
      Write(10,'(a)') "list_vars = list(f.keys())"
      Write(10,'(a)') "times = list(f['phi'].keys())"
      Write(10,'(a)') "time_series1 = []"
      Write(10,'(a)') "time_series2 = []"
      Write(10,'(a)') "for item in times:"
      Write(10,'(a)') "    time_series1.append(float(item))"
      Write(10,'(a)') "time_series1 = np.sort(np.array(time_series1))"
      Write(10,'(a)') "time_series2 = ['%.6f' % x for x in time_series1]"
      Write(10,'(a)') "times = []"
      Write(10,'(a)') "for item in time_series2:"
      Write(10,'(a)') "    if item[0] == '0':"
      Write(10,'(a)') "       item = item[1:]"
      Write(10,'(a)') "    times.append(item)"
      Write(10,'(a)') "nx, ny, nz = f['u'][times[0]].shape"
      Write(10,'(a)') "nz = f.attrs['nz'][0]"
      Write(10,'(a)') "ny = f.attrs['ny'][0]"
      Write(10,'(a)') "nx = f.attrs['nx'][0]"
      Write(10,'(a)') "dz = f.attrs['dx'][0]"
      Write(10,'(a)') "dy = f.attrs['dy'][0]"
      Write(10,'(a)') "dx = f.attrs['dz'][0]"
      Write(10,'(a)') "x = np.arange(nx) * dx + dx * 0.5"
      Write(10,'(a)') "y = np.arange(ny) * dy + dy * 0.5"
      Write(10,'(a)') "z = np.arange(nz) * dz + dz * 0.5"
      Write(10,'(a)') "X, Y, Z = np.mgrid[dz*0.5:nz*dx+dz*0.5:dz,"
      Write(10,'(a)') "                   dy*0.5:ny*dy+dy*0.5:dy,"
      Write(10,'(a)') "                   dx*0.5:nx*dz+dx*0.5:dx"
      Write(10,'(a)') "                   ]"
      Write(10,'(a)') "num_xyz = nx * ny * nz"
      Write(10,'(a)') "if (len(sys.argv) ==2):"
      Write(10,'(a)') "    new_dir = sys.argv[1]"
      Write(10,'(a)') "else:"
      Write(10,'(a)') "    new_dir = 'tecplot'"
      Write(10,'(a)') "os.system('mkdir -p ' + new_dir)"
      Write(10,'(a)') "f2 = h5py.File(new_dir + '/tecplot.h5','w')"
      Write(10,'(a)') "datacount = 1"
      Write(10,'(a)') "X, Y, Z = np.meshgrid(x,y,z)"
      Write(10,'(a)') "for group in times:"
      Write(10,'(a)') "    dk = str(datacount).zfill(6)"
      Write(10,'(a)') "    f2.create_group(dk)"
      Write(10,'(a)') "    grp = f2[dk]"
      Write(10,'(a)') "    for dat in list_vars:"
      Write(10,'(a)') "        grp.create_dataset(dat, data=f[dat][group])"
      Write(10,'(a)') "    grp.create_dataset('X', data=X)"
      Write(10,'(a)') "    grp.create_dataset('Y', data=Y)"
      Write(10,'(a)') "    grp.create_dataset('Z', data=Z)"
      Write(10,'(a)') "    datacount+=1"
      Write(10,'(a)') "    print('Write to data group ' + dk)"
      Write(10,'(a)') "f2.close()"
      Write(10,'(a)') "print('Write to    ' + new_dir + '/tecplot.h5      success')"
      close(10)
      Call system ('python gen_tecplot.py')
    End If

  End Subroutine H5ToTecplot

  Subroutine H5TOParaview
    Implicit None
    If (myid .eq. 0) Then
      open(10, file='gen_paraview.py',status='unknown')
      Write(10,'(a)') "from pyevtk.vtk import VtkFile, VtkStructuredGrid"
      Write(10,'(a)') "import numpy as np"
      Write(10,'(a)') "from pyevtk.vtk import VtkGroup"
      Write(10,'(a)') "import h5py"
      Write(10,'(a)') "import os"
      Write(10,'(a)') "import sys"
      Write(10,'(a)') "f = h5py.File('output.h5','r')"
      Write(10,'(a)') "vars = list(f.keys())"
      Write(10,'(a)') "scalar_vars = []"
      Write(10,'(a)') "scalar_vars_string = []"
      Write(10,'(a)') "for item in vars:"
      Write(10,'(a)') "    if (item != 'u' and item != 'v' and item != 'w'):"
      Write(10,'(a)') "        scalar_vars.append(item)"
      Write(10,'(a)') "        scalar_vars_string.append(str(item))"
      Write(10,'(a)') "times = list(f['phi'].keys())"
      Write(10,'(a)') "nz = f.attrs['nx'][0]"
      Write(10,'(a)') "ny = f.attrs['ny'][0]"
      Write(10,'(a)') "nx = f.attrs['nz'][0]"
      Write(10,'(a)') "dz = f.attrs['dx'][0]"
      Write(10,'(a)') "dy = f.attrs['dy'][0]"
      Write(10,'(a)') "dx = f.attrs['dz'][0]"
      Write(10,'(a)') "x = np.arange(nx) * dx + dx * 0.5"
      Write(10,'(a)') "y = np.arange(ny) * dy + dy * 0.5"
      Write(10,'(a)') "z = np.arange(nz) * dz + dz * 0.5"
      Write(10,'(a)') "X, Y, Z = np.mgrid[dz*0.5:nz*dx+dz*0.5:dz,"
      Write(10,'(a)') "                   dy*0.5:ny*dy+dy*0.5:dy,"
      Write(10,'(a)') "                   dx*0.5:nx*dz+dx*0.5:dx"
      Write(10,'(a)') "                   ]"
      Write(10,'(a)') "start, end = (1,1,1), (nx, ny, nz) #Modify 0->1"
      Write(10,'(a)') "if (len(sys.argv) ==2):"
      Write(10,'(a)') "    new_dir = sys.argv[1]"
      Write(10,'(a)') "else:"
      Write(10,'(a)') "    new_dir = 'paraview'"
      Write(10,'(a)') "new_data = new_dir + '/data'"
      Write(10,'(a)') "os.system('mkdir -p ' + new_dir)"
      Write(10,'(a)') "os.system('mkdir -p ' + new_data)"
      Write(10,'(a)') "for step in times:"
      Write(10,'(a)') "    filename=new_data + '/'+ str(step)"
      Write(10,'(a)') "    w = VtkFile(filename, VtkStructuredGrid) #evtk_test0"
      Write(10,'(a)') "    w.openGrid(start = start, end = end)"
      Write(10,'(a)') "    w.openPiece( start = start, end = end)"
      Write(10,'(a)') "    w.openData('Point', scalars = scalar_vars_string, vectors = 'Velocity')"
      Write(10,'(a)') "    for key in scalar_vars:"
      Write(10,'(a)') "        w.addData(str(key),np.array(f[key][step]))"
      Write(10,'(a)') "    vx = np.array(f['u'][step])"
      Write(10,'(a)') "    vy = np.array(f['v'][step])"
      Write(10,'(a)') "    vz = np.array(f['w'][step])"
      Write(10,'(a)') "    w.addData('Velocity', (vx,vy,vz))"
      Write(10,'(a)') "    w.closeData('Point')"
      Write(10,'(a)') "    w.openElement('Points')"
      Write(10,'(a)') "    w.addData('points', (X, Y, Z))"
      Write(10,'(a)') "    w.closeElement('Points')"
      Write(10,'(a)') "    w.closePiece()"
      Write(10,'(a)') "    w.closeGrid()"
      Write(10,'(a)') "    for key in scalar_vars:"
      Write(10,'(a)') "        w.appendData(data = np.array(f[key][step]))"
      Write(10,'(a)') "    w.appendData(data = (vx,vy,vz))"
      Write(10,'(a)') "    w.appendData((X, Y, Z))"
      Write(10,'(a)') "    w.save()"
      Write(10,'(a)') "    print('file: '+filename+' added')"
      Write(10,'(a)') "g = VtkGroup(new_dir+'/group')"
      Write(10,'(a)') "for step in times:"
      Write(10,'(a)') "    g.addFile(filepath = new_data + '/' + str(step)+'.vts', sim_time = float(step))"
      Write(10,'(a)') "g.save()"
      Write(10,'(a)') "print('group file: group.pvd added')"
      close(10)
      Call system ('python gen_paraview.py')
    End If

  End Subroutine h5TOParaview

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
    logical :: reorder = .true.
    !
    call MPI_CART_CREATE( MPI_COMM_WORLD, 3, dims, &
        periods, reorder, comm_cart, ierr)
    call MPI_CART_COORDS( comm_cart, myid, 3, coord, ierr)
    !
    call MPI_CART_SHIFT(comm_cart,0,1,left,right,ierr)
    call MPI_CART_SHIFT(comm_cart,1,1,front,back,ierr)
    call MPI_CART_SHIFT(comm_cart,2,1,bottom,top,ierr)
    !
    nl(1) = n(1) / dims(1)
    nl(2) = n(2) / dims(2)
    nl(3) = n(3) / dims(3)
    ntx = nl(1)+2*nglevel
    nty = nl(2)+2*nglevel
    ntz = nl(3)+2*nglevel
    ! Block
    !   Integer :: i
    !   Do i = 1, nproc
    !     if (myid .eq. i-1 ) Then
    !       print *,'=============='
    !       Print *, "myid:", myid
    !       Print *, "coord:", coord
    !       Print *, "neighbout:", left, right, front, back, bottom, top
    !     end if
    !     Call mpi_barrier(MPI_comm_world, ierr)
    !   end do
    !   call MPI_FINALIZE(ierr)
    !   stop
    ! End Block
    !
    ! Definitions of datatypes for velocity and pressure b.c.'s
    ! Note: array(i,j,k) is basically a 1-dim array;
    !       k is the outer and i is the inner loop counter =>
    !         * for fixed i, (j1+1)*(k1+1) blocks of 1 element,
    !           with (i1+1) elements between start and end
    !         * for fixed j, (k1+1) blocks of (i1+1) elements,
    !           with (i1+1)*(j1+1) elements between start and end
    !
    call MPI_TYPE_VECTOR(nty*ntz,        nglevel,            ntx,MPI_REAL_SP,xhalo,ierr)
    call MPI_TYPE_VECTOR(    ntz,    nglevel*ntx,        ntx*nty,MPI_REAL_SP,yhalo,ierr)
    call MPI_TYPE_VECTOR(      1,nglevel*ntx*nty,nglevel*ntx*nty,MPI_REAL_SP,zhalo,ierr)
    call MPI_TYPE_COMMIT(xhalo,ierr)
    call MPI_TYPE_COMMIT(yhalo,ierr)
    call MPI_TYPE_COMMIT(zhalo,ierr)
    return
  End Subroutine Initmpi

  Subroutine InitShape
    Implicit None
    ! Allocate variables
    Allocate(phi(0:nl(1)+1,0:nl(2)+1,0:nl(3)+1))
    Allocate(  p(0:nl(1)+1,0:nl(2)+1,0:nl(3)+1))
    Allocate(  u(0:nl(1)+1,0:nl(2)+1,0:nl(3)+1))
    Allocate(  v(0:nl(1)+1,0:nl(2)+1,0:nl(3)+1))
    Allocate(  w(0:nl(1)+1,0:nl(2)+1,0:nl(3)+1))
    Allocate( cx(0:nl(1)+1,0:nl(2)+1,0:nl(3)+1))
    Allocate( cy(0:nl(1)+1,0:nl(2)+1,0:nl(3)+1))
    Allocate( cz(0:nl(1)+1,0:nl(2)+1,0:nl(3)+1))
    phi = 0.0_sp
    p   = 0.0_sp
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
    integer , dimension(3), intent(in) :: n
    integer , intent(in) :: idir
    Real(sp), Intent(InOut) :: f(0:,0:,0:)
    ! Real(sp), Intent(InOut) :: f(0:nl(1)+1,0:nl(2)+1,0:nl(3)+1)
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
    case(3) ! z direction
      call MPI_SENDRECV(f(0,0,     1),1,zhalo,bottom,0, &
                        f(0,0,n(3)+1),1,zhalo,top ,0, &
                        comm_cart,status,ierr)
      call MPI_SENDRECV(f(0,0,n(3)),1,zhalo,top ,0, &
                        f(0,0,   0),1,zhalo,bottom,0, &
                        comm_cart,status,ierr)
      ! call MPI_IRECV(p(0,0,n(3)+1),1,zhalo,bottom ,0, &
                    ! comm_cart,requests(1),ierr)
      ! call MPI_IRECV(p(0,0     ,0),1,zhalo,top,1, &
      !               ! comm_cart,requests(2),ierror)
      ! call MPI_ISSEND(p(0,0   ,1),1,zhalo,top,0, &
      !               comm_cart,requests(3),ierr)
      ! call MPI_ISSEND(p(0,0,n(3)),1,zhalo,bottom ,1, &
      !               comm_cart,requests(4),ierror)
      ! call MPI_WAITALL(4, requests, statuses, ierr)
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
    phi_bc%ghost_flag(:,:) = .true.

    u_bc%name = 'u'
    u_bc%lohi(:,1) = 0
    u_bc%lohi(:,2) = nl(:)+1
    u_bc%lohi(1,2) = nl(1)
    u_bc%ghost_flag(:,:) = .true.
    u_bc%ghost_flag(1,:) = .false.

    v_bc%name = 'v'
    v_bc%lohi = 0
    v_bc%lohi(:,2) = nl(:)+1
    v_bc%lohi(2,2) = nl(2)
    v_bc%ghost_flag(:,:) = .true.
    v_bc%ghost_flag(2,:) = .false.

    w_bc%name = 'w'
    w_bc%lohi = 0
    w_bc%lohi(:,2) = nl(:)+1
    w_bc%lohi(3,2) = nl(3)
    w_bc%ghost_flag(:,:) = .true.
    w_bc%ghost_flag(3,:) = .false.

    p_bc%name = 'p'
    p_bc%lohi(:,1) = 0
    p_bc%lohi(:,2) = nl(:)+1
    p_bc%ghost_flag(3,:) = .true.

    ! At present, the values are not imported from file, but assigned directly here.
    ! For VOF problem, always set to 0 Nuemann
    phi_bc%bound_type(:,:) = 2
    phi_bc%bound_value(:,:) = 0

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
    Real(sp), Intent(InOut) :: f(0:nl(1)+1,0:nl(2)+1,0:nl(3)+1)
    Call updthalo(nl,1,f)
    Call updthalo(nl,2,f)
    Call updthalo(nl,3,f)
    If(left   .eq. MPI_PROC_NULL) Call self%setBC(1, -1, f)
    If(right  .eq. MPI_PROC_NULL) Call self%setBC(1,  1, f)
    If(front  .eq. MPI_PROC_NULL) Call self%setBC(2, -1, f)
    If(back   .eq. MPI_PROC_NULL) Call self%setBC(2,  1, f)
    If(bottom .eq. MPI_PROC_NULL) Call self%setBC(3, -1, f)
    If(top    .eq. MPI_PROC_NULL) Call self%setBC(3,  1, f)
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
        If (self%ghost_flag(idir,hl)) Then
          f(ind,:,:) = 2.0_sp * self%bound_value(idir,hl) - f(ind,:,:)
        Else
          f(ind,:,:) = self%bound_value(idir,hl)
        End If
      Case(2)
        f(ind,:,:) = f(ind-isign,:,:) + &
            isign * f(ind-isign,:,:) * self%bound_value(idir,hl) / dl(1)
      Case(3)
        f(ind,:,:) = 2.0_sp * f(ind-isign,:,:) - f(ind-isign-isign,:,:)
      End Select
    Case(2)
      Select Case(self%bound_type(idir,hl))
      Case(1)
        If (self%ghost_flag(idir,hl)) Then
          f(:,ind,:) = 2.0_sp * self%bound_value(idir,hl) - f(:,ind,:)
        Else
          f(:,ind,:) = self%bound_value(idir,hl)
      End If
      Case(2)
        f(:,ind,:) = f(:,ind-isign,:) + &
            isign * f(:,ind-isign,:) * self%bound_value(idir,hl) / dl(2)
      Case(3)
        f(:,ind,:) = 2.0_sp * f(:,ind-isign,:) - f(:,ind-isign-isign,:)
      End Select
    Case(3)
      Select Case(self%bound_type(idir,hl))
      Case(1)
        If (self%ghost_flag(idir,hl)) Then
          f(:,:,ind) = 2.0_sp * self%bound_value(idir,hl) - f(:,:,ind)
        Else
          f(:,:,ind) = self%bound_value(idir,hl)
        EndIf
      Case(2)
        f(:,:,ind) = f(:,:,ind-isign) + &
            isign * f(:,:,ind-isign) * self%bound_value(idir,hl) / dl(3)
      Case(3)
        f(:,:,ind) = 2.0_sp * f(:,:,ind-isign) - f(:,:,ind-isign-isign)
      End Select
    End Select
  End Subroutine SetBC

  ! Subroutine Periodzbound(var)
    ! Implicit None
    ! Real(sp), Intent(InOut) :: var(0:nl(1)+1,0:nl(2)+1,0:nl(3)+1)
    ! var(0:nl(1)+1,0:nl(2)+1,0)       = var(0:nl(1)+1,0:nl(2)+1,nl(3))
    ! var(0:nl(1)+1,0:nl(2)+1,nl(3)+1) = var(0:nl(1)+1,0:nl(2)+1,1)
  ! End Subroutine Periodzbound

  Subroutine Finalize
    Implicit None
    Call WriteFieldData
    ! Output directory
    If (myid .eq. 0) Then
      If (output_type .eq. 1) Then
        Call H5toParaview
      Else If (output_type .eq. 2) Then
        Call H5toTecplot
      End If
    End If

    Call MPI_Finalize(ierr)

  End Subroutine Finalize

End Module ModGlobal
