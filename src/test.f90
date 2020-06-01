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
! test1: HDF5 parallel input
! Object:
!    (1) Determine the i,j,k rank. Especially between n and dims, fortran and python
!    (2) for a given 2D or 3D shape, able to write correctly
!===============
Subroutine test1
  use ModGlobal
  Implicit None
  real(SP),allocatable :: p(:,:,:)
  Character(80) :: data_name

  Call MPI_INIT(ierr)
  call MPI_COMM_RANK(MPI_COMM_WORLD,myid,ierr)
  call MPI_COMM_SIZE(MPI_COMM_WORLD,nproc,ierr)

  ! ! case 1
  n = (/4,8,2/)
  dims = (/2,2/)
  ! ! case 2
  ! n = (/4,8,2/)
  ! dims = (/4,1/)
  ! ! case 3
  ! n = (/4,8,1/)
  ! dims = (/1,1/)


  Call InitMPI(n)


  Allocate(p(0:nl(1)+1,0:nl(2)+1,0:nl(3)+1))

  p = ToFloat(myid)


  ! do i = 1,4
  !   if (myid .eq. i-1) then
  !     print *, myid
  !     print *,p
  !   end if
  !   call MPI_Barrier(MPI_COMM_WORLD, ierr)
  ! end do


  Call H5Init(inputfield=.false.)
  Allocate(h5_output_field(1))
  h5_output_field%groupname = 'phi'
  print *, h5_output_field(1)%groupname
  Call HDF5CreateGroup(h5_output, h5_output_field(1))

  data_name = 'test'
  Call HDF5WriteData(h5_output_field(1), p, data_name)

  Call MPI_FINALIZE(ierr)

End Subroutine test1

!===============
! test2: HDF5 parallel input
! Object:
!    (1) Read attributes from namelist
!    (2) Read field from h5file
!    (3) Assign correct values to each processor
!===============
Subroutine test2
  use ModGlobal
  Implicit None
  real(SP),allocatable :: p(:,:,:)
  Character(80) :: data_name

  Call MPI_INIT(ierr)
  call MPI_COMM_RANK(MPI_COMM_WORLD,myid,ierr)
  call MPI_COMM_SIZE(MPI_COMM_WORLD,nproc,ierr)


  Call Init_param


  if (myid .eq. 0) print *, '======nx, ny, nz======'
  Call MPI_Barrier(MPI_COMM_WORLD, ierr)
  print *, myid, n
  Call MPI_Barrier(MPI_COMM_WORLD, ierr)

  if (myid .eq. 0) print *, '======dx, dy, dz======'
  Call MPI_Barrier(MPI_COMM_WORLD, ierr)
  print *, myid, dl
  Call MPI_Barrier(MPI_COMM_WORLD, ierr)

  if (myid .eq. 0) print *, '======dt======'
  Call MPI_Barrier(MPI_COMM_WORLD, ierr)
  print *, myid, dt
  Call MPI_Barrier(MPI_COMM_WORLD, ierr)

  if (myid .eq. 0) print *, '======px, py======'
  Call MPI_Barrier(MPI_COMM_WORLD, ierr)
  print *, myid, dims
  Call MPI_Barrier(MPI_COMM_WORLD, ierr)


  Call InitMPI(n)

  Allocate(p(0:nl(1)+1,0:nl(2)+1,0:nl(3)+1))

  if (myid .eq. 0) print *, '======local computaional grid======'
  Call MPI_Barrier(MPI_COMM_WORLD, ierr)
  print *, myid, nl
  Call MPI_Barrier(MPI_COMM_WORLD, ierr)

  if (myid .eq. 0) print *, '======local array size======'
  Call MPI_Barrier(MPI_COMM_WORLD, ierr)
  print *, myid, shape(p)
  Call MPI_Barrier(MPI_COMM_WORLD, ierr)

  Call H5Init(inputfield=.true.)
  Call HDF5OpenGroup(h5_input, h5_input_field(1))

  data_name = 'test'
  Call HDF5ReadData(h5_input_field(1), p, data_name)

  block
    real(sp) :: out_p(nl(1),nl(2),nl(3))
    integer :: i
    out_p(1:nl(1),1:nl(2),1:nl(3)) = p(1:nl(1),1:nl(2),1:nl(3))
    do i = 1,4
      if (myid .eq. i-1) then
        print *, myid
        print *,out_p
      end if
      call MPI_Barrier(MPI_COMM_WORLD, ierr)
    end do
  end block

  ! data_name = 'test'
  ! Call HDF5WriteData(h5_output_phi, p, data_name)

  Call MPI_FINALIZE(ierr)

End Subroutine test2

!===============
! test3: Initialization for all
! Object:
!    (1) Read everything correctly
!===============
Subroutine test3
  use ModGlobal
  Implicit None

  Call Init(inputfield=.true.)

  block
    real(sp) :: out_p(nl(1),nl(2),nl(3))
    integer :: i
    ! out_p(1:nl(1),1:nl(2),1:nl(3)) = phi(1:nl(1),1:nl(2),1:nl(3))
    ! out_p(1:nl(1),1:nl(2),1:nl(3)) = u(1:nl(1),1:nl(2),1:nl(3))
    out_p(1:nl(1),1:nl(2),1:nl(3)) = v%values(1:nl(1),1:nl(2),1:nl(3))
    do i = 1,4
      if (myid .eq. i-1) then
        print *, myid
        print *,out_p
      end if
      call MPI_Barrier(MPI_COMM_WORLD, ierr)
    end do
  end block

  Call MPI_FINALIZE(ierr)

end Subroutine test3


!===============
! test4: boundary conditions
! Object:
!     (1) MPI inner boundary
!     (2) outer boundary
!===============
Subroutine test4
  use ModGlobal
  Implicit None
  Integer :: nexch(2)

 Call Init(inputfield=.true.)

 nexch = nl(1:2)
 Call ExchMPI(nexch,u%values)

  block
    integer :: j,k
    do k = 1,4
      if (myid .eq. k-1) then
        print *, myid
        do j = 0, nl(2)+1
          print *, u%values(:,j,1)
        end do
      end if
      call MPI_Barrier(MPI_COMM_WORLD, ierr)
    end do
  end block

  Call MPI_FINALIZE(ierr)

end Subroutine test4
