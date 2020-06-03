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
    out_p(1:nl(1),1:nl(2),1:nl(3)) = v(1:nl(1),1:nl(2),1:nl(3))
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
  ! Integer :: nexch(2)

 Call Init(inputfield=.true.)

 ! After testing, this line is added to Init
 ! Call InitBound

 ! ! Test Dilichelet 
 ! if all comment, it will use the default Nuemann
 phi_bc%bound_type(:,:) = 1
 phi_bc%bound_value(:,:) = 10
 cx_bc%bound_type(:,:) = 1
 cx_bc%bound_value(:,:) = 10
 cy_bc%bound_type(:,:) = 1
 cy_bc%bound_value(:,:) = 10
 cy_bc%bound_type(:,:) = 1
 cz_bc%bound_value(:,:) = 10
 u_bc%bound_type(:,:) = 1
 u_bc%bound_value(:,:) = 10
 v_bc%bound_type(:,:) = 1
 v_bc%bound_value(:,:) = 10
 w_bc%bound_type(:,:) = 1
 w_bc%bound_value(:,:) = 10

 ! ! Do this test first, then apply the SetBCS fubction
 ! nexch = nl(1:2)
 ! call updthalo(nexch,1,u)
 ! call updthalo(nexch,2,u)

 ! if(left .eq.MPI_PROC_NULL)  Call u_bc%setBC(1, -1, u)
 ! if(right .eq.MPI_PROC_NULL) Call u_bc%setBC(1, 1, u)
 ! if(front .eq.MPI_PROC_NULL) Call u_bc%setBC(2, -1, u)
 ! if(back .eq.MPI_PROC_NULL)  Call u_bc%setBC(2, 1, u)
 ! Call u_bc%setBC(3, -1, u)
 ! Call u_bc%setBC(3, 1, u)

 ! After testing, we should be able to use a simplified Version
 Call u_bc%setBCS(u)
 Call v_bc%setBCS(v)
 Call w_bc%setBCS(w)
 Call phi_bc%setBCS(phi)

  block
    integer :: j,k
    ! If (myid .eq. 0) Print *, 'Dir 1: bc_node: ', phi_bc%lohi(1,1), phi_bc%lohi(1,2)
    ! If (myid .eq. 0) Print *, 'Dir 2: bc_node: ', phi_bc%lohi(2,1), phi_bc%lohi(2,2)
    ! If (myid .eq. 0) Print *, 'Dir 3: bc_node: ', phi_bc%lohi(3,1), phi_bc%lohi(3,2)
    ! do k = 1,4
    !   if (myid .eq. k-1) then
    !     print *, myid
    !     do j = 0, nl(2)+1
    !       print *, phi(:,j,2)
    !       print *, phi(:,j,0)
    !       print *, phi(:,j,nl(3)+1)
    !       print *, phi(:,j,0) - phi(:,j,1)
    !       print *, phi(:,j,nl(3)+1) - phi(:,j,nl(3))
    !     end do
    !   end if
    ! call MPI_Barrier(MPI_COMM_WORLD, ierr)
    ! end do
    ! If (myid .eq. 0) Print *, 'Dir 1: bc_node: ', u_bc%lohi(1,1), u_bc%lohi(1,2)
    ! If (myid .eq. 0) Print *, 'Dir 2: bc_node: ', u_bc%lohi(2,1), u_bc%lohi(2,2)
    ! If (myid .eq. 0) Print *, 'Dir 3: bc_node: ', u_bc%lohi(3,1), u_bc%lohi(3,2)
    ! do k = 1,4
    !   if (myid .eq. k-1) then
    !     print *, myid
    !     do j = 0, nl(2)+1
    !       print *, u(:,j,2)
    !       print *, u(:,j,0)
    !       print *, u(:,j,nl(3)+1)
    !       print *, u(:,j,0) - u(:,j,1)
    !       print *, u(:,j,nl(3)+1) - u(:,j,nl(3))
    !     end do
    !   end if
    !   call MPI_Barrier(MPI_COMM_WORLD, ierr)
    ! end do
    ! If (myid .eq. 0) Print *, 'Dir 1: bc_node: ', v_bc%lohi(1,1), v_bc%lohi(1,2)
    ! If (myid .eq. 0) Print *, 'Dir 2: bc_node: ', v_bc%lohi(2,1), v_bc%lohi(2,2)
    ! If (myid .eq. 0) Print *, 'Dir 3: bc_node: ', v_bc%lohi(3,1), v_bc%lohi(3,2)
    ! do k = 1,4
    !   if (myid .eq. k-1) then
    !     print *, myid
    !     do j = 0, nl(2)+1
    !       print *, v(:,j,2)
    !       print *, v(:,j,0)
    !       print *, v(:,j,nl(3)+1)
    !       print *, v(:,j,0) - v(:,j,1)
    !       print *, v(:,j,nl(3)+1) - v(:,j,nl(3))
    !     end do
    !   end if
    !   call MPI_Barrier(MPI_COMM_WORLD, ierr)
    ! end do
    If (myid .eq. 0) Print *, 'Dir 1: bc_node: ', w_bc%lohi(1,1), w_bc%lohi(1,2)
    If (myid .eq. 0) Print *, 'Dir 2: bc_node: ', w_bc%lohi(2,1), w_bc%lohi(2,2)
    If (myid .eq. 0) Print *, 'Dir 3: bc_node: ', w_bc%lohi(3,1), w_bc%lohi(3,2)
    do k = 1,4
      if (myid .eq. k-1) then
        print *, myid
        do j = 0, nl(2)+1
          print *, w(:,j,2)
          print *, w(:,j,0)
          print *, w(:,j,nl(3))
          print *, w(:,j,nl(3)+1)
          print *, w(:,j,0) - w(:,j,1)
          print *, w(:,j,nl(3)+1) - w(:,j,nl(3))
        end do
      end if
      call MPI_Barrier(MPI_COMM_WORLD, ierr)
    end do
  end block

  Call MPI_FINALIZE(ierr)

end Subroutine test4
