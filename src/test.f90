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
! test1: IO Primary test
!    (1) Read data set correctly
!    (2) Write data set correctly (see jupyter)
!    (3) Quick visualization test
!===============
Subroutine test1
  Use ModGlobal
  Use ModTools
  Use ModVOF
  Implicit None
  Character(80) :: data_name

  ! (1)
  Call Init(inputfield=.true.)

  ! Check if phi is loaded correctly
   ! print *, phi(6,5,6)
  if ( myid .eq. 0) then
    print *, nl
    print *, dl
    print *, dt
    print *, time
    print *, tend
  endif

  ! block
  !   integer :: j
  ! if (myid .eq. 0) then
  !     do j = 5, 10
  !       ! print *, data_out(:,j,1)
  !       print *, u(5:10,j,7)
  !     end do
  ! end if
  ! end block


  ! (2)
  data_name = 'final'
  Call HDF5WriteData(h5_output_field(1), phi, data_name)
  Call HDF5WriteData(h5_output_field(2),   u, data_name)
  Call HDF5WriteData(h5_output_field(3),   v, data_name)
  Call HDF5WriteData(h5_output_field(4),   w, data_name)

  ! (3)
  ! Call Visual3DContour(f1=phi)
  w = 0
  ! ! for single processor test
  ! w(15:18,15:18,15:18) = 1.0_sp
  ! for 2 2 block mpi test
  if (myid .eq. 3) w(0:3,0:3,0:3) = 1.0_sp
  Call Visual3DContour(f1=phi, f2=w)

  ! ! Test orientation
  ! phi = 0
  ! phi(10,:,:) = 1.0_sp
  ! Call Visual3DContour(f1=phi)
  ! phi = 0
  ! phi(:,10,:) = 1.0_sp
  ! Call Visual3DContour(f1=phi)
  ! phi = 0
  ! phi(:,:,10) = 1.0_sp
  ! Call Visual3DContour(f1=phi)

  Call MPI_FINALIZE(ierr)


end Subroutine test1

!===============
! VOF-PLIC test
!     simple cube advection
!-------------
! setup:
!     20 30 16 grids
!     u=v=w=1.0
!     dx=dy=dz=1.0
!     dx=dy=dz=1.0
!     dt = 0.2
!     tend = 2
!     phi[6:10,6:10,6:10] = 1
!===============
Subroutine test2
  Use ModGlobal
  Use ModVOF
  Use Modtools
  Implicit None
  Real(sp), allocatable, Dimension(:,:,:) :: f_beg
  Real(sp), allocatable, Dimension(:,:,:) :: f_end
  Integer :: nexch(2)

  call init(inputfield=.true.)
  u = 0.0
  v = 0.0
  w = 0.0
  u(1:nl(1),1:nl(2),1:nl(3)) = 1.0
  v(1:nl(1),1:nl(2),1:nl(3)) = 1.0
  w(1:nl(1),1:nl(2),1:nl(3)) = 1.0

  block
    integer :: j
    if (myid .eq. 0) then
      ! do j = 1, nl(2)+1
      ! print *, data_out(:,j,1)
      print *, u(:,5,7)
      ! end do
    end if
  end block

  ! block
  !   integer :: j
  !   if (myid .eq. 0) then
  !     do j = 1, nl(2)+1
  !       ! print *, data_out(:,j,1)
  !       print *, u(5:10,j,7)
  !     end do
  !   end if
  ! end block
  nexch = nl(1:2)
  call updthalo(nexch,1,u)
  call updthalo(nexch,2,u)
  call updthalo(nexch,1,v)
  call updthalo(nexch,2,v)
  call updthalo(nexch,1,w)
  call updthalo(nexch,2,w)
  ! Call u_bc%SetBCS(u)
  ! Call v_bc%SetBCS(v)
  ! Call w_bc%SetBCS(w)
  Allocate(f_beg(0:nl(1)+1,0:nl(2)+1,0:nl(3)+1))
  Allocate(f_end(0:nl(1)+1,0:nl(2)+1,0:nl(3)+1))
  f_beg = phi

  block
    integer :: j
  if (myid .eq. 0) then
      ! do j = 1, nl(2)+1
        ! print *, data_out(:,j,1)
        print *, u(:,5,7)
      ! end do
  end if
  end block

  ! print *, shape(phi)
  ! print *, nl

  ! VOF advection
  Do While (time < tend)
    Call APPLIC(Phi, u, v, w, nl, dl, dt)
    time =  time + dt
    print *, time
  End Do

  f_end = phi

  Call Visual3DContour(f1=f_beg, f2=f_end)
  Call MPI_FINALIZE(ierr)

End Subroutine test2
