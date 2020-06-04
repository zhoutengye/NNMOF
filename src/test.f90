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
! test1: Primary test
!    (1) Read data set correctly
!    (2) Write data set correctly (see jupyter)
!    (3) cube advection
! setup:
!     20 30 16 grids
!     u=v=w=1.0
!     dx=dy=dz=1.0
!     dx=dy=dz=1.0
!     dt = 0.2
!     tend = 2
!     phi[6:10,6:10,6:10] = 1
!     
!===============
Subroutine test1
  Use ModGlobal
  Use ModTools
  Use ModVOF
  Implicit None
  Character(80) :: data_name

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

  Call Visual3D(phi)

  ! VOF advection
  Do While (time < tend)
    ! call AdvSZ(u, phi, nl, dl, dt, 1)
    ! call AdvSZ(v, phi, nl, dl, dt, 2)
    ! call AdvSZ(w, phi, nl, dl, dt, 3)
    Call APPLIC(Phi, u, v, w, nl, dl, dt)
    time =  time + dt
    print *, time
  End Do

  block
    integer :: j
    Real(sp) :: data_out(nl(1),nl(2),nl(3))
    data_out(:,:,:) = u(1:nl(1),1:nl(2),1:nl(3))
  if (myid .eq. 0) then
      do j = 5, 10
        ! print *, data_out(:,j,1)
        print *, u(5:10,j,7)
      end do
  end if
  end block


  data_name = 'final'
  Call HDF5WriteData(h5_output_field(1), phi, data_name)
  Call HDF5WriteData(h5_output_field(2),   u, data_name)
  Call HDF5WriteData(h5_output_field(3),   v, data_name)
  Call HDF5WriteData(h5_output_field(4),   w, data_name)

  Call MPI_FINALIZE(ierr)

end Subroutine test1
