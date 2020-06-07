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
! test1: MOF reconstruction
! Given centroid and volume and obtains the normal vecor
!===============
Subroutine test1
  Use ModGlobal
  Use ModTools
  Use ModMOF
  Use Mod_VOF
  Implicit None
  Real(sp) :: f
  Real(sp) :: c(3)
  Real(sp) :: norm(3)

  Call Init(inputfield=.true.)

  Call MOF_Init

  f = 1.0_sp/3.0_sp
  c(1) = 1.0/4.0
  c(2) = 1.0/4.0
  c(3) = 1.0/4.0

  f = 0.8_sp
  c(1) = 0.1
  c(2) = 0.0
  c(3) = 0.0

  f = 0.8_sp
  c(1) = 0.0
  c(2) = 0.1
  c(3) = 0.0

  Call NormSussmanMOF(f,c,norm)

  if (myid .eq. 0) Print *,norm

  Call MPI_FINALIZE(ierr)


end Subroutine test1


!===============
! test1: MOF reconstruction
! Given centroid and volume and obtains the normal vecor
!===============
Subroutine test2
  Use ModGlobal
  Use ModTools
  Use ModMOF
  Use Mod_VOF
  Implicit None
  Real(sp), allocatable, Dimension(:,:,:) :: f_beg
  Real(sp), allocatable, Dimension(:,:,:) :: f_end
  Real(sp) :: v1, v2
  Real(sp) :: v11, v12
  Integer :: nexch(2)

  Call Init(inputfield=.true.)

  Call MOF_Init
  u = 0.0
  v = 0.0
  w = 0.0
  u(1:nl(1),1:nl(2),1:nl(3)) = 1.0
  v(1:nl(1),1:nl(2),1:nl(3)) = 1.0
  w(1:nl(1),1:nl(2),1:nl(3)) = 1.0

  nexch = nl(1:2)
  Call u_bc%SetBCS(u)
  Call v_bc%SetBCS(v)
  Call w_bc%SetBCS(w)
  Call phi_bc%SetBCS(phi)
  Call phi_bc%SetBCS(cx)
  Call phi_bc%SetBCS(cy)
  Call phi_bc%SetBCS(cz)
  Allocate(f_beg(0:nl(1)+1,0:nl(2)+1,0:nl(3)+1))
  Allocate(f_end(0:nl(1)+1,0:nl(2)+1,0:nl(3)+1))
  f_beg = phi

  block
    integer :: j
  if (myid .eq. 0) then
      do j = 5, 10
        ! print *, data_out(:,j,1)
        print *, u(5:10,j,7)
      end do
  end if
  end block
  

  ! VOF advection
  ! Do While (time < tend)
    ! Call VOFCIAM(Phi, u, v, w, nl, dl, dt)
    ! Call MOFCIAM(Phi, cx, cy, cz, u, v, w, nl, dl, dt)
    Call MOFWY(Phi, cx, cy, cz, u, v, w, nl, dl, dt)
    time =  time + dt
  ! End Do

  f_end = phi

  v1 = sum(f_beg(1:nl(1),1:nl(2),1:nl(3))) 
  v2 = sum(f_end(1:nl(1),1:nl(2),1:nl(3))) 

  Call MPI_Reduce(v1, v11, 1, MPI_REAL_SP, MPI_SUM, MPI_COMM_WORLD, 0, ierr)
  Call MPI_Reduce(v1, v12, 1, MPI_REAL_SP, MPI_SUM, MPI_COMM_WORLD, 0, ierr)
  if (myid .eq.0) then
    print *, v11
    print *, v12
  endif

  Call Visual3DContour(f1=f_end)
  ! Call Visual3DContour(f1=f_beg, f2=f_end)


  Call MPI_FINALIZE(ierr)


end Subroutine test2

