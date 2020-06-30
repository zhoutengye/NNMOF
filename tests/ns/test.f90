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
! test1: MISC:
!    (1) Norm2Angle and Angle2Norm
!    (2) Flood_BackwardC
!    (3) FindCentroid
!    (4) MOFZY
!===============
Subroutine test1
  Use ModGlobal
  Use ModTools
  Use ModVOF
  Use ModNavierStokes
  Implicit None
  Integer  :: nn = 0
  Integer :: rank

  Call Init(inputfield=.true.)


  Call InitNavierStokes

  ! VOF advection
  ! Do While (time < tend)
  !   nn = nn + 1
  !   rank = mod(nn+1,3)
  !   Call VOFAdvection(Phi, u, v, w, nl, dl, dt,rank, cx, cy, cz)
  !   time =  time + dt
  ! End Do
  ! print *, nn


  ! Call Visual3DContour(f1=f_end)


  Call MPI_FINALIZE(ierr)


end Subroutine test1

