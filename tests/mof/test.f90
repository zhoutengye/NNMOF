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
  use mod_cg3_polyhedron
  use MODSussman
  Implicit None
  Real(sp) :: f
  Real(sp) :: c(3)
  Real(sp) :: norm(3)
  Real(sp) :: init_norm(3)
  Real(sp) :: angle(2)
  type(t_polyhedron) :: poly
  Real(sp) :: ddx(3)

  Call Init(inputfield=.false.)

  ! ! ==================Test norm2angle and angle2norm
  ! f = 1.0/6.0_sp
  ! norm = (/ 1.0/3.0_sp, 1.0/3.0_sp, 1.0/3.0_sp /)
  ! norm = (/ -1.0/3.0_sp, 1.0/3.0_sp, 1.0/3.0_sp /)
  ! norm = (/ 1.0/3.0_sp, -1.0/3.0_sp, 1.0/3.0_sp /)
  ! norm = (/ 1.0/3.0_sp, 1.0/3.0_sp, -1.0/3.0_sp /)
  ! norm = (/ 1.0_sp, 0.0/3.0_sp, 0.0/3.0_sp /)
  ! norm = (/ 0.0_sp, 1.0/1.0_sp, 0.0/3.0_sp /)
  ! norm = (/ 0.0_sp, 0.0/3.0_sp, 1.0/1.0_sp /)
  ! ! norm = (/ 0.0_sp, 1.0/3.0_sp, 2.0/3.0_sp /)
  ! ! norm = (/ 0.5_sp, 0.0/1.0_sp, 1.0/2.0_sp /)
  ! norm = (/ -0.1_sp, 0.9/1.0_sp, 1.0/1.0_sp /)
  ! Call Normalization2(norm)
  ! if (myid .eq. 0) Print *, norm
  ! Call norm2angle(angle, norm)
  ! if (myid .eq. 0) Print *, angle
  ! Call angle2norm(angle, norm)
  ! if (myid .eq. 0) Print *, norm

  !  ! ==================Test backward flood algorithm for centroid
  ! f = 1.0/6.0_sp
  ! norm = (/ 1.0/3.0_sp, 1.0/3.0_sp, 1.0/3.0_sp /)
  ! norm = (/ -1.0/3.0_sp, -1.0/3.0_sp, -1.0/3.0_sp /)
  ! norm = (/ 1.0/3.0_sp, -1.0/3.0_sp, -1.0/3.0_sp /)
  ! norm = (/ -1.0/3.0_sp, 1.0/3.0_sp, -1.0/3.0_sp /)
  ! norm = (/ 0.0/1.0_sp, 1.0/1.0_sp, -0.0/3.0_sp /)
  ! norm = (/ 1.0/1.0_sp, 0.0/3.0_sp, -0.0/3.0_sp /)
  ! f = 1.0/2.0_sp
  ! norm = (/ 1.0/2.0_sp, 1.0/2.0_sp, -0.0/3.0_sp /)
  ! norm = (/ -1.0/2.0_sp, -1.0/2.0_sp, -0.0/3.0_sp /)

  ! f = 2.0/3.0_sp
  ! norm = (/ 1.0/3.0_sp, 1.0/3.0_sp, 1.0/3.0_sp /)
  ! ! norm = (/ 1.0/1.0_sp, 0.0/3.0_sp, -0.0/3.0_sp /)
  ! ! norm = (/ -1.0/1.0_sp, 0.0/3.0_sp, -0.0/3.0_sp /)

  ! Call FloodSZ_BackwardC(norm, f, c)
  ! print *, c

  ! ! ==========Test FindCentroid subroutine
  ! f = 1.0/6.0_sp
  ! ! norm = (/ 1.0/3.0_sp, 1.0/3.0_sp, 1.0/3.0_sp /)
  ! norm = (/ 1.0/1.0_sp, 0.0/3.0_sp, -0.0/3.0_sp /)
  ! Call Normalization2(norm)
  ! Call Norm2Angle(angle, norm)
  ! Call Angle2Norm(angle, norm)
  ! Call FindCentroid(angle, f, c)
  ! print *, c

  ! ===============Test MOF reconstruction
  ! f = 2.0/27.0_sp
  ! c = (/ 1.0/8.0_sp, 1.0/6.0_sp, 1.0/4.0_sp /)
  ! init_norm = (/ 1.0/3.0_sp, 1.0/3.0_sp, 1.0/4.0_sp /)

  ! f = 1.0/12.0_sp
  ! c = (/ 1.0/4.0_sp, 1.0/4.0_sp, 1.0/8.0_sp /)
  ! init_norm = (/ 1.0/1.0_sp, 1.0/8.0_sp, 1.0/8.0_sp /)

  f = 1.0/12.0_sp
  c = (/ 1.0/4.0_sp, 1.0/4.0_sp, 1.0/8.0_sp /)

  !  --------------Different initial guess
  ! init_norm = (/ -1.0/3.0_sp, -1.0/3.0_sp, -1.0/3.0_sp /)
  Call Normalization2(init_norm)
  ! Call MOFLemoine(f,c,norm, init_norm)
  ! print *, mof_niter
  ! Call MOFZY(f,c,norm)

  block
    integer :: ii
    integer :: nnn = 100000
    real(8) :: tt1, tt2
    real(8) :: ttt1, ttt2
    real(8) :: nnn1(2), nnn2(2)

    !! ZY
    MOFNorm => MOFZY
    Call cpu_time(ttt1)
    init_norm = (/ 1.0/2.0_sp, 1.0/2.0_sp,  1.0/3.0_sp /)
    Do ii = 1, nnn
      init_norm = init_norm + 0.001
      ! Call MOFZY(f,c,norm)
      Call MOFZY(f,c,norm,init_norm)
      nnn2 = nnn2 + mof_niter
    End Do
    Call cpu_time(ttt2)
    print *, "===ZY======================"
    print *, "Norm", norm
    print *, "CPU time", ttt2-ttt1
    print *, "Averaged iteration", dble(nnn2/dble(nnn))
    print *, ""

    !! Lenoine Gauss Newton
    MOFNorm => MOFLemoine_GaussNewton
    Call cpu_time(tt1)
    init_norm = (/ 1.0/2.0_sp, 1.0/2.0_sp,  1.0/3.0_sp /)
    Do ii = 1, nnn
      init_norm = init_norm + 0.001
      ! Call MOFLemoine_GaussNewton(f,c,norm)
      Call MOFLemoine_GaussNewton(f,c,norm, init_norm)
      nnn1 = nnn1 + mof_niter
    End Do
    Call cpu_time(tt2)
    print *, "===Lemoine Gauss Newton====="
    print *, "Norm", norm
    print *, "CPU time", tt2-tt1
    print *, "Averaged iteration", dble(nnn1/dble(nnn))
    print *, ""


    !! Lenoine BFGS
    MOFNorm => MOFLemoine_BFGS
    ddx = 1.0_sp
    Call Lemoine_create_cuboid(ddx, LemoinePoly)

    Call cpu_time(ttt1)
    nnn1 = 0
    init_norm = (/ 1.0/2.0_sp, 1.0/2.0_sp,  1.0/3.0_sp /)
    Do ii = 1, nnn
      init_norm = init_norm + 0.001
      ! Call MOFLemoine_BFGS(poly,f,c,norm)
      Call MOFLemoine_BFGS(f,c,norm,init_norm)
      nnn1 = nnn1 + mof_niter
    End Do
    Call cpu_time(ttt2)
    print *, "===Lemoine BFGS================="
    print *, "Norm", norm
    print *, "CPU time", ttt2-ttt1
    print *, "Averaged iteration", dble(nnn1/dble(nnn))
    print *, ""


    !! Sussman Gauss Newton
    MOFNorm => MOFSussmanGaussNewton
    Call MOFInit3d
    Call cpu_time(ttt1)
    nnn1 = 0
    ! sussman%mof_verbose    = 1          ! 0 or 1
    init_norm = (/ 1.0/2.0_sp, 1.0/2.0_sp,  1.0/3.0_sp /)
    Do ii = 1, nnn
      init_norm = init_norm + 0.001
      Call MOFSussmanGaussNewton(f,c,norm)
      ! Call MOFSussman(f,c,norm,init_norm)
      nnn1 = nnn1 + mof_niter
    End Do
    Call cpu_time(ttt2)
    Call Normalization1(norm)
    print *, "===Sussman GaussNewton================="
    print *, "Norm", norm
    print *, "CPU time", ttt2-ttt1
    print *, "Averaged iteration", dble(nnn1/dble(nnn))
    print *, ""

  End block

  ! Call MOFLemoine(f,c,norm, init_norm)
  ! print *, norm
  ! c = c - 0.5_sp
  ! Call MOFZY(f,c,norm)
  ! print *, norm
  
  Call Normalization2(norm)

  ! if (myid .eq. 0) Print *,norm

  Call MPI_FINALIZE(ierr)

end Subroutine test1


