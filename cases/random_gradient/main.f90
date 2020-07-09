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

  f = 9.9370462e-01
  norm = (/2.9198408e-01,   3.0561590e-01,   4.0240002e-01/)

  ! f = 9.7274085e-01
  ! norm = (/6.9810552e-01,   2.0122110e-01,   1.0067338e-01/)

  ! f= 1.9577624e-02
  ! norm = (/5.3086428e-01,   3.0702386e-01,   1.6211186e-01/)

  ! f = 9.8987215e-01
  ! norm = (/4.7135715e-01,   1.8905713e-02,   5.0973713e-01/)


  Call FloodSZ_BackwardC(norm,f,c)

  ! print *, c

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
    nnn2 = 0
    Do ii = 1, nnn
      init_norm = init_norm + 0.001
      Call MOFZY(f,c,norm)
      ! Call MOFZY(f,c,norm,init_norm)
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
    nnn1 = 0
    Do ii = 1, nnn
      init_norm = init_norm + 0.001
      Call MOFLemoine_GaussNewton(f,c,norm)
      ! Call MOFLemoine_GaussNewton(f,c,norm, init_norm)
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
      Call MOFLemoine_BFGS(f,c,norm)
      ! Call MOFLemoine_BFGS(f,c,norm,init_norm)
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
