#include "param.h"
# if defined TEST1
# define test_case gen_data
# elif defined TEST2
# define test_case compare
# endif
program compare
  Use ModGlobal
  Use ModTools
  Use ModVOF
  use mod_cg3_polyhedron
  use MODSussman
  use variables_mof
  use ModGlobal
  Implicit None

  Integer,parameter :: num_sampling = 170
  Real(sp) :: ddx(3)
  Real(sp) :: f, c(3), norm(3)
  Real(sp) :: ang(2)
  Real(sp) :: cc(3)
  Real(sp) :: norm_e(3)

  mof_tol = 1d-8
  mof3d_tol_derivative = 1d-8
  mof_tol_derivative = 1d-8
  GAUSSNEWTONTOL = 1.0e-8
  delta_theta = MOF_Pi / 180.0_sp  ! 1 degree=pi/180
  delta_theta = 1e-8

  Call Init(inputfield=.false.)

  f = 0.495985203706686
  c(1) = 0.550044028498869
  c(2) = 0.498369754183450
  c(3) = 0.255453397251650
  f = 0.997321933384063
  c(1) = 0.499996001865528
  c(2) = 0.501322010079336
  c(3) = 0.501134503838346
  f = 0.999958169986159
  c(1) = 0.500020832766753
  c(2) = 0.499982261797367
  c(3) = 0.499980171014034
  f = 8.091207358831330E-002
  c(1) = 0.590416624604458
  c(2) = 0.496889598162534
  c(3) = 0.955570448654066
  f = 6.090275906283370E-002
  c(1) = 0.499858052819548
  c(2) = 0.434412006242029
  c(3) = 3.202333031923170E-002
  norm_e(1) = 9.898470959310000E-005
  norm_e(2) = 4.573679092635140E-002
  norm_e(3) = 0.954164224364056
  call normalization2(norm_e)
  print *, norm_e

  MOFNorm => MOFZY
  MOFNorm => MOFLemoine_BFGS
  ddx = 1.0_sp
  mof3d_internal_is_analytic_gradient_enabled = .true.
  Call Lemoine_create_cuboid(ddx, LemoinePoly)

  Call MOFNorm(f,c,norm)
  Call Normalization2(Norm)
  Call norm2angle(ang,norm)
  Call angle2norm(ang,Norm_e)
  print *, ang

  cc(2) = c(1)
  cc(3) = c(2)
  cc(1) = c(3)

  Call MOFNorm(f,cc,norm)
  Call Normalization2(Norm)
  Call norm2angle(ang,norm)
  Call angle2norm(ang,Norm_e)
  print *, ang

  cc(3) = c(1)
  cc(1) = c(2)
  cc(2) = c(3)

  Call MOFNorm(f,cc,norm)
  Call Normalization2(Norm)
  Call norm2angle(ang,norm)
  Call angle2norm(ang,Norm_e)
  print *, ang

  ! c(1:3) = ( 0.5_sp - c(1:3) * f ) / ( 1.0_sp -f )
  ! f = 1.0_sp - f
  ! print *, f
  ! print *, c


  ! Call Initial_Guess(c-0.5_sp, f, norm)

  MOFNorm => MOFLemoine_GaussNewton
  Call MOFNorm(f,c,norm)

  block
    real(sp) :: angle_exact(2)
    Call Norm2Angle(angle_exact,norm)
    Call Normalization2(Norm)
    print *, norm 
    print *, angle_exact
  end block

  read(*,*)

  ! f = 1.0_sp / 48.0_sp
  ! c(1) = 1.0 / 8.0_sp
  ! c(2) = 1.0 / 8.0_sp
  ! c(3) = 1.0 / 8.0_sp


  MOFITERMAX = 10
  mof_tol = 1d-12
  mof3d_tol_derivative = 1d-8
  mof_tol_derivative = 1d-8
  GAUSSNEWTONTOL = 1.0e-8
  delta_theta = 1e-8
  delta_theta_max = 100.0_sp * MOF_Pi / 180.0_sp  ! 10 degrees
  ! ZY
  MOFNorm => MOFZY
  Call MOFNorm(f,c,norm)
  MOFNorm => MOFZY2
  Call MOFNorm(f,c,norm)
  print *,'========================================'
  print *, "mof_ZY"

  print *, 'norm:', norm
  print *, 'num_iter:', mof_niter
  print *,'========================================'

  ! Lemoine GN with analytic gradient
  MOFNorm => MOFLemoine_GaussNewton
  ddx = 1.0_sp
  mof3d_internal_is_analytic_gradient_enabled = .true.
  mof_use_symmetric_reconstruction = .false.
  Call MOFNorm(f,c,norm)
  print *,'========================================'
  print *, "mof_GN Lemoine"

  print *, 'norm:', norm
  print *, 'num_iter:', mof_niter
  print *,'========================================'

  ! Lemoine BFGS with analytic gradient
  MOFNorm => MOFLemoine_BFGS
  ddx = 1.0_sp
  mof3d_internal_is_analytic_gradient_enabled = .true.
  Call Lemoine_create_cuboid(ddx, LemoinePoly)
  Call MOFNorm(f,c,norm)
  print *,'========================================'
  print *, "mof_BFGS Lemoine analytic"

  print *, 'norm:', norm
  print *, 'num_iter:', mof_niter
  print *,'========================================'

!   f = 0.4916324916188766
!   c(1) = 0.6814008872523908
!   c(2) = 0.3493231324273239
!   c(3) =   0.5442845204284332
!   Call MOFNorm(f,c,norm)
!   f = 0.4817198738472945
!   c(1) = 0.5025450849361643
!   c(2) = 0.5146004749267803
!   c(3) = 0.2414947993762771
!   Call MOFNorm(f,c,norm)

! ! 0.0422517135424312      -0.4763812535244112       0.4813670329331576

!   f = 0.4360043067190202
!   c(1) = 0.4842581570824622
!   c(2) = 0.6867550443706089
!   c(3) =   0.3100005077490405
!   Call MOFNorm(f,c,norm)
! ! 0.5839208091039154       0.0041858354255728      -0.4118933554705118

!   f = 0.4987337488593603
!   c(1) = 0.2909413249110092
!   c(2) = 0.4988022168667937
!   c(3) = 0.6178639061746680
!   Call MOFNorm(f,c,norm)
!   print *, 'norm:', norm

  ! Lemoine BFGS with numerical gradient
  MOFNorm => MOFLemoine_BFGS
  ddx = 1.0_sp
  mof3d_internal_is_analytic_gradient_enabled = .false.
  mof_use_symmetric_reconstruction = .false.
  mof3d_use_optimized_centroid = .false.
  Call MOFNorm(f,c,norm)
  print *,'========================================'
  print *, "mof_BFGS Lemoine numerical"

  print *, 'norm:', norm
  print *, 'num_iter:', mof_niter
  print *,'========================================'

  ! Sussman Gauss Newton with numerical gradient
  Call MOFInit3d
  MOFNorm => MOFSussmanGaussNewton

end program compare
