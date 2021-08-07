#include "param.h"
# if defined TEST1
# define test_case gen_data
# elif defined TEST2
# define test_case compare
# elif defined TEST3
# define test_case compare_old_initial
# endif
program test
  use ModGlobal
  Implicit None

  call test_case

end program test

Subroutine gen_data
  Use ModVOFFunc
  Implicit None
  Integer, Parameter  :: num_sampling = 1000000
  Real(sp) :: n1(num_sampling)
  Real(sp) :: n2(num_sampling)
  Real(sp) :: n3(num_sampling)
  Real(sp) :: f(num_sampling)
  Real(sp) :: data(7,num_sampling)
  Real(sp) :: xc0(3), nr(3)
  Integer :: i
  Call random_number(n1)
  Call random_number(n2)
  Call random_number(n3)
  Call random_number(f)
  n1 = (n1 - 0.5_sp) * 2.0_sp
  n2 = (n2 - 0.5_sp) * 2.0_sp
  n3 = (n3 - 0.5_sp) * 2.0_sp
  open(10,file='fdata.dat',status='old')
  Do i = 1, num_sampling
    read(10,*) f(i)
    nr = (/n1(i),   n2(i),   n3(i)/)
    Call Normalization1(nr)
    Call FloodSZ_BackwardC(nr,f(i),xc0)
    data(1,i) = nr(1)
    data(2,i) = nr(2)
    data(3,i) = nr(3)
    data(4,i) = f(i)
    data(5,i) = xc0(1)
    data(6,i) = xc0(2)
    data(7,i) = xc0(3)
  End Do

  open(10,file='data.dat')
  Do i = 1, num_sampling
    Write(10,'(7F)')data(:,i)
  End Do
  close(10)


End subroutine gen_data

Subroutine compare
  Use ModGlobal
  Use ModTools
  Use ModVOF
  use mod_cg3_polyhedron
  use MODSussman
  use variables_mof
  Implicit None
  Integer,parameter :: num_sampling = 1000000
  Real(sp) :: data(7,num_sampling)
  Real(sp) :: norm_exact(3,num_sampling)
  Real(sp) :: norm_ZY(3,num_sampling)
  Real(sp) :: norm_BFGS1(3,num_sampling)
  Real(sp) :: norm_BFGS2(3,num_sampling)
  Real(sp) :: norm_GN(3,num_sampling)
  Real(sp) :: norm_Sussman(3,num_sampling)
  Real(sp) :: ddx(3)
  Real(sp) :: tt(5)
  Real(sp) :: num_iters(5,2)
  Real(sp) :: error(5,2)
  Integer :: sum_iter(2)
  Integer :: i
  Character(80) :: method


  Call Init(inputfield=.false.)

  ! Load data
  open(10,file='data.dat',status='old')
  Do i = 1,num_sampling
    Read(10,*) data(:,i)
  End Do

  norm_exact = data(1:3,:)

  MOFITERMAX = 100
  mof_tol = 1d-8
  mof3d_tol_derivative = 1d-8
  mof_tol_derivative = 1d-8
  GAUSSNEWTONTOL = 1.0e-8
  delta_theta = MOF_Pi / 180.0_sp  ! 1 degree=pi/180
  delta_theta = 1e-8

  ! ZY
  sum_iter = 0
  Call Initialize_NN
  MOFNorm => MOFNN
  method = 'MOFNN'

  ! ddx = 1.0_sp
  ! ! mof3d_internal_is_analytic_gradient_enabled = .false.
  ! ! mof3d_use_optimized_centroid = .false.
  ! mof3d_internal_is_analytic_gradient_enabled = .true.
  ! mof3d_use_optimized_centroid = .true.
  ! mof_use_symmetric_reconstruction = .false.
  ! Call Lemoine_create_cuboid(ddx, LemoinePoly)
  ! MOFNorm => MOFLemoine_BFGS

  block
    real(8) :: angle_init(2)
    real(8) :: f
    real(8) :: c(3)
    real(8) :: norm(3)
    real(8) :: angle_exact(2)
    real(8) :: err_temp
    real(8) :: tt1, tt2

    real(8) :: norm_cal(3,num_sampling)
    real(8) :: cen_cal(3)
    real(8) :: error_c

    Call cpu_time(tt1)
  Do i = 1, num_sampling
    f = data(4,i)
    c = data(5:7,i)
    Call Norm2Angle(angle_exact,  data(1:3,i))
    Call MOFNorm(f, c, norm_cal(1:3,i))
    ! Call Initial_Guess2(c-0.5_sp, f, angle_init,err_temp)
  End Do
    ! print *,angle_exact
    ! print *,f,c-0.5
    Call cpu_time(tt2)

    error_c = 0.0_sp
  Do i = 1, num_sampling
       Call FloodSZ_BackwardC(norm_cal(1:3,i),data(4,i),cen_cal)
       error_c = error_c + norm2(data(5:7,i)-cen_cal)
  End Do


  print *, error_c / dble(num_sampling)
    print *,tt2-tt1
  end block

end Subroutine compare

