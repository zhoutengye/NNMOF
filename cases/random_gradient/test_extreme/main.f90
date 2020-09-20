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
  ! delta_theta_max = 5.0_sp * MOF_Pi / 180.0_sp  ! 10 degrees

  ! ZY
  sum_iter = 0
  Call Initialize_NN
  MOFNorm => MOFNN
  method = 'MOFZY'
  Call onemethod(method,num_sampling, data, Norm_ZY, tt(1), num_iters(1,:), error(1,:))

  ! Lemoine GN with analytic gradient
  sum_iter = 0
  MOFNorm => MOFLemoine_GaussNewton
  ddx = 1.0_sp
  mof3d_internal_is_analytic_gradient_enabled = .true.
  mof_use_symmetric_reconstruction = .false.
  Call Lemoine_create_cuboid(ddx, LemoinePoly)
  method = 'Lemoine, Gauss-Newton'
  Call onemethod(method,num_sampling, data, Norm_GN, tt(2), num_iters(2,:), error(2,:))

  ! Lemoine BFGS with analytic gradient
  sum_iter = 0
  MOFNorm => MOFLemoine_BFGS
  ddx = 1.0_sp
  mof3d_internal_is_analytic_gradient_enabled = .true.
  mof_use_symmetric_reconstruction = .true.
  method = 'Lemoine, BFGS, analytic'
  Call onemethod(method,num_sampling, data, Norm_BFGS1, tt(3), num_iters(3,:), error(3,:))

  ! Lemoine BFGS with numerical gradient
  MOFNorm => MOFLemoine_BFGS
  ddx = 1.0_sp
  mof_use_symmetric_reconstruction = .false.
  mof3d_internal_is_analytic_gradient_enabled = .false.
  mof3d_use_optimized_centroid = .false.
  method = 'Lemoine, BFGS, numerical'
  Call onemethod(method, num_sampling, data, Norm_BFGS2, tt(4), num_iters(4,:), error(4,:))

  ! Sussman Gauss Newton with numerical gradient
  ! Call MOFInit3d
  ! MOFNorm => MOFSussmanGaussNewton
  ! method = 'Sussman, Gauss-Newton'
  ! Call onemethod(method, num_sampling, data, Norm_Sussman, tt(5), num_iters(5,:), error(5,:))

  ! Call MPI_FINALIZE(ierr)

  Open(10,file='norm_BFGS1.dat')
  Do i = 1, num_sampling
    Write(10,'(3F)') norm_BFGS1(:,i)
  End Do
  Close(10)

  Open(10,file='norm_BFGS2.dat')
  Do i = 1, num_sampling
    Write(10,'(3F)') norm_BFGS2(:,i)
  End Do
  Close(10)

  Open(10,file='norm_GN.dat')
  Do i = 1, num_sampling
    Write(10,'(3F)') norm_GN(:,i)
  End Do
  Close(10)

  Open(10,file='norm_ZY.dat')
  Do i = 1, num_sampling
    Write(10,'(3F)') norm_ZY(:,i)
  End Do
  Close(10)

  Open(10,file='norm_Sussman.dat')
  Do i = 1, num_sampling
    Write(10,'(3F)') norm_Sussman(:,i)
  End Do
  Close(10)

  Open(10,file='norm_analytic.dat')
  Do i = 1, num_sampling
    Write(10,'(3F)') norm_exact(:,i)
  End Do
  Close(10)

end Subroutine compare

Subroutine onemethod(method, num_sampling, data_in, norm_out, times, iter, error)
  Use ModGlobal
  Use ModTools
  Use ModVOF
  Use ModVOFFunc
  use mod_cg3_polyhedron
  use MODSussman
  use variables_mof
  Implicit None
  Character(80), Intent(In) :: method
  Integer, Intent(In) :: num_sampling
  Real(sp), Intent(In) :: data_in(7,num_sampling)
  Real(sp), Intent(Out) :: Norm_out(3,num_sampling)
  ! Real(sp), Intent(Out) :: angle_out(2,num_sampling)
  ! Real(sp), Intent(Out) :: cen_out(2,num_sampling)
  Real(sp), Intent(Out) :: times
  Real(sp), Intent(Out) :: iter(2)
  Real(sp), Intent(Out) :: error(2)
  Real(sp) :: angle_out(2,num_sampling)
  Real(sp) :: cen_out(2,num_sampling)
  Real(sp) :: c_error
  Real(sp) :: angle_error
  Real(sp) :: tot_angle_error = 0.0_sp
  Real(sp) :: Norm_analytic(3,num_sampling)
  Real(sp) :: Norm_analytic1(3)
  Real(sp) :: Angle_analytic(2)
  Integer :: sum_iter(2)
  Integer :: i
  Real(sp) :: f
  Real(sp) :: c(3)
  Real(sp) :: c_cal(3)
  Real(sp) :: norm(3)
  Real(sp) :: angle_init(2), norm_init(3)
  Real(sp) :: angle_error_init(3)
  Real(sp) :: norm_error_init(3)
  Real(sp) :: cen_error_init(3)
  Real(sp) :: f_error_init
  Real(sp) :: t1,t2
  Real(sp) :: max_error
  Real(sp) :: norm_error_max(3)
  Real(sp) :: norm_c_error_max(3)
  Real(sp) :: f_error_max
  Real(sp) :: err_temp

  max_error = 0.0_sp

  sum_iter = 0
  Call cpu_time(t1)
  Do i = 1, num_sampling
    f = data_in(4,i)
    c = data_in(5:7,i)
    Call MOFNorm(f,c,norm)
    norm_out(:,i) = norm
    sum_iter = sum_iter + mof_niter
  End Do
  Call cpu_time(t2)
  times = t2 - t1
  iter = dble(sum_iter) / dble(num_sampling)
  error = 0.0_sp
  norm_analytic = data_in(1:3,:)
  angle_error = 0.0_sp
  c_error = 0.0_sp
  angle_error_init = 0.0_sp
  Do i = 1, num_sampling
    f = data_in(4,i)
    c = data_in(5:7,i)
    Call Norm2Angle(angle_out(:,i),  norm_out(:,i))
    Call FloodSZ_BackwardC(norm_out(:,i),data_in(4,i),c_cal)
    f = data_in(4,i)
    c = data_in(5:7,i)
    Norm_analytic1 = norm_analytic(:,i)
    error(1) = error(1) + sum(abs(norm_out(:,i)-norm_analytic(:,i)))
    error(2) = error(2) + sum((norm_out(:,i)-norm_analytic(:,i)) * (norm_out(:,i)-norm_analytic(:,i)))
    c_error = c_error + sum(abs(c_cal-data_in(5:7,i)))
    Call Initial_Guess(c-0.5_sp, f, angle_init,err_temp)
    Call Normalization2(Norm_analytic1)
    Call Norm2Angle(angle_analytic,  norm_analytic1)
    angle_error = sum(abs(angle_out(:,i)-angle_analytic))
    tot_angle_error = tot_angle_error + angle_error
    if (f>0.5_sp) cycle
    if ( angle_error > max_error) then
      max_error = angle_error
      norm_error_max = norm_analytic1
      call Angle2Norm(angle_out(:,i),norm_c_error_max)
      f_error_max = f
    end if
    if ( angle_error_init(1) < abs(angle_init(1)-angle_analytic(1)) ) then
      ! angle_error_init(1) = abs(angle_out(1,i)-angle_analytic(1))
      angle_error_init(1) = abs(angle_init(1)-angle_analytic(1))
    endif
    if ( angle_error_init(2) < abs(angle_init(2)-angle_analytic(2)) ) then
      ! angle_error_init(2) = abs(angle_out(2,i)-angle_analytic(2))
      angle_error_init(2) = abs(angle_init(2)-angle_analytic(2))
    endif
    if ( angle_error_init(3) < abs(angle_init(1)-angle_analytic(1))+ abs(angle_init(2)-angle_analytic(2)) ) then
      ! angle_error_init(2) = abs(angle_out(2,i)-angle_analytic(2))
      angle_error_init(3) = abs(angle_init(1)-angle_analytic(1)) + abs(angle_init(2)-angle_analytic(2))
      norm_error_init = norm_analytic(:,i)
    ! Call Norm2Angle(angle_init,  norm_init)
    !  norm_error_init = norm_init
      cen_error_init  = c
      f_error_init  = f
    endif

  End Do

  error = error / dble(num_sampling)
  c_error = c_error / dble(num_sampling)
  tot_angle_error = tot_angle_error / dble(num_sampling)


  print *,'========================================'
  print *,trim(method)
  print *, 'time:', times
  print *, 'num_iter:', iter
  print *, 'error:', error
  print *, 'angle_error:', tot_angle_error
  print *, 'centroid_error:', c_error
  print *, 'initial_guess_error_max:', angle_error_init
  print *, 'initial_norm_error_max:', norm_error_init
  print *, 'initial_cen_error_max:', cen_error_init
  print *, 'f_cen_error_max:', f_error_init
  ! print *,'----------------------------------------'
  ! print *, 'max_error:', max_error
  ! print *, 'vof:', f_error_max
  ! print *, 'exact norm:', norm_error_max
  ! print *, 'calculated norm:', norm_c_error_max
  ! print *,'----------------------------------------'
  print *,'========================================'
  open(10, file = trim(method)//'_info.dat')
  Write(10,*) trim(method)
  Write(10,*)  'time=', times
  Write(10,*)  'num_iter1=', iter(1)
  Write(10,*)  'num_iter2=', iter(2)
  Write(10,*)  'L1_error=', error(1)
  Write(10,*)  'L2_error=', error(2)
  Write(10,*)  'angle_error=', tot_angle_error
  Write(10,*)  'centroid_error=', c_error
  close(10)
End Subroutine onemethod


Subroutine compare_old_initial
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
  ! delta_theta_max = 5.0_sp * MOF_Pi / 180.0_sp  ! 10 degrees

  ! ZY
  ! sum_iter = -4
  ! MOFNorm => MOFZY
  ! method = 'MOFZY'
  ! Call onemethod_old(method,num_sampling, data, Norm_ZY, tt(1), num_iters(1,:), error(1,:))

  ! Lemoine GN with analytic gradient
  sum_iter = 0
  MOFNorm => MOFLemoine_GaussNewton
  ddx = 1.0_sp
  mof3d_internal_is_analytic_gradient_enabled = .true.
  mof_use_symmetric_reconstruction = .false.
  Call Lemoine_create_cuboid(ddx, LemoinePoly)
  method = 'Lemoine, Gauss-Newton'
  Call onemethod_old(method,num_sampling, data, Norm_GN, tt(2), num_iters(2,:), error(2,:))

  ! Lemoine BFGS with analytic gradient
  sum_iter = 0
  MOFNorm => MOFLemoine_BFGS
  ddx = 1.0_sp
  mof3d_internal_is_analytic_gradient_enabled = .true.
  mof_use_symmetric_reconstruction = .true.
  method = 'Lemoine, BFGS, analytic'
  Call onemethod_old(method,num_sampling, data, Norm_BFGS1, tt(3), num_iters(3,:), error(3,:))

  ! Lemoine BFGS with numerical gradient
  MOFNorm => MOFLemoine_BFGS
  ddx = 1.0_sp
  mof_use_symmetric_reconstruction = .false.
  mof3d_internal_is_analytic_gradient_enabled = .false.
  mof3d_use_optimized_centroid = .false.
  method = 'Lemoine, BFGS, numerical'
  Call onemethod_old(method, num_sampling, data, Norm_BFGS2, tt(4), num_iters(4,:), error(4,:))

  ! Sussman Gauss Newton with numerical gradient
  ! Call MOFInit3d
  ! MOFNorm => MOFSussmanGaussNewton
  ! method = 'Sussman, Gauss-Newton'
  ! Call onemethod_old(method, num_sampling, data, Norm_Sussman, tt(5), num_iters(5,:), error(5,:))

  Call MPI_FINALIZE(ierr)

  Open(10,file='norm_BFGS1.dat')
  Do i = 1, num_sampling
    Write(10,'(3F)') norm_BFGS1(:,i)
  End Do
  Close(10)

  Open(10,file='norm_BFGS2.dat')
  Do i = 1, num_sampling
    Write(10,'(3F)') norm_BFGS2(:,i)
  End Do
  Close(10)

  Open(10,file='norm_GN.dat')
  Do i = 1, num_sampling
    Write(10,'(3F)') norm_GN(:,i)
  End Do
  Close(10)

  Open(10,file='norm_ZY.dat')
  Do i = 1, num_sampling
    Write(10,'(3F)') norm_ZY(:,i)
  End Do
  Close(10)

  Open(10,file='norm_Sussman.dat')
  Do i = 1, num_sampling
    Write(10,'(3F)') norm_Sussman(:,i)
  End Do
  Close(10)

  Open(10,file='norm_analytic.dat')
  Do i = 1, num_sampling
    Write(10,'(3F)') norm_exact(:,i)
  End Do
  Close(10)

end Subroutine compare_old_initial

Subroutine onemethod_old(method, num_sampling, data_in, norm_out, times, iter, error)
  Use ModGlobal
  Use ModTools
  Use ModVOF
  Use ModVOFFunc
  use mod_cg3_polyhedron
  use MODSussman
  use variables_mof
  Implicit None
  Character(80), Intent(In) :: method
  Integer, Intent(In) :: num_sampling
  Real(sp), Intent(In) :: data_in(7,num_sampling)
  Real(sp), Intent(Out) :: Norm_out(3,num_sampling)
  ! Real(sp), Intent(Out) :: angle_out(2,num_sampling)
  ! Real(sp), Intent(Out) :: cen_out(2,num_sampling)
  Real(sp), Intent(Out) :: times
  Real(sp), Intent(Out) :: iter(2)
  Real(sp), Intent(Out) :: error(2)
  Real(sp) :: angle_out(2,num_sampling)
  Real(sp) :: cen_out(2,num_sampling)
  Real(sp) :: c_error
  Real(sp) :: angle_error
  Real(sp) :: tot_angle_error = 0.0_sp
  Real(sp) :: Norm_analytic(3,num_sampling)
  Real(sp) :: Angle_analytic(2)
  Integer :: sum_iter(2)
  Integer :: i
  Real(sp) :: f
  Real(sp) :: c(3), c2(3)
  Real(sp) :: c_cal(3)
  Real(sp) :: norm(3), init_norm(3)
  Real(sp) :: t1,t2
  Real(sp) :: max_error
  Real(sp) :: norm_error_max(3)
  Real(sp) :: norm_c_error_max(3)
  Real(sp) :: f_error_max

  max_error = 0.0_sp

  sum_iter = 0
  Call cpu_time(t1)
  Do i = 1, num_sampling
    f = data_in(4,i)
    c = data_in(5:7,i)
    if (f<0.5_sp) then
      init_norm = 0.5_sp-c
    else
      c2 = ( 0.5 - c * f ) / (1.0_sp - f)
      init_norm = 0.5_sp-c2
    end if
    Call MOFNorm(f,c,norm,init_norm)
    ! Call MOFNorm(f,c,norm)
    norm_out(:,i) = norm
    sum_iter = sum_iter + mof_niter
  End Do
  Call cpu_time(t2)
  times = t2 - t1
  iter = dble(sum_iter) / dble(num_sampling)
  error = 0.0_sp
  norm_analytic = data_in(1:3,:)
  angle_error = 0.0_sp
  c_error = 0.0_sp
  Do i = 1, num_sampling
    Call Norm2Angle(angle_out(:,i),  norm_out(:,i))
    Call FloodSZ_BackwardC(norm_out(:,i),data_in(4,i),c_cal)
    f = data_in(4,i)
    c = data_in(5:7,i)
    error(1) = error(1) + sum(abs(norm_out(:,i)-norm_analytic(:,i)))
    error(2) = error(2) + sum((norm_out(:,i)-norm_analytic(:,i)) * (norm_out(:,i)-norm_analytic(:,i)))
    c_error = c_error + sum(abs(c_cal-data_in(5:7,i)))
    Call Norm2Angle(angle_analytic,  norm_analytic(:,i))
    angle_error = sum(abs(angle_out(:,i)-angle_analytic))
    tot_angle_error = tot_angle_error + angle_error
    if ( angle_error > max_error) then
      max_error = angle_error
      norm_error_max = norm_analytic(:,i)
      ! norm_error_max = data_in(5:7,i)
      call Angle2Norm(angle_out(:,i),norm_c_error_max)
      f_error_max = f
    end if
  End Do

  error = error / dble(num_sampling)
  c_error = c_error / dble(num_sampling)
  tot_angle_error = tot_angle_error / dble(num_sampling)


  print *,'========================================'
  print *,trim(method)
  print *, 'time:', times
  print *, 'num_iter:', iter
  print *, 'error:', error
  print *, 'angle_error:', tot_angle_error
  print *, 'centroid_error:', c_error
  ! print *,'----------------------------------------'
  ! print *, 'max_error:', max_error
  ! print *, 'vof:', f_error_max
  ! print *, 'exact norm:', norm_error_max
  ! print *, 'calculated norm:', norm_c_error_max
  ! print *,'----------------------------------------'
  print *,'========================================'
  open(10, file = trim(method)//'_old_info.dat')
  Write(10,*) trim(method)
  Write(10,*)  'time=', times
  Write(10,*)  'num_iter1=', iter(1)
  Write(10,*)  'num_iter2=', iter(2)
  Write(10,*)  'L1_error=', error(1)
  Write(10,*)  'L2_error=', error(2)
  Write(10,*)  'angle_error=', tot_angle_error
  Write(10,*)  'centroid_error=', c_error
  close(10)
End Subroutine onemethod_old
