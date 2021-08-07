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
  Integer, Parameter  :: num_sampling = 100000
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
  ! n1 = (n1 - 0.5_sp) * 2.0_sp
  ! n2 = (n2 - 0.5_sp) * 2.0_sp
  ! n3 = (n3 - 0.5_sp) * 2.0_sp
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
  Integer,parameter :: num_sampling = 100000
  Real(sp) :: data(7,num_sampling)
  Real(sp) :: norm_init(3), norm_analytic(3)
  Real(sp) :: angle_init(2), angle_analytic(2)
  Real(sp) :: f, c(3)
  Real(sp) :: f2, c2(3)
  Real(sp) :: error,error2,err
  Integer :: sum_iter(2)
  Integer :: i
  Character(80) :: method


  Call Init(inputfield=.false.)

  ! Load data
  open(10,file='data.dat',status='old')
  Do i = 1,num_sampling
    Read(10,*) data(:,i)
  End Do
  close(10)

  Open(10,file='error_1.dat')
  Open(20,file='error_2.dat')

  Do i = 1, num_sampling
    f = data(4,i)
    c = data(5:7,i)
    f2 = 1-f
    c2 = ( 0.5_sp - c * f ) /f2
    norm_analytic = data(1:3,i)
    Call Normalization2(norm_analytic)
    Call Norm2Angle(angle_analytic,  norm_analytic)
    Call Initial_Guess(c-0.5_sp, f, angle_init, err)
    error = abs(angle_init(1)-angle_analytic(1)) + &
        abs(angle_init(2)-angle_analytic(2))
    norm_init = 0.5_sp - c
    Call Norm2Angle(angle_init,  norm_init)
    error2 = abs(angle_init(1)-angle_analytic(1)) + &
        abs(angle_init(2)-angle_analytic(2))
    Write(10,'(4F)') c,error
    Write(10,'(4F)') c2,error
    Write(20,'(4F)') c,error2
    Write(20,'(4F)') c2,error2
  End Do
  close(10)
  close(20)


end Subroutine compare
