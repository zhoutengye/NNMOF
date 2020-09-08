#include "param.h"
# if defined TEST1
# define test_case gen_data
# elif defined TEST2
# define test_case compare
# elif defined TEST3
# define test_case single_test
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
  Real(sp) :: data(10,num_sampling)
  Real(sp) :: xc0(3), nr(3)
  Real(sp) :: angle_init(2), angle_exact(2), delta_angle(2)
  Real(sp) :: err_temp
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
    Call Normalization2(nr)
    Call Norm2Angle(angle_exact,nr)
    Call Initial_Guess3(xc0-0.5_sp, f(i), angle_init,err_temp)
    delta_angle = angle_exact - angle_init
    data(1,i) = xc0(1)
    data(2,i) = xc0(2)
    data(3,i) = xc0(3)
    data(4,i) = f(i)
    data(5,i) = angle_exact(1)
    data(6,i) = angle_exact(2)
    data(7,i) = angle_init(1)
    data(8,i) = angle_init(2)
    data(9,i) = delta_angle(1)
    data(10,i) = delta_angle(2)
  End Do

  close(10)

  open(10,file='exact_centroid.dat',status='unknown')
  open(11,file='exact_f.dat',status='unknown')
  open(12,file='exact_angle.dat',status='unknown')
  open(13,file='initial_angle.dat',status='unknown')
  open(14,file='delta_angle.dat',status='unknown')
  Do i = 1, num_sampling
    Write(10,'(3F)')data(1:3,i)
    Write(11,'(F)')data(4,i)
    Write(12,'(2F)')data(5:6,i)
    Write(13,'(2F)')data(7:8,i)
    Write(14,'(2F)')data(9:10,i)
  End Do
  close(10)

End subroutine gen_data

Subroutine single_test
  Use ModVOFFunc
  Implicit None
  Integer, Parameter  :: num_sampling = 1000000
  Real(sp) :: n1(num_sampling)
  Real(sp) :: n2(num_sampling)
  Real(sp) :: n3(num_sampling)
  Real(sp) :: f(num_sampling)
  Real(sp) :: data(10,num_sampling)
  Real(sp) :: xc0(3), nr(3),nr1(3),xc1(3), nr2(3)
  Real(sp) :: angle_init(2), angle_exact(2), delta_angle(2)
  Real(sp) :: err_temp
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
  End Do

  close(10)

  Do i = 1, 1
    nr = (/n1(i),   n2(i),   n3(i)/)
    Call Normalization1(nr)
    Call FloodSZ_BackwardC(nr,f(i),xc0)
    Call Normalization2(nr)
    Call Norm2Angle(angle_exact,nr)
    Call Initial_Guess3(xc0-0.5_sp, f(i), angle_init,err_temp)
    delta_angle = angle_exact - angle_init
    data(1,i) = xc0(1)
    data(2,i) = xc0(2)
    data(3,i) = xc0(3)
    data(4,i) = f(i)
    data(5,i) = angle_exact(1)
    data(6,i) = angle_exact(2)
    data(7,i) = angle_init(1)
    data(8,i) = angle_init(2)
    data(9,i) = delta_angle(1)
    data(10,i) = delta_angle(2)
    Call Angle2Norm(angle_init,nr1)
    Call FloodSZ_BackwardC(nr1,f(i),xc1)
    print *, f(i)
    print *, xc0-0.5_sp
    print *, xc1-0.5_sp
    print *, angle_exact
    Call Angle2Norm(angle_exact,nr2)
    print *, nr2
    print *, nr
  End Do


End subroutine single_test
