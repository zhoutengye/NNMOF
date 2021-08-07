#include "param.h"
program test
  use ModGlobal
  Use ModTools
  Use ModVOF
  Use ModSussman
  Implicit None
  Real(sp) :: nr(3)
  Real(sp) :: alpha
  Real(sp) :: x0(3), dx(3)
  Real(sp) :: f
  Real(sp) :: xc(3)
  Real(sp) :: cpu1, cpu2
  Integer :: i

  MOFNorm => MOFSussmanGaussNewton
  Call MOFInit3d

  x0 = 0.0_sp
  dx = 1.0_sp

  nr = 1.0_sp/3.0_sp
  alpha = 1.0/3.0_sp
  dx(3) = 0.5_sp
  nr(3) = 2.0_sp/3.0_sp
  alpha = 1.0/4.0_sp
  Call normalization1(nr)

  Call cpu_time(cpu1)
  Do i = 1, 1000000
    ! Call FloodSussman_forwardC(nr,alpha,x0,dx,f,xc)
    Call FloodSZ_forwardC(nr,alpha,x0,dx,f,xc)
  End Do
  Call cpu_time(cpu2)

  ! Call MOFSussmanGaussNewton(vof,alpha,x0,dx,f,xc)

  print *, f
  print *, xc
  print *, cpu2-cpu1

end program test

