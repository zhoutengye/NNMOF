Program Main
  Implicit None
  Integer, Parameter  :: sp = 8
  Real(sp) :: f
  Real(sp) :: c(3)
  Real(sp) :: norm(3)
  Real(sp) :: t1, t2 
  Integer  :: i

  ! alpha = FloodSZ_Backward(norm,vof)
  ! Call FloodSZ_ForwardC(norm,alpha,x0,deltax,vof,c3)

  f = 1.0_sp/6.0_sp
  c(1) = 33.0/48.0
  c(2) = 33.0/48.0
  c(3) = 33.0/48.0

  ! f = 1.0_sp/24.0_sp
  ! c(3) = 1.0/8.0
  ! c(2) = 2.0/4.0
  ! c(1) = 2.0/4.0

  ! f = 0.6_sp
  ! c(1) = 0.5_sp
  ! c(2) = 0.3_sp
  ! c(3) = 0.5_sp

  c = c - 0.5_sp

  Call cpu_time(t1)
  ! Do i = 1, 100000
    ! f = f + 1e-6
    Call MOFZY(f,c,norm)
  ! End Do
  Call cpu_time(t2)
  print *, t2-t1

End Program Main
