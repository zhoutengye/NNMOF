Program lid_driven_cavity
  use ModGlobal
  Use ModTools
  Use ModNavierStokes
  Implicit None
  Integer :: i,j,k
  Character(80) :: testname

  Call Init(inputfield=.false.)
  phi = 1.0_sp
  U = 0.0_sp
  V = 0.0_sp
  W = 0.0_sp
  P = 0.0_sp
  Call InitNavierStokes(U, V, W, Phi, P)

  ! set tup velocity
  u_bc%ghost_flag(2,1) = .false.
  u_bc%bound_type(2,1) = 1
  u_bc%bound_value(2,1) = 1.0_sp

  time = tstart
  ! Begin loop
  Do While (time < tend)
    Call SinglePhaseFlow(U, V, W, P)
    Call Monitor(U,V,W)
    Call WriteFieldData
    time = time + dt
    if (myid .eq. 0) print*, "time=", time, "n_iter=", n_iter,"div_max=", div_max
  End Do

  Call Finalize

End Program lid_driven_cavity
