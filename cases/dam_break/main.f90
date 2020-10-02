Program lid_driven_cavity
  use ModGlobal
  Use ModTools
  Use ModNavierStokes
  Implicit None
  Integer :: i,j,k
  Character(80) :: testname

  Type(CDSRK3) :: solver

  Call Init(inputfield=.false.)
  Phi = 0.0_sp
  print *, shape(phi)
  ! phi(0:20,0:2,0:40) = 1.0_sp
  ! phi(0:20,0:40,0:2) = 1.0_sp
  phi(0:20,0:40,0:11) = 1.0_sp
  ! phi(0:2,0:20,0:40) = 1.0_sp
  ! phi(0:10,0:10,0:10) = 1.0_sp
  ! phi(0:81,0:40,0:2) = 1.0_sp
  U = 0.0_sp
  V = 0.0_sp
  W = 0.0_sp
  P = 0.0_sp
  ! CSF => CSF_VOSET
  Call Solver%Initialize(U, V, W, Phi, P)

  time = tstart
  ! Begin loop
  Do While (time < tend)
    time = time + dt
    ! Call Solver%TwoPhaseFlow(U, V, W, Phi, P)
    Call Solver%TwoPhaseFlow(U, V, W, Phi, P)
    Call Monitor(U,V,W)
    Call WriteFieldData
    if (myid .eq. 0) print*, "time=", time, "n_iter=", n_iter,"div_max=", div_max, &
      "vol=", sum(phi(1:nl(1),1:nl(2),1:nl(3)))

    ! block 
    !   real(sp) :: div_out(0:nl(1)+1,0:nl(2)+1,0:nl(3)+1)
    !   div_out = 0.0_sp
    !   div_out(1:nl(1),1:nl(2),1:nl(3)) = div
    !   if (time >= 0.008_sp) then
    !     Call visual3dcontour(div_out)
    !   endif
    ! end block

  End Do

    ! Call visual2dcontour(phi,slice_dir=3,slice_coord=1)
    ! Call visual2dcontour(rho,slice_dir=3,slice_coord=1)
    ! Call visual2dcontour(rhox,slice_dir=3,slice_coord=1)
    ! Call visual2dcontour(rhoy,slice_dir=3,slice_coord=1)
    ! Call visual2dcontour(rhoz,slice_dir=3,slice_coord=1)

    Call visual3dcontour(phi)


  Call Finalize

End Program lid_driven_cavity
