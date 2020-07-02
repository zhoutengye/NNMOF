#include "param.h"
# if defined TEST1
# define test_case test1
# elif defined TEST2
# define test_case test2
# endif
program test
  use ModGlobal
  Implicit None

  call test_case

end program test

!===============
! test1:
! (1) Make sure everything reads correctly
! (2) Hydrostatic problem with uniform rho
! (3) Hydrostatic problem with variable rho
!===============
Subroutine test1
  Use ModGlobal
  Use ModTools
  Use ModNavierStokes
  Implicit None
  Integer :: i

  Call Init(inputfield=.true.)
  Call InitNavierStokes(U, V, W, Phi, P)

  ! (1) Test read parameters
  Do i = 1, 1
    If ( myid .eq. i-1) then
      print *, "=====id: ", myid, '=========='
      print *, "rho_l: ", rho_l
      print *, "rho_g: ", rho_g
      print *, "mu_l: ", mu_l
      print *, "mu_g: ", mu_g
      print *, "body_force: ", body_force
      print *, "rk_order: ", rk_order
      print *, "phi_bc_types: ", phi_bc%bound_type
      print *, "phi_bc_values: ", phi_bc%bound_value
      print *, "u_bc_types: ", u_bc%bound_type
      print *, "u_bc_values: ", u_bc%bound_value
      print *, "v_bc_types: ", v_bc%bound_type
      print *, "v_bc_values: ", v_bc%bound_value
      print *, "w_bc_types: ", w_bc%bound_type
      print *, "w_bc_values: ", w_bc%bound_value
    End If
  End Do

  ! (2) Test uniform density
  ! phi = 1.0_sp
  ! body_force = 0.0_sp
  ! body_force(1:3) = - 9.8_sp
  ! Call UpdtRhoMu(Phi)
  ! Call TwoPhaseFlow(U, V, W, Phi, P)
  ! Call Visual3DContour(P)

  ! (3) test variable density
  ! For simplicity only one processor
  ! phi = 0.0_sp
  ! phi(0:10,:,:) = 1.0_sp
  ! phi(:,:,0:10) = 1.0_sp
  ! body_force = 0.0_sp
  ! body_force(3) = - 9.8_sp
  ! Call UpdtRhoMu(Phi)
  ! Call TwoPhaseFlow(U, V, W, Phi, P)
  ! Call Visual3DContour(P)
  
  ! (3) Test variable density
  ! For simplicity only one processor
  phi = 0.0_sp
  cx = 0.0_sp
  cy = 0.0_sp
  cz = 0.0_sp
  body_force = 0.0_sp
  body_force(3) = - 9.8_sp
  if (myid .eq.0) then
    phi(:,:,1:10) = 1.0_sp
     cx(:,:,1:10) = 0.5_sp
     cy(:,:,1:10) = 0.5_sp
     cz(:,:,1:10) = 0.5_sp
  endif
  Call UpdtRhoMu(Phi)
  Do i = 1, 10
    Call TwoPhaseFlow(U, V, W, Phi, P, cx, cy, cz)
  End Do

  Call Visual3DContour(Phi)


  Call MPI_FINALIZE(ierr)


end Subroutine test1

!===============
! test2:
! (1) DAM-BREAK like case
!===============
Subroutine test2
  Use ModGlobal
  Use ModTools
  Use ModNavierStokes
  Implicit None
  Integer :: i, j, k
  Real(sp) :: vol1, vol2

  Call Init(inputfield=.false.)
  Call InitNavierStokes(U, V, W, Phi, P)

  phi = 0.0_sp
  U = 0.0_sp
  V = 0.0_sp
  W = 0.0_sp
  cx = 0.0_sp
  cy = 0.0_sp
  cz = 0.0_sp
  body_force = 0.0_sp
  body_force(3) = - 9.8_sp
  if (myid .eq.0) then
    phi(1:10,1:10,1:10) = 1.0_sp
    cx(1:10,1:10,1:10) = 0.5_sp
    cy(1:10,1:10,1:10) = 0.5_sp
    cz(1:10,1:10,1:10) = 0.5_sp
  endif
  Call UpdtRhoMu(Phi)
  vol1 = 0.0_sp
  Do k = 1,nl(3)
    Do j = 1,nl(2)
      Do i = 1,nl(2)
        vol1 = vol1+phi(i,j,k)
      End Do
    End Do
  End Do

  time = tstart
  Do While (time < tend)
    Call TwoPhaseFlow(U, V, W, Phi, P, cx, cy, cz)
    time = time + dt
    if (myid .eq. 0) print*, "time=", time, "n_iter=", n_iter
  End Do

  vol2 = 0.0_sp
  Do k = 1,nl(3)
    Do j = 1,nl(2)
      Do i = 1,nl(2)
        vol2 = vol2+phi(i,j,k)
      End Do
    End Do
  End Do

  Call Visual3DContour(Phi)

  print *, vol1, vol2

  Call MPI_FINALIZE(ierr)
End Subroutine Test2
