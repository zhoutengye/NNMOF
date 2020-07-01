#include "param.h"
# if defined TEST1
# define test_case test1
# elif defined TEST2
# define test_case test2
# elif defined TEST3
# define test_case test3
# elif defined TEST4
# define test_case test4
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
  phi(6:15,0:10,6:15) = 1.0_sp
  Call UpdtRhoMu(Phi)
  Do i = 1, 2
    Call TwoPhaseFlow(U, V, W, Phi, P)
  End Do

  Call Visual3DQuiver(U,V,W)
  
  ! Call Visual3DContour(P)
  ! Call Visual3DContour(Phi)
  ! Call Visual3DContour(u)
  ! Call Visual3DContour(v)
  ! Call Visual3DContour(w)

  Call MPI_FINALIZE(ierr)


end Subroutine test1

