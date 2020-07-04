#include "param.h"
# if defined TEST1
# define test_case test1
# elif defined TEST2
# define test_case test2
# elif defined TEST3
# define test_case test3
# endif
program test
  use ModGlobal
  Implicit None
  External test_case

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

  Call Init(inputfield=.false.)
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
  Real(sp) :: vol1, vol2, vol11, vol21

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
    Call TwoPhaseFlow(U, V, W, Phi, P)
    time = time + dt
    if (myid .eq. 0) print*, "time=", time, "n_iter=", n_iter
  End Do
  Call Visual3DContour(phi)

  vol2 = 0.0_sp
  Do k = 1,nl(3)
    Do j = 1,nl(2)
      Do i = 1,nl(2)
        vol2 = vol2+phi(i,j,k)
      End Do
    End Do
  End Do


  Call MPI_Reduce(vol1, vol11, 1, MPI_REAL_SP, MPI_SUM, MPI_COMM_WORLD, 0, ierr)
  Call MPI_Reduce(vol2, vol21, 1, MPI_REAL_SP, MPI_SUM, MPI_COMM_WORLD, 0, ierr)

  print *, vol1, vol2

  Call MPI_FINALIZE(ierr)
End Subroutine Test2

Subroutine test3
  Use ModGlobal
  Use ModTools
  Use ModNavierStokes
  Implicit None
  Real(sp),allocatable :: uvw(:,:,:)
  Integer :: i,j,k

  Call Init(inputfield=.false.)
  phi = 1.0_sp
  U = 0.0_sp
  V = 0.0_sp
  W = 0.0_sp
  P = 0.0_sp
  body_force = 0.0_sp
  Allocate(uvw(0:nl(1)+1,0:nl(2)+1,0:nl(3)+1))
  Call InitNavierStokes(U, V, W, Phi, P)
  Call BC_UVW(U,V,W)

  ! ! Ucase 1
  ! u_bc%ghost_flag(2,1) = .false.
  ! u_bc%bound_type(2,1) = 1
  ! u_bc%bound_value(2,1) = 1.0_sp
  ! ! Ucase 2
  ! u_bc%ghost_flag(2,2) = .false.
  ! u_bc%bound_type(2,2) = 1
  ! u_bc%bound_value(2,2) = 1.0_sp
  ! ! Ucase 3
  ! u_bc%ghost_flag(3,1) = .false.
  ! u_bc%bound_type(3,1) = 1
  ! u_bc%bound_value(3,1) = 1.0_sp
  ! Ucase 4
  ! u_bc%ghost_flag(3,2) = .false.
  ! u_bc%bound_type(3,2) = 1
  ! u_bc%bound_value(3,2) = 1.0_sp

  ! ! Vcase 1
  ! v_bc%ghost_flag(1,1) = .false.
  ! v_bc%bound_type(1,1) = 1
  ! v_bc%bound_value(1,1) = 1.0_sp
  ! ! Vcase 2
  ! v_bc%ghost_flag(1,2) = .false.
  ! v_bc%bound_type(1,2) = 1
  ! v_bc%bound_value(1,2) = 1.0_sp
  ! ! Vcase 3
  ! v_bc%ghost_flag(3,1) = .false.
  ! v_bc%bound_type(3,1) = 1
  ! v_bc%bound_value(3,1) = 1.0_sp
  ! Vcase 4
  ! v_bc%ghost_flag(3,2) = .false.
  ! v_bc%bound_type(3,2) = 1
  ! v_bc%bound_value(3,2) = 1.0_sp

  ! ! Wcase 1
  ! w_bc%ghost_flag(1,1) = .false.
  ! w_bc%bound_type(1,1) = 1
  ! w_bc%bound_value(1,1) = 1.0_sp
  ! ! Wcase 2
  ! w_bc%ghost_flag(1,2) = .false.
  ! w_bc%bound_type(1,2) = 1
  ! w_bc%bound_value(1,2) = 1.0_sp
  ! ! Wcase 3
  ! w_bc%ghost_flag(2,1) = .false.
  ! w_bc%bound_type(2,1) = 1
  ! w_bc%bound_value(2,1) = 1.0_sp
  !! Wcase 4
  w_bc%ghost_flag(2,2) = .false.
  w_bc%bound_type(2,2) = 1
  w_bc%bound_value(2,2) = 1.0_sp

  time = tstart
  Do While (time < tend)
    Call SinglePhaseFlow(U, V, W, P)
    ! Call TwoPhaseFlow(U, V, W, Phi, P)
    Call Monitor(U,V,W)
    time = time + dt
    if (myid .eq. 0) print*, "time=", time, "n_iter=", n_iter,"div_max=", div_max
    ! Call Visual3DQUIVER(v,u,w)
  End Do
  ! Call Visual3DContour(u)
  Do k=1,nl(3)
    Do j=1,nl(2)
      Do i=1,nl(1)
        uvw(i,j,k) = sqrt(u(i,j,k)**2+v(i,j,k)**2+w(i,j,k)**2)
      End Do
    End Do
  End Do
  Call Visual3DCONTOUR(uvw)

  Call MPI_FINALIZE(ierr)
End Subroutine Test3
