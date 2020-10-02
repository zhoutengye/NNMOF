!===============================================
! Contains several solver.
!    1. solvers are defined as derived types.
!    1. common subroutines are in file ns_common.f90
! Includes:
!     (1) CDSRK3: Central Difference with Runge-Kutta3 on staggered grid
!       (1-1) Two Phase Flow solver
!       (1-2) Single Phase Flow solver
!       (1-3) RK3 for advection and diffusion
!       (1-4) Update Rho
!       (1-5) boundary conditions for u, v, w
!       (1-6) Boudnary conditions for p
!       (1-7) Boudnary conditions for rho and mu
!       (1-8) Initialize solver
!     (2) VSIAM: Multi-moment solver on multi-moment grid
!------------------------------------------
! Author: Zhouteng Ye (yzt9zju@gmail.com)
!------------------------------------------
! Note:
!  1) Some of very frequently used variables are not passed via the
!  interface of the subroutine, but used from ModGlobal.
!  2) The field variables from the MofGlobal will not be used from module,
!     but passed from the interface of the subtuine.
!  3) Other field variables such as divergence, density, viscosity are
!     defined in this module and there is no need to pass those field
!     variables to via the interface of the subroutine.
!     but passed from the interface of the subtuine.
!  4) When calling TwoPhaseFlow
!   If cx, cy, cz appear,
!        VOFAdvection calls MOF
!   Otherwise,
!        VOFAdvection calls VOF
!
! The following variables are not passed from the subroutine interface,
!  but comes from the ModGlobal
!     mpi variables
!     dt: time step
!     nl: number of total grids, with rank 3
!     dl: grid size, with rank 3
!
!  2)  List of solvers:
!    In input.namelist, set hypre_solver:
!       1: SMG
!       2: PFMG
!       3: BicgSTAB
!       4: GMRES
!    List of pre_conditioner:
!       0: No preconditioner
!       1: SMG
!       2: FPMG
! Working pairs are
!     |------------|-----|------|----------|-------|
!     |            | SMG | FPMG | BicgSTAB | GMRES |
!     |------------|-----|------|----------|-------|
!     |   no pre:  | Yes | Yes  | Yes      | Yes   |
!     |   SMG:     |     |      | Yes      | Yes   |
!     |   FPMG:    |     |      | Yes      | Yes   |
!     |------------|-----|------|----------|-------|
!    SMG  seems work best
!    FPMG seems not working properly with multi-phase flow
!===============================================
#include "param.h"
Module ModNavierStokes
  Use mpi
  Use ModGlobal, only : sp
  Use ModGlobal, only : dl, nl, dt, dt0
  Use ModGlobal, only : Phi_bc, P_bc, U_bc, V_bc, W_bc
  Use ModGlobal, only : h5_input
  Use ModGlobal, only : myid, nproc, ierr
  Use ModGlobal, only : left, right, top, bottom, front, back
  Use ModGlobal, only : Field
  Use ModTools
  Use HypreStruct
  Use ModVOF
  Implicit None

  public
  ! private
  ! Public :: TwoPhaseFlow, SinglePhaseFlow, InitNavierStokes, Monitor
  ! Public :: CSF, CSF_VOSET, body_force, bc_uvw
  ! Public :: n_iter, div_max
  ! Module variables will not pass through the subroutine interface
  Real(sp), Allocatable :: Rho(:,:,:), Rhox(:,:,:), Rhoy(:,:,:), Rhoz(:,:,:), PP(:,:,:)
  Real(sp), Allocatable :: mu(:,:,:), muxy(:,:,:), muxz(:,:,:), muyz(:,:,:)
  Real(sp), Allocatable :: U_sx(:,:,:), U_sy(:,:,:), U_sz(:,:,:), U_v(:,:,:)
  Real(sp), Allocatable :: V_sx(:,:,:), V_sy(:,:,:), V_sz(:,:,:), V_v(:,:,:)
  Real(sp), Allocatable :: W_sx(:,:,:), W_sy(:,:,:), W_sz(:,:,:), W_v(:,:,:)
  Real(sp), Allocatable :: ls(:,:,:)
  Real(sp), Allocatable :: Div(:,:,:)
  Integer, Allocatable :: flag(:,:,:)
  Type(Field) :: ls_bc
  Real(sp) :: body_force(3)
  Real(sp) :: rho_l, rho_g
  Real(sp) ::  mu_l, mu_g
  Real(sp) ::  sigma
  Integer :: rk_order
  Real(sp), Allocatable :: rkcoef(:,:)
  Integer :: n_iter
  Real(sp) :: cfl, dt_min
  Real(sp) :: iter_tolerance
  Real(sp) :: p_residual
  Real(sp) :: div_max
  Integer :: step_rank = 0
  Integer :: out_nn = 0
  Real(sp) :: time_out
  integer :: jacobiflag = 0
  Real(sp) :: tec_m = 0.95

  ! PROCEDURE(InterfaceCSF), POINTER :: CSF => CSF_NULL
  ! Interface
  !   Subroutine InterfaceCSF(dudt, dvdt, dwdt, Phi)
  !     import
  !     Implicit None
  !     real(sp), intent(in), dimension(0:,0:,0:) :: Phi
  !     real(sp), intent(inout), dimension(:,:,:) :: dudt, dvdt, dwdt
  !   End Subroutine InterfaceCSF
  ! End Interface

  Type VSIAM
    Type(field) :: BC_Uv
    Type(field) :: BC_Usx
    Type(field) :: BC_Usy
    Type(field) :: BC_Usz
    Type(field) :: BC_Vv
    Type(field) :: BC_Vsx
    Type(field) :: BC_Vsy
    Type(field) :: BC_Vsz
    Type(field) :: BC_Wv
    Type(field) :: BC_Wsx
    Type(field) :: BC_Wsy
    Type(field) :: BC_Wsz
  Contains
    Procedure :: Initialize => InitNavierStokes_VSIAM
  !   Procedure :: SinglePhaseFlow => SinglePhaseflow_VSIAM
    Procedure :: TwoPhaseFlow => TwoPhaseflow_VSIAM
    Procedure :: AdvDiff => AdvDiffVSIAM
    Procedure :: UpdtRhoMu => UpdtRhoMu_VSIAM
    Procedure :: BC_U => BC_U_VSIAM
    Procedure :: BC_V => BC_V_VSIAM
    Procedure :: BC_W => BC_W_VSIAM
    Procedure :: BC_P => BC_P_VSIAM
    Procedure :: BC_RhoMu => BC_RhoMu_VSIAM
  End Type VSIAM

  Type CDSRK3
  Contains
    Procedure :: Initialize => InitNavierStokes_CDSRK3 
    Procedure :: SinglePhaseFlow => SinglePhaseflow_CDSRK3
    Procedure :: TwoPhaseFlow => TwoPhaseflow_CDSRK3
    Procedure :: RK3Single => RungeKuttaMomentumCDS
    Procedure :: RK3TwoPhase => RungeKuttaMomentumUPCDS
    Procedure :: UpdtRhoMu => UpdtRhoMu_CDSRK3
    Procedure :: BC_UVW => BC_UVW_CDSRK3
    Procedure :: BC_P => BC_P_CDSRK3
    Procedure :: BC_RhoMu => BC_RhoMu_CDSRK3
  End Type CDSRK3

  ! Type AdamsCN
  ! Contains
  !   Procedure :: SinglePhaseFlow => SinglePhaseflow_CDSRK3
  !   Procedure :: TwoPhaseFlow => TRungeKuttaMomentumUPCDS
  ! End Type AdamsCN

Contains

  !==========================
  ! (1-1) Two phase flow solver
  ! Key steps:
  !  (1) Runge Kutta time advance for Navier-Stokes equations
  !  (2) Free surface using VOF/MOF
  !--------------------------
  ! Inputs & outputs:
  !   U, V, W: velocity conponents
  !   Phi: volume fraction
  !   P: pressure
  !   cx, cy, cz: centroid of the volume fraction (OPTIONAL)
  ! Note: If cx, cy, cz appear,
  !        VOFAdvection calls MOF
  !   Otherwise,
  !        VOFAdvection calls VOF
  !==========================
  Subroutine TwoPhaseFlow_CDSRK3(self,U, V, W, Phi, P, cx, cy, cz)
    Implicit None
    Class(CDSRK3) :: self
    real(sp), Intent(inout), Dimension(0:,0:,0:) :: u, v, w, p, phi
    real(sp), Intent(inout), Dimension(0:,0:,0:), optional :: cx, cy, cz
    real(sp), dimension(nl(1),nl(2),nl(3))    :: dudtrko,dvdtrko,dwdtrko
    real(sp), dimension(0:nl(1)+1,0:nl(2)+1,0:nl(3)+1)    :: up,vp,wp
    Integer :: irk
    External Divergence, DuDtPartialP, Adjustdt

    ! Runge-Kutta time integration
    Do irk=1,rk_order
      ! Runge-Kutta step for momentum
      call self%RK3TwoPhase(rkcoef(:,irk), u, v, w, p, dudtrko,dvdtrko,dwdtrko, up, vp, wp,phi)

      ! Boundary condition for intermediate velocity
      Call self%BC_UVW(UP, VP, WP)
      ! Calculate divergence
      Call Divergence(up, vp, wp, nl, dl, div)

      ! Solve pressure Poisson equation
      Call Hypre_Poisson(PP, Rhox, Rhoy, Rhoz, Div, flag, n_iter, p_residual)
      ! Call Jacobi(PP)
      Call self%BC_P(PP)

      ! Preject the velocity to the divergence-free field
      Call DuDtPartialP(up, vp, wp, Rhox, Rhoy, Rhoz, pp, dt, nl, dl, u, v, w)
      ! Boundary condition for divergence-free velocity field
      Call self%BC_UVW(U, V, W)
      ! Update p
      p(1:nl(1),1:nl(2),1:nl(3)) = p(1:nl(1),1:nl(2),1:nl(3)) + pp(1:nl(1),1:nl(2),1:nl(3))
      p = p - sum(P(1:nl(1),1:nl(2),1:nl(3))) / (nl(1)) / (nl(2)) / (nl(3))
      Call self%BC_P(P)
    End Do

    ! VOF/MOF advection
    step_rank = step_rank+1
    If (present(cx)) Then ! MOF
      Call VOFAdvection(Phi, u, v, w, nl, dl, dt, mod(step_rank,6), cx, cy, cz)
    Else !VOF
      Call VOFAdvection(Phi, u, v, w, nl, dl, dt, mod(step_rank,6))
    End If

    ! Update Rho, Mu, and adjust time step
    Call self%UpdtRhoMu(Phi)
    Call Adjustdt(U, V, W, nl, dl, dt, dt0, dt_min, cfl)

  End Subroutine TwoPhaseFlow_CDSRK3

  !==========================
  ! (1-2) Single phase flow solver
  ! Key steps:
  !  (1) Runge Kutta time advance for Navier-Stokes equations
  !  (2) Free surface using VOF/MOF
  !--------------------------
  ! Inputs & outputs:
  !   U, V, W: velocity conponents
  !   Phi: volume fraction
  !   P: pressure
  !==========================
  Subroutine SinglePhaseFlow_CDSRK3(self,U, V, W, P)
    Implicit None
    Class(CDSRK3) :: self
    real(sp), Intent(inout), Dimension(0:,0:,0:) :: u, v, w, p
    real(sp), dimension(nl(1),nl(2),nl(3))    :: dudtrko,dvdtrko,dwdtrko
    real(sp), dimension(0:nl(1)+1,0:nl(2)+1,0:nl(3)+1)    :: up,vp,wp
    Integer :: irk
    External Divergence, DuDtPartialP, Adjustdt

    ! Runge-Kutta time integration
    Do irk=1,rk_order
      ! Runge-Kutta step for momentum
      call self%RK3Single(rkcoef(:,irk), u, v, w, p, dudtrko,dvdtrko,dwdtrko, up, vp, wp)
      ! Boundary condition for intermediate velocity
      Call self%BC_UVW(UP, VP, WP)
      ! Calculate divergence
      Call Divergence(up, vp,wp, nl, dl, div)
      ! Solve pressure Poisson equation
      Call Hypre_Poisson(PP, Rhox, Rhoy, Rhoz, Div, flag, n_iter, p_residual)
        ! Call Jacobi(PP)
      Call self%BC_P(PP)
      ! Preject the velocity to the divergence-free field
      Call DuDtPartialP(up, vp, wp, Rhox, Rhoy, Rhoz, PP, dt, nl, dl, u, v, w)
      ! Boundary condition for divergence-free velocity field
      Call self%BC_UVW(U, V, W)
      ! Update p
      p(1:nl(1),1:nl(2),1:nl(3)) = p(1:nl(1),1:nl(2),1:nl(3)) + pp(1:nl(1),1:nl(2),1:nl(3))
      Call self%BC_P(P)
    End Do

    Call Adjustdt(U, V, W, nl, dl, dt, dt0, dt_min, cfl)

  End Subroutine SinglePhaseFlow_CDSRK3

  !==========================
  ! (1-3) Runge Kutta time integration 
  ! Runge Kutta step for
  ! U^star = -1 / \rho * \nabla p^k-1/2 + \alpha * AD^k-1/2 + \beta AD^k-1/2
  ! For DNS
  !--------------------------
  ! Inputs:
  !   rkpar: runge-kutta parameters
  !   U, V, W: velocity conponents for time step k-1/2
  !   P: pressure for time step k-1/2
  !   dudtrko, dvdtrko, dwdtrko: du/dt, dv/dt, dw/dt for time step k-1/2
  ! Outpurts
  !   UP, VP, WP: u*star, v_star, w_star
  !   dudtrko, dvdtrko, dwdtrko: du/dt, dv/dt, dw/dt for time step k+1/2
  !==========================
  Subroutine RungeKuttaMomentumCDS(self, rkpar, u, v, w, p, dudtrko, dvdtrko, dwdtrko, up, vp, wp, Phi)
    Implicit None
    Class(CDSRK3) :: self
    real(sp), Intent(in), Dimension(0:,0:,0:) :: u,v,w, p
    real(sp), Intent(inout), Dimension(:,:,:) :: dudtrko,dvdtrko,dwdtrko
    real(sp), Intent(in), Dimension(2) :: rkpar
    real(sp), Intent(out), Dimension(0:,0:,0:)  :: up,vp,wp
    real(sp), Intent(in), Dimension(0:,0:,0:),optional  :: phi
    real(sp), Dimension(nl(1),nl(2),nl(3)) :: dudtrk,dvdtrk,dwdtrk
    real(sp) :: factor1,factor2,factor12
    Integer :: i, j, k
    External DuDtPartialP, AdvDiffU_CDS, AdvDiffV_CDS, AdvDiffW_CDS, BodyForce

    factor1 = rkpar(1)*dt
    factor2 = rkpar(2)*dt
    factor12 = factor1 + factor2

    ! Calculate -1 / \tho * \nabla p^k-1/2
    Call DuDtPartialP(u, v, w, rhox, rhoy, rhoz, P, factor12, nl, dl, up, vp, wp)

    ! Calculate AD^k+1/2
    dudtrk = 0.0_sp
    dvdtrk = 0.0_sp
    dwdtrk = 0.0_sp
    Call AdvDiffU_CDS(u, v, w, dl, nl, mu, muxy, muxz, rhox, dudtrk)
    Call AdvDiffV_CDS(u, v, w, dl, nl, mu, muxy, muyz, rhoy, dvdtrk)
    Call AdvDiffW_CDS(u, v, w, dl, nl, mu, muxz, muyz, rhoz, dwdtrk)
    ! Calculate Body force k+1/2 (Included in AD^k+1/2 term)
    Call BodyForce(dudtrk, dvdtrk, dwdtrk, nl, body_force)

    ! if ( present(phi) ) Call CSF(dudtrk, dvdtrk, dwdtrk, Phi)

    ! UP - U^start
    ! dUdtro = U^start - U^k-1/2
    do k=1,nl(3)
      do j=1,nl(2)
        do i=1,nl(1)
          up(i,j,k) = up(i,j,k) + factor1*dudtrk(i,j,k) + factor2*dudtrko(i,j,k)
          vp(i,j,k) = vp(i,j,k) + factor1*dvdtrk(i,j,k) + factor2*dvdtrko(i,j,k)
          wp(i,j,k) = wp(i,j,k) + factor1*dwdtrk(i,j,k) + factor2*dwdtrko(i,j,k)
          dudtrko(i,j,k) = dudtrk(i,j,k)
          dvdtrko(i,j,k) = dvdtrk(i,j,k)
          dwdtrko(i,j,k) = dwdtrk(i,j,k)
        enddo
      enddo
    enddo

  End Subroutine RungeKuttaMomentumCDS

  !==========================
  ! (1-3-1) Runge Kutta time integration 
  ! Runge Kutta step for
  ! U^star = -1 / \rho * \nabla p^k-1/2 + \alpha * AD^k-1/2 + \beta AD^k-1/2
  ! Upwind for advection, central difference for diffusion
  ! For multiphase flow
  !--------------------------
  ! Inputs:
  !   rkpar: runge-kutta parameters
  !   U, V, W: velocity conponents for time step k-1/2
  !   P: pressure for time step k-1/2
  !   dudtrko, dvdtrko, dwdtrko: du/dt, dv/dt, dw/dt for time step k-1/2
  ! Outpurts
  !   UP, VP, WP: u*star, v_star, w_star
  !   dudtrko, dvdtrko, dwdtrko: du/dt, dv/dt, dw/dt for time step k+1/2
  !==========================
  Subroutine RungeKuttaMomentumUPCDS(self, rkpar, u, v, w, p, dudtrko, dvdtrko, dwdtrko, up, vp, wp, Phi)
    Implicit None
    Class(CDSRK3) :: self
    real(sp), Intent(in), Dimension(0:,0:,0:) :: u,v,w, p
    real(sp), Intent(inout), Dimension(:,:,:) :: dudtrko,dvdtrko,dwdtrko
    real(sp), Intent(in), Dimension(2) :: rkpar
    real(sp), Intent(out), Dimension(0:,0:,0:)  :: up,vp,wp
    real(sp), Intent(in), Dimension(0:,0:,0:),optional  :: phi
    real(sp), Dimension(nl(1),nl(2),nl(3)) :: dudtrk,dvdtrk,dwdtrk
    real(sp) :: factor1,factor2,factor12
    Integer :: i, j, k
    External DuDtPartialP, AdvDiffU_CDS, AdvDiffV_CDS, AdvDiffW_CDS, BodyForce

    factor1 = rkpar(1)*dt
    factor2 = rkpar(2)*dt
    factor12 = factor1 + factor2

    ! Calculate -1 / \tho * \nabla p^k-1/2
    Call DuDtPartialP(u, v, w, rhox, rhoy, rhoz, P, factor12, nl, dl, up, vp, wp)

    ! Calculate AD^k+1/2
    dudtrk = 0.0_sp
    dvdtrk = 0.0_sp
    dwdtrk = 0.0_sp
    Call AdvDiffU_UPCDS(u, v, w, dl, nl, mu, muxy, muxz, rhox, dudtrk)
    Call AdvDiffV_UPCDS(u, v, w, dl, nl, mu, muxy, muyz, rhoy, dvdtrk)
    Call AdvDiffW_UPCDS(u, v, w, dl, nl, mu, muxz, muyz, rhoz, dwdtrk)
    ! Calculate Body force k+1/2 (Included in AD^k+1/2 term)
    Call BodyForce(dudtrk, dvdtrk, dwdtrk, nl, body_force)

    ! if ( present(phi) ) Call CSF(dudtrk, dvdtrk, dwdtrk, Phi)

    ! UP - U^start
    ! dUdtro = U^start - U^k-1/2
    do k=1,nl(3)
      do j=1,nl(2)
        do i=1,nl(1)
          up(i,j,k) = up(i,j,k) + factor1*dudtrk(i,j,k) + factor2*dudtrko(i,j,k)
          vp(i,j,k) = vp(i,j,k) + factor1*dvdtrk(i,j,k) + factor2*dvdtrko(i,j,k)
          wp(i,j,k) = wp(i,j,k) + factor1*dwdtrk(i,j,k) + factor2*dwdtrko(i,j,k)
          dudtrko(i,j,k) = dudtrk(i,j,k)
          dvdtrko(i,j,k) = dvdtrk(i,j,k)
          dwdtrko(i,j,k) = dwdtrk(i,j,k)
        enddo
      enddo
    enddo

  End Subroutine RungeKuttaMomentumUPCDS



  !==========================
  ! (1-4) Update Rho and Mu
  !    key steps:
  !    (1) Calculate the density and viscosity using volume fraction
  !    (2) Interpolate to get 1/Rho at cell face center
  !    (3) Interpolate to get mu at cell edge
  !--------------------------
  ! Inputs:
  !    phi: volume fraction
  ! Note:
  !    phi and mu related variables are updated in this subroutine
  !    and can be used in this module.
  !==========================
  Subroutine UpdtRhoMu_CDSRK3(self,Phi)
    Implicit None
    Class(CDSRK3) :: self
    real(sp), intent(out), dimension(0:,0:,0:) :: Phi
    Integer :: i,j,k

    Call Phi_bc%SetBCS(Phi)

    Do k=1,nl(3)
      Do j=1,nl(2)
        Do i=1,nl(1)
          Rho(i,j,k) = Phi(i,j,k) * rho_l + ( 1.0_sp - Phi(i,j,k) ) * rho_g
          Mu(i,j,k)  = Phi(i,j,k) *  mu_l + ( 1.0_sp - Phi(i,j,k) ) *  mu_g
        enddo
      enddo
    enddo

    Call self%BC_RHOMU(Rho, Mu)

    Do k=0,nl(3)
      Do j=0,nl(2)
        Do i=0,nl(1)
          Rhox(i,j,k) = 2.0_sp / (Rho(i,j,k) + Rho(i+1,j,k))
          Rhoy(i,j,k) = 2.0_sp / (Rho(i,j,k) + Rho(i,j+1,k))
          Rhoz(i,j,k) = 2.0_sp / (Rho(i,j,k) + Rho(i,j,k+1))
          Muxy(i,j,k) = 0.25_sp * (Mu(i,j,k) + Mu(i+1,  j,k) + Mu(i,j+1,  k) + Mu(i+1,j+1,  k))
          Muyz(i,j,k) = 0.25_sp * (Mu(i,j,k) + Mu(  i,j+1,k) + Mu(i,  j,k+1) + Mu(  i,j+1,k+1))
          Muxz(i,j,k) = 0.25_sp * (Mu(i,j,k) + Mu(i+1,  j,k) + Mu(i,  j,k+1) + Mu(i+1,  j,k+1))
        EndDo
      EndDo
    EndDo

  End Subroutine UpdtRhoMu_CDSRK3

  !==========================
  ! (1-5) Boudnary conditions for u, v, w
  !==========================
  Subroutine BC_UVW_CDSRK3(self,U, V, W)
    Implicit None
    Class(CDSRK3) :: self
    real(sp), intent(inout) , dimension(0:,0:,0:) :: U, V, W
    Call U_bc%SetBCS(U)
    Call V_bc%SetBCS(V)
    Call W_bc%SetBCS(W)
  End Subroutine BC_UVW_CDSRK3

  !==========================
  ! (1-6) Boudnary conditions for p
  !==========================
  Subroutine BC_P_CDSRK3(self,P)
    Implicit None
    Class(CDSRK3) :: self
    real(sp), intent(inout) , dimension(0:,0:,0:) :: P
    Call Phi_bc%SetBCS(P)
  End Subroutine BC_P_CDSRK3

  !==========================
  ! (1-7) Boudnary conditions for rho and mu
  !==========================
  Subroutine BC_RHOMU_CDSRK3(self,Rho, Mu)
    Implicit None
    Class(CDSRK3) :: self
    real(sp), intent(inout) , dimension(0:,0:,0:) :: Rho, Mu
    Call Phi_bc%SetBCS(Rho)
    Call Phi_bc%SetBCS(Mu)
  End Subroutine BC_RHOMU_CDSRK3

  !========================
  ! (1-8) Initialize the Navier-Stokes solver after common initialization
  ! Note:
  !   The bc types and values are initialized in basics.f90 for only VOF
  !   Here will overwrite some of the value.
  !==========================
  Subroutine InitNavierStokes_CDSRK3(self,u, v, w, phi, p)
    Implicit None
    Class(CDSRK3) :: self
    real(sp), intent(inout) , dimension(0:,0:,0:) :: p,u,v,w,phi
    Integer :: iter_max
    Character(80) :: file_name
    Character(80) :: input_name
    Integer :: bc_left, bc_right, bc_bottom, bc_top, bc_front, bc_back
    Integer :: bc_pair(2)
    Integer :: hypre_solver, hypre_PreConditioner
    Integer :: i

    !  variables for physics
    namelist /ns_physics/ rho_l, rho_g, mu_l, mu_g, sigma, body_force, cfl, dt_min
    namelist /ns_solver/ iter_tolerance, iter_max, rk_order, hypre_solver, Hypre_PreConditioner
    !  variables for boundary conditions
    namelist /ns_bc/ &
        bc_left, bc_right, &
        bc_back, bc_front, &
        bc_top, bc_bottom

    Call getarg(1,input_name)
    if (INPUT_NAME .eq. '') Then
      file_name = 'input.namelist'
      h5_input%filename = "input.h5"
    Else
      file_name = trim(input_name)//'.namelist'
      h5_input%filename = trim(input_name)//'.h5'
    endif

    Do i = 1, nproc
      if (myid .eq. i-1) then
        Open(10, file=file_name)
        Read(10, nml = ns_physics)
        Close(10)
        Open(10, file=file_name)
        Read(10, nml = ns_solver)
        Close(10)
        Open(10, file=file_name)
        Read(10, nml = ns_bc)
        Close(10)
      Endif
      Call MPI_barrier(MPI_COMM_WORLD, ierr)
    End Do

    ! Allocate variables
    Allocate( Rho(0:nl(1)+1, 0:nl(2)+1, 0:nl(3)+1))
    Allocate(Rhox(0:nl(1)+1, 0:nl(2)+1, 0:nl(3)+1))
    Allocate(Rhoy(0:nl(1)+1, 0:nl(2)+1, 0:nl(3)+1))
    Allocate(Rhoz(0:nl(1)+1, 0:nl(2)+1, 0:nl(3)+1))
    Allocate(  Mu(0:nl(1)+1, 0:nl(2)+1, 0:nl(3)+1))
    Allocate(Muxy(0:nl(1)+1, 0:nl(2)+1, 0:nl(3)+1))
    Allocate(Muxz(0:nl(1)+1, 0:nl(2)+1, 0:nl(3)+1))
    Allocate(Muyz(0:nl(1)+1, 0:nl(2)+1, 0:nl(3)+1))
    Allocate(  LS(0:nl(1)+1, 0:nl(2)+1, 0:nl(3)+1))
    Allocate(Flag(0:nl(1)+1, 0:nl(2)+1, 0:nl(3)+1))
    Allocate(PP(0:nl(1)+1, 0:nl(2)+1, 0:nl(3)+1))
    Allocate( Div(nl(1), nl(2), nl(3)))
    Allocate(rkcoef(2,rk_order))

    Rho  = 0.0_sp
    Rhox = 0.0_sp
    Rhoy = 0.0_sp
    Rhoz = 0.0_sp
    Mu   = 0.0_sp
    Div  = 0.0_sp
    LS   = 0.0_sp
    Flag = 1
    rkcoef   = 0.0_sp

    ls_bc%name = 'ls'
    ls_bc%lohi(:,1) = 0
    ls_bc%lohi(:,2) = nl(:)+1
    phi_bc%bound_type(:,:) = 3
    phi_bc%bound_value(:,:) = 0

    ! Boudnary conditions:
    ! 1. No-slip wall
    ! 2. Slip wall/ Symmetric

    ! Left and right
    bc_pair(1) = bc_left
    bc_pair(2) = bc_right
    Do i = 1, 2
      Select case (bc_pair(i))
      Case(1)
        u_bc%bound_type(1,i)  = 1
        u_bc%bound_value(1,i) = 0.0_sp
        v_bc%bound_type(1,i)  = 1
        v_bc%bound_value(1,i) = 0.0_sp
        w_bc%bound_type(1,i)  = 1
        w_bc%bound_value(1,i) = 0.0_sp
        p_bc%bound_type(1,i)  = 2
        p_bc%bound_value(1,i) = 0.0_sp
        phi_bc%bound_type(1,i)  = 2
        phi_bc%bound_value(1,i) = 0.0_sp
      Case(2)
        u_bc%bound_type(1,i)  = 1
        u_bc%bound_value(1,i) = 0.0_sp
        v_bc%bound_type(1,i)  = 2
        v_bc%bound_value(1,i) = 0.0_sp
        w_bc%bound_type(1,i)  = 2
        w_bc%bound_value(1,i) = 0.0_sp
        p_bc%bound_type(1,i)  = 2
        p_bc%bound_value(1,i) = 0.0_sp
        phi_bc%bound_type(1,i)  = 2
        phi_bc%bound_value(1,i) = 0.0_sp
      Case Default
        if (myid .eq.0) then
          print *, "======Fatal Error=============================="
          print *, "Wrong bc value, please check the namelist file"
          print *, "==============================================="
        end if
        Call MPI_Finalize(ierr)
        stop
      End Select
    End Do

    ! Front and back
    bc_pair(1) = bc_back
    bc_pair(2) = bc_front
    Do i = 1, 2
      Select case (bc_pair(i))
      Case(1)
        u_bc%bound_type(2,i)  = 1
        u_bc%bound_value(2,i) = 0.0_sp
        v_bc%bound_type(2,i)  = 1
        v_bc%bound_value(2,i) = 0.0_sp
        w_bc%bound_type(2,i)  = 1
        w_bc%bound_value(2,i) = 0.0_sp
        p_bc%bound_type(2,i)  = 2
        p_bc%bound_value(2,i) = 0.0_sp
        phi_bc%bound_type(2,i)  = 2
        phi_bc%bound_value(2,i) = 0.0_sp
      Case(2)
        u_bc%bound_type(2,i)  = 2
        u_bc%bound_value(2,i) = 0.0_sp
        v_bc%bound_type(2,i)  = 1
        v_bc%bound_value(2,i) = 0.0_sp
        w_bc%bound_type(2,i)  = 2
        w_bc%bound_value(2,i) = 0.0_sp
        p_bc%bound_type(2,i)  = 2
        p_bc%bound_value(2,i) = 0.0_sp
        phi_bc%bound_type(2,i)  = 2
        phi_bc%bound_value(2,i) = 0.0_sp
      End Select
    End Do

    ! Top and bottom
    bc_pair(1) = bc_bottom
    bc_pair(2) = bc_top
    Do i = 1, 2
      Select case (bc_pair(i))
      Case(1)
        u_bc%bound_type(3,i)  = 1
        u_bc%bound_value(3,i) = 0.0_sp
        v_bc%bound_type(3,i)  = 1
        v_bc%bound_value(3,i) = 0.0_sp
        w_bc%bound_type(3,i)  = 1
        w_bc%bound_value(3,i) = 0.0_sp
        p_bc%bound_type(3,i)  = 2
        p_bc%bound_value(3,i) = 0.0_sp
        phi_bc%bound_type(3,i)  = 2
        phi_bc%bound_value(3,i) = 0.0_sp
      Case(2)
        u_bc%bound_type(3,i)  = 2
        u_bc%bound_value(3,i) = 0.0_sp
        v_bc%bound_type(3,i)  = 2
        v_bc%bound_value(3,i) = 0.0_sp
        w_bc%bound_type(3,i)  = 1
        w_bc%bound_value(3,i) = 0.0_sp
        p_bc%bound_type(3,i)  = 2
        p_bc%bound_value(3,i) = 0.0_sp
        phi_bc%bound_type(3,i)  = 2
        phi_bc%bound_value(3,i) = 0.0_sp
      End Select
    End Do

    ! Assign RungeKutta
    If (rk_order .eq. 1) Then
      rkcoef(1,1) =  1.0_sp
      rkcoef(2,1) =  0.0_sp
    ElseIf (rk_order .eq. 2) Then
      rkcoef(1,1) =  0.5_sp
      rkcoef(2,1) =  0.0_sp
      rkcoef(1,2) =  1.0_sp
      rkcoef(2,2) =  -0.5_sp
    ElseIf (rk_order .eq. 3) Then
      rkcoef(1,1) =  32.0_sp / 60.0_sp
      rkcoef(2,1) =  0.0_sp  / 60.0_sp
      rkcoef(1,2) =  25.0_sp / 60.0_sp
      rkcoef(2,2) = -17.0_sp / 60.0_sp
      rkcoef(1,3) =  45.0_sp / 60.0_sp
      rkcoef(2,3) = -25.0_sp / 60.0_sp
    Else
      If ( myid .eq. 0 ) then
        print *, "======Fatal Error=============================="
        print *, "Wong rk order, should be 1 or 2 or 3"
        print *, "==============================================="
      End If
      Call MPI_Finalize(ierr)
      stop
    End If

    ! Initialize hypre
    Call Hypre_Initialize(iter_tolerance, iter_max, hypre_solver, Hypre_PreConditioner)

    ! Initialize the Rho and Mu field
    Call self%UpdtRhoMu(Phi)

    ! Initialize boundary condition
    Call self%BC_UVW(U,V,W)
    Call self%BC_P(P)

  End Subroutine InitNavierStokes_CDSRK3

  Subroutine TwoPhaseFlow_VSIAM(self,U, V, W, Phi, P, cx, cy, cz)
    Implicit None
    Class(VSIAM) :: self
    real(sp), Intent(inout), Dimension(0:,0:,0:) :: u, v, w, p, phi
    real(sp), Intent(inout), Dimension(0:,0:,0:), optional :: cx, cy, cz
    real(sp), dimension(0:nl(1)+1,0:nl(2)+1,0:nl(3)+1)    :: Un_sx, Vn_sy, Wn_sz
    real(sp), dimension(0:nl(1)+1,0:nl(2)+1,0:nl(3)+1)    :: Un_v, Vn_v, Wn_v
    External Divergence, DuDtPartialP, Adjustdt

    U_sx = U
    V_sy = V
    W_sz = W

    ! VSIAM Advection and Diffusion 
    Call self%AdvDiff

    Call Divergence(U_sx, V_sy, W_sz, nl, dl, div)

    ! Solve pressure Poisson equation
    Call Hypre_Poisson(PP, Rhox, Rhoy, Rhoz, Div, flag, n_iter, p_residual)
    ! Call Jacobi(PP)
    Call self%BC_P(PP)

    Un_sx = U_sx
    Vn_sy = V_sy
    Wn_sz = W_sz
    Un_v  = U_v
    Vn_v  = V_v
    Wn_v  = W_v

    ! Preject the velocity to the divergence-free field
    Call DuDtPartialP(Un_sx, Vn_sy, Wn_sz, Rhox, Rhoy, Rhoz, pp, dt, nl, dl, U_sx, V_sy, W_sz)
    ! Boundary condition for divergence-free velocity field
    ! Update p
    P = PP
    ! p(1:nl(1),1:nl(2),1:nl(3)) = p(1:nl(1),1:nl(2),1:nl(3)) + pp(1:nl(1),1:nl(2),1:nl(3))
    ! p = p - sum(P(1:nl(1),1:nl(2),1:nl(3))) / (nl(1)) / (nl(2)) / (nl(3))
    ! Call self%BC_P(P)

    ! Update VIA by TEC
    Call Tec_VIA(U_v, U_sx, Un_sx, nl, 1)
    Call Tec_VIA(V_v, V_sy, Vn_sy, nl, 2)
    Call Tec_VIA(W_v, W_sz, Wn_sz, nl, 3)

    ! Update SIA_X, SIA_Y, SIA_Z by TEC
    Call Tec_SIA(U_sy, U_v, Un_v, nl, tec_m, 2)
    Call Tec_SIA(U_sz, U_v, Un_v, nl, tec_m, 3)
    Call Tec_SIA(V_sx, V_v, Vn_v, nl, tec_m, 1)
    Call Tec_SIA(V_sz, V_v, Vn_v, nl, tec_m, 3)
    Call Tec_SIA(W_sx, W_v, Wn_v, nl, tec_m, 1)
    Call Tec_SIA(W_sy, W_v, Wn_v, nl, tec_m, 2)

    ! Boundary Conditions
    Call self%BC_V(V_sx, V_sy, V_sz, V_v)
    Call self%BC_W(W_sx, W_sy, W_sz, W_v)
    Call self%BC_U(U_sx, U_sy, U_sz, U_v)

    ! VOF/MOF advection
    step_rank = step_rank+1
    If (present(cx)) Then ! MOF
      Call VOFAdvection(Phi, u_sx, v_sy, w_sz, nl, dl, dt, mod(step_rank,6), cx, cy, cz)
    Else !VOF
      Call VOFAdvection(Phi, u_sx, v_sy, w_sz, nl, dl, dt, mod(step_rank,6))
    End If

    ! Update Rho, Mu, and adjust time step
    Call self%UpdtRhoMu(Phi)
    Call Adjustdt(U_v, V_v, W_v, nl, dl, dt, dt0, dt_min, cfl)

    U = U_sx
    V = V_sy
    W = W_sz

  End Subroutine TwoPhaseFlow_VSIAM

  Subroutine AdvDiffVSIAM(self)
    Implicit None
    Class(VSIAM) :: self
    real(sp), dimension(0:nl(1)+1,0:nl(2)+1,0:nl(3)+1)    :: Un_sx, Vn_sy, Wn_sz
    real(sp), dimension(0:nl(1)+1,0:nl(2)+1,0:nl(3)+1)    :: Un_v, Vn_v, Wn_v
    real(sp), dimension(nl(1),nl(2),nl(3))    :: RHS_U, RHS_V, RHS_W
    Integer :: i,j,k
    External CIPCSLR_1D, Tec_SIA, Central_Difference_RHS

    Un_sx = U_sx
    Vn_sy = V_sy
    Wn_sz = W_sz

    ! Calculation of U momentum
    Un_v = U_v
    Call CIPCSLR_1D(U_v, U_sx, Un_sx, nl, dl, dt, 1)
    Call Tec_SIA(U_sy, U_v, Un_v, nl, tec_m, 2)
    Call Tec_SIA(U_sz, U_v, Un_v, nl, tec_m, 3)
    Call self%BC_U(U_sx, U_sy, U_sz, U_v)
    Un_v = U_v
    Call CIPCSLR_1D(U_v, U_sy, Vn_sy, nl, dl, dt, 2)
    Call Tec_SIA(U_sx, U_v, Un_v, nl, tec_m, 1)
    Call Tec_SIA(U_sz, U_v, Un_v, nl, tec_m, 3)
    Call self%BC_U(U_sx, U_sy, U_sz, U_v)
    Un_v = U_v
    Call CIPCSLR_1D(U_v, U_sz, Wn_sz, nl, dl, dt, 3)
    Call Tec_SIA(U_sx, U_v, Un_v, nl, tec_m, 1)
    Call Tec_SIA(U_sy, U_v, Un_v, nl, tec_m, 2)
    Call self%BC_U(U_sx, U_sy, U_sz, U_v)

    ! Calculation of U momentum
    Vn_v = V_v
    Call CIPCSLR_1D(V_v, V_sx, Un_sx, nl, dl, dt, 1)
    Call Tec_SIA(V_sy, V_v, Vn_v, nl, tec_m, 2)
    Call Tec_SIA(V_sz, V_v, Vn_v, nl, tec_m, 3)
    Call self%BC_V(V_sx, V_sy, V_sz, V_v)
    Vn_v = V_v
    Call CIPCSLR_1D(V_v, V_sy, Vn_sy, nl, dl, dt, 2)
    Call Tec_SIA(V_sx, V_v, Vn_v, nl, tec_m, 1)
    Call Tec_SIA(V_sz, V_v, Vn_v, nl, tec_m, 3)
    Call self%BC_V(V_sx, V_sy, V_sz, V_v)
    Vn_v = V_v
    Call CIPCSLR_1D(V_v, V_sz, Wn_sz, nl, dl, dt, 3)
    Call Tec_SIA(V_sx, V_v, Vn_v, nl, tec_m, 1)
    Call Tec_SIA(V_sy, V_v, Vn_v, nl, tec_m, 2)
    Call self%BC_V(V_sx, V_sy, V_sz, V_v)

    ! Calculation of U momentum
    Wn_v = W_v
    Call CIPCSLR_1D(W_v, W_sx, Un_sx, nl, dl, dt, 1)
    Call Tec_SIA(W_sy, W_v, Wn_v, nl, tec_m, 2)
    Call Tec_SIA(W_sz, W_v, Wn_v, nl, tec_m, 3)
    Call self%BC_W(W_sx, W_sy, W_sz, W_v)
    Wn_v = W_v
    Call CIPCSLR_1D(W_v, W_sy, Vn_sy, nl, dl, dt, 2)
    Call Tec_SIA(W_sx, W_v, Wn_v, nl, tec_m, 1)
    Call Tec_SIA(W_sz, W_v, Wn_v, nl, tec_m, 3)
    Call self%BC_W(W_sx, W_sy, W_sz, W_v)
    Wn_v = W_v
    Call CIPCSLR_1D(W_v, W_sz, Wn_sz, nl, dl, dt, 3)
    Call Tec_SIA(W_sx, W_v, Wn_v, nl, tec_m, 1)
    Call Tec_SIA(W_sy, W_v, Wn_v, nl, tec_m, 2)
    Call self%BC_W(W_sx, W_sy, W_sz, W_v)

    Un_v = U_v
    Vn_v = V_v
    Wn_v = W_v 

    ! Calculation of the right-hand-side of diffusion term 
    Call Central_Difference_RHS(RHS_U,RHS_V,RHS_W,U_v, V_v, W_v, mu, rho, nl, dl)

    ! Update VIA
    Do k=1,nl(3)
      Do j=1,nl(2)
        Do i=1,nl(1)
      U_v(i,j,k) = Un_v(i,j,k) + RHS_U(i,j,k) * dt
      V_v(i,j,k) = Vn_v(i,j,k) + RHS_V(i,j,k) * dt
      W_v(i,j,k) = Wn_v(i,j,k) + RHS_W(i,j,k) * dt
        EndDo
        EndDo
    EndDo

    ! Update SIA_X, SIA_Y, SIA_Z by TEC
    Call Tec_SIA(U_sx, U_v, Un_v, nl, tec_m, 1)
    Call Tec_SIA(U_sy, U_v, Un_v, nl, tec_m, 2)
    Call Tec_SIA(U_sz, U_v, Un_v, nl, tec_m, 3)
    Call Tec_SIA(V_sx, V_v, Vn_v, nl, tec_m, 1)
    Call Tec_SIA(V_sy, V_v, Vn_v, nl, tec_m, 2)
    Call Tec_SIA(V_sz, V_v, Vn_v, nl, tec_m, 3)
    Call Tec_SIA(W_sx, W_v, Wn_v, nl, tec_m, 1)
    Call Tec_SIA(W_sy, W_v, Wn_v, nl, tec_m, 2)
    Call Tec_SIA(W_sz, W_v, Wn_v, nl, tec_m, 3)

    ! Boundary Conditions
    Call self%BC_V(V_sx, V_sy, V_sz, V_v)
    Call self%BC_W(W_sx, W_sy, W_sz, W_v)
    Call self%BC_U(U_sx, U_sy, U_sz, U_v)

  End Subroutine AdvDiffVSIAM

  !==========================
  ! (1-4) Update Rho and Mu
  !    key steps:
  !    (1) Calculate the density and viscosity using volume fraction
  !    (2) Interpolate to get 1/Rho at cell face center
  !    (3) Interpolate to get mu at cell edge
  !--------------------------
  ! Inputs:
  !    phi: volume fraction
  ! Note:
  !    phi and mu related variables are updated in this subroutine
  !    and can be used in this module.
  !==========================
  Subroutine UpdtRhoMu_VSIAM(self,Phi)
    Implicit None
    Class(VSIAM) :: self
    real(sp), intent(out), dimension(0:,0:,0:) :: Phi
    Integer :: i,j,k

    Call Phi_bc%SetBCS(Phi)

    Do k=1,nl(3)
      Do j=1,nl(2)
        Do i=1,nl(1)
          Rho(i,j,k) = Phi(i,j,k) * rho_l + ( 1.0_sp - Phi(i,j,k) ) * rho_g
          Mu(i,j,k)  = Phi(i,j,k) *  mu_l + ( 1.0_sp - Phi(i,j,k) ) *  mu_g
        enddo
      enddo
    enddo

    ! Call self%BC_RHOMU(Rho, Mu)

    Do k=0,nl(3)
      Do j=0,nl(2)
        Do i=0,nl(1)
          Rhox(i,j,k) = 2.0_sp / (Rho(i,j,k) + Rho(i+1,j,k))
          Rhoy(i,j,k) = 2.0_sp / (Rho(i,j,k) + Rho(i,j+1,k))
          Rhoz(i,j,k) = 2.0_sp / (Rho(i,j,k) + Rho(i,j,k+1))
        EndDo
      EndDo
    EndDo

  End Subroutine UpdtRhoMu_VSIAM

  !==========================
  ! (1-5) Boudnary conditions for u
  !==========================
  Subroutine BC_U_VSIAM(self, Usx, Usy, Usz, Uv)
    Implicit None
    Class(VSIAM) :: self
    real(sp), intent(inout) , dimension(0:,0:,0:) :: Uv, Usx, Usy, Usz
    Call self%BC_Usx%SetBCS(Usx)
    Call self%BC_Usy%SetBCS(Usy)
    Call self%BC_Usz%SetBCS(Usz)
    Call self%BC_Uv%SetBCS(Uv)
  End Subroutine BC_U_VSIAM

  !==========================
  ! (1-5) Boudnary conditions for v
  !==========================
  Subroutine BC_V_VSIAM(self, Vsx, Vsy, Vsz, Vv)
    Implicit None
    Class(VSIAM) :: self
    real(sp), intent(inout) , dimension(0:,0:,0:) :: Vsx, Vsy, Vsz, Vv
    Call self%BC_Vsx%SetBCS(Vsx)
    Call self%BC_Vsy%SetBCS(Vsy)
    Call self%BC_Vsz%SetBCS(Vsz)
    Call self%BC_Vv%SetBCS(Vv)
  End Subroutine BC_V_VSIAM

  !==========================
  ! (1-5) Boudnary conditions for w
  !==========================
  Subroutine BC_W_VSIAM(self, Wsx, Wsy, Wsz, Wv)
    Implicit None
    Class(VSIAM) :: self
    real(sp), intent(inout) , dimension(0:,0:,0:) :: Wsx, Wsy, Wsz, Wv
    Call self%BC_Wsx%SetBCS(Wsx)
    Call self%BC_Wsy%SetBCS(Wsy)
    Call self%BC_Wsz%SetBCS(Wsz)
    Call self%BC_Wv%SetBCS(Wv)
  End Subroutine BC_W_VSIAM

  !==========================
  ! (1-6) Boudnary conditions for p
  !==========================
  Subroutine BC_P_VSIAM(self,P)
    Implicit None
    Class(VSIAM) :: self
    real(sp), intent(inout) , dimension(0:,0:,0:) :: P
    Call Phi_bc%SetBCS(P)
  End Subroutine BC_P_VSIAM

  !==========================
  ! (1-7) Boudnary conditions for rho and mu
  !==========================
  Subroutine BC_RHOMU_VSIAM(self,Rho, Mu)
    Implicit None
    Class(VSIAM) :: self
    real(sp), intent(inout) , dimension(0:,0:,0:) :: Rho, Mu
    Call Phi_bc%SetBCS(Rho)
    Call Phi_bc%SetBCS(Mu)
  End Subroutine BC_RHOMU_VSIAM


  Subroutine InitNavierStokes_VSIAM(self,u, v, w, phi, p)
    Implicit None
    Class(VSIAM) :: self
    real(sp), intent(inout) , dimension(0:,0:,0:) :: p,u,v,w,phi
    Integer :: iter_max
    Character(80) :: file_name
    Character(80) :: input_name
    Integer :: bc_left, bc_right, bc_bottom, bc_top, bc_front, bc_back
    Integer :: bc_pair(2)
    Integer :: hypre_solver, hypre_PreConditioner
    Integer :: i

    !  variables for physics
    namelist /ns_physics/ rho_l, rho_g, mu_l, mu_g, sigma, body_force, cfl, dt_min
    namelist /ns_solver/ iter_tolerance, iter_max, rk_order, hypre_solver, Hypre_PreConditioner
    !  variables for boundary conditions
    namelist /ns_bc/ &
        bc_left, bc_right, &
        bc_back, bc_front, &
        bc_top, bc_bottom

    Call getarg(1,input_name)
    if (INPUT_NAME .eq. '') Then
      file_name = 'input.namelist'
      h5_input%filename = "input.h5"
    Else
      file_name = trim(input_name)//'.namelist'
      h5_input%filename = trim(input_name)//'.h5'
    endif

    Do i = 1, nproc
      if (myid .eq. i-1) then
        Open(10, file=file_name)
        Read(10, nml = ns_physics)
        Close(10)
        Open(10, file=file_name)
        Read(10, nml = ns_solver)
        Close(10)
        Open(10, file=file_name)
        Read(10, nml = ns_bc)
        Close(10)
      Endif
      Call MPI_barrier(MPI_COMM_WORLD, ierr)
    End Do

    ! Allocate variables
    Allocate( Rho(0:nl(1)+1, 0:nl(2)+1, 0:nl(3)+1))
    Allocate(Rhox(0:nl(1)+1, 0:nl(2)+1, 0:nl(3)+1))
    Allocate(Rhoy(0:nl(1)+1, 0:nl(2)+1, 0:nl(3)+1))
    Allocate(Rhoz(0:nl(1)+1, 0:nl(2)+1, 0:nl(3)+1))
    Allocate( U_v(0:nl(1)+1, 0:nl(2)+1, 0:nl(3)+1))
    Allocate(U_sx(0:nl(1)+1, 0:nl(2)+1, 0:nl(3)+1))
    Allocate(U_sy(0:nl(1)+1, 0:nl(2)+1, 0:nl(3)+1))
    Allocate(U_sz(0:nl(1)+1, 0:nl(2)+1, 0:nl(3)+1))
    Allocate( V_v(0:nl(1)+1, 0:nl(2)+1, 0:nl(3)+1))
    Allocate(V_sx(0:nl(1)+1, 0:nl(2)+1, 0:nl(3)+1))
    Allocate(V_sy(0:nl(1)+1, 0:nl(2)+1, 0:nl(3)+1))
    Allocate(V_sz(0:nl(1)+1, 0:nl(2)+1, 0:nl(3)+1))
    Allocate( W_v(0:nl(1)+1, 0:nl(2)+1, 0:nl(3)+1))
    Allocate(W_sx(0:nl(1)+1, 0:nl(2)+1, 0:nl(3)+1))
    Allocate(W_sy(0:nl(1)+1, 0:nl(2)+1, 0:nl(3)+1))
    Allocate(W_sz(0:nl(1)+1, 0:nl(2)+1, 0:nl(3)+1))
    Allocate(  Mu(0:nl(1)+1, 0:nl(2)+1, 0:nl(3)+1))
    Allocate(  LS(0:nl(1)+1, 0:nl(2)+1, 0:nl(3)+1))
    Allocate(Flag(0:nl(1)+1, 0:nl(2)+1, 0:nl(3)+1))
    Allocate(PP(0:nl(1)+1, 0:nl(2)+1, 0:nl(3)+1))
    Allocate( Div(nl(1), nl(2), nl(3)))
    Allocate(rkcoef(2,rk_order))

    Rho  = 0.0_sp
    Rhox = 0.0_sp
    Rhoy = 0.0_sp
    Rhoz = 0.0_sp
    Mu   = 0.0_sp
    Div  = 0.0_sp
    LS   = 0.0_sp
    U_v   = 0.0_sp
    U_sx   = 0.0_sp
    U_sy  = 0.0_sp
    U_sz  = 0.0_sp
    V_v   = 0.0_sp
    V_sx  = 0.0_sp
    V_sy  = 0.0_sp
    V_sz  = 0.0_sp
    W_v   = 0.0_sp
    W_sx  = 0.0_sp
    W_sy  = 0.0_sp
    W_sz  = 0.0_sp
    Flag = 1
    rkcoef   = 0.0_sp

    ls_bc%name = 'ls'
    ls_bc%lohi(:,1) = 0
    ls_bc%lohi(:,2) = nl(:)+1
    ls_bc%bound_type(:,:) = 3
    ls_bc%bound_value(:,:) = 0

    ! Boundary locations for via
    self%BC_Uv%name = 'u_v'
    self%BC_Uv%lohi(:,1) = 0
    self%BC_Uv%lohi(:,2) = nl(:)+1
    self%bc_Uv%ghost_flag(:,:) = .true.

    self%BC_Vv%name = 'v_v'
    self%BC_Vv%lohi(:,1) = 0
    self%BC_Vv%lohi(:,2) = nl(:)+1
    self%bc_Vv%ghost_flag(:,:) = .true.

    self%BC_Wv%name = 'w_v'
    self%BC_Wv%lohi(:,1) = 0
    self%BC_Wv%lohi(:,2) = nl(:)+1
    self%bc_Wv%ghost_flag(:,:) = .true.

    ! Boundary locations for sia x
    self%BC_Usx%name = 'u_sx'
    self%BC_Usx%lohi(:,1) = 0
    self%BC_Usx%lohi(:,2) = nl(:)+1
    self%BC_Usx%lohi(1,2) = nl(1)
    self%bc_Usx%ghost_flag(:,:) = .true.
    self%bc_Usx%ghost_flag(1,:) = .false.

    self%BC_Vsx%name = 'v_sx'
    self%BC_Vsx%lohi(:,1) = 0
    self%BC_Vsx%lohi(:,2) = nl(:)+1
    self%BC_Vsx%lohi(1,2) = nl(1)
    self%bc_Vsx%ghost_flag(:,:) = .true.
    self%bc_Vsx%ghost_flag(1,:) = .false.

    self%BC_Wsx%name = 'w_sx'
    self%BC_Wsx%lohi(:,1) = 0
    self%BC_Wsx%lohi(:,2) = nl(:)+1
    self%BC_Wsx%lohi(1,2) = nl(1)
    self%bc_Wsx%ghost_flag(:,:) = .true.
    self%bc_Wsx%ghost_flag(1,:) = .false.

    ! Boundary locations for sia y
    self%BC_Usy%name = 'u_sy'
    self%BC_Usy%lohi(:,1) = 0
    self%BC_Usy%lohi(:,2) = nl(:)+1
    self%BC_Usy%lohi(2,2) = nl(2)
    self%bc_Usy%ghost_flag(:,:) = .true.
    self%bc_Usy%ghost_flag(2,:) = .false.

    self%BC_Vsy%name = 'v_sy'
    self%BC_Vsy%lohi(:,1) = 0
    self%BC_Vsy%lohi(:,2) = nl(:)+1
    self%BC_Vsy%lohi(2,2) = nl(2)
    self%bc_Vsy%ghost_flag(:,:) = .true.
    self%bc_Vsy%ghost_flag(2,:) = .false.

    self%BC_Wsy%name = 'w_sy'
    self%BC_Wsy%lohi(:,1) = 0
    self%BC_Wsy%lohi(:,2) = nl(:)+1
    self%BC_Wsy%lohi(2,2) = nl(2)
    self%bc_Wsy%ghost_flag(:,:) = .true.
    self%bc_Wsy%ghost_flag(2,:) = .false.

    ! Boundary locations for sia z
    self%BC_Usz%name = 'u_sz'
    self%BC_Usz%lohi(:,1) = 0
    self%BC_Usz%lohi(:,2) = nl(:)+1
    self%BC_Usz%lohi(3,2) = nl(3)
    self%bc_Usz%ghost_flag(:,:) = .true.
    self%bc_Usz%ghost_flag(3,:) = .false.

    self%BC_Vsz%name = 'v_sz'
    self%BC_Vsz%lohi(:,1) = 0
    self%BC_Vsz%lohi(:,2) = nl(:)+1
    self%BC_Vsz%lohi(3,2) = nl(3)
    self%bc_Vsz%ghost_flag(:,:) = .true.
    self%bc_Vsz%ghost_flag(3,:) = .false.

    self%BC_Wsz%name = 'w_sz'
    self%BC_Wsz%lohi(:,1) = 0
    self%BC_Wsz%lohi(:,2) = nl(:)+1
    self%BC_Wsz%lohi(3,2) = nl(3)
    self%bc_Wsz%ghost_flag(:,:) = .true.
    self%bc_Wsz%ghost_flag(3,:) = .false.

    ! Boudnary conditions:
    ! 1. No-slip wall
    ! 2. Slip wall/ Symmetric

    ! Left and right
    bc_pair(1) = bc_left
    bc_pair(2) = bc_right
    Do i = 1, 2
      Select case (bc_pair(i))
      Case(1)
        p_bc%bound_type(1,i)  = 2
        p_bc%bound_value(1,i) = 0.0_sp
        phi_bc%bound_type(1,i)  = 2
        phi_bc%bound_value(1,i) = 0.0_sp

        self%BC_Uv%bound_type(1,i) = 1
        self%BC_Uv%bound_value(1,i) = 0.0_sp
        self%BC_Vv%bound_type(1,i) = 1
        self%BC_Vv%bound_value(1,i) = 0.0_sp
        self%BC_Wv%bound_type(1,i) = 1
        self%BC_Wv%bound_value(1,i) = 0.0_sp

        self%BC_Usx%bound_type(1,i) = 1
        self%BC_Usx%bound_value(1,i) = 0.0_sp
        self%BC_Vsx%bound_type(1,i) = 1
        self%BC_Vsx%bound_value(1,i) = 0.0_sp
        self%BC_Wsx%bound_type(1,i) = 1
        self%BC_Wsx%bound_value(1,i) = 0.0_sp

        self%BC_Usy%bound_type(1,i) = 1
        self%BC_Usy%bound_value(1,i) = 0.0_sp
        self%BC_Vsy%bound_type(1,i) = 1
        self%BC_Vsy%bound_value(1,i) = 0.0_sp
        self%BC_Wsy%bound_type(1,i) = 1
        self%BC_Wsy%bound_value(1,i) = 0.0_sp

        self%BC_Usz%bound_type(1,i) = 1
        self%BC_Usz%bound_value(1,i) = 0.0_sp
        self%BC_Vsz%bound_type(1,i) = 1
        self%BC_Vsz%bound_value(1,i) = 0.0_sp
        self%BC_Wsz%bound_type(1,i) = 1
        self%BC_Wsz%bound_value(1,i) = 0.0_sp
      Case(2)
        p_bc%bound_type(1,i)  = 2
        p_bc%bound_value(1,i) = 0.0_sp
        phi_bc%bound_type(1,i)  = 2
        phi_bc%bound_value(1,i) = 0.0_sp

        self%BC_Uv%bound_type(1,i) = 1
        self%BC_Uv%bound_value(1,i) = 0.0_sp
        self%BC_Vv%bound_type(1,i) = 2
        self%BC_Vv%bound_value(1,i) = 0.0_sp
        self%BC_Wv%bound_type(1,i) = 2
        self%BC_Wv%bound_value(1,i) = 0.0_sp

        self%BC_Usx%bound_type(1,i) = 1
        self%BC_Usx%bound_value(1,i) = 0.0_sp
        self%BC_Vsx%bound_type(1,i) = 2
        self%BC_Vsx%bound_value(1,i) = 0.0_sp
        self%BC_Wsx%bound_type(1,i) = 2
        self%BC_Wsx%bound_value(1,i) = 0.0_sp

        self%BC_Usy%bound_type(1,i) = 1
        self%BC_Usy%bound_value(1,i) = 0.0_sp
        self%BC_Vsy%bound_type(1,i) = 2
        self%BC_Vsy%bound_value(1,i) = 0.0_sp
        self%BC_Wsy%bound_type(1,i) = 2
        self%BC_Wsy%bound_value(1,i) = 0.0_sp

        self%BC_Usz%bound_type(1,i) = 1
        self%BC_Usz%bound_value(1,i) = 0.0_sp
        self%BC_Vsz%bound_type(1,i) = 2
        self%BC_Vsz%bound_value(1,i) = 0.0_sp
        self%BC_Wsz%bound_type(1,i) = 2
        self%BC_Wsz%bound_value(1,i) = 0.0_sp
      Case Default
        if (myid .eq.0) then
          print *, "======Fatal Error=============================="
          print *, "Wrong bc value, please check the namelist file"
          print *, "==============================================="
        end if
        Call MPI_Finalize(ierr)
        stop
      End Select
    End Do

    ! Front and back
    bc_pair(1) = bc_back
    bc_pair(2) = bc_front
    Do i = 1, 2
      Select case (bc_pair(i))
      Case(1)
        p_bc%bound_type(2,i)  = 2
        p_bc%bound_value(2,i) = 0.0_sp
        phi_bc%bound_type(2,i)  = 2
        phi_bc%bound_value(2,i) = 0.0_sp

        self%BC_Uv%bound_type(2,i) = 1
        self%BC_Uv%bound_value(2,i) = 0.0_sp
        self%BC_Vv%bound_type(2,i) = 1
        self%BC_Vv%bound_value(2,i) = 0.0_sp
        self%BC_Wv%bound_type(2,i) = 1
        self%BC_Wv%bound_value(2,i) = 0.0_sp

        self%BC_Usx%bound_type(2,i) = 1
        self%BC_Usx%bound_value(2,i) = 0.0_sp
        self%BC_Vsx%bound_type(2,i) = 1
        self%BC_Vsx%bound_value(2,i) = 0.0_sp
        self%BC_Wsx%bound_type(2,i) = 1
        self%BC_Wsx%bound_value(2,i) = 0.0_sp

        self%BC_Usy%bound_type(2,i) = 1
        self%BC_Usy%bound_value(2,i) = 0.0_sp
        self%BC_Vsy%bound_type(2,i) = 1
        self%BC_Vsy%bound_value(2,i) = 0.0_sp
        self%BC_Wsy%bound_type(2,i) = 1
        self%BC_Wsy%bound_value(2,i) = 0.0_sp

        self%BC_Usz%bound_type(2,i) = 1
        self%BC_Usz%bound_value(2,i) = 0.0_sp
        self%BC_Vsz%bound_type(2,i) = 1
        self%BC_Vsz%bound_value(2,i) = 0.0_sp
        self%BC_Wsz%bound_type(2,i) = 1
        self%BC_Wsz%bound_value(2,i) = 0.0_sp
      Case(2)
        p_bc%bound_type(2,i)  = 2
        p_bc%bound_value(2,i) = 0.0_sp
        phi_bc%bound_type(2,i)  = 2
        phi_bc%bound_value(2,i) = 0.0_sp

        self%BC_Uv%bound_type(2,i) = 2
        self%BC_Uv%bound_value(2,i) = 0.0_sp
        self%BC_Vv%bound_type(2,i) = 1
        self%BC_Vv%bound_value(2,i) = 0.0_sp
        self%BC_Wv%bound_type(2,i) = 1
        self%BC_Wv%bound_value(2,i) = 0.0_sp

        self%BC_Usx%bound_type(2,i) = 2
        self%BC_Usx%bound_value(2,i) = 0.0_sp
        self%BC_Vsx%bound_type(2,i) = 1
        self%BC_Vsx%bound_value(2,i) = 0.0_sp
        self%BC_Wsx%bound_type(2,i) = 2
        self%BC_Wsx%bound_value(2,i) = 0.0_sp

        self%BC_Usy%bound_type(2,i) = 2
        self%BC_Usy%bound_value(2,i) = 0.0_sp
        self%BC_Vsy%bound_type(2,i) = 1
        self%BC_Vsy%bound_value(2,i) = 0.0_sp
        self%BC_Wsy%bound_type(2,i) = 2
        self%BC_Wsy%bound_value(2,i) = 0.0_sp

        self%BC_Usz%bound_type(2,i) = 2
        self%BC_Usz%bound_value(2,i) = 0.0_sp
        self%BC_Vsz%bound_type(2,i) = 1
        self%BC_Vsz%bound_value(2,i) = 0.0_sp
        self%BC_Wsz%bound_type(2,i) = 2
        self%BC_Wsz%bound_value(2,i) = 0.0_sp
      End Select
    End Do

    ! Top and bottom
    bc_pair(1) = bc_bottom
    bc_pair(2) = bc_top
    Do i = 1, 2
      Select case (bc_pair(i))
      Case(1)
        p_bc%bound_type(3,i)  = 2
        p_bc%bound_value(3,i) = 0.0_sp
        phi_bc%bound_type(3,i)  = 2
        phi_bc%bound_value(3,i) = 0.0_sp

        self%BC_Uv%bound_type(3,i) = 1
        self%BC_Uv%bound_value(3,i) = 0.0_sp
        self%BC_Vv%bound_type(3,i) = 1
        self%BC_Vv%bound_value(3,i) = 0.0_sp
        self%BC_Wv%bound_type(3,i) = 1
        self%BC_Wv%bound_value(3,i) = 0.0_sp

        self%BC_Usx%bound_type(3,i) = 1
        self%BC_Usx%bound_value(3,i) = 0.0_sp
        self%BC_Vsx%bound_type(3,i) = 1
        self%BC_Vsx%bound_value(3,i) = 0.0_sp
        self%BC_Wsx%bound_type(3,i) = 1
        self%BC_Wsx%bound_value(3,i) = 0.0_sp

        self%BC_Usy%bound_type(3,i) = 1
        self%BC_Usy%bound_value(3,i) = 0.0_sp
        self%BC_Vsy%bound_type(3,i) = 1
        self%BC_Vsy%bound_value(3,i) = 0.0_sp
        self%BC_Wsy%bound_type(3,i) = 1
        self%BC_Wsy%bound_value(3,i) = 0.0_sp

        self%BC_Usz%bound_type(3,i) = 1
        self%BC_Usz%bound_value(3,i) = 0.0_sp
        self%BC_Vsz%bound_type(3,i) = 1
        self%BC_Vsz%bound_value(3,i) = 0.0_sp
        self%BC_Wsz%bound_type(3,i) = 1
        self%BC_Wsz%bound_value(3,i) = 0.0_sp
      Case(2)
        p_bc%bound_type(3,i)  = 2
        p_bc%bound_value(3,i) = 0.0_sp
        phi_bc%bound_type(3,i)  = 2
        phi_bc%bound_value(3,i) = 0.0_sp

        self%BC_Uv%bound_type(3,i) = 2
        self%BC_Uv%bound_value(3,i) = 0.0_sp
        self%BC_Vv%bound_type(3,i) = 2
        self%BC_Vv%bound_value(3,i) = 0.0_sp
        self%BC_Wv%bound_type(3,i) = 1
        self%BC_Wv%bound_value(3,i) = 0.0_sp

        self%BC_Usx%bound_type(3,i) = 2
        self%BC_Usx%bound_value(3,i) = 0.0_sp
        self%BC_Vsx%bound_type(3,i) = 2
        self%BC_Vsx%bound_value(3,i) = 0.0_sp
        self%BC_Wsx%bound_type(3,i) = 1
        self%BC_Wsx%bound_value(3,i) = 0.0_sp

        self%BC_Usy%bound_type(3,i) = 2
        self%BC_Usy%bound_value(3,i) = 0.0_sp
        self%BC_Vsy%bound_type(3,i) = 2
        self%BC_Vsy%bound_value(3,i) = 0.0_sp
        self%BC_Wsy%bound_type(3,i) = 1
        self%BC_Wsy%bound_value(3,i) = 0.0_sp

        self%BC_Usz%bound_type(3,i) = 2
        self%BC_Usz%bound_value(3,i) = 0.0_sp
        self%BC_Vsz%bound_type(3,i) = 2
        self%BC_Vsz%bound_value(3,i) = 0.0_sp
        self%BC_Wsz%bound_type(3,i) = 1
        self%BC_Wsz%bound_value(3,i) = 0.0_sp
      End Select
    End Do

    ! Assign RungeKutta
    If (rk_order .eq. 1) Then
      rkcoef(1,1) =  1.0_sp
      rkcoef(2,1) =  0.0_sp
    ElseIf (rk_order .eq. 2) Then
      rkcoef(1,1) =  0.5_sp
      rkcoef(2,1) =  0.0_sp
      rkcoef(1,2) =  1.0_sp
      rkcoef(2,2) =  -0.5_sp
    ElseIf (rk_order .eq. 3) Then
      rkcoef(1,1) =  32.0_sp / 60.0_sp
      rkcoef(2,1) =  0.0_sp  / 60.0_sp
      rkcoef(1,2) =  25.0_sp / 60.0_sp
      rkcoef(2,2) = -17.0_sp / 60.0_sp
      rkcoef(1,3) =  45.0_sp / 60.0_sp
      rkcoef(2,3) = -25.0_sp / 60.0_sp
    Else
      If ( myid .eq. 0 ) then
        print *, "======Fatal Error=============================="
        print *, "Wong rk order, should be 1 or 2 or 3"
        print *, "==============================================="
      End If
      Call MPI_Finalize(ierr)
      stop
    End If

    ! Initialize hypre
    Call Hypre_Initialize(iter_tolerance, iter_max, hypre_solver, Hypre_PreConditioner)

    ! Initialize the Rho and Mu field
    ! Call self%UpdtRhoMu(Phi)

    ! Initialize boundary condition
    ! Call self%BC_UVW(U,V,W)
    ! Call self%BC_P(P)

  End Subroutine InitNavierStokes_VSIAM

 Subroutine Monitor(U,V,W)
    Implicit None
    Real(sp), Intent(In), Dimension(0:,0:,0:) :: U, V, W
    Integer :: i, j, k
    External Divergence
    Call Divergence(U,V,W,nl,dl,div)
    div_max = 0.0_sp
    Do k = 1, nl(3)
      Do j = 1, nl(2)
        Do i = 1, nl(1)
          div_max = max(div_max,div(i,j,k))
        End Do
      End Do
    End Do
  End Subroutine Monitor

End Module ModNavierStokes
