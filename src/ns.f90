!===============================================
! Includes:
!     (1) Time integration
!       (1-1) Two phase flow solver
!       (1-2) Runge Kutta for momentum
!     (2) Solution procedures
!       (2-1) Advection-Diffusion
!          (2-1-1) Advection-Diffusion for U
!          (2-1-1) Advection-Diffusion for V
!          (2-1-1) Advection-Diffusion for W
!       (2-2) Partial P
!       (2-3) Divergence
!       (2-4) Projection
!       (2-5) Update Rho and Mu
!       (2-6) Continuum surface tension
!     (3) Initial and Boundary conditions
!       (3-1) Boundary condition for u, v, w
!       (3-2) Boundary condition for p
!       (3-3) Boundary condition for rho, mu
!       (3-4) Initialization
!     (4) Misc
!       (4-2) Adjust_dt
!       (4-1) Jacobi iteration
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
  Use ModTools
  Use HypreStruct
  Use ModVOF
  Implicit None

  private
  Public :: TwoPhaseFlow, SinglePhaseFlow, InitNavierStokes, Monitor
  Public :: CSF, CSF_VOSET
  ! Module variables will not pass through the subroutine interface
  Real(sp), Allocatable :: Rho(:,:,:), Rhox(:,:,:), Rhoy(:,:,:), Rhoz(:,:,:), PP(:,:,:)
  Real(sp), Allocatable :: mu(:,:,:), muxy(:,:,:), muxz(:,:,:), muyz(:,:,:)
  Real(sp), Allocatable :: Div(:,:,:)
  Integer, Allocatable :: flag(:,:,:)
  Real(sp) :: body_force(3)
  Real(sp) :: rho_l, rho_g
  Real(sp) ::  mu_l, mu_g
  Real(sp) ::  sigma
  Integer :: rk_order
  Real(sp), Allocatable :: rkcoef(:,:)
  Integer :: n_iter
  Real(sp) :: cfl, dt_min
  Real(sp) :: p_residual
  Real(sp) :: div_max
  Integer :: step_rank = 0
  Integer :: out_nn = 0
  Real(sp) :: time_out
  integer :: jacobiflag = 0

  PROCEDURE(InterfaceCSF), POINTER :: CSF => CSF_NULL
  Interface
    Subroutine InterfaceCSF(dudt, dvdt, dwdt, Phi)
      import
      Implicit None
      real(sp), intent(in), dimension(0:,0:,0:) :: Phi
      real(sp), intent(out), dimension(:,:,:) :: dudt, dvdt, dwdt
    End Subroutine InterfaceCSF
  End Interface

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
  Subroutine TwoPhaseFlow(U, V, W, Phi, P, cx, cy, cz)
    Implicit None
    real(sp), Intent(inout), Dimension(0:,0:,0:) :: u, v, w, p, phi
    real(sp), Intent(inout), Dimension(0:,0:,0:), optional :: cx, cy, cz
    real(sp), dimension(nl(1),nl(2),nl(3))    :: dudtrko,dvdtrko,dwdtrko
    real(sp), dimension(0:nl(1)+1,0:nl(2)+1,0:nl(3)+1)    :: up,vp,wp
    Integer :: irk

    ! Runge-Kutta time integration
    Do irk=1,rk_order
      ! Runge-Kutta step for momentum
      call RungeKuttaMomentum(rkcoef(:,irk), u, v, w, p, dudtrko,dvdtrko,dwdtrko, up, vp, wp,phi)
      Call Visual3DQuiver(VP,UP,WP)
      ! Boundary condition for intermediate velocity
      Call BC_UVW(UP, VP, WP)
      ! Calculate divergence
      Call Divergence(up, vp, wp, Div)
      ! Solve pressure Poisson equation
      Call Hypre_Poisson(PP, Rhox, Rhoy, Rhoz, Div, flag, n_iter, p_residual)
      Call BC_P(PP)
      ! Preject the velocity to the divergence-free field
      Call Projection(PP, u, v, w, up, vp, wp)
      ! Boundary condition for divergence-free velocity field
      Call BC_UVW(U, V, W)
      p(1:nl(1),1:nl(2),1:nl(3)) = p(1:nl(1),1:nl(2),1:nl(3)) + pp(1:nl(1),1:nl(2),1:nl(3))
      Call BC_P(P)
    End Do

    ! VOF/MOF advection
    step_rank = step_rank+1
    If (present(cx)) Then ! MOF
      Call VOFAdvection(Phi, u, v, w, nl, dl, dt, mod(step_rank,3), cx, cy, cz)
    Else !VOF
      Call VOFAdvection(Phi, u, v, w, nl, dl, dt, mod(step_rank,3))
    End If

    ! Update Rho, Mu, and adjust time step
    Call UpdtRhoMu(Phi)
    Call Adjustdt(U, V, W, cfl)

  End Subroutine TwoPhaseFlow

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
  !==========================
  Subroutine SinglePhaseFlow(U, V, W, P)
    Implicit None
    real(sp), Intent(inout), Dimension(0:,0:,0:) :: u, v, w, p
    real(sp), dimension(nl(1),nl(2),nl(3))    :: dudtrko,dvdtrko,dwdtrko
    real(sp), dimension(0:nl(1)+1,0:nl(2)+1,0:nl(3)+1)    :: up,vp,wp
    Integer :: irk

    ! Runge-Kutta time integration
    Do irk=1,rk_order
      ! Runge-Kutta step for momentum
      call RungeKuttaMomentum(rkcoef(:,irk), u, v, w, p, dudtrko,dvdtrko,dwdtrko, up, vp, wp)
      ! Boundary condition for intermediate velocity
      Call BC_UVW(UP, VP, WP)
      ! Calculate divergence
      Call Divergence(up, vp, wp, Div)
      ! Solve pressure Poisson equation
      ! if (jacobiflag .eq.0) then
      Call Hypre_Poisson(PP, Rhox, Rhoy, Rhoz, Div, flag, n_iter, p_residual)
        ! jacobiflag = 1
      ! else
        ! Call Jacobi(PP)
      ! endif
      Call BC_P(PP)
      ! Preject the velocity to the divergence-free field
      Call Projection(PP, u, v, w, up, vp, wp)
      ! Boundary condition for divergence-free velocity field
      Call BC_UVW(U, V, W)
      p(1:nl(1),1:nl(2),1:nl(3)) = p(1:nl(1),1:nl(2),1:nl(3)) + pp(1:nl(1),1:nl(2),1:nl(3))
      Call BC_P(P)
    End Do

    Call Adjustdt(U, V, W, cfl)

  End Subroutine SinglePhaseFlow

  !==========================
  ! (1-2) Runge Kutta time integration 
  ! Runge Kutta step for
  ! U^star = -1 / \rho * \nabla p^k-1/2 + \alpha * AD^k-1/2 + \beta AD^k-1/2
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
  Subroutine RungeKuttaMomentum(rkpar, u, v, w, p, dudtrko, dvdtrko, dwdtrko, up, vp, wp, Phi)
    Implicit None
    real(sp), Intent(in), Dimension(0:,0:,0:) :: u,v,w, p
    real(sp), Intent(inout), Dimension(:,:,:) :: dudtrko,dvdtrko,dwdtrko
    real(sp), Intent(in), Dimension(2) :: rkpar
    real(sp), Intent(out), Dimension(0:,0:,0:)  :: up,vp,wp
    real(sp), Intent(in), Dimension(0:,0:,0:),optional  :: phi
    real(sp), Dimension(nl(1),nl(2),nl(3)) :: dudtrk,dvdtrk,dwdtrk
    real(sp) :: factor1,factor2,factor12
    Integer :: i, j, k

    factor1 = rkpar(1)*dt
    factor2 = rkpar(2)*dt
    factor12 = factor1 + factor2

    ! Calculate -1 / \tho * \nabla p^k-1/2
    Call PartialP(P, Rhox, dudtrk, 1)
    Call PartialP(P, Rhoy, dvdtrk, 2)
    Call PartialP(P, Rhoz, dwdtrk, 3)
    do k=1,nl(3)
      do j=1,nl(2)
        do i=1,nl(1)
          up(i,j,k) = u(i,j,k) + factor12*dudtrk(i,j,k)
          vp(i,j,k) = v(i,j,k) + factor12*dvdtrk(i,j,k)
          wp(i,j,k) = w(i,j,k) + factor12*dwdtrk(i,j,k)
        enddo
      enddo
    enddo

    ! Calculate AD^k+1/2
    Call AdvDiffU(u, v, w, dudtrk)
    Call AdvDiffV(u, v, w, dvdtrk)
    Call AdvDiffW(u, v, w, dwdtrk)
    ! Calculate Body force k+1/2 (Included in AD^k+1/2 term)
    Call BodyForce(dudtrk, dvdtrk, dwdtrk)

    if ( present(phi) ) Call CSF(dudtrk, dvdtrk, dwdtrk, Phi)

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

  End Subroutine RungeKuttaMomentum

  !==========================
  ! (2-1-1) Calculate advection and diffusion
  ! Using central difference
  !-------------------------
  !    d(u)/dt_i+1/2 =
  !     ( (uu)_i+1 - (uu)_i + ( (mu*dudx)_i+1   - (mu*dudx)_i   ) / rhox_i+1/2 ) / dx
  !     ( (uv)_j+1 - (uv)_j + ( (mu*dudy)_j+1/2 - (mu*dudy)_j/2 ) / rhoy_j+1/2 ) / dy
  !     ( (uw)_k+1 - (uw)_k + ( (mu*dudz)_k+1/2 - (mu*dudz)_k/2 ) / rhoz_k+1/2 ) / dz
  !  where:
  !      (uu)_i = 0.25 * (u_i+1/2 + u_i-1/2)^2
  !      (uv)_j = 0.25 * (u_i+1/2 + u_i-1/2) * (v_j+1/2,v_j-1/2)
  !      (uw)_k = 0.25 * (u_i+1/2 + u_i-1/2) * (w_k+1/2,w_k-1/2)
  !      (dudx)_i = (u_i+1/2 + u_i-1/2) / dx
  !      (dudy)_j+1/2 = (u_j+1   + u_j) / dy
  !      (dudz)_k+1/2 = (u_k+1   + u_k) / dz
  !--------------------------
  ! Inputs:
  !    u, v, w: velocity conponent
  !    rho_dir: the reciprocal of rho, depends on direction
  !    mu_c: mu at cell center
  !    nl: dimension of grid
  !    dl: grid size dx, dy, dz
  ! Outputs:
  !    dudt: time integration of advection and diffusion
  !--------------------------
  ! Inputs:
  !   U, V, W: velocity conponents
  ! Outpurts
  !   dudtr: du/dt
  !==========================
  Subroutine AdvDiffU(u, v, w, dudt)
    Implicit None
    Real(sp), Dimension(0:,0:,0:), Intent(in)  :: u, v, w
    Real(sp), Dimension(:,:,:), Intent(out) :: dudt
    integer :: i,j,k
    Real(sp) :: uuip,uuim,uvjp,uvjm,uwkp,uwkm
    Real(sp) :: dudxp,dudxm,dudyp,dudym,dudzp,dudzm

    ! Calculate the advection-diffusion spatial term
    Do k=1,nl(3)
      Do j=1,nl(2)
        Do i=1,nl(1)
          ! uu, uv, uw
          uuip = 0.25_sp * (u(i+1,  j,  k) + u(i,j,k)) * (u(i+1,  j,  k) + u(i,  j,  k))
          uuim = 0.25_sp * (u(i-1,  j,  k) + u(i,j,k)) * (u(i-1,  j,  k) + u(i,  j,  k))
          uvjp = 0.25_sp * (u(  i,j+1,  k) + u(i,j,k)) * (v(i+1,  j,  k) + v(i,  j,  k))
          uvjm = 0.25_sp * (u(  i,j-1,  k) + u(i,j,k)) * (v(i+1,j-1,  k) + v(i,j-1,  k))
          uwkp = 0.25_sp * (u(  i,  j,k+1) + u(i,j,k)) * (w(i+1,  j,  k) + w(i,  j,  k))
          uwkm = 0.25_sp * (u(  i,  j,k-1) + u(i,j,k)) * (w(i+1,  j,k-1) + w(i,  j,k-1))
          ! dudx, dudy, dudz
          dudxp =   mu(i+1,j,k) * (u(i+1,  j,  k) - u(  i,  j,  k)) / dl(1)
          dudxm =   mu(i,j,k)   * (u(  i,  j,  k) - u(i-1,  j,  k)) / dl(1)
          dudyp = muxz(i,j,k)   * (u(  i,j+1,  k) - u(  i,  j,  k)) / dl(2)
          dudym = muxz(i,j-1,k) * (u(  i,  j,  k) - u(  i,j-1,  k)) / dl(2)
          dudzp = muxy(i,j,k)   * (u(  i,  j,k+1) - u(  i,  j,  k)) / dl(3)
          dudzm = muxy(i,j,k-1) * (u(  i,  j,  k) - u(  i,  j,k-1)) / dl(3)
          ! Momentum balance
          dudt(i,j,k) = &
              ((-uuip + uuim) + (dudxp-dudxm)) * rhox(i,j,k) / dl(1) + &
              ((-uvjp + uvjm) + (dudyp-dudym)) * rhoy(i,j,k) / dl(2) + &
              ((-uwkp + uwkm) + (dudzp-dudzm)) * rhoz(i,j,k) / dl(3)
        Enddo
      Enddo
    Enddo

  End Subroutine AdvDiffU

  !==========================
  ! (2-1-2) Calculate advection and diffusion for V
  ! similar with  AdvDiffU
  !--------------------------
  ! Inputs:
  !   U, V, W: velocity conponents
  ! Outpurts
  !   dvdt: dv/dt
  !==========================
  Subroutine AdvDiffV(u, v, w, dvdt)
    Implicit None
    Real(sp), Dimension(0:,0:,0:), Intent(in)  :: u, v, w
    Real(sp), Dimension(:,:,:), Intent(out) :: dvdt
    integer :: i,j,k
    real(sp) :: uvip,uvim,vvjp,vvjm,wvkp,wvkm
    real(sp) :: dvdxp,dvdxm,dvdyp,dvdym,dvdzp,dvdzm

    ! Calculate the advection-diffusion spatial term
    Do k=1,nl(3)
      Do j=1,nl(2)
        Do i=1,nl(1)
          ! uu, uv, uw
          uvip = 0.25_sp * (u(  i,  j,  k) + u(  i,j+1,  k)) * (v(i,j,k) + v(i+1,  j,  k))
          uvim = 0.25_sp * (u(i-1,  j,  k) + u(i-1,j+1,  k)) * (v(i,j,k) + v(i-1,  j,  k))
          vvjp = 0.25_sp * (v(  i,  j,  k) + v(  i,j+1,  k)) * (v(i,j,k) + v(  i,j+1,  k))
          vvjm = 0.25_sp * (v(  i,  j,  k) + v(  i,j-1,  k)) * (v(i,j,k) + v(  i,j-1,  k))
          wvkp = 0.25_sp * (w(  i,  j,  k) + w(  i,j+1,  k)) * (v(i,j,k) + v(  i,  j,k+1))
          wvkm = 0.25_sp * (w(  i,  j,k-1) + w(  i,j+1,k-1)) * (v(i,j,k) + v(  i,  j,k-1))
          ! dudx, dudy, dudz
          dvdxp = muyz(i,j,k)   * (v(i+1,  j,  k) - v(  i,  j,  k)) / dl(1)
          dvdxm = muyz(i-1,j,k) * (v(  i,  j,  k) - v(i-1,  j,  k)) / dl(1)
          dvdyp =   mu(i,j+1,k) * (v(  i,j+1,  k) - v(  i,  j,  k)) / dl(2)
          dvdym =   mu(i,j,k)   * (v(  i,  j,  k) - v(  i,j-1,  k)) / dl(2)
          dvdzp = muxy(i,j,k)   * (v(  i,  j,k+1) - v(  i,  j,  k)) / dl(3)
          dvdzm = muxy(i,j,k-1) * (v(  i,  j,  k) - v(  i,  j,k-1)) / dl(3)
          ! Momentum balance
          dvdt(i,j,k) = &
              ((-uvip + uvim) + (dvdxp-dvdxm)) * rhox(i,j,k) / dl(1) + &
              ((-vvjp + vvjm) + (dvdyp-dvdym)) * rhoy(i,j,k) / dl(2) + &
              ((-wvkp + wvkm) + (dvdzp-dvdzm)) * rhoz(i,j,k) / dl(3)
        Enddo
      Enddo
    Enddo

  End Subroutine AdvDiffV

  !==========================
  ! (2-1-3) Calculate advection and diffusion for V
  ! similar with  AdvDiffU
  !--------------------------
  ! Inputs:
  !   U, V, W: velocity conponents
  ! Outpurts
  !   dvdt: dv/dt
  !==========================
  Subroutine AdvDiffW(u, v, w, dwdt)
    Implicit None
    Real(sp), Dimension(0:,0:,0:), Intent(in)  :: u, v, w
    Real(sp), Dimension(:,:,:), Intent(out) :: dwdt
    integer :: i,j,k
    real(sp) :: uwip,uwim,vwjp,vwjm,wwkp,wwkm
    real(sp) :: dwdxp,dwdxm,dwdyp,dwdym,dwdzp,dwdzm

    ! Calculate the advection-diffusion spatial term
    Do k=1,nl(3)
      Do j=1,nl(2)
        Do i=1,nl(1)
          ! uu, uv, uw
          uwip = 0.25_sp*(w(i,j,k) + w(i+1,  j,  k)) * (u(i  ,j  ,k) + u(  i,  j,k+1))
          uwim = 0.25_sp*(w(i,j,k) + w(i-1,  j,  k)) * (u(i-1,j  ,k) + u(i-1,  j,k+1))
          vwjp = 0.25_sp*(w(i,j,k) + w(  i,j+1,  k)) * (v(i  ,j  ,k) + v(  i,  j,k+1))
          vwjm = 0.25_sp*(w(i,j,k) + w(  i,j-1,  k)) * (v(i  ,j-1,k) + v(  i,j-1,k+1))
          wwkp = 0.25_sp*(w(i,j,k) + w(  i,  j,k+1)) * (w(i  ,j  ,k) + w(  i,  j,k+1))
          wwkm = 0.25_sp*(w(i,j,k) + w(  i,  j,k-1)) * (w(i  ,j  ,k) + w(  i,  j,k-1))
          ! dudx, dudy, dudz
          dwdxp = muyz(i,j,  k) * (w(i+1,  j,  k) - w(  i,  j,  k)) / dl(1)
          dwdxm = muyz(i-1,j,k) * (w(  i,  j,  k) - w(i-1,  j,  k)) / dl(1)
          dwdyp = muxz(i,j,  k) * (w(  i,j+1,  k) - w(  i,  j,  k)) / dl(2)
          dwdym = muxz(i,j-1,k) * (w(  i,  j,  k) - w(  i,j-1,  k)) / dl(2)
          dwdzp =   mu(i,j,k+1) * (w(  i,  j,k+1) - w(  i,  j,  k)) / dl(3)
          dwdzm =   mu(i,j,  k) * (w(  i,  j,  k) - w(  i,  j,k-1)) / dl(3)
          ! Momentum balance
          dwdt(i,j,k) = &
              ((-uwip + uwim) + (dwdxp-dwdxm)) * rhox(i,j,k) / dl(1) + &
              ((-vwjp + vwjm) + (dwdyp-dwdym)) * rhoy(i,j,k) / dl(2) + &
              ((-wwkp + wwkm) + (dwdzp-dwdzm)) * rhoz(i,j,k) / dl(3)
        Enddo
      Enddo
    Enddo

  End Subroutine AdvDiffW

  !==========================
  ! (2-1) Calculate body force
  !    only support uniformly body force, for example, gravity
  !--------------------------
  ! Inputs & outputs:
  !    dudt, dvdt, dwdt: du/dt, dv/dt, dwdt
  !==========================
  Subroutine BodyForce(dudt, dvdt, dwdt)
    Implicit None
    Real(sp), Dimension(:,:,:), Intent(inout)  :: dudt, dvdt,dwdt
    Integer :: i, j, k
    Do k=1,nl(3)
      Do j=1,nl(2)
        Do i=1,nl(1)
          dudt(i,j,k) = dudt(i,j,k) + body_force(1)
          dvdt(i,j,k) = dvdt(i,j,k) + body_force(2)
          dwdt(i,j,k) = dwdt(i,j,k) + body_force(3)
        EndDo
      Enddo
    Enddo
  End Subroutine BodyForce

  !==========================
  ! (2-2) Calculate partial P by
  ! 1D view
  ! du/dt_i = - (p_i+1/2 - p_i-1/2) / dx
  !--------------------------
  ! Inputs:
  !    dir: direction
  !    p: pressure
  ! Outputs:
  !    dudt: time integration of pressure
  !==========================
  Subroutine PartialP(P, Rho_dir, dudt, dir)
    Integer, intent(in) :: dir
    real(sp), dimension(0:,0:,0:), intent(in) :: p
    real(sp), dimension(0:,0:,0:), intent(in) :: rho_dir
    real(sp), dimension(:,:,:), intent(out) :: dudt
    Integer :: i, j, k
    Integer :: ip, jp, kp

    If (dir .eq. 1) Then
      ip = 1; jp = 0; kp = 0
    ElseIf (dir .eq. 2) Then
      ip = 0; jp = 1; kp = 0
    Else
      ip = 0; jp = 0; kp = 1
    End If

    do k=1, nl(3)
      do j=1, nl(2)
        do i=1, nl(1)
          dudt(i,j,k) = - rho_dir(i,j,k) *  (p(i+ip,j+jp,k+kp) - p(i,j,k)) / dl(dir)
        enddo
      enddo
    enddo

  End Subroutine PartialP

  !==========================
  ! (2-3) Calculate divergence by
  ! div =  ( u_i+1/2 - u_i-1/2 ) / dx
  !        ( v_j+1/2 - v_j-1/2 ) / dy
  !        ( w_k+1/2 - w_k-1/2 ) / dz
  !--------------------------
  ! Inputs:
  !    up, vp, wp: velocity field
  ! Outputs:
  !    div: Divergence
  !==========================
  Subroutine Divergence(up, vp, wp, div)
    Implicit None
    real(sp), intent(in) , dimension(0:,0:,0:) :: up, vp, wp
    real(sp), intent(out) , dimension(:,:,:) :: div
    Integer :: i, j, k
    Do k = 1, nl(3)
      Do j = 1, nl(2)
        Do i = 1, nl(1)
          Div(i,j,k) = &
              ( up(i,j,k) - up(i-1,j,k) ) / dl(1) + &
              ( vp(i,j,k) - vp(i,j-1,k) ) / dl(2) + &
              ( wp(i,j,k) - wp(i,j,k-1) ) / dl(3)
        End Do
      End Do
    End Do
  End Subroutine Divergence

  !==========================
  ! (2-4) Projection method of Chorin
  ! Corrects the velocity field by pressure gradient
  ! 1-D view:
  ! u_i = up_i - (p_i+1/2 - p_i-1/2) * dt / dx
  !--------------------------
  ! Inputs:
  !    p: Pressure from Poisson equatiosn
  !    up, vp, wp: Intermediate velocity field
  ! Outputs:
  !    u, v, w: Divergence-free velocity field
  !==========================
  Subroutine Projection(P, u, v, w, up, vp, wp)
    Implicit None
    real(sp), intent(in) , dimension(0:,0:,0:) :: p,up,vp,wp
    real(sp), intent(out), dimension(0:,0:,0:) :: u,v,w
    Integer :: i, j, k

    do k=1,nl(3)
      do j=1,nl(2)
        do i=1,nl(1)
          u(i,j,k) = up(i,j,k) - rhox(i,j,k) * (p(i+1,  j,  k) - p(i,j,k)) * dt / dl(1)
          v(i,j,k) = vp(i,j,k) - rhoy(i,j,k) * (p(  i,j+1,  k) - p(i,j,k)) * dt / dl(2)
          w(i,j,k) = wp(i,j,k) - rhoz(i,j,k) * (p(  i,  j,k+1) - p(i,j,k)) * dt / dl(3)
        enddo
      enddo
    enddo

  End Subroutine Projection

  !==========================
  ! (2-5) Update Rho and Mu
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
  Subroutine UpdtRhoMu(Phi)
    Implicit None
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

    Call BC_RHOMU(Rho, Mu)

    Do k=0,nl(3)
      Do j=0,nl(2)
        Do i=0,nl(1)
          Rhox(i,j,k) = 2.0_sp / (Rho(i,j,k) + Rho(i+1,j,k))
          Rhoy(i,j,k) = 2.0_sp / (Rho(i,j,k) + Rho(i,j+1,k))
          Rhoz(i,j,k) = 2.0_sp / (Rho(i,j,k) + Rho(i,j,k+1))
          Muxy(i,j,k) = 0.25_sp * (Mu(i,j,k) + Mu(i+1,  j,k  ) + Mu(i,j+1,  k) + Mu(i+1,j+1,  k))
          Muyz(i,j,k) = 0.25_sp * (Mu(i,j,k) + Mu(  i,j+1,k  ) + Mu(i,  j,k+1) + Mu(  i,j+1,k+1))
          Muxz(i,j,k) = 0.25_sp * (Mu(i,j,k) + Mu(i+1,  j,  k) + Mu(i,  j,k+1) + Mu(i+1,  j,k+1))
        EndDo
      EndDo
    EndDo

  End Subroutine UpdtRhoMu

  !===============================================================
  !  Density-scaled Continuum Surface Force (CSF)
  !  with balanced force formation
  !  Key steps:
  !    (1) Scale volume function
  !    (2) Convert the VOF function to level set
  !    (3) calculate curvature at grid center
  !    (4) Interpolate the curvature to grid center
  !    (5) Calculte density-scaled balanced surface tension
  !
  !---------------------------------------------------------------
  ! Inputs:
  !   Phi: volume fraction
  !   sigma: surface tension coefficient
  ! Outputs:
  !   dudt, dvdt, dwdt: du/dt, dv/dt, dw/dt
  !===============================================================
  Subroutine CSF_VOSET(dudt, dvdt, dwdt, Phi)
    Implicit None
    real(sp), intent(in), dimension(0:,0:,0:) :: Phi
    real(sp), intent(out), dimension(:,:,:) :: dudt, dvdt, dwdt
    real(sp), dimension(0:nl(1)+1,0:nl(2)+1,0:nl(3)+1) :: Phi_DS
    real(sp), dimension(0:nl(1)+1,0:nl(2)+1,0:nl(3)+1) :: Ls
    real(sp), dimension(0:nl(1)+1,0:nl(2)+1,0:nl(3)+1) :: kappa_v
    real(sp), dimension(0:nl(1)+1,0:nl(2)+1,0:nl(3)+1) :: kappa_sx, kappa_sy, kappa_sz
    Integer :: i, j, k

    !    (1) Scale volume function
    Call Smooth_Phi(Phi, Phi_DS, nl(1), nl(2), nl(3))
    Call Phi_bc%SetBCS(Phi_DS)

    !    (2) Convert the VOF function to level set
    Call Get_LS_3D_init(Phi_DS, Ls, dl(1), dl(2), dl(3), nl(1), nl(2), nl(3))

    Call Phi_bc%SetBCS(LS)
    !    (3) calculate curvature at grid center
    Kappa_v = 0.0_sp
    Do k=0,nl(3)
      Do j=0,nl(2)
        Do i=0,nl(1)
          Kappa_v(i,j,k) = &
              (    ( Ls(i,j,k) - 2.d0 * Ls(i,j,k) + Ls(i+1,j,k) ) / dl(1) / dl(1) &
              &  + ( Ls(i,j,k) - 2.d0 * Ls(i,j,k) + Ls(i,j+1,k) ) / dl(2) / dl(2) &
              &  + ( Ls(i,j,k) - 2.d0 * Ls(i,j,k) + Ls(i,j,k+1) ) / dl(3) / dl(3) )
        EndDo
      EndDo
    EndDo
    Call Phi_bc%SetBCS(kappa_v)

    !    (4) Interpolate the curvature to grid center
    Do k=0,nl(3)
      Do j=0,nl(2)
        Do i=0,nl(1)
          Kappa_sx(i,j,k) = 0.5_sp * ( Kappa_v(i,j,k) + Kappa_v(i+1,  j,  k) )
          Kappa_sy(i,j,k) = 0.5_sp * ( Kappa_v(i,j,k) + Kappa_v(  i,j+1,  k) )
          Kappa_sz(i,j,k) = 0.5_sp * ( Kappa_v(i,j,k) + Kappa_v(  i,  j,k+1) )
        EndDo
      EndDo
    EndDo

    !    (5) Calculte density-scaled balanced surface tension
    Do k=1,nl(3)
      Do j=1,nl(2)
        Do i=1,nl(1)
          dudt(i,j,k) = dudt(i,j,k) - sigma * Kappa_sx(i,j,k) * ( Phi_DS(i+1,  j,  k) - Phi_DS(i,j,k) ) * rhox(i,j,k) / dl(1) * dt
          dvdt(i,j,k) = dvdt(i,j,k) - sigma * Kappa_sy(i,j,k) * ( Phi_DS(  i,j+1,  k) - Phi_DS(i,j,k) ) * rhoy(i,j,k) / dl(2) * dt
          dwdt(i,j,k) = dvdt(i,j,k) - sigma * Kappa_sz(i,j,k) * ( Phi_DS(  i,  j,k+1) - Phi_DS(i,j,k) ) * rhoz(i,j,k) / dl(3) * dt
        EndDo
      EndDo
  EndDo

  End Subroutine CSF_VOSET

  Subroutine CSF_NULL(dudt, dvdt, dwdt, Phi)
    Implicit None
    real(sp), intent(in), dimension(0:,0:,0:) :: Phi
    real(sp), intent(out), dimension(:,:,:) :: dudt, dvdt, dwdt
  End Subroutine CSF_NULL

  !==========================
  ! (3-1) Boudnary conditions for u, v, w
  !==========================
  Subroutine BC_UVW(U, V, W)
    Implicit None
    real(sp), intent(inout) , dimension(0:,0:,0:) :: U, V, W
    Call U_bc%SetBCS(U)
    Call V_bc%SetBCS(V)
    Call W_bc%SetBCS(W)

    ! If(left   .eq. MPI_PROC_NULL) U(    0,      :,      :) = 0.0_sp
    ! If(right  .eq. MPI_PROC_NULL) U(nl(1),      :,      :) = 0.0_sp
    ! If(front  .eq. MPI_PROC_NULL) U(    :,      0,      :) = - U(:,1,:)
    ! If(back   .eq. MPI_PROC_NULL) U(    :,nl(2)+1,      :) = 1.0_sp
    ! If(bottom .eq. MPI_PROC_NULL) U(    :,      :,      0) = - U(:,:,1)
    ! If(top    .eq. MPI_PROC_NULL) U(    :,      :,nl(3)+1) = - U(:,:,nl(3))

    ! If(left   .eq. MPI_PROC_NULL) V(      0,    :,      :) = - V(1,:,:)
    ! If(right  .eq. MPI_PROC_NULL) V(nl(1)+1,    :,      :) = - V(nl(1),:,:)
    ! If(front  .eq. MPI_PROC_NULL) V(      :,    0,      :) = 0.0_sp
    ! If(back   .eq. MPI_PROC_NULL) V(      :,nl(2),      :) = 0.0_sp
    ! If(bottom .eq. MPI_PROC_NULL) V(      :,    :,      0) = - V(:,:,1)
    ! If(top    .eq. MPI_PROC_NULL) V(      :,    :,nl(3)+1) = - V(:,:,nl(3))

    ! If(left   .eq. MPI_PROC_NULL) W(      0,      :,    :) = - W(1,:,:)
    ! If(right  .eq. MPI_PROC_NULL) W(nl(1)+1,      :,    :) = - W(nl(1),:,:)
    ! If(front  .eq. MPI_PROC_NULL) W(      :,      0,    :) = - W(:,1,:)
    ! If(back   .eq. MPI_PROC_NULL) W(      :,nl(2)+1,    :) = - W(:,nl(2),:)
    ! If(bottom .eq. MPI_PROC_NULL) W(      :,      :,    0) = 0.0_sp
    ! If(top    .eq. MPI_PROC_NULL) W(      :,      :,nl(3)) = 0.0_sp

    ! U(    0,      :,      :) = 0.0_sp
    ! U(nl(1),      :,      :) = 0.0_sp
    ! U(    :,      0,      :) = - U(:,1,:)
    ! U(    :,nl(2)+1,      :) = 1.0_sp
    ! U(    :,      :,      0) = - U(:,:,1)
    ! U(    :,      :,nl(3)+1) = - U(:,:,nl(3))

    ! V(      0,    :,      :) = - V(1,:,:)
    ! V(nl(1)+1,    :,      :) = - V(nl(1),:,:)
    ! V(      :,    0,      :) = 0.0_sp
    ! V(      :,nl(2),      :) = 0.0_sp
    ! V(      :,    :,      0) = - V(:,:,1)
    ! V(      :,    :,nl(3)+1) = - V(:,:,nl(3))

    ! W(      0,      :,    :) = - W(1,:,:)
    ! W(nl(1)+1,      :,    :) = - W(nl(1),:,:)
    ! W(      :,      0,    :) = - W(:,1,:)
    ! W(      :,nl(2)+1,    :) = W(:,nl(2),:)
    ! W(      :,      :,    0) = 0.0_sp
    ! W(      :,      :,nl(3)) = 0.0_sp

  End Subroutine BC_UVW

  Subroutine BC_UVW2(U, V, W)
    Implicit None
    real(sp), intent(inout) , dimension(0:,0:,0:) :: U, V, W
    Call U_bc%SetBCS(U)
    Call V_bc%SetBCS(V)
    Call W_bc%SetBCS(W)
  End Subroutine BC_UVW2

  !==========================
  ! (3-2) Boudnary conditions for p
  !==========================
  Subroutine BC_P(P)
    Implicit None
    real(sp), intent(inout) , dimension(0:,0:,0:) :: P
    Call Phi_bc%SetBCS(P)
  End Subroutine BC_P

  !==========================
  ! (3-3) Boudnary conditions for rho and mu
  !==========================
  Subroutine BC_RHOMU(Rho, Mu)
    Implicit None
    real(sp), intent(inout) , dimension(0:,0:,0:) :: Rho, Mu
    Call Phi_bc%SetBCS(Rho)
    Call Phi_bc%SetBCS(Mu)
  End Subroutine BC_RHOMU

  !==========================
  ! Initialize the Navier-Stokes solver after common initialization
  ! Note:
  !   The bc types and values are initialized in basics.f90 for only VOF
  !   Here will overwrite some of the value.
  !==========================
  Subroutine InitNavierStokes(u, v, w, phi, p)
    Implicit None
    real(sp), intent(inout) , dimension(0:,0:,0:) :: p,u,v,w,phi
    Real(sp) :: iter_tolerance
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
    Flag = 1
    rkcoef   = 0.0_sp

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
    Call UpdtRhoMu(Phi)

    ! Initialize hypre
    Call u_bc%SetBCS(u)
    Call v_bc%SetBCS(v)
    Call w_bc%SetBCS(w)
    Call phi_bc%SetBCS(phi)
    Call phi_bc%SetBCS(p)
    jacobiflag = 0

  End Subroutine InitNavierStokes

  !==========================
  !  (4-1) Adjust dt
  !==========================
  Subroutine Adjustdt(U, V, W, cfl)
    Implicit None
    Real(sp), Intent(In), Dimension(0:,0:,0:) :: U, V, W
    Real(sp) :: uvw, uvw2, dlmin
    Real(sp), Intent(In) :: cfl
    Integer :: i, j, k
    uvw = 0.0_sp
    Do k = 1, nl(3)
      Do j = 1, nl(2)
        Do i = 1, nl(1)
          uvw = max(uvw,abs(u(i,j,k))+abs(v(i,j,k))+abs(w(i,j,k)))
        End Do
      End Do
    End Do

    dlmin = min(min(nl(1),nl(2)),nl(3))

    Call MPI_Reduce(uvw, uvw2, 1, MPI_REAL_SP, MPI_MAX, MPI_COMM_WORLD, 0, ierr)

    dt = cfl * sqrt(3.0_sp) * dlmin / uvw

    dt = min(dt0, dt)
    dt = max(dt, dt_min)

    Call MPI_Barrier(MPI_COMM_WORLD, ierr)

  End Subroutine Adjustdt

  !==========================
  !  (4-2) Jacobi iteration
  !==========================
  Subroutine Jacobi(P)
    Implicit None
    Real(sp), Intent(InOut) :: P(0:,0:,0:)
    Real(sp) :: PP(0:nl(1)+1,0:nl(2)+1,0:nl(3)+1)
    Integer :: i, j, k, ll

    Do ll = 1, 1000
      Do k = 1, nl(3)
        Do j = 1, nl(2)
          Do i = 1, nl(1)
            PP(i,j,k) = - ( DIV(i,j,k) / dt &
                &  - rhox(i-1,j,k) * P(i-1,j,k) - rhox(i,j,k) * P(i+1,j,k)    &
                &  - rhoy(i,j-1,k) * P(i,j-1,k) - rhoy(i,j,k) * P(i,j+1,k)  &
                &  - rhoz(i,j,k-1) * P(i,j,k-1) - rhoz(i,j,k) * P(i,j,k+1) )  &
                &  / ( rhox(i,j,k) + rhox(i-1,j,k) + rhoy(i,j,k) + rhoy(i,j-1,k) + rhoz(i,j,k) + rhoz(i,j,k-1) )
          End Do
        End Do
      End Do

      Do k = 1, nl(3)
        Do j = 1, nl(2)
          Do i = 1, nl(1)
            P(i,j,k) = PP(i,j,k)
          End Do
        End Do
      End Do
      Call P_bc%SetBCS(P)
    End Do

  End Subroutine Jacobi

  Subroutine Monitor(U,V,W)
    Implicit None
    Real(sp), Intent(In), Dimension(0:,0:,0:) :: U, V, W
    Integer :: i, j, k
    Call Divergence(U,V,W,Div)
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
