Module ModNavierStokes
  Use ModGlobal
  Use ModTools
  Use HypreStruct
  Implicit None
  Real(sp), Allocatable :: Rho(:,:,:), Rhox(:,:,:), Rhoy(:,:,:), Rhoz(:,:,:)
  Real(sp), Allocatable :: mu(:,:,:), Div(:,:,:)
  Integer, Allocatable :: flag(:,:,:)
  Real(sp), Allocatable :: bforce(:)
  Real(sp) :: rho_l, rho_g
  Real(sp) ::  mu_l, mu_g
  Integer :: rk_order
  Real(sp), Allocatable :: rkcoef(:,:)


Contains

  Subroutine TwoPhaseFlow
    Implicit None
    real(sp), dimension(nl(1),nl(2),nl(3))    :: dudtrko,dvdtrko,dwdtrko
    real(sp), dimension(0:nl(1)+1,0:nl(2)+1,0:nl(3)+1)    :: up,vp,wp
    Integer :: irk

    Do irk=1,rk_order

      call RungeKuttaMomentum(rkcoef(:,irk), u, v, w, p, dudtrko,dvdtrko,dwdtrko, up, vp, wp)

      Call BC_UVW

      Call Divergence(u, v, w, Div)

      Call Hypre_Poisson(P, Div, Rhox, Rhoy, Rhoz, flag)

      Call BC_P

      Call Projection(P, u, v, w, up, vp, wp)

      Call BC_UVW
    End Do

    ! Call MOF
    Call UpdtRhoMu

  End Subroutine TwoPhaseFlow


  !==========================
  ! Runge Kutta time integration
  !   Key steps:
  ! Depending on the actual variable of u1, u2, u3,
  ! it will solve the corresponding momentum equation
  !    (1) Solve u: u1=u, u2=v, u3=w
  !    (2) Solve v: u1=v, u2=w, u3=u
  !    (3) Solve w: u1=w, u2=u, u3=v
  ! 1D view for u
  !    d(u)/dt_i =
  !     ( (uu)_i+1/2 - (uu)_i-1/2 + ( (mu*dudx)_i+1/2 - (mu*dudx)_i-1/2 ) / rho ) / dx
  !     ( (uv)_j+1/2 - (uv)_j-1/2 + ( (mu*dudy)_j+1/2 - (mu*dudy)_j-1/2 ) / rho ) / dy
  !     ( (uw)_k+1/2 - (uw)_k-1/2 + ( (mu*dudz)_k+1/2 - (mu*dudz)_k-1/2 ) / rho ) / dz
  !==========================
  Subroutine RungeKuttaMomentum(rkpar, u, v, w, p, dudtrko, dvdtrko, dwdtrko, up, vp, wp)
    Implicit None
    real(sp), intent(in), dimension(0:,0:,0:) :: u,v,w, p
    real(sp), intent(inout), dimension(:,:,:) :: dudtrko,dvdtrko,dwdtrko
    real(sp), intent(in), dimension(2) :: rkpar
    real(sp), dimension(nl(1),nl(2),nl(3)) :: dudtrk,dvdtrk,dwdtrk
    real(sp), intent(out), dimension(0:,0:,0:)    :: up,vp,wp
    real(sp) :: factor1,factor2,factor12
    Integer :: i, j, k

    factor1 = rkpar(1)*dt
    factor2 = rkpar(2)*dt
    factor12 = factor1 + factor2

    ! Calculate \partial P \partial x
    Call PartialP(P, dudtrk, dir=1)
    Call PartialP(P, dvdtrk, dir=2)
    Call PartialP(P, dwdtrk, dir=3)
    do k=1,nl(3)
      do j=1,nl(2)
        do i=1,nl(1)
          up(i,j,k) = u(i,j,k) + factor12*dudtrk(i,j,k)
          vp(i,j,k) = v(i,j,k) + factor12*dvdtrk(i,j,k)
          wp(i,j,k) = w(i,j,k) + factor12*dwdtrk(i,j,k)
        enddo
      enddo
    enddo

    ! Calculate Advection, diffusion
    ! Call AdvDiff(u, v, w, dudtrk, dir=1)
    ! Call AdvDiff(v, w, u, dvdtrk, dir=2)
    ! Call AdvDiff(u, u, v, dwdtrk, dir=3)
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
  ! Calculate Advection and diffusion
  ! 1D view for u
  !    d(u)/dt_i+1/2 =
  !     ( (uu)_i+1 - (uu)_i + ( (mu*dudx)_i+1 - (mu*dudx)_i ) / rho_i+1/2 ) / dx
  !     ( (uv)_i+1 - (uv)_i + ( (mu*dudy)_i+1 - (mu*dudy)_i ) / rho_i+1/2 ) / dy
  !     ( (uw)_i+1 - (uw)_i + ( (mu*dudz)_i+1 - (mu*dudz)_i ) / rho_i+1/2 ) / dz
  !  where:
  !      (uu)_i = 0.25 * (u_i+1/2 + u_i-1/2)^2
  !      (uv)_i = 0.25 * (u_i+1/2 + u_i-1/2) * (v_j+1/2,v_j-1/2)
  !      (uw)_i = 0.25 * (u_i+1/2 + u_i-1/2) * (w_k+1/2,w_k-1/2)
  !--------------------------
  ! Inputs:
  !    u1, u2, u3: velocity conponent
  !    rho_dir: the reciprocal of rho, depends on direction
  !    mu_c: mu at cell center
  ! Outputs:
  !    dudt: time integration of advection and diffusion
  !==========================
  Subroutine AdvDiffU(u, v, w, rho, mu_c, bforce, dudt)
    Implicit None
    Real(sp), Dimension(0:,0:,0:), Intent(in)  :: u, v, w
    Real(sp), Dimension(0:,0:,0:), Intent(in)  :: rho, mu_c
    Real(sp), Intent(in)  :: bforce
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
          dudxp = mu(i,j,k) * (u(i+1,  j,  k) - u(  i,  j,  k)) / dl(1)
          dudxm = mu(i,j,k) * (u(  i,  j,  k) - u(i-1,  j,  k)) / dl(1)
          dudyp = mu(i,j,k) * (u(  i,j+1,  k) - u(  i,  j,  k)) / dl(2)
          dudym = mu(i,j,k) * (u(  i,  j,  k) - u(  i,j-1,  k)) / dl(2)
          dudzp = mu(i,j,k) * (u(  i,  j,k+1) - u(  i,  j,  k)) / dl(3)
          dudzm = mu(i,j,k) * (u(  i,  j,  k) - u(  i,  j,k-1)) / dl(3)
          ! Momentum balance
          dudt(i,j,k) = &
              ((-uuip + uuim) + (dudxp-dudxm)) / rho(i,j,k) / dl(1) + &
              ((-uvjp + uvjm) + (dudyp-dudym)) / rho(i,j,k) / dl(2) + &
              ((-uwkp + uwkm) + (dudzp-dudzm)) / rho(i,j,k) / dl(3) + &
              bforce
        Enddo
      Enddo
    Enddo

  End Subroutine AdvDiffU

  !==========================
  ! Calculate divergence by
  ! 1D view
  ! du/dt_i = - (p_i+1/2 - p_i-1/2) / dx
  !--------------------------
  ! Inputs:
  !    dir: direction
  !    p: pressure
  ! Outputs:
  !    dudt: time integration of pressure
  !==========================
  Subroutine PartialP(P, dudt, dir)
    Integer, intent(in) :: dir
    real(sp), dimension(0:,0:,0:), intent(in) :: p
    real(sp), dimension(:,:,:), intent(out) :: dudt
    Integer :: i, j, k
    Integer :: ip, jp, kp

    If (dir .eq. 1) Then
      ip =-1; jp = 0; kp = 0
    ElseIf (dir .eq. 2) Then
      ip = 0; jp =-1; kp = 0
    Else
      ip = 0; jp = 0; kp =-1
    End If

    do k=1, nl(3)
      do j=1, nl(2)
        do i=1, nl(1)
          dudt(i,j,k) = - (p(i,j,k) - p(i+ip,j+jp,k+jp)) / dl(dir)
        enddo
      enddo
    enddo

  End Subroutine PartialP

  !==========================
  ! Calculate divergence by
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
              ( up(i,j,k) - up(i-1,  j,  k) ) / dl(1) + &
              ( vp(i,j,k) - vp(  i,j-1,  k) ) / dl(2) + &
              ( wp(i,j,k) - wp(  i,  j,k-1) ) / dl(3)
        End Do
      End Do
    End Do
  End Subroutine Divergence

  !==========================
  ! Projection procedure
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
          u(i,j,k) = up(i,j,k) - (p(i,j,k) - p(i-1,  j,  k)) * dt / dl(1)
          v(i,j,k) = vp(i,j,k) - (p(i,j,k) - p(  i,j-1,  k)) * dt / dl(2)
          w(i,j,k) = wp(i,j,k) - (p(i,j,k) - p(  i,  j,k-1)) * dt / dl(3)
        enddo
      enddo
    enddo

  End Subroutine Projection

  !==========================
  ! Boudnary conditions for u, v, w
  !==========================
  Subroutine UpdtRhoMu
    Implicit None
    Integer :: i,j,k
    Do k=1,nl(3)
      Do j=1,nl(2)
        Do i=1,nl(1)
          Rho(i,j,k) = Phi(i,j,k) * rho_l + Phi(i,j,k) * rho_g
          Mu(i,j,k) = Phi(i,j,k) *  mu_l + Phi(i,j,k) *  mu_g
        enddo
      enddo
    enddo

    Call Phi_bc%SetBCS(Phi)
    Call BC_RHOMU

    Do k=1,nl(3)
      Do j=1,nl(2)
        Do i=1,nl(1)
          Rhox(i,j,k) = 0.5_sp * ( 1.0_sp / Rho(i,j,k) + 1.0_sp / Rho(i+1,j,k) )
          Rhoy(i,j,k) = 0.5_sp * ( 1.0_sp / Rho(i,j,k) + 1.0_sp / Rho(i,j+1,k) )
          Rhoz(i,j,k) = 0.5_sp * ( 1.0_sp / Rho(i,j,k) + 1.0_sp / Rho(i,j,k+1) )
        EndDo
      EndDo
    EndDo
  End Subroutine UpdtRhoMu

  !==========================
  ! Boudnary conditions for u, v, w
  !==========================
  Subroutine BC_UVW
    Implicit None
    Call U_bc%SetBCS(U)
    Call V_bc%SetBCS(V)
    Call W_bc%SetBCS(W)
  End Subroutine BC_UVW

  !==========================
  ! Boudnary conditions for p
  !==========================
  Subroutine BC_P
    Implicit None
    Call Phi_bc%SetBCS(P)
  End Subroutine BC_P

  !==========================
  ! Boudnary conditions for rho and mu
  !==========================
  Subroutine BC_RHOMU
    Implicit None
    Call Phi_bc%SetBCS(RHO)
    Call Phi_bc%SetBCS(MU)
  End Subroutine BC_RHOMU

  !==========================
  ! Initialize the Navier-Stokes solver after common initialization
  ! Note:
  !   The bc types and values are initialized in basics.f90 for only VOF
  !   Here will overwrite some of the value.
  !==========================
  Subroutine InitNavierStokes
    Implicit None
    Real(sp) :: iter_tolerance
    Integer :: iter_max
    Character(80) :: file_name
    Character(80) :: input_name
    Real(sp), Dimension(3) :: phi_bc_value_lo, phi_bc_value_hi
    Real(sp), Dimension(3) ::   p_bc_value_lo,   p_bc_value_hi
    Real(sp), Dimension(3) ::   u_bc_value_lo,   u_bc_value_hi
    Real(sp), Dimension(3) ::   v_bc_value_lo,   v_bc_value_hi
    Real(sp), Dimension(3) ::   w_bc_value_lo,   w_bc_value_hi
    Integer,  Dimension(3) :: phi_bc_types_lo, phi_bc_types_hi
    Integer,  Dimension(3) ::   p_bc_types_lo,   p_bc_types_hi
    Integer,  Dimension(3) ::   u_bc_types_lo,   u_bc_types_hi
    Integer,  Dimension(3) ::   v_bc_types_lo,   v_bc_types_hi
    Integer,  Dimension(3) ::   w_bc_types_lo,   w_bc_types_hi
    Integer :: i

    !  variables for physics
    namelist /ns/ rho_l, rho_g, mu_l, mu_g, iter_tolerance, iter_max, rk_order
    !  variables for boundary conditions
    namelist /ns_bc/ &
        phi_bc_value_lo, phi_bc_value_hi, &
        p_bc_value_lo,   p_bc_value_hi, &
        u_bc_value_lo,   u_bc_value_hi, &
        v_bc_value_lo,   v_bc_value_hi, &
        w_bc_value_lo,   w_bc_value_hi, &
        phi_bc_types_lo, phi_bc_types_hi, &
        p_bc_types_lo,   p_bc_types_hi, &
        u_bc_types_lo,   u_bc_types_hi, &
        v_bc_types_lo,   v_bc_types_hi, &
        w_bc_types_lo,   w_bc_types_hi

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
        Read(10, nml = ns)
        Close(10)
        if (myid .eq. i-1) then
          Open(10, file=file_name)
          Read(10, nml = ns_bc)
          Close(10)
        Endif
      Endif
      Call MPI_barrier(MPI_COMM_WORLD, ierr)
    End Do

    ! Allocate variables
    Allocate( Rho(0:nl(1)+1, 0:nl(2)+1, 0:nl(3)+1))
    Allocate(Rhox(0:nl(1)+1, 0:nl(2)+1, 0:nl(3)+1))
    Allocate(Rhoy(0:nl(1)+1, 0:nl(2)+1, 0:nl(3)+1))
    Allocate(Rhoz(0:nl(1)+1, 0:nl(2)+1, 0:nl(3)+1))
    Allocate(  Mu(0:nl(1)+1, 0:nl(2)+1, 0:nl(3)+1))
    Allocate(Flag(0:nl(1)+1, 0:nl(2)+1, 0:nl(3)+1))
    Allocate( Div(nl(1), nl(2), nl(3)))
    Allocate(rkcoef(2,rk_order))

    ! Set boundary conditions for phi
    phi_bc%bound_type(:,1)  = phi_bc_types_lo
    phi_bc%bound_type(:,2)  = phi_bc_types_hi
    phi_bc%bound_value(:,1) = phi_bc_value_lo
    phi_bc%bound_value(:,2) = phi_bc_value_hi
    ! Set boundary conditions for p
    p_bc%bound_type(:,1)    = p_bc_types_lo
    p_bc%bound_type(:,2)    = p_bc_types_hi
    p_bc%bound_value(:,1)   = p_bc_value_lo
    p_bc%bound_value(:,2)   = p_bc_value_hi
    ! Set boundary conditions for u
    u_bc%bound_type(:,1)    = u_bc_types_lo
    u_bc%bound_type(:,2)    = u_bc_types_hi
    u_bc%bound_value(:,1)   = u_bc_value_lo
    u_bc%bound_value(:,2)   = u_bc_value_hi
    ! Set boundary conditions for v
    v_bc%bound_type(:,1)    = v_bc_types_lo
    v_bc%bound_type(:,2)    = v_bc_types_hi
    v_bc%bound_value(:,1)   = v_bc_value_lo
    v_bc%bound_value(:,2)   = v_bc_value_hi
    ! Set boundary conditions for w
    w_bc%bound_type(:,1)    = w_bc_types_lo
    w_bc%bound_type(:,2)    = w_bc_types_hi
    w_bc%bound_value(:,1)   = w_bc_value_lo
    w_bc%bound_value(:,2)   = w_bc_value_hi

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
    End If

    ! Initialize hypre
    Call Hypre_Initialize(iter_tolerance, iter_max)

    ! Initialize the Rho and Mu field
    Call UpdtRhoMu

    ! Initialize hypre
    Call u_bc%SetBCS(u)
    Call v_bc%SetBCS(v)
    Call w_bc%SetBCS(w)
    Call phi_bc%SetBCS(phi)
    Call phi_bc%SetBCS(p)

  End Subroutine InitNavierStokes

End Module ModNavierStokes
