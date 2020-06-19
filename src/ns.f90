Module ModNavierStokes
  Use ModGlobal
  Use HypreStruct
  Implicit None
  Real(sp), Allocatable :: Rho(:,:,:), Rhox(:,:,:), Rhoy(:,:,:), Rhoz(:,:,:)
  Real(sp), Allocatable :: mu(:,:,:), P(:,:,:), Div(:,:,:)
  Integer, Allocatable :: flag(:,:,:)
  Real(sp), Allocatable :: bforce(:)
  Real(sp) :: rho_l, rho_g
  Real(sp) ::  mu_l, mu_g
  Integer :: rk_order
  Real(sp), Allocatable :: rkcoef(:,:)


Contains

  Subroutine InitNavierStokes
    Implicit None
    Real(sp) :: iter_tolerance
    Integer :: iter_max
    Character(80) :: file_name
    Character(80) :: input_name
    Integer :: i

    namelist /ns/ rho_l, rho_g, mu_l, mu_g, iter_tolerance, iter_max, rk_order

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
      Endif
      Call MPI_barrier(MPI_COMM_WORLD, ierr)
    End Do

    ! Allocate variables
    Allocate( Rho(0:nl(1)+1, 0:nl(2)+1, 0:nl(3)+1))
    Allocate(Rhox(0:nl(1)+1, 0:nl(2)+1, 0:nl(3)+1))
    Allocate(Rhoy(0:nl(1)+1, 0:nl(2)+1, 0:nl(3)+1))
    Allocate(Rhoz(0:nl(1)+1, 0:nl(2)+1, 0:nl(3)+1))
    Allocate(  Mu(0:nl(1)+1, 0:nl(2)+1, 0:nl(3)+1))
    Allocate(   P(0:nl(1)+1, 0:nl(2)+1, 0:nl(3)+1))
    Allocate(Flag(0:nl(1)+1, 0:nl(2)+1, 0:nl(3)+1))
    Allocate( Div(nl(1), nl(2), nl(3)))
    Allocate(rkcoef(2,rk_order))

    ! Assign RungeKutta
    If (rk_order .eq. 1) Then
      rkcoef(1,1) =  1.0_sp
      rkcoef(2,1) =  0.0_sp
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

    Call UpdtRhoMu

    ! Initialize hypre
    Call u_bc%SetBCS(u)
    Call v_bc%SetBCS(v)
    Call w_bc%SetBCS(w)
    Call phi_bc%SetBCS(phi)
    Call phi_bc%SetBCS(p)

  End Subroutine InitNavierStokes

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
    Call AdvDiff(u, v, w, dudtrk, dir=1)
    Call AdvDiff(v, w, u, dvdtrk, dir=2)
    Call AdvDiff(u, u, v, dwdtrk, dir=3)
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

  Subroutine AdvDiff(u1, u2, u3, dudt, dir)
    Implicit None
    Integer , Intent(in) :: dir
    Real(sp), Dimension(0:,0:,0:), Intent(in)  :: u1, u2, u3
    Real(sp), Dimension(:,:,:), Intent(out) :: dudt
    integer :: im,ip,jm,jp,km,kp,i,j,k
    Real(sp) :: uuip,uuim,uvjp,uvjm,uwkp,uwkm
    Real(sp) :: dudxp,dudxm,dudyp,dudym,dudzp,dudzm

    ! Calculate the advection-diffusion spatial term
    Do k=1,nl(3)
      kp = k + 1
      km = k - 1
      Do j=1,nl(2)
        jp = j + 1
        jm = j - 1
        Do i=1,nl(1)
          ip = i + 1
          im = i - 1
          uuip  = 0.25_sp * (u1(ip, j, k) + u1(i,j,k)) * (u1(ip,j ,k ) + u1(i,j ,k ))
          uuim  = 0.25_sp * (u1(im, j, k) + u1(i,j,k)) * (u1(im,j ,k ) + u1(i,j ,k ))
          uvjp  = 0.25_sp * (u1( i,jp, k) + u1(i,j,k)) * (u2(ip,j ,k ) + u2(i,j ,k ))
          uvjm  = 0.25_sp * (u1( i,jm, k) + u1(i,j,k)) * (u2(ip,jm,k ) + u2(i,jm,k ))
          uwkp  = 0.25_sp * (u1( i, j,kp) + u1(i,j,k)) * (u3(ip,j ,k ) + u3(i,j ,k ))
          uwkm  = 0.25_sp * (u1( i, j,km) + u1(i,j,k)) * (u3(ip,j ,km) + u3(i,j ,km))
          dudxp = mu(i,j,k) * (u1(ip, j, k) - u1( i, j, k)) / dl(1)
          dudxm = mu(i,j,k) * (u1( i, j, k) - u1(im, j, k)) / dl(1)
          dudyp = mu(i,j,k) * (u1( i,jp, k) - u1( i, j, k)) / dl(2)
          dudym = mu(i,j,k) * (u1( i, j, k) - u1( i,jm, k)) / dl(2)
          dudzp = mu(i,j,k) * (u1( i, j,kp) - u1( i, j, k)) / dl(3)
          dudzm = mu(i,j,k) * (u1( i, j, k) - u1( i, j,km)) / dl(3)
          !
          ! Momentum balance
          !
          dudt(i,j,k) = ((-uuip + uuim) + (dudxp-dudxm)) / rho(i,j,k) / dl(1) + &
                        ((-uvjp + uvjm) + (dudyp-dudym)) / rho(i,j,k) / dl(2) + &
                        ((-uwkp + uwkm) + (dudzp-dudzm)) / rho(i,j,k) / dl(3) + &
                        bforce(dir)
        Enddo
      Enddo
    Enddo

  End Subroutine AdvDiff

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

  Subroutine Divergence(up, vp, wp, div)
    Implicit None
    real(sp), intent(in) , dimension(0:,0:,0:) :: up, vp, wp
    real(sp), intent(out) , dimension(:,:,:) :: div
    Integer :: i, j, k
    Do k = 1, nl(3)
      Do j = 1, nl(2)
        Do i = 1, nl(1)
          Div(i,j,k) = &
              ( up(i+1,  j,  k) - up(i,j,k) ) / dl(1) + &
              ( vp(  i,j+1,  k) - vp(i,j,k) ) / dl(2) + &
              ( wp(  i,  j,k+1) - wp(i,j,k) ) / dl(3)
        End Do
      End Do
    End Do
  End Subroutine Divergence

  Subroutine Projection(P, u, v, w, up, vp, wp)
    Implicit None
    real(sp), intent(in) , dimension(0:,0:,0:) :: p,up,vp,wp
    real(sp), intent(out), dimension(0:,0:,0:) :: u,v,w
    Integer :: i, j, k

    do k=1,nl(3)
      do j=1,nl(2)
        do i=1,nl(1)
          u(i,j,k) = up(i,j,k) - (p(i+1,  j,  k) - p(i,j,k)) * dt / dl(1)
          v(i,j,k) = vp(i,j,k) - (p(  i,j+1,  k) - p(i,j,k)) * dt / dl(2)
          w(i,j,k) = wp(i,j,k) - (p(  i,  j,k+1) - p(i,j,k)) * dt / dl(3)
        enddo
      enddo
    enddo

  End Subroutine Projection

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
          Rhoy(i,j,k) = 0.5_sp * ( 1.0_sp / Rho(i,j,k) + 1.0_sp / Rho(i+1,j,k) )
          Rhoz(i,j,k) = 0.5_sp * ( 1.0_sp / Rho(i,j,k) + 1.0_sp / Rho(i+1,j,k) )
        EndDo
      EndDo
    EndDo
  End Subroutine UpdtRhoMu

  Subroutine BC_UVW
    Implicit None
    Call U_bc%SetBCS(U)
    Call V_bc%SetBCS(V)
    Call W_bc%SetBCS(W)
  End Subroutine BC_UVW

  Subroutine BC_P
    Implicit None
    Call Phi_bc%SetBCS(P)
  End Subroutine BC_P

  Subroutine BC_RHOMU
    Implicit None
    Call Phi_bc%SetBCS(RHO)
    Call Phi_bc%SetBCS(MU)
  End Subroutine BC_RHOMU

End Module ModNavierStokes
