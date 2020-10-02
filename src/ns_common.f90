#include "param.h"

!==========================
! (1-1) Calculate advection and diffusion
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
Subroutine AdvDiffU_CDS(u, v, w, dl, nl, mu, muxy, muxz, rhox, dudt)
  Use ModGlobal, only : sp
  Implicit None
  Real(sp), Intent(In) :: dl(3)
  Integer, Intent(In) :: nl(3)
  Real(sp), Dimension(0:nl(1)+1,0:nl(2)+1,0:nl(3)+1), Intent(in)  :: u, v, w
  Real(sp), Dimension(0:nl(1)+1,0:nl(2)+1,0:nl(3)+1), Intent(in)  :: mu, muxy, muxz, rhox
  Real(sp), Dimension(nl(1),nl(2),nl(3)), Intent(out) :: dudt
  integer :: i,j,k
  Real(sp),Dimension(0:nl(1)+1,0:nl(2)+1,0:nl(3)+1) :: uu, uv, uw, dudx, dudy, dudz

  ! Calculate the advection-diffusion spatial term
  Do k=0,nl(3)
    Do j=0,nl(2)
      Do i=1,nl(1) + 1
        uu(i,j,k) = 0.25_sp * (u(i-1,  j,  k) + u(i,j,k)) * (u(i-1,j,k)+ u(i,j,k))
        uv(i,j,k) = 0.25_sp * (u(  i,j+1,  k) + u(i,j,k)) * (v(i-1,j,k)+ v(i,j,k))
        uw(i,j,k) = 0.25_sp * (u(  i,  j,k+1) + u(i,j,k)) * (w(i-1,j,k)+ w(i,j,k))
        dudx(i,j,k) =   mu(i,j,k) * ( u(i,  j,  k) - u(i-1,j,k) ) / dl(1)
        dudy(i,j,k) = muxz(i,j,k) * ( u(i,j+1,  k) - u(  i,j,k) ) / dl(2)
        dudz(i,j,k) = muxy(i,j,k) * ( u(i,  j,k+1) - u(  i,j,k) ) / dl(3)
      Enddo
    Enddo
  Enddo

  Do k=1,nl(3)
    Do j=1,nl(2)
      Do i=1,nl(1)
        dudt(i,j,k) = &
            ((-uu(i+1,j,k) + uu(i,  j,  k)) + (dudx(i+1,j,k)-dudx(i,  j,  k)) * rhox(i,j,k) ) / dl(1) + &
            ((-uv(  i,j,k) + uv(i,j-1,  k)) + (dudy(  i,j,k)-dudy(i,j-1,  k)) * rhox(i,j,k) ) / dl(2) + &
            ((-uw(  i,j,k) + uw(i,  j,k-1)) + (dudz(  i,j,k)-dudz(i,  j,k-1)) * rhox(i,j,k) ) / dl(3)
      Enddo
    Enddo
  Enddo

End Subroutine AdvDiffU_CDS

!==========================
! (1-2) Calculate advection and diffusion for V
! similar with  AdvDiffU_CDS
!--------------------------
! Inputs:
!   U, V, W: velocity conponents
! Outpurts
!   dvdt: dv/dt
!==========================
Subroutine AdvDiffV_CDS(u, v, w, dl, nl, mu, muxy, muyz, rhoy, dvdt)
  Use ModGlobal, only : sp
  Implicit None
  Real(sp), Intent(In) :: dl(3)
  Integer, Intent(In) :: nl(3)
  Real(sp), Dimension(0:nl(1)+1,0:nl(2)+1,0:nl(3)+1), Intent(in)  :: u, v, w
  Real(sp), Dimension(0:nl(1)+1,0:nl(2)+1,0:nl(3)+1), Intent(in)  :: mu, muyz, muxy, rhoy
  Real(sp), Dimension(nl(1),nl(2),nl(3)), Intent(out) :: dvdt
  integer :: i,j,k
  Real(sp),Dimension(0:nl(1)+1,0:nl(2)+1,0:nl(3)+1) :: uv, vv, vw, dvdx, dvdy, dvdz

  ! Calculate the advection-diffusion spatial term
  Do k=0,nl(3)
    Do j=1,nl(2) + 1
      Do i=0,nl(1) 
        uv(i,j,k) = 0.25_sp * (v(i+1,  j,  k) + v(i,j,k)) * (u(i,j-1,k)+ u(i,j,k))
        vv(i,j,k) = 0.25_sp * (v(  i,j-1,  k) + v(i,j,k)) * (v(i,j-1,k)+ v(i,j,k))
        vw(i,j,k) = 0.25_sp * (v(  i,  j,k+1) + v(i,j,k)) * (w(i,j-1,k)+ w(i,j,k))
        dvdx(i,j,k) = muyz(i,j,k) * ( v(i+1,j,  k) - v(i,  j,k) ) / dl(1)
        dvdy(i,j,k) =   mu(i,j,k) * ( v(  i,j,  k) - v(i,j-1,k) ) / dl(2)
        dvdz(i,j,k) = muxy(i,j,k) * ( v(  i,j,k+1) - v(i,  j,k) ) / dl(3)
      Enddo
    Enddo
  Enddo

  Do k=1,nl(3)
    Do j=1,nl(2)
      Do i=1,nl(1)
        dvdt(i,j,k) = &
             ((-uv(i,  j,k) + uv(i-1,j,  k)) + (dvdx(i,  j,k)-dvdx(i-1,j,  k)) * rhoy(i,j,k) ) / dl(1) + &
             ((-vv(i,j+1,k) + vv(  i,j,  k)) + (dvdy(i,j+1,k)-dvdy(  i,j,  k)) * rhoy(i,j,k) ) / dl(2) + &
             ((-vw(i,  j,k) + vw(  i,j,k-1)) + (dvdz(i,  j,k)-dvdz(  i,j,k-1)) * rhoy(i,j,k) ) / dl(3)
      Enddo
    Enddo
  Enddo

End Subroutine AdvDiffV_CDS

!==========================
! (1-3) Calculate advection and diffusion for V
! similar with  AdvDiffU_CDS
!--------------------------
! Inputs:
!   U, V, W: velocity conponents
! Outpurts
!   dvdt: dv/dt
!==========================
Subroutine AdvDiffW_CDS(u, v, w, dl, nl, mu, muxz, muyz, rhoz, dwdt)
  Use ModGlobal, only : sp
  Implicit None
  Real(sp), Intent(In) :: dl(3)
  Integer, Intent(In) :: nl(3)
  Real(sp), Dimension(0:nl(1)+1,0:nl(2)+1,0:nl(3)+1), Intent(in)  :: u, v, w
  Real(sp), Dimension(0:nl(1)+1,0:nl(2)+1,0:nl(3)+1), Intent(in)  :: mu, muxz, muyz, rhoz
  Real(sp), Dimension(nl(1),nl(2),nl(3)), Intent(out) :: dwdt
  integer :: i,j,k
  Real(sp),Dimension(0:nl(1)+1,0:nl(2)+1,0:nl(3)+1) :: uw, vw, ww, dwdx, dwdy, dwdz

  ! Calculate the advection-diffusion spatial term
  Do k=1,nl(3) + 1
    Do j=0,nl(2)
      Do i=0,nl(1)
        uw(i,j,k) = 0.25_sp * (u(i,j,k-1) + u(i,j,k)) * (w(i+1,  j,  k)+ w(i,j,k))
        vw(i,j,k) = 0.25_sp * (v(i,j,k-1) + v(i,j,k)) * (w(  i,j+1,  k)+ w(i,j,k))
        ww(i,j,k) = 0.25_sp * (w(i,j,k-1) + w(i,j,k)) * (w(  i,  j,k-1)+ w(i,j,k))
        dwdx(i,j,k) = muyz(i,j,k) * ( w(i+1,  j,k) - w(i,j,  k) ) / dl(1)
        dwdy(i,j,k) = muxz(i,j,k) * ( w(  i,j+1,k) - w(i,j,  k) ) / dl(2)
        dwdz(i,j,k) =   mu(i,j,k) * ( w(  i,  j,k) - w(i,j,k-1) ) / dl(3)
      Enddo
    Enddo
  Enddo

  Do k=1,nl(3)
    Do j=1,nl(2)
      Do i=1,nl(1)
        dwdt(i,j,k) = &
            ((-uw(i,j,  k) + uw(i-1,  j,k)) + (dwdx(i,j,  k)-dwdx(i-1,  j,k)) * rhoz(i,j,k) ) / dl(1) + &
            ((-vw(i,j,  k) + vw(  i,j-1,k)) + (dwdy(i,j,  k)-dwdy(  i,j-1,k)) * rhoz(i,j,k) ) / dl(2) + &
            ((-ww(i,j,k+1) + ww(  i,  j,k)) + (dwdz(i,j,k+1)-dwdz(  i,  j,k)) * rhoz(i,j,k) ) / dl(3)
      Enddo
    Enddo
  Enddo

End Subroutine AdvDiffW_CDS

!==========================
! (1-4) CIPCSLR scheme
!--------------------------
! Inputs:
!   U, V, W: velocity conponents
! Outpurts
!   dvdt: dv/dt
!==========================
Subroutine CIPCSLR_1D(Via, Sia, Usx, nl, dl, dt1, dir)
  Use ModGlobal, only : sp
  Implicit None
  Integer, Intent(In) :: nl(3)
  Real(sp), Intent(In) :: dl(3)
  Integer,  Intent(In) :: dir
  Real(sp), Intent(In) :: dt1
  Real(sp), Intent(InOut) :: Via(0:nl(1)+1, 0:nl(2)+1, 0:nl(3)+1)
  Real(sp), Intent(InOut) :: Sia(0:nl(1)+1, 0:nl(2)+1, 0:nl(3)+1)
  Real(sp), Intent(In) :: Usx(0:nl(1)+1, 0:nl(2)+1, 0:nl(3)+1)

  Real(sp), Dimension(0:nl(1)+1, 0:nl(2)+1, 0:nl(3)+1) :: Flux, Fn_sx

  Real(8) :: dd, xx
  Real(8) :: a1, b1, c1
  Real(8) :: ep = 1.d-10

  Integer :: i, j, k
  Integer :: ii, jj, kk
  Integer :: ip, jp, kp
  Integer :: iv, jv, kv

  if (dir == 1) then
    ii = 1; jj = 0; kk = 0
  else if (dir == 1) then
    ii = 0; jj = 1; kk = 0
  else if (dir == 3) then
    ii = 0; jj = 0; kk = 1
  end if

  Do k = 1, nl(3)
    Do j = 1, nl(2)
      Do i = 1, nl(1)
        xx = - Usx(i,j,k) * dt1
        if ( xx .lt. 0.0_sp) then
          ip = i - ii
          jp = j - jj
          kp = k - kk
          iv = i 
          jv = j 
          kv = k 
          dd = - dl(dir)
        else
          ip = i + ii
          jp = j + jj
          kp = k + kk
          iv = i + ii
          jv = j + jj
          kv = k + kk
          dd = dl(dir)
        end if

        a1  = ( dAbs( ( Via(iv,jv,kv) - Sia(i,j,k) + ep ) / ( Sia(ip,jp,kp) - Via(iv,jv,kv) + ep ) )-1.0_sp ) / dd
        b1  = a1 * Via(iv,jv,kv) + ( Via(iv,jv,kv) - Sia(i,j,k) ) / dd
        c1  = 1.d0 + a1 * xx            
        Fn_sx(i,j,k) =  ( ( a1 * b1 * xx + 2.d0 * b1 ) * xx + Sia(i,j,k) ) / c1 / c1
        Flux(i,j,k) =  -xx * ( b1 * xx + Sia(i,j,k) ) / c1
      End Do
    End Do
  End Do

  Call BC_Flux(flux, Usx, dt1,dir)

  Do k = 1, nl(3)
    Do j = 1, nl(2)
      Do i = 1, nl(1)
        ! Update F_v by flux
        Via(i,j,k)  =  Via(i,j,k) - ( Flux(i,j,k) - Flux(i-ii,j-jj,k-kk) ) / Dl(dir)
        Sia(i,j,k) =  Fn_sx(i,j,k)
      End Do
    End Do
  End Do

End Subroutine CIPCSLR_1D

Subroutine BC_Flux(flux, U, dt, dir)
  Use mpi
  Use ModGlobal, only : sp
  Use ModGlobal, only : nl
  Use ModGlobal, only : left, right, front, back, bottom, top
  Use ModGlobal, only : updthalo
  Implicit None
  Real(sp), Intent(InOut) :: flux(0:nl(1)+1,0:nl(2)+1,0:nl(3)+1)
  Real(sp), Intent(In) :: u(0:nl(1)+1,0:nl(2)+1,0:nl(3)+1)
  Real(sp), Intent(In) :: dt
  Integer, Intent(In) :: dir

  Call updthalo(nl,1,flux)
  Call updthalo(nl,2,flux)
  Call updthalo(nl,3,flux)

  If (dir == 1) Then
    If(left   .eq. MPI_PROC_NULL) flux(    0,0:nl(2)+1,0:nl(3)+1) = U(    0,0:nl(2)+1,0:nl(3)+1) * dt
    If(right  .eq. MPI_PROC_NULL) flux(nl(1),0:nl(2)+1,0:nl(3)+1) = U(nl(1),0:nl(2)+1,0:nl(3)+1) * dt
  Else If (dir == 2) Then
    If(front  .eq. MPI_PROC_NULL) flux(0:nl(1)+1,    0,0:nl(3)+1) = U(0:nl(1)+1,    0,0:nl(3)+1) * dt
    If(back   .eq. MPI_PROC_NULL) flux(0:nl(1)+1,nl(2),0:nl(3)+1) = U(0:nl(1)+1,nl(2),0:nl(3)+1) * dt
  Else If (dir == 3) Then
    If(bottom .eq. MPI_PROC_NULL) flux(0:nl(1)+1,0:nl(2)+1,    0) = U(0:nl(1)+1,0:nl(2)+1,    0) * dt
    If(top    .eq. MPI_PROC_NULL) flux(0:nl(1)+1,0:nl(2)+1,nl(3)) = U(0:nl(1)+1,0:nl(2)+1,nl(3)) * dt
  End If

End Subroutine BC_Flux

Subroutine Central_Difference_RHS(RHS_X,RHS_Y,RHS_Z,&
                                  U_v, V_V, W_v, &
                                  mu, rho, nl, dl)
  Use ModGlobal, only : sp
  Implicit None
  Integer, Intent(In) :: nl(3)
  Real(sp),  Intent(In) :: dl(3)
  Real(sp), Intent(In) :: U_v(0:nl(1)+1,0:nl(2)+1,0:nl(3)+1)
  Real(sp), Intent(In) :: V_v(0:nl(1)+1,0:nl(2)+1,0:nl(3)+1)
  Real(sp), Intent(In) :: W_v(0:nl(1)+1,0:nl(2)+1,0:nl(3)+1)
  Real(sp), Intent(Out) :: RHS_X(nl(1),nl(2),nl(3))
  Real(sp), Intent(Out) :: RHS_Y(nl(1),nl(2),nl(3))
  Real(sp), Intent(Out) :: RHS_Z(nl(1),nl(2),nl(3))
  Real(sp), Intent(In) :: Rho(0:nl(1)+1,0:nl(2)+1,0:nl(3)+1)
  Real(sp), Intent(In) :: Mu(0:nl(1)+1,0:nl(2)+1,0:nl(3)+1)

  Real(8) :: mu_L, mu_R
  Real(8) :: mu_U, mu_D
  Real(8) :: mu_F, mu_B

  Real(8) :: du_dx_sx_R, du_dx_sx_L
  Real(8) :: du_dy_sy_U, du_dy_sy_D
  Real(8) :: du_dz_sz_B, du_dz_sz_F

  Real(8) :: dv_dx_sx_R, dv_dx_sx_L
  Real(8) :: dv_dy_sy_U, dv_dy_sy_D
  Real(8) :: dv_dz_sz_B, dv_dz_sz_F

  Real(8) :: dw_dx_sx_R, dw_dx_sx_L
  Real(8) :: dw_dy_sy_U, dw_dy_sy_D
  Real(8) :: dw_dz_sz_B, dw_dz_sz_F

  Real(8) :: Mu_PPu, Mu_PPv, Mu_ppw
  Integer :: i , j , k

  RHS_X = 0.d0
  RHS_Y = 0.d0
  RHS_Z = 0.d0


  Do k =  1, nl(3)
  Do j =  1, nl(2)
  Do i =  1, nl(1)

  ! Interpolation to six points around SIA_x , SIA_Y, SIA_z position
    mu_R =  0.5_sp * ( Mu(i,j,k) + Mu(i+1,j,k) )
    mu_L =  0.5_sp * ( Mu(i,j,k) + Mu(i-1,j,k) )
    mu_U =  0.5_sp * ( Mu(i,j,k) + Mu(i,j+1,k) )
    mu_D =  0.5_sp * ( Mu(i,j,k) + Mu(i,j-1,k) )
    mu_B =  0.5_sp * ( Mu(i,j,k) + Mu(i,j,k+1) )
    mu_F =  0.5_sp * ( Mu(i,j,k) + Mu(i,j,k-1) )

  ! Compute partial U at VIA
    du_dx_sx_R = ( U_v(i+1,  j,  k) - U_v(  i,  j,  k) ) /  dl(1)
    du_dx_sx_L = ( U_v(  i,  j,  k) - U_v(i-1,  j,  k) ) /  dl(1)
    du_dy_sy_U = ( U_v(  i,j+1,  k) - U_v(  i,  j,  k) ) /  dl(2)
    du_dy_sy_D = ( U_v(  i,  j,  k) - U_v(  i,j-1,  k) ) /  dl(2)
    du_dz_sz_B = ( U_v(  i,  j,k+1) - U_v(  i,  j,  k) ) /  dl(3)
    du_dz_sz_F = ( U_v(  i,  j,  k) - U_v(  i,  j,k-1) ) /  dl(3)

  ! Mu_PPu: Mu times partial partial U at VIA 
    Mu_PPu = ( mu_R * du_dx_sx_R - mu_L * du_dx_sx_L ) / dl(1) &
      &    + ( mu_U * du_dy_sy_U - mu_D * du_dy_sy_D ) / dl(2) &
      &    + ( mu_B * du_dz_sz_B - mu_F * du_dz_sz_F ) / dl(3)


  ! Compute partial V at VIA
    RHS_X(i,j,k) = Mu_PPu / Rho(i,j,k)

  ! Compute partial V at VIA
    dv_dx_sx_R = ( V_v(i+1,  j,  k) - V_v(  i,  j,  k) ) /  dl(1)
    dv_dx_sx_L = ( V_v(  i,  j,  k) - V_v(i-1,  j,  k) ) /  dl(1)
    dv_dy_sy_U = ( V_v(  i,j+1,  k) - V_v(  i,  j,  k) ) /  dl(2)
    dv_dy_sy_D = ( V_v(  i,  j,  k) - V_v(  i,j-1,  k) ) /  dl(2)
    dv_dz_sz_B = ( V_v(  i,  j,k+1) - V_v(  i,  j,  k) ) /  dl(3)
    dv_dz_sz_F = ( V_v(  i,  j,  k) - V_v(  i,  j,k-1) ) /  dl(3)

  ! Mu_PPu: Mu times partial partial V at VIA 
    Mu_PPv = ( mu_R * dv_dx_sx_R - mu_L * dv_dx_sx_L ) / dl(1) &
      &    + ( mu_U * dv_dy_sy_U - mu_D * dv_dy_sy_D ) / dl(2) &
      &    + ( mu_B * dv_dz_sz_B - mu_F * dv_dz_sz_F ) / dl(3)

  ! Compute partial U at VIA
    RHS_Y(i,j,k) = Mu_PPv / Rho(i,j,k)

    ! Compute partial W at VIA
    dw_dx_sx_R = ( W_v(i+1,  j,  k) - W_v(  i,  j,  k) ) /  dl(1)
    dw_dx_sx_L = ( W_v(  i,  j,  k) - W_v(i-1,  j,  k) ) /  dl(1)
    dw_dy_sy_U = ( W_v(  i,j+1,  k) - W_v(  i,  j,  k) ) /  dl(2)
    dw_dy_sy_D = ( W_v(  i,  j,  k) - W_v(  i,j-1,  k) ) /  dl(2)
    dw_dz_sz_B = ( W_v(  i,  j,k+1) - W_v(  i,  j,  k) ) /  dl(3)
    dw_dz_sz_F = ( W_v(  i,  j,  k) - W_v(  i,  j,k-1) ) /  dl(3)

  ! Mu_PPu: Mu times partial partial W at VIA 
    Mu_PPw = ( mu_R * dw_dx_sx_R - mu_L * dw_dx_sx_L ) / dl(1) &
      &    + ( mu_U * dw_dy_sy_U - mu_D * dw_dy_sy_D ) / dl(2) &
      &    + ( mu_B * dw_dz_sz_B - mu_F * dw_dz_sz_F ) / dl(3)

  ! Compute partial U at VIA
    RHS_Z(i,j,k) =  Mu_PPw / Rho(i,j,k)

  EndDo
  EndDo
  EndDo

End Subroutine Central_Difference_RHS


Subroutine Tec_SIA(SIA, VIA, VIA_n, nl, tec_m, dir)
  Use ModGlobal, only : sp
  Implicit None
  Integer, Intent(In) :: nl(3)
  Integer,  Intent(In) :: dir
  Real(sp), Intent(In) :: tec_m
  Real(sp), Intent(InOut) :: SIA(0:nl(1)+1,0:nl(2)+1,0:nl(3)+1)
  Real(sp), Intent(In)    :: VIA(0:nl(1)+1,0:nl(2)+1,0:nl(3)+1)
  Real(sp), Intent(In)    :: VIA_n(0:nl(1)+1,0:nl(2)+1,0:nl(3)+1)

  Integer :: i,j,k
  Integer :: ii, jj, kk

  if (dir == 1) then
    ii = 1; jj = 0; kk = 0
  else if (dir == 1) then
    ii = 0; jj = 1; kk = 0
  else if (dir == 3) then
    ii = 0; jj = 0; kk = 1
  end if

  Do k = 1, nl(3)
    Do j = 1, nl(2)
      Do i = 1, nl(1)
      SIA(i,j,k) = tec_m * ( SIA(i,j,k) + &
                   0.5_sp * ( VIA(i,j,k) - VIA_n(i,j,k) + VIA(i+i,j+jj,k+kk) - VIA_n(i+ii,j+jj,k+kk) )  ) &
                             + 0.5_sp * (1.d0 - tec_m) * (  VIA(i,j,k) + VIA(i+ii,j+jj,k+kk) ) 
      EndDo
    EndDo
  EndDo

End Subroutine Tec_SIA

Subroutine Tec_VIA(VIA, SIA, SIA_n, nl, dir)

  Use ModGlobal, only : sp
  Implicit None
  Integer, Intent(In) :: nl(3)
  Integer,  Intent(In) :: dir
  Real(sp), Intent(InOut) :: VIA(0:nl(1)+1,0:nl(2)+1,0:nl(3)+1)
  Real(sp), Intent(In)    :: SIA(0:nl(1)+1,0:nl(2)+1,0:nl(3)+1)
  Real(sp), Intent(In)    :: SIA_n(0:nl(1)+1,0:nl(2)+1,0:nl(3)+1)

  Integer :: i,j,k
  Integer :: ii, jj, kk

  if (dir == 1) then
    ii = 1; jj = 0; kk = 0
  else if (dir == 1) then
    ii = 0; jj = 1; kk = 0
  else if (dir == 3) then
    ii = 0; jj = 0; kk = 1
  end if

  Do k = 1, nl(3)+1
    Do j = 1, nl(2)+1
      Do i = 1, nl(1)+1
        VIA(i,j,k) =  VIA(i,j,k) + 0.5d0 * ( SIA(i,j,k) - SIA_n(i,j,k) + SIA(i-ii,j-jj,k-kk) - SIA_n(i-ii,j-jj,k-kk) )
      EndDo
    EndDo
  EndDo

End Subroutine Tec_VIA



!==========================
! (2-1) Calculate body force
!    only support uniformly body force, for example, gravity
!--------------------------
! Inputs & outputs:
!    dudt, dvdt, dwdt: du/dt, dv/dt, dwdt
!==========================
Subroutine BodyForce(dudt, dvdt, dwdt, nl, body_force)
  Use ModGlobal, only : sp
  Implicit None
  Integer :: nl(3)
  Real(sp), Dimension(1:nl(1),1:nl(2),1:nl(3)), Intent(inout)  :: dudt, dvdt,dwdt
  Real(sp), Dimension(3), Intent(in)  :: body_force
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
! (2-2) Calculate advance of partial P by
! 1D view
! u_i = up_i - (p_i+1/2 - p_i-1/2) * dt / dx
!--------------------------
! Inputs:
!    dir: direction
!    p: pressure
! Outputs:
!    dudt: time integration of pressure
!==========================
Subroutine DuDtPartialP(U1, V1, W1, Rhox, Rhoy, Rhoz, P, dt1, nl, dl, U2, V2, W2)
  Use ModGlobal, only : sp
  Real(sp), Intent(In) :: dt1
  Integer, Intent(In) :: nl(3)
  Real(sp), Intent(In) :: dl(3)
  Real(sp), Dimension(0:nl(1)+1,0:nl(2)+1,0:nl(3)+1), intent(in) :: U1, V1, W1
  Real(sp), Dimension(0:nl(1)+1,0:nl(2)+1,0:nl(3)+1), intent(in) :: rhox, rhoy, rhoz
  Real(sp), Dimension(0:nl(1)+1,0:nl(2)+1,0:nl(3)+1), intent(in) :: p
  Real(sp), Dimension(0:nl(1)+1,0:nl(2)+1,0:nl(3)+1), intent(out) :: U2, V2, W2
  Integer :: i, j, k

  do k=1,nl(3)
    do j=1,nl(2)
      do i=1,nl(1)
        U2(i,j,k) = U1(i,j,k) - rhox(i,j,k) * (p(i+1,  j,  k) - p(i,j,k)) * dt1 / dl(1)
        V2(i,j,k) = V1(i,j,k) - rhoy(i,j,k) * (p(  i,j+1,  k) - p(i,j,k)) * dt1 / dl(2)
        W2(i,j,k) = W1(i,j,k) - rhoz(i,j,k) * (p(  i,  j,k+1) - p(i,j,k)) * dt1 / dl(3)
      enddo
    enddo
  enddo

End Subroutine DuDtPartialP

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
Subroutine Divergence(up, vp, wp, nl, dl, div)
  Use ModGlobal, only : sp
  Implicit None
  Integer :: nl(3)
  Real(sp) :: dl(3)
  real(sp), intent(in) , dimension(0:nl(1)+1,0:nl(2)+1,0:nl(3)+1) :: up, vp, wp
  real(sp), intent(out) , dimension(1:nl(1),1:nl(2),1:nl(3)) :: div
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
! Subroutine CSF_VOSET(dudt, dvdt, dwdt, Phi)
  ! Implicit None
  ! Integer :: nl(3)
  ! Real(sp) :: dl(3)
  ! real(sp), intent(in), dimension(0:,0:,0:) :: Phi
  ! real(sp), intent(inout), dimension(:,:,:) :: dudt, dvdt, dwdt
  ! real(sp), dimension(0:nl(1)+1,0:nl(2)+1,0:nl(3)+1) :: Phi_DS
  ! real(sp), dimension(0:nl(1)+1,0:nl(2)+1,0:nl(3)+1) :: kappa_v
  ! real(sp), dimension(0:nl(1)+1,0:nl(2)+1,0:nl(3)+1) :: kappa_sx, kappa_sy, kappa_sz
  ! Integer :: i, j, k
! 
    !  (1) Scale volume function
  ! Call Smooth_Phi(Phi, Phi_DS, nl(1), nl(2), nl(3))
  ! Call Phi_bc%SetBCS(Phi_DS)
! 
    !  (2) Convert the VOF function to level set
  ! Call Get_LS_3D_init(Phi_DS, Ls, dl(1), dl(2), dl(3), nl(1), nl(2), nl(3))
! 
  ! Call ls_bc%SetBCS(LS)
    !  (3) calculate curvature at grid center
  ! Kappa_v = 0.0_sp
  ! Do k=1,nl(3)
    ! Do j=1,nl(2)
      ! Do i=1,nl(1)
        ! If (Phi_DS(i,j,k) < 1.0_sp .or. Phi_DS(i,j,k) > 0.0_sp) Then
          ! Kappa_v(i,j,k) = &
              ! (    ( Ls(i-1,  j,  k) - 2.d0 * Ls(i,j,k) + Ls(i+1,  j,  k) ) / dl(1) / dl(1) &
              ! &  + ( Ls(  i,j-1,  k) - 2.d0 * Ls(i,j,k) + Ls(  i,j+1,  k) ) / dl(2) / dl(2) &
              ! &  + ( Ls(  i,  j,k-1) - 2.d0 * Ls(i,j,k) + Ls(  i,  j,k+1) ) / dl(3) / dl(3) )
        ! Else
          ! Kappa_v(i,j,k) = 0.0_sp
        ! End If
      ! EndDo
    ! EndDo
  ! EndDo
  ! Call Phi_bc%SetBCS(kappa_v)
! 
    !  (4) Interpolate the curvature to grid center
  ! Do k=0,nl(3)
    ! Do j=0,nl(2)
      ! Do i=0,nl(1)
        ! Kappa_sx(i,j,k) = 0.5_sp * ( Kappa_v(i,j,k) + Kappa_v(i+1,  j,  k) )
        ! Kappa_sy(i,j,k) = 0.5_sp * ( Kappa_v(i,j,k) + Kappa_v(  i,j+1,  k) )
        ! Kappa_sz(i,j,k) = 0.5_sp * ( Kappa_v(i,j,k) + Kappa_v(  i,  j,k+1) )
      ! EndDo
    ! EndDo
  ! EndDo
! 
    !  (5) Calculte density-scaled balanced surface tension
  ! Do k=1,nl(3)
    ! Do j=1,nl(2)
      ! Do i=1,nl(1)
        ! dudt(i,j,k) = dudt(i,j,k) - sigma * Kappa_sx(i,j,k) * ( Phi_DS(i+1,  j,  k) - Phi_DS(i,j,k) ) * rhox(i,j,k) / dl(1)
        ! dvdt(i,j,k) = dvdt(i,j,k) - sigma * Kappa_sy(i,j,k) * ( Phi_DS(  i,j+1,  k) - Phi_DS(i,j,k) ) * rhoy(i,j,k) / dl(2)
        ! dwdt(i,j,k) = dwdt(i,j,k) - sigma * Kappa_sz(i,j,k) * ( Phi_DS(  i,  j,k+1) - Phi_DS(i,j,k) ) * rhoz(i,j,k) / dl(3)
      ! EndDo
    ! EndDo
  ! EndDo
! 
! End Subroutine CSF_VOSET

! Subroutine CSF_NULL(dudt, dvdt, dwdt, Phi)
!   Implicit None
!   real(sp), intent(in), dimension(0:,0:,0:) :: Phi
!   real(sp), intent(inout), dimension(:,:,:) :: dudt, dvdt, dwdt
! End Subroutine CSF_NULL

!==========================
!  (4-1) Adjust dt
!==========================
Subroutine Adjustdt(U, V, W, nl, dl, dt, dt0, dt_min, cfl)
  use mpi
  Use ModGlobal, only : sp
  Implicit None
  Integer :: nl(3)
  Real(sp) :: dl(3)
  Real(sp) :: dt, dt0, dt_min
  Real(sp), Intent(In), Dimension(0:nl(1)+1,0:nl(2)+1,0:nl(3)+1) :: U, V, W
  Real(sp) :: uvw, uvw2, dlmin
  Real(sp), Intent(In) :: cfl
  Integer :: i, j, k
  Integer :: ierr
  uvw = 0.0_sp
  Do k = 1, nl(3)
    Do j = 1, nl(2)
      Do i = 1, nl(1)
        uvw = max(uvw,abs(u(i,j,k))+abs(v(i,j,k))+abs(w(i,j,k)))
      End Do
    End Do
  End Do

  dlmin = min(min(dl(1),dl(2)),dl(3))

  Call MPI_Reduce(uvw, uvw2, 1, MPI_REAL_SP, MPI_MAX, MPI_COMM_WORLD, 0, ierr)

  dt = cfl * sqrt(3.0_sp) * dlmin / uvw

  dt = min(dt0, dt)
  dt = max(dt, dt_min)

  Call MPI_Barrier(MPI_COMM_WORLD, ierr)

End Subroutine Adjustdt

  !==========================
  !  (4-2) Jacobi iteration
  !==========================
  ! Subroutine Jacobi(P)
  !   Use ModGlobal, only : dl, nl, dt, dt0
  !   Implicit None
  !   Real(sp), Intent(InOut) :: P(0:,0:,0:)
  !   Real(sp) :: res(0:nl(1)+1,0:nl(2)+1,0:nl(3)+1)
  !   Real(sp) :: res_max, res_max2
  !   Real(sp) :: rho_xx_L, rho_xx_R
  !   Real(sp) :: rho_yy_U, rho_yy_D
  !   Real(sp) :: rho_zz_F, rho_zz_B
  !   Real(sp) :: aa
  !   Integer :: i, j, k, ll

  !   Do ll = 1, 10000
  !     Do k = 1, nl(3)
  !       Do j = 1, nl(2)
  !         Do i = 1, nl(1)

  !           rho_xx_R =  Rhox(  i,j,k) / dl(1) / dl(1)
  !           rho_xx_L =  Rhox(i-1,j,k) / dl(1) / dl(1)
  !           rho_yy_U =  Rhoy(i,  j,k) / dl(2) / dl(2)
  !           rho_yy_D =  Rhoy(i,j-1,k) / dl(2) / dl(2)
  !           rho_zz_B =  Rhoz(i,j,  k) / dl(3) / dl(3)
  !           rho_zz_F =  Rhoz(i,j,k-1) / dl(3) / dl(3)

  !           Aa = - ( rho_xx_R + rho_xx_L   &
  !           &     + rho_yy_U + rho_yy_D   &
  !           &     + rho_zz_B + rho_zz_F ) 

  !           Res(i,j,k) = ( DIV(i,j,k) / dt &
  !            &  - rho_xx_L * P(i-1,j,k) - rho_xx_R * P(i+1,j,k)    &
  !            &  - rho_yy_D * P(i,j-1,k) - rho_yy_U * P(i,j+1,k)    &
  !            &  - rho_zz_F * P(i,j,k-1) - rho_zz_B * P(i,j,k+1) )  &
  !            &  / Aa - P(i,j,k)
  !         End Do
  !       End Do
  !     End Do

  !     res_max = 0.0_sp
  !     Do k = 1, nl(3)
  !       Do j = 1, nl(2)
  !         Do i = 1, nl(1)
  !           if (abs(res(i,j,k)) > res_max ) then
  !             res_max = res(i,j,k)
  !           end if
  !             P(i,j,k) = P(i,j,k) + res(i,j,k)
  !         End Do
  !       End Do
  !     End Do
  !     Call P_bc%SetBCS(P)
  !     Call MPI_AllReduce(res_max, res_max2, 1, MPI_REAL_SP, MPI_MAX, MPI_COMM_WORLD, ierr)
  !     if (res_max2 < iter_tolerance) exit
  !   End Do

  !   n_iter = ll

  !   if (myid .eq.0 ) print *, myid,ll, res_max2

! Subroutine AdvDiffU2(u, v, w, dudt)
!     Implicit None
!     Real(sp), Dimension(0:,0:,0:), Intent(in)  :: u, v, w
!     Real(sp), Dimension(:,:,:), Intent(out) :: dudt
!     integer :: i,j,k
!     Real(sp),Dimension(0:nl(1)+1,0:nl(2)+1,0:nl(3)+1) :: uu, uv, uw, dudx, dudy, dudz

!     ! Calculate the advection-diffusion spatial term
!     Do k=0,nl(3)
!       Do j=0,nl(2)
!         Do i=1,nl(1) + 1
!           dudx(i,j,k) =   mu(i,j,k) * ( u(i,  j,  k) - u(i-1,j,k) ) / dl(1)
!           dudy(i,j,k) = muxz(i,j,k) * ( u(i,j+1,  k) - u(  i,j,k) ) / dl(2)
!           dudz(i,j,k) = muxy(i,j,k) * ( u(i,  j,k+1) - u(  i,j,k) ) / dl(3)
!         Enddo
!       Enddo
!     Enddo

!     Do k=1,nl(3)
!       Do j=1,nl(2)
!         Do i=1,nl(1) 
!           uu(i,j,k) = u(i,j,k)
!           uv(i,j,k) = 0.25_sp * ( v(i,j,k) + v(i,j-1,k) + v(i+1,j,k) + v(i+1,j-1,k) )
!           uw(i,j,k) = 0.25_sp * ( w(i,j,k) + w(i,j,k-1) + w(i+1,j,k) + w(i+1,j,k-1))
!         Enddo
!       Enddo
!     Enddo

!     Do k=1,nl(3)
!       Do j=1,nl(2)
!         Do i=1,nl(1)
!           dudt(i,j,k) = &
!               (0.5_sp * uu(i,j,k) * (u(i-1,j,k)-u(i+1,j,k)) + (dudx(i+1,j,k)-dudx(i,  j,  k)) * rhox(i,j,k) ) / dl(1) + &
!               (0.5_sp * uv(i,j,k) * (u(i,j-1,k)-u(i,j+1,k)) + (dudy(  i,j,k)-dudy(i,j-1,  k)) * rhox(i,j,k) ) / dl(2) + &
!               (0.5_sp * uw(i,j,k) * (u(i,j,k-1)-u(i,j,k+1)) + (dudz(  i,j,k)-dudz(i,  j,k-1)) * rhox(i,j,k) ) / dl(3)
!           ! dudt(i,j,k) = &
!           !     (-uu(i+1,j,k) + uu(i,  j,  k)) / dl(1) + &
!           !     (-uv(  i,j,k) + uv(i,j-1,  k)) / dl(2) + &
!           !     (-uw(  i,j,k) + uw(i,  j,k-1)) / dl(3)
!         Enddo
!       Enddo
!     Enddo

!   End Subroutine AdvDiffU2

! Subroutine AdvDiffV2(u, v, w, dvdt)
!     Implicit None
!     Real(sp), Dimension(0:,0:,0:), Intent(in)  :: u, v, w
!     Real(sp), Dimension(:,:,:), Intent(out) :: dvdt
!     integer :: i,j,k
!     Real(sp),Dimension(0:nl(1)+1,0:nl(2)+1,0:nl(3)+1) :: uv, vv, vw, dvdx, dvdy, dvdz

!     ! Calculate the advection-diffusion spatial term
!     Do k=0,nl(3)
!       Do j=1,nl(2) + 1
!         Do i=0,nl(1)
!           dvdx(i,j,k) = muyz(i,j,k) * ( v(i+1,j,  k) - v(i,  j,k) ) / dl(1)
!           dvdy(i,j,k) =   mu(i,j,k) * ( v(  i,j,  k) - v(i,j-1,k) ) / dl(2)
!           dvdz(i,j,k) = muxy(i,j,k) * ( v(  i,j,k+1) - v(i,  j,k) ) / dl(3)
!         Enddo
!       Enddo
!     Enddo

!     Do k=1,nl(3)
!       Do j=1,nl(2)
!         Do i=1,nl(1) 
!           uv(i,j,k) = 0.25_sp * ( u(i,j,k) + u(i,j+1,k) + u(i-1,j,k) + u(i-1,j+1,k) )
!           vv(i,j,k) = v(i,j,k)
!           vw(i,j,k) = 0.25_sp * ( w(i,j,k) + w(i,j,k-1) + w(i,j+1,k) + w(i,j+1,k-1) )
!         Enddo
!       Enddo
!     Enddo

!     Do k=1,nl(3)
!       Do j=1,nl(2)
!         Do i=1,nl(1)
!           dvdt(i,j,k) = &
!               (0.5_sp * uv(i,j,k) * (v(i-1,j,k)-v(i+1,j,k)) + (dvdx(i,  j,k)-dvdx(i-1,j,  k)) * rhoy(i,j,k) ) / dl(1) + &
!               (0.5_sp * vv(i,j,k) * (v(i,j-1,k)-v(i,j+1,k)) + (dvdy(i,j+1,k)-dvdy(  i,j,  k)) * rhoy(i,j,k) ) / dl(2) + &
!               (0.5_sp * vw(i,j,k) * (v(i,j,k-1)-v(i,j,k+1)) + (dvdz(i,  j,k)-dvdz(  i,j,k-1)) * rhoy(i,j,k) ) / dl(3)
!           ! dudt(i,j,k) = &
!           !     (-uu(i+1,j,k) + uu(i,  j,  k)) / dl(1) + &
!           !     (-uv(  i,j,k) + uv(i,j-1,  k)) / dl(2) + &
!           !     (-uw(  i,j,k) + uw(i,  j,k-1)) / dl(3)
!         Enddo
!       Enddo
!     Enddo

!   End Subroutine AdvDiffV2

!   Subroutine AdvDiffW2(u, v, w, dwdt)
!     Implicit None
!     Real(sp), Dimension(0:,0:,0:), Intent(in)  :: u, v, w
!     Real(sp), Dimension(:,:,:), Intent(out) :: dwdt
!     integer :: i,j,k
!     Real(sp),Dimension(0:nl(1)+1,0:nl(2)+1,0:nl(3)+1) :: uw, vw, ww, dwdx, dwdy, dwdz

!     ! Calculate the advection-diffusion spatial term
!     Do k=1,nl(3) + 1
!       Do j=0,nl(2)
!         Do i=0,nl(1)
!           uw(i,j,k) = 0.25_sp * ( u(i,j,k) + u(i+1,j,k) + u(i,j,k-1) + u(i+1,j,k-1) )
!           vw(i,j,k) = 0.25_sp * ( v(i,j,k) + v(i,j,k-1) + v(i,j+1,k) + v(i,j+1,k-1))
!           ww(i,j,k) = w(i,j,k) 
!           uw(i,j,k) = 0.25_sp * (u(i,j,k-1) + u(i,j,k)) * (w(i+1,  j,  k)+ w(i,j,k))
!           vw(i,j,k) = 0.25_sp * (v(i,j,k-1) + v(i,j,k)) * (w(  i,j+1,  k)+ w(i,j,k))
!           ww(i,j,k) = 0.25_sp * (w(i,j,k-1) + w(i,j,k)) * (w(  i,  j,k-1)+ w(i,j,k))
!           dwdx(i,j,k) = muyz(i,j,k) * ( w(i+1,  j,k) - w(i,j,  k) ) / dl(1)
!           dwdy(i,j,k) = muxz(i,j,k) * ( w(  i,j+1,k) - w(i,j,  k) ) / dl(2)
!           dwdz(i,j,k) =   mu(i,j,k) * ( w(  i,  j,k) - w(i,j,k-1) ) / dl(3)
!         Enddo
!       Enddo
!     Enddo

!     Do k=1,nl(3)
!       Do j=1,nl(2)
!         Do i=1,nl(1) 
!           uw(i,j,k) = 0.25_sp * ( u(i,j,k) + u(i,j,k+1) + u(i-1,j,k) + u(i-1,j,k+1) )
!           vw(i,j,k) = 0.25_sp * ( v(i,j,k) + v(i,j,k+1) + v(i,j-1,k) + w(i,j-1,k+1) )
!           ww(i,j,k) = w(i,j,k)
!         Enddo
!       Enddo
!     Enddo

!     Do k=1,nl(3)
!       Do j=1,nl(2)
!         Do i=1,nl(1)
!           dwdt(i,j,k) = &
!               (0.5_sp * uw(i,j,k) * (w(i-1,j,k)-w(i+1,j,k)) + (dwdx(i,j,  k)-dwdx(i-1,  j,k)) * rhoz(i,j,k) ) / dl(1) + &
!               (0.5_sp * vw(i,j,k) * (w(i,j-1,k)-w(i,j+1,k)) + (dwdy(i,j,  k)-dwdy(  i,j-1,k)) * rhoz(i,j,k) ) / dl(2) + &
!               (0.5_sp * ww(i,j,k) * (w(i,j,k-1)-w(i,j,k+1)) + (dwdz(i,j,k+1)-dwdz(  i,  j,k)) * rhoz(i,j,k) ) / dl(3)
!               ! (0.5_sp * uw(i,j,k) * (w(i-1,j,k)-w(i+1,j,k)) + (dvdx(i,  j,k)-dvdx(i-1,j,  k)) * rhoz(i,j,k) ) / dl(1) + &
!               ! (0.5_sp * vw(i,j,k) * (w(i,j-1,k)-w(i,j+1,k)) + (dvdy(i,j+1,k)-dvdy(  i,j,  k)) * rhoz(i,j,k) ) / dl(2) + &
!               ! (0.5_sp * ww(i,j,k) * (w(i,j,k-1)-w(i,j,k+1)) + (dvdz(i,  j,k)-dvdz(  i,j,k-1)) * rhoz(i,j,k) ) / dl(3)
!               ! (-uw(i,j,  k) + uw(i-1,  j,k)) / dl(1) + &
!               ! (-vw(i,j,  k) + vw(  i,j-1,k)) / dl(2) + &
!               ! (-ww(i,j,k+1) + ww(  i,  j,k)) / dl(3)
!         Enddo
!       Enddo
!     Enddo

!   End Subroutine AdvDiffW2

Subroutine AdvDiffU_UPCDS(u, v, w, dl, nl, mu, muxy, muxz, rhox, dudt)
  Use ModGlobal, only : sp
    Implicit None
  Real(sp), Intent(In) :: dl(3)
  Integer, Intent(In) :: nl(3)
  Real(sp), Dimension(0:nl(1)+1,0:nl(2)+1,0:nl(3)+1), Intent(in)  :: u, v, w
  Real(sp), Dimension(0:nl(1)+1,0:nl(2)+1,0:nl(3)+1), Intent(in)  :: mu, muxy, muxz, rhox
  Real(sp), Dimension(nl(1),nl(2),nl(3)), Intent(out) :: dudt
  integer :: i,j,k
  Real(sp),Dimension(0:nl(1)+1,0:nl(2)+1,0:nl(3)+1) :: uu, uv, uw, dudx, dudy, dudz
  Real(sp),Dimension(0:nl(1)+1,0:nl(2)+1,0:nl(3)+1) :: dudxup, dudyup, dudzup

  ! Calculate the advection-diffusion spatial term
  Do k=0,nl(3)
    Do j=0,nl(2)
      Do i=1,nl(1) + 1
        dudx(i,j,k) =   mu(i,j,k) * ( u(i,  j,  k) - u(i-1,j,k) ) / dl(1)
        dudy(i,j,k) = muxz(i,j,k) * ( u(i,j+1,  k) - u(  i,j,k) ) / dl(2)
        dudz(i,j,k) = muxy(i,j,k) * ( u(i,  j,k+1) - u(  i,j,k) ) / dl(3)
      Enddo
    Enddo
  Enddo

  Do k=1,nl(3)
    Do j=1,nl(2)
      Do i=1,nl(1) 
        uu(i,j,k) = u(i,j,k)
        uv(i,j,k) = 0.25_sp * ( v(i,j,k) + v(i,j-1,k) + v(i+1,j,k) + v(i+1,j-1,k) )
        uw(i,j,k) = 0.25_sp * ( w(i,j,k) + w(i,j,k-1) + w(i+1,j,k) + w(i+1,j,k-1))
        if( uu(i,j,k) > 0) then
          dudxup(i,j,k) = u(i-1,j,k) - u(i,j,k)
        else
          dudxup(i,j,k) = u(i,j,k) - u(i+1,j,k)
        endif
        if( uv(i,j,k) > 0) then
          dudyup(i,j,k) = u(i,j-1,k) - u(i,j,k)
        else
          dudyup(i,j,k) = u(i,j,k) - u(i,j+1,k)
        endif
        if( uw(i,j,k) > 0) then
          dudzup(i,j,k) = u(i,j,k-1) - u(i,j,k)
        else
          dudzup(i,j,k) = u(i,j,k) - u(i,j,k+1)
        endif
      Enddo
    Enddo
  Enddo

  Do k=1,nl(3)
    Do j=1,nl(2)
      Do i=1,nl(1)
        dudt(i,j,k) = &
            ( uu(i,j,k) * dudxup(i,j,k) + (dudx(i+1,j,k)-dudx(i,  j,  k)) * rhox(i,j,k) ) / dl(1) + &
            ( uv(i,j,k) * dudyup(i,j,k) + (dudy(  i,j,k)-dudy(i,j-1,  k)) * rhox(i,j,k) ) / dl(2) + &
            ( uw(i,j,k) * dudzup(i,j,k) + (dudz(  i,j,k)-dudz(i,  j,k-1)) * rhox(i,j,k) ) / dl(3)
      Enddo
    Enddo
  Enddo

End Subroutine AdvDiffU_UPCDS

Subroutine AdvDiffV_UPCDS(u, v, w, dl, nl, mu, muxy, muyz, rhoy, dvdt)
  Use ModGlobal, only : sp
  Implicit None
  Real(sp), Intent(In) :: dl(3)
  Integer, Intent(In) :: nl(3)
  Real(sp), Dimension(0:nl(1)+1,0:nl(2)+1,0:nl(3)+1), Intent(in)  :: u, v, w
  Real(sp), Dimension(0:nl(1)+1,0:nl(2)+1,0:nl(3)+1), Intent(in)  :: mu, muyz, muxy, rhoy
  Real(sp), Dimension(nl(1),nl(2),nl(3)), Intent(out) :: dvdt
  integer :: i,j,k
  Real(sp),Dimension(0:nl(1)+1,0:nl(2)+1,0:nl(3)+1) :: uv, vv, vw, dvdx, dvdy, dvdz
  Real(sp),Dimension(0:nl(1)+1,0:nl(2)+1,0:nl(3)+1) :: dvdxup, dvdyup, dvdzup

  ! Calculate the advection-diffusion spatial term
  Do k=0,nl(3)
    Do j=1,nl(2) + 1
      Do i=0,nl(1)
        dvdx(i,j,k) = muyz(i,j,k) * ( v(i+1,j,  k) - v(i,  j,k) ) / dl(1)
        dvdy(i,j,k) =   mu(i,j,k) * ( v(  i,j,  k) - v(i,j-1,k) ) / dl(2)
        dvdz(i,j,k) = muxy(i,j,k) * ( v(  i,j,k+1) - v(i,  j,k) ) / dl(3)
      Enddo
    Enddo
  Enddo

  Do k=1,nl(3)
    Do j=1,nl(2)
      Do i=1,nl(1) 
        uv(i,j,k) = 0.25_sp * ( u(i,j,k) + u(i,j+1,k) + u(i-1,j,k) + u(i-1,j+1,k) )
        vv(i,j,k) = v(i,j,k)
        vw(i,j,k) = 0.25_sp * ( w(i,j,k) + w(i,j,k-1) + w(i,j+1,k) + w(i,j+1,k-1) )
        if( uv(i,j,k) > 0) then
          dvdxup(i,j,k) = v(i-1,j,k) - v(i,j,k)
        else
          dvdxup(i,j,k) = v(i,j,k) - v(i+1,j,k)
        endif
        if( vv(i,j,k) > 0) then
          dvdyup(i,j,k) = v(i,j-1,k) - v(i,j,k)
        else
          dvdyup(i,j,k) = v(i,j,k) - v(i,j+1,k)
        endif
        if( vw(i,j,k) > 0) then
          dvdzup(i,j,k) = v(i,j,k-1) - v(i,j,k)
        else
          dvdzup(i,j,k) = v(i,j,k) - v(i,j,k+1)
        endif
      Enddo
    Enddo
  Enddo

  Do k=1,nl(3)
    Do j=1,nl(2)
      Do i=1,nl(1)
        dvdt(i,j,k) = &
            (uv(i,j,k) * dvdxup(i,j,k) + (dvdx(i,  j,k)-dvdx(i-1,j,  k)) * rhoy(i,j,k) ) / dl(1) + &
            (vv(i,j,k) * dvdyup(i,j,k) + (dvdy(i,j+1,k)-dvdy(  i,j,  k)) * rhoy(i,j,k) ) / dl(2) + &
            (vw(i,j,k) * dvdzup(i,j,k) + (dvdz(i,  j,k)-dvdz(  i,j,k-1)) * rhoy(i,j,k) ) / dl(3)
      Enddo
    Enddo
  Enddo

End Subroutine AdvDiffV_UPCDS

Subroutine AdvDiffW_UPCDS(u, v, w, dl, nl, mu, muxz, muyz, rhoz, dwdt)
  Use ModGlobal, only : sp
  Implicit None
  Real(sp), Intent(In) :: dl(3)
  Integer, Intent(In) :: nl(3)
  Real(sp), Dimension(0:nl(1)+1,0:nl(2)+1,0:nl(3)+1), Intent(in)  :: u, v, w
  Real(sp), Dimension(0:nl(1)+1,0:nl(2)+1,0:nl(3)+1), Intent(in)  :: mu, muxz, muyz, rhoz
  Real(sp), Dimension(nl(1),nl(2),nl(3)), Intent(out) :: dwdt
  integer :: i,j,k
  Real(sp),Dimension(0:nl(1)+1,0:nl(2)+1,0:nl(3)+1) :: uw, vw, ww, dwdx, dwdy, dwdz
  Real(sp),Dimension(0:nl(1)+1,0:nl(2)+1,0:nl(3)+1) :: dwdxup, dwdyup, dwdzup

  ! Calculate the advection-diffusion spatial term
  Do k=1,nl(3) + 1
    Do j=0,nl(2)
      Do i=0,nl(1)
        dwdx(i,j,k) = muyz(i,j,k) * ( w(i+1,  j,k) - w(i,j,  k) ) / dl(1)
        dwdy(i,j,k) = muxz(i,j,k) * ( w(  i,j+1,k) - w(i,j,  k) ) / dl(2)
        dwdz(i,j,k) =   mu(i,j,k) * ( w(  i,  j,k) - w(i,j,k-1) ) / dl(3)
      Enddo
    Enddo
  Enddo

  Do k=1,nl(3)
    Do j=1,nl(2)
      Do i=1,nl(1) 
        uw(i,j,k) = 0.25_sp * ( u(i,j,k) + u(i,j,k+1) + u(i-1,j,k) + u(i-1,j,k+1) )
        vw(i,j,k) = 0.25_sp * ( v(i,j,k) + v(i,j,k+1) + v(i,j-1,k) + w(i,j-1,k+1) )
        ww(i,j,k) = w(i,j,k)
        if( uw(i,j,k) > 0) then
          dwdxup(i,j,k) = w(i-1,j,k) - w(i,j,k)
        else
          dwdxup(i,j,k) = w(i,j,k) - w(i+1,j,k)
        endif
        if( vw(i,j,k) > 0) then
          dwdyup(i,j,k) = w(i,j-1,k) - w(i,j,k)
        else
          dwdyup(i,j,k) = w(i,j,k) - w(i,j+1,k)
        endif
        if( ww(i,j,k) > 0) then
          dwdzup(i,j,k) = w(i,j,k-1) - w(i,j,k)
        else
          dwdzup(i,j,k) = w(i,j,k) - w(i,j,k+1)
        endif
      Enddo
    Enddo
  Enddo

  Do k=1,nl(3)
    Do j=1,nl(2)
      Do i=1,nl(1)
        dwdt(i,j,k) = &
            (uw(i,j,k) * dwdxup(i,j,k) + (dwdx(i,j,  k)-dwdx(i-1,  j,k)) * rhoz(i,j,k) ) / dl(1) + &
            (vw(i,j,k) * dwdyup(i,j,k) + (dwdy(i,j,  k)-dwdy(  i,j-1,k)) * rhoz(i,j,k) ) / dl(2) + &
            (ww(i,j,k) * dwdzup(i,j,k) + (dwdz(i,j,k+1)-dwdz(  i,  j,k)) * rhoz(i,j,k) ) / dl(3)
      Enddo
    Enddo
  Enddo

End Subroutine AdvDiffW_UPCDS