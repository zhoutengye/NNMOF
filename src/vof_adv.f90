!=================
! Includes:
!     (1) VOF/MOF Advection
!       (1.1) VOF-CIAM (Lagrangian)
!       (1.2) VOF-WY (Eulerian)
!       (1.3) MOF-CIAM (Lanrangian)
!       (1.4) MOF-WY (Eulerian) (!!! NOT WORKING)
!     (2) Directional split VOF/MOF advection
!       (1.1) VOF-CIAM (Lagrangian)
!       (1.2) VOF-WY (Eulerian)
!       (1.3) MOF-CIAM (Lanrangian)
!       (1.4) MOF-WY (Eulerian) (!!! NOT WORKING)
!     (3) Advection of centroid
!       (3.1) Langrangian advection
!       (3.2) Eulerian advection
!-----------------
! Author: Zhouteng Ye (yzt9zju@gmail.com)
!=================
Module ModVOF
  Use ModGlobal, only : sp
  Use ModGlobal, only : updthalo
  Use ModGlobal, only : setBCS, phi_bc
  Use ModVOFFunc
#if defined(VISUALIZE)
  Use ModTools
#endif

Contains

  Subroutine VOFCIAM(Phi, u, v, w, nl, dl, dt)
    Implicit None
    Real(sp) :: dt
    Integer,Intent(In)     :: nl(3)
    Real(sp),Intent(In)    :: dl(3)
    Real(sp),Intent(InOut)    :: Phi(0:nl(1)+1,0:nl(2)+1,0:nl(3)+1)
    Real(sp),Intent(In)       :: u(0:nl(1)+1,0:nl(2)+1,0:nl(3)+1)
    Real(sp),Intent(In)       :: v(0:nl(1)+1,0:nl(2)+1,0:nl(3)+1)
    Real(sp),Intent(In)       :: w(0:nl(1)+1,0:nl(2)+1,0:nl(3)+1)
    call AdvCIAM(u, phi, nl, dl, dt, 1)
    call AdvCIAM(v, phi, nl, dl, dt, 2)
    call AdvCIAM(w, phi, nl, dl, dt, 3)
  End Subroutine VOFCIAM

  Subroutine VOFWY(Phi, u, v, w, nl, dl, dt)
    Implicit None
    Real(sp) :: dt
    Integer,Intent(In)     :: nl(3)
    Real(sp),Intent(In)    :: dl(3)
    Real(sp),Intent(InOut)    :: Phi(0:nl(1)+1,0:nl(2)+1,0:nl(3)+1)
    Real(sp),Intent(In)       :: u(0:nl(1)+1,0:nl(2)+1,0:nl(3)+1)
    Real(sp),Intent(In)       :: v(0:nl(1)+1,0:nl(2)+1,0:nl(3)+1)
    Real(sp),Intent(In)       :: w(0:nl(1)+1,0:nl(2)+1,0:nl(3)+1)
    call AdvWY(u, phi, nl, dl, dt, 1)
    call AdvWY(v, phi, nl, dl, dt, 2)
    call AdvWY(w, phi, nl, dl, dt, 3)
  End Subroutine VOFWY

  Subroutine MOFCIAM(Phi, cx, cy, cz, u, v, w, nl, dl, dt)
    Implicit None
    Real(sp) :: dt
    Integer,Intent(In)     :: nl(3)
    Real(sp),Intent(In)    :: dl(3)
    Real(sp),Intent(InOut)    :: Phi(0:nl(1)+1,0:nl(2)+1,0:nl(3)+1)
    Real(sp),Intent(InOut)    :: cx(0:nl(1)+1,0:nl(2)+1,0:nl(3)+1)
    Real(sp),Intent(InOut)    :: cy(0:nl(1)+1,0:nl(2)+1,0:nl(3)+1)
    Real(sp),Intent(InOut)    :: cz(0:nl(1)+1,0:nl(2)+1,0:nl(3)+1)
    Real(sp),Intent(In)       :: u(0:nl(1)+1,0:nl(2)+1,0:nl(3)+1)
    Real(sp),Intent(In)       :: v(0:nl(1)+1,0:nl(2)+1,0:nl(3)+1)
    Real(sp),Intent(In)       :: w(0:nl(1)+1,0:nl(2)+1,0:nl(3)+1)
    call AdvCIAM_MOF(u, cx, cy, cz, phi, nl, dl, dt, 1)
    call AdvCIAM_MOF(v, cx, cy, cz, phi, nl, dl, dt, 2)
    call AdvCIAM_MOF(w, cx, cy, cz, phi, nl, dl, dt, 3)
  End Subroutine MOFCIAM

  Subroutine MOFCIAM2(Phi, cx, cy, cz, u, v, w, nl, dl, dt,rank)
    Implicit None
    Real(sp) :: dt
    Integer,Intent(In)     :: nl(3)
    Real(sp),Intent(In)    :: dl(3)
    Real(sp),Intent(InOut)    :: Phi(0:nl(1)+1,0:nl(2)+1,0:nl(3)+1)
    Real(sp),Intent(InOut)    :: cx(0:nl(1)+1,0:nl(2)+1,0:nl(3)+1)
    Real(sp),Intent(InOut)    :: cy(0:nl(1)+1,0:nl(2)+1,0:nl(3)+1)
    Real(sp),Intent(InOut)    :: cz(0:nl(1)+1,0:nl(2)+1,0:nl(3)+1)
    Real(sp),Intent(In)       :: u(0:nl(1)+1,0:nl(2)+1,0:nl(3)+1)
    Real(sp),Intent(In)       :: v(0:nl(1)+1,0:nl(2)+1,0:nl(3)+1)
    Real(sp),Intent(In)       :: w(0:nl(1)+1,0:nl(2)+1,0:nl(3)+1)
    Integer :: rank
    if (rank == 1) Then
      call AdvCIAM_MOF(u, cx, cy, cz, phi, nl, dl, dt, 1)
      call AdvCIAM_MOF(v, cx, cy, cz, phi, nl, dl, dt, 2)
      call AdvCIAM_MOF(w, cx, cy, cz, phi, nl, dl, dt, 3)
    else  if (rank == 2) Then
      call AdvCIAM_MOF(w, cx, cy, cz, phi, nl, dl, dt, 3)
      call AdvCIAM_MOF(u, cx, cy, cz, phi, nl, dl, dt, 1)
      call AdvCIAM_MOF(v, cx, cy, cz, phi, nl, dl, dt, 2)
    else  if (rank == 0) Then
      call AdvCIAM_MOF(v, cx, cy, cz, phi, nl, dl, dt, 2)
      call AdvCIAM_MOF(w, cx, cy, cz, phi, nl, dl, dt, 3)
      call AdvCIAM_MOF(u, cx, cy, cz, phi, nl, dl, dt, 1)
    end if
  End Subroutine MOFCIAM2

  Subroutine MOFWY(Phi, cx, cy, cz, u, v, w, nl, dl, dt)
    Implicit None
    Real(sp) :: dt
    Integer,Intent(In)     :: nl(3)
    Real(sp),Intent(In)    :: dl(3)
    Real(sp),Intent(InOut)    :: Phi(0:nl(1)+1,0:nl(2)+1,0:nl(3)+1)
    Real(sp),Intent(InOut)    :: cx(0:nl(1)+1,0:nl(2)+1,0:nl(3)+1)
    Real(sp),Intent(InOut)    :: cy(0:nl(1)+1,0:nl(2)+1,0:nl(3)+1)
    Real(sp),Intent(InOut)    :: cz(0:nl(1)+1,0:nl(2)+1,0:nl(3)+1)
    Real(sp),Intent(In)       :: u(0:nl(1)+1,0:nl(2)+1,0:nl(3)+1)
    Real(sp),Intent(In)       :: v(0:nl(1)+1,0:nl(2)+1,0:nl(3)+1)
    Real(sp),Intent(In)       :: w(0:nl(1)+1,0:nl(2)+1,0:nl(3)+1)
    call AdvWY_MOF(u, cx, cy, cz, phi, nl, dl, dt, 1)
    call AdvWY_MOF(v, cx, cy, cz, phi, nl, dl, dt, 2)
    call AdvWY_MOF(w, cx, cy, cz, phi, nl, dl, dt, 3)
  End Subroutine MOFWY



!=======================================================
! Advection algorithm of CIAM PLIC
!=======================================================
Subroutine AdvCIAM(us, f, nl, dl, dt, dir)
  ! Use ModVOFFunc
  Implicit None
  Integer :: dir
  Real(sp) :: dt
  Integer,Intent(In)     :: nl(3)
  Real(sp),Intent(In)    :: dl(3)
  Real(sp),Intent(In)    :: us(0:nl(1)+1,0:nl(2)+1,0:nl(3)+1)
  Real(sp),Intent(InOut) :: f(0:nl(1)+1,0:nl(2)+1,0:nl(3)+1)
  Integer  :: i, j, k, ii, jj, kk, ic1, ic2, ic3
  Real(sp) :: a1,a2,alpha
  Real(sp) :: x0(3), deltax(3)
  Real(sp),dimension(0:nl(1)+1,0:nl(2)+1,0:nl(3)+1) :: vof1,vof2,vof3
  Real(sp) :: norm(3)
  Real(sp) :: f_block(3,3,3)
  ! Real(sp) :: FloodSZ_Backward, FloodSZ_Forward
  Real(sp) :: EPSC = 1.0e-12
  Integer :: nexch(2)
  nexch = nl(1:2)

  ! default, f = 0
  vof1 = 0.0_sp
  vof2 = 0.0_sp
  vof3 = 0.0_sp

  if (dir.eq.1)  then
    ii  = 1; ic1 = 1; ic2 = 2; ic3 = 3
  else
    ii = 0
  end if
  if (dir.eq.2)  then
    jj  = 1; ic1 = 2; ic2 = 3; ic3 = 1
  else; jj = 0
  end if
  if (dir.eq.3)  then
    kk  = 1; ic1 = 3; ic2 = 1; ic3 = 2
  else; kk = 0
  end if

  !***
  do k=1,nl(3)
    do j=1,nl(2)
      do i=1,nl(1)
        a1 = us(i,j,k) *dt/dl(1)
        a2 = us(i+ii,j+jj,k+kk) *dt/dl(dir)

        ! f = 1
        if (f(i,j,k) .GE. 1.0_sp-epsc) then
          vof1(i,j,k) = DMAX1(-a1,0.0_sp)
          vof2(i,j,k) = 1.0_sp - DMAX1(a1,0.0_sp) + DMIN1(a2,0.0_sp)
          vof3(i,j,k) = DMAX1(a2,0.0_sp)

        ! 0 < f < 1
        else if (f(i,j,k) .GT. 0.0_sp) then
          !*(1)* normal vector
          f_block = f(i-1:i+1,j-1:j+1,k-1:k+1)
          Call NormMYCS(f_block, norm)
          Call Normalization1(norm)
          !*(2) get alpha;
          alpha = FloodSZ_Backward(norm,f(i,j,k))
          norm(dir) = norm(dir)/(1.0_sp - a1 + a2)
          alpha = alpha + norm(dir)*a1
          !*(3) get fluxes
          x0=0.0_sp; deltax=1.0_sp
          if (a1 .LT. 0.0_sp) then
            x0(dir)=a1;
            deltax(dir)=-a1
            vof1(i,j,k) = FloodSZ_Forward(norm,alpha,x0,deltax)
          end if
          if (a2 .GT. 0.0_sp) then
            x0(dir)=1.0_sp;
            deltax(dir)=a2
            vof3(i,j,k) = FloodSZ_Forward(norm,alpha,x0,deltax)
          end if
          x0(dir) = DMAX1(a1,0.0_sp)
          deltax(dir) = 1.0_sp - x0(dir) + DMIN1(0.0_sp,a2)
          vof2(i,j,k) = FloodSZ_Forward(norm,alpha,x0,deltax)
        endif
      enddo
    enddo
  enddo

  ! apply proper boundary conditions to vof1, vof2, vof3
  call phi_bc%setBCS(vof1)
  call phi_bc%setBCS(vof2)
  call phi_bc%setBCS(vof3)
  !new values of c and  clip it: 0. <= c <= 1.
  do k=1,nl(3)
    do j=1,nl(2)
      do i=1,nl(1)
        f(i,j,k) = vof3(i-ii,j-jj,k-kk) + vof2(i,j,k) + vof1(i+ii,j+jj,k+kk)
        if (f(i,j,k) < EPSC) then
          f(i,j,k) = 0.0_sp
        elseif ( f(i,j,k) >  (1.d0 - EPSC)) then
          f(i,j,k) = 1.0_sp
        endif
      enddo
    enddo
  enddo
  ! apply proper boundary conditions to c
  call phi_bc%setBCS(f)
  !***
  return
End Subroutine AdvCIAM

!=======================================================
! Advection algorithm of Dick-Weymouth PLIC
!=======================================================
Subroutine AdvWY(us, f, nl, dl, dt, dir)
  Implicit None
  Integer :: dir
  Real(sp) :: dt
  Integer,Intent(In)     :: nl(3)
  Real(sp),Intent(In)    :: dl(3)
  Real(sp),Intent(In)    :: us(0:nl(1)+1,0:nl(2)+1,0:nl(3)+1)
  Real(sp),Intent(InOut) :: f(0:nl(1)+1,0:nl(2)+1,0:nl(3)+1)
  Integer  :: i, j, k, ii, jj, kk, ic1, ic2, ic3
  Real(sp) :: a1,a2,alpha
  Real(sp) :: x0(3), deltax(3)
  Real(sp),dimension(0:nl(1)+1,0:nl(2)+1,0:nl(3)+1) :: vof1,vof3
  Real(sp) :: norm(3)
  Real(sp) :: f_block(3,3,3)
  ! Real(sp) :: FloodSZ_Backward, FloodSZ_Forward
  Real(sp) :: EPSC = 1.0e-12
  Integer :: nexch(2)
  nexch = nl(1:2)

  ! default, f = 0
  vof1 = 0.0_sp
  vof3 = 0.0_sp

  if (dir.eq.1)  then
    ii  = 1; ic1 = 1; ic2 = 2; ic3 = 3
  else
    ii = 0
  end if
  if (dir.eq.2)  then
    jj  = 1; ic1 = 2; ic2 = 3; ic3 = 1
  else; jj = 0
  end if
  if (dir.eq.3)  then
    kk  = 1; ic1 = 3; ic2 = 1; ic3 = 2
  else; kk = 0
  end if

  !***
  do k=1,nl(3)
    do j=1,nl(2)
      do i=1,nl(1)
        a1 = us(i,j,k) *dt/dl(1)
        a2 = us(i+ii,j+jj,k+kk) *dt/dl(dir)

        ! f = 1
        if (f(i,j,k) .GE. 1.0_sp-epsc) then
          vof1(i,j,k) = DMAX1(-a1,0.0_sp)
          vof3(i,j,k) = DMAX1(a2,0.0_sp)

        ! 0 < f < 1
        else if (f(i,j,k) .GT. 0.0_sp) then
          !*(1)* normal vector
          f_block = f(i-1:i+1,j-1:j+1,k-1:k+1)
          Call NormMYCS(f_block, norm)
          Call Normalization1(norm)
          !*(2) get alpha;
          alpha = FloodSZ_Backward(norm,f(i,j,k))
          !*(3) get fluxes
          x0=0.0_sp; deltax=1.0_sp
          if (a1 .LT. 0.0_sp) then
            deltax(dir)=-a1
            vof1(i,j,k) = FloodSZ_Forward(norm,alpha,x0,deltax)
          end if
          if (a2 .GT. 0.0_sp) then
            x0(dir)=1.0_sp-a2;
            deltax(dir)=a2
            vof3(i,j,k) = FloodSZ_Forward(norm,alpha,x0,deltax)
          end if
        endif
      enddo
    enddo
  enddo

  ! apply proper boundary conditions to vof1, vof3
  call phi_bc%setBCS(vof1)
  call phi_bc%setBCS(vof3)
  !new values of f and  clip it: 0. <= f <= 1.
  do k=1,nl(3)
    do j=1,nl(2)
      do i=1,nl(1)
        a1 = us(i,j,k)*dt/dl(dir)
        a2 = us(i+ii,j+jj,k+kk)*dt/dl(dir)
        f(i,j,k) = f(i,j,k) - (vof3(i,j,k) - vof1(i+ii,j+jj,k+kk)) + & 
            (vof3(i-ii,j-jj,k-kk) - vof1(i,j,k)) !+ vof2(i,j,k)*(a2-a1);
        if (f(i,j,k) < EPSC) then
          f(i,j,k) = 0.d0
        elseif (f(i,j,k) >  (1.0_sp - EPSC)) then
          f(i,j,k) = 1.0_sp
        endif
      enddo
    enddo
  enddo
  ! apply proper boundary conditions to c
  call phi_bc%setBCS(f)
  !***
  return
End Subroutine AdvWY

!=======================================================
! Advection algorithm of CIAM MOF
!=======================================================
Subroutine AdvCIAM_MOF(us, cx, cy, cz, f, nl, dl, dt, dir)
  Implicit None
  Integer :: dir
  Real(sp) :: dt
  Integer,Intent(In)     :: nl(3)
  Real(sp),Intent(In)    :: dl(3)
  Real(sp),Intent(In)    :: us(0:nl(1)+1,0:nl(2)+1,0:nl(3)+1)
  Real(sp),Intent(InOut) :: cx(0:nl(1)+1,0:nl(2)+1,0:nl(3)+1)
  Real(sp),Intent(InOut) :: cy(0:nl(1)+1,0:nl(2)+1,0:nl(3)+1)
  Real(sp),Intent(InOut) :: cz(0:nl(1)+1,0:nl(2)+1,0:nl(3)+1)
  Real(sp),Intent(InOut) :: f(0:nl(1)+1,0:nl(2)+1,0:nl(3)+1)
  Integer  :: i, j, k, ii, jj, kk, ic1, ic2, ic3
  Real(sp) :: a1,a2,alpha
  Real(sp) :: mx, my, mz
  Real(sp) :: x0(3), deltax(3)
  Real(sp) :: c1xyz(3), c2xyz(3), c3xyz(3)
  Real(sp),dimension(0:nl(1)+1,0:nl(2)+1,0:nl(3)+1) :: vof1,vof2,vof3
  Real(sp),dimension(0:nl(1)+1,0:nl(2)+1,0:nl(3)+1) :: c1x, c1y, c1z
  Real(sp),dimension(0:nl(1)+1,0:nl(2)+1,0:nl(3)+1) :: c2x, c2y, c2z
  Real(sp),dimension(0:nl(1)+1,0:nl(2)+1,0:nl(3)+1) :: c3x, c3y, c3z
  Real(sp) :: f_block(3,3,3), init_norm(3)
  Real(sp) :: norm(3)
  Real(sp) :: c3(3)
  ! Real(sp) :: FloodSZ_Backward
  Real(sp) :: EPSC = 1.0e-12
  Integer :: nexch(2)
  nexch = nl(1:2)

  ! default, f = 0
  vof1 = 0.0_sp
  vof2 = 0.0_sp
  vof3 = 0.0_sp

  if (dir.eq.1)  then
    ii  = 1; ic1 = 1; ic2 = 2; ic3 = 3
  else
    ii = 0
  end if
  if (dir.eq.2)  then
    jj  = 1; ic1 = 2; ic2 = 3; ic3 = 1
  else; jj = 0
  end if
  if (dir.eq.3)  then
    kk  = 1; ic1 = 3; ic2 = 1; ic3 = 2
  else; kk = 0
  end if

  !***
  do k=1,nl(3)
    do j=1,nl(2)
      do i=1,nl(1)
        a1 = us(i,j,k) *dt/dl(1)
        a2 = us(i+ii,j+jj,k+kk) *dt/dl(dir)

        c1xyz = 0.0_sp
        c2xyz = 0.0_sp
        c3xyz = 0.0_sp
        ! f = 1
        if (f(i,j,k) .GE. 1.0_sp-epsc) then
          vof1(i,j,k) = DMAX1(-a1,0.0_sp)
          vof2(i,j,k) = 1.0_sp - DMAX1(a1,0.0_sp) + DMIN1(a2,0.0_sp)
          vof3(i,j,k) = DMAX1(a2,0.0_sp)
          c1xyz = 0.5_sp
          c2xyz = 0.5_sp
          c3xyz = 0.5_sp
          c1xyz(dir) = - DMAX1(-a1/2.0_sp,0.0_sp)
          c2xyz(dir) = 0.5_sp + DMAX1(a1/2.0_sp,0.0_sp) - DMIN1(a2/2.0_sp,0.0_sp)
          c3xyz(dir) = 1.0_sp + DMAX1(a2/2.0_sp,0.0_sp)

        ! 0 < f < 1
        else if (f(i,j,k) .GT. 0.0_sp) then
          !*(1)* normal vector
          c3(1) = cx(i,j,k)
          c3(2) = cy(i,j,k)
          c3(3) = cz(i,j,k)
          f_block = f(i-1:i+1,j-1:j+1,k-1:k+1)
          Call NormMYCS(f_block, init_norm)
          ! Call NormMOF(f(i,j,k), c3, norm, init_norm=init_norm)
          Call NormMOF(f(i,j,k), c3, norm)
          Call Normalization1(norm)
          !*(2) get alpha;
          alpha = FloodSZ_Backward(norm,f(i,j,k))
          norm(dir) = norm(dir)/(1.0_sp - a1 + a2)
          alpha = alpha + norm(dir)*a1
          !*(3) get fluxes
          x0=0.0_sp; deltax=1.0_sp
          if (a1 .LT. 0.0_sp) then
            x0(dir)=a1;
            deltax(dir)=-a1
            Call FloodSZ_ForwardC(norm,alpha,x0,deltax,vof1(i,j,k),c1xyz)
          end if
          if (a2 .GT. 0.0_sp) then
            x0(dir)=1.0_sp;
            deltax(dir)=a2
            Call FloodSZ_ForwardC(norm,alpha,x0,deltax,vof3(i,j,k),c3xyz)
          end if
          x0(dir) = DMAX1(a1,0.0_sp)
          deltax(dir) = 1.0_sp - x0(dir) + DMIN1(0.0_sp,a2)
          Call FloodSZ_ForwardC(norm,alpha,x0,deltax,vof2(i,j,k),c2xyz)
        endif
        ! Call Centroid_Eulerian_Adv(c1xyz, a1, a2, 0.0_sp, 1.0_sp, dir)
        ! Call Centroid_Eulerian_Adv(c2xyz, a1, a2, 0.0_sp, 1.0_sp, dir)
        ! Call Centroid_Eulerian_Adv(c3xyz, a1, a2, 0.0_sp, 1.0_sp, dir)
        ! If (vof1(i,j,k) .gt. epsc) Then
        !   Call Centroid_Lagrangian_Adv(c1xyz, -a1, -a2, 0.0_sp, 1.0_sp, dir)
        ! EndIf
        ! If (vof2(i,j,k) .gt. epsc) Then
        !   Call Centroid_Lagrangian_Adv(c2xyz, -a1, -a2, 0.0_sp, 1.0_sp, dir)
        ! EndIf
        ! If (vof3(i,j,k) .gt. epsc) Then
        !   Call Centroid_Lagrangian_Adv(c3xyz, -a1, -a2, 0.0_sp, 1.0_sp, dir)
        ! EndIf
        c1xyz(dir) = c1xyz(dir) + 1.0_sp
        c3xyz(dir) = c3xyz(dir) - 1.0_sp

        c1x(i,j,k) = c1xyz(1)
        c1y(i,j,k) = c1xyz(2)
        c1z(i,j,k) = c1xyz(3)
        c3x(i,j,k) = c3xyz(1)
        c3y(i,j,k) = c3xyz(2)
        c3z(i,j,k) = c3xyz(3)
        c2x(i,j,k) = c2xyz(1)
        c2y(i,j,k) = c2xyz(2)
        c2z(i,j,k) = c2xyz(3)

      enddo
    enddo
  enddo

  ! apply proper boundary conditions to vof1, vof2, vof3
  call phi_bc%setBCS(vof1)
  call phi_bc%setBCS(vof2)
  call phi_bc%setBCS(vof3)
  call phi_bc%setBCS(c1x)
  call phi_bc%setBCS(c2x)
  call phi_bc%setBCS(c3x)
  call phi_bc%setBCS(c1y)
  call phi_bc%setBCS(c2y)
  call phi_bc%setBCS(c3y)
  call phi_bc%setBCS(c1z)
  call phi_bc%setBCS(c2z)
  call phi_bc%setBCS(c3z)
  !new values of c and  clip it: 0. <= c <= 1.
  do k=1,nl(3)
    do j=1,nl(2)
      do i=1,nl(1)
        mx = vof3(i-ii,j-jj,k-kk) * c3x(i-ii,j-jj,k-kk) + &
            vof2(i,j,k) * c2x(i,j,k) + &
            vof1(i+ii,j+jj,k+kk) * c1x(i+ii,j+jj,k+kk)
        my = vof3(i-ii,j-jj,k-kk) * c3y(i-ii,j-jj,k-kk) +&
            vof2(i,j,k) * c2y(i,j,k) + &
            vof1(i+ii,j+jj,k+kk) * c1y(i+ii,j+jj,k+kk)
        ! mz = vof3(i-ii,j-jj,k-kk) * c3z(i-ii,j-j,k-kk) +&
        !     vof2(i,j,k) * c2z(i,j,k) + &
        !     vof1(i+ii,j+jj,k+kk) * c1z(i+ii,j+jj,k+kk)
        mz = vof3(i-ii,j-jj,k-kk) * c3z(i-ii,j-jj,k-kk) +&
            vof2(i,j,k) * c2z(i,j,k) + &
            vof1(i+ii,j+jj,k+kk) * c1z(i+ii,j+jj,k+kk)

        f(i,j,k) = vof3(i-ii,j-jj,k-kk) + &
            vof2(i,j,k) + &
            vof1(i+ii,j+jj,k+kk)
        cx(i,j,k) = mx / f(i,j,k)
        cy(i,j,k) = my / f(i,j,k)
        cz(i,j,k) = mz / f(i,j,k)
        if (f(i,j,k) < EPSC) then
          f(i,j,k) = 0.0_sp
          cx(i,j,k) = 0.0_sp
          cy(i,j,k) = 0.0_sp
          cz(i,j,k) = 0.0_sp
        elseif ( f(i,j,k) >  (1.d0 - EPSC)) then
          f(i,j,k) = 1.0_sp
          cx(i,j,k) = 0.5_sp
          cy(i,j,k) = 0.5_sp
          cz(i,j,k) = 0.5_sp
        endif
        if (cx(i,j,k) < 0.0_sp) then
          cx(i,j,k) = 0.0_sp
          print *, 'cx<1', i,j,k
        elseif (cx(i,j,k) >  (1.0_sp )) then
          cx(i,j,k) = 1.0_sp
          print *, 'cx>1', i,j,k
        endif
        if (cy(i,j,k) < 0.0_sp) then
          cy(i,j,k) = 0.0_sp
          print *, 'cy<1', i,j,k
        elseif (cy(i,j,k) >  (1.0_sp )) then
          cy(i,j,k) = 1.0_sp
          print *, 'cy>1', i,j,k
        endif
        if (cz(i,j,k) < 0.0_sp) then
          cz(i,j,k) = 0.0_sp
          print *, 'cz<0', i,j,k
        elseif (cz(i,j,k) >  (1.0_sp )) then
          cz(i,j,k) = 1.0_sp
          print *, 'cz>1', i,j,k
        endif
      enddo
    enddo
  enddo

  ! apply proper boundary conditions to c
  call phi_bc%setBCS(f)
  call phi_bc%setBCS(cx)
  call phi_bc%setBCS(cy)
  call phi_bc%setBCS(cz)
  !***
  return
End Subroutine AdvCIAM_MOF

Subroutine AdvWY_MOF(us, cx, cy, cz, f, nl, dl, dt, dir)
  Implicit None
  Integer :: dir
  Real(sp) :: dt
  Integer,Intent(In)     :: nl(3)
  Real(sp),Intent(In)    :: dl(3)
  Real(sp),Intent(In)    :: us(0:nl(1)+1,0:nl(2)+1,0:nl(3)+1)
  Real(sp),Intent(InOut) :: f(0:nl(1)+1,0:nl(2)+1,0:nl(3)+1)
  Real(sp),Intent(InOut) :: cx(0:nl(1)+1,0:nl(2)+1,0:nl(3)+1)
  Real(sp),Intent(InOut) :: cy(0:nl(1)+1,0:nl(2)+1,0:nl(3)+1)
  Real(sp),Intent(InOut) :: cz(0:nl(1)+1,0:nl(2)+1,0:nl(3)+1)
  Integer  :: i, j, k, ii, jj, kk, ic1, ic2, ic3
  Real(sp) :: a1,a2,alpha
  Real(sp) :: mx, my, mz
  Real(sp) :: x0(3), deltax(3)
  Real(sp) :: c1xyz(3), c2xyz(3), c3xyz(3)
  Real(sp),dimension(0:nl(1)+1,0:nl(2)+1,0:nl(3)+1) :: vof1,vof3
  Real(sp),dimension(0:nl(1)+1,0:nl(2)+1,0:nl(3)+1) :: c1x, c1y, c1z
  Real(sp),dimension(0:nl(1)+1,0:nl(2)+1,0:nl(3)+1) :: c2x, c2y, c2z
  Real(sp),dimension(0:nl(1)+1,0:nl(2)+1,0:nl(3)+1) :: c3x, c3y, c3z
  Real(sp) :: norm(3)
  Real(sp) :: c3(3)
  ! Real(sp) :: FloodSZ_Backward
  Real(sp) :: EPSC = 1.0e-12
  Integer :: nexch(2)
  nexch = nl(1:2)

  ! default, f = 0
  vof1 = 0.0_sp
  vof3 = 0.0_sp

  if (dir.eq.1)  then
    ii  = 1; ic1 = 1; ic2 = 2; ic3 = 3
  else
    ii = 0
  end if
  if (dir.eq.2)  then
    jj  = 1; ic1 = 2; ic2 = 3; ic3 = 1
  else; jj = 0
  end if
  if (dir.eq.3)  then
    kk  = 1; ic1 = 3; ic2 = 1; ic3 = 2
  else; kk = 0
  end if

  !***
  do k=1,nl(3)
    do j=1,nl(2)
      do i=1,nl(1)
        a1 = us(i,j,k) *dt/dl(1)
        a2 = us(i+ii,j+jj,k+kk) *dt/dl(dir)

        c1xyz = 0.0_sp
        c3xyz = 0.0_sp
        ! f = 1
        if (f(i,j,k) .GE. 1.0_sp-epsc) then
          vof1(i,j,k) = DMAX1(-a1,0.0_sp)
          vof3(i,j,k) = DMAX1(a2,0.0_sp)
          c1xyz = 0.5_sp
          c3xyz = 0.5_sp
          c1xyz(dir) = DMAX1(-a1/2,0.0_sp)
          c3xyz(dir) = 1.0_sp - DMAX1(a2/2.0_sp,0.0_sp)

        ! 0 < f < 1
        else if (f(i,j,k) .GT. 0.0_sp) then
          !*(1)* normal vector
          c3(1)  = cx(i,j,k)
          c3(2)  = cy(i,j,k)
          c3(3)  = cz(i,j,k)
          Call NormMOF(f(i,j,k), c3, norm)
          Call Normalization1(norm)
          ! f_block = f(i-1:i+1,j-1:j+1,k-1:k+1)
          ! Call NormMYCS(f_block, norm)
          !*(2) get alpha;
          alpha = FloodSZ_Backward(norm,f(i,j,k))
          !*(3) get fluxes
          x0=0.0_sp; deltax=1.0_sp
          if (a1 .LT. 0.0_sp) then
            deltax(dir)=-a1
            Call FloodSZ_ForwardC(norm,alpha,x0,deltax,vof1(i,j,k),c1xyz)
          end if
          if (a2 .GT. 0.0_sp) then
            x0(dir)=1.0_sp-a2;
            deltax(dir)=a2
            Call FloodSZ_ForwardC(norm,alpha,x0,deltax,vof3(i,j,k),c3xyz)
          end if
        endif
        c2xyz(1) = cx(i,j,k); c2xyz(2) = cy(i,j,k); c2xyz(3) = cz(i,j,k);
        ! Call Centroid_Eulerian_Adv(c1xyz, a1, a2, 0.0_sp, 1.0_sp, dir)
        ! Call Centroid_Eulerian_Adv(c2xyz, a1, a2, 0.0_sp, 1.0_sp, dir)
        ! Call Centroid_Eulerian_Adv(c3xyz, a1, a2, 0.0_sp, 1.0_sp, dir)
        Call Centroid_Lagrangian_Adv(c1xyz, a1, a2, 0.0_sp, 1.0_sp, dir)
        Call Centroid_Lagrangian_Adv(c2xyz, a1, a2, 0.0_sp, 1.0_sp, dir)
        Call Centroid_Lagrangian_Adv(c3xyz, a1, a2, 0.0_sp, 1.0_sp, dir)
        ! Call Visual3DContour(f1=cz)
        c1xyz(dir)  = c1xyz(dir) + 1.0_sp
        c3xyz(dir)  = c3xyz(dir) - 1.0_sp
        c1x(i,j,k) = c1xyz(1)
        c1y(i,j,k) = c1xyz(2)
        c1z(i,j,k) = c1xyz(3)
        c3x(i,j,k) = c3xyz(1)
        c3y(i,j,k) = c3xyz(2)
        c3z(i,j,k) = c3xyz(3)
        c2x(i,j,k) = c2xyz(1)
        c2y(i,j,k) = c2xyz(2)
        c2z(i,j,k) = c2xyz(3)
      enddo
    enddo
  enddo

  ! apply proper boundary conditions to vof1, vof3
  call phi_bc%setBCS(vof1)
  call phi_bc%setBCS(vof3)
  call phi_bc%setBCS(c1x)
  call phi_bc%setBCS(c2x)
  call phi_bc%setBCS(c3x)
  call phi_bc%setBCS(c1y)
  call phi_bc%setBCS(c2y)
  call phi_bc%setBCS(c3y)
  call phi_bc%setBCS(c1z)
  call phi_bc%setBCS(c2z)
  call phi_bc%setBCS(c3z)
  !new values of f and  clip it: 0. <= f <= 1.
  print * , f(6,8,8), cx(6,8,8), cy(6,8,8), cz(6,8,8)
  do k=1,nl(3)
    do j=1,nl(2)
      do i=1,nl(1)
        mx = f(i,j,k) * c2x(i,j,k) &
            - vof1(i,j,k) * c1x(i,j,k) &
            - vof3(i,j,k) * c3x(i,j,k) &
            + vof1(i+ii,j+jj,k+kk) * c1x(i+ii,j+jj,k+kk) &
            + vof3(i-ii,j-jj,k-kk) * c3x(i-ii,j-jj,k-kk)
        my = f(i,j,k) * c2y(i,j,k) &
            - vof1(i,j,k) * c1y(i,j,k) &
            - vof3(i,j,k) * c3y(i,j,k) &
            + vof1(i+ii,j+jj,k+kk) * c1y(i+ii,j+jj,k+kk) &
            + vof3(i-ii,j-jj,k-kk) * c3y(i-ii,j-jj,k-kk)
        mz = f(i,j,k) * c2z(i,j,k) &
            - vof1(i,j,k) * c1z(i,j,k) &
            - vof3(i,j,k) * c3z(i,j,k) &
            + vof1(i+ii,j+jj,k+kk) * c1z(i+ii,j+jj,k+kk) &
            + vof3(i-ii,j-jj,k-kk) * c3z(i-ii,j-jj,k-kk)
        f(i,j,k) = f(i,j,k) &
            - vof1(i,j,k) &
            - vof3(i,j,k) &
            + vof1(i+ii,j+jj,k+kk) &
            + vof3(i-ii,j-jj,k-kk)
        cx(i,j,k) = mx / f(i,j,k)
        cy(i,j,k) = my / f(i,j,k)
        cz(i,j,k) = mz / f(i,j,k)
        if (f(i,j,k) < EPSC) then
          f(i,j,k) = 0.0_sp
          cx(i,j,k) = 0.0_sp
          cy(i,j,k) = 0.0_sp
          cz(i,j,k) = 0.0_sp
        elseif (f(i,j,k) >  (1.0_sp - EPSC)) then
          f(i,j,k) = 1.0_sp
          cx(i,j,k) = 0.5_sp
          cy(i,j,k) = 0.5_sp
          cz(i,j,k) = 0.5_sp
        endif
        if (cx(i,j,k) < 0.0_sp) then
          cx(i,j,k) = 0.0_sp
          print *, 'cx<1', i,j,k
        elseif (cx(i,j,k) >  (1.0_sp )) then
          cx(i,j,k) = 1.0_sp
          print *, 'cx>1', i,j,k
        endif
        if (cy(i,j,k) < 0.0_sp) then
          cy(i,j,k) = 0.0_sp
          print *, 'cy<1', i,j,k
        elseif (cy(i,j,k) >  (1.0_sp )) then
          cy(i,j,k) = 1.0_sp
          print *, 'cy>1', i,j,k
        endif
        if (cz(i,j,k) < 0.0_sp) then
          cz(i,j,k) = 0.0_sp
          print *, 'cz<0', i,j,k
        elseif (cz(i,j,k) >  (1.0_sp )) then
          cz(i,j,k) = 1.0_sp
          print *, 'cz>1', i,j,k
        endif
      enddo
    enddo
  enddo

  ! Call Visual3DContour(f1=cx)
  ! Call Visual3DContour(f1=cy)
  ! Call Visual3DContour(f1=cz)

  print * , f(6,8,8), cx(6,8,8), cy(6,8,8), cz(6,8,8)
  print * , f(11,8,8), cx(11,8,8), cy(11,8,8), cz(11,8,8)
  print * , c1x(10,8,8), c1y(10,8,8), c1z(10,8,8)
  print * , c3x(10,8,8), c3y(10,8,8), c3z(10,8,8)
  ! apply proper boundary conditions to c
  ! apply proper boundary conditions to c
  call phi_bc%setBCS(f)
  call phi_bc%setBCS(cx)
  call phi_bc%setBCS(cy)
  call phi_bc%setBCS(cz)
  !***
  return
End Subroutine AdvWY_MOF

Subroutine Centroid_Lagrangian_Adv(c, ul, ur, xl, xr, dir)
  Implicit None
  Real(sp), Intent(Inout) :: c(3)
  Real(sp), Intent(In) :: ul, ur
  Real(sp), Intent(In) :: xl, xr
  Integer,  Intent(In) :: dir

  c(dir) = (1.0_sp + ur - ul) * c(dir) - (ur * xl - ul * xr)

End Subroutine Centroid_Lagrangian_Adv

Subroutine Centroid_Eulerian_Adv(c, ul, ur, xl, xr, dir)
  Implicit None
  Real(sp), Intent(Inout) :: c(3)
  Real(sp), Intent(In) :: ul, ur
  Real(sp), Intent(In) :: xl, xr
  Integer,  Intent(In) :: dir

  c(dir) =  (c(dir) - (ur * xl - ul * xr)) / (1.0_sp - ur + ul)

End Subroutine Centroid_Eulerian_Adv

End Module ModVOF
