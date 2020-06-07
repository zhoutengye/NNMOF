# if defined NORM_PY
# define NormCalc NormParkerYoungs
# elif defined NORM_CS
# define NormCalc NormCS
# elif defined NORM_MYCS
# define NormCalc NormMYCS
# endif

# define FloodForward FloodSZ_Forward
# define FloodBackward FloodSZ_Backward

Module Mod_VOF
  Use ModGlobal, only : sp
  Use ModGlobal, only : updthalo
  Use ModGlobal, only : setBCS, phi_bc

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
  Real(sp) :: norm(3), abs_norm(3)
  Real(sp) :: f_block(3,3,3)
  Real(sp) :: FloodBackward, FloodForward
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
        if (f(i,j,k) .EQ. 1.0_sp) then
          vof1(i,j,k) = DMAX1(-a1,0.0_sp)
          vof2(i,j,k) = 1.0_sp - DMAX1(a1,0.0_sp) + DMIN1(a2,0.0_sp)
          vof3(i,j,k) = DMAX1(a2,0.0_sp)

        ! 0 < f < 1
        else if (f(i,j,k) .GT. 0.0_sp) then
          !*(1)* normal vector
          f_block = f(i-1:i+1,j-1:j+1,k-1:k+1)
          Call NormCalc(f_block, norm, abs_norm)
          !*(2) get alpha;
          alpha = FloodBackward(norm,f(i,j,k))
          norm(dir) = norm(dir)/(1.0_sp - a1 + a2)
          alpha = alpha + norm(dir)*a1
          !*(3) get fluxes
          x0=0.0_sp; deltax=1.0_sp
          if (a1 .LT. 0.0_sp) then
            x0(dir)=a1;
            deltax(dir)=-a1
            vof1(i,j,k) = FloodForward(norm,alpha,x0,deltax)
          end if
          if (a2 .GT. 0.0_sp) then
            x0(dir)=1.0_sp;
            deltax(dir)=a2
            vof3(i,j,k) = FloodForward(norm,alpha,x0,deltax)
          end if
          x0(dir) = DMAX1(a1,0.0_sp)
          deltax(dir) = 1.0_sp - x0(dir) + DMIN1(0.0_sp,a2)
          vof2(i,j,k) = FloodForward(norm,alpha,x0,deltax)
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
  Real(sp),dimension(0:nl(1)+1,0:nl(2)+1,0:nl(3)+1) :: vof1,vof2, vof3
  Real(sp) :: norm(3), abs_norm(3)
  Real(sp) :: f_block(3,3,3)
  Real(sp) :: FloodBackward, FloodForward
  Real(sp) :: EPSC = 1.0e-12
  Integer :: nexch(2)
  nexch = nl(1:2)

  ! default, f = 0
  vof1 = 0.0_sp
  vof2 = f
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
        if (f(i,j,k) .EQ. 1.0_sp) then
          vof1(i,j,k) = DMAX1(-a1,0.0_sp)
          vof3(i,j,k) = DMAX1(a2,0.0_sp)

        ! 0 < f < 1
        else if (f(i,j,k) .GT. 0.0_sp) then
          !*(1)* normal vector
          f_block = f(i-1:i+1,j-1:j+1,k-1:k+1)
          Call NormCalc(f_block, norm, abs_norm)
          !*(2) get alpha;
          alpha = FloodBackward(norm,f(i,j,k))
          norm(dir) = norm(dir)/(1.0_sp - a1 + a2)
          alpha = alpha + norm(dir)*a1
          !*(3) get fluxes
          x0=0.0_sp; deltax=1.0_sp
          if (a1 .LT. 0.0_sp) then
            deltax(dir)=-a1
            vof1(i,j,k) = FloodForward(norm,alpha,x0,deltax)
          end if
          if (a2 .GT. 0.0_sp) then
            x0(dir)=1.0_sp-a2;
            deltax(dir)=a2
            vof3(i,j,k) = FloodForward(norm,alpha,x0,deltax)
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
  Real(sp) :: norm(3), abs_norm(3)
  Real(sp) :: c3(3)
  Real(sp) :: FloodBackward!, FloodForward
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
        if (f(i,j,k) .EQ. 1.0_sp) then
          vof1(i,j,k) = DMAX1(-a1,0.0_sp)
          vof2(i,j,k) = 1.0_sp - DMAX1(a1,0.0_sp) + DMIN1(a2,0.0_sp)
          vof3(i,j,k) = DMAX1(a2,0.0_sp)
          c1xyz = 0.5
          c2xyz = 0.5
          c3xyz = 0.5
          c1xyz(dir) = DMAX1(-a1/2.0_sp,0.0_sp)
          c2xyz(dir) = 0.5_sp - DMAX1(-a1/2.0_sp,0.0_sp) - Dmin1(a2/2.0_sp,0.0_sp)
          c3xyz(dir) = DMAX1(a1/2.0_sp,0.0_sp)

        ! 0 < f < 1
        else if (f(i,j,k) .GT. 0.0_sp) then
          !*(1)* normal vector
          c3(1) = cx(i,j,k)
          c3(2) = cy(i,j,k)
          c3(3) = cz(i,j,k)
          Call NormMOF(f, c3, norm, abs_norm)
          !*(2) get alpha;
          alpha = FloodBackward(norm,f(i,j,k))
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
          c1x(i,j,k) = c1xyz(1);
          c1y(i,j,k) = c1xyz(2);
          c1z(i,j,k) = c1xyz(3);
          c2x(i,j,k) = c2xyz(1);
          c2y(i,j,k) = c2xyz(2);
          c2z(i,j,k) = c2xyz(3);
          c3x(i,j,k) = c3xyz(1);
          c3y(i,j,k) = c3xyz(2);
          c3z(i,j,k) = c3xyz(3);
        endif
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
        mx = vof3(i-ii,j-jj,k-kk) * c3x(i,j,k) &
            + vof2(i,j,k) * c2x(i,j,k) + vof1(i+ii,j+jj,k+kk) * c1x(i,j,k)
        my = vof3(i-ii,j-jj,k-kk) * c3y(i,j,k) &
            + vof2(i,j,k) * c2y(i,j,k) + vof1(i+ii,j+jj,k+kk) * c1y(i,j,k)
        mz = vof3(i-ii,j-jj,k-kk) * c3z(i,j,k) &
            + vof2(i,j,k) * c2z(i,j,k) + vof1(i+ii,j+jj,k+kk) * c1z(i,j,k)
        f(i,j,k) = vof3(i-ii,j-jj,k-kk) + vof2(i,j,k) + vof1(i+ii,j+jj,k+kk)
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
      enddo
    enddo
  enddo
  ! apply proper boundary conditions to c
  call phi_bc%setBCS(f)
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
  Real(sp) :: c1xyz(3), c3xyz(3)
  Real(sp),dimension(0:nl(1)+1,0:nl(2)+1,0:nl(3)+1) :: vof1,vof2, vof3
  Real(sp),dimension(0:nl(1)+1,0:nl(2)+1,0:nl(3)+1) :: c1x, c1y, c1z
  Real(sp),dimension(0:nl(1)+1,0:nl(2)+1,0:nl(3)+1) :: c3x, c3y, c3z
  Real(sp) :: norm(3), abs_norm(3)
  Real(sp) :: c3(3)
  Real(sp) :: FloodBackward!, FloodForward
  Real(sp) :: EPSC = 1.0e-12
  Integer :: nexch(2)
  nexch = nl(1:2)

  ! default, f = 0
  vof1 = 0.0_sp
  vof2 = f
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
        if (f(i,j,k) .EQ. 1.0_sp) then
          vof1(i,j,k) = DMAX1(-a1,0.0_sp)
          vof3(i,j,k) = DMAX1(a2,0.0_sp)
          c1xyz = 0.5_sp
          c3xyz = 0.5_sp
          c1xyz(dir) = DMAX1(-a1/2.0_sp,0.0_sp)
          c3xyz(dir) = DMAX1(a2/2.0_sp,0.0_sp)

        ! 0 < f < 1
        else if (f(i,j,k) .GT. 0.0_sp) then
          !*(1)* normal vector
          c3(1)  = cx(i,j,k)
          c3(2)  = cy(i,j,k)
          c3(3)  = cz(i,j,k)
          Call NormMOF(f(i,j,k), norm, abs_norm)
          !*(2) get alpha;
          alpha = FloodBackward(norm,f(i,j,k))
          norm(dir) = norm(dir)/(1.0_sp - a1 + a2)
          alpha = alpha + norm(dir)*a1
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
          c1x(i,j,k) = c1xyz(1);
          c1y(i,j,k) = c1xyz(2);
          c1z(i,j,k) = c1xyz(3);
          c3x(i,j,k) = c3xyz(1);
          c3y(i,j,k) = c3xyz(2);
          c3z(i,j,k) = c3xyz(3);
        endif
      enddo
    enddo
  enddo

  ! apply proper boundary conditions to vof1, vof3
  call phi_bc%setBCS(vof1)
  call phi_bc%setBCS(vof3)
  call phi_bc%setBCS(c1x)
  call phi_bc%setBCS(c3x)
  call phi_bc%setBCS(c1y)
  call phi_bc%setBCS(c3y)
  call phi_bc%setBCS(c1z)
  call phi_bc%setBCS(c3z)
  !new values of f and  clip it: 0. <= f <= 1.
  call phi_bc%setBCS(vof3)
  !new values of f and  clip it: 0. <= f <= 1.
  !new values of f and  clip it: 0. <= f <= 1.
  do k=1,nl(3)
    do j=1,nl(2)
      do i=1,nl(1)
        a1 = us(i,j,k)*dt/dl(dir)
        a2 = us(i+ii,j+jj,k+kk)*dt/dl(dir)
        mx = f(i,j,k) * cx(i,j,k) - &
            (vof3(i,j,k) * c3x(i,j,k) - vof1(i+ii,j+jj,k+kk) * c1x(i+ii,j+jj,k+kk)) + &
            (vof3(i-ii,j-jj,k-kk) * c3x(i-ii,j-jj,k-kk) - vof1(i,j,k) * c1x(i,j,k))
        my = f(i,j,k) * cy(i,j,k) - &
            (vof3(i,j,k) * c3y(i,j,k) - vof1(i+ii,j+jj,k+kk) * c1y(i+ii,j+jj,k+kk)) + &
            (vof3(i-ii,j-jj,k-kk) * c3y(i-ii,j-jj,k-kk) - vof1(i,j,k) * c1y(i,j,k))
        mz = f(i,j,k) * cz(i,j,k) - &
            (vof3(i,j,k) * c3z(i,j,k) - vof1(i+ii,j+jj,k+kk) * c1z(i+ii,j+jj,k+kk)) + &
            (vof3(i-ii,j-jj,k-kk) * c3z(i-ii,j-jj,k-kk) - vof1(i,j,k) * c1z(i,j,k))
        f(i,j,k) = f(i,j,k) - (vof3(i,j,k) - vof1(i+ii,j+jj,k+kk)) + & 
            (vof3(i-ii,j-jj,k-kk) - vof1(i,j,k)) !+ vof2(i,j,k)*(a2-a1);
        cx(i,j,k) = mx
        cy(i,j,k) = my
        cz(i,j,k) = mz
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
      enddo
    enddo
  enddo
  ! apply proper boundary conditions to c
  call phi_bc%setBCS(f)
  !***
  return
End Subroutine AdvWY_MOF

End Module MOD_VOF
