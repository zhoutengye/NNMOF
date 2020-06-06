# if defined NORM_PY
# define NormCalc NormParkerYoungs
# elif defined NORM_CS
# define NormCalc NormCS
# elif defined NORM_MYCS
# define NormCalc NormMYCS
# endif

# define FloodForward FloodSZ_Forward
# define FloodBackward FloodSZ_Backward

!=======================================================
! Advection algorithm of CIAM PLIC
!=======================================================
Subroutine AdvCIAM(us, f, nl, dl, dt, dir)
  Use ModGlobal, only : sp
  Use ModGlobal, only : setBCS, phi_bc
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
  Use ModGlobal, only : sp
  Use ModGlobal, only : setBCS, phi_bc
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
          alpha = FloodBackward(norm,c(i,j,k))
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

Subroutine AdvCIAM_MOF(us, f, c, nl, dl, dt, dir)
  Use ModGlobal, only : sp
  Use ModGlobal, only : setBCS, phi_bc
  Implicit None
  Integer :: dir
  Real(sp) :: dt
  Integer,Intent(In)     :: nl(3)
  Real(sp),Intent(In)    :: dl(3)
  Real(sp),Intent(In)    :: us(0:nl(1)+1,0:nl(2)+1,0:nl(3)+1)
  Real(sp),Intent(InOut) :: f(0:nl(1)+1,0:nl(2)+1,0:nl(3)+1)
  Real(sp),Intent(InOut) :: c(0:nl(1)+1,0:nl(2)+1,0:nl(3)+1)
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
            Call FloodSZ_forwardC(norm, alpha, x0, 
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
End Subroutine AdvCIAM_MOF
