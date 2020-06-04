!=======================================================
! Advection algorithm of SZ PLIC
!=======================================================
SUBROUTINE AdvSZ(us, c, nl, dl, dt, dir)
  Use ModGlobal, only : sp
  Use ModGlobal, only : setBCS, phi_bc
  Implicit None
  Integer :: dir
  Real(sp) :: dt
  Integer,Intent(In)     :: nl(3)
  Real(sp),Intent(In)    :: dl(3)
  Real(sp),Intent(In)    :: us(0:nl(1)+1,0:nl(2)+1,0:nl(3)+1)
  Real(sp),Intent(InOut) :: c(0:nl(1)+1,0:nl(2)+1,0:nl(3)+1)
  Real(sp) :: mm1,mm2
  Integer  :: i, j, k, ii, jj, kk, ic1, ic2, ic3
  Real(sp) :: a1,a2,alpha
  Real(sp),dimension(0:nl(1)+1,0:nl(2)+1,0:nl(3)+1) :: vof1,vof2,vof3
  Real(sp) :: norm(3), abs_norm(3)
  Real(sp) :: c_block(3,3,3)
  Real(sp) :: FloodAPPLIC_Backward, FloodAPPLIC_Forward
  Integer :: nexch(2)
  nexch = nl(1:2)

  vof1 = 0.0_sp
  vof2 = 0.0_sp
  vof3 = 0.0_sp

  if (dir.eq.1)  then
    ii  = 1
    ic1 = 1
    ic2 = 2
    ic3 = 3
  else
    ii = 0
  end if
  if (dir.eq.2)  then
    jj  = 1
    ic1 = 2
    ic2 = 3
    ic3 = 1
  else
    jj = 0
  end if
  if (dir.eq.3)  then
    kk  = 1
    ic1 = 3
    ic2 = 1
    ic3 = 2
  else
    kk = 0
  end if

  !***
  do k=1,nl(3)
    do j=1,nl(2)
      do i=1,nl(1)
        a1 = us(i,j,k) *dt/dl(1)
        a2 = us(i+ii,j+jj,k+kk) *dt/dl(dir)

        !***
        !     3 cases: 1: DEFAULT (c=0. and fluxes=0.); 2: c=1.; 3:c>0.
        !***

        if (c(i,j,k) .EQ. 1.0_sp) then
          vof1(i,j,k) = DMAX1(-a1,0.0_sp)
          vof2(i,j,k) = 1.0_sp - DMAX1(a1,0.0_sp) + DMIN1(a2,0.0_sp)
          vof3(i,j,k) = DMAX1(a2,0.0_sp)

        else if (c(i,j,k) .GT. 0.0_sp) then
          !***
          !     ;    (6) get fluxes

          !*(1)* normal vector
          c_block = c(i-1:i+1,j-1:j+1,k-1:k+1)
          Call NormParkerYoungs(c_block, norm)
          ! Call NormAPPLIC(c_block, norm)

          !*(2) mx,my,mz>0. and mx+my+mz = 1.;
          Call Normalization1(norm, abs_norm)

          !*(3) get alpha;
          alpha = FloodAPPLIC_Backward(abs_norm(1),abs_norm(2),abs_norm(3),c(i,j,k))

          !*(4) back to original plane;
          alpha = alpha + DMIN1(0.0_sp,norm(1)) + DMIN1(0.0_sp,norm(2)) + &
          &                 DMIN1(0.0_sp,norm(3))

          !*(5) semi-lagrangian advection
          mm1 = DMAX1(a1,0.0_sp)
          mm2 = 1.0_sp - mm1 + DMIN1(0.0_sp,a2)

          !*(6) get fluxes
          norm(dir) = norm(dir)/(1.0_sp - a1 + a2)
          alpha = alpha + norm(dir)*a1
          if (a1 .LT. 0.0_sp) vof1(i,j,k) = FloodAPPLIC_Forward(norm(ic1),norm(ic2),norm(ic3),alpha,a1  ,-a1)
          if (a2 .GT. 0.0_sp) vof3(i,j,k) = FloodAPPLIC_Forward(norm(ic1),norm(ic2),norm(ic3),alpha,1.0_sp,a2)
          vof2(i,j,k) = FloodAPPLIC_Forward(norm(ic1),norm(ic2),norm(ic3),alpha,mm1,mm2)
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
        c(i,j,k) = vof3(i-ii,j-jj,k-kk) + vof2(i,j,k) + vof1(i+ii,j+jj,k+kk)
        c(i,j,k) = DMAX1(0.0_sp,DMIN1(1.0_sp,c(i,j,k)))
      enddo
    enddo
  enddo
  ! apply proper boundary conditions to c
  call phi_bc%setBCS(c)
  !***
  return
End Subroutine AdvSZ
