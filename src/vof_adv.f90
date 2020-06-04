!=======================================================
! Advection algorithm of SZ PLIC
!=======================================================
SUBROUTINE AdvSZ(us, c, nl, dl, dt, dir)
  Use ModGlobal, only : sp
  Implicit None
  Integer :: dir
  Real(sp) :: dt
  Integer,Intent(In)     :: nl(3)
  Real(sp),Intent(In)    :: dl(3)
  Real(sp),Intent(In)    :: us(0:nl(1)+1,0:nl(2)+1,0:nl(3)+1)
  Real(sp),Intent(InOut) :: c(0:nl(1)+1,0:nl(2)+1,0:nl(3)+1)
  Real(sp) :: mx,my,mz,mm1,mm2
  Integer  :: i,j,k,invx,invy,invz
  Real(sp) :: a1,a2,alpha
  Real(sp),dimension(0:nl(1)+1,0:nl(2)+1,0:nl(3)+1) :: vof1,vof2,vof3
  Real(sp) :: norm(3)
  Real(sp) :: c_block(3,3,3)
  Real(sp) :: FloodAPPLIC_Backward, FloodAPPLIC_Forward

  vof1 = 0.0d0
  vof2 = 0.0d0
  vof3 = 0.0d0
  !***
  do k=1,nl(3)
    do j=1,nl(2)
      do i=1,nl(1)
        a1 = us(i,j,k) *dt/dl(1)
        if (dir.eq.1) then
          a2 = us(i+1,j,k) *dt/dl(2)
        elseif (dir.eq.2) then
          a2 = us(i,j+1,k) *dt/dl(3)
        elseif (dir.eq.3) then
          a2 = us(i,j,k+1) *dt/dl(4)
        endif

        !***
        !     3 cases: 1: DEFAULT (c=0. and fluxes=0.); 2: c=1.; 3:c>0.
        !***

        if (c(i,j,k) .EQ. 1.0d0) then
          vof1(i,j,k) = DMAX1(-a1,0.d0)
          vof2(i,j,k) = 1.d0 - DMAX1(a1,0.d0) + DMIN1(a2,0.d0)
          vof3(i,j,k) = DMAX1(a2,0.d0)

        else if (c(i,j,k) .GT. 0.d0) then
          !***
          !     (1) normal vector: mx,my,mz; (2) mx,my,mz>0. and mx+my+mz = 1.;
          !     (3) get alpha;               (4) back to original plane;
          !     (5) lagrangian advection;    (6) get fluxes

          !*(1)*
          c_block = c(i-1:i+1,j-1:j+1,k-1:k+1)
          ! Call NormParkerYoungs(c, norm)
                  mm1 = c(i-1,j-1,k-1)+c(i-1,j-1,k+1)+c(i-1,j+1,k-1) & 
                      +c(i-1,j+1,k+1)+2.0d0*(c(i-1,j-1,k)+c(i-1,j+1,k)&
                      +c(i-1,j,k-1)+c(i-1,j,k+1))+4.0d0*c(i-1,j,k)
                  mm2 = c(i+1,j-1,k-1)+c(i+1,j-1,k+1)+c(i+1,j+1,k-1) &
                      +c(i+1,j+1,k+1)+2.0d0*(c(i+1,j-1,k)+c(i+1,j+1,k) &
                      +c(i+1,j,k-1)+c(i+1,j,k+1))+4.0d0*c(i+1,j,k)
                  mx = mm1 - mm2
                  
                  mm1 = c(i-1,j-1,k-1)+c(i-1,j-1,k+1)+c(i+1,j-1,k-1) &
                      +c(i+1,j-1,k+1)+2.0d0*(c(i-1,j-1,k)+c(i+1,j-1,k) &
                      +c(i,j-1,k-1)+c(i,j-1,k+1))+4.0d0*c(i,j-1,k)
                  mm2 = c(i-1,j+1,k-1)+c(i-1,j+1,k+1)+c(i+1,j+1,k-1) &
                      +c(i+1,j+1,k+1)+2.0d0*(c(i-1,j+1,k)+c(i+1,j+1,k) &
                      +c(i,j+1,k-1)+c(i,j+1,k+1))+4.0d0*c(i,j+1,k)
                  my = mm1 - mm2
                  
                  mm1 = c(i-1,j-1,k-1)+c(i-1,j+1,k-1)+c(i+1,j-1,k-1) &
                      +c(i+1,j+1,k-1)+2.0d0*(c(i-1,j,k-1)+c(i+1,j,k-1) &
                      +c(i,j-1,k-1)+c(i,j+1,k-1))+4.0d0*c(i,j,k-1)
                  mm2 = c(i-1,j-1,k+1)+c(i-1,j+1,k+1)+c(i+1,j-1,k+1) &
                      +c(i+1,j+1,k+1)+2.0d0*(c(i-1,j,k+1)+c(i+1,j,k+1) &
                      +c(i,j-1,k+1)+c(i,j+1,k+1))+4.0d0*c(i,j,k+1)
                  mz = mm1 - mm2
          ! mx = norm(1)
          ! my = norm(2)
          ! mz = norm(3)

          !*(2)*
          invx = 1
          invy = 1
          invz = 1
          if (mx .LT. 0.0d0) then
            mx = -mx
            invx = -1
          endif
          if (my .LT. 0.0d0) then
            my = -my
            invy = -1
          endif
          if (mz .LT. 0.0d0) then
            mz = -mz
            invz = -1
          endif
          mm2 = mx+my+mz
          mx = mx/mm2
          my = my/mm2
          mz = mz/mm2

          !*(3)*  
          alpha = FloodAPPLIC_Backward(mx,my,mz,c(i,j,k))

          !*(4)*  
          mx = invx*mx
          my = invy*my
          mz = invz*mz
          alpha = alpha + DMIN1(0.d0,mx) + DMIN1(0.d0,my) + &
          &                 DMIN1(0.d0,mz)

          !*(5)*
          mm1 = DMAX1(a1,0.0d0)
          mm2 = 1.d0 - mm1 + DMIN1(0.d0,a2)

          if (dir.eq.1) then
            mx = mx/(1.0d0 - a1 + a2)
            alpha = alpha + mx*a1
            if (a1 .LT. 0.d0) vof1(i,j,k) = FloodAPPLIC_Forward(mx,my,mz,alpha,a1  ,-a1)
            if (a2 .GT. 0.d0) vof3(i,j,k) = FloodAPPLIC_Forward(mx,my,mz,alpha,1.d0,a2)
            vof2(i,j,k) = FloodAPPLIC_Forward(mx,my,mz,alpha,mm1,mm2)
          elseif (dir.eq.2) then
            my = my/(1.0d0 - a1 + a2)
            alpha = alpha + my*a1
            if (a1 .LT. 0.d0) vof1(i,j,k) = FloodAPPLIC_Forward(my,mz,mx,alpha,a1  ,-a1)
            if (a2 .GT. 0.d0) vof3(i,j,k) = FloodAPPLIC_Forward(my,mz,mx,alpha,1.d0,a2)
            vof2(i,j,k) = FloodAPPLIC_Forward(my,mz,mx,alpha,mm1,mm2)
          elseif (dir.eq.3) then
            mz = mz/(1.0d0 - a1 + a2)
            alpha = alpha + mz*a1
            if (a1 .LT. 0.d0) vof1(i,j,k) = FloodAPPLIC_Forward(mz,mx,my,alpha,a1  ,-a1)
            if (a2 .GT. 0.d0) vof3(i,j,k) = FloodAPPLIC_Forward(mz,mx,my,alpha,1.d0,a2)
            vof2(i,j,k) = FloodAPPLIC_Forward(mz,mx,my,alpha,mm1,mm2)
          endif
        endif
      enddo
    enddo
  enddo
  !***
  !     (1) apply proper boundary conditions to fluxes
  !     (2) new values of c and  clip it: 0. <= c <= 1.
  !     (3) apply proper boundary conditions to c
  !*(1)*
  ! call bc_flux(vof1,vof3,us,nx,ny,nz,d)
 !*(2)*
  do k=1,nl(3)
    do j=1,nl(2)
      do i=1,nl(1)
        if (dir.eq.1) then
          c(i,j,k) = vof3(i-1,j,k) + vof2(i,j,k) + vof1(i+1,j,k)
        elseif (dir.eq.2) then
          c(i,j,k) = vof3(i,j-1,k) + vof2(i,j,k) + vof1(i,j+1,k)
        elseif (dir.eq.3) then
          c(i,j,k) = vof3(i,j,k-1) + vof2(i,j,k) + vof1(i,j,k+1)
        endif
        c(i,j,k) = DMAX1(0.0d0,DMIN1(1.0d0,c(i,j,k)))
      enddo
    enddo
  enddo
  !*(3)*
  ! call bc_c(c,nx,ny,nz)
  !***
  return
End Subroutine AdvSZ
