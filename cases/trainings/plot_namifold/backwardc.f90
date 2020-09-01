 !=======================================================
  ! (2-1) Backward flooding algorithm of SZ finding alpha
  ! Calculate the cutting plane from the normal vector and
  ! volume function in UNIT CUBE. Calculate alpha in
  !    m1 x1 + m2 x2 + m3 x3 = alpha
  !
  ! Adopted from Paris simulator by Zaleski
  ! Used in:
  !     VOF reconstruction
  !     MOF reconstruction
  !-------------------------------------------------------
  ! Input:
  !      nr: normal vector(nx, ny, nz)
  !      cc: vof function
  ! Output:
  !      alpha
  !=======================================================
  Real(8) Function FloodSZ_Backward(nr,cc)
    Implicit None
    Integer, Parameter  :: sp = 8
    REAL(8), INTENT(IN):: nr(3),cc
    REAL(8) :: cch,c01,c02,c03,np1,np2,np3
    REAL(8) :: m1,m2,m3,m12,numer,denom,p,pst,q,arc,csarc, alpha
    REAL(8), PARAMETER :: athird=1.d0/3.d0,eps0=1.d-50
    INTRINSIC DABS,DMIN1,DMAX1,DSQRT,DACOS,DCOS

    np1 = DABS(nr(1))                                 ! need positive coefficients
    np2 = DABS(nr(2))
    np3 = DABS(nr(3))
    m1 = DMIN1(np1,np2)                              ! order positive coefficients
    m3 = DMAX1(np1,np2)
    if (np3 < m1) then
      m2 = m1
      m1 = np3
    else if (np3 >= m3) then
      m2 = m3
      m3 = np3
    else
      m2 = np3
    endif

    denom = DMAX1(6.d0*m1*m2*m3,eps0)                           
    cch = DMIN1(cc,1.d0-cc)                              ! limit to: 0 < cch < 1/2
    c01 = m1*m1*m1/denom                                          ! get cch ranges
    c02  = c01 + 0.5d0*(m2-m1)/m3
    m12 = m1 + m2
    if (m12 <= m3) then
      c03 = 0.5d0*m12/m3
    else
      numer = m3*m3*(3.d0*m12-m3) + m1*m1*(m1-3.d0*m3) + m2*m2*(m2-3.d0*m3)
      c03 = numer/denom
    endif

    ! 1: C<=C1; 2: C1<=C<=C2; 3: C2<=C<=C3; 4: C3<=C<=1/2 (a: m12<=m3; b: m3<m12)) 
    if (cch <= c01) then
      alpha = (denom*cch)**athird                                       ! case (1)
    else if (cch <= c02) then 
      alpha = 0.5d0*(m1 + DSQRT(m1*m1 + 8.d0*m2*m3*(cch-c01)))          ! case (2)
    else if (cch <= c03) then
      p = 2.d0*m1*m2
      q = 1.5d0*m1*m2*(m12 - 2.d0*m3*cch)
      pst = DSQRT(p)
      arc = athird*DACOS(q/(p*pst))
      csarc = DCOS(arc)
      alpha = pst*(DSQRT(3.d0*(1.d0-csarc*csarc)) - csarc) + m12        ! case (3)
    else if (m12 <= m3) then
      alpha = m3*cch + 0.5d0*m12                                       ! case (4a)
    else
      p = m12*m3 + m1*m2 - 0.25d0
      q = 1.5d0*m1*m2*m3*(0.5d0-cch)
      pst = DSQRT(p)
      arc = athird*DACOS(q/(p*pst))
      csarc = DCOS(arc)
      alpha = pst*(DSQRT(3.d0*(1.d0-csarc*csarc)) - csarc) + 0.5d0     ! case (4b)
    endif

    if (cc > 0.5d0)  alpha = 1.d0 - alpha

    FloodSZ_backward = alpha
    FloodSZ_backward = alpha + DMIN1(0.0_sp,nr(1)) + DMIN1(0.0_sp,nr(2)) + &
        &                 DMIN1(0.0_sp,nr(3))
    return
  end FUNCTION FloodSZ_Backward

  !=======================================================
  ! (2.2) Backward flooding algorithm of SZ finding centroid
  ! Calculate the cutting plane from the normal vector and
  ! volume function in UNIT CUBE.
  ! Calculate centroid of the polehedron by cutting plane
  !    m1 x1 + m2 x2 + m3 x3 = alpha
  !
  ! Adopted from Paris simulator by Zaleski
  ! Used in: MOF reconstruction
  !-------------------------------------------------------
  ! Input:
  !      nr: normal vector(nx, ny, nz)
  !      cc: vof function
  ! Output:
  !      xc0: centroid (cx,cy,cz)
  !=======================================================
  SUBROUTINE FloodSZ_BackwardC(nr,cc,xc0)

    IMPLICIT NONE
    REAL(8), INTENT(IN) :: nr(3),cc
    REAL(8), INTENT(OUT) :: xc0(3)
    REAL(8) :: ctd0(3)
    REAL(8) :: cch,ccr,alh,c01,c02,c03,np1,np2,np3
    REAL(8) :: m1,m2,m3,m12,numer,denom,p,pst,q,arc,csarc
    REAL(8) :: tmp1,tmp2,tmp3,tmp4,a2,bot1
    INTEGER :: ind(3)
    REAL(8), PARAMETER :: athird=1.d0/3.d0,eps0=1.d-50 
    INTRINSIC DMAX1,DMIN1,DSQRT,DACOS,DCOS

    np1 = DABS(nr(1))                                 ! need positive coefficients
    np2 = DABS(nr(2))
    np3 = DABS(nr(3))
    ! coefficients in ascending order
    if (np1 <= np2) then                             
      m1 = np1
      m3 = np2
      ind(1) = 1                                            
      ind(3) = 2
    else
      m1 = np2
      m3 = np1
      ind(1) = 2
      ind(3) = 1
    endif

    if (np3 < m1) then
      m2 = m1
      m1 = np3
      ind(2) = ind(1)
      ind(1) = 3
    else if (np3 >= m3) then
      m2 = m3
      m3 = np3
      ind(2) = ind(3)
      ind(3) = 3
    else
      m2 = np3
      ind(2) = 3
    endif

    denom = DMAX1(6.d0*m1*m2*m3,eps0)                           
    cch = DMIN1(cc,1.d0-cc)                              ! limit to: 0 < cch < 1/2
    ccr = DMAX1(cch,eps0)                                     ! avoid cch = m1 = 0
    c01 = m1*m1*m1/denom                                          ! get cch ranges 
    c02  = c01 + 0.5d0*(m2-m1)/m3
    m12 = m1 + m2
    if (m12 <= m3) then
      c03 = 0.5d0*m12/m3
    else
      numer = m3*m3*(3.d0*m12-m3) + m1*m1*(m1-3.d0*m3) + m2*m2*(m2-3.d0*m3)
      c03 = numer/denom
    endif

    ! 1: C<=C1; 2: C1<=C<=C2; 3: C2<=C<=C3; 4: C3<=C<=1/2 (a: m12<=m3; b: m3<m12)) 
    if (ccr <= c01) then
      alh = 0.25d0*(denom*cch)**athird 
      ctd0(1) = alh/m1 
      ctd0(2) = alh/m2
      ctd0(3) = alh/m3                                                ! case (1)
    else if (ccr <= c02) then 
      alh = 0.5d0*(m1 + DSQRT(m1*m1 + 8.d0*m2*m3*(cch - c01)))    
      tmp1 = (2.d0*alh-m1)
      tmp2 = (2.d0*alh*alh - 2.d0*alh*m1 + m1*m1)
      bot1 = DMAX1(4.d0*(3.d0*alh*alh - 3.d0*alh*m1 + m1*m1),eps0)
      tmp2 = tmp1*tmp2/bot1
      ctd0(1) = 0.5d0 - m1*tmp1/bot1 
      ctd0(2) = tmp2/m2
      ctd0(3) = tmp2/m3                                                ! case (2)
    else if (ccr <= c03) then
      p = 2.d0*m1*m2                                                   
      q = 1.5d0*m1*m2*(m12 - 2.d0*m3*cch)
      pst = DSQRT(p)
      arc = athird*DACOS(q/(p*pst))
      csarc = DCOS(arc)
      alh = pst*(DSQRT(3.d0*(1.d0-csarc*csarc)) - csarc) + m12
      a2 = alh*alh
      tmp1 = (alh-m1)*(alh-m1)*(alh-m1)
      tmp2 = (alh-m2)*(alh-m2)*(alh-m2)
      bot1 = 4.d0*(a2*alh - tmp1 - tmp2)
      ctd0(1) = (a2*a2 - tmp1*(alh+3.d0*m1) - tmp2*(alh-m2))/(m1*bot1) 
      ctd0(2) = (a2*a2 - tmp1*(alh-m1) - tmp2*(alh+3.d0*m2))/(m2*bot1) 
      ctd0(3) = (a2*a2 - tmp1*(alh-m1) - tmp2*(alh-m2))/(m3*bot1)      ! case (3)
    else if (m12 <= m3) then
      alh = m3*cch + 0.5d0*m12                                         
      bot1 = DMAX1((2.d0*alh - m12),eps0)
      ctd0(1) = 0.5d0 - m1/(6.d0*bot1)
      ctd0(2) = 0.5d0 - m2/(6.d0*bot1)                                ! case (4a) 
      ctd0(3) = ((3.d0*alh - 2.d0*m12)*bot1 + alh*m12 - m1*m2)/(6.d0*m3*bot1)
    else
      p = m12*m3 + m1*m2 - 0.25d0                                     
      q = 1.5d0*m1*m2*m3*(0.5d0-cch)
      pst = DSQRT(p)
      arc = athird*DACOS(q/(p*pst))
      csarc = DCOS(arc)
      alh = pst*(DSQRT(3.d0*(1.d0-csarc*csarc)) - csarc) + 0.5d0 
      tmp1 = m1 + m2 + m3
      tmp2 = m1*m1 + m2*m2 + m3*m3
      tmp3 = m1*m1*m1 + m2*m2*m2 + m3*m3*m3
      tmp4 = m1*m1*m1*m1 + m2*m2*m2*m2 + m3*m3*m3*m3
      a2 = alh*alh
      bot1 = 4.d0*(2.d0*a2*alh - 3.d0*a2*tmp1 + 3.d0*alh*tmp2 - tmp3)
      tmp1 = 2.d0*a2*a2 - 4.d0*a2*alh*tmp1 + 6.d0*a2*tmp2 - 4.d0*alh*tmp3 + tmp4
      ctd0(1) = (tmp1 + 4.d0*m1*(alh-m1)*(alh-m1)*(alh-m1))/(m1*bot1)
      ctd0(2) = (tmp1 + 4.d0*m2*(alh-m2)*(alh-m2)*(alh-m2))/(m2*bot1)
      ctd0(3) = (tmp1 + 4.d0*m3*(alh-m3)*(alh-m3)*(alh-m3))/(m3*bot1) ! case (4b)
    endif

    ! back to actual c value, but the centroid of a full cell is the cell center
    if (cc > 0.5d0) then
      ctd0(1) = (0.5d0 - (1.d0 - cc)*(1.d0-ctd0(1)))/cc
      ctd0(2) = (0.5d0 - (1.d0 - cc)*(1.d0-ctd0(2)))/cc
      ctd0(3) = (0.5d0 - (1.d0 - cc)*(1.d0-ctd0(3)))/cc
    endif

    ! get correct indexing
    xc0(ind(1)) = ctd0(1)                              
    xc0(ind(2)) = ctd0(2)
    xc0(ind(3)) = ctd0(3)

    ! take care of negative coefficients
    if (nr(1) < 0.d0) xc0(1) = 1.d0 - xc0(1) 
    if (nr(2) < 0.d0) xc0(2) = 1.d0 - xc0(2)
    if (nr(3) < 0.d0) xc0(3) = 1.d0 - xc0(3)

  END SUBROUTINE FloodSZ_BackwardC

  !=======================================================
  ! (2-3) Forward flooding algorithm of SZ finding f
  ! FIND THE "CUT VOLUME" V0
  ! for GIVEN
  !    r0, dr0
  !  and
  !    m1 x1 + m2 x2 + m3 x3 = alpha
  !
  ! Adopted from Paris simulator by Zaleski
  ! Used in:
  !     VOF advection
  !     MOF advection
  !-------------------------------------------------------
  ! Input:
  !      nr: normal vector(nx, ny, nz)
  !      alpha: alpha
  !      x0: the origin (x0,y0,z0)
  !      dx: dimension of the cell (dx0,dy0,dz0)
  ! Output:
  !      vof function
  !=======================================================
  Real(8) Function FloodSZ_Forward(nr,alpha,x0,dx)
    IMPLICIT NONE
    REAL(8), INTENT(IN):: nr(3),x0(3),dx(3),alpha
    REAL(8) :: al,almax,alh,np1,np2,np3,m1,m2,m3,m12,mm,denom,frac,top
    REAL(8), PARAMETER :: eps0=1.d-50
    INTRINSIC DMAX1,DMIN1,DABS

    ! move origin to x0
    al = alpha - nr(1)*x0(1) - nr(2)*x0(2) - nr(3)*x0(3)
    ! reflect the figure when negative coefficients
    al = al + DMAX1(0.d0,-nr(1)*dx(1)) + DMAX1(0.d0,-nr(2)*dx(2)) &
        + DMAX1(0.d0,-nr(3)*dx(3))
    np1 = DABS(nr(1))                                 ! need positive coefficients
    np2 = DABS(nr(2))
    np3 = DABS(nr(3))
    almax = np1*dx(1) + np2*dx(2) + np3*dx(3)
    al = DMAX1(0.d0,DMIN1(1.d0,al/almax))           !get new al within safe limits
    alh = DMIN1(al,1.d0-al)                              ! limit to: 0 < alh < 1/2


    ! normalized equation: m1*y1 + m2*y2 + m3*y3 = alh, with 0 <= m1 <= m2 <= m3
    ! the problem is then solved again in the unit cube
    np1 = np1/almax;
    np2 = np2/almax;
    np3 = np3/almax;
    m1 = DMIN1(np1*dx(1),np2*dx(2))                           ! order coefficients
    m3 = DMAX1(np1*dx(1),np2*dx(2))
    top = np3*dx(3)
    if (top < m1) then
      m2 = m1
      m1 = top
    else if (top >= m3) then
      m2 = m3
      m3 = top
    else
      m2 = top
    endif

    m12 = m1 + m2
    mm = DMIN1(m12,m3)
    denom = DMAX1(6.d0*m1*m2*m3,eps0)

    ! 1: al<=m1; 2: m1<=al<=m2; 3: m2<=al<=mm; 4: mm<=al<=1/2 (a:m12<=m3; b:m3<m12)) 
    if (alh <= m1) then
      frac = alh*alh*alh/denom                                         ! case (1)
    else if (alh <= m2) then
      frac = 0.5d0*alh*(alh-m1)/(m2*m3) +  m1*m1*m1/denom              ! case (2)
    else if (alh <= mm) then
      top = alh*alh*(3.d0*m12-alh) + m1*m1*(m1-3.d0*alh)              
      frac = (top + m2*m2*(m2-3.d0*alh))/denom                         ! case (3)
    else if (m12 <= m3) then
      frac = (alh - 0.5d0*m12)/m3                                     ! case (4a)
    else
      top = alh*alh*(3.d0-2.d0*alh) + m1*m1*(m1-3.d0*alh)             
      frac = (top + m2*m2*(m2-3.d0*alh) + m3*m3*(m3-3.d0*alh))/denom  ! case (4b)
    endif

    top = dx(1)*dx(2)*dx(3)
    if (al <= 0.5d0) then
      FloodSZ_forward = frac*top
    else
      FloodSZ_forward = (1.d0-frac)*top
    endif

  end FUNCTION FloodSZ_Forward

  !=======================================================
  ! (2-4) Forward flooding algorithm of SZ finding f and centroid
  ! FIND THE "CUT VOLUME" V0
  ! for GIVEN
  !    r0, dr0
  !  and
  !    m1 x1 + m2 x2 + m3 x3 = alpha
  !
  ! Adopted from Paris simulator by Zaleski
  ! Used in:
  !     MOF advection
  !-------------------------------------------------------
  ! Input:
  !      nr: normal vector(nx, ny, nz)
  !      alpha: alpha
  !      x0: the origin (x0,y0,z0)
  !      dx: dimension of the cell (dx0,dy0,dz0)
  ! Output:
  !      f: vof function
  !      xc0: centroid (cx, cy, cz)
  !=======================================================
  SUBROUTINE FloodSZ_forwardC(nr,alpha,x0,dx,f,xc0)

    IMPLICIT NONE
    REAL(8), INTENT(IN) :: nr(3),x0(3),dx(3),alpha
    REAL(8), INTENT(OUT) :: f
    REAL(8), INTENT(OUT) :: xc0(3)
    REAL(8) :: ctd0(3)
    REAL(8) :: al,almax,alh,alr,np1,np2,np3,m1,m2,m3,m12,mm,denom
    REAL(8) :: tmp1,tmp2,tmp3,tmp4,a2,bot1,frac
    REAL(8), PARAMETER :: eps0=1.d-50
    INTEGER :: ind(3)
    INTRINSIC DMAX1,DMIN1,DABS

    ! move origin to x0
    al = alpha - nr(1)*x0(1) - nr(2)*x0(2) - nr(3)*x0(3)
    ! reflect the figure when negative coefficients
    al = al + DMAX1(0.d0,-nr(1)*dx(1)) + DMAX1(0.d0,-nr(2)*dx(2)) &
        + DMAX1(0.d0,-nr(3)*dx(3))
    np1 = DABS(nr(1))                                 ! need positive coefficients
    np2 = DABS(nr(2))
    np3 = DABS(nr(3))
    almax = np1*dx(1) + np2*dx(2) + np3*dx(3)
    al = DMAX1(0.d0,DMIN1(1.d0,al/almax))           !get new al within safe limits
    alh = DMIN1(al,1.d0-al)                              ! limit to: 0 < alh < 1/2
    alr = DMAX1(alh,eps0)                                     ! avoid alh = m1 = 0

    ! normalized equation: m1*y1 + m2*y2 + m3*y3 = alh, with 0 <= m1 <= m2 <= m3
    ! the problem is then solved again in the unit cube
    np1 = np1*dx(1)/almax;
    np2 = np2*dx(2)/almax;
    np3 = np3*dx(3)/almax;
    ! coefficients in ascending order
    if (np1 <= np2) then
      m1 = np1
      m3 = np2
      ind(1) = 1
      ind(3) = 2
    else
      m1 = np2
      m3 = np1
      ind(1) = 2
      ind(3) = 1
    endif

    if (np3 < m1) then
      m2 = m1
      m1 = np3
      ind(2) = ind(1)
      ind(1) = 3
    else if (np3 >= m3) then
      m2 = m3
      m3 = np3
      ind(2) = ind(3)
      ind(3) = 3
    else
      m2 = np3
      ind(2) = 3
    endif

    m12 = m1 + m2
    mm = DMIN1(m12,m3)
    denom = DMAX1(6.d0*m1*m2*m3,eps0)

    ! 1: al<=m1; 2: m1<=al<=m2; 3: m2<=al<=mm; 4: mm<=al<=1/2 (a:m12<=m3; b:m3<m12)) 
    if (alr <= m1) then
      tmp1 = 0.25d0*alh
      ctd0(1) = tmp1/m1 
      ctd0(2) = tmp1/m2
      ctd0(3) = tmp1/m3                                                
      frac = alh*alh*alh/denom                                          ! case (1)
    else if (alr <= m2) then
      tmp1 = (2.d0*alh-m1)
      tmp2 = (2.d0*alh*alh - 2.d0*alh*m1 + m1*m1)
      bot1 = 4.d0*(3.d0*alr*alr - 3.d0*alh*m1 + m1*m1)
      tmp2 = tmp1*tmp2/bot1
      ctd0(1) = 0.5d0 - m1*tmp1/bot1 
      ctd0(2) = tmp2/m2
      ctd0(3) = tmp2/m3                                                
      frac = 0.5d0*alh*(alh-m1)/(m2*m3) +  m1*m1*m1/denom               ! case (2)
    else if (alr <= mm) then
      a2 = alh*alh
      tmp1 = (alh-m1)*(alh-m1)*(alh-m1)
      tmp2 = (alh-m2)*(alh-m2)*(alh-m2)
      bot1 = 4.d0*(a2*alh - tmp1 - tmp2)
      ctd0(1) = (a2*a2 - tmp1*(alh+3.d0*m1) - tmp2*(alh-m2))/(m1*bot1) 
      ctd0(2) = (a2*a2 - tmp1*(alh-m1) - tmp2*(alh+3.d0*m2))/(m2*bot1) 
      ctd0(3) = (a2*a2 - tmp1*(alh-m1) - tmp2*(alh-m2))/(m3*bot1)      
      tmp3 = alh*alh*(3.d0*m12-alh) + m1*m1*(m1-3.d0*alh)              
      frac = (tmp3 + m2*m2*(m2-3.d0*alh))/denom                         ! case (3)
    else if (m12 <= m3) then
      bot1 = DMAX1((2.d0*alh - m12),eps0)
      ctd0(1) = 0.5d0 - m1/(6.d0*bot1)
      ctd0(2) = 0.5d0 - m2/(6.d0*bot1)                                
      ctd0(3) = ((3.d0*alh - 2.d0*m12)*bot1 + alh*m12 - m1*m2)/(6.d0*m3*bot1)
      frac = (alh - 0.5d0*m12)/m3                                      ! case (4a)
    else
      tmp1 = m1 + m2 + m3
      tmp2 = m1*m1 + m2*m2 + m3*m3
      tmp3 = m1*m1*m1 + m2*m2*m2 + m3*m3*m3
      tmp4 = m1*m1*m1*m1 + m2*m2*m2*m2 + m3*m3*m3*m3
      a2 = alh*alh
      bot1 = 4.d0*(2.d0*a2*alh - 3.d0*a2*tmp1 + 3.d0*alh*tmp2 - tmp3)
      tmp1 = 2.d0*a2*a2 - 4.d0*a2*alh*tmp1 + 6.d0*a2*tmp2 - 4.d0*alh*tmp3 + tmp4
      ctd0(1) = (tmp1 + 4.d0*m1*(alh-m1)*(alh-m1)*(alh-m1))/(m1*bot1)
      ctd0(2) = (tmp1 + 4.d0*m2*(alh-m2)*(alh-m2)*(alh-m2))/(m2*bot1)
      ctd0(3) = (tmp1 + 4.d0*m3*(alh-m3)*(alh-m3)*(alh-m3))/(m3*bot1) 
      tmp1 = alh*alh*(3.d0-2.d0*alh) + m1*m1*(m1-3.d0*alh)             
      frac = (tmp1 + m2*m2*(m2-3.d0*alh) + m3*m3*(m3-3.d0*alh))/denom  ! case (4b)
    endif

    ! back to actual al value, but the centroid of a full cell is the cell center
    ! must be careful that frac = cc when al < 0.5, otherwise frac = 1 - cc
    if (al > 0.5d0) then
      ctd0(1) = (0.5d0 - frac*(1.d0-ctd0(1)))/(1.d0-frac)
      ctd0(2) = (0.5d0 - frac*(1.d0-ctd0(2)))/(1.d0-frac)
      ctd0(3) = (0.5d0 - frac*(1.d0-ctd0(3)))/(1.d0-frac)
    endif

    ! get correct indexing
    xc0(ind(1)) = ctd0(1)                              
    xc0(ind(2)) = ctd0(2)
    xc0(ind(3)) = ctd0(3)

    ! take care of negative coefficients
    if (nr(1) < 0.d0) xc0(1) = 1.d0 - xc0(1) 
    if (nr(2) < 0.d0) xc0(2) = 1.d0 - xc0(2)
    if (nr(3) < 0.d0) xc0(3) = 1.d0 - xc0(3)

    ! final position of the centroid with respect to the cell origin
    ! and by considering the actual sides of the hexahedron 
    xc0(1) = x0(1) + xc0(1)*dx(1) 
    xc0(2) = x0(2) + xc0(2)*dx(2) 
    xc0(3) = x0(3) + xc0(3)*dx(3) 

    if (al <= 0.5d0) then
      f = frac * dx(1)*dx(2)*dx(3)
    else
      f = (1.d0-frac) * dx(1)*dx(2)*dx(3)
    endif

  End Subroutine FloodSZ_forwardC

    Subroutine Normalization1(norm)
    Implicit None
    Real(8), Intent(Inout) :: norm(3)
    Real(8) :: abs_norm(3)
    Real(8) :: aa

    abs_norm = abs(norm)
    aa = abs_norm(1) + abs_norm(2) + abs_norm(3)
    If( aa .gt. 1.d-10) Then
      norm = norm / aa
    End If
  End Subroutine Normalization1

  Subroutine Norm2Angle(angle, norm)
    Implicit None
    Integer, Parameter  :: sp = 8
    Real(sp), Intent(Out)  :: angle(2)
    Real(sp), Intent(In) :: norm(3)
    Real(sp) :: epsc = 1e-10
    Real(sp) :: MOF_Pi = 3.141592653589793
    angle(2) = acos(norm(3))
    If (norm(1) > epsc) Then
      angle(1) = atan(norm(2) / norm(1) )
    Else If (norm(1) < -epsc) Then
      If (norm(2) .ge. epsc) Then
        angle(1) = atan(norm(2) / norm(1) ) + MOF_Pi
      Else
        angle(1) = atan(norm(2) / norm(1) ) - MOF_Pi
      EndIf
    Else
      If (norm(2) .gt. epsc) Then
        angle(1) = MOF_Pi / 2.0_sp
      Else If (norm(2) .lt. -epsc) Then
        angle(1) = - MOF_Pi / 2.0_sp
      Else
        angle(1) = 0.0_sp
      End If
    End If

  End Subroutine Norm2Angle

  !=======================================================
  ! (4-4) Convert the spherical angle to to Caetesian Norm
  !-------------------------------------------------------
  ! Input:  norm (normalized normal vector, nx, ny, nz)
  ! Output: angle (angle, theta, phi)
  !=======================================================
  Subroutine Angle2Norm(angle, norm)
    Implicit None
    Integer, Parameter  :: sp = 8
    Real(sp), Intent(In)  :: angle(2)
    Real(sp), Intent(Out) :: norm(3)
    norm(3) = cos(angle(2))
    norm(1) = sin(angle(2)) * cos(angle(1))
    norm(2) = sin(angle(2)) * sin(angle(1))
  End Subroutine Angle2Norm



