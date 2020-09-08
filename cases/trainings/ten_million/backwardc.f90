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


  Subroutine Normalization2(norm)
    Implicit None
    Real(8), Intent(Inout) :: norm(3)
    Real(8) :: aa
    aa = Sqrt( norm(1) * norm(1) + &
        &       norm(2) * norm(2) + &
        &       norm(3) * norm(3) )
    If( aa .gt. 1.d-10) Then
      norm(1) = norm(1) / aa
      norm(2) = norm(2) / aa
      norm(3) = norm(3) / aa
    End If
  End Subroutine Normalization2
