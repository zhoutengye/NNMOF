!=================
! Note: Some of the subroutines are in this module,
!       som eare available in the extension module (with [])
! Includes:
!     (1) Norm Calculation
!       [1-1] Norm calculation Parker and Youngs
!       (1-2) Norm calculation Central difference
!       [1-3] Norm calculation Mixed Central difference and Youngs (MYC)
!       [1-3] Norm calculation ELVIRA
!     (2) Flood algorithm (VOF Reconstruction)
!       (2-1) Flood algorithm backward finding alpha
!       (2-2) Flood algorithm backward finding centroid
!       (2-3) Flood algorithm forward finding vof
!       (2-4) Flood algorithm forward finding vof and centroid
!       [2-5] Flood algorithm backward of THINC/SW finding jump center
!       [2-6] Flood algorithm forward of THINC/SW finding vof
!       [2-7] The forward algorithm from Sussman's code
!     (3) MOF Reconstruction
!       (3-1) MOF reconstruction using quick centroid and Gauss-Newton iteration
!       [3-1-1] MOF reconstruction using quick centroid and Gauss-Newton iteration (Old version)
!       [3-2] MOF reconstruction using Lemoine's analytic gradient and Gauss-Newton iteration
!       [3-3] MOF reconstruction using Lemoine's analytic gradient and BFGS iteration
!       [3-4] MOF reconstruction using Sussman's computational geometry and Gauss-Newton iteration
!       [3-5-1] MOF reconstruction using machine learning (NNMOF)
!       [3-5-2] MOF reconstruction using machine learning-Stablized (NNMOF-Stab)
!       (3-6) Find centroid for given norm and volume function
!     (4) Misc
!       (4-1) Normalization vector 1 (nx+ny+nz=1)
!       (4-2) Normalization vector 2 (nx**2+ny**2+nz**2=1)
!       (4-3) Cartesian norm to Spherical angle
!       (4-4) Spherical angle to Cartesian norm
!       (4-5) Anvance of angle
!       (4-6-1) Initial_Guess (JCP note)
!       [4-6-2] Initial_GuessOld (Lemoine)
!       [4-6-3] Initial_GuessNN (Nueral Network)
!       [4-7] Initialize NNMOF
!       [4-8] Initialize Lemoine's algorithm
!     [5] Level set related
!       [5-1] Smooth vof function
!       [5-2] VOF to level set function near interface
!       [5-3] Fast sweeping method (FSM) for ls function for far field
!       [5-4] permulation function for FSM
!-----------------
! Author: Zhouteng Ye (yzt9zju@gmail.com)
!-----------------
! Note:
!    The origin and cell size for each function/subtourines are
!         |+++++++++++++++++++|++++++++++++++++++|+++++++++++++++++++|
!         |                   | Origin           |  cell size        |
!         |+++++++++++++++++++|++++++++++++++++++|+++++++++++++++++++|
!         | Backward flooidng | 0,0,0            | 1,1,1             |
!         | Forward flooidng  | x0(1),x0(2),x0(3)| dx(1),dx(2),dx(3) |
!         | MOFZY             ! -0.5,-0.5,-0.5   | 1,1,1             |
!         | MOFSussman        ! -0.5,-0.5,-0.5   | 1,1,1             |
!         | MOFLemoine-BFGS   ! 0,0,0            | 1,1,1             |
!         | MOFLemoine-GN     ! 0,0,0            | 1,1,1             |
!         |+++++++++++++++++++|++++++++++++++++++|+++++++++++++++++++|
!    While calling those functions, be carefully about the grid origin and size
!
!   A Mixed Central difference and Youngs (MYCorocefure pointer MOFNorm is used to determine which MOF to choose,
!   by default, it is pointed to MOFZY. To change the target function, for example
!   to MOFSussmanGaussNewton, use the following sentence in code
!        MOFNorm => MOFSussmanGaussNewton
!=================
Module ModVOFFunc
  Use ModGlobal, only : sp
  Implicit None
  public

  ! MOF iteration parameters
  Integer :: MOFITERMAX = 10
  Real(sp) :: mof_tol = 1.0d-8
  Real(sp) :: mof_tol_dangle = 1.0d-8
  Real(sp) :: det_lim = 1.0d-30
  Real(sp),parameter :: MOF_Pi = 3.1415926535897932d0
  Real(sp) :: epsc = 1.0d-15
  Integer :: mof_niter(2)
  Real(sp) :: delta_theta =1d-5
  Real(sp) :: delta_theta_max = 3.0_sp * MOF_Pi / 180.0_sp  ! 10 degrees

  PROCEDURE(InterfaceMOF), POINTER :: MOFNorm => MOFZY
  PROCEDURE(InterfaceVOF), POINTER :: VOFNorm => NormParkerYoungs

  Interface
    Subroutine InterfaceVOF(f,norm)
      Implicit None
      Real(8), Intent(In) :: f(3,3,3)
      Real(8), Intent(Out) :: norm(3)
    End Subroutine InterfaceVOF
  End Interface

  Interface
    Subroutine InterfaceMOF(f,c,norm, init_norm)
      Implicit None
      Real(8), Intent(In) :: f
      Real(8), Intent(In) :: c(3)
      Real(8), Intent(Out) :: norm(3)
      Real(8), Intent(In), optional :: init_norm(3)
    End Subroutine InterfaceMOF
  End Interface

Contains
  !=======================================================
  ! (1-1) Normal vector Youngs (1992)
  !-------------------------------------------------------
  ! Input: f (vof function, 3*3*3 stencil)
  ! Output: norm (normal vector)
  !=======================================================
  Subroutine NormParkerYoungs(f, norm)
    Implicit None
    Real(sp), Intent(In)   :: F(3,3,3)
    Real(sp), Intent(Out)  :: norm(3)
    Real(sp) :: mm1, mm2
    norm(3) = 0.0_sp
    mm1 = f(1,1,1)+f(1,1,3)+f(1,3,1) &
        +f(1,3,3)+2.0_sp*(f(1,1,2)+f(1,3,2)&
        +f(1,2,1)+f(1,2,3))+4.0_sp*f(1,2,2)
    mm2 = f(3,1,1)+f(3,1,3)+f(3,3,1) &
        +f(3,3,3)+2.0_sp*(f(3,1,2)+f(3,3,2) &
        +f(3,2,1)+f(3,2,3))+4.0_sp*f(3,2,2)
    norm(1) = mm1 - mm2

    mm1 = f(1,1,1)+f(1,1,3)+f(3,1,1) &
        +f(3,1,3)+2.0_sp*(f(1,1,2)+f(3,1,2) &
        +f(2,1,1)+f(2,1,3))+4.0_sp*f(2,1,2)
    mm2 = f(1,3,1)+f(1,3,3)+f(3,3,1) &
        +f(3,3,3)+2.0_sp*(f(1,3,2)+f(3,3,2) &
        +f(2,3,1)+f(2,3,3))+4.0_sp*f(2,3,2)
    norm(2) = mm1 - mm2

    mm1 = f(1,1,1)+f(1,3,1)+f(3,1,1) &
        +f(3,3,1)+2.0_sp*(f(1,2,1)+f(3,2,1) &
        +f(2,1,1)+f(2,3,1))+4.0_sp*f(2,2,1)
    mm2 = f(1,1,3)+f(1,3,3)+f(3,1,3) &
        +f(3,3,3)+2.0_sp*(f(1,2,3)+f(3,2,3) &
        +f(2,1,3)+f(2,3,3))+4.0_sp*f(2,2,3)
    norm(3) = mm1 - mm2

  End Subroutine NormParkerYoungs

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
  Real(sp) Function FloodSZ_Backward(nr,cc)
    Implicit None
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
  End Function FloodSZ_Backward

  !=======================================================
  ! (2-2) Backward flooding algorithm of SZ finding centroid
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
  Real(sp) Function FloodSZ_Forward(nr,alpha,x0,dx)
    IMPLICIT NONE
    REAL(sp), INTENT(IN):: nr(3),x0(3),dx(3),alpha
    REAL(sp) :: al,almax,alh,np1,np2,np3,m1,m2,m3,m12,mm,denom,frac,top
    REAL(sp), PARAMETER :: eps0=1.d-50
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

  !=======================================================
  ! (3-1) MOF reconstruction using Gauss-Newton iteration
  !     For given vof and centroid, find the best norm
  ! Key steps:
  !   - shift the angle positive/negative in each direction (4 values)
  !   - Find the centroid corresponding with shifted angle
  !   - Solve Jacobi Matrix
  !          - if det != 0, determine the shift angle for next step
  !          0 if det = 0, choose the best shift angle
  !-------------------------------------------------------
  ! Input:
  !      f: vof function
  !      c: centroid (cx, cy, cz)
  ! Optional Input:
  !      init_norm: Initial guess of the normal vector
  ! Output:
  !      norm: vof function
  !=======================================================
  Subroutine MOFZY(f, c, norm, Init_Norm)
    ! Use ModOptimizer
    Implicit None
    Real(sp), Intent(In)  :: f
    Real(sp), Intent(In)  :: c(3)
    Real(sp), Intent(Out) :: norm(3)
    Real(sp), Optional, Intent(In) :: Init_Norm(3)
    Real(sp)  :: vof
    Real(sp)  :: c_mof(3) ,c_mof_sym(3)
    Real(sp), Dimension(3)   :: norm_2
    Real(sp), Dimension(2)   :: delangle, angle_init,  new_angle
    Real(sp), Dimension(2)   :: angle_base, angle_plus, angle_minus
    Real(sp), Dimension(3)   :: cenopt, cenp, cenm, cen_init
    Real(sp), Dimension(3)   :: cenopt_sym
    Real(sp), Dimension(3,2) :: cen_plus, cen_minus
    Real(sp), Dimension(3,2) :: Jacobian
    Real(sp), Dimension(2,2) :: Hessian, HessianT
    Real(sp), Dimension(2)   :: gradient
    Real(sp), Dimension(3)   :: c_diff, c_diff_sym
    Real(sp)   :: det
    Real(sp), Dimension(2,MOFITERMAX+1) :: angle_array
    Real(sp), Dimension(MOFITERMAX+1)   :: err_array

    Real(sp) :: err
    Real(sp) :: mof_tol2
    Integer :: singular_flag = 0
    Integer :: i_angle, j_angle, dir, iter, i, j

    c_mof = c - 0.5_sp

    if (f .ge. 0.5_sp) then
      vof = 1.0_sp - f
      c_mof =  - c_mof * f / vof
    else
      vof = f
      c_mof = c_mof
    endif
    c_mof_sym = - c_mof * f / (1.0_sp - f)

    ! Initialize angle
    If(present(Init_Norm))then
      norm_2 = Init_Norm
      Call Normalization2(Norm_2)
      Call Norm2Angle(angle_init,norm_2)
      Call FindCentroid(angle_init,vof,cen_init)
      c_diff = cen_init - c_mof
      err = dot_product(c_diff,c_diff)
    Else
      Call Initial_Guess(c_mof, vof, angle_init, err)
    EndIf

    ! Initialize other data
    do dir=1,2
      angle_array(dir,1)= angle_init(dir)
    enddo
    do dir=1,3
      err_array(1) = err
    enddo

    iter = 0
    err = err_array(1)
    mof_niter = 0

    mof_tol2 = mof_tol * (1.0_sp-f)**2
    mof_tol2 = mof_tol

    Do While ((iter.lt.MOFITERMAX).and. &
        (err.gt.mof_tol2))

      Do i_angle=1,2
        angle_base(i_angle) = angle_array(i_angle, iter+1)
      End Do

      ! theta phi shift angle
      ! to calculate the numerical gradient
      Do i_angle=1,2
        Do j_angle=1,2
          angle_plus(j_angle)=angle_base(j_angle)
          angle_minus(j_angle)=angle_base(j_angle)
        End Do
        angle_plus(i_angle)=angle_plus(i_angle)+delta_theta
        angle_minus(i_angle)=angle_minus(i_angle)-delta_theta
        Call FindCentroid(angle_plus,vof,cenp)
        do dir=1,3
          cen_plus(dir,i_angle)=cenp(dir)
        end do
        Call FindCentroid(angle_minus,vof,cenm)
        do dir=1,3
          cen_minus(dir,i_angle)=cenm(dir)
        enddo

        ! jacobian matrix (numerical partial gradient):
        do dir=1,3
          Jacobian(dir,i_angle)=(cenp(dir)-cenm(dir))/(2.0_sp*delta_theta)
        enddo

      End Do ! end theta phi shift angle

      Call FindCentroid(angle_base,vof,cenopt)
      cenopt_sym = - cenopt * f / (1.0_sp - f)
      c_diff = cenopt - c_mof
      c_diff_sym = cenopt_sym - c_mof_sym

      Do i=1,2
        Do j=1,2
          Hessian(i,j) = dot_product(Jacobian(:,i), Jacobian(:,j))
        EndDo
        gradient(i) = dot_product(c_diff, Jacobian(:,i))
      EndDo
      ! err = dot_product(c_diff, c_diff) &
      !     + dot_product(c_diff_sym, c_diff_sym)
      err = dot_product(c_diff, c_diff)

      det = Hessian(1,1)*Hessian(2,2) - Hessian(1,2)*Hessian(2,1)
      ! Calculate the inverse of the matrix
      HessianT(1,1) = +Hessian(2,2) / det
      HessianT(2,1) = -Hessian(2,1) / det
      HessianT(1,2) = -Hessian(1,2) / det
      HessianT(2,2) = +Hessian(1,1) / det

      ! Call matinv2(Hessian, HessianT, det)
      If (det .lt. det_lim) Then
        Singular_flag = 1
      End If

      ! Delta angle
      Do i=1,2
        delangle(i) = - dot_product(HessianT(:,i), gradient)
      End Do

      If (singular_flag.eq.1) then
        delangle(1:2)=0.0_sp
        err = 0.0_sp
      Else
        Do i_angle = 1,2
          If (delangle(i_angle).gt.delta_theta_max) then
            delangle(i_angle)=delta_theta_max
          Else If (delangle(i_angle).lt.-delta_theta_max) then
            delangle(i_angle)=-delta_theta_max
          End If
        End Do
      End If


      call advance_angle(angle_base,delangle)
      do i_angle=1,2
        ! call advance_angle(angle_base(i_angle),delangle(i_angle))
        angle_array(i_angle,iter+2)=angle_base(i_angle)
      enddo

      err_array(iter+2)=err
      iter=iter+1
      if ( dot_product(delangle,delangle) .lt. mof_tol_dangle) exit

    End Do

    mof_niter(1) = iter + 1
#if defined(DEBUG)
    Block
      Integer :: ii
      Do ii = 1, mof_niter(1)
        print *, '=====step',ii-1,'=========='
        ! print *, cen_array(:,ii)
        ! print *, angle_array(:,ii)
        ! print *, err_array(ii)
      End Do
    End Block
#endif

    do dir=1,2
      new_angle(dir)=angle_array(dir,iter+1)
    enddo

    ! Convert angle to norm
    call Angle2Norm(new_angle,norm)
    Call Normalization1(norm)

    if (f .ge. 0.5_sp) norm = -norm
#if defined(DEBUG)
    print *, '=====MOF Reconc norm=========='
    print *, 'Normal vector:', norm
    print *, '=============================='
#endif
  End Subroutine MOFZY

  !=======================================================
  ! (3-7) Find the centroid for MOF iteration
  !  shift the angle to Cartesian norm, then calculate
  !  the norm, finally shift the origin by -0.5 in each
  !  direction
  !-------------------------------------------------------
  ! Input:
  !      angle: angle (theta, phi)
  !      vof: volume fraction
  ! Output:
  !      cen_mof: centroid of mof
  !=======================================================
  Subroutine FindCentroid(angle, vof, cen_mof)
    Implicit None
    Integer, Parameter  :: sp = 8
    Real(sp), Intent(In) :: angle(2)
    Real(sp), Intent(In) :: vof
    Real(sp), Intent(Out) :: cen_mof(3)
    Real(sp) :: cen_sz(3)
    Real(sp) :: norm(3)

    Call Angle2Norm(angle, norm)
    Call Normalization1(norm)
    Call FloodSZ_BackwardC(norm,vof,cen_sz)
    cen_mof = cen_sz - 0.5_sp
  End Subroutine FindCentroid

  !=======================================================
  ! (4-1) Normalize the normal vector that satisfies
  !    sigma norm(i) = 1
  ! and
  !    all norm(i) > 0
  !-------------------------------------------------------
  ! Input: norm (normal vector)
  ! Output: norm (normalized normal vector)
  !=======================================================
  Subroutine Normalization1(norm)
    Implicit None
    Real(sp), Intent(Inout) :: norm(3)
    Real(sp) :: abs_norm(3)
    Real(sp) :: aa

    abs_norm = abs(norm)
    aa = abs_norm(1) + abs_norm(2) + abs_norm(3)
    norm = norm / (aa + 1d-30)
  End Subroutine Normalization1

  !=======================================================
  ! (4-2) Normalize the normal vector that is a unit vector
  !     sigma norm(i)^2 = 1
  !-------------------------------------------------------
  ! Input: norm (normal vector)
  ! Output: norm (normalized normal vector)
  !=======================================================
  Subroutine Normalization2(norm)
    Implicit None
    Real(sp), Intent(Inout) :: norm(3)
    Real(sp) :: aa
    aa = Sqrt( norm(1) * norm(1) + &
        &       norm(2) * norm(2) + &
        &       norm(3) * norm(3) )
    norm(1) = norm(1) / ( aa + 1d-30)
    norm(2) = norm(2) / ( aa + 1d-30)
    norm(3) = norm(3) / ( aa + 1d-30)
  End Subroutine Normalization2

  !=======================================================
  ! (4-3) Convert the Cartesian norm to spherical angle
  !-------------------------------------------------------
  ! Input:  angle (angle, theta, phi)
  ! Output: norm (normalized normal vector, nx, ny, nz)
  !=======================================================

  Subroutine Norm2Angle(angle, norm)
    Implicit None
    Integer, Parameter  :: sp = 8
    Real(sp), Intent(Out)  :: angle(2)
    Real(sp), Intent(In) :: norm(3)
    Real(sp) :: norm0(3)
    norm0 = norm
    Call Normalization2(norm0)
    angle(2) = dacos(norm0(3))
    angle(1) = datan2(norm0(2), norm0(1))
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
    norm(3) = dcos(angle(2))
    norm(1) = dsin(angle(2)) * dcos(angle(1))
    norm(2) = dsin(angle(2)) * dsin(angle(1))
  End Subroutine Angle2Norm

  !=======================================================
  ! (4-5) Advance angle with delta angle
  !       limit in between - pi and pi
  !-------------------------------------------------------
  ! Input:  norm (normalized normal vector, nx, ny, nz)
  ! Output: angle (angle, theta, phi)
  !=======================================================
  Subroutine advance_angle(angle,delangle)
    Implicit None
    Integer, Parameter  :: sp = 8
    REAL(sp) angle(2),delangle(2)
    ! REAL(sp) angle,delangle
    angle=angle+delangle
    if (angle(1).lt.-MOF_Pi) then
      angle(1)=angle(1)+2.0_sp*MOF_Pi
    endif
    if (angle(1).gt.MOF_Pi) then
      angle(1)=angle(1)-2.0_sp*MOF_Pi
    end if
    if (angle(2).lt.0.0_sp) then
      angle(2)=-MOF_Pi
    endif
    if (angle(2).gt.MOF_Pi) then
      angle(2)=2.0_sp * MOF_Pi - angle(2)
    end if
  End Subroutine advance_angle

  !=======================================================
  ! (4-6-1) Find the initial guess
  !-------------------------------------------------------
  ! Input:  
  !      c_norm: centroid ranging from [0,1] (cx,cy,cz)
  !      vof: volume fraction ranging from [0,1]
  ! Output: 
  !      angle: initial guess of angle (theta, phi)
  !      error: initial error angle (angle, theta, phi)
  !=======================================================
  Subroutine Initial_Guess(C_mof, vof, angle, error)
    Implicit None
    Real(sp), Intent(In) :: C_mof(3)
    Real(sp), Intent(In) :: vof
    Real(sp), Intent(Out) :: angle(2)
    Real(sp), Intent(Out) :: error
    Real(sp) :: norm_1(3), norm_2(3)
    Real(sp) :: cen1(3), cen2(3)
    Real(sp) :: angle1(2), angle2(2)
    Real(sp) :: err1, err2
    Real(sp) :: permutation(3)
    Real(sp) :: relative_c(3)
    Integer :: dir

    Norm_1(1) = - c_mof(1)
    Norm_1(2) = - c_mof(2)
    Norm_1(3) = - c_mof(3)

    relative_c = 1.0_sp
    Do dir = 1, 3
      If ( c_mof(dir) .gt. 0.0_sp) Then
        permutation(dir) = -1.0_sp
        relative_c(dir) = 0.5-c_mof(dir)
      Else
        permutation(dir) = 1.0_sp
        relative_c(dir) = 0.5+c_mof(dir)
      End If
      norm_2(dir) = 1.0_sp / ( relative_c(dir) + 1d-20)
    End Do
    norm_2(:) = norm_2(:) * permutation(:)

    call normalization2(norm_1)
    call norm2angle(angle1,norm_1)
    call findcentroid(angle1,vof,cen1)
    call normalization2(norm_2)
    call norm2angle(angle2,norm_2)
    call findcentroid(angle2,vof,cen2)

    err1 = norm2(cen1-c_mof)
    err2 = norm2(cen2-c_mof)
    if (err1 < err2) then
      angle = angle1
      error = err1
    else
      angle = angle2
      error = err2
    endif

  End Subroutine Initial_Guess

End Module ModVOFFunc

