#define MLMethod Decision_Tree
program test
  Implicit None

  call err_map_angle
  call err_map_xyz
  call err_map_angle_training

end program test

Subroutine err_map_xyz
  Implicit None
  Integer,parameter :: len = 4000000
  Integer,parameter :: len1 = 100
  Integer, parameter :: SP = 8
  Real(SP) :: pi
  Real(SP) :: f(len)
  Real(SP) :: f1(len/2)
  Real(SP) :: f2(len/2)
  Real(SP) :: p(len)
  Real(SP) :: t(len)
  Real(SP) :: nx(len)
  Real(SP) :: ny(len)
  Real(SP) :: nz(len)
  Real(SP) :: x(len)
  Real(SP) :: y(len)
  Real(SP) :: z(len)
  Real(SP) :: xc(len1)
  Real(SP) :: yc(len1)
  Real(SP) :: zc(len1)
  Real(SP) :: dp(len)
  Real(SP) :: dt(len)
  Real(SP) :: pt(2), pt2(2), norm(3), c(3)
  Real(SP) :: xi, err
  Real(SP) :: ff(len1,len1,len1)
  Integer  :: nf(len1,len1,len1)
  Real(SP) :: pp(len1,len1,len1)
  Integer  :: np(len1,len1,len1)
  Real(SP) :: tt(len1,len1,len1)
  Integer  :: nt(len1,len1,len1)

  Integer :: i,j,k,l

  pi = 4.0_sp * atan (1.0_sp)

  Call random_number(nx)
  Call random_number(ny)
  Call random_number(nz)
  Call random_number(p)
  Call random_number(t)
  Call random_number(f1)
  Call random_number(f2)
  nx = (nx - 0.5_sp) * 2.0_sp
  ny = (ny - 0.5_sp) * 2.0_sp
  nz = (nz - 0.5_sp) * 2.0_sp
  f1 = f1 / 20.0_sp
  p =  (p -0.5 ) * 2.0 * pi 
  t =  t * pi 
  f(1:len/2) = f1
  f(len/2+1:len) = f2

  Do i = 1, len
    ! norm = (/nx(i), ny(i), nz(i)/)
    ! Call Normalization2(norm)
    ! Call norm2angle(pt, norm)
    pt = (/p(i), t(i)/)
    Call Angle2Norm(pt, norm)
    Call Normalization1(norm)
    Call FloodSZ_BackwardC(norm, f(i), c)
    c = c - 0.5_sp
    Call Initial_GuessOld(c, f(i), pt2, err)
    ! Call MOFZY(f(k),c,norm)
    dp(i) = pt2(1) - pt(1)
    dt(i) = pt2(2) - pt(2)
    x(i) = c(1)
    y(i) = c(2)
    z(i) = c(3)
  End Do

  nf = 0
  ff = 0.0_sp
  np = 0
  pp = 0.0_sp
  nt = 0
  tt = 0.0_sp
  Do l = 1, len
    i = floor( (x(l)+0.5) / 0.01_sp) + 1
    j = floor( (y(l)+0.5) / 0.01_sp) + 1
    k = floor( (z(l)+0.5) / 0.01_sp) + 1
    ff(i,j,k) = ff(i,j,k) + f(l)
    nf(i,j,k) = nf(i,j,k) + 1
    pp(i,j,k) = pp(i,j,k) + dp(l)
    np(i,j,k) = np(i,j,k) + 1
    tt(i,j,k) = tt(i,j,k) + dt(l)
    nt(i,j,k) = nt(i,j,k) + 1
  End Do

  open(10,file='pts3.dat')
  ! Do k = 1, len1
  !   Do j = 1, len1
  !     Do i = 1, len1
  !       Write(10,'(6F10.5)') xc(i), yc(j), zc(k), pp(i,j,k), tt(i,j,k), ff(i,j,k)
  !     End Do
  !   End Do
  ! End Do
  Do i = 1, len
    Write(10,'(6F10.5)') x(i), y(i), z(i), dp(i), dt(i), f(i)
  End Do

end Subroutine err_map_xyz

Subroutine err_map_angle
  Implicit None
  Integer,parameter :: len = 100
  Integer,parameter :: len2 = 100
  Integer, parameter :: SP = 8
  Real(SP) :: pi
  Real(SP) :: p(len)
  Real(SP) :: t(len)
  Real(SP) :: f(len2)
  Real(SP) :: x(len,len,len2)
  Real(SP) :: y(len,len,len2)
  Real(SP) :: z(len,len,len2)
  Real(SP) :: dp(len,len,len2)
  Real(SP) :: dt(len,len,len2)
  Real(SP) :: dp2(len,len,len2)
  Real(SP) :: dt2(len,len,len2)
  Real(SP) :: phi(len,len,len2)
  Real(SP) :: theta(len,len,len2)
  Real(SP) :: pt(2), pt2(2), norm(3), c(3)
  Real(SP) :: xi, err

  Integer :: i, j, k

  pi = 4.0_sp * atan (1.0_sp)

  Do i = 1, len
    xi = dble(i) / len - 0.5_sp / len
    p(i) = xi * 2.0_sp * pi - pi
    t(i) = xi * pi
  End Do
  Do i = 1, len2
    xi = dble(i) / len - 0.5_sp / len
    f(i) = xi
  End Do

  Do k = 1, len2
    Do j = 1, len
      Do i = 1, len
        pt = (/p(i),t(j)/)
        Call Angle2Norm(pt,norm)
        Call Normalization1(norm)
        Call FloodSZ_BackwardC(norm, f(k), c)
        c = c - 0.5_sp
        Call Initial_GuessOld(c, f(k), pt2, err)
        ! Call MOFZY(f(k),c,norm)
        dp(i,j,k) = pt2(1) - pt(1)
        dt(i,j,k) = pt2(2) - pt(2)
        x(i,j,k) = c(1)
        y(i,j,k) = c(2)
        z(i,j,k) = c(3)
      End Do
    End Do
  End Do

  open(10,file='pts2.dat')
  Do k = 1, len2
    Do j = 1, len
      Do i = 1, len
        Write(10,'(8F10.5)') p(i), t(j), f(k), dp(i,j,k), dt(i,j,k), x(i,j,k), y(i,j,k), z(i,j,k)
      End Do
    End Do
  End Do

end Subroutine err_map_angle

Subroutine err_map_angle_training
  Use DecisionTree
  Use RandomForest
  Use NueralNetwork
  Implicit None
  Integer,parameter :: len = 100
  Integer,parameter :: len2 = 100
  Integer, parameter :: SP = 8
  Real(SP) :: pi
  Real(SP) :: p(len)
  Real(SP) :: t(len)
  Real(SP) :: f(len2)
  Real(SP) :: x(len,len,len2)
  Real(SP) :: y(len,len,len2)
  Real(SP) :: z(len,len,len2)
  Real(SP) :: dp(len,len,len2)
  Real(SP) :: dt(len,len,len2)
  Real(SP) :: dp2(len,len,len2)
  Real(SP) :: dt2(len,len,len2)
  Real(SP) :: phi(len,len,len2)
  Real(SP) :: theta(len,len,len2)
  Real(SP) :: pt(2), pt2(2), norm(3), c(3)
  Real(SP) :: xi, err
  Type(MLMethod) :: ml

  Integer :: i, j, k

  Call ml%Initialization

  pi = 4.0_sp * atan (1.0_sp)

  Do i = 1, len
    xi = dble(i) / len - 0.5_sp / len
    p(i) = xi * 2.0_sp * pi - pi
    t(i) = xi * pi
  End Do
  Do i = 1, len2
    xi = dble(i) / len - 0.5_sp / len
    f(i) = xi
  End Do

  Do k = 1, len2
    Do j = 1, len
      Do i = 1, len
        pt = (/p(i),t(j)/)
        Call Angle2Norm(pt,norm)
        Call Normalization1(norm)
        Call FloodSZ_BackwardC(norm, f(k), c)
        c = c - 0.5_sp
        Call Initial_GuessOld(c, f(k), pt2, err)
        Call Initial_GuessNN(ml, c, f(k), pt, err)
        dp(i,j,k) = pt2(1) - pt(1)
        dt(i,j,k) = pt2(2) - pt(2)
        x(i,j,k) = c(1)
        y(i,j,k) = c(2)
        z(i,j,k) = c(3)
      End Do
    End Do
  End Do

  open(10,file='pts1.dat')
  Do k = 1, len2
    Do j = 1, len
      Do i = 1, len
        Write(10,'(8F10.5)') p(i), t(j), f(k), dp(i,j,k), dt(i,j,k), x(i,j,k), y(i,j,k), z(i,j,k)
      End Do
    End Do
  End Do

end Subroutine err_map_angle_training


Subroutine FloodSZ_BackwardC(nr,cc,xc0)

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

End Subroutine FloodSZ_BackwardC

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
  angle(2) = dacos(norm(3))
  angle(1) = datan2(norm(2), norm(1))
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

Subroutine Initial_GuessOld(C_mof, vof, angle, error)
  Implicit None
  Integer, Parameter  :: sp = 8
  Real(sp), Intent(In) :: C_mof(3)
  Real(sp), Intent(In) :: vof
  Real(sp), Intent(Out) :: angle(2)
  Real(sp), Intent(Out) :: error
  Real(sp) :: norm(3)
  Real(sp) :: cen(3)

  Norm(1) = - c_mof(1)
  Norm(2) = - c_mof(2)
  Norm(3) = - c_mof(3)
  call normalization2(norm)
  call norm2angle(angle,norm)
  error = norm2(cen-c_mof)

End Subroutine Initial_GuessOld

  Subroutine Initial_GuessNN(ml, C_mof, vof, angle, error)
    Use DecisionTree
    Use RandomForest
    Use NueralNetwork
    Implicit None
    Type(MLMethod) :: ml
    Integer, Parameter :: sp = 8
    Real(sp), Intent(In) :: C_mof(3)
    Real(sp), Intent(In) :: vof
    Real(sp), Intent(Out) :: angle(2)
    Real(sp), Intent(Out) :: error
    Real(sp) :: norm_1(3), norm_2(3)
    Real(sp) :: cen1(3), cen2(3), cen3(3)
    Real(sp) :: angle1(2), angle2(2), angle3(2)
    Real(sp) :: err1, err2, err3
    Real(sp) :: permutation(3)
    Real(sp) :: relative_c(3)
    Real(sp) :: inputs(3)
    Integer :: dir

    Norm_1(1) = - c_mof(1)
    Norm_1(2) = - c_mof(2)
    Norm_1(3) = - c_mof(3)

    call normalization2(norm_1)
    call norm2angle(angle1,norm_1)

    inputs(1) = angle1(1)
    inputs(2) = angle1(2)
    inputs(3) = vof
    angle = angle1 + ml%predict(inputs)

    error = 1.0_sp

  End Subroutine Initial_GuessNN


  Subroutine MOFZY(f, c, norm)
    ! Use ModOptimizer
    Use DecisionTree
    Implicit None
    Integer,  parameter :: sp=8
  Integer,parameter :: MOFITERMAX = 10
  Real(sp) :: mof_tol = 1.0d-8
  Real(sp) :: mof_tol_dangle = 1.0d-8
  Real(sp) :: det_lim = 1.0d-30
  Real(sp),parameter :: MOF_Pi = 3.1415926535897932d0
  Real(sp) :: epsc = 1.0d-15
  Integer :: mof_niter(2)
  Real(sp) :: delta_theta =1d-5
  Real(sp) :: delta_theta_max = 3.0_sp * MOF_Pi / 180.0_sp  ! 10 degrees

    Real(sp), Intent(In)  :: f
    Real(sp), Intent(In)  :: c(3)
    Real(sp), Intent(Out) :: norm(3)
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
      Call Initial_Guess(c_mof, vof, angle_init, err)

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


    do dir=1,2
      new_angle(dir)=angle_array(dir,iter+1)
    enddo

    ! Convert angle to norm
    call Angle2Norm(new_angle,norm)
    Call Normalization1(norm)

    if (f .ge. 0.5_sp) norm = -norm
  End Subroutine MOFZY


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

  Subroutine Initial_Guess(C_mof, vof, angle, error)
    Implicit None
    Integer, Parameter  :: sp = 8
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
    ! angle = angle2

  End Subroutine Initial_Guess

  Subroutine advance_angle(angle,delangle)
    Implicit None
    Integer, Parameter  :: sp = 8
    REAL(sp) angle(2),delangle(2)
    REAL(sp) mof_pi
    mof_pi = 4.0_sp * atan (1.0_sp)
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

