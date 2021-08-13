!=================
! Extension module for ModVOFFunc in vof_func.f90
! Functions/subroutines are listed tems listed in vof.func.f90 with []
!=================
Module ModVOFFuncExt
  Use ModGlobal, only: sp
  Use ModVOFFunc
  Use NueralNetwork
  Use DecisionTree
  Use RandomForest
  use mod_cg3_polyhedron
  use mod_cg3_complete_polyhedron_structure
  Implicit None

  type(t_polyhedron) :: LemoinePoly
  ! type(Neural_Network) :: ML
  ! type(Decision_Tree) :: ML
  type(Random_Forest) :: ML
Contains

  !=======================================================
  ! [1-2] Normal vector Central Scheme
  ! * returns normal normalized so that |mx|+|my|+|mz| = 1*
  ! Adopted from Paris Simulator by Zakeski
  !-------------------------------------------------------
  ! Input: c (vof function, 3*3*3 stencil)
  ! Output: norm (normal vector)
  !=======================================================
  Subroutine NormCS(c,norm)
    !***
    Implicit None
    Real(8), Intent(In)  :: c(3,3,3)
    Real(8), Intent(Out) :: norm(3)
    Real(8) m1,m2,m(0:3,0:2),t0,t1,t2
    Integer cn

    ! write the plane as: sgn(mx) X =  my Y +  mz Z + alpha 
    !                           m00 X = m01 Y + m02 Z + alpha 

    m1 = c(1,2,1) + c(1,2,3) + c(1,1,2) + c(1,3,2) + &
        c(1,2,2)
    m2 = c(3,2,1) + c(3,2,3) + c(3,1,2) + c(3,3,2) + &
        c(3,2,2)

    if(m1>m2) then
      m(0,0) = 1.
    else
      m(0,0) = -1.
    end if

    m1 = c(1,1,2)+ c(3,1,2)+ c(2,1,2)
    m2 = c(1,3,2)+ c(3,3,2)+ c(2,3,2)
    m(0,1) = 0.5*(m1-m2)

    m1 = c(1,2,1)+ c(3,2,1)+ c(2,2,1)
    m2 = c(1,2,3)+ c(3,2,3)+ c(2,2,3)
    m(0,2) = 0.5*(m1-m2)

    ! write the plane as: sgn(my) Y =  mx X +  mz Z + alpha, 
    !                          m11 Y = m10 X + m12 Z + alpha.
    m1 = c(1,1,2) + c(1,3,2) + c(1,2,2)
    m2 = c(3,1,2) + c(3,3,2) + c(3,2,2)
    m(1,0) = 0.5*(m1-m2)

    m1 = c(2,1,1) + c(2,1,3) + c(3,1,2) + c(1,1,2) +&
        c(2,1,2)
    m2 = c(2,3,1) + c(2,3,3) + c(3,3,2) + c(1,3,2) +&
        c(2,3,2)

    if(m1>m2) then
      m(1,1) = 1.
    else
      m(1,1) = -1.
    end if

    m1 = c(2,1,1)+ c(2,2,1)+ c(2,3,1)
    m2 = c(2,1,3)+ c(2,2,3)+ c(2,3,3)
    m(1,2) = 0.5*(m1-m2)

    m1 = c(1,2,1)+ c(1,2,3)+ c(1,2,2)
    m2 = c(3,2,1)+ c(3,2,3)+ c(3,2,2)
    m(2,0) = 0.5*(m1-m2)

    m1 = c(2,1,1)+ c(2,1,3)+ c(2,1,2)
    m2 = c(2,3,1)+ c(2,3,3)+ c(2,3,2)
    m(2,1) = 0.5*(m1-m2)

    m1 = c(1,2,1) + c(3,2,1) + c(2,1,1) + c(2,3,1) +&
        c(2,2,1)
    m2 = c(1,2,3) + c(3,2,3) + c(2,1,3) + c(2,3,3) +&
        c(2,2,3)

    if(m1>m2) then
      m(2,2) = 1.
    else
      m(2,2) = -1.
    end if

    ! normalize each set (mx,my,mz): |mx|+|my|+|mz| = 1

    t0 = DABS(m(0,0)) + DABS(m(0,1)) + DABS(m(0,2))
    m(0,0) = m(0,0)/t0
    m(0,1) = m(0,1)/t0
    m(0,2) = m(0,2)/t0

    t0 = DABS(m(1,0)) + DABS(m(1,1)) + DABS(m(1,2))
    m(1,0) = m(1,0)/t0
    m(1,1) = m(1,1)/t0
    m(1,2) = m(1,2)/t0

    t0 = DABS(m(2,0)) + DABS(m(2,1)) + DABS(m(2,2))
    m(2,0) = m(2,0)/t0
    m(2,1) = m(2,1)/t0
    m(2,2) = m(2,2)/t0

    ! choose among the three central schemes */ 
    t0 = DABS(m(0,0))
    t1 = DABS(m(1,1))
    t2 = DABS(m(2,2))

    cn = 0
    if (t1 > t0) then
      t0 = t1
      cn = 1
    endif

    if (t2 > t0) cn = 2

    ! components of the normal vector */
    norm(1) = m(cn,0)
    norm(2) = m(cn,1)
    norm(3) = m(cn,2)

    return
  End Subroutine NormCS

  !=======================================================
  ! [1-3] Mixed Youngs and Central difference
  ! * returns normal normalized so that |mx|+|my|+|mz| = 1*
  ! Adopted from Paris Simulator by Zakeski
  !-------------------------------------------------------
  ! Input: c (vof function, 3*3*3 stencil)
  ! Output: norm (normal vector)
  !=======================================================
  Subroutine NormMYCS(c,norm)
    !***
    Implicit None
    Real(8), Intent(In)  :: c(3,3,3)
    Real(8), Intent(Out) :: norm(3)
    Real(8) m1,m2,m(0:3,0:2),t0,t1,t2
    Real(8) mm(3)
    Integer cn
    ! write the plane as: sgn(mx) X =  my Y +  mz Z + alpha 
    !                           m00 X = m01 Y + m02 Z + alpha 

    m1 = c(1,2,1) + c(1,2,3) + c(1,1,2) + c(1,3,2) + &
        c(1,2,2)
    m2 = c(3,2,1) + c(3,2,3) + c(3,1,2) + c(3,3,2) + &
        c(3,2,2)

    if(m1>m2) then
      m(0,0) = 1.
    else
      m(0,0) = -1.
    end if

    m1 = c(1,1,2)+ c(3,1,2)+ c(2,1,2)
    m2 = c(1,3,2)+ c(3,3,2)+ c(2,3,2)
    m(0,1) = 0.5*(m1-m2)

    m1 = c(1,2,1)+ c(3,2,1)+ c(2,2,1)
    m2 = c(1,2,3)+ c(3,2,3)+ c(2,2,3)
    m(0,2) = 0.5*(m1-m2)

    ! write the plane as: sgn(my) Y =  mx X +  mz Z + alpha, 
    !                          m11 Y = m10 X + m12 Z + alpha.
    m1 = c(1,1,2) + c(1,3,2) + c(1,2,2)
    m2 = c(3,1,2) + c(3,3,2) + c(3,2,2)
    m(1,0) = 0.5*(m1-m2)

    m1 = c(2,1,1) + c(2,1,3) + c(3,1,2) + c(1,1,2) +&
        c(2,1,2)
    m2 = c(2,3,1) + c(2,3,3) + c(3,3,2) + c(1,3,2) +&
        c(2,3,2)

    if(m1>m2) then
      m(1,1) = 1.
    else
      m(1,1) = -1.
    end if

    m1 = c(2,1,1)+ c(2,2,1)+ c(2,3,1)
    m2 = c(2,1,3)+ c(2,2,3)+ c(2,3,3)
    m(1,2) = 0.5*(m1-m2)

    m1 = c(1,2,1)+ c(1,2,3)+ c(1,2,2)
    m2 = c(3,2,1)+ c(3,2,3)+ c(3,2,2)
    m(2,0) = 0.5*(m1-m2)

    m1 = c(2,1,1)+ c(2,1,3)+ c(2,1,2)
    m2 = c(2,3,1)+ c(2,3,3)+ c(2,3,2)
    m(2,1) = 0.5*(m1-m2)

    m1 = c(1,2,1) + c(3,2,1) + c(2,1,1) + c(2,3,1) +&
        c(2,2,1)
    m2 = c(1,2,3) + c(3,2,3) + c(2,1,3) + c(2,3,3) +&
        c(2,2,3)

    if(m1>m2) then
      m(2,2) = 1.
    else
      m(2,2) = -1.
    end if

    ! normalize each set (mx,my,mz): |mx|+|my|+|mz| = 1

    t0 = DABS(m(0,0)) + DABS(m(0,1)) + DABS(m(0,2))
    m(0,0) = m(0,0)/t0
    m(0,1) = m(0,1)/t0
    m(0,2) = m(0,2)/t0

    t0 = DABS(m(1,0)) + DABS(m(1,1)) + DABS(m(1,2))
    m(1,0) = m(1,0)/t0
    m(1,1) = m(1,1)/t0
    m(1,2) = m(1,2)/t0

    t0 = DABS(m(2,0)) + DABS(m(2,1)) + DABS(m(2,2))
    m(2,0) = m(2,0)/t0
    m(2,1) = m(2,1)/t0
    m(2,2) = m(2,2)/t0

    ! choose among the three central schemes */ 
    t0 = DABS(m(0,0))
    t1 = DABS(m(1,1))
    t2 = DABS(m(2,2))

    cn = 0
    if (t1 > t0) then
      t0 = t1
      cn = 1
    endif

    if (t2 > t0) cn = 2

    ! Youngs-CIAM scheme */
    Call NormParkerYoungs(c,mm)
    m(3,0:2) = mm(1:3)
    t0 = DABS(m(3,0)) + DABS(m(3,1)) + DABS(m(3,2)) + 1d-20
    m(3,0) = m(3,0)/t0
    m(3,1) = m(3,1)/t0
    m(3,2) = m(3,2)/t0

    t0 = DABS (m(3,0))
    t1 = DABS (m(3,1))
    t2 = DABS (m(3,2))
    if (t1 > t0)  t0 = t1
    if (t2 > t0)  t0 = t2

    ! ! choose between the previous choice and Youngs-CIAM 
    if (DABS(m(cn,cn)) > t0)  cn = 3

    ! components of the normal vector */
    norm(1:3) = m(cn,0:2)

    return
  End Subroutine NormMYCS

  !=======================================================
  ! [1-4] Efficient Least-square volume averaging
  ! * returns normal normalized so that |mx|+|my|+|mz| = 1*
  ! In the ELVIRA reference paper
  !    A Conservative Three-Dimensional Eulerian Method for Coupled Solid–Fluid Shock Capturing
  ! by G.H Miller and P.Colella, 5*5*5 stencil is used,
  ! However, here we use 3*3*3 stencil.
  ! 
  ! For each direction, every 24 posible triplets that contains
  ! The central cell are evaluated.
  ! The symmetric and flip cases are estimated as well.
  ! It is 24 * 6 = 144 in total
  !-------------------------------------------------------
  ! Input: c (vof function, 3*3*3 stencil)
  ! Output: norm (normal vector)
  !=======================================================
  Subroutine NormELVIRA(c,norm)
    Implicit None
    Real(8), Intent(In) :: c(3,3,3)
    Real(8), Intent(Out) :: norm(3)
    Real(sp) :: HVolume(3,3)
    Integer, Dimension(24) :: p1is, p1js, p3is, p3js
    Integer :: p2i = 2
    Integer :: p2j = 2
    Integer :: p1i, p1j, p3i, p3j
    Real(sp) :: block_vol
    Real(sp) :: vol_e(3,3,3)
    Real(sp) :: vec_1(3), vec_2(3)
    Real(sp) :: dxl(3), x0(3)
    Real(sp) :: alpha, norm_e(3)
    Real(sp) :: err, err_min
    Integer :: i, j, k, dir, nelvira, ii, jj

    p1is =(/1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,2,2,2,2,2,2,3,3,3/)
    p1js =(/1,1,1,1,1,1,2,2,2,2,2,3,3,3,3,1,1,1,3,3,3,1,1,2/)
    p3is =(/1,1,2,2,3,3,1,2,2,3,3,2,2,3,3,3,3,3,3,3,3,3,3,3/)
    p3js =(/2,3,1,3,2,1,3,1,3,1,3,1,3,2,3,1,2,3,1,2,3,2,3,3/)

    block_vol = sum(c) / 27.0_sp

    norm = 1.0/3.0_sp
    err_min = 100.0_sp

    Do Dir = 1, 3
      If (dir == 1) Then
        Do j = 1,3
          Do i = 1,3
            HVolume(i,j) = sum(c(:,i,j))
          End Do
        End Do
      Else If (dir == 2) Then
        Do j = 1,3
          Do i = 1,3
            HVolume(i,j) = sum(c(i,:,j))
          End Do
        End Do
      Else
        Do j = 1,3
          Do i = 1,3
            HVolume(i,j) = sum(c(i,j,:))
          End Do
        End Do
      End If
      ! Find the minimum normal voector at the 24 candidates
      ! Note that the flip case of the 24 candidates should be calculated as well
      ! So, actually 48 candidates
      ! When times by 3, corresponds with the 72 or 144 values in ELVIRA paper
      Do nelvira = 1, 24
        p1i = p1is(nelvira)
        p1j = p1js(nelvira)
        p3i = p3is(nelvira)
        p3j = p3js(nelvira)
        If (dir == 1) Then
          vec_1 = (/Hvolume(p1i,p1j)-Hvolume(p2i,p2j), dble(p1i-p2i), dble(p1j-p2j)/)
          vec_2 = (/Hvolume(p3i,p3j)-Hvolume(p2i,p2j), dble(p3i-p2i), dble(p3j-p2j)/)
        Else If (dir == 2) Then
          vec_1 = (/dble(p1i-p2i), Hvolume(p1i,p1j)-Hvolume(p2i,p2j), dble(p1j-p2j)/)
          vec_2 = (/dble(p3i-p2i), Hvolume(p3i,p3j)-Hvolume(p2i,p2j), dble(p3j-p2j)/)
        Else
          vec_1 = (/dble(p1i-p2i), dble(p1j-p2j), Hvolume(p1i,p1j)-Hvolume(p2i,p2j)/)
          vec_2 = (/dble(p3i-p2i), dble(p3j-p2j), Hvolume(p3i,p3j)-Hvolume(p2i,p2j)/)
        End If
        norm_e(1) = vec_1(2) * vec_2(3) - vec_1(3) * vec_2(2)
        norm_e(2) = vec_1(3) * vec_2(1) - vec_1(1) * vec_2(3)
        norm_e(3) = vec_1(1) * vec_2(2) - vec_1(2) * vec_2(1)
        Call Normalization1(norm_e)

        Do ii = 1, 2
          norm_e = - norm_e
        Do jj = 1, 2
          norm_e(dir) = - norm_e(dir)
          alpha =  FloodSZ_backward(norm_e,block_vol)
          dxl = 1.0/3.0_sp
          Do k = 1, 3
            Do j = 1, 3
              Do i = 1, 3
                x0(1) = dble(i-1) / 3.0_sp
                x0(2) = dble(j-1) / 3.0_sp
                x0(3) = dble(k-1) / 3.0_sp
                vol_e(i,j,k) = FloodSZ_Forward(norm_e,alpha,x0,dxl)
              End Do
            End Do
          End Do
          err = norm2( (vol_e-c))
          If ( err < err_min) Then
            norm = norm_e
            err_min = err
          End If
        End Do
        End Do

      End Do
    End Do

    return
  End Subroutine NormELVIRA

  !=======================================================
  ! [2-5] Forward flooding algorithm 
  ! FIND THE "CUT VOLUME" V0
  ! for GIVEN
  !    r0, dr0
  !  and
  !    m1 x1 + m2 x2 + m3 x3 = alpha
  !
  ! Adopted from Sussan's Legacy code
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
  SUBROUTINE FloodSussman_forwardC(nr,alpha,x0,dx,f,xc0)
    Use ModSussman
    Implicit None

    REAL(8), INTENT(IN) :: nr(3),x0(3),dx(3),alpha
    REAL(8), INTENT(OUT) :: f
    REAL(8), INTENT(OUT) :: xc0(3)
    REAL(8) :: nr0(3)

    Real(8) :: cen_zero(3), cen_zero2(3)
    Real(8) :: intercept_sussman
    Real(8) :: f_sussman
    Real(8) :: cen_sussman(3)
    Real(8) :: angle_sussman(2)
    Integer :: default_one = 1
    Real(8) :: maxdx
    Real(8) :: refcentroid(3)
    Real(8) :: dx_scale(3)
    Real(8) :: xsten0_scale(-3:3,3)
    Real(8) :: refcentroid_scale(3)
    Real(8) :: dx0(3)
    Real(8) :: nr2(3)

    f_sussman = 0.0_sp
    cen_sussman = 0.0_sp
    xc0 = 0.0_sp
    f =  FloodSZ_Forward(nr,alpha,x0,dx)

    ! Only one component of dx <1, other two =1
    f_sussman = f / dx(1) / dx(2) / dx(3)
    nr2 = -nr * dx
    Call Normalization2(nr2)
    Call Norm2Angle(angle_sussman,nr2)

    dx0 = 1.0_sp

    ! Evaluate in a unit cell [-0.5,0,5]^3
    if (f_sussman>1.0-epsc) then
      xc0 = dx / 2.0_sp + x0
    elseif (f_sussman<epsc) then
      xc0 = 0.0_sp
    else
      call multi_rotatefunc( &
        sussman%bfact, &
        dx0, &
        sussman%xsten0, &
        sussman%nhalf0, &
        sussman%xtetlist_vof, &
        default_one, &
        sussman%xtetlist_cen, &
        default_one, &
        sussman%nmax, &
        refcentroid_scale,  &
        f_sussman, &
        sussman%continuous_mof, &
        angle_sussman, &
        levelrz, &
        cen_zero2, & 
        intercept_sussman, &
        cen_sussman, & 
        default_one, &
        default_one, &
        sussman%sdim)
      xc0 = ( cen_sussman + 0.5_sp )* dx + x0
    endif

  End Subroutine FloodSussman_forwardC


  !=======================================================
  ! [2-6] Backward flooding algorithm of THINC/SW
  ! FIND THE "jump center" the sigmoid/tanh function
  ! for given vof function
  !    f
  !  and THINC parameter
  !    betagamma2
  !
  !-------------------------------------------------------
  ! Input:
  !      f: vof function
  !      batagamma2: 2 * beta * gamma
  ! Output:
  !      THINC1DForward: coordinate of jump center
  !=======================================================
  Real(8) Function THINC1DBackward(f, betagamma2)

    Implicit None
    REAL(8) :: f, betagamma2
    Real(8) :: a0, a1, a2

    a0 = dexp(betagamma2 * f)
    a1 = dexp(betagamma2) - a0
    a2 = a0 - 1.0_sp + 1.0d-20
    THINC1DBackward = 1.0d0 /  betagamma2 * dlog(dabs(a1 / a2))

  End Function THINC1DBackward

  !=======================================================
  ! [2-7] Forward flooding algorithm of THINC/SW
  ! FIND THE "CUT VOLUME" V0
  ! for given doner region
  !    x0, deltax
  !  and THINC parameter
  !    x_center, betagamma2
  !
  !-------------------------------------------------------
  ! Input:
  !      x_center: center of the jump for sigmoid/Tanh
  !      batagamma2: 2 * beta * gamma
  !      x0: the origin of the direction
  !      deltax: u*dt/dx
  ! Output:
  !      THINC1DForward: vof function
  !=======================================================
  Real(8) Function THINC1DForward(x_center, betagamma2, x0, deltax)

    Implicit None
    REAL(8) :: x_center, betagamma2, x0, deltax
    Real(8) :: a0, a1

    a1 = dexp(betagamma2 * ( x_center - ( x0 + deltax) ) ) + 1.0d0
    a0 = dexp(betagamma2 * ( x_center - x0 ) ) + 1.0d0 + 1.0d-20

    THINC1DForward = deltax + 1.0_sp / betagamma2 * dlog( a1 / a0 )

  End Function THINC1DForward

  !=======================================================
  ! [3-1-1] MOF reconstruction using Gauss-Newton iteration
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
  Subroutine MOFZY2(f, c, norm, Init_Norm)
    ! Use ModOptimizer
    Implicit None
    Real(sp), Intent(In)  :: f
    Real(sp), Intent(In)  :: c(3)
    Real(sp), Intent(Out) :: norm(3)
    Real(sp), Optional, Intent(In) :: Init_Norm(3)

    Real(sp)  :: c_mof(3)
    Real(sp), Dimension(3)   :: norm_2
    Real(sp), Dimension(2)   :: delangle, angle_init, angle_previous, new_angle
    Real(sp), Dimension(2)   :: angle_base, angle_plus, angle_minus
    Real(sp), Dimension(2)   :: err_plus, err_minus
    Real(sp), Dimension(3)   :: dbase, dopt, dp, dm
    Real(sp), Dimension(3,2) :: d_plus, d_minus
    Real(sp), Dimension(3)   :: cenopt, cenp, cenm, cen_init
    Real(sp), Dimension(3,2) :: cen_plus, cen_minus
    Real(sp), Dimension(3,2) :: Jacobian
    Real(sp), Dimension(2,2) :: Hessian, HessianT
    Real(sp), Dimension(2)   :: gradient
    Real(sp)   :: det
    Real(sp), Dimension(3,MOFITERMAX+1) :: d_array, cen_array
    Real(sp), Dimension(2,MOFITERMAX+1) :: angle_array
    Real(sp), Dimension(MOFITERMAX+1)   :: err_array

    Real(sp) :: err, err_local_min
    Integer :: singular_flag
    Integer :: i_angle, j_angle, dir, iter, i, j

    c_mof = c - 0.5_sp

    ! Initialize angle
    If(present(Init_Norm))then
      norm_2 = Init_Norm
    Else
      ! Call Initial_Guess(c_mof-0.5_sp, f, angle_init, c_diff)
      norm_2 = - c_mof
    EndIf

    Call Normalization2(Norm_2)
    Call Norm2Angle(angle_init,norm_2)
    ! Initialize other data
    do dir=1,2
      angle_array(dir,1)= angle_init(dir)
    enddo
    Call FindCentroid(angle_init,f,cen_init)
    do dir=1,3
      d_array(dir,1) = cen_init(dir) - c_mof(dir)
      cen_array(dir,1)=cen_init(dir)
      err_array(1)=err_array(1) + d_array(dir,1) * d_array(dir,1)
    enddo
    ! err_array(1) = sqrt(dot_product(d_array(:,1),d_array(:,1)))
    err_array(1) = dot_product(d_array(:,1),d_array(:,1))

    iter = 0
    err = err_array(1)
    err_local_min = err_array(1)
    mof_niter = 0
    Do While ((iter.lt.MOFITERMAX).and. &
        (err.gt.mof_tol).and. &
        (err_local_min.gt.mof_tol))

      do dir=1,3
        dbase(dir)=d_array(dir,iter+1)
      enddo

      do i_angle=1,2
        angle_base(i_angle)=angle_array(i_angle,iter+1)
      enddo

      ! theta phi shift angle
      ! to calculate the numerical gradient
      Do i_angle=1,2
        Do j_angle=1,2
          angle_plus(j_angle)=angle_base(j_angle)
          angle_minus(j_angle)=angle_base(j_angle)
        End Do
        angle_plus(i_angle)=angle_plus(i_angle)+delta_theta
        angle_minus(i_angle)=angle_minus(i_angle)-delta_theta
        Call FindCentroid(angle_plus,f,cenp)
        err_plus(i_angle)= 0.0_sp
        do dir=1,3
          dp(dir) = c_mof(dir) - cenp(dir)
          err_plus(i_angle)=err_plus(i_angle)+dp(dir)**2
          d_plus(dir,i_angle)=dp(dir)
          cen_plus(dir,i_angle)=cenp(dir)
        end do
        ! err_plus(i_angle)=sqrt(err_plus(i_angle))
        Call FindCentroid(angle_minus,f,cenm)
        err_minus(i_angle)= 0.0_sp
        do dir=1,3
          dm(dir) = c_mof(dir) - cenm(dir)
          err_minus(i_angle)=err_minus(i_angle)+dm(dir)**2
          d_minus(dir,i_angle)=dm(dir)
          cen_minus(dir,i_angle)=cenm(dir)
        enddo
        ! err_minus(i_angle)=sqrt(err_minus(i_angle))

        ! jacobian matrix (numerical partial gradient):
        do dir=1,3
          Jacobian(dir,i_angle)=(cenp(dir)-cenm(dir))/(2.0_sp*delta_theta)
        enddo

      End Do ! end theta phi shift angle

      ! JT * X
      Do i_angle=1,2
        gradient(i_angle) = dot_product(dbase, Jacobian(:,i_angle))
      EndDo


      !! Begin Gauss Newton
      Singular_flag = 0
      Do i=1,2
        Do j=1,2
          Hessian(i,j) = dot_product(Jacobian(:,i), Jacobian(:,j))
        EndDo
      EndDo


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

      Do i_angle = 1,2
        If (delangle(i_angle).gt.delta_theta_max) then
          delangle(i_angle)=delta_theta_max
        Else If (delangle(i_angle).lt.-delta_theta_max) then
          delangle(i_angle)=-delta_theta_max
        End If
      End Do

      !! End Gauss Newton
      ! call GaussNewton(Jacobian, gradient, 3, 2, delangle, Singular_flag)

      ! Find delta angle
      if (singular_flag.eq.0) then
        ! -pi<angle<pi
        call advance_angle(angle_base,delangle)
        do i_angle=1,2
          angle_previous(i_angle)=angle_base(i_angle)
          ! call advance_angle(angle_base(i_angle),delangle(i_angle))
        enddo
        Call FindCentroid(angle_base,f, cenopt)
        do dir=1,3
          dopt(dir) = cenopt(dir) - c_mof(dir)
        enddo

      else if (singular_flag.eq.1) then
        err_local_min= 0.0_sp

        do dir=1,3
          dopt(dir)=dbase(dir)
        enddo
        do i_angle=1,2
          ! angle_base(i_angle)=angle_array(i_angle,iter+1)
          angle_array(i_angle,iter+1)=angle_base(i_angle)
        enddo
      else
        print *,"singular_flag invalid"
        stop
      End If

      err= 0.0_sp
      do dir=1,3
        err=err+dopt(dir)**2
      enddo
      ! err=sqrt(err)

      if (singular_flag.eq.1) then
        ! do nothing
      else if (singular_flag.eq.0) then
        do i_angle=1,2
          do j_angle=1,2
            angle_plus(j_angle)=angle_previous(j_angle)
            angle_minus(j_angle)=angle_previous(j_angle)
          enddo
          angle_plus(i_angle)=angle_plus(i_angle)+delta_theta
          angle_minus(i_angle)=angle_minus(i_angle)-delta_theta

          if ((err.le.err_plus(i_angle)).and. &
              (err.le.err_minus(i_angle))) then
            ! do nothing
          else

            if (err.ge.err_plus(i_angle)) then
              err=err_plus(i_angle)
              do dir=1,3
                dopt(dir)=d_plus(dir,i_angle)
                cenopt(dir)=cen_plus(dir,i_angle)
              enddo
              do j_angle=1,2
                angle_base(j_angle)=angle_plus(j_angle)
              enddo
            endif

            if (err.ge.err_minus(i_angle)) then
              err=err_minus(i_angle)
              do dir=1,3
                dopt(dir)=d_minus(dir,i_angle)
                cenopt(dir)=cen_minus(dir,i_angle)
              enddo
              do j_angle=1,2
                angle_base(j_angle)=angle_minus(j_angle)
              enddo
            endif

          endif ! use safe update

        enddo ! i_angle     if (singular_flag.eq.1) then
      else
        print *,"singular_flag invalid"
        stop
      endif

      do dir=1,3
        d_array(dir,iter+2)=dopt(dir)
        cen_array(dir,iter+2)=cenopt(dir)
      enddo

      do i_angle=1,2
        angle_array(i_angle,iter+2)=angle_base(i_angle)
      enddo
      err_array(iter+2)=err / 2.0_sp
      iter=iter+1
    End Do

    mof_niter(1) = iter + 1
#if defined(DEBUG)
    Block
      Integer :: ii
      Do ii = 1, mof_niter(1)
        print *, '=====step',ii-1,'=========='
        print *, cen_array(:,ii)
        print *, angle_array(:,ii)
        print *, err_array(ii)
        print *, err, err_local_min
      End Do
    End Block
#endif

    ! print *, niter

    do dir=1,2
      new_angle(dir)=angle_array(dir,iter+1)
    enddo

    ! Convert angle to norm
    call Angle2Norm(new_angle,norm)
    Call Normalization1(norm)
#if defined(DEBUG)
    print *, '=====MOF Reconc norm=========='
    print *, 'Normal vector:', norm
    print *, '=============================='
#endif
  End Subroutine MOFZY2


  !=======================================================
  ! [3-2] MOF reconstruction using Analytic gradient and Gauss-Newton iteration
  !-------------------------------------------------------
  ! The algorithm uses analytic slope by Lemoine (2020, JCP)
  ! Gauss-Newton iteration is used
  ! Input:
  !      f: vof function
  !      c: centroid (cx, cy, cz)
  ! Optional Input:
  !      init_norm: Initial guess of the normal vector
  ! Output:
  !      norm: vof function
  !=======================================================
  Subroutine MOFLemoine_GaussNewton(f, c, norm, Init_Norm)
    ! Use ModOptimizer
    Use mod_mof3d_analytic_centroid
    Implicit None
    Real(sp), Intent(In)  :: f
    Real(sp), Intent(In)  :: c(3)
    Real(sp), Intent(Out) :: norm(3)
    Real(sp), Optional, Intent(In) :: Init_Norm(3)

    Real(sp)  :: vof
    Real(sp)  :: c_mof(3)
    Real(sp), Dimension(3)   :: norm_2
    Real(sp), Dimension(2)   :: delangle, angle_init, new_angle
    Real(sp), Dimension(2)   :: angle_base
    Real(sp), Dimension(2,MOFITERMAX+1) :: angle_array
    Real(sp), Dimension(MOFITERMAX+1)   :: err_array

    Real(sp) :: err
    Integer :: singular_flag
    Integer :: i_angle, dir, iter, i, j
    real(sp) :: dxs(3)
    Real(sp) :: c_diff(3)
    Real(sp) :: gradient(2)
    Real(sp) :: Jacobian(3,2)
    Real(sp) :: Hessian(2,2)
    Real(sp) :: HessianT(2,2)
    Real(sp) :: det

    dxs = 1.0_sp

    if (f .ge. 0.5_sp) then
      vof = 1.0_sp - f
      c_mof = ( 0.5 - c * f ) / vof
    else
      vof = f
      c_mof = c
    endif

    ! Initialize angle
    If(present(Init_Norm))then
      norm_2 = Init_Norm
      Call Normalization2(Norm_2)
      Call Norm2Angle(angle_init,norm_2)
      Call mof3d_compute_analytic_gradient_GN(angle_init, c_mof, vof, dxs, c_diff, Jacobian)
      err = norm2(c_diff)
    Else
      Call Initial_Guess(c_mof-0.5_sp, vof, angle_init, err)
    EndIf

    do dir=1,2
      angle_array(dir,1)= angle_init(dir)
    enddo

    do dir=1,3
      err_array(1) = err
    enddo

    iter = 0
    err = err_array(1)
    mof_niter = 0
    Do While ((iter.lt.MOFITERMAX) .and. (err.gt.mof_tol))

      Do i_angle=1,2
        angle_base(i_angle) = angle_array(i_angle, iter+1)
      End Do

      Call mof3d_compute_analytic_gradient_GN(angle_base, c_mof, vof, dxs, c_diff, Jacobian)
      gradient = 0.0_sp
      hessian = 0.0_sp
      Singular_flag = 0
      Do i=1,2
        Do j=1,2
          Hessian(i,j) = dot_product(Jacobian(:,i), Jacobian(:,j))
        EndDo
        gradient(i) = dot_product(c_diff, Jacobian(:,i))
      EndDo
      err = dot_product(c_diff, c_diff)

      det = Hessian(1,1)*Hessian(2,2) - Hessian(1,2)*Hessian(2,1)
      ! Calculate the inverse of the matrix
      HessianT(1,1) = + Hessian(2,2) / det
      HessianT(2,1) = - Hessian(2,1) / det
      HessianT(1,2) = - Hessian(1,2) / det
      HessianT(2,2) = + Hessian(1,1) / det

      If (det .lt. det_lim) Then
        Singular_flag = 1
      End If

      ! Delta angle
      Do i=1,2
        delangle(i) = - dot_product(HessianT(:,i), gradient)
      End Do
      !! End Gauss Newton

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

      call Angle2Norm(angle_base,norm)
      Call Normalization1(norm)
      ! print *, norm
      err_array(iter+2)=err
      iter=iter+1

      if ( dot_product(delangle,delangle) .lt. mof_tol_dangle) exit
    End Do
    mof_niter(1) = iter+1 

    do dir=1,2
      new_angle(dir)=angle_array(dir,iter+1)
    enddo
    call Angle2Norm(new_angle,norm)
    Call Normalization1(norm)

    if (f .ge. 0.5_sp) norm = -norm

  End Subroutine MOFLemoine_GaussNewton


  !=======================================================
  ! [3-3] MOF reconstruction using analytic gradient and BFGS 
  !     For given vof and centroid, find the best norm
  !-------------------------------------------------------
  ! Input:
  !      f: vof function
  !      c: centroid (cx, cy, cz)
  ! Optional Input:
  !      init_norm: Initial guess of the normal vector
  ! Output:
  !      norm: vof function
  !=======================================================
  Subroutine MOFLemoine_BFGS(f, c, norm, Init_Norm)
    ! Use ModOptimizer
   use variables_mof
   Use mod_cg3_polyhedron
    Use mod_mof3d_bfgs
    Implicit None
    Real(sp), Intent(In)  :: f
    Real(sp), Intent(In)  :: c(3)
    Real(sp), Intent(Out) :: norm(3)
    Real(sp), Optional, Intent(In) :: Init_Norm(3)

    Real(sp) :: vof
    Real(sp) :: c_mof(3), c_mof_sym(3)
    Real(sp), Dimension(3)   :: norm_2
    Real(sp), Dimension(2)   :: angle_init
    Real(sp) :: err

    Real(sp) :: dxs(3)
    Integer :: nstat(2)
    Real(sp) :: residual(2)

    dxs = 1.0_sp

    if (f .ge. 0.5_sp) then
      vof = 1.0_sp - f
      c_mof = ( 0.5 - c * f ) / vof
    else
      vof = f
      c_mof = c
    endif
    c_mof_sym = ( 0.5 - c_mof * vof ) / (1.0_sp - vof)

    ! Initialize angle
    If(present(Init_Norm))then
      norm_2 = Init_Norm
      Call Normalization2(Norm_2)
      Call Norm2Angle(angle_init,norm_2)
    Else
      Call Initial_GuessOld(c_mof-0.5_sp, vof, angle_init, err)
    EndIf

    call mof3d_bfgs(LemoinePoly, c_mof, c_mof_sym, vof, angle_init, norm, nstat, residual)

    mof_niter = nstat 

    Call Normalization1(norm)

    if (f .ge. 0.5_sp) norm = -norm

  End Subroutine MOFLemoine_BFGS

  !=======================================================
  ! [3-4] MOF reconstruction using numerical gradient and BFGS 
  !     For given vof and centroid, find the best norm
  !    Sussman's original MOF version
  !    uses computational geometry algorithm to find the centroid
  !-------------------------------------------------------
  ! Input:
  !      vof: vof function
  !      centroid: centroid (cx, cy, cz)
  ! Optional Input:
  !      init_norm: Initial guess of the normal vector
  ! Output:
  !      norm: vof function
  !=======================================================
  Subroutine MOFSussmanGaussNewton(vof, centroid, norm, init_norm)
    Use ModSussman
    Implicit None
    Real(8), Intent(In)  :: vof
    Real(8), Intent(In)  :: centroid(3)
    Real(8), Intent(Out) :: norm(3)
    Real(8), Intent(In), optional :: init_norm(3)

    Real(8)  :: c_mof(3)
    Real(8), Dimension(-1:1,-1:1,-1:1,1) :: ls_mof, lsnormal
    Real(8) :: norm_2(3)
    Integer :: lsnormal_valid(1)
    Real(8) :: npredict(3)
    Real(8) :: angle_init(2), err
    Real(8) :: intercept
    Integer :: default_one = 1

    c_mof = centroid - 0.5_sp

    If(present(Init_Norm))then
      norm_2 = Init_Norm
      Call Normalization2(Norm_2)
    Else
      Call Initial_Guess(c_mof-0.5_sp, vof, angle_init, err)
      Call Angle2Norm(angle_init,norm_2)
    EndIf

    Norm_2 = -norm_2

    lsnormal_valid = 1

    npredict = Norm_2
    mof_niter = 0

    call find_cut_geom_slope( &
        ls_mof, &
        lsnormal, &
        lsnormal_valid, &
        sussman%bfact, &
        sussman%dx, &
        sussman%xsten0, &
        sussman%nhalf0, &
        c_mof,&
        vof, &
        levelrz, &
        npredict, &
        sussman%continuous_mof, &
        norm, &
        intercept, &
        sussman%xtetlist_vof, &
        default_one, &
        sussman%xtetlist_cen, &
        default_one, &
        sussman%multi_centroidA, &
        sussman%nmax, &
        default_one, &
        default_one, &
        default_one, &
        sussman%sdim)

    norm(1:3) = -norm
    mof_niter(1) = nn_iter + 1
    mof_niter(2) = 0

  End Subroutine MOFSussmanGaussNewton

  !=======================================================
  ! [3-5-1] MOF reconstruction using Machine Learning
  !     For given vof and centroid, find the best norm
  !-------------------------------------------------------
  ! Input:
  !      f: vof function
  !      c: centroid (cx, cy, cz)
  ! Optional Input:
  !      init_norm: Initial guess of the normal vector
  ! Output:
  !      norm: vof function
  !=======================================================
  Subroutine MOFNN(f, c, norm, Init_Norm)
    use variables_mof
    Use mod_cg3_polyhedron
    Use mod_mof3d_bfgs
    Implicit None
    Real(sp), Intent(In)  :: f
    Real(sp), Intent(In)  :: c(3)
    Real(sp), Intent(Out) :: norm(3)
    Real(sp), Optional, Intent(In) :: Init_Norm(3)

    Real(sp) :: vof
    Real(sp) :: c_mof(3)
    Real(sp), Dimension(2)   :: angle_init
    Real(sp), Dimension(2)   :: delta_angle_correct
    Real(sp), Dimension(2)   :: angle_correct
    Real(sp), Dimension(3)   :: inputs
    Real(sp) :: err
    vof = f
    c_mof = c

    Call Initial_GuessNN(c_mof-0.5_sp, vof, angle_init, err)
    Call Angle2Norm(angle_init, norm)

    Call Normalization1(norm)

  End Subroutine MOFNN

  !=======================================================
  ! [3-5-1] MOF reconstruction using Nueral Network- stablized
  !  Using NN value as initial guess and do 1 iteration.
  !-------------------------------------------------------
  !      f: vof function
  !      c: centroid (cx, cy, cz)
  ! Optional Input:
  !      init_norm: Initial guess of the normal vector
  ! Output:
  !      norm: vof function
  !=======================================================
  Subroutine MOFNNStab(f, c, norm, Init_Norm)
    ! Use ModOptimizer
    Use mod_mof3d_analytic_centroid
    Implicit None
    Real(sp), Intent(In)  :: f
    Real(sp), Intent(In)  :: c(3)
    Real(sp), Intent(Out) :: norm(3)
    Real(sp), Optional, Intent(In) :: Init_Norm(3)

    Real(sp)  :: vof
    Real(sp)  :: c_mof(3)
    Real(sp), Dimension(3)   :: norm_2
    Real(sp), Dimension(2)   :: delangle, angle_init, new_angle
    Real(sp), Dimension(2)   :: angle_base
    Real(sp), Dimension(2,MOFITERMAX+1) :: angle_array
    Real(sp), Dimension(MOFITERMAX+1)   :: err_array

    Real(sp) :: err
    Integer :: singular_flag
    Integer :: i_angle, dir, iter, i, j
    real(sp) :: dxs(3)
    Real(sp) :: c_diff(3)
    Real(sp) :: gradient(2)
    Real(sp) :: Jacobian(3,2)
    Real(sp) :: Hessian(2,2)
    Real(sp) :: HessianT(2,2)
    Real(sp) :: det

    dxs = 1.0_sp

    if (f .ge. 0.5_sp) then
      vof = 1.0_sp - f
      c_mof = ( 0.5 - c * f ) / vof
    else
      vof = f
      c_mof = c
    endif

    ! Initialize angle
    If(present(Init_Norm))then
      norm_2 = Init_Norm
      Call Normalization2(Norm_2)
      Call Norm2Angle(angle_init,norm_2)
      Call mof3d_compute_analytic_gradient_GN(angle_init, c_mof, vof, dxs, c_diff, Jacobian)
      err = norm2(c_diff)
    Else
      Call Initial_GuessNN(c_mof-0.5_sp, vof, angle_init, err)
    EndIf

    do dir=1,2
      angle_array(dir,1)= angle_init(dir)
    enddo

    do dir=1,3
      err_array(1) = err
    enddo

    iter = 0
    err = err_array(1)
    mof_niter = 0
    mofitermax = 1
    Do While ((iter.lt.MOFITERMAX) .and. (err.gt.mof_tol))

      Do i_angle=1,2
        angle_base(i_angle) = angle_array(i_angle, iter+1)
      End Do

      Call mof3d_compute_analytic_gradient_GN(angle_base, c_mof, vof, dxs, c_diff, Jacobian)
      gradient = 0.0_sp
      hessian = 0.0_sp
      Singular_flag = 0
      Do i=1,2
        Do j=1,2
          Hessian(i,j) = dot_product(Jacobian(:,i), Jacobian(:,j))
        EndDo
        gradient(i) = dot_product(c_diff, Jacobian(:,i))
      EndDo
      err = dot_product(c_diff, c_diff)

      det = Hessian(1,1)*Hessian(2,2) - Hessian(1,2)*Hessian(2,1)
      ! Calculate the inverse of the matrix
      HessianT(1,1) = + Hessian(2,2) / det
      HessianT(2,1) = - Hessian(2,1) / det
      HessianT(1,2) = - Hessian(1,2) / det
      HessianT(2,2) = + Hessian(1,1) / det

      If (det .lt. det_lim) Then
        Singular_flag = 1
      End If

      ! Delta angle
      Do i=1,2
        delangle(i) = - dot_product(HessianT(:,i), gradient)
      End Do
      !! End Gauss Newton

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

      call Angle2Norm(angle_base,norm)
      Call Normalization1(norm)
      ! print *, norm
      err_array(iter+2)=err
      iter=iter+1

      if ( dot_product(delangle,delangle) .lt. mof_tol_dangle) exit
    End Do
    mof_niter(1) = iter+1 

    do dir=1,2
      new_angle(dir)=angle_array(dir,iter+1)
    enddo
    call Angle2Norm(new_angle,norm)
    Call Normalization1(norm)

    if (f .ge. 0.5_sp) norm = -norm

  End Subroutine MOFNNStab

  ! subroutine direction_to_spherical_angles(direction, angles)
  !   double precision, dimension(3), intent(in) :: direction
  !   double precision, dimension(2), intent(out) :: angles

  !   double precision, parameter :: PI = 2d0*acos(0d0)

  !   if ((abs(direction(3)) - 1d0) < epsilon(1d0)) then
  !     angles = [atan2(direction(2), direction(1)), acos(direction(3))]
  !   else
  !     if (direction(3) > 0d0) then
  !       angles = [0d0, 0d0]
  !     else
  !       angles = [0d0, PI]
  !     end if
  !   end if
  ! end subroutine direction_to_spherical_angles

  !=======================================================
  ! [4-6-2] Find the initial guess (Lemoine)
  !-------------------------------------------------------
  ! Input:  
  !      c_norm: centroid ranging from [0,1] (cx,cy,cz)
  !      vof: volume fraction ranging from [0,1]
  ! Output: 
  !      angle: initial guess of angle (theta, phi)
  !      error: initial error angle (angle, theta, phi)
  !=======================================================
  Subroutine Initial_GuessOld(C_mof, vof, angle, error)
    Implicit None
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
    call findcentroid(angle,vof,cen)
    error = norm2(cen-c_mof)

  End Subroutine Initial_GuessOld

  !=======================================================
  ! [4-6-3] Find the initial guess (NN)
  !-------------------------------------------------------
  ! Input:  
  !      c_norm: centroid ranging from [0,1] (cx,cy,cz)
  !      vof: volume fraction ranging from [0,1]
  ! Output: 
  !      angle: initial guess of angle (theta, phi)
  !      error: initial error angle (angle, theta, phi)
  !=======================================================
  Subroutine Initial_GuessNN(C_mof, vof, angle, error)
    Implicit None
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

    inputs(1:2) = angle1
    inputs(3) = vof
    angle = angle1 + ml%predict(inputs)

    error = 1.0_sp

  End Subroutine Initial_GuessNN




  !=======================================================
  ! (4-7) Initilize Nueral Network
  !=======================================================
  Subroutine Initialize_NN
    implicit none
    Call ML%Initialization
  End Subroutine Initialize_NN  

  !=======================================================
  ! (4-8) Initilize Lemoine's algorithm
  !=======================================================
  Subroutine Lemoine_create_cuboid(c, poly)
    double precision, dimension(3), intent(in) :: c
    type(t_polyhedron), intent(out) :: poly
    integer :: error_id

    poly%nb_points = 8
    poly%nb_faces = 6

    allocate(poly%point(3, poly%nb_points))

    poly%point(:,1) = [0.0d0, 0.0d0, 0.0d0]
    poly%point(:,2) = [0.0d0, 0.0d0, c(3) ]
    poly%point(:,3) = [0.0d0, c(2) , 0.0d0]
    poly%point(:,4) = [0.0d0, c(2) , c(3) ]
    poly%point(:,5) = [c(1) , 0.0d0, 0.0d0]
    poly%point(:,6) = [c(1) , 0.0d0, c(3) ]
    poly%point(:,7) = [c(1) , c(2) , 0.0d0]
    poly%point(:,8) = [c(1) , c(2) , c(3) ]

    allocate(poly%face(poly%nb_faces))

    poly%face(1)%size = 4
    allocate(poly%face(1)%id(poly%face(1)%size))
    poly%face(1)%id(1:4) = [8, 4, 2, 6]

    poly%face(2)%size = 4
    allocate(poly%face(2)%id(poly%face(2)%size))
    poly%face(2)%id(1:4) = [8, 6, 5, 7]

    poly%face(3)%size = 4
    allocate(poly%face(3)%id(poly%face(3)%size))
    poly%face(3)%id(1:4) = [8, 7, 3, 4]

    poly%face(4)%size = 4
    allocate(poly%face(4)%id(poly%face(4)%size))
    poly%face(4)%id(1:4) = [4, 3, 1, 2]

    poly%face(5)%size = 4
    allocate(poly%face(5)%id(poly%face(5)%size))
    poly%face(5)%id(1:4) = [1, 3, 7, 5]

    poly%face(6)%size = 4
    allocate(poly%face(6)%id(poly%face(6)%size))
    poly%face(6)%id(1:4) = [2, 1, 5, 6]

    call cg3_complete_polyhedron_structure(poly, error_id)
  end subroutine Lemoine_create_cuboid

  !=======================================================
  ! (5-1) Smooth vof function
  !=======================================================
  Subroutine Smooth_Phi(Phi1A, Phi_DS, nx, ny, nz)
    Implicit None
    real(sp), intent(in), dimension(0:,0:,0:) :: Phi1A
    real(sp), intent(out), dimension(0:,0:,0:) :: Phi_DS
    Integer, Intent(In) :: nx, ny, nz
    Integer :: i,j,k

    Do k=1,nz
      Do j=1,ny
        Do i=1,nx
          PHI_DS(i,j,k) = 27.d0/64.d0 * PHI1A(i,j,k)+ &
              9.d0/128.d0  * ( PHI1A(i,j,k+1) + PHI1A(i,j,k-1) + PHI1A(i,j+1,k) + PHI1A(i,j-1,k) + &
              PHI1A(i+1,j,k) + PHI1A(i-1,j,k) ) + &
              3.d0/256.d0  * ( PHI1A(i,j+1,k+1) + PHI1A(i,j+1,k-1) + PHI1A(i,j-1,k+1) + &
              PHI1A(i,j-1,k-1) + &
              PHI1A(i+1,j,k+1) + PHI1A(i+1,j,k-1) + PHI1A(i-1,j,k+1) + PHI1A(i-1,j,k-1) + &
              PHI1A(i+1,j+1,k) + PHI1A(i+1,j-1,k) + PHI1A(i-1,j+1,k) + PHI1A(i-1,j-1,k)  )  + &
              1.d0/512.d0 * ( PHI1A(i+1,j+1,k+1) + PHI1A(i+1,j+1,k-1) + PHI1A(i+1,j-1,k+1) + &
              PHI1A(i+1,j-1,k-1) + PHI1A(i-1,j+1,k+1) + PHI1A(i-1,j+1,k-1) + PHI1A(i-1,j-1,k+1) + &
              PHI1A(i-1,j-1,k-1) )
        Enddo
      Enddo
    Enddo

  End Subroutine Smooth_Phi

  !===========================================
  ! (5-3) Fast sweeping method (FSM) for ls function for far field
  !===========================================
  SUBROUTINE Get_LS_3D_init(f,ls,dx,dy,dz,nx,ny,nz)
    !============================================================
    use ModTools, only : Heaviside, Dirac
    Use ModGlobal, only : Phi_bc
    Implicit None
    Real(8), Intent(In) :: dx, dy, dz
    integer , Intent(In)  ::  nx,ny,nz
    REAL(8), intent(in) :: f(0:nx+1,0:ny+1,0:nz+1)
    REAL(8), intent(out) :: ls(0:nx+1,0:ny+1,0:nz+1)
    Integer :: s(nx,ny,nz)
    REAL(8) :: ls_v1(0:nx+1,0:ny+1,0:nz+1)
    INTEGER :: I,J,K,L
    INTEGER   ::  Istart, Istep, Iend
    integer   ::  Jstart, Jstep, Jend
    integer   ::  Kstart, Kstep, Kend
    real(8) :: hs, dr, h, a, b, c
    real(8) :: dnew
    !============================================================
    !ÖÃÁã¡ª¡ª½çÃæÍø¸ñ×´Ì¬¡¢LSº¯ÊýºÍÐéÄâLSº¯Êý

    s     = 0.d0
    ls    = 0.d0
    ls_v1 = 0.d0
    !============================================================
    !³õÊ¼»¯ÐéÄâLSº¯Êý
    DO K=1,NZ
      DO I=1,NX
        DO J=1,NY
          IF (f(I,J,K) .ge. 0.5d0) THEN
            ls_v1(I,J,K) = (NX*dx+NY*dy+NZ*dZ)
          ELSE
            ls_v1(I,J,K) = -(NX*dx+NY*dy+NZ*NZ)
          ENDIF
        ENDDO
      ENDDO
    ENDDO
    !Boundary
    Call Phi_bc%setbcs(ls_v1)
    !============================================================
    !¸ù¾ÝÐéÄâLSº¯ÊýÉèÖÃÍø¸ñ½çÃæ×´Ì¬
    do k=1,nz
      DO I=1,NX
        DO J=1,NY
          IF ( ls_v1(I,J,k)*ls_v1(I-1,J,k)   .le. 0.0d0 .OR. &
              &         ls_v1(I,J,k)*ls_v1(I+1,J,k)   .le. 0.0d0 .OR. &
              &         ls_v1(I,J,k)*ls_v1(I,J-1,k)   .le. 0.0d0 .OR. &
              &         ls_v1(I,J,k)*ls_v1(I,J+1,k)   .le. 0.0d0 .OR. &
              &         ls_v1(I,J,k)*ls_v1(I-1,J-1,k) .le. 0.0d0 .OR. &
              &         ls_v1(I,J,k)*ls_v1(I-1,J+1,k) .le. 0.0d0 .OR. &
              &         ls_v1(I,J,k)*ls_v1(I+1,J-1,k) .le. 0.0d0 .OR. &
              &         ls_v1(I,J,k)*ls_v1(I+1,J+1,k) .le. 0.0d0 .OR. &
              &         ls_v1(I,J,k)*ls_v1(I-1,J,k+1)   .le. 0.0d0 .OR. &
              &         ls_v1(I,J,k)*ls_v1(I+1,J,k+1)   .le. 0.0d0 .OR. &
              &         ls_v1(I,J,k)*ls_v1(I,J-1,k+1)   .le. 0.0d0 .OR. &
              &         ls_v1(I,J,k)*ls_v1(I,J+1,k+1)   .le. 0.0d0 .OR. &
              &         ls_v1(I,J,k)*ls_v1(I-1,J-1,k+1) .le. 0.0d0 .OR. &
              &         ls_v1(I,J,k)*ls_v1(I-1,J+1,k+1) .le. 0.0d0 .OR. &
              &         ls_v1(I,J,k)*ls_v1(I+1,J-1,k+1) .le. 0.0d0 .OR. &
              &         ls_v1(I,J,k)*ls_v1(I+1,J+1,k+1) .le. 0.0d0 .OR. &
              &         ls_v1(I,J,k)*ls_v1(I-1,J,k-1)   .le. 0.0d0 .OR. &
              &         ls_v1(I,J,k)*ls_v1(I+1,J,k-1)   .le. 0.0d0 .OR. &
              &         ls_v1(I,J,k)*ls_v1(I,J-1,k-1)   .le. 0.0d0 .OR. &
              &         ls_v1(I,J,k)*ls_v1(I,J+1,k-1)   .le. 0.0d0 .OR. &
              &         ls_v1(I,J,k)*ls_v1(I-1,J-1,k-1) .le. 0.0d0 .OR. &
              &         ls_v1(I,J,k)*ls_v1(I-1,J+1,k-1) .le. 0.0d0 .OR. &
              &         ls_v1(I,J,k)*ls_v1(I+1,J-1,k-1) .le. 0.0d0 .OR. &
              &         ls_v1(I,J,k)*ls_v1(I+1,J+1,k-1) .le. 0.0d0 .OR. &
              &         ls_v1(I,J,k)*ls_v1(I,J,k+1) .le. 0.0d0 .OR. &
              &         ls_v1(I,J,k)*ls_v1(I,J,k-1) .le. 0.0d0  ) THEN
            s(I,J,k) = 1
          ELSE
            s(I,J,k) = 0
          ENDIF
        ENDDO
      ENDDO
    enddo
    !============================================================
    !¸üÐÂ½çÃæÍø¸ñ×´Ì¬Îª"1"(½çÃæÍø¸ñ)µÄÐéÄâLSº¯Êý
    do k=1,nz
      DO I=1,NX
        DO J=1,NY

          IF ( s(I,J,k) .eq. 1 ) THEN
            ls_v1(I,J,k) = 0.0d0
            hs = Heaviside(ls_v1(I,J,k),dsqrt(2.0d0)*dx)
            DO WHILE (dabs(hs-f(I,J,k)) .gt. 1.0d-15)
              hs = Heaviside(ls_v1(I,J,k),dsqrt(2.0d0)*dx)
              dr = Dirac(ls_v1(I,J,k),dsqrt(2.0d0)*dx)
              ls_v1(I,J,k) = ls_v1(I,J,k)-(hs-f(I,J,k))/(dr+1.0d-15)
            ENDDO
          ENDIF
        ENDDO
      ENDDO
    enddo


    Call Phi_bc%setbcs(ls_v1)
    !============================================================
    !¸üÐÂ½çÃæÍø¸ñ×´Ì¬Îª"0"(·Ç½çÃæÍø¸ñ)µÄÐéÄâLSº¯Êý
    h = dx
    !-------------------------------------
    !X+,Y+·½ÏòÉ¨Ãè
    DO l=1,8

      if (l.eq.1) then
        !X+,Y+,Z+·½ÏòÉ¨Ãè
        Istart = 1
        Istep  = 1
        Iend   = nx

        Jstart = 1
        Jstep  = 1
        Jend   = ny

        Kstart = 1
        Kstep  = 1
        Kend   = nz
      elseif (l.eq.2) then
        !X-,Y+,Z+·½ÏòÉ¨Ãè
        Istart = nx
        Istep  = -1
        Iend   = 1

        Jstart = 1
        Jstep  = 1
        Jend   = ny

        Kstart = 1
        Kstep  = 1
        Kend   = nz
      elseif (l.eq.3) then
        !X+,Y-,Z+·½ÏòÉ¨Ãè
        Istart = 1
        Istep  = 1
        Iend   = nx

        Jstart = ny
        Jstep  = -1
        Jend   = 1

        Kstart = 1
        Kstep  = 1
        Kend   = nz
      elseif  (l.eq.4) then
        !X-,Y-,Z+·½ÏòÉ¨Ãè
        Istart = nx
        Istep  = -1
        Iend   = 1

        Jstart = ny
        Jstep  = -1
        Jend   = 1

        Kstart = 1
        Kstep  = 1
        Kend   = nz
      ELSEif (l.eq.5) then
        !X+,Y+,Z-·½ÏòÉ¨Ãè
        Istart = 1
        Istep  = 1
        Iend   = nx

        Jstart = 1
        Jstep  = 1
        Jend   = ny

        Kstart = nz
        Kstep  = -1
        Kend   = 1
      elseif (l.eq.6) then
        !X-,Y+,Z-·½ÏòÉ¨Ãè
        Istart = nx
        Istep  = -1
        Iend   = 1

        Jstart = 1
        Jstep  = 1
        Jend   = ny

        Kstart = nz
        Kstep  = -1
        Kend   = 1
      elseif (l.eq.7) then
        !X+,Y-,Z-·½ÏòÉ¨Ãè
        Istart = 1
        Istep  = 1
        Iend   = nx

        Jstart = ny
        Jstep  = -1
        Jend   = 1

        Kstart = nz
        Kstep  = -1
        Kend   = 1
      else 
        !X-,Y-,Z-·½ÏòÉ¨Ãè
        Istart = nx
        Istep  = -1
        Iend   = 1

        Jstart = ny
        Jstep  = -1
        Jend   = 1

        Kstart = nz
        Kstep  = -1
        Kend   = 1
      endif

      DO K=KSTART,KEND,KSTEP
        DO J=Jstart,Jend,Jstep
          DO I=Istart,Iend,Istep

            !ls>=0
            IF ( ls_v1(I,J,K) .ge. 0.0d0 ) THEN
              a = dmin1(ls_v1(I-1,J,K),ls_v1(I+1,J,K))
              b = dmin1(ls_v1(I,J-1,K),ls_v1(I,J+1,K))
              c = dmin1(ls_v1(I,J,K-1),LS_V1(I,J,K+1))
              call rank3(a,b,c)

              if ( (a+h).lt.b ) then
                dnew = a + h 
              else 
                dnew = (a+b+dsqrt(2.0*h*h-(a-b)**2))/2.0d0
                IF ( dnew .gt. c ) THEN
                  dnew = (a+b+c+dsqrt(3.0*h*h-(a-b)**2-(a-c)**2-(b-c)**2))/3.0d0

                endif
              endif

              IF ( s(I,J,k) .eq. 0 ) THEN
                ls_v1(I,J,k) = dmin1(ls_v1(I,J,k),dnew)
              ENDIF
              !ls<0
            ELSE
              a = dmax1(ls_v1(I-1,J,k),ls_v1(I+1,J,k))
              b = dmax1(ls_v1(I,J-1,k),ls_v1(I,J+1,k))
              c = dmax1(ls_v1(I,J,k-1),ls_v1(I,J,k+1))
              call rank3(a,b,c)

              if( (c-h).gt.b ) then
                dnew = c - h
              else
                dnew = (b+c-dsqrt(2.0*h*h-(b-c)**2))/2.0d0
                IF ( dnew .lt. a ) THEN
                  dnew = (a+b+c-dsqrt(3.0*h*h-(a-b)**2-(a-c)**2-(b-c)**2))/3.0d0
                endif
              endif

              IF ( s(I,J,k) .eq. 0 ) THEN
                ls_v1(I,J,k) = dmax1(ls_v1(I,J,k),dnew)
              ENDIF
            ENDIF
          ENDDO
        ENDDO
      ENDDO
      Call Phi_bc%setbcs(ls_v1)
    ENDDO




    !============================================================
    !¸üÐÂLSº¯Êý
    ls(1:NX,1:NY,1:nz)=LS_V1(1:NX,1:NY,1:nz)
    !============================================================
    RETURN

  END subroutine GET_LS_3D_INIT

  !===========================================
  ! (5-3) Fast sweeping method (FSM) for ls function for far field
  !===========================================
  SUBROUTINE VoSET(f33,dx,dy,dz,nx,ny,nzls)
    !============================================================
    use ModTools, only : Heaviside, Dirac
    Use ModGlobal, only : Phi_bc
    Implicit None
    Real(8), Intent(In) :: dx, dy, dz
    integer , Intent(In)  ::  nx,ny,nz
    REAL(8), intent(in) :: f(0:nx+1,0:ny+1,0:nz+1)
    REAL(8), intent(out) :: ls(0:nx+1,0:ny+1,0:nz+1)
    Integer :: s(nx,ny,nz)
    REAL(8) :: ls_v1(0:nx+1,0:ny+1,0:nz+1)
    INTEGER :: I,J,K,L
    INTEGER   ::  Istart, Istep, Iend
    integer   ::  Jstart, Jstep, Jend
    integer   ::  Kstart, Kstep, Kend
    real(8) :: hs, dr, h, a, b, c
    real(8) :: dnew
    !============================================================
    !ÖÃÁã¡ª¡ª½çÃæÍø¸ñ×´Ì¬¡¢LSº¯ÊýºÍÐéÄâLSº¯Êý

    Do 
    s     = 0.d0
    ls    = 0.d0
    ls_v1 = 0.d0
    !============================================================
    !³õÊ¼»¯ÐéÄâLSº¯Êý
    DO K=1,NZ
      DO I=1,NX
        DO J=1,NY
          IF (f(I,J,K) .ge. 0.5d0) THEN
            ls_v1(I,J,K) = (NX*dx+NY*dy+NZ*dZ)
          ELSE
            ls_v1(I,J,K) = -(NX*dx+NY*dy+NZ*NZ)
          ENDIF
        ENDDO
      ENDDO
    ENDDO
    !Boundary
    Call Phi_bc%setbcs(ls_v1)
    !============================================================
    !¸ù¾ÝÐéÄâLSº¯ÊýÉèÖÃÍø¸ñ½çÃæ×´Ì¬
    do k=1,nz
      DO I=1,NX
        DO J=1,NY
          IF ( ls_v1(I,J,k)*ls_v1(I-1,J,k)   .le. 0.0d0 .OR. &
              &         ls_v1(I,J,k)*ls_v1(I+1,J,k)   .le. 0.0d0 .OR. &
              &         ls_v1(I,J,k)*ls_v1(I,J-1,k)   .le. 0.0d0 .OR. &
              &         ls_v1(I,J,k)*ls_v1(I,J+1,k)   .le. 0.0d0 .OR. &
              &         ls_v1(I,J,k)*ls_v1(I-1,J-1,k) .le. 0.0d0 .OR. &
              &         ls_v1(I,J,k)*ls_v1(I-1,J+1,k) .le. 0.0d0 .OR. &
              &         ls_v1(I,J,k)*ls_v1(I+1,J-1,k) .le. 0.0d0 .OR. &
              &         ls_v1(I,J,k)*ls_v1(I+1,J+1,k) .le. 0.0d0 .OR. &
              &         ls_v1(I,J,k)*ls_v1(I-1,J,k+1)   .le. 0.0d0 .OR. &
              &         ls_v1(I,J,k)*ls_v1(I+1,J,k+1)   .le. 0.0d0 .OR. &
              &         ls_v1(I,J,k)*ls_v1(I,J-1,k+1)   .le. 0.0d0 .OR. &
              &         ls_v1(I,J,k)*ls_v1(I,J+1,k+1)   .le. 0.0d0 .OR. &
              &         ls_v1(I,J,k)*ls_v1(I-1,J-1,k+1) .le. 0.0d0 .OR. &
              &         ls_v1(I,J,k)*ls_v1(I-1,J+1,k+1) .le. 0.0d0 .OR. &
              &         ls_v1(I,J,k)*ls_v1(I+1,J-1,k+1) .le. 0.0d0 .OR. &
              &         ls_v1(I,J,k)*ls_v1(I+1,J+1,k+1) .le. 0.0d0 .OR. &
              &         ls_v1(I,J,k)*ls_v1(I-1,J,k-1)   .le. 0.0d0 .OR. &
              &         ls_v1(I,J,k)*ls_v1(I+1,J,k-1)   .le. 0.0d0 .OR. &
              &         ls_v1(I,J,k)*ls_v1(I,J-1,k-1)   .le. 0.0d0 .OR. &
              &         ls_v1(I,J,k)*ls_v1(I,J+1,k-1)   .le. 0.0d0 .OR. &
              &         ls_v1(I,J,k)*ls_v1(I-1,J-1,k-1) .le. 0.0d0 .OR. &
              &         ls_v1(I,J,k)*ls_v1(I-1,J+1,k-1) .le. 0.0d0 .OR. &
              &         ls_v1(I,J,k)*ls_v1(I+1,J-1,k-1) .le. 0.0d0 .OR. &
              &         ls_v1(I,J,k)*ls_v1(I+1,J+1,k-1) .le. 0.0d0 .OR. &
              &         ls_v1(I,J,k)*ls_v1(I,J,k+1) .le. 0.0d0 .OR. &
              &         ls_v1(I,J,k)*ls_v1(I,J,k-1) .le. 0.0d0  ) THEN
            s(I,J,k) = 1
          ELSE
            s(I,J,k) = 0
          ENDIF
        ENDDO
      ENDDO
    enddo
    !============================================================
    !¸üÐÂ½çÃæÍø¸ñ×´Ì¬Îª"1"(½çÃæÍø¸ñ)µÄÐéÄâLSº¯Êý
    do k=1,nz
      DO I=1,NX
        DO J=1,NY

          IF ( s(I,J,k) .eq. 1 ) THEN
            ls_v1(I,J,k) = 0.0d0
            hs = Heaviside(ls_v1(I,J,k),dsqrt(2.0d0)*dx)
            DO WHILE (dabs(hs-f(I,J,k)) .gt. 1.0d-15)
              hs = Heaviside(ls_v1(I,J,k),dsqrt(2.0d0)*dx)
              dr = Dirac(ls_v1(I,J,k),dsqrt(2.0d0)*dx)
              ls_v1(I,J,k) = ls_v1(I,J,k)-(hs-f(I,J,k))/(dr+1.0d-15)
            ENDDO
          ENDIF
        ENDDO
      ENDDO
    enddo


    Call Phi_bc%setbcs(ls_v1)
    !============================================================
    !¸üÐÂ½çÃæÍø¸ñ×´Ì¬Îª"0"(·Ç½çÃæÍø¸ñ)µÄÐéÄâLSº¯Êý
    h = dx
    !-------------------------------------
    !X+,Y+·½ÏòÉ¨Ãè
    DO l=1,8

      if (l.eq.1) then
        !X+,Y+,Z+·½ÏòÉ¨Ãè
        Istart = 1
        Istep  = 1
        Iend   = nx

        Jstart = 1
        Jstep  = 1
        Jend   = ny

        Kstart = 1
        Kstep  = 1
        Kend   = nz
      elseif (l.eq.2) then
        !X-,Y+,Z+·½ÏòÉ¨Ãè
        Istart = nx
        Istep  = -1
        Iend   = 1

        Jstart = 1
        Jstep  = 1
        Jend   = ny

        Kstart = 1
        Kstep  = 1
        Kend   = nz
      elseif (l.eq.3) then
        !X+,Y-,Z+·½ÏòÉ¨Ãè
        Istart = 1
        Istep  = 1
        Iend   = nx

        Jstart = ny
        Jstep  = -1
        Jend   = 1

        Kstart = 1
        Kstep  = 1
        Kend   = nz
      elseif  (l.eq.4) then
        !X-,Y-,Z+·½ÏòÉ¨Ãè
        Istart = nx
        Istep  = -1
        Iend   = 1

        Jstart = ny
        Jstep  = -1
        Jend   = 1

        Kstart = 1
        Kstep  = 1
        Kend   = nz
      ELSEif (l.eq.5) then
        !X+,Y+,Z-·½ÏòÉ¨Ãè
        Istart = 1
        Istep  = 1
        Iend   = nx

        Jstart = 1
        Jstep  = 1
        Jend   = ny

        Kstart = nz
        Kstep  = -1
        Kend   = 1
      elseif (l.eq.6) then
        !X-,Y+,Z-·½ÏòÉ¨Ãè
        Istart = nx
        Istep  = -1
        Iend   = 1

        Jstart = 1
        Jstep  = 1
        Jend   = ny

        Kstart = nz
        Kstep  = -1
        Kend   = 1
      elseif (l.eq.7) then
        !X+,Y-,Z-·½ÏòÉ¨Ãè
        Istart = 1
        Istep  = 1
        Iend   = nx

        Jstart = ny
        Jstep  = -1
        Jend   = 1

        Kstart = nz
        Kstep  = -1
        Kend   = 1
      else 
        !X-,Y-,Z-·½ÏòÉ¨Ãè
        Istart = nx
        Istep  = -1
        Iend   = 1

        Jstart = ny
        Jstep  = -1
        Jend   = 1

        Kstart = nz
        Kstep  = -1
        Kend   = 1
      endif

      DO K=KSTART,KEND,KSTEP
        DO J=Jstart,Jend,Jstep
          DO I=Istart,Iend,Istep

            !ls>=0
            IF ( ls_v1(I,J,K) .ge. 0.0d0 ) THEN
              a = dmin1(ls_v1(I-1,J,K),ls_v1(I+1,J,K))
              b = dmin1(ls_v1(I,J-1,K),ls_v1(I,J+1,K))
              c = dmin1(ls_v1(I,J,K-1),LS_V1(I,J,K+1))
              call rank3(a,b,c)

              if ( (a+h).lt.b ) then
                dnew = a + h 
              else 
                dnew = (a+b+dsqrt(2.0*h*h-(a-b)**2))/2.0d0
                IF ( dnew .gt. c ) THEN
                  dnew = (a+b+c+dsqrt(3.0*h*h-(a-b)**2-(a-c)**2-(b-c)**2))/3.0d0

                endif
              endif

              IF ( s(I,J,k) .eq. 0 ) THEN
                ls_v1(I,J,k) = dmin1(ls_v1(I,J,k),dnew)
              ENDIF
              !ls<0
            ELSE
              a = dmax1(ls_v1(I-1,J,k),ls_v1(I+1,J,k))
              b = dmax1(ls_v1(I,J-1,k),ls_v1(I,J+1,k))
              c = dmax1(ls_v1(I,J,k-1),ls_v1(I,J,k+1))
              call rank3(a,b,c)

              if( (c-h).gt.b ) then
                dnew = c - h
              else
                dnew = (b+c-dsqrt(2.0*h*h-(b-c)**2))/2.0d0
                IF ( dnew .lt. a ) THEN
                  dnew = (a+b+c-dsqrt(3.0*h*h-(a-b)**2-(a-c)**2-(b-c)**2))/3.0d0
                endif
              endif

              IF ( s(I,J,k) .eq. 0 ) THEN
                ls_v1(I,J,k) = dmax1(ls_v1(I,J,k),dnew)
              ENDIF
            ENDIF
          ENDDO
        ENDDO
      ENDDO
      Call Phi_bc%setbcs(ls_v1)
    ENDDO




    !============================================================
    !¸üÐÂLSº¯Êý
    ls(1:NX,1:NY,1:nz)=LS_V1(1:NX,1:NY,1:nz)
    !============================================================
    RETURN

  END subroutine GET_LS_3D_INIT

  !==========================
  ! (5-4) permulation function for FSM
  !==========================

  subroutine rank3(a,b,c)

    implicit none

    real(8)   ::  a,b,c
    real(8)   ::  d

    if(a.gt.b) then
      d=a
      a=b
      b=d
    endif

    if(a.gt.c) then
      d=a
      a=c
      c=d
    endif

    if(b.gt.c) then
      d=b
      b=c
      c=d
    endif

  end subroutine rank3

End Module ModVOFFuncExt