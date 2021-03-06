!=================
! Includes:
!     (1) Norm Calculation
!       (1-1) Norm calculation Parker and Youngs
!       (1-2) Norm calculation Central difference
!       (1-3) Norm calculation Mixed Central difference and Youngs (MYC)
!       (1-3) Norm calculation ELVIRA
!     (2) Flood algorithm (VOF Reconstruction)
!       (2-1) Flood algorithm backward finding alpha
!       (2-2) Flood algorithm backward finding centroid
!       (2-3) Flood algorithm forward finding vof
!       (2-4) Flood algorithm forward finding vof and centroid
!       (2-5) Flood algorithm backward of THINC/SW finding jump center
!       (2-6) Flood algorithm forward of THINC/SW finding vof
!     (3) MOF Reconstruction
!       (3-1) MOF reconstruction using quick centroid and Gauss-Newton iteration
!       (3-2) MOF reconstruction using Lemoine's analytic gradient and Gauss-Newton iteration
!       (3-3) MOF reconstruction using Lemoine's analytic gradient and BFGS iteration
!       (3-4) MOF reconstruction using Sussman's computational geometry and Gauss-Newton iteration
!       (3-5) Find centroid for given norm and volume function
!     (4) Misc
!       (4-1) Normalization vector 1 (nx+ny+nz=1)
!       (4-2) Normalization vector 2 (nx**2+ny**2+nz**2=1)
!       (4-3) Cartesian norm to Spherical angle
!       (4-4) Spherical angle to Cartesian norm
!       (4-5) Anvance of angle
!     (5) Level set related
!       (5-1) Smooth vof function
!       (5-2) VOF to level set function near interface
!       (5-3) Fast sweeping method (FSM) for ls function for far field
!       (5-4) permulation function for FSM
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
!   A pMixed Central difference and Youngs (MYCorocefure pointer MOFNorm is used to determine which MOF to choose,
!   by default, it is pointed to MOFZY. To change the target function, for example
!   to MOFSussmanGaussNewton, use the following sentence in code
!        MOFNorm => MOFSussmanGaussNewton
!=================
Module ModVOFFunc
  Use ModGlobal, only : sp
  use mod_cg3_polyhedron
  use mod_cg3_complete_polyhedron_structure
  Implicit None
  public

  ! MOF iteration parameters
  Integer, Parameter :: MOFITERMAX = 10
  Real(sp), Parameter :: tol = 1.0e-8
  Real(sp), Parameter :: local_tol = 1.0e-8
  ! Real(sp), Parameter :: GaussNewtonTol = 1e-13
  Real(sp), Parameter :: MOF_Pi = 3.1415926535897932d0
  Real(sp), Parameter :: epsc = 1.0e-12
  Real(sp) :: mof_niter(2)

  type(t_polyhedron) :: LemoinePoly

  PROCEDURE(InterfaceMOF), POINTER :: MOFNorm => MOFZY
  PROCEDURE(InterfaceVOF), POINTER :: VOFNorm => NormMYCS

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
  ! (1-2) Normal vector Central Scheme
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
  ! (1-3) Mixed Youngs and Central difference
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

    ! ! choose between the previous choice and Youngs-CIAM 
    if (DABS(m(cn,cn)) > maxval(abs(mm)))  cn = 3

    ! components of the normal vector */
    norm(1:3) = m(cn,0:2)

    return
  End Subroutine NormMYCS

  !=======================================================
  ! (1-3) Efficient Least-square volume averaging
  ! * returns normal normalized so that |mx|+|my|+|mz| = 1*
  ! In the ELVIRA reference paper, 5*5*5 stencil is used,
  ! However, here 3*3*3 stencil is use.
  ! For each direction, every 24 posible triplets that contains
  ! The central cell are evaluated.
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
    Integer :: i, j, k, dir, nelvira

    p1is =(/1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,2,2,2,2,2,2,3,3,3/)
    p1js =(/1,1,1,1,1,1,2,2,2,2,2,3,3,3,3,1,1,1,3,3,3,1,1,2/)
    p3is =(/1,1,2,2,3,3,1,2,2,3,3,2,2,3,3,3,3,3,3,3,3,3,3,3/)
    p3js =(/2,3,1,3,2,3,3,1,3,1,3,1,3,2,3,1,2,3,1,2,3,2,3,3/)

    block_vol = sum(c) / 27.0_sp

    norm = 1.0/3.0_sp
    err_min = 100.0_sp

    Do Dir = 1, 3
      If (dir == 1) Then
        Do j = 1,3
          Do i = 1,3
            HVolume(i,j) = sum(c(i,j,:))
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
            HVolume(i,j) = sum(c(:,i,j))
          End Do
        End Do
      End If
      ! Find the minimum normal voector at the 24 candidates
      ! Note that the flip case of the 24 candidates should be calculated as well
      ! So, actually 48 candidates
      Do nelvira = 1, 24
        p1i = p1is(nelvira)
        p1j = p1js(nelvira)
        p3i = p3is(nelvira)
        p3j = p3js(nelvira)
        If (dir == 1) Then
          vec_1 = (/dble(p1i-p2i), dble(p1j-p2j), Hvolume(p1i,p1j)-Hvolume(p2i,p2j)/)
          vec_2 = (/dble(p2i-p2i), dble(p3j-p2j), Hvolume(p3i,p3j)-Hvolume(p2i,p2j)/)
        Else If (dir == 2) Then
          vec_1 = (/dble(p1i-p2i), Hvolume(p1i,p1j)-Hvolume(p2i,p2j), dble(p1j-p2j)/)
          vec_2 = (/dble(p2i-p2i), Hvolume(p3i,p3j)-Hvolume(p2i,p2j), dble(p3j-p2j)/)
        Else
          vec_1 = (/Hvolume(p1i,p1j)-Hvolume(p2i,p2j), dble(p1i-p2i), dble(p1j-p2j)/)
          vec_2 = (/Hvolume(p3i,p3j)-Hvolume(p2i,p2j), dble(p2i-p2i), dble(p3j-p2j)/)
        End If
        norm_e(1) = vec_1(2) * vec_2(3) - vec_1(3) * vec_2(2)
        norm_e(2) = vec_1(3) * vec_2(1) - vec_1(1) * vec_2(3)
        norm_e(3) = vec_1(1) * vec_2(2) - vec_1(2) * vec_2(1)
        Call Normalization1(norm_e)

        alpha =  3.0 * FloodSZ_backward(norm_e,block_vol)
        dxl = 1.0_sp
        Do k = 1, 3
          Do j = 1, 3
            Do i = 1, 3
              x0(1) = dble(i-1)
              x0(2) = dble(j-1)
              x0(3) = dble(k-1)
              vol_e(i,j,k) = FloodSZ_Forward(norm_e,alpha,x0,dxl)
            End Do
          End Do
        End Do
        err = sum( ((1.0_sp - vol_e) - c) * ((1.0_sp - vol_e) - c) )
        If ( err < err_min) Then
          norm = norm_e
          err_min = err
        End If
        err = sum( (vol_e-c) * (vol_e-c) )
        If ( err < err_min) Then
          norm = norm_e
          err_min = err
        End If
        norm_e(3) = - norm_e(3)
        Call Normalization1(norm_e)

        alpha =  3.0 * FloodSZ_backward(norm_e,block_vol)
        dxl = 1.0_sp
        Do k = 1, 3
          Do j = 1, 3
            Do i = 1, 3
              x0(1) = dble(i-1)
              x0(2) = dble(j-1)
              x0(3) = dble(k-1)
              vol_e(i,j,k) = FloodSZ_Forward(norm_e,alpha,x0,dxl)
            End Do
          End Do
        End Do
        err = sum( ((1.0_sp - vol_e) - c) * ((1.0_sp - vol_e) - c) )
        If ( err < err_min) Then
          norm = norm_e
          err_min = err
        End If
        err = sum( (vol_e-c) * (vol_e-c) )
        If ( err < err_min) Then
          norm = norm_e
          err_min = err
        End If

      End Do
    End Do

    return
  End Subroutine NormELVIRA

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
  ! (2-5) Backward flooding algorithm of THINC/SW
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
    a2 = a0 - 1.0_sp
    THINC1DBackward = 1.0d0 / betagamma2 * dlog(abs(a1 / a2))

  End Function THINC1DBackward

  !=======================================================
  ! (2-6) Forward flooding algorithm of THINC/SW
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
    a0 = dexp(betagamma2 * ( x_center - x0 ) ) + 1.0d0

    THINC1DForward = deltax + 1.0_sp / betagamma2 * dlog( a1 / a0 )

  End Function THINC1DForward


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

    Real(sp)  :: c_mof(3)
    Real(sp), Dimension(3)   :: norm_2
    Real(sp) :: delta_theta
    Real(sp) :: delta_theta_max
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
      Norm_2(1) = - c_mof(1)
      Norm_2(2) = - c_mof(2)
      Norm_2(3) = - c_mof(3)
    EndIf

    Call Normalization2(Norm_2)
    Call Norm2Angle(angle_init,norm_2)
    ! Initialize other data
    do dir=1,2
      angle_array(dir,1)= angle_init(dir)
    enddo
    Call FindCentroid(angle_init,f,cen_init)
    do dir=1,3
      d_array(dir,1) = c_mof(dir) - cen_init(dir)
      cen_array(dir,1)=cen_init(dir)
      err_array(1)=err_array(1) + d_array(dir,1) * d_array(dir,1)
    enddo
    err_array(1) = sqrt(dot_product(d_array(:,1),d_array(:,1)))
    delta_theta = MOF_Pi / 1800.0_sp  ! 1 degree=pi/180
    delta_theta_max = 10.0_sp * MOF_Pi / 180.0_sp  ! 10 degrees

    iter = 0
    err = err_array(1)
    err_local_min = err_array(1)
    mof_niter = 0
    Do While ((iter.lt.MOFITERMAX).and. &
        (err.gt.tol).and. &
        (err_local_min.gt.local_tol))

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
        err_plus(i_angle)=sqrt(err_plus(i_angle))
        Call FindCentroid(angle_minus,f,cenm)
        err_minus(i_angle)= 0.0_sp
        do dir=1,3
          dm(dir) = c_mof(dir) - cenm(dir)
          err_minus(i_angle)=err_minus(i_angle)+dm(dir)**2
          d_minus(dir,i_angle)=dm(dir)
          cen_minus(dir,i_angle)=cenm(dir)
        enddo
        err_minus(i_angle)=sqrt(err_minus(i_angle))

        ! jacobian matrix (numerical partial gradient):
        do dir=1,3
          Jacobian(dir,i_angle)=(dp(dir)-dm(dir))/(2.0_sp*delta_theta)
        enddo

      End Do ! end theta phi shift angle

      ! JT * X
      Do i_angle=1,2
        gradient(i_angle) = 2.0_sp * dot_product(dbase, Jacobian(:,i_angle))
      EndDo


      !! Begin Gauss Newton
      Singular_flag = 0
      Do i=1,2
        Do j=1,2
          Hessian(i,j) = 2.0_sp * dot_product(Jacobian(:,i), Jacobian(:,j))
        EndDo
      EndDo

      det = Hessian(1,1)*Hessian(2,2) - Hessian(1,2)*Hessian(2,1)
      ! Calculate the inverse of the matrix
      HessianT(1,1) = +Hessian(2,2) / det
      HessianT(2,1) = -Hessian(2,1) / det
      HessianT(1,2) = -Hessian(1,2) / det
      HessianT(2,2) = +Hessian(1,1) / det

      ! Call matinv2(Hessian, HessianT, det)
      If (det .lt. 1.0e-20) Then
        Singular_flag = 1
      End If

      ! Delta angle
      Do i=1,2
        delangle(i) = - dot_product(HessianT(:,i), gradient)
      End Do

      !! End Gauss Newton
      ! call GaussNewton(Jacobian, gradient, 3, 2, delangle, Singular_flag)

      ! Find delta angle
      if (singular_flag.eq.0) then
        ! -pi<angle<pi
        do i_angle=1,2
          angle_previous(i_angle)=angle_base(i_angle)
          call advance_angle(angle_base(i_angle),delangle(i_angle))
        enddo
        Call FindCentroid(angle_base,f, cenopt)
        do dir=1,3
          dopt(dir) = c_mof(dir) - cenopt(dir)
        enddo

      else if (singular_flag.eq.1) then
        err_local_min= 0.0_sp

        do dir=1,3
          dopt(dir)=dbase(dir)
        enddo
        do i_angle=1,2
          angle_base(i_angle)=angle_array(i_angle,iter+1)
        enddo
      else
        print *,"singular_flag invalid"
        stop
      End If

      err= 0.0_sp
      do dir=1,3
        err=err+dopt(dir)**2
      enddo
      err=sqrt(err)

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

    mof_niter(1) = iter
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
  End Subroutine MOFZY

  !=======================================================
  ! (3-2) MOF reconstruction using Analytic gradient and Gauss-Newton iteration
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
    Real(sp) :: delta_theta_max
    Real(sp), Dimension(2)   :: delangle, angle_init, new_angle
    Real(sp), Dimension(2)   :: angle_base
    Real(sp), Dimension(2,MOFITERMAX+1) :: angle_array
    Real(sp), Dimension(MOFITERMAX+1)   :: err_array

    Real(sp) :: err
    Integer :: singular_flag
    Integer :: i_angle, dir, iter, i, j
    Real(sp) :: dxs(3)
    Real(sp) :: c_diff(3)
    Real(sp) :: gradient(2)
    Real(sp) :: Jacobian(3,2)
    Real(sp) :: Hessian(2,2)
    Real(sp) :: HessianT(2,2)
    Real(sp) :: det

    dxs = 1.0_sp
    delta_theta_max = 3.0_sp * MOF_Pi / 180.0_sp  ! 10 degrees

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
    Else
      Norm_2(1) = 0.5_sp - c_mof(1)
      Norm_2(2) = 0.5_sp - c_mof(2)
      Norm_2(3) = 0.5_sp - c_mof(3)
    EndIf

    Call Normalization2(Norm_2)
    Call Norm2Angle(angle_init,norm_2)
    do dir=1,2
      angle_array(dir,1)= angle_init(dir)
    enddo
    Call mof3d_compute_analytic_gradient_GN(angle_init, c_mof, vof, dxs, c_diff, Jacobian)
    err = dot_product(c_diff, c_diff)
    do dir=1,3
      err_array(1) = err
    enddo

    iter = 0
    err = err_array(1)
    ! print *, err
    ! print *, iter, err, err_local_min
    mof_niter = 0
    Do While ((iter.lt.MOFITERMAX) .and. (err.gt.tol))

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

      If (det .lt. 1.0e-30) Then
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

      do i_angle=1,2
        call advance_angle(angle_base(i_angle),delangle(i_angle))
        angle_array(i_angle,iter+2)=angle_base(i_angle)
      enddo

      call Angle2Norm(angle_base,norm)
      Call Normalization1(norm)
      ! print *, norm
      err_array(iter+2)=err
      iter=iter+1
    End Do
    mof_niter(1) = iter

    do dir=1,2
      new_angle(dir)=angle_array(dir,iter+1)
    enddo
    call Angle2Norm(new_angle,norm)
    Call Normalization1(norm)

    if (f .ge. 0.5_sp) norm = -norm

  End Subroutine MOFLemoine_GaussNewton


  !=======================================================
  ! (3-3) MOF reconstruction using analytic gradient and BFGS 
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
    Real(sp) :: c_mof(3)
    Real(sp), Dimension(3)   :: norm_2
    Real(sp) :: delta_theta_max
    Real(sp), Dimension(2)   :: angle_init

    Real(sp) :: dxs(3)
    Integer :: nstat(2)
    Real(sp) :: residual(2)

    dxs = 1.0_sp
    delta_theta_max = 10.0_sp * MOF_Pi / 180.0_sp  ! 10 degrees

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
    Else
      Norm_2(1) = 0.5_sp - c_mof(1)
      Norm_2(2) = 0.5_sp - c_mof(2)
      Norm_2(3) = 0.5_sp - c_mof(3)
    EndIf

    ! print *, c_mof
    ! print *, vof
    Call Normalization2(Norm_2)
    ! norm_2 = - norm_2
    Call Norm2Angle(angle_init,norm_2)
    ! Call direction_to_spherical_angles(norm_2, angle_init)
    ! print *, norm_2

    ! volume = vol
    ! ref_centroid = ref_c
    ! angles = ang
    ! print *, c, f, angle_init
    call mof3d_bfgs(LemoinePoly, c_mof, c_mof, vof, angle_init, norm, nstat, residual)

    ! print *, nstat
    ! print *, residual

    mof_niter = nstat


    Call Normalization1(norm)
    ! call Angle2Norm(new_angle,norm)
    ! print *, norm

    if (f .ge. 0.5_sp) norm = -norm

  End Subroutine MOFLemoine_BFGS

  !=======================================================
  ! (3-4) MOF reconstruction using numerical gradient and BFGS 
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
    Real(8) :: intercept
    Integer :: default_one = 1

    c_mof = centroid - 0.5_sp

    If(present(Init_Norm))then
      norm_2 = Init_Norm
    Else
        Norm_2(1) = c_mof(1)
        Norm_2(2) = c_mof(2)
        Norm_2(3) = c_mof(3)
    EndIf

    norm_2 = norm_2 / norm2(norm_2)

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
    mof_niter(1) = nn_iter

  End Subroutine MOFSussmanGaussNewton

  !=======================================================
  ! (3-5) MOF reconstruction using Machine Learning
  !     For given vof and centroid, find the best norm
  ! Key steps:
  !   - Convert to local region
  !   - Calculate the gradient using machine learning
  !   - Flip the normal vector
  !-------------------------------------------------------
  ! Input:
  !      f: vof function
  !      c: centroid (cx, cy, cz)
  ! Optional Input:
  !      init_norm: Initial guess of the normal vector
  ! Output:
  !      norm: vof function
  !=======================================================
  ! Subroutine MOFNN(f, c, norm, Init_Norm)
  !   Implicit None
  !   Real(sp) :: flip(3)
  !   Real(sp) :: symmetric
  !   Real(sp) :: c_mof(3)

  !   c = c - 0.5_sp
  !   If ( f .ge. 0.5_sp ) Then
  !     symmetric = - 1.0_sp
  !   Else
  !     symmetric = 1.0_sp
  !   endif
  !   Do i = 1,3
  !     If ( c(i) < 0.0_sp ) Then
  !       flip(i) = .false.
  !       c_mof(i) = c(i)
  !     Else
  !       flip(i) = .true.
  !       c_mof(i) = - c(i)
  !     End Do
  !   End Do

  !   Do i = 1,3
  !     Norm(i) = Norm(i) * flip(i) * symmetric
  !   End Do

  ! End Subroutine MOFNN


  !=======================================================
  ! (3-5) Find the centroid for MOF iteration
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
  ! (3-5) Find the centroid for MOF iteration
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
    If( aa .gt. 1.d-10) Then
      norm = norm / aa
    End If
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
    If( aa .gt. 1.d-10) Then
      norm(1) = norm(1) / aa
      norm(2) = norm(2) / aa
      norm(3) = norm(3) / aa
    End If
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
      Else If (norm(2) .gt. epsc) Then
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

  subroutine direction_to_spherical_angles(direction, angles)
    double precision, dimension(3), intent(in) :: direction
    double precision, dimension(2), intent(out) :: angles

    double precision, parameter :: PI = 2d0*acos(0d0)

    if ((abs(direction(3)) - 1d0) < epsilon(1d0)) then
      angles = [atan2(direction(2), direction(1)), acos(direction(3))]
    else
      if (direction(3) > 0d0) then
        angles = [0d0, 0d0]
      else
        angles = [0d0, PI]
      end if
    end if
  end subroutine direction_to_spherical_angles

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
    REAL(sp) angle,delangle
    angle=angle+delangle
    if (angle.lt.-MOF_Pi) then
      angle=angle+2.0_sp*MOF_Pi
    endif
    if (angle.gt.MOF_Pi) then
      angle=angle-2.0_sp*MOF_Pi
    end if
  End Subroutine advance_angle

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

End Module ModVOFFunc

