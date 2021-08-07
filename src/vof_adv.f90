!=================
! Includes:
!     (1) VOF/MOF Advection
!       (1-1) VOF-WY (Explicit Implicit with mass correction)
!       [1-2] VOF-LE (Lagrangian Explicit)
!       [1-3] VOF-EI (Eulerian Implicit)
!       [1-4-1] VOF-EILE3DS (Eulerian Explicit Lanrangian Ixplicit simplified)
!       [1-4-2] VOF-EILE3D (Eulerian Explicit Lanrangian Ixplicit)
!       [1-4-3] VOF-EIEALE (Mixed Eulerian Lagrangian with Algebraic step)
!       [1-4-4] VOF-EILE2D (Mixed Eulerian Lagrangian in 2D)
!       (1-5) MOF-WY (Explicit Implicit with mass correction)
!       [1-6] MOF-LE (Lanrangian Explicit)
!       [1-7] MOF-EI (Eulerian Impllicit)
!       [1-8-1] MOF-EILE3DS (Eulerian Explicit Lanrangian Ixplicit)
!       [1-8-2] MOF-EILE3D (Eulerian Explicit Lanrangian Ixplicit)
!       [1-8-3] MOF-EIEALE (Mixed Eulerian Lagrangian with Algebraic step)
!       [1-8-4] MOF-EILE2D (Mixed Eulerian Lagrangian in 2D)
!       [1-9]  THINC-EI (Eulerian Implicit)
!       [1-10] THINC-LE (Lagrangian Explicit)
!       [1-11] THINC-WY (Explicit Implicit with mass correction)
!       [1-12-1] THINC-EILE2D (Mixed Eulerian Lagrangian in 2D)
!     (2) Directional split VOF/MOF advection
!       (2-1)  VOF-LE (Lagrangian Explicit)
!       (2-2)  VOF-EI (Eulerian Implicit)
!       (2-3)  VOF-WY (Weymouth-Yue)
!       [2-4]  VOF-EA (Eulerian Algebaric)
!       (2-5)  MOF-LE (Lanrangian Explicit)
!       (2-6)  MOF-EI (Eulerian Implicit)
!       (2-7)  MOF-WY (Eulerian Implicit)
!       [2-8]  MOF-EA (Eularian Algebaric)
!       [2-9]  THINC-EI (Eularian Implicit)
!       [2-10] THINC-LE (Lagrangian Explicit)
!       [2-10] THINC-WY (Weymouth-Yue)
!     (3) Misc
!       [3-1] Velocity decomposition for EILE3D
!       [3-2] Centroid advection (Eulerian)
!       [3-3] Centroid advection (Lagrangian)
!
!  Module procedure:
!    The vof subroutine can be called with VOFAdvection,
!        if centroid appears, MOF_EILE3DS is used
!        if centroid not appear, VOF_EILE3DS is used
!-----------------
! Author: Zhouteng Ye (yzt9zju@gmail.com)
!=================
Module ModVOF
  Use ModGlobal, only: sp
  Use ModGlobal, only: updthalo
  Use ModGlobal, only: setBCS, phi_bc
#if defined(EXTENDS)
  Use ModVOFFuncExt
#else
  Use ModVOFFunc
#endif
#if defined(VISUALIZE)
  Use ModTools
#endif

  Interface VOFAdvection
    Module Procedure VOF_WY
    Module Procedure MOF_WY
  End Interface VOFAdvection

  PROCEDURE(InterfaceForwardC), POINTER :: ForwardC => FloodSZ_ForwardC

  Interface
    Subroutine InterfaceForwardC(nr,alpha,x0,dx,f,xc0)
      import
      Implicit None
      REAL(8), INTENT(IN) :: nr(3),x0(3),dx(3),alpha
      REAL(8), INTENT(OUT) :: f
      REAL(8), INTENT(OUT) :: xc0(3)
    End Subroutine InterfaceForwardC
  End Interface

  real(sp) :: epsc_vof = 1.0d-14

  Integer :: num_iter(2) = 0
  Integer :: grid_iter = 0

  Real(sp) :: t_rec
  Real(sp) :: t_adv
  Real(sp) :: t_cor
  Contains

  !=======================================================
  ! (1-1) VOF-WY 
  ! VOF advection with Weymouth-Yue's strategy
  !-------------------------------------
  ! Sweep sequence:
  !   1-2-3
  !   2-3-1
  !   3-1-2
  ! All: WY-WY-WY
  !=======================================================
  Subroutine VOF_WY(Phi, u, v, w, nl, dl, dt, rank, phi0)
    Implicit None
    Real(sp) :: dt
    Integer, Intent(In)            :: nl(3)
    Real(sp), Intent(In)           :: dl(3)
    Real(sp), Intent(InOut)        :: Phi(0:nl(1) + 1, 0:nl(2) + 1, 0:nl(3) + 1)
    Real(sp), Intent(In), optional :: Phi0(0:nl(1) + 1, 0:nl(2) + 1, 0:nl(3) + 1)
    Real(sp), Intent(In)           :: u(0:nl(1) + 1, 0:nl(2) + 1, 0:nl(3) + 1)
    Real(sp), Intent(In)           :: v(0:nl(1) + 1, 0:nl(2) + 1, 0:nl(3) + 1)
    Real(sp), Intent(In)           :: w(0:nl(1) + 1, 0:nl(2) + 1, 0:nl(3) + 1)
    Integer :: rank
    num_iter = 0
    grid_iter = 0
    t_rec = 0.0_sp
    t_adv = 0.0_sp
    t_cor = 0.0_sp
    if (rank == 1 .or. rank == 4) Then
      call VOFAdv_WY(u, phi0, phi, nl, dl, dt, 1)
      call VOFAdv_WY(v, phi0, phi, nl, dl, dt, 2)
      call VOFAdv_WY(w, phi0, phi, nl, dl, dt, 3)
    else if (rank == 2 .or. rank == 5) Then
      call VOFAdv_WY(v, phi0, phi, nl, dl, dt, 2)
      call VOFAdv_WY(w, phi0, phi, nl, dl, dt, 3)
      call VOFAdv_WY(u, phi0, phi, nl, dl, dt, 1)
    else if (rank == 0 .or. rank == 3) Then
      call VOFAdv_WY(w, phi0, phi, nl, dl, dt, 3)
      call VOFAdv_WY(u, phi0, phi, nl, dl, dt, 1)
      call VOFAdv_WY(v, phi0, phi, nl, dl, dt, 2)
    endif
  End Subroutine VOF_WY

  !=======================================================
  ! (1-5) VOF-WY 
  ! MOF advection with Weymouth-Yue's strategy
  !-------------------------------------
  ! Sweep sequence:
  !   1-2-3
  !   2-3-1
  !   3-1-2
  ! All: WY-WY-WY
  !=======================================================
  Subroutine MOF_WY(Phi, u, v, w, nl, dl, dt, rank, cx, cy, cz, phi0)
    Implicit None
    Real(sp) :: dt
    Integer, Intent(In)            :: nl(3)
    Real(sp), Intent(In)           :: dl(3)
    Real(sp), Intent(InOut)        :: Phi(0:nl(1) + 1, 0:nl(2) + 1, 0:nl(3) + 1)
    Real(sp), Intent(In), Optional :: Phi0(0:nl(1) + 1, 0:nl(2) + 1, 0:nl(3) + 1)
    Real(sp), Intent(InOut)        :: cx(0:nl(1) + 1, 0:nl(2) + 1, 0:nl(3) + 1)
    Real(sp), Intent(InOut)        :: cy(0:nl(1) + 1, 0:nl(2) + 1, 0:nl(3) + 1)
    Real(sp), Intent(InOut)        :: cz(0:nl(1) + 1, 0:nl(2) + 1, 0:nl(3) + 1)
    Real(sp), Intent(In)           :: u(0:nl(1) + 1, 0:nl(2) + 1, 0:nl(3) + 1)
    Real(sp), Intent(In)           :: v(0:nl(1) + 1, 0:nl(2) + 1, 0:nl(3) + 1)
    Real(sp), Intent(In)           :: w(0:nl(1) + 1, 0:nl(2) + 1, 0:nl(3) + 1)
    Integer :: rank
    num_iter = 0
    grid_iter = 0
    t_rec = 0.0_sp
    t_adv = 0.0_sp
    t_cor = 0.0_sp
    if (rank == 0 .or. rank == 3) Then
      call MOFAdv_WY(u, phi0, cx, cy, cz, phi, nl, dl, dt, 1)
      call MOFAdv_WY(v, phi0, cx, cy, cz, phi, nl, dl, dt, 2)
      call MOFAdv_WY(w, phi0, cx, cy, cz, phi, nl, dl, dt, 3)
    else if (rank == 1 .or. rank == 4) Then
      call MOFAdv_WY(w, phi0, cx, cy, cz, phi, nl, dl, dt, 3)
      call MOFAdv_WY(u, phi0, cx, cy, cz, phi, nl, dl, dt, 1)
      call MOFAdv_WY(v, phi0, cx, cy, cz, phi, nl, dl, dt, 2)
    else if (rank == 2 .or. rank == 5) Then
      call MOFAdv_WY(v, phi0, cx, cy, cz, phi, nl, dl, dt, 2)
      call MOFAdv_WY(w, phi0, cx, cy, cz, phi, nl, dl, dt, 3)
      call MOFAdv_WY(u, phi0, cx, cy, cz, phi, nl, dl, dt, 1)
    end if
  End Subroutine MOF_WY

  !=======================================================
  ! (2-1) Advection algorithm of CIAM PLIC (LE)
  !=======================================================
  Subroutine VOFAdv_LE(us, f, nl, dl, dt, dir)
    ! Use ModVOFFunc
    Implicit None
    Integer :: dir
    Real(sp) :: dt
    Integer, Intent(In)     :: nl(3)
    Real(sp), Intent(In)    :: dl(3)
    Real(sp), Intent(In)    :: us(0:nl(1) + 1, 0:nl(2) + 1, 0:nl(3) + 1)
    Real(sp), Intent(InOut) :: f(0:nl(1) + 1, 0:nl(2) + 1, 0:nl(3) + 1)
    Integer  :: i, j, k, ii, jj, kk
    Real(sp) :: a1, a2, alpha
    Real(sp) :: x0(3), deltax(3)
    Real(sp), dimension(0:nl(1) + 1, 0:nl(2) + 1, 0:nl(3) + 1) :: vof1, vof2, vof3
    Real(sp), dimension(0:nl(1) + 1, 0:nl(2) + 1, 0:nl(3) + 1) :: n1, n2, n3
    Integer,  dimension(0:nl(1) + 1, 0:nl(2) + 1, 0:nl(3) + 1) :: flag
    Real(sp) :: tt1, tt2

    ii = 0; jj=0; kk=0
    If (dir .eq. 1) ii = 1
    If (dir .eq. 2) jj = 1
    If (dir .eq. 3) kk = 1

    ! Normal vector
    Call VOFNormalVectors(f, nl, flag, n1, n2, n3)

    ! CAlculate flux
    Call VOF_LagrangianFlux(us, f, nl, dl, dt, dir, &
                              n1, n2, n3, flag, &
                              vof1, vof2, vof3)

    Call CPU_TIME(tt1)  !!! cpu time for correction term
    !new values of c and  clip it: 0. <= c <= 1.
    Do k = 1, nl(3)
      Do j = 1, nl(2)
        Do i = 1, nl(1)
          f(i, j, k) = vof3(i - ii, j - jj, k - kk) + vof2(i, j, k) + vof1(i + ii, j + jj, k + kk)
          If (f(i, j, k) < EPSC_VOF) Then
            f(i, j, k) = 0.0_sp
          ElseIf (f(i, j, k) > (1.d0 - EPSC_VOF)) Then
            f(i, j, k) = 1.0_sp
          EndIf
        EndDo
      EndDo
    EndDo

    Call CPU_TIME(tt2)  !!! cpu time for correction term
    t_cor = t_cor + tt2 - tt1

    ! apply proper boundary conditions to c
    call phi_bc%setBCS(f)
    !***
    return
  End Subroutine VOFAdv_LE

  !=======================================================
  ! (2-2) Advection algorithm of Dick-Weymouth PLIC (EI)
  !=======================================================
  Subroutine VOFAdv_EI(us, f, nl, dl, dt, dir)
    Implicit None
    Integer :: dir
    Real(sp) :: dt
    Integer, Intent(In)     :: nl(3)
    Real(sp), Intent(In)    :: dl(3)
    Real(sp), Intent(In)    :: us(0:nl(1) + 1, 0:nl(2) + 1, 0:nl(3) + 1)
    Real(sp), Intent(InOut) :: f(0:nl(1) + 1, 0:nl(2) + 1, 0:nl(3) + 1)
    Integer  :: i, j, k, ii, jj, kk
    Real(sp) :: a1, a2, alpha
    Real(sp) :: x0(3), deltax(3)
    Real(sp), dimension(0:nl(1) + 1, 0:nl(2) + 1, 0:nl(3) + 1) :: vof1, vof2, vof3
    Real(sp), dimension(0:nl(1) + 1, 0:nl(2) + 1, 0:nl(3) + 1) :: n1, n2, n3
    Integer, dimension(0:nl(1) + 1, 0:nl(2) + 1, 0:nl(3) + 1) :: flag
    Real(sp) :: norm(3)
    Real(sp) :: f_block(3, 3, 3)
    Real(sp) :: tt1, tt2

    ii = 0; jj=0; kk=0
    If (dir .eq. 1) ii = 1
    If (dir .eq. 2) jj = 1
    If (dir .eq. 3) kk = 1

    Call VOFNormalVectors(f, nl, flag, n1, n2, n3)

    ! Calculate flux
    Call VOF_EulerianFlux(us, f, nl, dl, dt, dir, &
                              n1, n2, n3, flag, &
                              vof1, vof2, vof3)

    Call CPU_TIME(tt1)  !!! cpu time for correction term

    !new values of f and  clip it: 0. <= f <= 1.
    Do k = 1, nl(3)
      Do j = 1, nl(2)
        Do i = 1, nl(1)
          a1 = us(i - ii, j - jj, k - kk)*dt/dl(dir)
          a2 = us(i, j, k)*dt/dl(dir)
          f(i, j, k) = (vof3(i - ii, j - jj, k - kk) + vof2(i, j, k) + vof1(i + ii, j + jj, k + kk))/(1.0_sp - a2 + a1)
          If (f(i, j, k) < EPSC_VOF) Then
            f(i, j, k) = 0.0_sp
          ElseIf (f(i, j, k) > (1.0_sp - EPSC_VOF)) Then
            f(i, j, k) = 1.0_sp
          EndIf
        EndDo
      EndDo
    EndDo
    Call CPU_TIME(tt2)  !!! cpu time for correction term
    t_cor = t_cor + tt2 - tt1

    ! apply proper boundary conditions to volume fraction
    call phi_bc%setBCS(f)

  End Subroutine VOFAdv_EI

  !=======================================================
  ! (2-3) Advection algorithm of Weymouth-Dick (EI)
  !-------------------------------------------------------
  ! Weymouth-Yue's mass correction strategy
  !=======================================================
  Subroutine VOFAdv_WY(us, f0, f, nl, dl, dt, dir)
    Implicit None
    Integer :: dir
    Real(sp) :: dt
    Integer, Intent(In)     :: nl(3)
    Real(sp), Intent(In)    :: dl(3)
    Real(sp), Intent(In)    :: us(0:nl(1) + 1, 0:nl(2) + 1, 0:nl(3) + 1)
    Real(sp), Intent(In)    :: f0(0:nl(1) + 1, 0:nl(2) + 1, 0:nl(3) + 1)
    Real(sp), Intent(InOut) :: f(0:nl(1) + 1, 0:nl(2) + 1, 0:nl(3) + 1)
    Integer  :: i, j, k, ii, jj, kk
    Real(sp) :: a1, a2, alpha
    Real(sp) :: b1, b2
    Real(sp) :: c1, c2
    Real(sp) :: x0(3), deltax(3)
    Real(sp), dimension(0:nl(1) + 1, 0:nl(2) + 1, 0:nl(3) + 1) :: vof1, vof2, vof3
    Real(sp), dimension(0:nl(1) + 1, 0:nl(2) + 1, 0:nl(3) + 1) :: n1, n2, n3
    Integer, dimension(0:nl(1) + 1, 0:nl(2) + 1, 0:nl(3) + 1) :: flag
    Real(sp) :: norm(3)
    Real(sp) :: f_block(3, 3, 3)
    Real(sp) :: exp_factor
    Real(sp) :: exp_factor1, exp_factor2
    Real(sp) :: tt1, tt2
    ! Real(sp) :: FloodSZ_Backward, FloodSZ_Forward

    ii = 0; jj=0; kk=0
    If (dir .eq. 1) ii = 1
    If (dir .eq. 2) jj = 1
    If (dir .eq. 3) kk = 1

    Call VOFNormalVectors(f, nl,flag, n1, n2, n3)

    Call VOF_EulerianFlux(us, f, nl, dl, dt, dir, &
                              n1, n2, n3, flag, &
                              vof1, vof2, vof3)

    Call CPU_TIME(tt1)  !!! cpu time for correction term

    !new values of f and  clip it: 0. <= f <= 1.
    do k = 1, nl(3)
      do j = 1, nl(2)
        do i = 1, nl(1)
          a1 = us(i - ii, j - jj, k - kk)*dt/dl(dir)
          a2 = us(i, j, k)*dt/dl(dir)
          f(i, j, k) =  vof3(i - ii, j - jj, k - kk) + &
                        vof1(i + ii, j + jj, k + kk) &
                        +vof2(i,j,k) 
          if ( f0(i,j,k) > 0.5_sp) f(i,j,k) = f(i,j,k) + a2 - a1
          if (f(i, j, k) < EPSC_VOF) then
            f(i, j, k) = 0.0_sp
          elseif (f(i, j, k) > (1.0_sp - EPSC_VOF)) then
            f(i, j, k) = 1.0_sp
          endif
        enddo
      enddo
    enddo
    Call CPU_TIME(tt2)  !!! cpu time for correction term
    t_cor = t_cor + tt2 - tt1
    ! apply proper boundary conditions to c
    call phi_bc%setBCS(f)
    !***
  End Subroutine VOFAdv_WY

  !=======================================================
  ! (2-5) Advection algorithm of CIAM MOF (LE)
  !=======================================================
  Subroutine MOFAdv_LE(us, cx, cy, cz, f, nl, dl, dt, dir)
    Implicit None
    Integer :: dir
    Real(sp) :: dt
    Integer, Intent(In)     :: nl(3)
    Real(sp), Intent(In)    :: dl(3)
    Real(sp), Intent(In)    :: us(0:nl(1) + 1, 0:nl(2) + 1, 0:nl(3) + 1)
    Real(sp), Intent(InOut) :: cx(0:nl(1) + 1, 0:nl(2) + 1, 0:nl(3) + 1)
    Real(sp), Intent(InOut) :: cy(0:nl(1) + 1, 0:nl(2) + 1, 0:nl(3) + 1)
    Real(sp), Intent(InOut) :: cz(0:nl(1) + 1, 0:nl(2) + 1, 0:nl(3) + 1)
    Real(sp), Intent(InOut) :: f(0:nl(1) + 1, 0:nl(2) + 1, 0:nl(3) + 1)
    Integer  :: i, j, k, ii, jj, kk
    Real(sp) :: a1, a2, alpha
    Real(sp) :: mx, my, mz
    Real(sp), dimension(0:nl(1) + 1, 0:nl(2) + 1, 0:nl(3) + 1) :: vof1, vof2, vof3
    Real(sp), dimension(0:nl(1) + 1, 0:nl(2) + 1, 0:nl(3) + 1) :: c1x, c1y, c1z
    Real(sp), dimension(0:nl(1) + 1, 0:nl(2) + 1, 0:nl(3) + 1) :: c2x, c2y, c2z
    Real(sp), dimension(0:nl(1) + 1, 0:nl(2) + 1, 0:nl(3) + 1) :: c3x, c3y, c3z
    Real(sp), dimension(0:nl(1) + 1, 0:nl(2) + 1, 0:nl(3) + 1) :: n1, n2, n3
    Integer, dimension(0:nl(1) + 1, 0:nl(2) + 1, 0:nl(3) + 1) :: flag
    Real(sp) :: norm(3)
    Real(sp) :: c3(3)
    Real(sp) :: tt1,tt2
    ! Real(sp) :: FloodSZ_Backward

    ! default, f = 0
    vof1 = 0.0_sp
    vof2 = 0.0_sp
    vof3 = 0.0_sp

    ii = 0; jj=0; kk=0
    If (dir .eq. 1) ii = 1
    If (dir .eq. 2) jj = 1
    If (dir .eq. 3) kk = 1

    Call MOFNormalVectors(f, cx, cy, cz, nl, flag, n1, n2, n3)

    Call MOF_LagrangianFlux(us, f, cx, cy, cz, nl, dl, dt, dir, &
                              n1, n2, n3, flag, &
                              c1x, c1y, c1z, &
                              c2x, c2y, c2z, &
                              c3x, c3y, c3z, &
                              vof1, vof2, vof3)

    Call CPU_TIME(tt1)  !!! cpu time for correction term
    ! new values of c and  clip it: 0. <= c <= 1.
    do k = 1, nl(3)
      do j = 1, nl(2)
        do i = 1, nl(1)

          f(i, j, k) =  vof3(i - ii, j - jj, k - kk) + &
                        vof2(i, j, k) + &
                        vof1(i + ii, j + jj, k + kk)
          if (f(i, j, k) < EPSC_VOF) then
            f(i, j, k)  = 0.0_sp
            cx(i, j, k) = 0.0_sp
            cy(i, j, k) = 0.0_sp
            cz(i, j, k) = 0.0_sp
          elseif (f(i, j, k) > (1.0_sp - EPSC_VOF)) then
            f(i, j, k) = 1.0_sp
            cx(i, j, k) = 0.5_sp
            cy(i, j, k) = 0.5_sp
            cz(i, j, k) = 0.5_sp
          else
            mx =  vof3(i - ii, j - jj, k - kk)*c3x(i - ii, j - jj, k - kk) + &
                  vof2(i, j, k)*c2x(i, j, k) + &
                  vof1(i + ii, j + jj, k + kk)*c1x(i + ii, j + jj, k + kk)
            my =  vof3(i - ii, j - jj, k - kk)*c3y(i - ii, j - jj, k - kk) + &
                  vof2(i, j, k)*c2y(i, j, k) + &
                  vof1(i + ii, j + jj, k + kk)*c1y(i + ii, j + jj, k + kk)
            mz =  vof3(i - ii, j - jj, k - kk)*c3z(i - ii, j - jj, k - kk) + &
                  vof2(i, j, k)*c2z(i, j, k) + &
                  vof1(i + ii, j + jj, k + kk)*c1z(i + ii, j + jj, k + kk)
            mx = min(max(mx, 0.0_sp), 0.5_sp)
            my = min(max(my, 0.0_sp), 0.5_sp)
            mz = min(max(mz, 0.0_sp), 0.5_sp)

            cx(i, j, k) = mx / f(i, j, k)
            cy(i, j, k) = my / f(i, j, k)
            cz(i, j, k) = mz / f(i, j, k)
            cx(i, j, k) = min(max(cx(i, j, k), 0.0_sp), 1.0_sp)
            cy(i, j, k) = min(max(cy(i, j, k), 0.0_sp), 1.0_sp)
            cz(i, j, k) = min(max(cz(i, j, k), 0.0_sp), 1.0_sp)
          endif
        enddo
      enddo
    enddo
    Call CPU_TIME(tt2)
    t_cor = t_cor + tt2 - tt1

    ! apply proper boundary conditions to c
    call phi_bc%setBCS(f)
    call phi_bc%setBCS(cx)
    call phi_bc%setBCS(cy)
    call phi_bc%setBCS(cz)
    !***
    return
  End Subroutine MOFAdv_LE

  !=======================================================
  ! (2-6) Advection algorithm of CIAM MOF (LE)
  !=======================================================
  Subroutine MOFAdv_EI(us, cx, cy, cz, f, nl, dl, dt, dir)
    Implicit None
    Integer :: dir
    Real(sp) :: dt
    Integer, Intent(In)     :: nl(3)
    Real(sp), Intent(In)    :: dl(3)
    Real(sp), Intent(In)    :: us(0:nl(1) + 1, 0:nl(2) + 1, 0:nl(3) + 1)
    Real(sp), Intent(InOut) :: f(0:nl(1) + 1, 0:nl(2) + 1, 0:nl(3) + 1)
    Real(sp), Intent(InOut) :: cx(0:nl(1) + 1, 0:nl(2) + 1, 0:nl(3) + 1)
    Real(sp), Intent(InOut) :: cy(0:nl(1) + 1, 0:nl(2) + 1, 0:nl(3) + 1)
    Real(sp), Intent(InOut) :: cz(0:nl(1) + 1, 0:nl(2) + 1, 0:nl(3) + 1)
    Integer  :: i, j, k, ii, jj, kk
    Real(sp) :: a1, a2, alpha
    Real(sp) :: mx, my, mz
    Real(sp) :: x0(3), deltax(3)
    Real(sp) :: cxyz(3), c1xyz(3), c2xyz(3), c3xyz(3)
    Real(sp), dimension(0:nl(1) + 1, 0:nl(2) + 1, 0:nl(3) + 1) :: vof1, vof2, vof3
    Real(sp), dimension(0:nl(1) + 1, 0:nl(2) + 1, 0:nl(3) + 1) :: c1x, c1y, c1z
    Real(sp), dimension(0:nl(1) + 1, 0:nl(2) + 1, 0:nl(3) + 1) :: c2x, c2y, c2z
    Real(sp), dimension(0:nl(1) + 1, 0:nl(2) + 1, 0:nl(3) + 1) :: c3x, c3y, c3z
    Real(sp), dimension(0:nl(1) + 1, 0:nl(2) + 1, 0:nl(3) + 1) :: n1, n2, n3
    Integer, dimension(0:nl(1) + 1, 0:nl(2) + 1, 0:nl(3) + 1) :: flag
    Real(sp) :: norm(3)
    Real(sp) :: f_block(3, 3, 3)
    Real(sp) :: exp_factor
    Real(sp) :: tt1, tt2

    ! default, f = 0
    vof1 = 0.0_sp
    vof3 = 0.0_sp
    ii = 0; jj=0; kk=0
    If (dir .eq. 1) ii = 1
    If (dir .eq. 2) jj = 1
    If (dir .eq. 3) kk = 1

    Call MOFNormalVectors(f, cx, cy, cz, nl, flag, n1, n2, n3)
    Call MOF_EulerianFlux(us, f, cx, cy, cz, nl, dl, dt, dir, &
                              n1, n2, n3, flag, &
                              c1x, c1y, c1z, &
                              c2x, c2y, c2z, &
                              c3x, c3y, c3z, &
                              vof1, vof2, vof3)

    Call CPU_TIME(tt1)  !!! cpu time for correction term

   !new values of f and  clip it: 0. <= f <= 1.
    do k = 1, nl(3)
      do j = 1, nl(2)
        do i = 1, nl(1)
          a1 = us(i - ii, j - jj, k - kk)*dt/dl(dir)
          a2 = us(i, j, k)*dt/dl(dir)
          exp_factor = 1.0_sp/(1.0_sp - a2 + a1)
          f(i, j, k) = (vof3(i - ii, j - jj, k - kk) + vof2(i, j, k) + vof1(i + ii, j + jj, k + kk))/(1.0_sp - a2 + a1)
          if (f(i, j, k) < EPSC_VOF) then
            f(i, j, k) = 0.0_sp
            cx(i, j, k) = 0.0_sp
            cy(i, j, k) = 0.0_sp
            cz(i, j, k) = 0.0_sp
          elseif (f(i, j, k) > (1.0_sp - EPSC_VOF)) then
            f(i, j, k) = 1.0_sp
            cx(i, j, k) = 0.5_sp
            cy(i, j, k) = 0.5_sp
            cz(i, j, k) = 0.5_sp
          else
            mx = vof3(i - ii, j - jj, k - kk)*c3x(i - ii, j - jj, k - kk) + &
                 vof2(i, j, k)*c2x(i, j, k) + &
                 vof1(i + ii, j + jj, k + kk)*c1x(i + ii, j + jj, k + kk)
            my = vof3(i - ii, j - jj, k - kk)*c3y(i - ii, j - jj, k - kk) + &
                 vof2(i, j, k)*c2y(i, j, k) + &
                 vof1(i + ii, j + jj, k + kk)*c1y(i + ii, j + jj, k + kk)
            mz = vof3(i - ii, j - jj, k - kk)*c3z(i - ii, j - jj, k - kk) + &
                 vof2(i, j, k)*c2z(i, j, k) + &
                 vof1(i + ii, j + jj, k + kk)*c1z(i + ii, j + jj, k + kk)
            mx = mx*exp_factor
            my = my*exp_factor
            mz = mz*exp_factor
            if (dir == 1) mx = (mx + f(i, j, k)*a1)*exp_factor
            if (dir == 2) my = (my + f(i, j, k)*a1)*exp_factor
            if (dir == 3) mz = (mz + f(i, j, k)*a1)*exp_factor
            mx = min(max(mx, 0.0_sp), 0.5_sp)
            my = min(max(my, 0.0_sp), 0.5_sp)
            mz = min(max(mz, 0.0_sp), 0.5_sp)
            cx(i, j, k) = mx/f(i, j, k)
            cy(i, j, k) = my/f(i, j, k)
            cz(i, j, k) = mz/f(i, j, k)
            cx(i, j, k) = min(max(cx(i, j, k), 0.0_sp), 1.0_sp)
            cy(i, j, k) = min(max(cy(i, j, k), 0.0_sp), 1.0_sp)
            cz(i, j, k) = min(max(cz(i, j, k), 0.0_sp), 1.0_sp)
          endif
        enddo
      enddo
    enddo

    Call CPU_TIME(tt2)  !!! cpu time for correction term
    t_cor = t_cor + tt2 - tt1

    ! apply proper boundary conditions to c
    call phi_bc%setBCS(f)
    !***
    return
  End Subroutine MOFAdv_EI

  !=======================================================
  ! (2-7) Advection algorithm of Weymouth-Dick correction
  !=======================================================
  Subroutine MOFAdv_WY(us, f0, cx, cy, cz, f, nl, dl, dt, dir)
    Implicit None
    Integer :: dir
    Real(sp) :: dt
    Integer, Intent(In)     :: nl(3)
    Real(sp), Intent(In)    :: dl(3)
    Real(sp), Intent(In)    :: us(0:nl(1) + 1, 0:nl(2) + 1, 0:nl(3) + 1)
    Real(sp), Intent(In)    :: f0(0:nl(1) + 1, 0:nl(2) + 1, 0:nl(3) + 1)
    Real(sp), Intent(InOut) :: f(0:nl(1) + 1, 0:nl(2) + 1, 0:nl(3) + 1)
    Real(sp), Intent(InOut) :: cx(0:nl(1) + 1, 0:nl(2) + 1, 0:nl(3) + 1)
    Real(sp), Intent(InOut) :: cy(0:nl(1) + 1, 0:nl(2) + 1, 0:nl(3) + 1)
    Real(sp), Intent(InOut) :: cz(0:nl(1) + 1, 0:nl(2) + 1, 0:nl(3) + 1)
    Integer  :: i, j, k
    Integer  :: ii, jj, kk
    Real(sp) :: a1, a2, alpha
    Real(sp) :: b1, b2
    Real(sp) :: c1, c2
    Real(sp) :: mx, my, mz
    Real(sp) :: x0(3), deltax(3)
    Real(sp) :: cxyz(3), c1xyz(3), c2xyz(3), c3xyz(3)
    Real(sp), dimension(0:nl(1) + 1, 0:nl(2) + 1, 0:nl(3) + 1) :: vof1, vof2, vof3
    Real(sp), dimension(0:nl(1) + 1, 0:nl(2) + 1, 0:nl(3) + 1) :: c1x, c1y, c1z
    Real(sp), dimension(0:nl(1) + 1, 0:nl(2) + 1, 0:nl(3) + 1) :: c2x, c2y, c2z
    Real(sp), dimension(0:nl(1) + 1, 0:nl(2) + 1, 0:nl(3) + 1) :: c3x, c3y, c3z
    Real(sp), dimension(0:nl(1) + 1, 0:nl(2) + 1, 0:nl(3) + 1) :: n1, n2, n3
    Integer, dimension(0:nl(1) + 1, 0:nl(2) + 1, 0:nl(3) + 1) :: flag
    Real(sp) :: norm(3)
    Real(sp) :: f_block(3, 3, 3)
    Real(sp) :: exp_factor
    Real(sp) :: tt1, tt2

    vof1 = 0.0_sp
    vof3 = 0.0_sp
    ii = 0; jj=0; kk=0
    If (dir .eq. 1) ii = 1
    If (dir .eq. 2) jj = 1
    If (dir .eq. 3) kk = 1

    ! default, f = 0
    Call MOFNormalVectors(f, cx, cy, cz, nl, flag, n1, n2, n3)
    Call MOF_EulerianFlux(us, f, cx, cy, cz, nl, dl, dt, dir, &
                              n1, n2, n3, flag, &
                              c1x, c1y, c1z, &
                              c2x, c2y, c2z, &
                              c3x, c3y, c3z, &
                              vof1, vof2, vof3)

    Call CPU_TIME(tt1)  !!! cpu time for correction term

    !new values of f and  clip it: 0. <= f <= 1.
    do k = 1, nl(3)
      do j = 1, nl(2)
        do i = 1, nl(1)
          a1 = us(i - ii, j - jj, k - kk)*dt/dl(dir)
          a2 = us(i, j, k)*dt/dl(dir)
          f(i, j, k) =  vof3(i - ii, j - jj, k - kk) + &
                        vof2(i, j, k)  + &
                        vof1(i + ii, j + jj, k + kk) 
          if ( f0(i,j,k) > 0.5_sp) f(i,j,k) = f(i,j,k) + a2 - a1
          exp_factor = 1.0_sp/(1.0_sp - a2 + a1)
          if (f(i, j, k) < EPSC_VOF) then
            f(i, j, k) = 0.0_sp
            cx(i, j, k) = 0.0_sp
            cy(i, j, k) = 0.0_sp
            cz(i, j, k) = 0.0_sp
          elseif (f(i, j, k) > (1.0_sp - EPSC_VOF)) then
            f(i, j, k) = 1.0_sp
            cx(i, j, k) = 0.5_sp
            cy(i, j, k) = 0.5_sp
            cz(i, j, k) = 0.5_sp
          else
            mx = vof3(i - ii, j - jj, k - kk)*c3x(i - ii, j - jj, k - kk) + &
                 vof2(i, j, k)*c2x(i, j, k) + &
                 vof1(i + ii, j + jj, k + kk)*c1x(i + ii, j + jj, k + kk)
            my = vof3(i - ii, j - jj, k - kk)*c3y(i - ii, j - jj, k - kk) + &
                 vof2(i, j, k)*c2y(i, j, k) + &
                 vof1(i + ii, j + jj, k + kk)*c1y(i + ii, j + jj, k + kk)
            mz = vof3(i - ii, j - jj, k - kk)*c3z(i - ii, j - jj, k - kk) + &
                 vof2(i, j, k)*c2z(i, j, k) + &
                 vof1(i + ii, j + jj, k + kk)*c1z(i + ii, j + jj, k + kk)
            mx = mx*exp_factor
            my = my*exp_factor
            mz = mz*exp_factor
            if (dir == 1) mx = (mx + f(i, j, k)*a1)*exp_factor
            if (dir == 2) my = (my + f(i, j, k)*a1)*exp_factor
            if (dir == 3) mz = (mz + f(i, j, k)*a1)*exp_factor
            mx = min(max(mx, 0.0_sp), 0.5_sp)
            my = min(max(my, 0.0_sp), 0.5_sp)
            mz = min(max(mz, 0.0_sp), 0.5_sp)
            cx(i, j, k) = mx/f(i, j, k)
            cy(i, j, k) = my/f(i, j, k)
            cz(i, j, k) = mz/f(i, j, k)
            cx(i, j, k) = min(max(cx(i, j, k), 0.0_sp), 1.0_sp)
            cy(i, j, k) = min(max(cy(i, j, k), 0.0_sp), 1.0_sp)
            cz(i, j, k) = min(max(cz(i, j, k), 0.0_sp), 1.0_sp)
            ! f(i,j,k) = min(max(f(i,j,k),0.0_sp),1.0_sp)
          endif
        enddo
      enddo
    enddo
    Call CPU_TIME(tt2)  !!! cpu time for correction term
    t_cor = t_cor + tt2 - tt1
    ! apply proper boundary conditions to c
    call phi_bc%setBCS(f)
    !***
    return
  End Subroutine MOFAdv_WY

  Subroutine VOFNormalVectors(f, nl, flag, n1, n2, n3)
    Implicit None
    Integer, Intent(In)     :: nl(3)
    Real(sp), dimension(0:nl(1) + 1, 0:nl(2) + 1, 0:nl(3) + 1), Intent(In) :: f
    Real(sp), dimension(0:nl(1) + 1, 0:nl(2) + 1, 0:nl(3) + 1), Intent(Out) :: n1, n2, n3
    Integer, dimension(0:nl(1) + 1, 0:nl(2) + 1, 0:nl(3) + 1), Intent(Out) :: flag
    Integer :: i, j, k
    Real(sp) :: f_block(3,3,3)
    Real(sp) :: norm(3)
    Real(sp) :: tt1, tt2

    Call CPU_TIME(tt1)  !!! cpu time for reconstruction

    Do k = 1, nl(3)
      Do j = 1, nl(2)
        Do i = 1, nl(1)
          If (f(i, j, k) .GE. 1.0_sp - epsc_vof) Then
            flag(i,j,k) = 0
          Else If (f(i, j, k) .GE. epsc) Then
            flag(i,j,k) = 1
            !*(1)* normal vector
            f_block = f(i - 1:i + 1, j - 1:j + 1, k - 1:k + 1)
            Call VOFNorm(f_block, norm)
            Call Normalization1(norm)
            n1(i,j,k) = norm(1)
            n2(i,j,k) = norm(2)
            n3(i,j,k) = norm(3)
            grid_iter = grid_iter + 1
          Else
            flag(i,j,k) = -1
          EndIf
        EndDo 
      EndDo 
    EndDo 

    Call CPU_TIME(tt2)  !!! cpu time for reconstruction
    t_rec = t_rec + tt2 - tt1

  End Subroutine 

  Subroutine MOFNormalVectors(f,cx, cy, cz, nl, flag, n1, n2, n3)
    Implicit None
    Integer, Intent(In)     :: nl(3)
    Real(sp), dimension(0:nl(1) + 1, 0:nl(2) + 1, 0:nl(3) + 1), Intent(In) :: f
    Real(sp), dimension(0:nl(1) + 1, 0:nl(2) + 1, 0:nl(3) + 1), Intent(In) :: cx, cy, cz
    Real(sp), dimension(0:nl(1) + 1, 0:nl(2) + 1, 0:nl(3) + 1), Intent(Out) :: n1, n2, n3
    Integer, dimension(0:nl(1) + 1, 0:nl(2) + 1, 0:nl(3) + 1), Intent(Out) :: flag
    Integer :: i, j, k
    Real(sp) :: f_block(3,3,3)
    Real(sp) :: norm(3)
    Real(sp) :: c3(3)
    Real(sp) :: tt1, tt2

    Call CPU_TIME(tt1)  !!! cpu time for reconstruction

    Do k = 1, nl(3)
      Do j = 1, nl(2)
        Do i = 1, nl(1)
          If (f(i, j, k) .GE. 1.0_sp - epsc_vof) Then
            flag(i,j,k) = 0
          Else If (f(i, j, k) .GE. epsc) Then
            flag(i,j,k) = 1
            !*(1)* normal vector
            c3(1) = cx(i, j, k)
            c3(2) = cy(i, j, k)
            c3(3) = cz(i, j, k)
            Call MOFNorm(f(i, j, k), c3, norm)
            Call Normalization1(norm)
            n1(i,j,k) = norm(1)
            n2(i,j,k) = norm(2)
            n3(i,j,k) = norm(3)
            num_iter = num_iter + mof_niter
            grid_iter = grid_iter + 1
          Else
            flag(i,j,k) = -1
          EndIf
        EndDo 
      EndDo 
    EndDo 

    Call CPU_TIME(tt2)  !!! cpu time for reconstruction
    t_rec = t_rec + tt2 - tt1

  End Subroutine 

  Subroutine VOF_EulerianFlux(us, f, nl, dl, dt, dir, &
                              n1, n2, n3, flag, &
                              vof1, vof2, vof3)
    Implicit None
    Integer :: dir
    Real(sp) :: dt
    Integer, Intent(In)     :: nl(3)
    Real(sp), Intent(In)    :: dl(3)
    Real(sp), Intent(In)    :: us(0:nl(1) + 1, 0:nl(2) + 1, 0:nl(3) + 1)
    Real(sp), Intent(In) :: f(0:nl(1) + 1, 0:nl(2) + 1, 0:nl(3) + 1)
    Real(sp), dimension(0:nl(1) + 1, 0:nl(2) + 1, 0:nl(3) + 1), Intent(In) :: n1, n2, n3
    Integer, dimension(0:nl(1) + 1, 0:nl(2) + 1, 0:nl(3) + 1), Intent(In) :: flag
    Real(sp), dimension(0:nl(1) + 1, 0:nl(2) + 1, 0:nl(3) + 1), Intent(Out) :: vof1, vof2, vof3
    Integer  :: i, j, k, ii, jj, kk
    Real(sp) :: a1, a2, alpha
    Real(sp) :: x0(3), deltax(3)
    Real(sp) :: norm(3)
    Real(sp) :: tt1, tt2

    ! default, f = 0
    vof1 = 0.0_sp
    vof3 = 0.0_sp

    ii = 0; jj=0; kk=0
    If (dir .eq. 1) ii = 1
    If (dir .eq. 2) jj = 1
    If (dir .eq. 3) kk = 1

    Call CPU_TIME(tt1)  !!! cpu time for reconstruction

    !***
    Do k = 1, nl(3)
      Do j = 1, nl(2)
        Do i = 1, nl(1)
          a1 = us(i - ii, j - jj, k - kk)*dt/dl(dir)
          a2 = us(i, j, k)*dt/dl(dir)

          vof1(i, j, k) = 0.d0
          vof3(i, j, k) = 0.d0
          vof2(i, j, k) = 0.d0
          ! f = 1
          If (flag(i, j, k) .EQ. 0) Then
            vof1(i, j, k) = DMAX1(-a1, 0.0_sp)
            vof3(i, j, k) = DMAX1(a2, 0.0_sp)

            ! 0 < f < 1
          Else If (flag(i, j, k) .EQ. 1) Then
            !*(2) get alpha;
            norm(1) = n1(i,j,k); norm(2) = n2(i,j,k); norm(3) = n3(i,j,k)
            alpha = FloodSZ_Backward(norm, f(i, j, k))
            !*(3) get fluxes
            x0 = 0.0_sp; deltax = 1.0_sp
            If (a1 .LT. 0.0_sp) Then
              deltax(dir) = -a1
              vof1(i, j, k) = FloodSZ_Forward(norm, alpha, x0, deltax)
            End If
            If (a2 .GT. 0.0_sp) Then
              x0(dir) = 1.0_sp - a2; 
              deltax(dir) = a2
              vof3(i, j, k) = FloodSZ_Forward(norm, alpha, x0, deltax)
            end If
          EndIf
          vof2(i, j, k) = f(i, j, k) - vof1(i, j, k) - vof3(i, j, k)
        EndDo
      EndDo
    EndDo

    Call CPU_TIME(tt2)  !!! cpu time for reconstruction
    t_adv = t_adv + tt2 - tt1

    ! apply proper boundary conditions to vof1, vof3
    Call phi_bc%setBCS(vof1)
    Call phi_bc%setBCS(vof3)

  End Subroutine VOF_EulerianFlux

  Subroutine VOF_LagrangianFlux(us, f, nl, dl, dt, dir, &
                              n1, n2, n3, flag, &
                              vof1, vof2, vof3)
    Implicit None
    Integer :: dir
    Real(sp) :: dt
    Integer, Intent(In)     :: nl(3)
    Real(sp), Intent(In)    :: dl(3)
    Real(sp), Intent(In)    :: us(0:nl(1) + 1, 0:nl(2) + 1, 0:nl(3) + 1)
    Real(sp), Intent(InOut) :: f(0:nl(1) + 1, 0:nl(2) + 1, 0:nl(3) + 1)
    Real(sp), dimension(0:nl(1) + 1, 0:nl(2) + 1, 0:nl(3) + 1), Intent(In) :: n1, n2, n3
    Integer, dimension(0:nl(1) + 1, 0:nl(2) + 1, 0:nl(3) + 1), Intent(In) :: flag
    Real(sp), dimension(0:nl(1) + 1, 0:nl(2) + 1, 0:nl(3) + 1), Intent(Out) :: vof1, vof2, vof3
    Integer  :: i, j, k, ii, jj, kk
    Real(sp) :: a1, a2, alpha
    Real(sp) :: x0(3), deltax(3)
    Real(sp) :: norm(3)
    Real(sp) :: tt1, tt2

    ! default, f = 0
    vof1 = 0.0_sp
    vof2 = 0.0_sp
    vof3 = 0.0_sp
    ii = 0; jj=0; kk=0
    if (dir .eq. 1) ii = 1
    if (dir .eq. 2) jj = 1
    if (dir .eq. 3) kk = 1

    Call CPU_TIME(tt1)  !!! cpu time for reconstruction

    !***
    do k = 1, nl(3)
      do j = 1, nl(2)
        do i = 1, nl(1)
          a1 = us(i - ii, j - jj, k - kk)*dt/dl(dir)
          a2 = us(i, j, k)*dt/dl(dir)

          ! f = 1
          if (flag(i, j, k) .eq. 0) then
            vof1(i, j, k) = DMAX1(-a1, 0.0_sp)
            vof2(i, j, k) = 1.0_sp - DMAX1(a1, 0.0_sp) + DMIN1(a2, 0.0_sp)
            vof3(i, j, k) = DMAX1(a2, 0.0_sp)

            ! 0 < f < 1
          else if (flag(i, j, k) .eq. 1) then
            !*(2) get alpha;
            norm(1) = n1(i,j,k); norm(2) = n2(i,j,k); norm(3) = n3(i,j,k)
            alpha = FloodSZ_Backward(norm, f(i, j, k))
            norm(dir) = norm(dir)/(1.0_sp - a1 + a2)
            alpha = alpha + norm(dir)*a1
            !*(3) get fluxes
            x0 = 0.0_sp; deltax = 1.0_sp
            if (a1 .LT. 0.0_sp) then
              x0(dir) = a1; 
              deltax(dir) = -a1
              vof1(i, j, k) = FloodSZ_Forward(norm, alpha, x0, deltax)
            end if
            if (a2 .GT. 0.0_sp) then
              x0(dir) = 1.0_sp; 
              deltax(dir) = a2
              vof3(i, j, k) = FloodSZ_Forward(norm, alpha, x0, deltax)
            end if
            x0(dir) = DMAX1(a1, 0.0_sp)
            deltax(dir) = 1.0_sp - x0(dir) + DMIN1(0.0_sp, a2)
            vof2(i, j, k) = FloodSZ_Forward(norm, alpha, x0, deltax)
          endif
        enddo
      enddo
    enddo

    Call CPU_TIME(tt2)  !!! cpu time for reconstruction
    t_adv = t_adv + tt2 - tt1

    ! apply proper boundary conditions to vof1, vof2, vof3
    call phi_bc%setBCS(vof1)
    call phi_bc%setBCS(vof3)

  End Subroutine VOF_LagrangianFlux

  Subroutine MOF_EulerianFlux(us, f, cx, cy, cz, nl, dl, dt, dir, &
                              n1, n2, n3, flag, &
                              c1x, c1y, c1z, &
                              c2x, c2y, c2z, &
                              c3x, c3y, c3z, &
                              vof1, vof2, vof3)
    Implicit None
    Integer :: dir
    Real(sp) :: dt
    Integer, Intent(In)     :: nl(3)
    Real(sp), Intent(In)    :: dl(3)
    Real(sp), Intent(In)    :: us(0:nl(1) + 1, 0:nl(2) + 1, 0:nl(3) + 1)
    Real(sp), Intent(InOut) :: f(0:nl(1) + 1, 0:nl(2) + 1, 0:nl(3) + 1)
    Real(sp), dimension(0:nl(1) + 1, 0:nl(2) + 1, 0:nl(3) + 1), Intent(In) :: n1, n2, n3
    Real(sp), dimension(0:nl(1) + 1, 0:nl(2) + 1, 0:nl(3) + 1), Intent(In) :: cx, cy, cz
    Integer, dimension(0:nl(1) + 1, 0:nl(2) + 1, 0:nl(3) + 1), Intent(In) :: flag
    Real(sp), dimension(0:nl(1) + 1, 0:nl(2) + 1, 0:nl(3) + 1), Intent(Out) :: vof1, vof2, vof3
    Real(sp), dimension(0:nl(1) + 1, 0:nl(2) + 1, 0:nl(3) + 1), Intent(Out) :: c1x, c1y, c1z
    Real(sp), dimension(0:nl(1) + 1, 0:nl(2) + 1, 0:nl(3) + 1), Intent(Out) :: c2x, c2y, c2z
    Real(sp), dimension(0:nl(1) + 1, 0:nl(2) + 1, 0:nl(3) + 1), Intent(Out) :: c3x, c3y, c3z
    Integer  :: i, j, k, ii, jj, kk
    Real(sp) :: a1, a2, alpha
    Real(sp) :: x0(3), deltax(3)
    Real(sp) :: cxyz(3), c1xyz(3), c2xyz(3), c3xyz(3)
    Real(sp) :: norm(3)
    Real(sp) :: tt1, tt2

    ! default, f = 0
    vof1 = 0.0_sp
    vof2 = 0.0_sp
    vof3 = 0.0_sp
    ii = 0; jj=0; kk=0
    If (dir .eq. 1) ii = 1
    If (dir .eq. 2) jj = 1
    If (dir .eq. 3) kk = 1

    Call CPU_TIME(tt1)  !!! cpu time for reconstruction

    !***
    do k = 1, nl(3)
      do j = 1, nl(2)
        do i = 1, nl(1)
          a1 = us(i - ii, j - jj, k - kk)*dt/dl(dir)
          a2 = us(i, j, k)*dt/dl(dir)

          vof1(i, j, k) = 0.d0
          vof3(i, j, k) = 0.d0
          vof2(i, j, k) = 0.d0
          c1xyz = 0.0_sp
          c2xyz = 0.0_sp
          c3xyz = 0.0_sp
          ! f = 1
          If (flag(i, j, k) .EQ. 0) Then
            vof1(i, j, k) = DMAX1(-a1, 0.0_sp)
            vof3(i, j, k) = DMAX1(a2, 0.0_sp)
            c1xyz = 0.5_sp
            c3xyz = 0.5_sp
            c1xyz(dir) = DMAX1(-a1/2.0_sp, 0.0_sp)
            c3xyz(dir) = 1.0_sp - DMAX1(a2/2.0_sp, 0.0_sp)

            vof2(i, j, k) = 1.0_sp - DMAX1(-a1, 0.0_sp) - DMAX1(a2, 0.0_sp)
            c2xyz = 0.5_sp
            c2xyz(dir) = 0.5_sp + DMAX1(-a1/2.0_sp, 0.0_sp) - DMAX1(a2/2.0_sp, 0.0_sp)

            ! 0 < f < 1
          Else if (flag(i, j, k) .EQ. 1) Then
            !*(1)* normal vector
            norm(1) = n1(i,j,k)
            norm(2) = n2(i,j,k)
            norm(3) = n3(i,j,k)
            !*(2) get alpha;
            alpha = FloodSZ_Backward(norm, f(i, j, k))
            !*(3) get fluxes
            x0 = 0.0_sp; deltax = 1.0_sp
            if (a1 .LT. 0.0_sp) then
              deltax(dir) = -a1
              Call ForwardC(norm, alpha, x0, deltax, vof1(i, j, k), c1xyz)
            end if
            if (a2 .GT. 0.0_sp) then
              x0(dir) = 1.0_sp - a2; 
              deltax(dir) = a2
              Call ForwardC(norm, alpha, x0, deltax, vof3(i, j, k), c3xyz)
            end if
            x0(dir) = DMAX1(-a1, 0.0_sp)
            deltax(dir) = 1.0_sp - x0(dir) - DMAX1(0.0_sp, a2)
            Call ForwardC(norm, alpha, x0, deltax, vof2(i, j, k), c2xyz)
          endif
          If (vof1(i, j, k) > epsc_vof) Then
            c1xyz(dir) = c1xyz(dir) + 1.0_sp
          End If
          If (vof3(i, j, k) > epsc_vof) Then
            c3xyz(dir) = c3xyz(dir) - 1.0_sp
          End If

          c1x(i, j, k) = c1xyz(1)
          c1y(i, j, k) = c1xyz(2)
          c1z(i, j, k) = c1xyz(3)
          c3x(i, j, k) = c3xyz(1)
          c3y(i, j, k) = c3xyz(2)
          c3z(i, j, k) = c3xyz(3)
          c2x(i, j, k) = c2xyz(1)
          c2y(i, j, k) = c2xyz(2)
          c2z(i, j, k) = c2xyz(3)

        enddo
      enddo
    enddo

    Call CPU_TIME(tt2)  !!! cpu time for reconstruction
    t_adv = t_adv + tt2 - tt1

    ! apply proper boundary conditions to vof1, vof3
    call phi_bc%setBCS(vof1)
    call phi_bc%setBCS(vof3)
    call phi_bc%setBCS(c1x)
    call phi_bc%setBCS(c3x)
    call phi_bc%setBCS(c1y)
    call phi_bc%setBCS(c3y)
    call phi_bc%setBCS(c1z)
    call phi_bc%setBCS(c3z)

  End Subroutine MOF_EulerianFlux

  Subroutine MOF_LagrangianFlux(us, f, cx, cy, cz, nl, dl, dt, dir, &
                              n1, n2, n3, flag, &
                              c1x, c1y, c1z, &
                              c2x, c2y, c2z, &
                              c3x, c3y, c3z, &
                              vof1, vof2, vof3)
    Implicit None
    Integer :: dir
    Real(sp) :: dt
    Integer, Intent(In)     :: nl(3)
    Real(sp), Intent(In)    :: dl(3)
    Real(sp), Intent(In)    :: us(0:nl(1) + 1, 0:nl(2) + 1, 0:nl(3) + 1)
    Real(sp), Intent(InOut) :: f(0:nl(1) + 1, 0:nl(2) + 1, 0:nl(3) + 1)
    Real(sp), dimension(0:nl(1) + 1, 0:nl(2) + 1, 0:nl(3) + 1), Intent(In) :: n1, n2, n3
    Real(sp), dimension(0:nl(1) + 1, 0:nl(2) + 1, 0:nl(3) + 1), Intent(In) :: cx, cy, cz
    Integer, dimension(0:nl(1) + 1, 0:nl(2) + 1, 0:nl(3) + 1), Intent(In) :: flag
    Real(sp), dimension(0:nl(1) + 1, 0:nl(2) + 1, 0:nl(3) + 1), Intent(Out) :: vof1, vof2, vof3
    Real(sp), dimension(0:nl(1) + 1, 0:nl(2) + 1, 0:nl(3) + 1), Intent(Out) :: c1x, c1y, c1z
    Real(sp), dimension(0:nl(1) + 1, 0:nl(2) + 1, 0:nl(3) + 1), Intent(Out) :: c2x, c2y, c2z
    Real(sp), dimension(0:nl(1) + 1, 0:nl(2) + 1, 0:nl(3) + 1), Intent(Out) :: c3x, c3y, c3z
    Integer  :: i, j, k, ii, jj, kk
    Real(sp) :: a1, a2, alpha
    Real(sp) :: x0(3), deltax(3)
    Real(sp) :: c1xyz(3), c2xyz(3), c3xyz(3)
    Real(sp) :: norm(3)
    Real(sp) :: tt1, tt2

    ! default, f = 0
    vof1 = 0.0_sp
    vof2 = 0.0_sp
    vof3 = 0.0_sp
    ii = 0; jj=0; kk=0
    If (dir .eq. 1) ii = 1
    If (dir .eq. 2) jj = 1
    If (dir .eq. 3) kk = 1

    Call CPU_TIME(tt1)  !!! cpu time for reconstruction

    Do k = 1, nl(3)
      Do j = 1, nl(2)
        Do i = 1, nl(1)

          a1 = us(i - ii, j - jj, k - kk)*dt/dl(dir)
          a2 = us(i, j, k)*dt/dl(dir)
          c1xyz = 0.0_sp
          c2xyz = 0.0_sp
          c3xyz = 0.0_sp
          ! f = 1
          If (flag(i, j, k) .EQ. 0) Then
            vof1(i, j, k) = DMAX1(-a1, 0.0_sp)
            vof2(i, j, k) = 1.0_sp - DMAX1(a1, 0.0_sp) + DMIN1(a2, 0.0_sp)
            vof3(i, j, k) = DMAX1(a2, 0.0_sp)
            c1xyz = 0.5_sp
            c2xyz = 0.5_sp
            c3xyz = 0.5_sp
            c1xyz(dir) = -DMAX1(-a1/2.0_sp, 0.0_sp)
            c2xyz(dir) = 0.5_sp + DMAX1(a1/2.0_sp, 0.0_sp) + DMIN1(a2/2.0_sp, 0.0_sp)
            c3xyz(dir) = 1.0_sp + DMAX1(a2/2.0_sp, 0.0_sp)

            ! 0 < f < 1
          Else if (flag(i, j, k) .EQ. 1) Then
            norm(1) = n1(i,j,k)
            norm(2) = n2(i,j,k)
            norm(3) = n3(i,j,k)
            !*(2) get alpha;
            alpha = FloodSZ_Backward(norm, f(i, j, k))
            norm(dir) = norm(dir)/(1.0_sp - a1 + a2)
            alpha = alpha + norm(dir)*a1
            !*(3) get fluxes
            x0 = 0.0_sp; deltax = 1.0_sp
            If (a1 .LT. 0.0_sp) Then
              x0(dir) = a1; 
              deltax(dir) = -a1
              Call ForwardC(norm, alpha, x0, deltax, vof1(i, j, k), c1xyz)
            End If
            If (a2 .GT. 0.0_sp) Then
              x0(dir) = 1.0_sp; 
              deltax(dir) = a2
              Call ForwardC(norm, alpha, x0, deltax, vof3(i, j, k), c3xyz)
            End If
            x0(dir) = DMAX1(a1, 0.0_sp)
            deltax(dir) = 1.0_sp - x0(dir) + DMIN1(0.0_sp, a2)
            Call ForwardC(norm, alpha, x0, deltax, vof2(i, j, k), c2xyz)
          EndIf

          If (vof1(i, j, k) .ge. epsc_vof) Then
            c1xyz(dir) = c1xyz(dir) + 1.0_sp
          End If
          If (vof3(i, j, k) .ge. epsc_vof) Then
            c3xyz(dir) = c3xyz(dir) - 1.0_sp
          End If

          c1x(i, j, k) = c1xyz(1)
          c1y(i, j, k) = c1xyz(2)
          c1z(i, j, k) = c1xyz(3)
          c3x(i, j, k) = c3xyz(1)
          c3y(i, j, k) = c3xyz(2)
          c3z(i, j, k) = c3xyz(3)
          c2x(i, j, k) = c2xyz(1)
          c2y(i, j, k) = c2xyz(2)
          c2z(i, j, k) = c2xyz(3)

        enddo
      enddo
    enddo

    Call CPU_TIME(tt2)  !!! cpu time for reconstruction
    t_adv = t_adv + tt2 - tt1

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

  End Subroutine MOF_LagrangianFlux

End Module ModVOF