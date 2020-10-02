!=================
! Includes:
!     (1) VOF/MOF Advection
!       (1-1) VOF-LE (Lagrangian Explicit)
!       (1-2) VOF-EI (Eulerian Implicit)
!       (1-3) VOF-EILE3DS (Eulerian Explicit Lanrangian Ixplicit)
!       (1-4) MOF-LE (Lanrangian Explicit)
!       (1-5) MOF-EI (Eulerian Impllicit)
!       (1-6) MOF-EILE3DS (Eulerian Explicit Lanrangian Ixplicit)
!       (1-7) THINC-EI (Eulerian Implicit)
!     (2) Directional split VOF/MOF advection
!       (2-1) VOF-LE (Lagrangian Explicit)
!       (2-2) VOF-EI (Eulerian Implicit)
!       (2-3) MOF-LE (Lanrangian Explicit)
!       (2-4) MOF-EI (Eulerian Implicit)
!       (2-3) THINC-EI (Eulerian Implicit)
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
  Use ModVOFFunc
#if defined(VISUALIZE)
  Use ModTools
#endif

  Interface VOFAdvection
    Module Procedure VOF_EILE3DS
    Module Procedure MOF_EILE3DS
  End Interface VOFAdvection

  real(sp) :: epsc_vof = 1.0d-12

  Integer(sp) :: num_iter(2) = 0
  Integer(sp) :: grid_iter = 0

  ! Originally, the names of the subroutines are
  ! AdvCIAM, AdvWY ...
  !! This part allows the old cases to call the modified subroutine name
  ! May delete it later
  PROCEDURE(InterfaceMOFTIME), POINTER :: MOFCIAM => MOF_LE
  PROCEDURE(InterfaceMOFTIME), POINTER :: MOFWY => MOF_EI
  PROCEDURE(InterfaceVOFTIME), POINTER :: VOFCIAM => VOF_LE
  PROCEDURE(InterfaceVOFTIME), POINTER :: VOFWY => VOF_EI
  PROCEDURE(InterfaceMOFDIR),  POINTER :: AdvCIAM_MOF => MOFAdv_LE
  PROCEDURE(InterfaceMOFDIR),  POINTER :: AdvWY_MOF => MOFAdv_EI
  PROCEDURE(InterfaceVOFDIR),  POINTER :: AdvCIAM => VOFAdv_LE
  PROCEDURE(InterfaceVOFDIR),  POINTER :: AdvWY => VOFAdv_EI

  Interface
    Subroutine InterfaceVOFTIME(Phi, u, v, w, nl, dl, dt, rank)
      import
      Implicit None
      Real(sp) :: dt
      Integer, Intent(In)     :: nl(3)
      Real(sp), Intent(In)    :: dl(3)
      Real(sp), Intent(InOut)    :: Phi(0:nl(1) + 1, 0:nl(2) + 1, 0:nl(3) + 1)
      Real(sp), Intent(In)       :: u(0:nl(1) + 1, 0:nl(2) + 1, 0:nl(3) + 1)
      Real(sp), Intent(In)       :: v(0:nl(1) + 1, 0:nl(2) + 1, 0:nl(3) + 1)
      Real(sp), Intent(In)       :: w(0:nl(1) + 1, 0:nl(2) + 1, 0:nl(3) + 1)
      Integer :: rank
    End Subroutine InterfaceVOFTIME

    Subroutine InterfaceMOFTIME(Phi, u, v, w, nl, dl, dt, rank, cx, cy, cz)
      import
      Implicit None
      Real(sp) :: dt
      Integer, Intent(In)     :: nl(3)
      Real(sp), Intent(In)    :: dl(3)
      Real(sp), Intent(InOut)    :: Phi(0:nl(1) + 1, 0:nl(2) + 1, 0:nl(3) + 1)
      Real(sp), Intent(InOut)    :: cx(0:nl(1) + 1, 0:nl(2) + 1, 0:nl(3) + 1)
      Real(sp), Intent(InOut)    :: cy(0:nl(1) + 1, 0:nl(2) + 1, 0:nl(3) + 1)
      Real(sp), Intent(InOut)    :: cz(0:nl(1) + 1, 0:nl(2) + 1, 0:nl(3) + 1)
      Real(sp), Intent(In)       :: u(0:nl(1) + 1, 0:nl(2) + 1, 0:nl(3) + 1)
      Real(sp), Intent(In)       :: v(0:nl(1) + 1, 0:nl(2) + 1, 0:nl(3) + 1)
      Real(sp), Intent(In)       :: w(0:nl(1) + 1, 0:nl(2) + 1, 0:nl(3) + 1)
      Integer :: rank
    End Subroutine InterfaceMOFTIME

    Subroutine InterfaceVOFDIR(us, f, nl, dl, dt, dir)
      import
      Implicit None
      Integer :: dir
      Real(sp) :: dt
      Integer, Intent(In)     :: nl(3)
      Real(sp), Intent(In)    :: dl(3)
      Real(sp), Intent(In)    :: us(0:nl(1) + 1, 0:nl(2) + 1, 0:nl(3) + 1)
      Real(sp), Intent(InOut) :: f(0:nl(1) + 1, 0:nl(2) + 1, 0:nl(3) + 1)
    End subroutine InterfaceVOFDIR

    Subroutine InterfaceMOFDIR(us, cx, cy, cz, f, nl, dl, dt, dir)
      import
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
    End Subroutine InterfaceMOFDIR

  End Interface

  Contains


! (1-1) VOF-LE (Lagrangian Explicit)
  Subroutine VOF_LE(Phi, u, v, w, nl, dl, dt, rank)
    Implicit None
    Real(sp) :: dt
    Integer, Intent(In)     :: nl(3)
    Real(sp), Intent(In)    :: dl(3)
    Real(sp), Intent(InOut)    :: Phi(0:nl(1) + 1, 0:nl(2) + 1, 0:nl(3) + 1)
    Real(sp), Intent(In)       :: u(0:nl(1) + 1, 0:nl(2) + 1, 0:nl(3) + 1)
    Real(sp), Intent(In)       :: v(0:nl(1) + 1, 0:nl(2) + 1, 0:nl(3) + 1)
    Real(sp), Intent(In)       :: w(0:nl(1) + 1, 0:nl(2) + 1, 0:nl(3) + 1)
    Integer :: rank
    if (rank == 0 .or. rank == 3) Then
      call VOFAdv_LE(u, phi, nl, dl, dt, 1)
      call VOFAdv_LE(v, phi, nl, dl, dt, 2)
      call VOFAdv_LE(w, phi, nl, dl, dt, 3)
    else if (rank == 1 .or. rank == 4) Then
      call VOFAdv_LE(v, phi, nl, dl, dt, 2)
      call VOFAdv_LE(w, phi, nl, dl, dt, 3)
      call VOFAdv_LE(u, phi, nl, dl, dt, 1)
    else if (rank == 2 .or. rank == 5) Then
      call VOFAdv_LE(w, phi, nl, dl, dt, 3)
      call VOFAdv_LE(u, phi, nl, dl, dt, 1)
      call VOFAdv_LE(v, phi, nl, dl, dt, 2)
    endif
  End Subroutine VOF_LE

!=======================================================
! (1-2) VOF-EI (Eulerian Implicit)
!=======================================================
  Subroutine VOF_EI(Phi, u, v, w, nl, dl, dt, rank)
    Implicit None
    Real(sp) :: dt
    Integer, Intent(In)     :: nl(3)
    Real(sp), Intent(In)    :: dl(3)
    Real(sp), Intent(InOut)    :: Phi(0:nl(1) + 1, 0:nl(2) + 1, 0:nl(3) + 1)
    Real(sp), Intent(In)       :: u(0:nl(1) + 1, 0:nl(2) + 1, 0:nl(3) + 1)
    Real(sp), Intent(In)       :: v(0:nl(1) + 1, 0:nl(2) + 1, 0:nl(3) + 1)
    Real(sp), Intent(In)       :: w(0:nl(1) + 1, 0:nl(2) + 1, 0:nl(3) + 1)
    Integer :: rank
    if (rank == 1 .or. rank == 4) Then
      call VOFAdv_EI(u, phi, nl, dl, dt, 1)
      call VOFAdv_EI(v, phi, nl, dl, dt, 2)
      call VOFAdv_EI(w, phi, nl, dl, dt, 3)
    else if (rank == 2 .or. rank == 5) Then
      call VOFAdv_EI(v, phi, nl, dl, dt, 2)
      call VOFAdv_EI(w, phi, nl, dl, dt, 3)
      call VOFAdv_EI(u, phi, nl, dl, dt, 1)
    else if (rank == 0 .or. rank == 3) Then
      call VOFAdv_EI(w, phi, nl, dl, dt, 3)
      call VOFAdv_EI(u, phi, nl, dl, dt, 1)
      call VOFAdv_EI(v, phi, nl, dl, dt, 2)
    endif
  End Subroutine VOF_EI

!=======================================================
! (1-3) VOF-EILE3DS (Eulerian Explicit Lanrangian Ixplicit)
!=======================================================
  Subroutine VOF_EILE3DS(Phi, u, v, w, nl, dl, dt, rank)
    Implicit None
    Real(sp) :: dt
    Integer, Intent(In)     :: nl(3)
    Real(sp), Intent(In)    :: dl(3)
    Real(sp), Intent(InOut)    :: Phi(0:nl(1) + 1, 0:nl(2) + 1, 0:nl(3) + 1)
    Real(sp), Intent(In)       :: u(0:nl(1) + 1, 0:nl(2) + 1, 0:nl(3) + 1)
    Real(sp), Intent(In)       :: v(0:nl(1) + 1, 0:nl(2) + 1, 0:nl(3) + 1)
    Real(sp), Intent(In)       :: w(0:nl(1) + 1, 0:nl(2) + 1, 0:nl(3) + 1)
    Real(sp)       :: u1(0:nl(1) + 1, 0:nl(2) + 1, 0:nl(3) + 1)
    Real(sp)       :: v1(0:nl(1) + 1, 0:nl(2) + 1, 0:nl(3) + 1)
    Real(sp)       :: w1(0:nl(1) + 1, 0:nl(2) + 1, 0:nl(3) + 1)
    Integer :: rank
    if (rank == 0 .or. rank == 2 .or. rank == 4) Then
      call VOFAdv_EI(u, phi, nl, dl, dt, 1)
      call VOFAdv_EI(v, phi, nl, dl, dt, 2)
      call VOFAdv_EI(w, phi, nl, dl, dt, 3)
    else if (rank == 1 .or. rank == 3 .or. rank == 5) Then
      call VOFAdv_LE(w, phi, nl, dl, dt, 3)
      call VOFAdv_LE(v, phi, nl, dl, dt, 2)
      call VOFAdv_LE(u, phi, nl, dl, dt, 1)
    endif

  End Subroutine VOF_EILE3DS

!=======================================================
! (1-4) MOF-LE (Lanrangian Explicit)
!=======================================================
  Subroutine MOF_LE(Phi, u, v, w, nl, dl, dt, rank, cx, cy, cz)
    Implicit None
    Real(sp) :: dt
    Integer, Intent(In)     :: nl(3)
    Real(sp), Intent(In)    :: dl(3)
    Real(sp), Intent(InOut)    :: Phi(0:nl(1) + 1, 0:nl(2) + 1, 0:nl(3) + 1)
    Real(sp), Intent(InOut)    :: cx(0:nl(1) + 1, 0:nl(2) + 1, 0:nl(3) + 1)
    Real(sp), Intent(InOut)    :: cy(0:nl(1) + 1, 0:nl(2) + 1, 0:nl(3) + 1)
    Real(sp), Intent(InOut)    :: cz(0:nl(1) + 1, 0:nl(2) + 1, 0:nl(3) + 1)
    Real(sp), Intent(In)       :: u(0:nl(1) + 1, 0:nl(2) + 1, 0:nl(3) + 1)
    Real(sp), Intent(In)       :: v(0:nl(1) + 1, 0:nl(2) + 1, 0:nl(3) + 1)
    Real(sp), Intent(In)       :: w(0:nl(1) + 1, 0:nl(2) + 1, 0:nl(3) + 1)
    Integer :: rank
    num_iter = 0
    if (rank == 0 .or. rank == 4) Then
      call MOFAdv_LE(u, cx, cy, cz, phi, nl, dl, dt, 1)
      call MOFAdv_LE(v, cx, cy, cz, phi, nl, dl, dt, 2)
      call MOFAdv_LE(w, cx, cy, cz, phi, nl, dl, dt, 3)
    else if (rank == 1 .or. rank == 5) Then
      call MOFAdv_LE(w, cx, cy, cz, phi, nl, dl, dt, 3)
      call MOFAdv_LE(u, cx, cy, cz, phi, nl, dl, dt, 1)
      call MOFAdv_LE(v, cx, cy, cz, phi, nl, dl, dt, 2)
    else if (rank == 2 .or. rank == 3) Then
      call MOFAdv_LE(v, cx, cy, cz, phi, nl, dl, dt, 2)
      call MOFAdv_LE(w, cx, cy, cz, phi, nl, dl, dt, 3)
      call MOFAdv_LE(u, cx, cy, cz, phi, nl, dl, dt, 1)
    end if
  End Subroutine MOF_LE

!=======================================================
! (1-5) MOF-EI (Eulerian Impllicit)
!=======================================================
  Subroutine MOF_EI(Phi, u, v, w, nl, dl, dt, rank, cx, cy, cz)
    Implicit None
    Real(sp) :: dt
    Integer, Intent(In)     :: nl(3)
    Real(sp), Intent(In)    :: dl(3)
    Real(sp), Intent(InOut)    :: Phi(0:nl(1) + 1, 0:nl(2) + 1, 0:nl(3) + 1)
    Real(sp), Intent(InOut)    :: cx(0:nl(1) + 1, 0:nl(2) + 1, 0:nl(3) + 1)
    Real(sp), Intent(InOut)    :: cy(0:nl(1) + 1, 0:nl(2) + 1, 0:nl(3) + 1)
    Real(sp), Intent(InOut)    :: cz(0:nl(1) + 1, 0:nl(2) + 1, 0:nl(3) + 1)
    Real(sp), Intent(In)       :: u(0:nl(1) + 1, 0:nl(2) + 1, 0:nl(3) + 1)
    Real(sp), Intent(In)       :: v(0:nl(1) + 1, 0:nl(2) + 1, 0:nl(3) + 1)
    Real(sp), Intent(In)       :: w(0:nl(1) + 1, 0:nl(2) + 1, 0:nl(3) + 1)
    Integer :: rank
    if (rank == 0 .or. rank == 4) Then
      call MOFAdv_EI(u, cx, cy, cz, phi, nl, dl, dt, 1)
      call MOFAdv_EI(v, cx, cy, cz, phi, nl, dl, dt, 2)
      call MOFAdv_EI(w, cx, cy, cz, phi, nl, dl, dt, 3)
    else if (rank == 1 .or. rank == 5) Then
      call MOFAdv_EI(w, cx, cy, cz, phi, nl, dl, dt, 3)
      call MOFAdv_EI(u, cx, cy, cz, phi, nl, dl, dt, 1)
      call MOFAdv_EI(v, cx, cy, cz, phi, nl, dl, dt, 2)
    else if (rank == 2 .or. rank == 3) Then
      call MOFAdv_EI(v, cx, cy, cz, phi, nl, dl, dt, 2)
      call MOFAdv_EI(w, cx, cy, cz, phi, nl, dl, dt, 3)
      call MOFAdv_EI(u, cx, cy, cz, phi, nl, dl, dt, 1)
    else
      call MOFAdv_EI(u, cx, cy, cz, phi, nl, dl, dt, 1)
      call MOFAdv_EI(v, cx, cy, cz, phi, nl, dl, dt, 2)
      call MOFAdv_EI(w, cx, cy, cz, phi, nl, dl, dt, 3)
    end if
  End Subroutine MOF_EI

!=======================================================
! (1-6) MOF-EILE3DS (Eulerian Explicit Lanrangian Ixplicit)
!=======================================================
  Subroutine MOF_EILE3DS(Phi, u, v, w, nl, dl, dt, rank, cx, cy, cz)
    Implicit None
    Real(sp) :: dt
    Integer, Intent(In)     :: nl(3)
    Real(sp), Intent(In)    :: dl(3)
    Real(sp), Intent(InOut)    :: Phi(0:nl(1) + 1, 0:nl(2) + 1, 0:nl(3) + 1)
    Real(sp), Intent(InOut)    :: cx(0:nl(1) + 1, 0:nl(2) + 1, 0:nl(3) + 1)
    Real(sp), Intent(InOut)    :: cy(0:nl(1) + 1, 0:nl(2) + 1, 0:nl(3) + 1)
    Real(sp), Intent(InOut)    :: cz(0:nl(1) + 1, 0:nl(2) + 1, 0:nl(3) + 1)
    Real(sp), Intent(In)       :: u(0:nl(1) + 1, 0:nl(2) + 1, 0:nl(3) + 1)
    Real(sp), Intent(In)       :: v(0:nl(1) + 1, 0:nl(2) + 1, 0:nl(3) + 1)
    Real(sp), Intent(In)       :: w(0:nl(1) + 1, 0:nl(2) + 1, 0:nl(3) + 1)
    Integer :: rank
    if (rank == 0 .or. rank == 2 .or. rank == 4) Then
      call MOFAdv_EI(u, cx, cy, cz, phi, nl, dl, dt, 1)
      call MOFAdv_LE(v, cx, cy, cz, phi, nl, dl, dt, 2)
      call MOFAdv_EI(w, cx, cy, cz, phi, nl, dl, dt, 3)
    else if (rank == 1 .or. rank == 3 .or. rank == 5) Then
      call MOFAdv_LE(w, cx, cy, cz, phi, nl, dl, dt, 3)
      call MOFAdv_EI(v, cx, cy, cz, phi, nl, dl, dt, 2)
      call MOFAdv_LE(u, cx, cy, cz, phi, nl, dl, dt, 1)
    endif
  End Subroutine MOF_EILE3DS

!=======================================================
! (1-7) THINC-EI (Eulerian Implicit)
!=======================================================
  Subroutine VOFTHINC_EI(Phi, u, v, w, nl, dl, dt, rank)
    Implicit None
    Real(sp) :: dt
    Integer, Intent(In)     :: nl(3)
    Real(sp), Intent(In)    :: dl(3)
    Real(sp), Intent(InOut)    :: Phi(0:nl(1) + 1, 0:nl(2) + 1, 0:nl(3) + 1)
    Real(sp), Intent(In)       :: u(0:nl(1) + 1, 0:nl(2) + 1, 0:nl(3) + 1)
    Real(sp), Intent(In)       :: v(0:nl(1) + 1, 0:nl(2) + 1, 0:nl(3) + 1)
    Real(sp), Intent(In)       :: w(0:nl(1) + 1, 0:nl(2) + 1, 0:nl(3) + 1)
    Integer :: rank
    if (rank == 0 .or. rank == 3) Then
      call VOFAdvTHINC_EI(u, phi, nl, dl, dt, 1)
      call VOFAdvTHINC_EI(v, phi, nl, dl, dt, 2)
      call VOFAdvTHINC_EI(w, phi, nl, dl, dt, 3)
    else if (rank == 1 .or. rank == 4) Then
      call VOFAdvTHINC_EI(v, phi, nl, dl, dt, 2)
      call VOFAdvTHINC_EI(w, phi, nl, dl, dt, 3)
      call VOFAdvTHINC_EI(u, phi, nl, dl, dt, 1)
    else if (rank == 2 .or. rank == 5) Then
      call VOFAdvTHINC_EI(w, phi, nl, dl, dt, 3)
      call VOFAdvTHINC_EI(u, phi, nl, dl, dt, 1)
      call VOFAdvTHINC_EI(v, phi, nl, dl, dt, 2)
    endif
  End Subroutine VOFTHINC_EI

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
    Integer  :: i, j, k, ii, jj, kk, ic1, ic2, ic3
    Real(sp) :: a1, a2, alpha
    Real(sp) :: x0(3), deltax(3)
    Real(sp), dimension(0:nl(1) + 1, 0:nl(2) + 1, 0:nl(3) + 1) :: vof1, vof2, vof3
    Real(sp) :: norm(3)
    Real(sp) :: f_block(3, 3, 3)
    ! Real(sp) :: FloodSZ_Backward, FloodSZ_Forward

    ! default, f = 0
    vof1 = 0.0_sp
    vof2 = 0.0_sp
    vof3 = 0.0_sp

    if (dir .eq. 1) then
      ii = 1; ic1 = 1; ic2 = 2; ic3 = 3
    else
      ii = 0
    end if
    if (dir .eq. 2) then
      jj = 1; ic1 = 2; ic2 = 3; ic3 = 1
    else; jj = 0
    end if
    if (dir .eq. 3) then
      kk = 1; ic1 = 3; ic2 = 1; ic3 = 2
    else; kk = 0
    end if

    !***
    do k = 1, nl(3)
      do j = 1, nl(2)
        do i = 1, nl(1)
          a1 = us(i - ii, j - jj, k - kk)*dt/dl(dir)
          a2 = us(i, j, k)*dt/dl(dir)

          ! f = 1
          if (f(i, j, k) .eq. 1.0_sp) then
            vof1(i, j, k) = DMAX1(-a1, 0.0_sp)
            vof2(i, j, k) = 1.0_sp - DMAX1(a1, 0.0_sp) + DMIN1(a2, 0.0_sp)
            vof3(i, j, k) = DMAX1(a2, 0.0_sp)

            ! 0 < f < 1
          else if (f(i, j, k) .GT. 0.0_sp) then
            !*(1)* normal vector
            f_block = f(i - 1:i + 1, j - 1:j + 1, k - 1:k + 1)
            Call VOFNorm(f_block, norm)
            Call Normalization1(norm)
            !*(2) get alpha;
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

    ! apply proper boundary conditions to vof1, vof2, vof3
    call phi_bc%setBCS(vof1)
    call phi_bc%setBCS(vof3)
    !new values of c and  clip it: 0. <= c <= 1.
    do k = 1, nl(3)
      do j = 1, nl(2)
        do i = 1, nl(1)
          f(i, j, k) = vof3(i - ii, j - jj, k - kk) + vof2(i, j, k) + vof1(i + ii, j + jj, k + kk)
          if (f(i, j, k) < EPSC_VOF) then
            f(i, j, k) = 0.0_sp
          elseif (f(i, j, k) > (1.d0 - EPSC_VOF)) then
            f(i, j, k) = 1.0_sp
          endif
        enddo
      enddo
    enddo
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
    Integer  :: i, j, k, ii, jj, kk, ic1, ic2, ic3
    Real(sp) :: a1, a2, alpha
    Real(sp) :: x0(3), deltax(3)
    Real(sp), dimension(0:nl(1) + 1, 0:nl(2) + 1, 0:nl(3) + 1) :: vof1, vof2, vof3
    Real(sp) :: norm(3)
    Real(sp) :: f_block(3, 3, 3)
    ! Real(sp) :: FloodSZ_Backward, FloodSZ_Forward

    ! default, f = 0
    vof1 = 0.0_sp
    vof3 = 0.0_sp

    if (dir .eq. 1) then
      ii = 1; ic1 = 1; ic2 = 2; ic3 = 3
    else
      ii = 0
    end if
    if (dir .eq. 2) then
      jj = 1; ic1 = 2; ic2 = 3; ic3 = 1
    else; jj = 0
    end if
    if (dir .eq. 3) then
      kk = 1; ic1 = 3; ic2 = 1; ic3 = 2
    else; kk = 0
    end if

    !***
    do k = 1, nl(3)
      do j = 1, nl(2)
        do i = 1, nl(1)
          a1 = us(i - ii, j - jj, k - kk)*dt/dl(dir)
          a2 = us(i, j, k)*dt/dl(dir)

          vof1(i, j, k) = 0.d0
          vof3(i, j, k) = 0.d0
          vof2(i, j, k) = 0.d0
          ! f = 1
          if (f(i, j, k) .GT. 1.0_sp - epsc_vof) then
            vof1(i, j, k) = DMAX1(-a1, 0.0_sp)
            vof3(i, j, k) = DMAX1(a2, 0.0_sp)

            ! 0 < f < 1
          else if (f(i, j, k) .GT. epsc_vof) then
            !*(1)* normal vector
            f_block = f(i - 1:i + 1, j - 1:j + 1, k - 1:k + 1)
            Call VOFNorm(f_block, norm)
            Call Normalization1(norm)
            !*(2) get alpha;
            alpha = FloodSZ_Backward(norm, f(i, j, k))
            ! alpha = AL3D(norm,f(i,j,k))
            !*(3) get fluxes
            x0 = 0.0_sp; deltax = 1.0_sp
            if (a1 .LT. 0.0_sp) then
              deltax(dir) = -a1
              vof1(i, j, k) = FloodSZ_Forward(norm, alpha, x0, deltax)
              ! vof1(i,j,k) = FL3D(norm,alpha,x0,deltax)
            end if
            if (a2 .GT. 0.0_sp) then
              x0(dir) = 1.0_sp - a2; 
              deltax(dir) = a2
              vof3(i, j, k) = FloodSZ_Forward(norm, alpha, x0, deltax)
              ! vof3(i,j,k) = FL3D(norm,alpha,x0,deltax)
            end if
          endif
          vof2(i, j, k) = f(i, j, k) - vof1(i, j, k) - vof3(i, j, k)
        enddo
      enddo
    enddo

    ! apply proper boundary conditions to vof1, vof3
    call phi_bc%setBCS(vof1)
    call phi_bc%setBCS(vof3)
    !new values of f and  clip it: 0. <= f <= 1.
    do k = 1, nl(3)
      do j = 1, nl(2)
        do i = 1, nl(1)
          a1 = us(i - ii, j - jj, k - kk)*dt/dl(dir)
          a2 = us(i, j, k)*dt/dl(dir)
          f(i, j, k) = (vof3(i - ii, j - jj, k - kk) + vof2(i, j, k) + vof1(i + ii, j + jj, k + kk))/(1.0_sp - a2 + a1)
          ! f(i,j,k) = ( vof3(i-ii,j-jj,k-kk) + vof2(i,j,k) + vof1(i+ii,j+jj,k+kk) ) * (1.0_sp+a2-a1)
          ! f(i,j,k) = ( vof3(i-ii,j-jj,k-kk) + vof2(i,j,k) + vof1(i+ii,j+jj,k+kk) )  &
          ! + f(i,j,k) * (a2-a1)
          ! f(i,j,k) = ( vof3(i-ii,j-jj,k-kk) + vof2(i,j,k) + vof1(i+ii,j+jj,k+kk) )
          if (f(i, j, k) < EPSC_VOF) then
            f(i, j, k) = 0.0_sp
          elseif (f(i, j, k) > (1.0_sp - EPSC_VOF)) then
            f(i, j, k) = 1.0_sp
          endif
        enddo
      enddo
    enddo
    ! apply proper boundary conditions to c
    call phi_bc%setBCS(f)
    !***
    return
  End Subroutine VOFAdv_EI

!=======================================================
! (2-4) Advection algorithm of CIAM MOF (LE)
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
    Integer  :: i, j, k, ii, jj, kk, ic1, ic2, ic3
    Real(sp) :: a1, a2, alpha
    Real(sp) :: mx, my, mz
    Real(sp) :: x0(3), deltax(3)
    Real(sp) :: c1xyz(3), c2xyz(3), c3xyz(3)
    Real(sp), dimension(0:nl(1) + 1, 0:nl(2) + 1, 0:nl(3) + 1) :: vof1, vof2, vof3
    Real(sp), dimension(0:nl(1) + 1, 0:nl(2) + 1, 0:nl(3) + 1) :: c1x, c1y, c1z
    Real(sp), dimension(0:nl(1) + 1, 0:nl(2) + 1, 0:nl(3) + 1) :: c2x, c2y, c2z
    Real(sp), dimension(0:nl(1) + 1, 0:nl(2) + 1, 0:nl(3) + 1) :: c3x, c3y, c3z
    ! Real(sp) :: f_block(3,3,3), init_norm(3)
    Real(sp) :: norm(3)
    Real(sp) :: c3(3)
    ! Real(sp) :: FloodSZ_Backward

    ! default, f = 0
    vof1 = 0.0_sp
    vof2 = 0.0_sp
    vof3 = 0.0_sp

    if (dir .eq. 1) then
      ii = 1; ic1 = 1; ic2 = 2; ic3 = 3
    else
      ii = 0
    end if
    if (dir .eq. 2) then
      jj = 1; ic1 = 2; ic2 = 3; ic3 = 1
    else; jj = 0
    end if
    if (dir .eq. 3) then
      kk = 1; ic1 = 3; ic2 = 1; ic3 = 2
    else; kk = 0
    end if

    !***
    do k = 1, nl(3)
      do j = 1, nl(2)
        do i = 1, nl(1)
          a1 = us(i - ii, j - jj, k - kk)*dt/dl(dir)
          a2 = us(i, j, k)*dt/dl(dir)

          c1xyz = 0.0_sp
          c2xyz = 0.0_sp
          c3xyz = 0.0_sp
          ! f = 1
          if (f(i, j, k) .GE. 1.0_sp - epsc_vof) then
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
          else if (f(i, j, k) .GT. 0.0_sp) then
            !*(1)* normal vector
            c3(1) = cx(i, j, k)
            c3(2) = cy(i, j, k)
            c3(3) = cz(i, j, k)
            ! f_block = f(i-1:i+1,j-1:j+1,k-1:k+1)
            ! Call NormMYCS(f_block, norm)
            ! Call NormMOF(f(i,j,k), c3, norm, init_norm=init_norm)
            Call MOFNorm(f(i, j, k), c3, norm)
            num_iter = num_iter + mof_niter
            grid_iter = grid_iter + 1
            !*(2) get alpha;
            alpha = FloodSZ_Backward(norm, f(i, j, k))
            norm(dir) = norm(dir)/(1.0_sp - a1 + a2)
            alpha = alpha + norm(dir)*a1
            !*(3) get fluxes
            x0 = 0.0_sp; deltax = 1.0_sp
            if (a1 .LT. 0.0_sp) then
              x0(dir) = a1; 
              deltax(dir) = -a1
              Call FloodSZ_ForwardC(norm, alpha, x0, deltax, vof1(i, j, k), c1xyz)
            end if
            if (a2 .GT. 0.0_sp) then
              x0(dir) = 1.0_sp; 
              deltax(dir) = a2
              Call FloodSZ_ForwardC(norm, alpha, x0, deltax, vof3(i, j, k), c3xyz)
            end if
            x0(dir) = DMAX1(a1, 0.0_sp)
            deltax(dir) = 1.0_sp - x0(dir) + DMIN1(0.0_sp, a2)
            Call FloodSZ_ForwardC(norm, alpha, x0, deltax, vof2(i, j, k), c2xyz)
          endif

          ! c2xyz(dir) = c2xyz(dir) + a1 - a2

          If (vof1(i, j, k) .ge. epsc_vof) Then
            ! c1xyz(dir) = c1xyz(dir) - a1
            ! Call Centroid_Eulerian_Adv(c1xyz, a1, a2, 0.0_sp, 1.0_sp, dir)
            c1xyz(dir) = c1xyz(dir) + 1.0_sp
          End If
          If (vof3(i, j, k) .ge. epsc_vof) Then
            ! c3xyz(dir) = c3xyz(dir) + a2
            ! Call Centroid_Eulerian_Adv(c3xyz, a1, a2, 0.0_sp, 1.0_sp, dir)
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
    !new values of c and  clip it: 0. <= c <= 1.
    do k = 1, nl(3)
      do j = 1, nl(2)
        do i = 1, nl(1)

          mx = vof3(i - ii, j - jj, k - kk)*c3x(i - ii, j - jj, k - kk) + &
               vof2(i, j, k)*c2x(i, j, k) + &
               vof1(i + ii, j + jj, k + kk)*c1x(i + ii, j + jj, k + kk)
          my = vof3(i - ii, j - jj, k - kk)*c3y(i - ii, j - jj, k - kk) + &
               vof2(i, j, k)*c2y(i, j, k) + &
               vof1(i + ii, j + jj, k + kk)*c1y(i + ii, j + jj, k + kk)
          mz = vof3(i - ii, j - jj, k - kk)*c3z(i - ii, j - jj, k - kk) + &
               vof2(i, j, k)*c2z(i, j, k) + &
               vof1(i + ii, j + jj, k + kk)*c1z(i + ii, j + jj, k + kk)

          f(i, j, k) = vof3(i - ii, j - jj, k - kk) + &
                       vof2(i, j, k) + &
                       vof1(i + ii, j + jj, k + kk)
          cx(i, j, k) = mx/(f(i, j, k) + 1e-30)
          cy(i, j, k) = my/(f(i, j, k) + 1e-30)
          cz(i, j, k) = mz/(f(i, j, k) + 1e-30)
          if (f(i, j, k) < EPSC_VOF) then
            f(i, j, k) = 0.0_sp
            cx(i, j, k) = 0.0_sp
            cy(i, j, k) = 0.0_sp
            cz(i, j, k) = 0.0_sp
          elseif (f(i, j, k) > (1.d0 - EPSC_VOF)) then
            f(i, j, k) = 1.0_sp
            cx(i, j, k) = 0.5_sp
            cy(i, j, k) = 0.5_sp
            cz(i, j, k) = 0.5_sp
          endif
          if (cx(i, j, k) < 0.0_sp) then
            cx(i, j, k) = 0.0_sp
            print *, 'cx<1', i, j, k
          elseif (cx(i, j, k) > (1.0_sp)) then
            cx(i, j, k) = 1.0_sp
            print *, 'cx>1', i, j, k
          endif
          if (cy(i, j, k) < 0.0_sp) then
            cy(i, j, k) = 0.0_sp
            print *, 'cy<1', i, j, k
          elseif (cy(i, j, k) > (1.0_sp)) then
            cy(i, j, k) = 1.0_sp
            print *, 'cy>1', i, j, k
          endif
          if (cz(i, j, k) < 0.0_sp) then
            cz(i, j, k) = 0.0_sp
            print *, 'cz<0', i, j, k
          elseif (cz(i, j, k) > (1.0_sp)) then
            cz(i, j, k) = 1.0_sp
            print *, 'cz>1', i, j, k
          endif
        enddo
      enddo
    enddo

    ! apply proper boundary conditions to c
    call phi_bc%setBCS(f)
    call phi_bc%setBCS(cx)
    call phi_bc%setBCS(cy)
    call phi_bc%setBCS(cz)
    !***
    return
  End Subroutine MOFAdv_LE

!=======================================================
! (2-5) Advection algorithm of CIAM MOF (LE)
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
    Integer  :: i, j, k, ii, jj, kk, ic1, ic2, ic3
    Real(sp) :: a1, a2, alpha
    Real(sp) :: mx, my, mz
    Real(sp) :: x0(3), deltax(3)
    Real(sp) :: cxyz(3), c1xyz(3), c2xyz(3), c3xyz(3)
    Real(sp), dimension(0:nl(1) + 1, 0:nl(2) + 1, 0:nl(3) + 1) :: vof1, vof2, vof3
    Real(sp), dimension(0:nl(1) + 1, 0:nl(2) + 1, 0:nl(3) + 1) :: c1x, c1y, c1z
    Real(sp), dimension(0:nl(1) + 1, 0:nl(2) + 1, 0:nl(3) + 1) :: c2x, c2y, c2z
    Real(sp), dimension(0:nl(1) + 1, 0:nl(2) + 1, 0:nl(3) + 1) :: c3x, c3y, c3z
    Real(sp) :: norm(3)
    Real(sp) :: f_block(3, 3, 3)
    Real(sp) :: exp_factor
    ! Real(sp) :: FloodSZ_Backward, FloodSZ_Forward

    ! default, f = 0
    vof1 = 0.0_sp
    vof3 = 0.0_sp

    if (dir .eq. 1) then
      ii = 1; ic1 = 1; ic2 = 2; ic3 = 3
    else
      ii = 0
    end if
    if (dir .eq. 2) then
      jj = 1; ic1 = 2; ic2 = 3; ic3 = 1
    else; jj = 0
    end if
    if (dir .eq. 3) then
      kk = 1; ic1 = 3; ic2 = 1; ic3 = 2
    else; kk = 0
    end if

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
          if (f(i, j, k) .GT. 1.0_sp - epsc_vof) then
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
          else if (f(i, j, k) .GT. epsc_vof) then
            !*(1)* normal vector
            f_block = f(i - 1:i + 1, j - 1:j + 1, k - 1:k + 1)
            cxyz(1) = cx(i, j, k)
            cxyz(2) = cy(i, j, k)
            cxyz(3) = cz(i, j, k)
            Call VOFNorm(f_block, norm)
            Call MOFNorm(f(i, j, k), cxyz, norm)
            Call Normalization1(norm)
            !*(2) get alpha;
            alpha = FloodSZ_Backward(norm, f(i, j, k))
            ! alpha = AL3D(norm,f(i,j,k))
            !*(3) get fluxes
            x0 = 0.0_sp; deltax = 1.0_sp
            if (a1 .LT. 0.0_sp) then
              deltax(dir) = -a1
              Call FloodSZ_ForwardC(norm, alpha, x0, deltax, vof1(i, j, k), c1xyz)
              ! vof1(i,j,k) = FL3D(norm,alpha,x0,deltax)
            end if
            if (a2 .GT. 0.0_sp) then
              x0(dir) = 1.0_sp - a2; 
              deltax(dir) = a2
              Call FloodSZ_ForwardC(norm, alpha, x0, deltax, vof3(i, j, k), c3xyz)
              ! vof3(i,j,k) = FL3D(norm,alpha,x0,deltax)
            end if
            x0(dir) = DMAX1(-a1, 0.0_sp)
            deltax(dir) = 1.0_sp - x0(dir) - DMAX1(0.0_sp, a2)
            Call FloodSZ_ForwardC(norm, alpha, x0, deltax, vof2(i, j, k), c2xyz)
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

    ! apply proper boundary conditions to vof1, vof3
    call phi_bc%setBCS(vof1)
    call phi_bc%setBCS(vof3)
    call phi_bc%setBCS(c1x)
    call phi_bc%setBCS(c3x)
    call phi_bc%setBCS(c1y)
    call phi_bc%setBCS(c3y)
    call phi_bc%setBCS(c1z)
    call phi_bc%setBCS(c3z)
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
            ! f(i,j,k) = min(max(f(i,j,k),0.0_sp),1.0_sp)
          endif
        enddo
      enddo
    enddo
    ! apply proper boundary conditions to c
    call phi_bc%setBCS(f)
    !***
    return
  End Subroutine MOFAdv_EI

!=======================================================
! (2-3) Advection algorithm of THINC (EI)
!=======================================================
  Subroutine VOFAdvTHINC_EI(us, f, nl, dl, dt, dir)
    ! Use ModVOFFunc
    Implicit None
    Integer :: dir
    Real(sp) :: dt
    Integer, Intent(In)     :: nl(3)
    Real(sp), Intent(In)    :: dl(3)
    Real(sp), Intent(In)    :: us(0:nl(1) + 1, 0:nl(2) + 1, 0:nl(3) + 1)
    Real(sp), Intent(InOut) :: f(0:nl(1) + 1, 0:nl(2) + 1, 0:nl(3) + 1)
    Integer  :: i, j, k, ii, jj, kk, ic1, ic2, ic3
    Real(sp) :: a1, a2, beta, gamma, betagamma2, x_center
    Real(sp) :: x0(3), deltax(3)
    Real(sp), dimension(0:nl(1) + 1, 0:nl(2) + 1, 0:nl(3) + 1) :: vof1, vof2, vof3
    Real(sp) :: norm(3)
    Real(sp) :: f_block(3, 3, 3)
    ! Real(sp) :: FloodSZ_Backward, FloodSZ_Forward

    ! default, f = 0
    vof1 = 0.0_sp
    vof2 = 0.0_sp
    vof3 = 0.0_sp

    if (dir .eq. 1) then
      ii = 1; ic1 = 1; ic2 = 2; ic3 = 3
    else
      ii = 0
    end if
    if (dir .eq. 2) then
      jj = 1; ic1 = 2; ic2 = 3; ic3 = 1
    else; jj = 0
    end if
    if (dir .eq. 3) then
      kk = 1; ic1 = 3; ic2 = 1; ic3 = 2
    else; kk = 0
    end if

    !***
    do k = 1, nl(3)
      do j = 1, nl(2)
        do i = 1, nl(1)
          a1 = us(i - ii, j - jj, k - kk)*dt/dl(dir)
          a2 = us(i, j, k)*dt/dl(dir)

          ! f = 1
          if (f(i, j, k) .GE. 1.0_sp - epsc_vof) then
            vof1(i, j, k) = DMAX1(-a1, 0.0_sp)
            vof3(i, j, k) = DMAX1(a2, 0.0_sp)

            ! 0 < f < 1
          else if (f(i, j, k) .GT. 0.0_sp) then
            !*(1)* normal vector
            f_block = f(i - 1:i + 1, j - 1:j + 1, k - 1:k + 1)
            Call VOFNorm(f_block, norm)
            Call Normalization1(norm)
            !*(2) determine gamma
            if (norm(dir) .gt. 0) Then
              gamma = -1.0_sp
            else
              gamma = 1.0_sp
            endif
            beta = 3.5_sp*dabs(norm(dir)) + 0.01_sp
            ! Use the formulation of THINC/SW
            betagamma2 = 2.0_sp*beta*gamma
            x_center = THINC1DBackward(f(i, j, k), betagamma2)
            x0 = 0.0_sp; deltax = 1.0_sp
            if (a1 .LT. 0.0_sp) then
              x0(dir) = 0; 
              deltax(dir) = -a1
              vof1(i, j, k) = THINC1DForward(x_center, betagamma2, x0(dir), deltax(dir))
            end if
            if (a2 .GT. 0.0_sp) then
              x0(dir) = 1.0_sp - a2; 
              deltax(dir) = a2
              vof3(i, j, k) = THINC1DForward(x_center, betagamma2, x0(dir), deltax(dir))
            end if
          endif
          vof2(i, j, k) = f(i, j, k) - vof1(i, j, k) - vof3(i, j, k)
        enddo
      enddo
    enddo

    ! apply proper boundary conditions to vof1, vof2, vof3
    call phi_bc%setBCS(vof1)
    call phi_bc%setBCS(vof3)
    !new values of c and  clip it: 0. <= c <= 1.
    do k = 1, nl(3)
      do j = 1, nl(2)
        do i = 1, nl(1)
          a1 = us(i - ii, j - jj, k - kk)*dt/dl(dir)
          a2 = us(i, j, k)*dt/dl(dir)
          f(i, j, k) = (vof3(i - ii, j - jj, k - kk) + vof2(i, j, k) + vof1(i + ii, j + jj, k + kk))/(1.0_sp - a2 + a1)
          if (f(i, j, k) < EPSC_VOF) then
            f(i, j, k) = 0.0_sp
          elseif (f(i, j, k) > (1.d0 - EPSC_VOF)) then
            f(i, j, k) = 1.0_sp
          endif
        enddo
      enddo
    enddo
    ! apply proper boundary conditions to c
    call phi_bc%setBCS(f)
    !***
    return
  End Subroutine VOFAdvTHINC_EI



! Subroutine Centroid_Lagrangian_Adv(c, ul, ur, xl, xr, dir)
!   Implicit None
!   Real(sp), Intent(InOut) :: c(3)
!   Real(sp), Intent(In)    :: ul, ur
!   Real(sp), Intent(In)    :: xl, xr
!   Integer,  Intent(In)    :: dir

!   c(dir) = (1.0_sp + ur - ul) * c(dir) - (ur * xl - ul * xr)

! End Subroutine Centroid_Lagrangian_Adv

! Subroutine Centroid_Eulerian_Adv(c, ul, ur, xl, xr, dir)
!   Implicit None
!   Real(sp), Intent(InOut) :: c(3)
!   Real(sp), Intent(In)    :: ul, ur
!   Real(sp), Intent(In)    :: xl, xr
!   Integer,  Intent(In)    :: dir

!   c(dir) =  (c(dir) - (ur * xl - ul * xr)) / (1.0_sp - ur + ul)

! End Subroutine Centroid_Eulerian_Adv

End Module ModVOF
