Module ModVOFExt
  Use ModVOFFuncExt
  Use ModVOF

  ! Originally, the names of the subroutines are
  ! AdvCIAM, AdvWY ...
  !! This part allows the old cases to call the modified subroutine name
  ! May delete it later

  Interface
    Subroutine InterfaceVOFTIME(Phi, u, v, w, nl, dl, dt, rank,Phi0)
      import
      Implicit None
      Real(sp) :: dt
      Integer, Intent(In)            :: nl(3)
      Real(sp), Intent(In)           :: dl(3)
      Real(sp), Intent(InOut)        ::  Phi(0:nl(1) + 1, 0:nl(2) + 1, 0:nl(3) + 1)
      Real(sp), Intent(In), optional :: Phi0(0:nl(1) + 1, 0:nl(2) + 1, 0:nl(3) + 1)
      Real(sp), Intent(In)           ::    u(0:nl(1) + 1, 0:nl(2) + 1, 0:nl(3) + 1)
      Real(sp), Intent(In)           ::    v(0:nl(1) + 1, 0:nl(2) + 1, 0:nl(3) + 1)
      Real(sp), Intent(In)           ::    w(0:nl(1) + 1, 0:nl(2) + 1, 0:nl(3) + 1)
      Integer :: rank
    End Subroutine InterfaceVOFTIME

    Subroutine InterfaceMOFTIME(Phi, u, v, w, nl, dl, dt, rank, cx, cy, cz, phi0)
      import
      Implicit None
      Real(sp) :: dt
      Integer, Intent(In)     :: nl(3)
      Real(sp), Intent(In)    :: dl(3)
      Real(sp), Intent(InOut)        ::  Phi(0:nl(1) + 1, 0:nl(2) + 1, 0:nl(3) + 1)
      Real(sp), Intent(In), Optional :: Phi0(0:nl(1) + 1, 0:nl(2) + 1, 0:nl(3) + 1)
      Real(sp), Intent(InOut)        ::   cx(0:nl(1) + 1, 0:nl(2) + 1, 0:nl(3) + 1)
      Real(sp), Intent(InOut)        ::   cy(0:nl(1) + 1, 0:nl(2) + 1, 0:nl(3) + 1)
      Real(sp), Intent(InOut)        ::   cz(0:nl(1) + 1, 0:nl(2) + 1, 0:nl(3) + 1)
      Real(sp), Intent(In)           ::    u(0:nl(1) + 1, 0:nl(2) + 1, 0:nl(3) + 1)
      Real(sp), Intent(In)           ::    v(0:nl(1) + 1, 0:nl(2) + 1, 0:nl(3) + 1)
      Real(sp), Intent(In)           ::    w(0:nl(1) + 1, 0:nl(2) + 1, 0:nl(3) + 1)
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

    Subroutine InterfaceTHINC_Beta(norm, f, dir, beta)
      import
      Real(SP), Intent(In) :: norm(3)
      Real(SP), Intent(In) :: f
      Integer , Intent(In) :: dir
      Real(SP), Intent(Out) :: beta
    End Subroutine InterfaceTHINC_Beta

  End Interface

  PROCEDURE(InterfaceMOFTIME), POINTER :: MOFCIAM => MOF_LE
  PROCEDURE(InterfaceMOFTIME), POINTER :: MOFWY => MOF_EI
  PROCEDURE(InterfaceVOFTIME), POINTER :: VOFCIAM => VOF_LE
  PROCEDURE(InterfaceVOFTIME), POINTER :: VOFWY => VOF_EI
  PROCEDURE(InterfaceMOFDIR),  POINTER :: AdvCIAM_MOF => MOFAdv_LE
  PROCEDURE(InterfaceMOFDIR),  POINTER :: AdvWY_MOF => MOFAdv_EI
  PROCEDURE(InterfaceVOFDIR),  POINTER :: AdvCIAM => VOFAdv_LE
  PROCEDURE(InterfaceVOFDIR),  POINTER :: AdvWY => VOFAdv_EI
  PROCEDURE(InterfaceTHINC_Beta),  POINTER :: THINC_Beta => THINC_Beta_SL

Contains 

  !=======================================================
  ! [1-2] VOF-LE (Lagrangian Explicit)
  !-------------------------------------------------------
  ! Sweep sequence:
  !   1-2-3
  !   2-3-1
  !   3-1-2
  ! All: LE-LE-LE
  !=======================================================
  Subroutine VOF_LE(Phi, u, v, w, nl, dl, dt, rank, Phi0)
    Implicit None
    Real(sp) :: dt
    Integer, Intent(In)            :: nl(3)
    Real(sp), Intent(In)           :: dl(3)
    Real(sp), Intent(InOut)        ::  Phi(0:nl(1) + 1, 0:nl(2) + 1, 0:nl(3) + 1)
    Real(sp), Intent(In), optional :: Phi0(0:nl(1) + 1, 0:nl(2) + 1, 0:nl(3) + 1)
    Real(sp), Intent(In)           ::    u(0:nl(1) + 1, 0:nl(2) + 1, 0:nl(3) + 1)
    Real(sp), Intent(In)           ::    v(0:nl(1) + 1, 0:nl(2) + 1, 0:nl(3) + 1)
    Real(sp), Intent(In)           ::    w(0:nl(1) + 1, 0:nl(2) + 1, 0:nl(3) + 1)
    Integer :: rank
    num_iter = 0
    grid_iter = 0
    t_rec = 0.0_sp
    t_adv = 0.0_sp
    t_cor = 0.0_sp
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
  ! [1-3] VOF-EI (Eulerian Implicit)
  !-------------------------------------
  ! Sweep sequence:
  !   1-2-3
  !   2-3-1
  !   3-1-2
  ! All: EI-EI-EI
  !=======================================================
  Subroutine VOF_EI(Phi, u, v, w, nl, dl, dt, rank, Phi0)
    Implicit None
    Real(sp) :: dt
    Integer, Intent(In)            :: nl(3)
    Real(sp), Intent(In)           :: dl(3)
    Real(sp), Intent(InOut)        ::  Phi(0:nl(1) + 1, 0:nl(2) + 1, 0:nl(3) + 1)
    Real(sp), Intent(In), optional :: Phi0(0:nl(1) + 1, 0:nl(2) + 1, 0:nl(3) + 1)
    Real(sp), Intent(In)           ::    u(0:nl(1) + 1, 0:nl(2) + 1, 0:nl(3) + 1)
    Real(sp), Intent(In)           ::    v(0:nl(1) + 1, 0:nl(2) + 1, 0:nl(3) + 1)
    Real(sp), Intent(In)           ::    w(0:nl(1) + 1, 0:nl(2) + 1, 0:nl(3) + 1)
    Integer :: rank
    num_iter = 0
    grid_iter = 0
    t_rec = 0.0_sp
    t_adv = 0.0_sp
    t_cor = 0.0_sp
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
  ! [1-4-1] VOF-EILE3DS (Eulerian Implicit Lanrangian Explicit)
  !-------------------------------------------------------
  ! Sweep sequence:
  !   1-2-3 (m1)
  !   1-2-3 (m2)
  !   2-3-1 (m1)
  !   2-3-1 (m2)
  !   3-1-2 (m1)
  !   3-1-2 (m2)
  ! m1: EI-LE-EI
  ! m2: LE-EI-LE
  !=======================================================
  Subroutine VOF_EILE3DS(Phi, u, v, w, nl, dl, dt, rank, Phi0)
    Implicit None
    Real(sp) :: dt
    Integer, Intent(In)            :: nl(3)
    Real(sp), Intent(In)           :: dl(3)
    Real(sp), Intent(InOut)        ::  Phi(0:nl(1) + 1, 0:nl(2) + 1, 0:nl(3) + 1)
    Real(sp), Intent(In), optional :: Phi0(0:nl(1) + 1, 0:nl(2) + 1, 0:nl(3) + 1)
    Real(sp), Intent(In)           ::    u(0:nl(1) + 1, 0:nl(2) + 1, 0:nl(3) + 1)
    Real(sp), Intent(In)           ::    v(0:nl(1) + 1, 0:nl(2) + 1, 0:nl(3) + 1)
    Real(sp), Intent(In)           ::    w(0:nl(1) + 1, 0:nl(2) + 1, 0:nl(3) + 1)
    Real(sp)       :: u1(0:nl(1) + 1, 0:nl(2) + 1, 0:nl(3) + 1)
    Real(sp)       :: v1(0:nl(1) + 1, 0:nl(2) + 1, 0:nl(3) + 1)
    Real(sp)       :: w1(0:nl(1) + 1, 0:nl(2) + 1, 0:nl(3) + 1)
    Integer :: rank
    num_iter = 0
    grid_iter = 0
    t_rec = 0.0_sp
    t_adv = 0.0_sp
    t_cor = 0.0_sp
    If (rank == 0) Then
      call VOFAdv_EI(u, phi, nl, dl, dt, 1)
      call VOFAdv_LE(v, phi, nl, dl, dt, 2)
      call VOFAdv_EI(w, phi, nl, dl, dt, 3)
    Else If (rank == 1) Then
      call VOFAdv_LE(v, phi, nl, dl, dt, 2)
      call VOFAdv_EI(u, phi, nl, dl, dt, 3)
      call VOFAdv_LE(w, phi, nl, dl, dt, 1)
    Else If (rank == 2) Then
      call VOFAdv_EI(w, phi, nl, dl, dt, 3)
      call VOFAdv_LE(u, phi, nl, dl, dt, 1)
      call VOFAdv_EI(v, phi, nl, dl, dt, 2)
    Else If (rank == 3) Then
      call VOFAdv_LE(u, phi, nl, dl, dt, 1)
      call VOFAdv_EI(v, phi, nl, dl, dt, 2)
      call VOFAdv_LE(w, phi, nl, dl, dt, 3)
    Else If (rank == 4) Then
      call VOFAdv_EI(v, phi, nl, dl, dt, 2)
      call VOFAdv_LE(w, phi, nl, dl, dt, 3)
      call VOFAdv_EI(u, phi, nl, dl, dt, 1)
    Else If (rank == 5) Then
      call VOFAdv_LE(w, phi, nl, dl, dt, 3)
      call VOFAdv_EI(u, phi, nl, dl, dt, 1)
      call VOFAdv_LE(v, phi, nl, dl, dt, 2)
    EndIf

  End Subroutine VOF_EILE3DS

  !=======================================================
  ! [1-4-2] VOF-EILE3DS (Eulerian Implicit Lanrangian Explicit)
  !-------------------------------------------------------
  ! Sweep sequence:
  !   1-2-3-2-3-1
  !   2-3-1-1-3-2
  !   3-1-2-1-2-3
  ! All: EI-LE-EI-LE-EI-LE
  !=======================================================
  Subroutine VOF_EILE3D(Phi, u, v, w, nl, dl, dt, rank, Phi0)
    Implicit None
    Real(sp) :: dt
    Integer, Intent(In)            :: nl(3)
    Real(sp), Intent(In)           :: dl(3)
    Real(sp), Intent(InOut)        ::  Phi(0:nl(1) + 1, 0:nl(2) + 1, 0:nl(3) + 1)
    Real(sp), Intent(In), optional :: Phi0(0:nl(1) + 1, 0:nl(2) + 1, 0:nl(3) + 1)
    Real(sp), Intent(In)           ::    u(0:nl(1) + 1, 0:nl(2) + 1, 0:nl(3) + 1)
    Real(sp), Intent(In)           ::    v(0:nl(1) + 1, 0:nl(2) + 1, 0:nl(3) + 1)
    Real(sp), Intent(In)           ::    w(0:nl(1) + 1, 0:nl(2) + 1, 0:nl(3) + 1)
    Real(sp), Dimension(0:nl(1) + 1, 0:nl(2) + 1, 0:nl(3) + 1) :: u1, u2, v1, v3, w2, w3
    Integer :: rank
    num_iter = 0
    grid_iter = 0
    t_rec = 0.0_sp
    t_adv = 0.0_sp
    t_cor = 0.0_sp

    Call Vel_Decomp(u,v,w,u1,v1,w2,u2,v3,w3,nl)
    If (rank == 0 .or. rank == 3) Then
      Call VOFAdv_EI(u1, phi, nl, dl, dt, 1)
      Call VOFAdv_LE(v1, phi, nl, dl, dt, 2)
      Call VOFAdv_EI(w2, phi, nl, dl, dt, 3)
      Call VOFAdv_LE(u2, phi, nl, dl, dt, 1)
      Call VOFAdv_EI(v3, phi, nl, dl, dt, 2)
      Call VOFAdv_LE(w3, phi, nl, dl, dt, 3)
    Else If (rank == 1 .or. rank == 4) Then
      Call VOFAdv_EI(w2, phi, nl, dl, dt, 3)
      Call VOFAdv_LE(u2, phi, nl, dl, dt, 1)
      Call VOFAdv_EI(v3, phi, nl, dl, dt, 2)
      Call VOFAdv_LE(w3, phi, nl, dl, dt, 3)
      Call VOFAdv_EI(u1, phi, nl, dl, dt, 1)
      Call VOFAdv_LE(v1, phi, nl, dl, dt, 2)
    Else If (rank == 2 .or. rank == 5) Then
      Call VOFAdv_EI(v3, phi, nl, dl, dt, 2)
      Call VOFAdv_LE(w3, phi, nl, dl, dt, 3)
      Call VOFAdv_EI(u1, phi, nl, dl, dt, 1)
      Call VOFAdv_LE(v1, phi, nl, dl, dt, 2)
      Call VOFAdv_EI(w2, phi, nl, dl, dt, 3)
      Call VOFAdv_LE(u2, phi, nl, dl, dt, 1)
    End If

  End Subroutine VOF_EILE3D

  !=======================================================
  ! [1-4-3] VOF-EIEALE 
  !     (Eulerian Implicit - 
  !      Eulerian Algebaric -
  !      Lanrangian Explicit)
  !-------------------------------------------------------
  ! Sweep sequence:
  !   1-2-3
  !   2-3-1
  !   3-1-2
  ! All: EI-EA-LE
  !=======================================================
  Subroutine VOF_EIEALE(Phi, u, v, w, nl, dl, dt, rank, Phi0)
    Implicit None
    Real(sp) :: dt
    Integer, Intent(In)            :: nl(3)
    Real(sp), Intent(In)           :: dl(3)
    Real(sp), Intent(InOut)        ::  Phi(0:nl(1) + 1, 0:nl(2) + 1, 0:nl(3) + 1)
    Real(sp), Intent(In), optional :: Phi0(0:nl(1) + 1, 0:nl(2) + 1, 0:nl(3) + 1)
    Real(sp), Intent(In)           ::    u(0:nl(1) + 1, 0:nl(2) + 1, 0:nl(3) + 1)
    Real(sp), Intent(In)           ::    v(0:nl(1) + 1, 0:nl(2) + 1, 0:nl(3) + 1)
    Real(sp), Intent(In)           ::    w(0:nl(1) + 1, 0:nl(2) + 1, 0:nl(3) + 1)
    Integer :: rank
    num_iter = 0
    grid_iter = 0
    t_rec = 0.0_sp
    t_adv = 0.0_sp
    t_cor = 0.0_sp
    if (rank == 0 .or. rank == 3) Then
      call VOFAdv_EI(u, phi, nl, dl, dt, 1)
      call VOFAdv_EA(v, u, w, phi, nl, dl, dt, 2)
      call VOFAdv_LE(w, phi, nl, dl, dt, 3)
    else if (rank == 1 .or. rank == 4) Then
      call VOFAdv_EI(v, phi, nl, dl, dt, 2)
      call VOFAdv_EA(w, v, u, phi, nl, dl, dt, 3)
      call VOFAdv_LE(u, phi, nl, dl, dt, 1)
    else if (rank == 2 .or. rank == 5) Then
      call VOFAdv_EI(w, phi, nl, dl, dt, 3)
      call VOFAdv_EA(u, w, v, phi, nl, dl, dt, 1)
      call VOFAdv_LE(v, phi, nl, dl, dt, 2)
    endif
  End Subroutine VOF_EIEALE

  !=======================================================
  ! [1-4-4] VOF-EILE2D
  !     (Eulerian Implicit - 
  !      Lanrangian Explicit)
  !-------------------------------------------------------
  ! Sweep sequence:
  !   1-2
  !   2-1
  !=======================================================
  Subroutine VOF_EILE2D(Phi, u, v, w, nl, dl, dt, rank, Phi0)
    Implicit None
    Real(sp) :: dt
    Integer, Intent(In)            :: nl(3)
    Real(sp), Intent(In)           :: dl(3)
    Real(sp), Intent(InOut)        ::  Phi(0:nl(1) + 1, 0:nl(2) + 1, 0:nl(3) + 1)
    Real(sp), Intent(In), optional :: Phi0(0:nl(1) + 1, 0:nl(2) + 1, 0:nl(3) + 1)
    Real(sp), Intent(In)           ::    u(0:nl(1) + 1, 0:nl(2) + 1, 0:nl(3) + 1)
    Real(sp), Intent(In)           ::    v(0:nl(1) + 1, 0:nl(2) + 1, 0:nl(3) + 1)
    Real(sp), Intent(In)           ::    w(0:nl(1) + 1, 0:nl(2) + 1, 0:nl(3) + 1)
    Integer :: rank
    num_iter = 0
    grid_iter = 0
    t_rec = 0.0_sp
    t_adv = 0.0_sp
    t_cor = 0.0_sp
    if (rank == 0 .or. rank == 2 .or. rank == 4) Then
      call VOFAdv_EI(u, phi, nl, dl, dt, 1)
      call VOFAdv_LE(v, phi, nl, dl, dt, 2)
    else if (rank == 1 .or. rank == 3 .or. rank == 5) Then
      call VOFAdv_EI(v, phi, nl, dl, dt, 2)
      call VOFAdv_LE(u, phi, nl, dl, dt, 1)
    endif
  End Subroutine VOF_EILE2D

  !=======================================================
  ! [1-6] MOF-LE (Lanrangian Explicit)
  !-------------------------------------------------------
  ! Sweep sequence:
  !   1-2-3
  !   2-3-1
  !   3-1-2
  ! All: LE-LE-LE
  !=======================================================
  Subroutine MOF_LE(Phi, u, v, w, nl, dl, dt, rank, cx, cy, cz, Phi0)
    Implicit None
    Real(sp) :: dt
    Integer, Intent(In)     :: nl(3)
    Real(sp), Intent(In)    :: dl(3)
    Real(sp), Intent(InOut)        ::  Phi(0:nl(1) + 1, 0:nl(2) + 1, 0:nl(3) + 1)
    Real(sp), Intent(In), Optional :: Phi0(0:nl(1) + 1, 0:nl(2) + 1, 0:nl(3) + 1)
    Real(sp), Intent(InOut)        ::   cx(0:nl(1) + 1, 0:nl(2) + 1, 0:nl(3) + 1)
    Real(sp), Intent(InOut)        ::   cy(0:nl(1) + 1, 0:nl(2) + 1, 0:nl(3) + 1)
    Real(sp), Intent(InOut)        ::   cz(0:nl(1) + 1, 0:nl(2) + 1, 0:nl(3) + 1)
    Real(sp), Intent(In)           ::    u(0:nl(1) + 1, 0:nl(2) + 1, 0:nl(3) + 1)
    Real(sp), Intent(In)           ::    v(0:nl(1) + 1, 0:nl(2) + 1, 0:nl(3) + 1)
    Real(sp), Intent(In)           ::    w(0:nl(1) + 1, 0:nl(2) + 1, 0:nl(3) + 1)
    Integer :: rank
    num_iter = 0
    grid_iter = 0
    t_rec = 0.0_sp
    t_adv = 0.0_sp
    t_cor = 0.0_sp
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
  ! [1-7] MOF-EI (Eulurian Implicit)
  !-------------------------------------------------------
  ! Sweep sequence:
  !   1-2-3
  !   2-3-1
  !   3-1-2
  ! All: LE-LE-LE
  !=======================================================
  Subroutine MOF_EI(Phi, u, v, w, nl, dl, dt, rank, cx, cy, cz, Phi0)
    Implicit None
    Real(sp) :: dt
    Integer, Intent(In)     :: nl(3)
    Real(sp), Intent(In)    :: dl(3)
    Real(sp), Intent(InOut)        ::  Phi(0:nl(1) + 1, 0:nl(2) + 1, 0:nl(3) + 1)
    Real(sp), Intent(In), Optional :: Phi0(0:nl(1) + 1, 0:nl(2) + 1, 0:nl(3) + 1)
    Real(sp), Intent(InOut)        ::   cx(0:nl(1) + 1, 0:nl(2) + 1, 0:nl(3) + 1)
    Real(sp), Intent(InOut)        ::   cy(0:nl(1) + 1, 0:nl(2) + 1, 0:nl(3) + 1)
    Real(sp), Intent(InOut)        ::   cz(0:nl(1) + 1, 0:nl(2) + 1, 0:nl(3) + 1)
    Real(sp), Intent(In)           ::    u(0:nl(1) + 1, 0:nl(2) + 1, 0:nl(3) + 1)
    Real(sp), Intent(In)           ::    v(0:nl(1) + 1, 0:nl(2) + 1, 0:nl(3) + 1)
    Real(sp), Intent(In)           ::    w(0:nl(1) + 1, 0:nl(2) + 1, 0:nl(3) + 1)
    Integer :: rank
    num_iter = 0
    grid_iter = 0
    t_rec = 0.0_sp
    t_adv = 0.0_sp
    t_cor = 0.0_sp
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
  ! [1-8-1] MOF-EILE3DS (Eulerian Implicit Lanrangian Explicit)
  !-------------------------------------------------------
  ! Sweep sequence:
  !   1-2-3 (m1)
  !   1-2-3 (m2)
  !   2-3-1 (m1)
  !   2-3-1 (m2)
  !   3-1-2 (m1)
  !   3-1-2 (m2)
  ! m1: EI-LE-EI
  ! m2: LE-EI-LE
  !=======================================================
  Subroutine MOF_EILE3DS(Phi, u, v, w, nl, dl, dt, rank, cx, cy, cz, Phi0)
    Implicit None
    Real(sp) :: dt
    Integer, Intent(In)     :: nl(3)
    Real(sp), Intent(In)    :: dl(3)
    Real(sp), Intent(InOut)        ::  Phi(0:nl(1) + 1, 0:nl(2) + 1, 0:nl(3) + 1)
    Real(sp), Intent(In), Optional :: Phi0(0:nl(1) + 1, 0:nl(2) + 1, 0:nl(3) + 1)
    Real(sp), Intent(InOut)        ::   cx(0:nl(1) + 1, 0:nl(2) + 1, 0:nl(3) + 1)
    Real(sp), Intent(InOut)        ::   cy(0:nl(1) + 1, 0:nl(2) + 1, 0:nl(3) + 1)
    Real(sp), Intent(InOut)        ::   cz(0:nl(1) + 1, 0:nl(2) + 1, 0:nl(3) + 1)
    Real(sp), Intent(In)           ::    u(0:nl(1) + 1, 0:nl(2) + 1, 0:nl(3) + 1)
    Real(sp), Intent(In)           ::    v(0:nl(1) + 1, 0:nl(2) + 1, 0:nl(3) + 1)
    Real(sp), Intent(In)           ::    w(0:nl(1) + 1, 0:nl(2) + 1, 0:nl(3) + 1)
    Integer :: rank
    num_iter = 0
    grid_iter = 0
    t_rec = 0.0_sp
    t_adv = 0.0_sp
    t_cor = 0.0_sp
    If (rank == 0) Then
      Call MOFAdv_EI(u, cx, cy, cz, phi, nl, dl, dt, 1)
      Call MOFAdv_LE(v, cx, cy, cz, phi, nl, dl, dt, 2)
      Call MOFAdv_EI(w, cx, cy, cz, phi, nl, dl, dt, 3)
    Else If (rank == 1) Then
      Call MOFAdv_LE(u, cx, cy, cz, phi, nl, dl, dt, 1)
      Call MOFAdv_EI(v, cx, cy, cz, phi, nl, dl, dt, 2)
      Call MOFAdv_LE(w, cx, cy, cz, phi, nl, dl, dt, 3)
    Else If (rank == 2) Then
      Call MOFAdv_EI(v, cx, cy, cz, phi, nl, dl, dt, 2)
      Call MOFAdv_LE(w, cx, cy, cz, phi, nl, dl, dt, 3)
      Call MOFAdv_EI(u, cx, cy, cz, phi, nl, dl, dt, 1)
    Else If (rank == 3) Then
      Call MOFAdv_LE(v, cx, cy, cz, phi, nl, dl, dt, 2)
      Call MOFAdv_EI(w, cx, cy, cz, phi, nl, dl, dt, 3)
      Call MOFAdv_LE(u, cx, cy, cz, phi, nl, dl, dt, 1)
    Else If (rank == 4) Then
      Call MOFAdv_EI(w, cx, cy, cz, phi, nl, dl, dt, 3)
      Call MOFAdv_LE(u, cx, cy, cz, phi, nl, dl, dt, 1)
      Call MOFAdv_EI(v, cx, cy, cz, phi, nl, dl, dt, 2)
    Else If (rank == 5) Then
      Call MOFAdv_LE(w, cx, cy, cz, phi, nl, dl, dt, 3)
      Call MOFAdv_EI(u, cx, cy, cz, phi, nl, dl, dt, 1)
      Call MOFAdv_LE(v, cx, cy, cz, phi, nl, dl, dt, 2)
    EndIf

  End Subroutine MOF_EILE3DS

  !=======================================================
  ! [1-8-2] MOF-EILE3DS (Eulerian Implicit Lanrangian Explicit)
  !-------------------------------------------------------
  ! Sweep sequence:
  !   1-2-3 (m1)
  !   1-2-3 (m2)
  !   2-3-1 (m1)
  !   2-3-1 (m2)
  !   3-1-2 (m1)
  !   3-1-2 (m2)
  ! m1: EI-LE-EI
  ! m2: LE-EI-LE
  !=======================================================
  Subroutine MOF_EILE3D(Phi, u, v, w, nl, dl, dt, rank, cx, cy, cz, Phi0)
    Implicit None
    Real(sp) :: dt
    Integer, Intent(In)     :: nl(3)
    Real(sp), Intent(In)    :: dl(3)
    Real(sp), Intent(InOut)        ::  Phi(0:nl(1) + 1, 0:nl(2) + 1, 0:nl(3) + 1)
    Real(sp), Intent(In), Optional :: Phi0(0:nl(1) + 1, 0:nl(2) + 1, 0:nl(3) + 1)
    Real(sp), Intent(InOut)        ::   cx(0:nl(1) + 1, 0:nl(2) + 1, 0:nl(3) + 1)
    Real(sp), Intent(InOut)        ::   cy(0:nl(1) + 1, 0:nl(2) + 1, 0:nl(3) + 1)
    Real(sp), Intent(InOut)        ::   cz(0:nl(1) + 1, 0:nl(2) + 1, 0:nl(3) + 1)
    Real(sp), Intent(In)           ::    u(0:nl(1) + 1, 0:nl(2) + 1, 0:nl(3) + 1)
    Real(sp), Intent(In)           ::    v(0:nl(1) + 1, 0:nl(2) + 1, 0:nl(3) + 1)
    Real(sp), Intent(In)           ::    w(0:nl(1) + 1, 0:nl(2) + 1, 0:nl(3) + 1)
    Real(sp), Dimension(0:nl(1) + 1, 0:nl(2) + 1, 0:nl(3) + 1) :: u1, u2, v1, v3, w2, w3
    Integer :: rank


    num_iter = 0
    grid_iter = 0
    t_rec = 0.0_sp
    t_adv = 0.0_sp
    t_cor = 0.0_sp
    Call Vel_Decomp(u,v,w,u1,v1,w2,u2,v3,w3,nl)
    Block
      Real(sp), Dimension(0:nl(1) + 1, 0:nl(2) + 1, 0:nl(3) + 1) :: div1, div2,div3
      Integer :: i, j, k
      print *, maxval(abs(u1))
      print *, maxval(abs(u2))
      print *, maxval(abs(v1))
      print *, maxval(abs(v3))
      print *, maxval(abs(w2)), maxloc(abs(w2))
      print *, maxval(abs(w3)), maxloc(abs(w3))
      print *, ''
      print *, maxval(abs(u))
      print *, maxval(abs(v))
      print *, maxval(abs(w))
      print *, ''
      Do k = 1, nl(3)
        Do j = 1, nl(2)
          Do i = 1, nl(1)
            div1(i,j,k) = u1(i,j,k) - u1(i-1,j,k) + v1(i,j,k) - v1(i,j-1,k)
            div2(i,j,k) = w2(i,j,k) - w2(i,j,k-1) + u2(i,j,k) - u2(i-1,j,k)
            div3(i,j,k) = v3(i,j,k) - v3(i,j-1,k) + w3(i,j,k) - w3(i,j,k-1)
          End Do
        End Do
      End Do
      print *, maxval(abs(div1))
      print *, maxval(abs(div2))
      print *, maxval(abs(div3))
      read(*,*)
    End Block
    If (rank == 0 .or. rank == 3) Then
      Call MOFAdv_EI(u1, cx, cy, cz, phi, nl, dl, dt, 1)
      Call MOFAdv_LE(v1, cx, cy, cz, phi, nl, dl, dt, 2)
      Call MOFAdv_EI(w2, cx, cy, cz, phi, nl, dl, dt, 3)
      Call MOFAdv_LE(u2, cx, cy, cz, phi, nl, dl, dt, 1)
      Call MOFAdv_EI(v3, cx, cy, cz, phi, nl, dl, dt, 2)
      Call MOFAdv_LE(w3, cx, cy, cz, phi, nl, dl, dt, 3)
    Else If (rank == 1 .or. rank == 4) Then
      Call MOFAdv_EI(w2, cx, cy, cz, phi, nl, dl, dt, 3)
      Call MOFAdv_LE(u2, cx, cy, cz, phi, nl, dl, dt, 1)
      Call MOFAdv_EI(v3, cx, cy, cz, phi, nl, dl, dt, 2)
      Call MOFAdv_LE(w3, cx, cy, cz, phi, nl, dl, dt, 3)
      Call MOFAdv_EI(u1, cx, cy, cz, phi, nl, dl, dt, 1)
      Call MOFAdv_LE(v1, cx, cy, cz, phi, nl, dl, dt, 2)
    Else If (rank == 2 .or. rank == 5) Then
      Call MOFAdv_EI(v3, cx, cy, cz, phi, nl, dl, dt, 2)
      Call MOFAdv_LE(w3, cx, cy, cz, phi, nl, dl, dt, 3)
      Call MOFAdv_EI(u1, cx, cy, cz, phi, nl, dl, dt, 1)
      Call MOFAdv_LE(v1, cx, cy, cz, phi, nl, dl, dt, 2)
      Call MOFAdv_EI(w2, cx, cy, cz, phi, nl, dl, dt, 3)
      Call MOFAdv_LE(u2, cx, cy, cz, phi, nl, dl, dt, 1)
    End If

  End Subroutine MOF_EILE3D

  !=======================================================
  ! [1-8-3] MOF-EIEALE 
  !     (Eulerian Implicit - 
  !      Eulerian Algebaric -
  !      Lanrangian Explicit)
  !-------------------------------------------------------
  ! Sweep sequence:
  !   1-2-3
  !   2-3-1
  !   3-1-2
  ! All: EI-EA-LE
  !=======================================================
  Subroutine MOF_EIEALE(Phi, u, v, w, nl, dl, dt, rank, cx, cy, cz, Phi0)
    Implicit None
    Real(sp) :: dt
    Integer, Intent(In)     :: nl(3)
    Real(sp), Intent(In)    :: dl(3)
    Real(sp), Intent(InOut)        ::  Phi(0:nl(1) + 1, 0:nl(2) + 1, 0:nl(3) + 1)
    Real(sp), Intent(In), Optional :: Phi0(0:nl(1) + 1, 0:nl(2) + 1, 0:nl(3) + 1)
    Real(sp), Intent(InOut)        ::   cx(0:nl(1) + 1, 0:nl(2) + 1, 0:nl(3) + 1)
    Real(sp), Intent(InOut)        ::   cy(0:nl(1) + 1, 0:nl(2) + 1, 0:nl(3) + 1)
    Real(sp), Intent(InOut)        ::   cz(0:nl(1) + 1, 0:nl(2) + 1, 0:nl(3) + 1)
    Real(sp), Intent(In)           ::    u(0:nl(1) + 1, 0:nl(2) + 1, 0:nl(3) + 1)
    Real(sp), Intent(In)           ::    v(0:nl(1) + 1, 0:nl(2) + 1, 0:nl(3) + 1)
    Real(sp), Intent(In)           ::    w(0:nl(1) + 1, 0:nl(2) + 1, 0:nl(3) + 1)
    Integer :: rank
    num_iter = 0
    grid_iter = 0
    t_rec = 0.0_sp
    t_adv = 0.0_sp
    t_cor = 0.0_sp
    If (rank == 0 .or. rank == 3) Then
      call MOFAdv_EI(u, cx, cy, cz, phi, nl, dl, dt, 1)
      call MOFAdv_EA(v, u, w, cx, cy, cz, phi, nl, dl, dt, 2)
      call MOFAdv_LE(w, cx, cy, cz, phi, nl, dl, dt, 3)
    Else If (rank == 1 .or. rank == 4) Then
      call MOFAdv_EI(w, cx, cy, cz, phi, nl, dl, dt, 3)
      call MOFAdv_EA(u, w, v, cx, cy, cz, phi, nl, dl, dt, 1)
      call MOFAdv_LE(v, cx, cy, cz, phi, nl, dl, dt, 2)
    Else If (rank == 2 .or. rank == 5) Then
      call MOFAdv_EI(v, cx, cy, cz, phi, nl, dl, dt, 2)
      call MOFAdv_EA(w, v, u, cx, cy, cz, phi, nl, dl, dt, 3)
      call MOFAdv_LE(u, cx, cy, cz, phi, nl, dl, dt, 1)
    EndIf
  End Subroutine MOF_EIEALE

  !=======================================================
  ! [1-8-4] MOF-EILE2D
  !     (Eulerian Implicit - 
  !      Lanrangian Explicit)
  !-------------------------------------------------------
  ! Sweep sequence:
  !   1-2
  !   2-1
  ! All: EI-LE
  !=======================================================
  Subroutine MOF_EILE2D(Phi, u, v, w, nl, dl, dt, rank, cx, cy, cz, Phi0)
    Implicit None
    Real(sp) :: dt
    Integer, Intent(In)     :: nl(3)
    Real(sp), Intent(In)    :: dl(3)
    Real(sp), Intent(InOut)        ::  Phi(0:nl(1) + 1, 0:nl(2) + 1, 0:nl(3) + 1)
    Real(sp), Intent(In), Optional :: Phi0(0:nl(1) + 1, 0:nl(2) + 1, 0:nl(3) + 1)
    Real(sp), Intent(InOut)        ::   cx(0:nl(1) + 1, 0:nl(2) + 1, 0:nl(3) + 1)
    Real(sp), Intent(InOut)        ::   cy(0:nl(1) + 1, 0:nl(2) + 1, 0:nl(3) + 1)
    Real(sp), Intent(InOut)        ::   cz(0:nl(1) + 1, 0:nl(2) + 1, 0:nl(3) + 1)
    Real(sp), Intent(In)           ::    u(0:nl(1) + 1, 0:nl(2) + 1, 0:nl(3) + 1)
    Real(sp), Intent(In)           ::    v(0:nl(1) + 1, 0:nl(2) + 1, 0:nl(3) + 1)
    Real(sp), Intent(In)           ::    w(0:nl(1) + 1, 0:nl(2) + 1, 0:nl(3) + 1)
    Integer :: rank
    num_iter = 0
    grid_iter = 0
    t_rec = 0.0_sp
    t_adv = 0.0_sp
    t_cor = 0.0_sp
    If (rank == 0 .or. rank == 2.or. rank == 4) Then
      call MOFAdv_EI(u, cx, cy, cz, phi, nl, dl, dt, 1)
      call MOFAdv_LE(w, cx, cy, cz, phi, nl, dl, dt, 2)
    Else If (rank == 1 .or. rank == 3.or. rank == 5) Then
      call MOFAdv_EI(v, cx, cy, cz, phi, nl, dl, dt, 2)
      call MOFAdv_LE(u, cx, cy, cz, phi, nl, dl, dt, 1)
    EndIf
  End Subroutine MOF_EILE2D

  !=======================================================
  ! (1-9) THINC-EI (Eulerian Implicit)
  !=======================================================
  Subroutine VOFTHINC_EI(Phi, u, v, w, nl, dl, dt, rank, Phi0)
    Implicit None
    Real(sp) :: dt
    Integer, Intent(In)            :: nl(3)
    Real(sp), Intent(In)           :: dl(3)
    Real(sp), Intent(InOut)        ::  Phi(0:nl(1) + 1, 0:nl(2) + 1, 0:nl(3) + 1)
    Real(sp), Intent(In), Optional :: Phi0(0:nl(1) + 1, 0:nl(2) + 1, 0:nl(3) + 1)
    Real(sp), Intent(In)           ::    u(0:nl(1) + 1, 0:nl(2) + 1, 0:nl(3) + 1)
    Real(sp), Intent(In)           ::    v(0:nl(1) + 1, 0:nl(2) + 1, 0:nl(3) + 1)
    Real(sp), Intent(In)           ::    w(0:nl(1) + 1, 0:nl(2) + 1, 0:nl(3) + 1)
    Integer :: rank
    num_iter = 0
    grid_iter = 0
    t_rec = 0.0_sp
    t_adv = 0.0_sp
    t_cor = 0.0_sp
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
  ! (1-10) THINC-LE (Lagrangian Explicit)
  !=======================================================
  Subroutine VOFTHINC_LE(Phi, u, v, w, nl, dl, dt, rank, Phi0)
    Implicit None
    Real(sp) :: dt
    Integer, Intent(In)            :: nl(3)
    Real(sp), Intent(In)           :: dl(3)
    Real(sp), Intent(InOut)        ::  Phi(0:nl(1) + 1, 0:nl(2) + 1, 0:nl(3) + 1)
    Real(sp), Intent(In), Optional :: Phi0(0:nl(1) + 1, 0:nl(2) + 1, 0:nl(3) + 1)
    Real(sp), Intent(In)           ::    u(0:nl(1) + 1, 0:nl(2) + 1, 0:nl(3) + 1)
    Real(sp), Intent(In)           ::    v(0:nl(1) + 1, 0:nl(2) + 1, 0:nl(3) + 1)
    Real(sp), Intent(In)           ::    w(0:nl(1) + 1, 0:nl(2) + 1, 0:nl(3) + 1)
    Integer :: rank
    num_iter = 0
    grid_iter = 0
    t_rec = 0.0_sp
    t_adv = 0.0_sp
    t_cor = 0.0_sp
    if (rank == 0 .or. rank == 3) Then
      call VOFAdvTHINC_LE(u, phi, nl, dl, dt, 1)
      call VOFAdvTHINC_LE(v, phi, nl, dl, dt, 2)
      call VOFAdvTHINC_LE(w, phi, nl, dl, dt, 3)
    else if (rank == 1 .or. rank == 4) Then
      call VOFAdvTHINC_LE(v, phi, nl, dl, dt, 2)
      call VOFAdvTHINC_LE(w, phi, nl, dl, dt, 3)
      call VOFAdvTHINC_LE(u, phi, nl, dl, dt, 1)
    else if (rank == 2 .or. rank == 5) Then
      call VOFAdvTHINC_LE(w, phi, nl, dl, dt, 3)
      call VOFAdvTHINC_LE(u, phi, nl, dl, dt, 1)
      call VOFAdvTHINC_LE(v, phi, nl, dl, dt, 2)
    endif
  End Subroutine VOFTHINC_LE

  !=======================================================
  ! (1-11) THINC-WY (Lagrangian Explicit)
  !=======================================================
  Subroutine VOFTHINC_WY(Phi, u, v, w, nl, dl, dt, rank, Phi0)
    Implicit None
    Real(sp) :: dt
    Integer, Intent(In)            :: nl(3)
    Real(sp), Intent(In)           :: dl(3)
    Real(sp), Intent(InOut)        ::  Phi(0:nl(1) + 1, 0:nl(2) + 1, 0:nl(3) + 1)
    Real(sp), Intent(In), Optional :: Phi0(0:nl(1) + 1, 0:nl(2) + 1, 0:nl(3) + 1)
    Real(sp), Intent(In)           ::    u(0:nl(1) + 1, 0:nl(2) + 1, 0:nl(3) + 1)
    Real(sp), Intent(In)           ::    v(0:nl(1) + 1, 0:nl(2) + 1, 0:nl(3) + 1)
    Real(sp), Intent(In)           ::    w(0:nl(1) + 1, 0:nl(2) + 1, 0:nl(3) + 1)
    Integer :: rank
    num_iter = 0
    grid_iter = 0
    t_rec = 0.0_sp
    t_adv = 0.0_sp
    t_cor = 0.0_sp
    if (rank == 0 .or. rank == 3) Then
      call VOFAdvTHINC_WY(u, phi0, phi, nl, dl, dt, 1)
      call VOFAdvTHINC_WY(v, phi0, phi, nl, dl, dt, 2)
      call VOFAdvTHINC_WY(w, phi0, phi, nl, dl, dt, 3)
    else if (rank == 1 .or. rank == 4) Then
      call VOFAdvTHINC_WY(v, phi0, phi, nl, dl, dt, 2)
      call VOFAdvTHINC_WY(w, phi0, phi, nl, dl, dt, 3)
      call VOFAdvTHINC_WY(u, phi0, phi, nl, dl, dt, 1)
    else if (rank == 2 .or. rank == 5) Then
      call VOFAdvTHINC_WY(w, phi0, phi, nl, dl, dt, 3)
      call VOFAdvTHINC_WY(u, phi0, phi, nl, dl, dt, 1)
      call VOFAdvTHINC_WY(v, phi0, phi, nl, dl, dt, 2)
    endif
  End Subroutine VOFTHINC_WY

  !=======================================================
  ! [1-12] THINC-EILE2D (Mixed Eulerian Implicit Lagrangian Explicit)
  !=======================================================
  Subroutine VOFTHINC_EILE2D(Phi, u, v, w, nl, dl, dt, rank, Phi0)
    Implicit None
    Real(sp) :: dt
    Integer, Intent(In)            :: nl(3)
    Real(sp), Intent(In)           :: dl(3)
    Real(sp), Intent(InOut)        ::  Phi(0:nl(1) + 1, 0:nl(2) + 1, 0:nl(3) + 1)
    Real(sp), Intent(In), Optional :: Phi0(0:nl(1) + 1, 0:nl(2) + 1, 0:nl(3) + 1)
    Real(sp), Intent(In)           ::    u(0:nl(1) + 1, 0:nl(2) + 1, 0:nl(3) + 1)
    Real(sp), Intent(In)           ::    v(0:nl(1) + 1, 0:nl(2) + 1, 0:nl(3) + 1)
    Real(sp), Intent(In)           ::    w(0:nl(1) + 1, 0:nl(2) + 1, 0:nl(3) + 1)
    Integer :: rank
    num_iter = 0
    grid_iter = 0
    t_rec = 0.0_sp
    t_adv = 0.0_sp
    t_cor = 0.0_sp
    if (rank == 0 .or. rank == 2.or. rank == 4) Then
      call VOFAdvTHINC_EI(u, phi, nl, dl, dt, 1)
      call VOFAdvTHINC_LE(v, phi, nl, dl, dt, 2)
    else if (rank == 1 .or. rank == 3.or. rank == 5) Then
      call VOFAdvTHINC_EI(v, phi, nl, dl, dt, 2)
      call VOFAdvTHINC_LE(u, phi, nl, dl, dt, 1)
    endif
  End Subroutine VOFTHINC_EILE2D

  !=======================================================
  ! [2-4] Advection algorithm of EA
  !-------------------------------------------------------
  ! Algebaric mapping
  !=======================================================
  Subroutine VOFAdv_EA(us, u1, u3, f, nl, dl, dt, dir)
    Implicit None
    Integer :: dir
    Real(sp) :: dt
    Integer, Intent(In)     :: nl(3)
    Real(sp), Intent(In)    :: dl(3)
    Real(sp), Intent(In)    :: us(0:nl(1) + 1, 0:nl(2) + 1, 0:nl(3) + 1)
    Real(sp), Intent(In)    :: u1(0:nl(1) + 1, 0:nl(2) + 1, 0:nl(3) + 1)
    Real(sp), Intent(In)    :: u3(0:nl(1) + 1, 0:nl(2) + 1, 0:nl(3) + 1)
    Real(sp), Intent(InOut) :: f(0:nl(1) + 1, 0:nl(2) + 1, 0:nl(3) + 1)
    Integer  :: i, j, k, ii, jj, kk
    Integer  :: ii1, jj1, kk1
    Integer  :: ii3, jj3, kk3
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

    vof1 = 0.0_sp
    vof3 = 0.0_sp
    ii1 = 0; jj1 = 0; kk1=0
    ii3 = 0; jj3 = 0; kk3=0
    ii  = 0; jj  = 0; kk =0
    if (dir .eq. 1) then
      ii = 1; kk1=1; jj3=1
    end if
    if (dir .eq. 2) then
      jj = 1; ii1=1; kk3 = 1
    end if
    if (dir .eq. 3) then
      kk = 1; jj1=1; ii3=1
    end if

    Call VOFNormalVectors(f, nl, flag, n1, n2, n3)

    ! Calculate flux
    Call VOF_EulerianFlux(us, f, nl, dl, dt, dir, &
                              n1, n2, n3, flag, &
                              vof1, vof2, vof3)

    Call CPU_TIME(tt1)  !!! cpu time for correction term

    !new values of f and  clip it: 0. <= f <= 1.
    do k = 1, nl(3)
      do j = 1, nl(2)
        do i = 1, nl(1)
          b1 = u1(i - ii1, j - jj1, k - kk1)*dt/dl(dir)
          b2 = u1(i, j, k)*dt/dl(dir)
          c1 = u3(i - ii3, j - jj3, k - kk3)*dt/dl(dir)
          c2 = u3(i, j, k)*dt/dl(dir)
          a1 = us(i - ii, j - jj, k - kk)*dt/dl(dir)
          a2 = us(i, j, k)*dt/dl(dir)
          exp_factor1 = 1.0_sp/(1.0_sp - b2 + b1)
          exp_factor2 = 1.0_sp + c2 - c1
          exp_factor = 1.0_sp / exp_factor1 / exp_factor2
          f(i, j, k) =  f(i, j, k) * exp_factor +&
                        (vof3(i - ii, j - jj, k - kk) + &
                        vof1(i + ii, j + jj, k + kk) &
                        -( vof1(i,j,k) + vof3(i,j,k) ) &
                        ) / exp_factor2
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
  End Subroutine VOFAdv_EA

  !=======================================================
  ! [2-8] Advection algorithm of EA MOF 
  !=======================================================
  Subroutine MOFAdv_EA(us, u1, u3, cx, cy, cz, f, nl, dl, dt, dir)
    Implicit None
    Integer :: dir
    Real(sp) :: dt
    Integer, Intent(In)     :: nl(3)
    Real(sp), Intent(In)    :: dl(3)
    Real(sp), Intent(In)    :: us(0:nl(1) + 1, 0:nl(2) + 1, 0:nl(3) + 1)
    Real(sp), Intent(In)    :: u1(0:nl(1) + 1, 0:nl(2) + 1, 0:nl(3) + 1)
    Real(sp), Intent(In)    :: u3(0:nl(1) + 1, 0:nl(2) + 1, 0:nl(3) + 1)
    Real(sp), Intent(InOut) :: f(0:nl(1) + 1, 0:nl(2) + 1, 0:nl(3) + 1)
    Real(sp), Intent(InOut) :: cx(0:nl(1) + 1, 0:nl(2) + 1, 0:nl(3) + 1)
    Real(sp), Intent(InOut) :: cy(0:nl(1) + 1, 0:nl(2) + 1, 0:nl(3) + 1)
    Real(sp), Intent(InOut) :: cz(0:nl(1) + 1, 0:nl(2) + 1, 0:nl(3) + 1)
    Integer  :: i, j, k
    Integer  :: ii, jj, kk
    Integer  :: ii1, jj1, kk1
    Integer  :: ii3, jj3, kk3
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
    Real(sp) :: exp_factor1, exp_factor2
    Real(sp) :: tt1, tt2

    vof1 = 0.0_sp
    vof3 = 0.0_sp
    ii1 = 0; jj1 = 0; kk1=0
    ii3 = 0; jj3 = 0; kk3=0
    ii  = 0; jj  = 0; kk =0
    if (dir .eq. 1) then
      ii = 1; kk1=1; jj3=1
    end if
    if (dir .eq. 2) then
      jj = 1; ii1=1; kk3 = 1
    end if
    if (dir .eq. 3) then
      kk = 1; jj1=1; ii3=1
    end if

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
          ! exp_factor = 1.0_sp/(1.0_sp - a2 + a1)
          b1 = u1(i - ii1, j - jj1, k - kk1)*dt/dl(dir)
          b2 = u1(i, j, k)*dt/dl(dir)
          c1 = u3(i - ii3, j - jj3, k - kk3)*dt/dl(dir)
          c2 = u3(i, j, k)*dt/dl(dir)
          exp_factor1 = 1.0_sp/(1.0_sp - b2 + b1)
          exp_factor2 = 1.0_sp + c2 - c1
          exp_factor = 1.0_sp / exp_factor1 / exp_factor2
          f(i, j, k) =  (vof3(i - ii, j - jj, k - kk) + &
                        f(i, j, k) / exp_factor1 + &
                        vof1(i + ii, j + jj, k + kk) &
                        -vof1(i,j,k) - vof3(i,j,k) &
                        ) / exp_factor2
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
  End Subroutine MOFAdv_EA

!=======================================================
! [2-9] Advection algorithm of THINC (EI)
!=======================================================
  Subroutine VOFAdvTHINC_EI(us, f, nl, dl, dt, dir)
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
    Call THINC_EulerianFlux(us, f, nl, dl, dt, dir, &
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

  End Subroutine VOFAdvTHINC_EI

!=======================================================
! [2-10] Advection algorithm of THINC (LE)
!=======================================================
  Subroutine VOFAdvTHINC_LE(us, f, nl, dl, dt, dir)
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
    Call THINC_LagrangianFlux(us, f, nl, dl, dt, dir, &
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
 
  End Subroutine VOFAdvTHINC_LE

!=======================================================
! [2-11] Advection algorithm of THINC (WY)
!=======================================================
  Subroutine VOFAdvTHINC_WY(us, f0, f, nl, dl, dt, dir)
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
    Call THINC_EulerianFlux(us, f, nl, dl, dt, dir, &
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

    ! apply proper boundary conditions to volume fraction
    call phi_bc%setBCS(f)

  End Subroutine VOFAdvTHINC_WY

  Subroutine THINC_EulerianFlux(us, f, nl, dl, dt, dir, &
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
    Real(sp) :: a1, a2, beta, gamma, betagamma2, x_center, alpha
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

          ! f = 1
          If (flag(i, j, k) .EQ. 0) Then
            vof1(i, j, k) = DMAX1(-a1, 0.0_sp)
            vof3(i, j, k) = DMAX1(a2, 0.0_sp)

            ! 0 < f < 1
          Else If (flag(i, j, k) .EQ. 1) Then
            !*(2) get gamma beta;
            norm(1) = n1(i,j,k); norm(2) = n2(i,j,k); norm(3) = n3(i,j,k)
            Call Normalization2(norm)
            if (norm(dir) .gt. 0.0_sp) Then
              gamma = -1.0_sp
            else
              gamma = 1.0_sp
            endif
            Call THINC_Beta(norm,f(i,j,k),dir,beta)
            ! Calculate beta
            betagamma2 = 2.0_sp*beta*gamma
            x_center = THINC1DBackward(f(i, j, k), betagamma2)
            !*(3) get fluxes
            x0 = 0.0_sp; deltax = 1.0_sp
            If (a1 .LT. 0.0_sp) Then
              deltax(dir) = -a1
              vof1(i, j, k) = THINC1DForward(x_center, betagamma2, x0(dir), deltax(dir))
            End If
            If (a2 .GT. 0.0_sp) Then
              x0(dir) = 1.0_sp - a2; 
              deltax(dir) = a2
              vof3(i, j, k) = THINC1DForward(x_center, betagamma2, x0(dir), deltax(dir))
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

  End Subroutine THINC_EulerianFlux

  Subroutine THINC_LagrangianFlux(us, f, nl, dl, dt, dir, &
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
    Real(sp) :: a1, a2, beta, gamma, betagamma2, x_center, alpha
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
            Call Normalization2(norm)
            if (norm(dir) .gt. 0.0_sp) Then
              gamma = -1.0_sp
            else
              gamma = 1.0_sp
            endif
            ! Calculate beta
            Call THINC_Beta(norm,f(i,j,k),dir,beta)
            betagamma2 = 2.0_sp*beta*gamma
            x_center = THINC1DBackward(f(i, j, k), betagamma2)
            betagamma2 = betagamma2 / (1.0_sp - a1 + a2)
            x_center = ( x_center ) * (1.0_sp - a1 + a2) + a1
            !*(3) get fluxes
            x0 = 0.0_sp; deltax = 1.0_sp
            if (a1 .LT. 0.0_sp) then
              x0(dir) = a1; 
              deltax(dir) = -a1
              vof1(i, j, k) = THINC1DForward(x_center, betagamma2, x0(dir), deltax(dir))
            end if
            if (a2 .GT. 0.0_sp) then
              x0(dir) = 1.0_sp; 
              deltax(dir) = a2
              vof3(i, j, k) = THINC1DForward(x_center, betagamma2, x0(dir), deltax(dir))
            end if
            x0(dir) = DMAX1(a1, 0.0_sp)
            deltax(dir) = 1.0_sp - x0(dir) + DMIN1(0.0_sp, a2)
            vof2(i, j, k) = THINC1DForward(x_center, betagamma2, x0(dir), deltax(dir))
          endif
        enddo
      enddo
    enddo

    Call CPU_TIME(tt2)  !!! cpu time for reconstruction
    t_adv = t_adv + tt2 - tt1

    ! apply proper boundary conditions to vof1, vof2, vof3
    call phi_bc%setBCS(vof1)
    call phi_bc%setBCS(vof3)

  End Subroutine THINC_LagrangianFlux

  !=======================================================
  ! [3-1] Velocity decomposition for EILE3D
  !=======================================================
  Subroutine Vel_Decomp(u,v,w,u1,v1,w1,u2,v2,w2,nl)
    Implicit None
    Integer,Intent(In) :: nl(3)
    Real(sp), Intent(In) , Dimension(0:nl(1) + 1, 0:nl(2) + 1, 0:nl(3) + 1) :: u, v, w
    Real(sp), Intent(Out), Dimension(0:nl(1) + 1, 0:nl(2) + 1, 0:nl(3) + 1) :: u1, u2, v1, v2, w1, w2

    Integer :: i,j,k
    u1 = u / 2.0_sp
    v1(:,0,:) = v(:,0,:) / 2.0_sp
    w1(:,:,0) = w(:,:,0) / 2.0_sp
    Do k = 1, nl(3)
      Do j = 1, nl(2)
        Do i = 1, nl(1)
          v1(i,j,k) = u1(i-1,j,k) - u1(i,j,k) + v1(i,j-1,  k)
          w1(i,j,k) = u1(i-1,j,k) - u1(i,j,k) + w1(i,  j,k-1)
        End Do
      End Do
    End Do
    u2 = u - u1
    v2 = v - v1
    w2 = w - w1
  End Subroutine Vel_Decomp

  !=======================================================
  ! [3-2] Centroid advection (Eulerian)
  !=======================================================
  Subroutine Centroid_Lagrangian_Adv(c, ul, ur, xl, xr, dir)
    Implicit None
    Real(sp), Intent(InOut) :: c(3)
    Real(sp), Intent(In)    :: ul, ur
    Real(sp), Intent(In)    :: xl, xr
    Integer,  Intent(In)    :: dir

    c(dir) = (1.0_sp + ur - ul) * c(dir) - (ur * xl - ul * xr)

  End Subroutine Centroid_Lagrangian_Adv

  !=======================================================
  ! [3-3] Centroid advection (Lagrangian)
  !=======================================================
  Subroutine Centroid_Eulerian_Adv(c, ul, ur, xl, xr, dir)
    Implicit None
    Real(sp), Intent(InOut) :: c(3)
    Real(sp), Intent(In)    :: ul, ur
    Real(sp), Intent(In)    :: xl, xr
    Integer,  Intent(In)    :: dir

    c(dir) =  (c(dir) - (ur * xl - ul * xr)) / (1.0_sp - ur + ul)

  End Subroutine Centroid_Eulerian_Adv

  Subroutine THINC_Beta_Constant(norm, f, dir, beta)
    Real(SP),Intent(In) :: norm(3)
    Real(SP),Intent(In) :: f
    Integer ,Intent(In) :: dir
    Real(SP),Intent(Out) :: beta
    beta = 3.5
  End Subroutine THINC_Beta_Constant

  Subroutine THINC_Beta_SW(norm, f, dir, beta)
    Real(SP),Intent(In) :: norm(3)
    Real(SP),Intent(In) :: f
    Integer ,Intent(In) :: dir
    Real(SP),Intent(Out) :: beta
    beta = 2.3_sp*dabs(norm(dir)) + 0.01_sp
  End Subroutine THINC_Beta_SW

  Subroutine THINC_Beta_SL(norm, f, dir, beta)
    Real(SP),Intent(In) :: norm(3)
    Real(SP),Intent(In) :: f
    Integer ,Intent(In) :: dir
    Real(SP),Intent(Out) :: beta
    beta = dabs(norm(dir) / dsqrt(1.0 - norm(dir)**2)+1.0d-15) * 2.0d0 * (dexp(dabs(f-0.5d0)))
    beta = min(max(1.0d-2,beta),1.0d2)
  End Subroutine THINC_Beta_SL

End Module ModVOFExt

