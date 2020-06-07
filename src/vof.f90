Module ModVOF
  Use ModGlobal, only : sp
  Implicit None
Contains
  !=================
  ! Norm: Parker & Youngs
  ! Reconstruction: THINC
  ! Advection: THINC
  !=================
  Subroutine THINCSW
    Implicit None
  End Subroutine THINCSW

  ! Subroutine VOFAdv(Phi, U, v, w, nl, dl, dt)
  !   Integer,Intent(In)           :: nl(3)
  !   Real(sp),Intent(In)          :: dl(3)
  !   Real(sp),Intent(InOut)       :: Phi(0:nl(1)+1,0:nl(2)+1,0:nl(3)+1)
  !   Real(sp),Intent(In)          :: u(0:nl(1)+1,0:nl(2)+1,0:nl(3)+1)
  !   Real(sp),Intent(In)          :: v(0:nl(1)+1,0:nl(2)+1,0:nl(3)+1)
  !   Real(sp),Intent(In)          :: w(0:nl(1)+1,0:nl(2)+1,0:nl(3)+1)
  !   Real(sp),Intent(In),optional :: u(0:nl(1)+1,0:nl(2)+1,0:nl(3)+1)
  !   Real(sp),Intent(In),optional :: v(0:nl(1)+1,0:nl(2)+1,0:nl(3)+1)
  !   Real(sp),Intent(In),optional :: w(0:nl(1)+1,0:nl(2)+1,0:nl(3)+1)
  !   Real(sp) :: dt
  ! End Subroutine VOFAdv

  !=================
  ! Norm: Parker & Youngs
  ! Reconstruction: APPLIC
  ! Advection: APPLIC
  !=================
  Subroutine VOFCIAM(Phi, u, v, w, nl, dl, dt)
    Use ModGlobal, only : updthalo
    Implicit None
    Real(sp) :: dt
    Integer,Intent(In)     :: nl(3)
    Real(sp),Intent(In)    :: dl(3)
    Real(sp),Intent(InOut)    :: Phi(0:nl(1)+1,0:nl(2)+1,0:nl(3)+1)
    Real(sp),Intent(In)       :: u(0:nl(1)+1,0:nl(2)+1,0:nl(3)+1)
    Real(sp),Intent(In)       :: v(0:nl(1)+1,0:nl(2)+1,0:nl(3)+1)
    Real(sp),Intent(In)       :: w(0:nl(1)+1,0:nl(2)+1,0:nl(3)+1)
    call AdvCIAM(u, phi, nl, dl, dt, 1)
    call AdvCIAM(v, phi, nl, dl, dt, 2)
    call AdvCIAM(w, phi, nl, dl, dt, 3)
  End Subroutine VOFCIAM

  Subroutine VOFWY(Phi, u, v, w, nl, dl, dt)
    Use ModGlobal, only : updthalo
    Implicit None
    Real(sp) :: dt
    Integer,Intent(In)     :: nl(3)
    Real(sp),Intent(In)    :: dl(3)
    Real(sp),Intent(InOut)    :: Phi(0:nl(1)+1,0:nl(2)+1,0:nl(3)+1)
    Real(sp),Intent(In)       :: u(0:nl(1)+1,0:nl(2)+1,0:nl(3)+1)
    Real(sp),Intent(In)       :: v(0:nl(1)+1,0:nl(2)+1,0:nl(3)+1)
    Real(sp),Intent(In)       :: w(0:nl(1)+1,0:nl(2)+1,0:nl(3)+1)
    call AdvWY(u, phi, nl, dl, dt, 1)
    call AdvWY(v, phi, nl, dl, dt, 2)
    call AdvWY(w, phi, nl, dl, dt, 3)
  End Subroutine VOFWY

  Subroutine MOFWY(Phi, u, v, w, nl, dl, dt)
    Use ModGlobal, only : updthalo
    Implicit None
    Real(sp) :: dt
    Integer,Intent(In)     :: nl(3)
    Real(sp),Intent(In)    :: dl(3)
    Real(sp),Intent(InOut)    :: Phi(0:nl(1)+1,0:nl(2)+1,0:nl(3)+1)
    Real(sp),Intent(In)       :: u(0:nl(1)+1,0:nl(2)+1,0:nl(3)+1)
    Real(sp),Intent(In)       :: v(0:nl(1)+1,0:nl(2)+1,0:nl(3)+1)
    Real(sp),Intent(In)       :: w(0:nl(1)+1,0:nl(2)+1,0:nl(3)+1)
    call AdvWY(u, phi, nl, dl, dt, 1)
    call AdvWY(v, phi, nl, dl, dt, 2)
    call AdvWY(w, phi, nl, dl, dt, 3)
  End Subroutine MOFWY


  !=================
  ! Norm: MOF
  ! Reconstruction: PLIC
  ! Advection: SZ with centroid
  !=================
  Subroutine MOFSussmanSZ
    Implicit None
  End Subroutine MOFSussmanSZ

  !=================
  ! Norm: DNN
  ! Reconstruction: APPLIc
  ! Advection: APPLIC with centroid
  !=================
  Subroutine NNMOF
    Implicit None
  End Subroutine NNMOF

End Module MODVOF
