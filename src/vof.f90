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

  !=================
  ! Norm: Parker & Youngs
  ! Reconstruction: APPLIC
  ! Advection: APPLIC
  !=================
  Subroutine APPLIC(Phi, U, v, w, nl, dl, dt)
    Use ModGlobal, only : updthalo
    Implicit None
    Real(sp) :: dt
    Integer,Intent(In)     :: nl(3)
    Real(sp),Intent(In)    :: dl(3)
    Real(sp),Intent(InOut)    :: Phi(0:nl(1)+1,0:nl(2)+1,0:nl(3)+1)
    Real(sp),Intent(In)       :: u(0:nl(1)+1,0:nl(2)+1,0:nl(3)+1)
    Real(sp),Intent(In)       :: v(0:nl(1)+1,0:nl(2)+1,0:nl(3)+1)
    Real(sp),Intent(In)       :: w(0:nl(1)+1,0:nl(2)+1,0:nl(3)+1)
    Integer :: nexch(2)
    nexch = nl(1:2)
    call AdvSZ(u, phi, nl, dl, dt, 1)
    call AdvSZ(v, phi, nl, dl, dt, 2)
    call AdvSZ(w, phi, nl, dl, dt, 3)
  End Subroutine APPLIC

  !=================
  ! Norm: Parker & Youngs
  ! Reconstruction: PLIC
  ! Advection: SZ
  !=================
  Subroutine SZPLIC(Phi, u, v, w, nl, dl, dt)
    Implicit None
    Real(sp) :: dt
    Integer,Intent(In)     :: nl(3)
    Real(sp),Intent(In)    :: dl(3)
    Real(sp),Intent(InOut)    :: Phi(0:nl(1)+1,0:nl(2)+1,0:nl(3)+1)
    Real(sp),Intent(In)       :: u(0:nl(1)+1,0:nl(2)+1,0:nl(3)+1)
    Real(sp),Intent(In)       :: v(0:nl(1)+1,0:nl(2)+1,0:nl(3)+1)
    Real(sp),Intent(In)       :: w(0:nl(1)+1,0:nl(2)+1,0:nl(3)+1)
      call AdvSZ(u, phi, nl, dl, dt, 1)
      call AdvSZ(v, phi, nl, dl, dt, 2)
      call AdvSZ(w, phi, nl, dl, dt, 3)
  End Subroutine SZPLIC

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
