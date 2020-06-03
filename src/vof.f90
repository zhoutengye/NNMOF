Module VOF
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
  Subroutine APPLIC
    Implicit None
  End Subroutine APPLIC

  !=================
  ! Norm: Parker & Youngs
  ! Reconstruction: PLIC
  ! Advection: SZ
  !=================
  Subroutine SZPLIC(Phi, u, v, w)
    Implicit None
      call AdvSZ(u, phi, nl, dl, dt, 1)
      call AdvSZ(v, phi, nl, dl, dt, 2)
      ! call AdvSZ(w, phi, nl, dl, dt, 3)
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

End Module VOF
