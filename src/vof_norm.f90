  !=======================================================
  ! Normalize the normal vector that satisfies
  !    sigma norm(i) = 1
  ! and
  !    all norm(i) > 0
  Subroutine Normalization1(norm)
    Use ModGlobal, only : sp
    Implicit None
    Real(sp), Intent(Inout) :: norm(3)
    Real(sp) :: aa
    norm = abs(norm)
    aa = norm(1) + norm(2) + norm(3)
    If( aa .gt. 1.d-10) Then
       norm(1) = norm(1) / aa
       norm(2) = norm(2) / aa
       norm(3) = norm(3) / aa
    End If
  End Subroutine Normalization1

  !=======================================================
  ! Normalize the normal vector that is a unit vector
  !     sigma norm(i)^2 = 1
  !=======================================================
  Subroutine Normalization2(norm)
    Use ModGlobal, only : sp
    Implicit None
    Real(sp), Intent(Inout) :: norm(3)
    Real(sp) :: aa
    aa = Sqrt( norm(1) * norm(1) + &
        &      norm(2) * norm(2) + &
        &      norm(3) * norm(3) )
    If( aa .gt. 1.d-10) Then
       norm(1) = norm(1) / aa
       norm(2) = norm(2) / aa
       norm(3) = norm(3) / aa
    End If
  End Subroutine Normalization2

  !=======================================================
  ! Normal vector APPLIC by akoi (2016)
  ! To be added
  !=======================================================
  Subroutine NormAPPLIC(f, norm)
    Use ModGlobal, only : sp
    Implicit None
    Real(sp), Intent(in)   :: F(3,3,3)
    Real(sp), Intent(out)  :: norm(3)
    norm(3) = 0
  End Subroutine NormAPPLIC

  !=======================================================
  ! Normal vector Parker and Youngs (1992)
  !=======================================================
  Subroutine NormParkaerYoungs(f, norm)
    Use ModGlobal, only : sp
    Implicit None
    Real(sp), Intent(in)   :: F(3,3,3)
    Real(sp), Intent(out)  :: norm(3)
    norm(1) =  ( F(3,1,2) + F(3,3,2) + 4.0_sp * F(3,2,2) + F(3,2,1) + F(3,2,3) ) &
        &    - ( F(1,1,2) + F(1,3,2) + 4.0_sp * F(1,2,2) + F(1,2,1) + F(1,2,3) )
    norm(2) =  ( F(1,3,2) + F(3,3,2) + 4.0_sp * F(2,3,2) + F(2,3,1) + F(2,3,3) ) &
        &    - ( F(1,1,2) + F(3,1,2) + 4.0_sp * F(2,1,2) + F(2,1,1) + F(2,1,3) )
    norm(3) =  ( F(1,2,3) + F(3,2,3) + 4.0_sp * F(2,2,3) + F(2,1,3) + F(2,3,3) ) &
        &    - ( F(1,2,1) + F(3,2,1) + 4.0_sp * F(2,2,1) + F(2,1,1) + F(2,3,1) )
  End Subroutine NormParkaerYoungs

!===============================================================
! Normal vector LIVERA
!===============================================================
Subroutine NormLVERA
  Use ModGlobal, only : sp
  Implicit None
End Subroutine NormLVERA

!===============================================================
! Normal vector MOF
!===============================================================
Subroutine NormMOF
  Use ModGlobal, only : sp
  Implicit None
End Subroutine NormMOF
