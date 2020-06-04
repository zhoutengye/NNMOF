  !=======================================================
  ! Normalize the normal vector that satisfies
  !    sigma norm(i) = 1
  ! and
  !    all norm(i) > 0
  Subroutine Normalization1(norm, abs_norm)
    Use ModGlobal, only : sp
    Implicit None
    Real(sp), Intent(Inout) :: norm(3)
    Real(sp), Intent(Out) :: abs_norm(3)
    Real(sp) :: aa

    abs_norm = abs(norm)
    aa = abs_norm(1) + abs_norm(2) + abs_norm(3)
    If( aa .gt. 1.d-10) Then
      norm = norm / aa
      abs_norm = abs_norm / aa
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
  End Subroutine NormAPPLIC

  !=======================================================
  ! Normal vector Parker and Youngs (1992)
  !=======================================================
  Subroutine NormParkerYoungs(f, norm)
    Use ModGlobal, only : sp
    Implicit None
    Real(sp), Intent(in)   :: F(3,3,3)
    Real(sp), Intent(out)  :: norm(3)
    norm(3) = 0.0_sp
    norm(1) =  ( F(1,1,2) + F(1,3,2) + 4.0_sp * F(1,2,2) + F(1,2,1) + F(1,2,3) ) &
        - ( F(3,1,2) + F(3,3,2) + 4.0_sp * F(3,2,2) + F(3,2,1) + F(3,2,3) ) 
        norm(2) =  ( F(1,1,2) + F(3,1,2) + 4.0_sp * F(2,1,2) + F(2,1,1) + F(2,1,3) ) &
        - ( F(1,3,2) + F(3,3,2) + 4.0_sp * F(2,3,2) + F(2,3,1) + F(2,3,3) ) 
        norm(3) =  ( F(1,2,1) + F(3,2,1) + 4.0_sp * F(2,2,1) + F(2,1,1) + F(2,3,1) ) &
        - ( F(1,2,3) + F(3,2,3) + 4.0_sp * F(2,2,3) + F(2,1,3) + F(2,3,3) )
  End Subroutine NormParkerYoungs

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
