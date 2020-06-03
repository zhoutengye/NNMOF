!=======================================================
! Normalize the normal vector that satisfies
!    sigma norm(i) = 1
! and
!    all norm(i) > 0
!=======================================================
Function FloodSZ(norm)
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
End Function FloodSZ
