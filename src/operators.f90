Module Operators
  Use ModGlobal, only : sp
  Implicit None
Contains

  !==========================
  ! Calculate the central difference of a variable
  !    Face => Center
  !--------------------------
  ! Inputs:
  !    f: value of variable, 3D array
  !    dl: grid size, dx, dy, dz
  !    nl: dimensions
  !    dir: direction
  ! Outputs:
  !    dfdt: the gradient of f1
  !==========================
  Subroutine CentralDifference_F2C(f, dfdl, dl, nl, dir)
    Implicit None
    Real(sp), Intent(In)  :: f(0:,0:,0:)
    Real(sp), Intent(Out) :: dfdl(0:,0:,0:)
    Real(sp), Intent(In)  :: dl(3)
    Integer,  Intent(In)  :: nl(3)
    Integer :: dir
    Integer :: i, j, k
    Integer :: ii, jj, kk
    If (dir .eq. 1) Then
      ii = 1; jj = 0; kk = 0
    Else If (dir .eq. 2) Then
      ii = 0; jj = 1; kk = 0
    Else
      ii = 0; jj = 0; kk = 1
    End If

    Do k = 1, nl(3)
      Do j = 1, nl(2)
        Do i = 1, nl(1)
          dfdl(i,j,k) = ( f(i,j,k) - f(i-ii,j-jj,k-kk) ) / dl(dir)
        End Do
      End Do
    End Do
  End Subroutine CentralDifference_F2C

  !==========================
  ! Calculate the central difference of a variable
  !    Face => Center
  !--------------------------
  ! Inputs:
  !    f: value of variable, 3D array
  !    dl: grid size, dx, dy, dz
  !    nl: dimensions
  !    dir: direction
  ! Outputs:
  !    fc: averaged value of f at cell center
  !==========================
  Subroutine Average_F2C(f, fc, dl, nl, dir)
    Implicit None
    Real(sp), Intent(In)  :: f(0:,0:,0:)
    Real(sp), Intent(Out) :: fc(0:,0:,0:)
    Real(sp), Intent(In)  :: dl(3)
    Integer,  Intent(In)  :: nl(3)
    Integer :: dir
    Integer :: i, j, k
    Integer :: ii, jj, kk
    If (dir .eq. 1) Then
      ii = 1; jj = 0; kk = 0
    Else If (dir .eq. 2) Then
      ii = 0; jj = 1; kk = 0
    Else
      ii = 0; jj = 0; kk = 1
    End If

    Do k = 1, nl(3)
      Do j = 1, nl(2)
        Do i = 1, nl(1)
          fc(i,j,k) = ( f(i,j,k) + f(i-ii,j-jj,k-kk) ) / 0.5_sp
        End Do
      End Do
    End Do
  End Subroutine Average_F2C

  !==========================
  ! Calculate the central difference of a variable
  !    Face => Edge
  !--------------------------
  ! Inputs:
  !    f: value of variable, 3D array
  !    dl: grid size, dx, dy, dz
  !    nl: dimensions
  !    dir: direction
  ! Outputs:
  !    fe: averaged value of f at cell edge
  !==========================
  Subroutine Average_F2E(f, fc, dl, nl, dir)
    Implicit None
    Real(sp), Intent(In)  :: f(0:,0:,0:)
    Real(sp), Intent(Out) :: fv(0:,0:,0:)
    Real(sp), Intent(In)  :: dl(3)
    Integer,  Intent(In)  :: nl(3)
    Integer :: dir
    Integer :: i, j, k
    Integer :: ii, jj, kk
    If (dir .eq. 1) Then
      ii = 1; jj = 0; kk = 0
    Else If (dir .eq. 2) Then
      ii = 0; jj = 1; kk = 0
    Else
      ii = 0; jj = 0; kk = 1
    End If

    Do k = 1, nl(3)
      Do j = 1, nl(2)
        Do i = 1, nl(1)
          fc(i,j,k) = ( f(i,j,k) + f(i-ii,j-jj,k-kk) ) / 0.5_sp
        End Do
      End Do
    End Do
  End Subroutine Average_F2E


End Module Operators
