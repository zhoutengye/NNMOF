!=======================================================
! Normalize the normal vector that satisfies
!    sigma norm(i) = 1
! and
!    all norm(i) > 0
!=======================================================
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
! Normal vector Youngs (1992)
! To be added
!=======================================================
Subroutine NormParkerYoungs(f, norm, abs_norm)
  Use ModGlobal, only : sp
  Implicit None
  Real(sp), Intent(in)   :: F(3,3,3)
  Real(sp), Intent(out)  :: norm(3)
  Real(sp), Intent(out)  :: abs_norm(3)
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

  !*(2) mx,my,mz>0. and mx+my+mz = 1.;
  Call Normalization1(norm, abs_norm)
End Subroutine NormParkerYoungs

!=======================================================
! Normal vector Central Scheme
! * returns normal normalized so that |mx|+|my|+|mz| = 1*
! Adopted from Paris Simulator by Zakeski
!=======================================================
Subroutine NormCS(c,norm,abs_norm)
  !***
  Implicit None
  Real(8) c(3,3,3)
  Real(8) norm(3), abs_norm(3)
  Real(8) m1,m2,m(0:3,0:2),t0,t1,t2
  Integer cn

  ! write the plane as: sgn(mx) X =  my Y +  mz Z + alpha 
  !                           m00 X = m01 Y + m02 Z + alpha 

  m1 = c(1,2,1) + c(1,2,3) + c(1,1,2) + c(1,3,2) + &
      c(1,2,2)
  m2 = c(3,2,1) + c(3,2,3) + c(3,1,2) + c(3,3,2) + &
      c(3,2,2)

  if(m1>m2) then
    m(0,0) = 1.
  else
    m(0,0) = -1.
  end if

  m1 = c(1,1,2)+ c(3,1,2)+ c(2,1,2)
  m2 = c(1,3,2)+ c(3,3,2)+ c(2,3,2)
  m(0,1) = 0.5*(m1-m2)

  m1 = c(1,2,1)+ c(3,2,1)+ c(2,2,1)
  m2 = c(1,2,3)+ c(3,2,3)+ c(2,2,3)
  m(0,2) = 0.5*(m1-m2)

  ! write the plane as: sgn(my) Y =  mx X +  mz Z + alpha, 
  !                          m11 Y = m10 X + m12 Z + alpha.
  m1 = c(1,1,2) + c(1,3,2) + c(1,2,2)
  m2 = c(3,1,2) + c(3,3,2) + c(3,2,2)
  m(1,0) = 0.5*(m1-m2)

  m1 = c(2,1,1) + c(2,1,3) + c(3,1,2) + c(1,1,2) +&
      c(2,1,2)
  m2 = c(2,3,1) + c(2,3,3) + c(3,3,2) + c(1,3,2) +&
      c(2,3,2)

  if(m1>m2) then
    m(1,1) = 1.
  else
    m(1,1) = -1.
  end if

  m1 = c(2,1,1)+ c(2,2,1)+ c(2,3,1)
  m2 = c(2,1,3)+ c(2,2,3)+ c(2,3,3)
  m(1,2) = 0.5*(m1-m2)

  m1 = c(1,2,1)+ c(1,2,3)+ c(1,2,2)
  m2 = c(3,2,1)+ c(3,2,3)+ c(3,2,2)
  m(2,0) = 0.5*(m1-m2)

  m1 = c(2,1,1)+ c(2,1,3)+ c(2,1,2)
  m2 = c(2,3,1)+ c(2,3,3)+ c(2,3,2)
  m(2,1) = 0.5*(m1-m2)

  m1 = c(1,2,1) + c(3,2,1) + c(2,1,1) + c(2,3,1) +&
      c(2,2,1)
  m2 = c(1,2,3) + c(3,2,3) + c(2,1,3) + c(2,3,3) +&
      c(2,2,3)

  if(m1>m2) then
    m(2,2) = 1.
  else
    m(2,2) = -1.
  end if

  ! normalize each set (mx,my,mz): |mx|+|my|+|mz| = 1

  t0 = DABS(m(0,0)) + DABS(m(0,1)) + DABS(m(0,2))
  m(0,0) = m(0,0)/t0
  m(0,1) = m(0,1)/t0
  m(0,2) = m(0,2)/t0

  t0 = DABS(m(1,0)) + DABS(m(1,1)) + DABS(m(1,2))
  m(1,0) = m(1,0)/t0
  m(1,1) = m(1,1)/t0
  m(1,2) = m(1,2)/t0

  t0 = DABS(m(2,0)) + DABS(m(2,1)) + DABS(m(2,2))
  m(2,0) = m(2,0)/t0
  m(2,1) = m(2,1)/t0
  m(2,2) = m(2,2)/t0

  ! choose among the three central schemes */ 
  t0 = DABS(m(0,0))
  t1 = DABS(m(1,1))
  t2 = DABS(m(2,2))

  cn = 0
  if (t1 > t0) then
    t0 = t1
    cn = 1
  endif

  if (t2 > t0) cn = 2

  ! components of the normal vector */
  norm(1) = m(cn,0)
  norm(2) = m(cn,1)
  norm(3) = m(cn,2)
  abs_norm(1) = abs(norm(1))
  abs_norm(2) = abs(norm(2))
  abs_norm(3) = abs(norm(3))

  return
End Subroutine NormCS

!=======================================================
! Mixed Youngs and Central difference
! * returns normal normalized so that |mx|+|my|+|mz| = 1*
! Adopted from Paris Simulator by Zakeski
!=======================================================
Subroutine NormMYCS(c,norm,abs_norm)
  !***
  Implicit None
  Real(8) c(3,3,3)
  Real(8) norm(3), abs_norm(3)
  Real(8) m1,m2,m(0:3,0:2),t0,t1,t2
  Real(8) mm(3), abs_mm(3)
  Integer cn
  ! write the plane as: sgn(mx) X =  my Y +  mz Z + alpha 
  !                           m00 X = m01 Y + m02 Z + alpha 

  m1 = c(1,2,1) + c(1,2,3) + c(1,1,2) + c(1,3,2) + &
      c(1,2,2)
  m2 = c(3,2,1) + c(3,2,3) + c(3,1,2) + c(3,3,2) + &
      c(3,2,2)

  if(m1>m2) then
    m(0,0) = 1.
  else
    m(0,0) = -1.
  end if

  m1 = c(1,1,2)+ c(3,1,2)+ c(2,1,2)
  m2 = c(1,3,2)+ c(3,3,2)+ c(2,3,2)
  m(0,1) = 0.5*(m1-m2)

  m1 = c(1,2,1)+ c(3,2,1)+ c(2,2,1)
  m2 = c(1,2,3)+ c(3,2,3)+ c(2,2,3)
  m(0,2) = 0.5*(m1-m2)

  ! write the plane as: sgn(my) Y =  mx X +  mz Z + alpha, 
  !                          m11 Y = m10 X + m12 Z + alpha.
  m1 = c(1,1,2) + c(1,3,2) + c(1,2,2)
  m2 = c(3,1,2) + c(3,3,2) + c(3,2,2)
  m(1,0) = 0.5*(m1-m2)

  m1 = c(2,1,1) + c(2,1,3) + c(3,1,2) + c(1,1,2) +&
      c(2,1,2)
  m2 = c(2,3,1) + c(2,3,3) + c(3,3,2) + c(1,3,2) +&
      c(2,3,2)

  if(m1>m2) then
    m(1,1) = 1.
  else
    m(1,1) = -1.
  end if

  m1 = c(2,1,1)+ c(2,2,1)+ c(2,3,1)
  m2 = c(2,1,3)+ c(2,2,3)+ c(2,3,3)
  m(1,2) = 0.5*(m1-m2)

  m1 = c(1,2,1)+ c(1,2,3)+ c(1,2,2)
  m2 = c(3,2,1)+ c(3,2,3)+ c(3,2,2)
  m(2,0) = 0.5*(m1-m2)

  m1 = c(2,1,1)+ c(2,1,3)+ c(2,1,2)
  m2 = c(2,3,1)+ c(2,3,3)+ c(2,3,2)
  m(2,1) = 0.5*(m1-m2)

  m1 = c(1,2,1) + c(3,2,1) + c(2,1,1) + c(2,3,1) +&
      c(2,2,1)
  m2 = c(1,2,3) + c(3,2,3) + c(2,1,3) + c(2,3,3) +&
      c(2,2,3)

  if(m1>m2) then
    m(2,2) = 1.
  else
    m(2,2) = -1.
  end if

  ! normalize each set (mx,my,mz): |mx|+|my|+|mz| = 1

  t0 = DABS(m(0,0)) + DABS(m(0,1)) + DABS(m(0,2))
  m(0,0) = m(0,0)/t0
  m(0,1) = m(0,1)/t0
  m(0,2) = m(0,2)/t0

  t0 = DABS(m(1,0)) + DABS(m(1,1)) + DABS(m(1,2))
  m(1,0) = m(1,0)/t0
  m(1,1) = m(1,1)/t0
  m(1,2) = m(1,2)/t0

  t0 = DABS(m(2,0)) + DABS(m(2,1)) + DABS(m(2,2))
  m(2,0) = m(2,0)/t0
  m(2,1) = m(2,1)/t0
  m(2,2) = m(2,2)/t0

  ! choose among the three central schemes */ 
  t0 = DABS(m(0,0))
  t1 = DABS(m(1,1))
  t2 = DABS(m(2,2))

  cn = 0
  if (t1 > t0) then
    t0 = t1
    cn = 1
  endif

  if (t2 > t0) cn = 2

  ! Youngs-CIAM scheme */
  Call NormParkerYoungs(c,mm,abs_mm)
  m(3,0:2) = mm(1:3)

  ! ! choose between the previous choice and Youngs-CIAM 
  if (DABS(m(cn,cn)) > maxval(abs_mm))  cn = 3

  ! components of the normal vector */
  norm(1:3) = m(cn,0:2)
  abs_norm = abs(norm)

  return
End Subroutine NormMYCS

Subroutine NormELVIRA(c, norm, abs_norm)
  Implicit None
  Real(8) c(3,3,3)
  Real(8) norm(3), abs_norm(3)
  Real(8) m1,m2,m(0:3,0:2),t0,t1,t2
  Real(8) mm(3), abs_mm(3)
  Integer cn

End Subroutine NormELVIRA

!===============================================================
! Normal vector MOF
! f: vof function
! c: centroid
!===============================================================
Subroutine NormMOF(f,c,norm,abs_norm)
  Use ModMOF
  Implicit None
  Real(8), Intent(In) :: f
  Real(8), Intent(In) :: c(3)
  Real(8), Intent(Out) :: norm(3), abs_norm(3)
  Call NormSussmanMOF(f,c,norm)
  Call Normalization1(norm, abs_norm)
End Subroutine NormMOF


