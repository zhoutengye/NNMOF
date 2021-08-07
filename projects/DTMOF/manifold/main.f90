#include "param.h"
# if defined TEST1
# define test_case gen_data
# elif defined TEST2
# define test_case compare
# elif defined TEST3
# define test_case compare_old_initial
# endif
program test
  use ModGlobal
  Implicit None

  call test_case

end program test

Subroutine gen_data
  Use ModVOFFunc
  Implicit None
  Integer, Parameter  :: num_sampling = 1000000
  Real(sp) :: n1(num_sampling)
  Real(sp) :: n2(num_sampling)
  Real(sp) :: n3(num_sampling)
  Real(sp) :: f(num_sampling)
  Real(sp) :: data(7,num_sampling)
  Real(sp) :: xc0(3), nr(3)
  Integer :: i
  Call random_number(n1)
  Call random_number(n2)
  Call random_number(n3)
  Call random_number(f)
  n1 = (n1 - 0.5_sp) * 2.0_sp
  n2 = (n2 - 0.5_sp) * 2.0_sp
  n3 = (n3 - 0.5_sp) * 2.0_sp
  open(10,file='fdata.dat',status='old')
  Do i = 1, num_sampling
    read(10,*) f(i)
    nr = (/n1(i),   n2(i),   n3(i)/)
    Call Normalization1(nr)
    Call FloodSZ_BackwardC(nr,f(i),xc0)
    data(1,i) = nr(1)
    data(2,i) = nr(2)
    data(3,i) = nr(3)
    data(4,i) = f(i)
    data(5,i) = xc0(1)
    data(6,i) = xc0(2)
    data(7,i) = xc0(3)
  End Do

  open(10,file='data.dat')
  Do i = 1, num_sampling
    Write(10,'(7F)')data(:,i)
  End Do
  close(10)


End subroutine gen_data

Subroutine compare
  Use ModGlobal
  Use ModTools
  Use ModVOF
  use mod_cg3_polyhedron
  use MODSussman
  use variables_mof
  Implicit None
  Integer,parameter :: num_sampling = 10000
  Real(sp) :: n1(num_sampling)
  Real(sp) :: n2(num_sampling)
  Real(sp) :: n3(num_sampling)
  Real(sp) :: p1(num_sampling)
  Real(sp) :: t1(num_sampling)
  Real(sp) :: ddx(3)
  Real(sp) :: num_iters(5,2)
  Real(sp) :: error(5,2)
  Real(sp) :: f, pt1(2), pt2(2), pt(2)
  Real(sp) :: xc0(3), nr(3)
  Integer :: sum_iter(2)
  Integer :: i
  Character(80) :: method
  Integer :: partition

  f = 0.3_SP

  Call Init(inputfield=.false.)

  Call random_number(n1)
  Call random_number(n2)
  Call random_number(n3)
  Call random_number(p1)
  Call random_number(t1)

  n1 = ( n1 - 0.5_SP ) * 2.0_SP
  n2 = ( n2 - 0.5_SP ) * 2.0_SP
  n3 = ( n3 - 0.5_SP ) * 2.0_SP
  p1 = p1 * 2.0 * Pi
  t1 = t1 * Pi

  Call Initialize_NN
  MOFNorm => MOFNNStab

  open(10,file='pt1.dat')
  open(11,file='pt2.dat')
  open(12,file='par.dat')
  Do i = 1, num_sampling
    nr = (/n1(i),   n2(i),   n3(i)/)
    pt = (/p1(i), t1(i)/)
    Call Angle2Norm(pt, nr)
    Call Norm2Angle(pt1, nr)
    Call Normalization1(nr)
    Call FloodSZ_BackwardC2(nr,f,xc0,partition)
    Call MOFNorm(f, xc0, nr)
    Call Norm2Angle(pt2, nr)
    Write(10,'(2F)') pt1
    Write(11,'(2F)') pt2
    Write(12,'(2F)') real(partition)
  End Do

end Subroutine compare

SUBROUTINE FloodSZ_BackwardC2(nr,cc,xc0,partition)

  IMPLICIT NONE
  REAL(8), INTENT(IN) :: nr(3),cc
  REAL(8), INTENT(OUT) :: xc0(3)
  Integer, INTENT(OUT) :: partition
  REAL(8) :: ctd0(3)
  REAL(8) :: cch,ccr,alh,c01,c02,c03,np1,np2,np3
  REAL(8) :: m1,m2,m3,m12,numer,denom,p,pst,q,arc,csarc
  REAL(8) :: tmp1,tmp2,tmp3,tmp4,a2,bot1
  INTEGER :: ind(3)
  REAL(8), PARAMETER :: athird=1.d0/3.d0,eps0=1.d-50 
  INTRINSIC DMAX1,DMIN1,DSQRT,DACOS,DCOS

  np1 = DABS(nr(1))                                 ! need positive coefficients
  np2 = DABS(nr(2))
  np3 = DABS(nr(3))
  ! coefficients in ascending order
  if (np1 <= np2) then                             
    m1 = np1
    m3 = np2
    ind(1) = 1                                            
    ind(3) = 2
  else
    m1 = np2
    m3 = np1
    ind(1) = 2
    ind(3) = 1
  endif

  if (np3 < m1) then
    m2 = m1
    m1 = np3
    ind(2) = ind(1)
    ind(1) = 3
  else if (np3 >= m3) then
    m2 = m3
    m3 = np3
    ind(2) = ind(3)
    ind(3) = 3
  else
    m2 = np3
    ind(2) = 3
  endif

  denom = DMAX1(6.d0*m1*m2*m3,eps0)                           
  cch = DMIN1(cc,1.d0-cc)                              ! limit to: 0 < cch < 1/2
  ccr = DMAX1(cch,eps0)                                     ! avoid cch = m1 = 0
  c01 = m1*m1*m1/denom                                          ! get cch ranges 
  c02  = c01 + 0.5d0*(m2-m1)/m3
  m12 = m1 + m2
  if (m12 <= m3) then
    c03 = 0.5d0*m12/m3
  else
    numer = m3*m3*(3.d0*m12-m3) + m1*m1*(m1-3.d0*m3) + m2*m2*(m2-3.d0*m3)
    c03 = numer/denom
  endif

  ! 1: C<=C1; 2: C1<=C<=C2; 3: C2<=C<=C3; 4: C3<=C<=1/2 (a: m12<=m3; b: m3<m12)) 
  if (ccr <= c01) then
    alh = 0.25d0*(denom*cch)**athird 
    ctd0(1) = alh/m1 
    ctd0(2) = alh/m2
    ctd0(3) = alh/m3                                                ! case (1)
    partition = 1
  else if (ccr <= c02) then 
    alh = 0.5d0*(m1 + DSQRT(m1*m1 + 8.d0*m2*m3*(cch - c01)))    
    tmp1 = (2.d0*alh-m1)
    tmp2 = (2.d0*alh*alh - 2.d0*alh*m1 + m1*m1)
    bot1 = DMAX1(4.d0*(3.d0*alh*alh - 3.d0*alh*m1 + m1*m1),eps0)
    tmp2 = tmp1*tmp2/bot1
    ctd0(1) = 0.5d0 - m1*tmp1/bot1 
    ctd0(2) = tmp2/m2
    ctd0(3) = tmp2/m3                                                ! case (2)
    partition = 2
  else if (ccr <= c03) then
    p = 2.d0*m1*m2                                                   
    q = 1.5d0*m1*m2*(m12 - 2.d0*m3*cch)
    pst = DSQRT(p)
    arc = athird*DACOS(q/(p*pst))
    csarc = DCOS(arc)
    alh = pst*(DSQRT(3.d0*(1.d0-csarc*csarc)) - csarc) + m12
    a2 = alh*alh
    tmp1 = (alh-m1)*(alh-m1)*(alh-m1)
    tmp2 = (alh-m2)*(alh-m2)*(alh-m2)
    bot1 = 4.d0*(a2*alh - tmp1 - tmp2)
    ctd0(1) = (a2*a2 - tmp1*(alh+3.d0*m1) - tmp2*(alh-m2))/(m1*bot1) 
    ctd0(2) = (a2*a2 - tmp1*(alh-m1) - tmp2*(alh+3.d0*m2))/(m2*bot1) 
    ctd0(3) = (a2*a2 - tmp1*(alh-m1) - tmp2*(alh-m2))/(m3*bot1)      ! case (3)
    partition = 3
  else if (m12 <= m3) then
    alh = m3*cch + 0.5d0*m12                                         
    bot1 = DMAX1((2.d0*alh - m12),eps0)
    ctd0(1) = 0.5d0 - m1/(6.d0*bot1)
    ctd0(2) = 0.5d0 - m2/(6.d0*bot1)                                ! case (4a) 
    ctd0(3) = ((3.d0*alh - 2.d0*m12)*bot1 + alh*m12 - m1*m2)/(6.d0*m3*bot1)
    partition = 4
  else
    p = m12*m3 + m1*m2 - 0.25d0                                     
    q = 1.5d0*m1*m2*m3*(0.5d0-cch)
    pst = DSQRT(p)
    arc = athird*DACOS(q/(p*pst))
    csarc = DCOS(arc)
    alh = pst*(DSQRT(3.d0*(1.d0-csarc*csarc)) - csarc) + 0.5d0 
    tmp1 = m1 + m2 + m3
    tmp2 = m1*m1 + m2*m2 + m3*m3
    tmp3 = m1*m1*m1 + m2*m2*m2 + m3*m3*m3
    tmp4 = m1*m1*m1*m1 + m2*m2*m2*m2 + m3*m3*m3*m3
    a2 = alh*alh
    bot1 = 4.d0*(2.d0*a2*alh - 3.d0*a2*tmp1 + 3.d0*alh*tmp2 - tmp3)
    tmp1 = 2.d0*a2*a2 - 4.d0*a2*alh*tmp1 + 6.d0*a2*tmp2 - 4.d0*alh*tmp3 + tmp4
    ctd0(1) = (tmp1 + 4.d0*m1*(alh-m1)*(alh-m1)*(alh-m1))/(m1*bot1)
    ctd0(2) = (tmp1 + 4.d0*m2*(alh-m2)*(alh-m2)*(alh-m2))/(m2*bot1)
    ctd0(3) = (tmp1 + 4.d0*m3*(alh-m3)*(alh-m3)*(alh-m3))/(m3*bot1) ! case (4b)
    partition = 5
  endif

  ! back to actual c value, but the centroid of a full cell is the cell center
  if (cc > 0.5d0) then
    ctd0(1) = (0.5d0 - (1.d0 - cc)*(1.d0-ctd0(1)))/cc
    ctd0(2) = (0.5d0 - (1.d0 - cc)*(1.d0-ctd0(2)))/cc
    ctd0(3) = (0.5d0 - (1.d0 - cc)*(1.d0-ctd0(3)))/cc
  endif

  ! get correct indexing
  xc0(ind(1)) = ctd0(1)                              
  xc0(ind(2)) = ctd0(2)
  xc0(ind(3)) = ctd0(3)

  ! take care of negative coefficients
  if (nr(1) < 0.d0) xc0(1) = 1.d0 - xc0(1) 
  if (nr(2) < 0.d0) xc0(2) = 1.d0 - xc0(2)
  if (nr(3) < 0.d0) xc0(3) = 1.d0 - xc0(3)

END SUBROUTINE FloodSZ_BackwardC2

