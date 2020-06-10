#include "param.h"
# if defined TEST1
# define test_case test1
# elif defined TEST2
# define test_case test2
# elif defined TEST3
# define test_case test3
# elif defined TEST4
# define test_case test4
# endif
program test
  use ModGlobal
  Implicit None

  call test_case

end program test

!===============
! test1: MOF reconstruction
! Given centroid and volume and obtains the normal vecor
!===============
Subroutine test1
  Use ModGlobal
  Use ModTools
  Use ModMOF
  Use Mod_VOF
  Implicit None
  Real(sp) :: f
  Real(sp) :: c(3)
  Real(sp) :: norm(3)
  Real(sp) :: abs_norm(3)

  Call Init(inputfield=.true.)

  ! Call MOF_Init

  ! f = 1.0_sp/3.0_sp
  ! c(1) = 1.0/4.0
  ! c(2) = 1.0/4.0
  ! c(3) = 1.0/4.0

  ! f = 2.0_sp/3.0_sp
  ! c(1) = 33.0/48.0
  ! c(2) = 33.0/48.0
  ! c(3) = 33.0/48.0

  f = 2.0_sp/3.0_sp
  c(1) = 33.0/48.0
  c(2) = 33.0/48.0
  c(3) = 33.0/48.0

  ! f = 0.8_sp
  ! c(1) = 0.5_sp
  ! c(2) = 0.4_sp
  ! c(3) = 0.5_sp

  ! f = 0.25_sp
  ! c(1) = 1.0_sp/3.0_sp
  ! c(2) = 1.0_sp/3.0_sp
  ! c(3) = 1.0_sp/2.0_sp

  ! f = 1.0_sp/24.0_sp
  ! c(1) = 1.0_sp/8.0_sp
  ! c(2) = 1.0_sp/8.0_sp
  ! c(3) = 1.0_sp/8.0_sp

  ! f = 0.5_sp
  ! c(1) = 0.2_sp
  ! c(2) = 0.49_sp
  ! c(3) = 0.49_sp

  ! f = 0.6_sp
  ! c(1) = 0.3
  ! c(2) = 0.5
  ! c(3) = 0.5

  ! Call NormMOF(f,c,norm,abs_norm)
  ! Call NormMOF(f,c,norm,abs_norm)

  c = c - 0.5_sp

  Call MOFZY(f,c,norm)

  print *, ''
  if (myid .eq. 0) Print *,norm

  Call MPI_FINALIZE(ierr)


end Subroutine test1

!===============
! test1: Centroid flooding
!===============
Subroutine test2
  Use ModGlobal
  Use ModTools
  Use ModMOF
  Use Mod_VOF
  Implicit None
  Real(sp) :: f
  Real(sp) :: alpha, x0(3), deltax(3)
  Real(sp) :: c(3)
  Real(sp) :: norm(3)

  Call Init(inputfield=.true.)

  ! Test whether centroid calculates correctly
  norm(1) = 1.0_sp / 3.0_sp
  norm(2) = 1.0_sp / 3.0_sp
  norm(3) = 1.0_sp / 3.0_sp
  f = 1.0_sp
  alpha = 1.0_sp / 3.0_sp
  x0 = 0.0_sp
  ! x0(1) = 0.5_sp
  x0(2) = 0.5_sp
  deltax = 1.0_sp
  Call FloodSZ_forwardC(norm, alpha, x0, deltax, f,c)
  if (myid .eq. 0) Print *, f, c

  ! Test whether centroid advection working
  ! Call Centroid_Lagrangian_Adv(c, 0.1_sp, 0.1_sp, 0.0_sp, 1.0_sp, 1)
  ! Call Centroid_Eulerian_Adv(c, 0.1_sp, 0.1_sp, 0.0_sp, 1.0_sp, 1)
  Call Centroid_Lagrangian_Adv(c, 0.1_sp, 0.2_sp, 0.0_sp, 1.0_sp, 1)
  if (myid .eq. 0) Print *, f, c
  Call Centroid_Eulerian_Adv(c, 0.1_sp, 0.2_sp, 0.0_sp, 1.0_sp, 2)
  if (myid .eq. 0) Print *, f, c

  Call MPI_FINALIZE(ierr)


end Subroutine test2


!===============
! test1: MOF reconstruction
! Given centroid and volume and obtains the normal vecor
!===============
Subroutine test3
  Use ModGlobal
  Use ModTools
  Use ModMOF
  Use Mod_VOF
  Implicit None
  Real(sp), allocatable, Dimension(:,:,:) :: f_beg
  Real(sp), allocatable, Dimension(:,:,:) :: f_end
  Real(sp) :: v1, v2
  Real(sp) :: v11, v12
  Integer  :: nexch(2)
  Integer  :: nn = 0
  Character(80) :: data_name
  Integer :: i, j, k

  Call Init(inputfield=.true.)

  Call MOF_Init
  u = 0.0
  v = 0.0
  w = 0.0
  u(1:nl(1),1:nl(2),1:nl(3)) = 1.0
  v(1:nl(1),1:nl(2),1:nl(3)) = 1.0
  w(1:nl(1),1:nl(2),1:nl(3)) = 1.0

  nexch = nl(1:2)
  Call u_bc%SetBCS(u)
  Call v_bc%SetBCS(v)
  Call w_bc%SetBCS(w)
  Call phi_bc%SetBCS(phi)
  Call phi_bc%SetBCS(cx)
  Call phi_bc%SetBCS(cy)
  Call phi_bc%SetBCS(cz)
  Allocate(f_beg(0:nl(1)+1,0:nl(2)+1,0:nl(3)+1))
  Allocate(f_end(0:nl(1)+1,0:nl(2)+1,0:nl(3)+1))
  f_beg = phi



  ! VOF advection
  ! Do While (time < tend)
    nn = nn + 1
    ! Call VOFCIAM(Phi, u, v, w, nl, dl, dt)
    ! Call VOFWY(Phi, u, v, w, nl, dl, dt)
    ! Call VOFHybrid(Phi, u, v, w, nl, dl, dt,nn)
    ! Call MOFCIAM(Phi, cx, cy, cz, u, v, w, nl, dl, dt)
    ! Call MOFWY(Phi, cx, cy, cz, u, v, w, nl, dl, dt)
    ! call AdvWY_MOF(u, cx, cy, cz, phi, nl, dl, dt, 1)
    ! Call Visual3DContour(f1=phi, slice_dir='x',slice_coord=8)
    ! call AdvWY_MOF(v, cx, cy, cz, phi, nl, dl, dt, 2)
    ! Call Visual3DContour(f1=phi, slice_dir='x',slice_coord=8)
    ! Call Visual3DContour(f1=phi)
    ! call AdvWY_MOF(w, cx, cy, cz, phi, nl, dl, dt, 3)
    ! Call Visual3DContour(f1=phi)
    ! call AdvCIAM_MOF(w, cx, cy, cz, phi, nl, dl, dt, 3)
    ! Call Visual3DContour(f1=phi)
    ! call AdvCIAM_MOF(v, cx, cy, cz, phi, nl, dl, dt, 2)
    ! Call Visual3DContour(f1=phi)
    ! call AdvCIAM_MOF(w, cx, cy, cz, phi, nl, dl, dt, 3)
    ! Call Visual3DContour(f1=phi)
    call AdvCIAM_MOF(u, cx, cy, cz, phi, nl, dl, dt, 1)
    ! Call Visual3DContour(f1=phi)
    ! call AdvCIAM_MOF(v, cx, cy, cz, phi, nl, dl, dt, 2)
    ! Call Visual3DContour(f1=phi)
    time =  time + dt
  ! End Do

    data_name = 'test'
  do nn = 1, n_vars
      Select Case(Trim(h5_output_field(nn)%groupname))
      Case('phi')
        Call HDF5WriteData(h5_output_field(nn), phi,data_name)
      Case('u')
        Call HDF5WriteData(h5_output_field(nn), u,data_name)
      Case('v')
        Call HDF5WriteData(h5_output_field(nn), v,data_name)
      Case('w')
        Call HDF5WriteData(h5_output_field(nn), w,data_name)
      Case('cx')
        Call HDF5WriteData(h5_output_field(nn), cx,data_name)
      Case('cy')
        Call HDF5WriteData(h5_output_field(nn), cy,data_name)
      Case('cz')
        Call HDF5WriteData(h5_output_field(nn), cz,data_name)
      End Select
  end do


  f_end = phi

  v1 = sum(f_beg(1:nl(1),1:nl(2),1:nl(3)))
  v2 = sum(f_end(1:nl(1),1:nl(2),1:nl(3)))

  Call MPI_Reduce(v1, v11, 1, MPI_REAL_SP, MPI_SUM, MPI_COMM_WORLD, 0, ierr)
  Call MPI_Reduce(v1, v12, 1, MPI_REAL_SP, MPI_SUM, MPI_COMM_WORLD, 0, ierr)
  if (myid .eq.0) then
    print *, v11
    print *, v12
  endif

  ! Call Visual3DContour(f1=f_end)
  Call Visual3DContour(f1=f_beg, f2=f_end)


  Call MPI_FINALIZE(ierr)


end Subroutine test3

