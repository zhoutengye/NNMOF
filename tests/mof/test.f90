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
! test1: MISC:
!    (1) Norm2Angle and Angle2Norm
!    (2) Flood_BackwardC
!    (3) FindCentroid
!    (4) MOFZY
!===============
Subroutine test1
  Use ModGlobal
  Use ModTools
  Use ModVOF
  Implicit None
  Real(sp) :: f
  Real(sp) :: c(3)
  Real(sp) :: norm(3)
  Real(sp) :: init_norm(3)
  Real(sp) :: angle(2)

  Call Init(inputfield=.false.)

  ! ! ==================Test norm2angle and angle2norm
  ! f = 1.0/6.0_sp
  ! norm = (/ 1.0/3.0_sp, 1.0/3.0_sp, 1.0/3.0_sp /)
  ! norm = (/ -1.0/3.0_sp, 1.0/3.0_sp, 1.0/3.0_sp /)
  ! norm = (/ 1.0/3.0_sp, -1.0/3.0_sp, 1.0/3.0_sp /)
  ! norm = (/ 1.0/3.0_sp, 1.0/3.0_sp, -1.0/3.0_sp /)
  ! norm = (/ 1.0_sp, 0.0/3.0_sp, 0.0/3.0_sp /)
  ! norm = (/ 0.0_sp, 1.0/1.0_sp, 0.0/3.0_sp /)
  ! norm = (/ 0.0_sp, 0.0/3.0_sp, 1.0/1.0_sp /)
  ! ! norm = (/ 0.0_sp, 1.0/3.0_sp, 2.0/3.0_sp /)
  ! ! norm = (/ 0.5_sp, 0.0/1.0_sp, 1.0/2.0_sp /)
  ! norm = (/ -0.1_sp, 0.9/1.0_sp, 1.0/1.0_sp /)
  ! Call Normalization2(norm)
  ! if (myid .eq. 0) Print *, norm
  ! Call norm2angle(angle, norm)
  ! if (myid .eq. 0) Print *, angle
  ! Call angle2norm(angle, norm)
  ! if (myid .eq. 0) Print *, norm

  !  ! ==================Test backward flood algorithm for centroid
  ! f = 1.0/6.0_sp
  ! norm = (/ 1.0/3.0_sp, 1.0/3.0_sp, 1.0/3.0_sp /)
  ! norm = (/ -1.0/3.0_sp, -1.0/3.0_sp, -1.0/3.0_sp /)
  ! norm = (/ 1.0/3.0_sp, -1.0/3.0_sp, -1.0/3.0_sp /)
  ! norm = (/ -1.0/3.0_sp, 1.0/3.0_sp, -1.0/3.0_sp /)
  ! norm = (/ 0.0/1.0_sp, 1.0/1.0_sp, -0.0/3.0_sp /)
  ! norm = (/ 1.0/1.0_sp, 0.0/3.0_sp, -0.0/3.0_sp /)
  ! f = 1.0/2.0_sp
  ! norm = (/ 1.0/2.0_sp, 1.0/2.0_sp, -0.0/3.0_sp /)
  ! norm = (/ -1.0/2.0_sp, -1.0/2.0_sp, -0.0/3.0_sp /)

  ! f = 2.0/3.0_sp
  ! norm = (/ 1.0/3.0_sp, 1.0/3.0_sp, 1.0/3.0_sp /)
  ! ! norm = (/ 1.0/1.0_sp, 0.0/3.0_sp, -0.0/3.0_sp /)
  ! ! norm = (/ -1.0/1.0_sp, 0.0/3.0_sp, -0.0/3.0_sp /)

  ! Call FloodSZ_BackwardC(norm, f, c)
  ! print *, c

  ! ! ==========Test FindCentroid subroutine
  ! f = 1.0/6.0_sp
  ! ! norm = (/ 1.0/3.0_sp, 1.0/3.0_sp, 1.0/3.0_sp /)
  ! norm = (/ 1.0/1.0_sp, 0.0/3.0_sp, -0.0/3.0_sp /)
  ! Call Normalization2(norm)
  ! Call Norm2Angle(angle, norm)
  ! Call Angle2Norm(angle, norm)
  ! Call FindCentroid(angle, f, c)
  ! print *, c

  ! ===============Test MOF reconstruction
  f = 1.0/48.0_sp
  c = (/ 1.0/8.0_sp, 1.0/8.0_sp, 1.0/8.0_sp /)
  init_norm = (/ 1.0/3.0_sp, 1.0/3.0_sp, 1.0/4.0_sp /)

  f = 1.0/6.0_sp
  c = (/ 1.0/4.0_sp, 1.0/4.0_sp, 1.0/4.0_sp /)
  ! init_norm = (/ 1.0/1.0_sp, 1.0/8.0_sp, 1.0/8.0_sp /)
  init_norm = (/ 1.0/3.0_sp, 1.0/3.0_sp, 1.0/4.0_sp /)

  !  --------------Different initial guess
  ! init_norm = (/ -1.0/3.0_sp, -1.0/3.0_sp, -1.0/3.0_sp /)
  Call Normalization2(init_norm)
  ! Call MOFLemoine(f,c,norm, init_norm)
  ! print *, mof_niter
  ! Call MOFZY(f,c,norm)
  print *, norm
  block
    integer :: ii
    integer :: nnn = 1
    real(8) :: tt1, tt2
    real(8) :: ttt1, ttt2
    real(8) :: nnn1, nnn2
    Call cpu_time(tt1)
    Do ii = 1, nnn
      init_norm = init_norm + 0.001
      Call MOFLemoine(f,c,norm, init_norm)
      nnn1 = nnn1 + mof_niter
    End Do
    Call cpu_time(tt2)
    print *, "===Lemoine================="
    print *, "Norm", norm
    print *, "CPU time", tt2-tt1
    print *, "Averaged iteration", dble(nnn1/dble(nnn))
    print *, ""
    Call cpu_time(ttt1)
    c = c - 0.5_sp
    Do ii = 1, nnn
      init_norm = init_norm + 0.001
      Call MOFZY(f,c,norm,init_norm)
      nnn2 = nnn2 + mof_niter
    End Do
    Call cpu_time(ttt2)
    print *, "===ZY======================"
    print *, "Norm", norm
    print *, "CPU time", ttt2-ttt1
    print *, "Averaged iteration", dble(nnn2/dble(nnn))
    print *, ""
  End block

  ! Call MOFLemoine(f,c,norm, init_norm)
  ! print *, norm
  ! c = c - 0.5_sp
  ! Call MOFZY(f,c,norm)
  ! print *, norm
  
  Call Normalization2(norm)

  ! if (myid .eq. 0) Print *,norm

  Call MPI_FINALIZE(ierr)


end Subroutine test1

Subroutine test3
  Use ModGlobal
  Use ModTools
  Use ModVOF
  Implicit None
  Real(sp), allocatable, Dimension(:,:,:) :: f_beg
  Real(sp), allocatable, Dimension(:,:,:) :: f_end
  Real(sp), allocatable, Dimension(:,:,:) :: f_exact
  Real(sp) :: v1, v2
  Real(sp) :: v11, v12
  Integer  :: nn = 0
  Character(80) :: data_name
  Integer :: i, j, k
  Real(sp) :: err
  Integer :: rank

  Call Init(inputfield=.true.)

  Call u_bc%SetBCS(u)
  Call v_bc%SetBCS(v)
  Call w_bc%SetBCS(w)
  Call phi_bc%SetBCS(phi)
  Call phi_bc%SetBCS(cx)
  Call phi_bc%SetBCS(cy)
  Call phi_bc%SetBCS(cz)
  Allocate(f_beg(0:nl(1)+1,0:nl(2)+1,0:nl(3)+1))
  Allocate(f_end(0:nl(1)+1,0:nl(2)+1,0:nl(3)+1))
  Allocate(f_exact(0:nl(1)+1,0:nl(2)+1,0:nl(3)+1))
  f_beg = phi

  f_exact = f_beg
  ! f_exact(11:15,6:10,6:10) = 1.0_sp

  ! VOF advection
  Do While (time < tend)
    nn = nn + 1
    ! Call VOFCIAM(Phi, u, v, w, nl, dl, dt)
    ! Call VOFWY(Phi, u, v, w, nl, dl, dt)
    rank = mod(nn+1,3)
    Call VOFAdvection(Phi, u, v, w, nl, dl, dt,rank, cx, cy, cz)
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
    ! call AdvCIAM_MOF(u, cx, cy, cz, phi, nl, dl, dt, 1)
    ! Call Visual3DContour(f1=phi)
    ! call AdvCIAM_MOF(v, cx, cy, cz, phi, nl, dl, dt, 2)
    ! Call Visual3DContour(f1=phi)
    time =  time + dt
  End Do
  print *, nn

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

  err = 0.0_sp
  Do k = 1, nl(3)
    Do j = 1, nl(2)
      Do i = 1, nl(1)
        err = err + abs(f_end(i,j,k)-f_exact(i,j,k))
      End Do
    End Do
  End Do


  v1 = sum(f_beg(1:nl(1),1:nl(2),1:nl(3)))
  v2 = sum(f_end(1:nl(1),1:nl(2),1:nl(3)))

  Call MPI_Reduce(v1, v11, 1, MPI_REAL_SP, MPI_SUM, MPI_COMM_WORLD, 0, ierr)
  Call MPI_Reduce(v1, v12, 1, MPI_REAL_SP, MPI_SUM, MPI_COMM_WORLD, 0, ierr)
  if (myid .eq.0) then
    print *, v11, v12
    print *, err
  endif

  ! Call Visual3DContour(f1=f_end)
  Call Visual3DContour(f1=f_beg, f2=f_end)
  Call Visual2DContour(f1=f_beg, f2=f_end, slice_dir=3, slice_coord=1)


  Call MPI_FINALIZE(ierr)


end Subroutine test3

