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
  Implicit None

  call test_case

end program test

Subroutine test_case
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
  Real(sp) :: tt1, tt2
  Integer :: rank
  Real(sp) :: error_rr
  Real(sp) :: error_r
  Real(sp) :: error_g
  Real(sp) :: error_m

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
  Call Visual3DContour(f1=f_beg)


  MOFNorm => MOFZY
  ! VOF advection
  mof_tol = 1.0e-8
  mofitermax = 10
  delta_theta = 0.1_sp * MOF_Pi / 180.0_sp  ! 10 degrees
  delta_theta = 1e-8
  delta_theta_max = 10.0_sp * MOF_Pi / 180.0_sp  ! 10 degrees
  Call CPU_Time(tt1)
  Do While (time < tend)
    nn = nn + 1
    ! Call VOFWY(Phi, u, v, w, nl, dl, dt)
    ! Call VOFCIAM(Phi, u, v, w, nl, dl, dt)
    if (myid .eq. 0) print *, 'step =', nn, 'mof_ier', mof_niter(1)
    rank = mod(nn+1,3)
    Call MOFCIAM(Phi, u, v, w, nl, dl, dt,rank, cx, cy, cz)
    ! Call VOFCIAM(Phi, u, v, w, nl, dl, dt,rank)
    time =  time + dt
  End Do
  Call CPU_Time(tt2)
  if (myid .eq.0 ) print *, nn

    data_name = 'init'
    Call HDF5WriteFrame(data_name)

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

  Do k = 1, nl(3)
    Do j = 1, nl(2)
      Do i = 1, nl(1)
        error_rr = error_rr + abs(f_beg(i,j,k) - f_end(i,j,k))
      End Do
    End Do
  End Do

  Call MPI_Reduce(v1, v11, 1, MPI_REAL_SP, MPI_SUM, MPI_COMM_WORLD, 0, ierr)
  Call MPI_Reduce(v1, v12, 1, MPI_REAL_SP, MPI_SUM, MPI_COMM_WORLD, 0, ierr)
  Call MPI_Reduce(error_rr, error_r, 1, MPI_REAL_SP, MPI_SUM, MPI_COMM_WORLD, 0, ierr)

  error_m = (v11 - v12) / v11
  error_g = (error_r) / dble(n(1)) / dble(n(2)) / dble(n(3))
  error_r = (error_r) / v11

  if (myid .eq.0) then
    print *, 'cpu_time = ', tt2-tt1
    print *, 'realtive distortion error = ', error_r
    print *, 'absolute error = ', error_g
    print *, 'conservation error = ', error_m
    open(10,file='errors.dat')
    write(10,*) , 'cpu_time = ', tt2-tt1
    write(10,*) , 'realtive distortion error = ', error_r
    write(10,*) , 'absolute error = ', error_g
    write(10,*) , 'conservation error = ', error_m
    close(10)
  endif

  ! Call Visual3DContour(f1=f_end)
  ! Call Visual3DContour(f1=f_beg, f2=f_end)
  Call Visual3DContour(f1=f_end)
  Call Visual2DContour(f1=f_beg, f2=f_end, slice_dir=3, slice_coord=25)

  Call MPI_FINALIZE(ierr)


end Subroutine test_case

