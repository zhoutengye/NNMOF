# include "param.h"
program main
  use ModGlobal
  Implicit None

  call zalesak

end program main

Subroutine zalesak
  Use ModGlobal
  Use ModTools
  Use ModVOF
  use MODSussman

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
  Real(sp) :: ddx(3)

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

  u = -u
  v = -v
  w = -w

  f_beg = phi
  f_exact = f_beg

 
  ! Call MOFInit3d
  ! MOFNorm => MOFSussmanGaussNewton

  ! MOFNorm => MOFLemoine_BFGS
  ! ddx = 1.0_sp
  ! Call Lemoine_create_cuboid(ddx, LemoinePoly)

  MOFNorm => MOFLemoine_GaussNewton

  ! MOFNorm => MOFZY

  ! VOF advection
  Call CPU_Time(tt1)
  Do While (time < tend)
    if (myid .eq. 0) print *, 'step =', nn
    rank = mod(nn+1,3)
    Call MOFWY(Phi, u, v, w, nl, dl, dt,rank, cx, cy, cz)
    ! Call VOFCIAM(Phi, u, v, w, nl, dl, dt,rank)
    nn = nn + 1
    time =  time + dt
  End Do
  Call CPU_Time(tt2)

  data_name = 'final'
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

  Call MPI_Reduce(v1, v11, 1, MPI_REAL_SP, MPI_SUM, MPI_COMM_WORLD, 0, ierr)
  Call MPI_Reduce(v1, v12, 1, MPI_REAL_SP, MPI_SUM, MPI_COMM_WORLD, 0, ierr)
  if (myid .eq.0) then
    print *, 'cpu_time =', tt2-tt1
    print *, 'Initial volume:', v11
    print *, 'Initial volume:', v12
    print *, err
  endif

  Call Visual3DContour(f1=f_end)
  Call Visual2DContour(f1=f_end, slice_dir=3, slice_coord=nl(3)/2)


  Call MPI_FINALIZE(ierr)


end Subroutine zalesak

