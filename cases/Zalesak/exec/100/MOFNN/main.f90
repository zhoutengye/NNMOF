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
  use variables_mof

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
  Real(sp) :: error_rr
  Real(sp) :: error_r
  Real(sp) :: error_g
  Real(sp) :: error_m
  Integer :: sum_iter(2) = 0

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

  u = - u
  v = - v
  w = - w
  ! u = u
  ! v = v
  ! w = w

  f_beg = phi
  f_exact = f_beg

 
  Call Initialize_NN
  MOFNorm => MOFNN

  !!! =====Method 5, My numerical gradient, Gauss Newton=========
  ! MOFNorm => MOFZY

  ! VOF advection
  Call CPU_Time(tt1)
  Do While (time < tend)
    if (myid .eq. 0) print *, 'step =', nn
    rank = mod(nn+1,6)
    ! Call MOFWY(Phi, u, v, w, nl, dl, dt,rank, cx, cy, cz)
    Call MOFCIAM(Phi, u, v, w, nl, dl, dt,rank, cx, cy, cz)
    sum_iter = sum_iter+mof_niter
    ! Call MOF_EI_LE(Phi, u, v, w, nl, dl, dt,rank, cx, cy, cz)
    ! Call VOFWY(Phi, u, v, w, nl, dl, dt,rank)
    ! Call VOFCIAM(Phi, u, v, w, nl, dl, dt,rank)
    ! Call VOFTHINC(Phi, u, v, w, nl, dl, dt,rank)
    ! Call VOF_EI_LE(Phi, u, v, w, nl, dl, dt,rank)
    nn = nn + 1
    time =  time + dt
    Call WriteFieldData
  End Do
  Call CPU_Time(tt2)

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
    print *, 'iter = ', sum_iter
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

  Call Visual3DContour(f1=f_end)
  Call Visual2DContour(f1=f_end, slice_dir=3, slice_coord=nl(3)/2)

  Call Finalize()


end Subroutine zalesak

