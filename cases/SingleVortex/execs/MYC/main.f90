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
  Real(sp), Allocatable :: xx(:), yy(:), zz(:)
  Real(sp) :: dxx
  Real(sp) :: tt
  Real(sp) :: tout, dout
  Real(sp) :: x_start, y_start

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

  Allocate(xx(nl(1)))
  Allocate(yy(nl(3)))
  Allocate(zz(nl(3)))

  print *, nl


  TT = 4.0_sp
  tout = 0.1_sp
  dxx = 1.0_sp / nl(1)
  x_start = dble(nl(1)) * dble(coord(2)) * dxx
  y_start = dble(nl(2)) * dble(coord(1)) * dxx
  Do i = 1, nl(1)
    xx(i) = x_start + dble(i-1) * dxx
  End Do
  Do j = 1, nl(2)
    yy(j) = y_start + dble(j-1) * dxx
  End Do
  Do k = 1, nl(3)
    zz(k) = dble(k-1) * dxx
  End Do

  tout = 0.0_sp
  dout = 0.5_sp

  ! VOF advection
  Call CPU_Time(tt1)
  Do While (time < tend)

    do k=1,nl(3)
      do j=1,nl(2)
        do i=1,nl(1)
          u(i,j,k)=2.0_sp*(sin(Pi*(xx(i)-0.5_sp))**2) &
              * sin(2.0_sp*Pi*xx(j)) &
              * sin(2*Pi*xx(k)) &
              * cos(Pi*time/TT)
          v(i,j,k)= -sin(2.0_sp*Pi*xx(i)) &
              * sin(Pi*xx(j))**2 &
              * sin(2*Pi*xx(k)) &
              * cos(Pi*time/TT)
          w(i,j,k)= -sin(2.0_sp*Pi*xx(i)) &
              * sin(2.0_sp*Pi*xx(j)) &
              * (sin(Pi*xx(k))**2) &
              * cos(Pi*time/TT)
        enddo
      enddo
    enddo
    Call u_bc%SetBCS(u)
    Call v_bc%SetBCS(v)
    Call w_bc%SetBCS(w)
    if (myid .eq. 0) print *, 'step =', nn
    rank = mod(nn+1,3)
    Call VOFWY(Phi, u, v, w, nl, dl, dt)
    nn = nn + 1
    time =  time + dt
    ! If (time .ge. tout) Then
      ! Call Visual3DContour(f1=phi)
      ! tout = tout + dout
    ! End If
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

