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
! test1: single cell
!===============
Subroutine test1
  ! Use ModGlobal
  Use ModTools
  Implicit None

  ! Call Init(inputfield=.false.)

  Type(Octree) :: gr

  gr%level = 1
  gr%maxlevel = 2
  TreeChild%xc(1) = 0.5_sp
  TreeChild%xc(2) = 0.5_sp
  TreeChild%xc(3) = 0.5_sp
  TreeChild%dx(1) = 1.0_sp
  TreeChild%dx(2) = 1.0_sp
  TreeChild%dx(3) = 1.0_sp
  ShapeLevelSet => ShapeTest

  Call VolumeCentroidQuadTree(gr)

  print *, gr%vof
  print *, gr%centroid
  print *, 'tt'

  Call MPI_FINALIZE(ierr)

end Subroutine test1

Real(sp) Function ShapeTest(x,y,z)
  Implicit None
  Real(sp) :: x, y, z
  ShapeTest = x+y+z-0.5_sp
End Function ShapeTest

! !===============
! ! test1: MOF reconstruction
! ! Given centroid and volume and obtains the normal vecor
! !===============
! Subroutine test2
!   Use ModGlobal
!   Use ModTools
!   Use ModVOF
!   Implicit None
!   Real(sp), allocatable, Dimension(:,:,:) :: f_beg
!   Real(sp), allocatable, Dimension(:,:,:) :: f_end
!   Real(sp), allocatable, Dimension(:,:,:) :: f_exact
!   Real(sp) :: v1, v2
!   Real(sp) :: v11, v12
!   Integer  :: nexch(2)
!   Integer  :: nn = 0
!   Character(80) :: data_name
!   Integer :: i, j, k
!   Real(sp) :: err
!   Integer :: rank

!   Call Init(inputfield=.true.)

!   u = 1.0_sp
!   v = 1.0_sp
!   w = 1.0_sp

!   nexch = nl(1:2)
!   Call u_bc%SetBCS(u)
!   Call v_bc%SetBCS(v)
!   Call w_bc%SetBCS(w)
!   Call phi_bc%SetBCS(phi)
!   Call phi_bc%SetBCS(cx)
!   Call phi_bc%SetBCS(cy)
!   Call phi_bc%SetBCS(cz)
!   Allocate(f_beg(0:nl(1)+1,0:nl(2)+1,0:nl(3)+1))
!   Allocate(f_end(0:nl(1)+1,0:nl(2)+1,0:nl(3)+1))
!   Allocate(f_exact(0:nl(1)+1,0:nl(2)+1,0:nl(3)+1))
!   f_beg = phi

!   f_exact = 0.0_sp
!   ! f_exact(11:15,6:10,6:10) = 1.0_sp
!   ! f_exact(11:15,11:15,11:15) = 1.0_sp
!   f_exact = phi

!   ! VOF advection
!   Do While (time < tend-dt)
!     nn = nn + 1
!     ! Call VOFCIAM(Phi, u, v, w, nl, dl, dt)
!     ! Call VOFWY(Phi, u, v, w, nl, dl, dt)
!     rank = mod(nn+1,3)
!     ! Call MOFCIAM(Phi, cx, cy, cz, u, v, w, nl, dl, dt)
!     Call MOFCIAM2(Phi, cx, cy, cz, u, v, w, nl, dl, dt,rank)
!     ! Call MOFWY(Phi, cx, cy, cz, u, v, w, nl, dl, dt)
!     ! call AdvWY_MOF(u, cx, cy, cz, phi, nl, dl, dt, 1)
!     ! Call Visual3DContour(f1=phi, slice_dir='x',slice_coord=8)
!     ! call AdvWY_MOF(v, cx, cy, cz, phi, nl, dl, dt, 2)
!     ! Call Visual3DContour(f1=phi, slice_dir='x',slice_coord=8)
!     ! Call Visual3DContour(f1=phi)
!     ! call AdvWY_MOF(w, cx, cy, cz, phi, nl, dl, dt, 3)
!     ! Call Visual3DContour(f1=phi)
!     ! call AdvCIAM_MOF(w, cx, cy, cz, phi, nl, dl, dt, 3)
!     ! Call Visual3DContour(f1=phi)
!     ! call AdvCIAM_MOF(v, cx, cy, cz, phi, nl, dl, dt, 2)
!     ! Call Visual3DContour(f1=phi)
!     ! call AdvCIAM_MOF(w, cx, cy, cz, phi, nl, dl, dt, 3)
!     ! Call Visual3DContour(f1=phi)
!     ! call AdvCIAM_MOF(u, cx, cy, cz, phi, nl, dl, dt, 1)
!     ! Call Visual3DContour(f1=phi)
!     ! call AdvCIAM_MOF(v, cx, cy, cz, phi, nl, dl, dt, 2)
!     ! Call Visual3DContour(f1=phi)
!     time =  time + dt
!   End Do
!   print *, nn
!   Call Visual3DContour(f1=f_beg)

!     data_name = 'test'
!   do nn = 1, n_vars
!       Select Case(Trim(h5_output_field(nn)%groupname))
!       Case('phi')
!         Call HDF5WriteData(h5_output_field(nn), phi,data_name)
!       Case('u')
!         Call HDF5WriteData(h5_output_field(nn), u,data_name)
!       Case('v')
!         Call HDF5WriteData(h5_output_field(nn), v,data_name)
!       Case('w')
!         Call HDF5WriteData(h5_output_field(nn), w,data_name)
!       Case('cx')
!         Call HDF5WriteData(h5_output_field(nn), cx,data_name)
!       Case('cy')
!         Call HDF5WriteData(h5_output_field(nn), cy,data_name)
!       Case('cz')
!         Call HDF5WriteData(h5_output_field(nn), cz,data_name)
!       End Select
!   end do

!   f_end = phi

!   err = 0.0_sp
!   Do k = 1, nl(3)
!     Do j = 1, nl(2)
!       Do i = 1, nl(1)
!         err = err + abs(f_end(i,j,k)-f_exact(i,j,k))
!       End Do
!     End Do
!   End Do


!   v1 = sum(f_beg(1:nl(1),1:nl(2),1:nl(3)))
!   v2 = sum(f_end(1:nl(1),1:nl(2),1:nl(3)))

!   Call MPI_Reduce(v1, v11, 1, MPI_REAL_SP, MPI_SUM, MPI_COMM_WORLD, 0, ierr)
!   Call MPI_Reduce(v1, v12, 1, MPI_REAL_SP, MPI_SUM, MPI_COMM_WORLD, 0, ierr)
!   if (myid .eq.0) then
!     print *, v11, v12
!     print *, err
!   endif

!   ! Call Visual3DContour(f1=f_end)
!   ! Call Visual3DContour(f1=f_beg, f2=f_end)
!   ! Call Visual2DContour(f1=f_beg, f2=f_end, slice_dir=2, slice_coord=10)
!   Call Visual3DContour(f1=f_end)
!   Call Visual2DContour(f1=f_end, slice_dir=2, slice_coord=10)


!   Call MPI_FINALIZE(ierr)


! end Subroutine test2

! Subroutine test3
!   Use ModGlobal
!   Use ModTools
!   Use ModVOF
!   Implicit None
!   Real(sp), allocatable, Dimension(:,:,:) :: f_beg
!   Real(sp), allocatable, Dimension(:,:,:) :: f_end
!   Real(sp), allocatable, Dimension(:,:,:) :: f_exact
!   Real(sp) :: v1, v2
!   Real(sp) :: v11, v12
!   Integer  :: nexch(2)
!   Integer  :: nn = 0
!   Character(80) :: data_name
!   Integer :: i, j, k
!   Real(sp) :: err
!   Integer :: rank

!   Call Init(inputfield=.true.)

!   nexch = nl(1:2)
!   Call u_bc%SetBCS(u)
!   Call v_bc%SetBCS(v)
!   Call w_bc%SetBCS(w)
!   Call phi_bc%SetBCS(phi)
!   Call phi_bc%SetBCS(cx)
!   Call phi_bc%SetBCS(cy)
!   Call phi_bc%SetBCS(cz)
!   Allocate(f_beg(0:nl(1)+1,0:nl(2)+1,0:nl(3)+1))
!   Allocate(f_end(0:nl(1)+1,0:nl(2)+1,0:nl(3)+1))
!   Allocate(f_exact(0:nl(1)+1,0:nl(2)+1,0:nl(3)+1))
!   f_beg = phi

!   f_exact = f_beg
!   ! f_exact(11:15,6:10,6:10) = 1.0_sp

!   ! VOF advection
!   Do While (time < tend)
!     nn = nn + 1
!     ! Call VOFCIAM(Phi, u, v, w, nl, dl, dt)
!     ! Call VOFWY(Phi, u, v, w, nl, dl, dt)
!     rank = mod(nn+1,3)
!     ! Call MOFCIAM(Phi, cx, cy, cz, u, v, w, nl, dl, dt)
!     Call MOFCIAM2(Phi, cx, cy, cz, u, v, w, nl, dl, dt,rank)
!     ! Call MOFWY(Phi, cx, cy, cz, u, v, w, nl, dl, dt)
!     ! call AdvWY_MOF(u, cx, cy, cz, phi, nl, dl, dt, 1)
!     ! Call Visual3DContour(f1=phi, slice_dir='x',slice_coord=8)
!     ! call AdvWY_MOF(v, cx, cy, cz, phi, nl, dl, dt, 2)
!     ! Call Visual3DContour(f1=phi, slice_dir='x',slice_coord=8)
!     ! Call Visual3DContour(f1=phi)
!     ! call AdvWY_MOF(w, cx, cy, cz, phi, nl, dl, dt, 3)
!     ! Call Visual3DContour(f1=phi)
!     ! call AdvCIAM_MOF(w, cx, cy, cz, phi, nl, dl, dt, 3)
!     ! Call Visual3DContour(f1=phi)
!     ! call AdvCIAM_MOF(v, cx, cy, cz, phi, nl, dl, dt, 2)
!     ! Call Visual3DContour(f1=phi)
!     ! call AdvCIAM_MOF(w, cx, cy, cz, phi, nl, dl, dt, 3)
!     ! Call Visual3DContour(f1=phi)
!     ! call AdvCIAM_MOF(u, cx, cy, cz, phi, nl, dl, dt, 1)
!     ! Call Visual3DContour(f1=phi)
!     ! call AdvCIAM_MOF(v, cx, cy, cz, phi, nl, dl, dt, 2)
!     ! Call Visual3DContour(f1=phi)
!     time =  time + dt
!   End Do
!   print *, nn

!     data_name = 'test'
!   do nn = 1, n_vars
!       Select Case(Trim(h5_output_field(nn)%groupname))
!       Case('phi')
!         Call HDF5WriteData(h5_output_field(nn), phi,data_name)
!       Case('u')
!         Call HDF5WriteData(h5_output_field(nn), u,data_name)
!       Case('v')
!         Call HDF5WriteData(h5_output_field(nn), v,data_name)
!       Case('w')
!         Call HDF5WriteData(h5_output_field(nn), w,data_name)
!       Case('cx')
!         Call HDF5WriteData(h5_output_field(nn), cx,data_name)
!       Case('cy')
!         Call HDF5WriteData(h5_output_field(nn), cy,data_name)
!       Case('cz')
!         Call HDF5WriteData(h5_output_field(nn), cz,data_name)
!       End Select
!   end do

!   f_end = phi

!   err = 0.0_sp
!   Do k = 1, nl(3)
!     Do j = 1, nl(2)
!       Do i = 1, nl(1)
!         err = err + abs(f_end(i,j,k)-f_exact(i,j,k))
!       End Do
!     End Do
!   End Do


!   v1 = sum(f_beg(1:nl(1),1:nl(2),1:nl(3)))
!   v2 = sum(f_end(1:nl(1),1:nl(2),1:nl(3)))

!   Call MPI_Reduce(v1, v11, 1, MPI_REAL_SP, MPI_SUM, MPI_COMM_WORLD, 0, ierr)
!   Call MPI_Reduce(v1, v12, 1, MPI_REAL_SP, MPI_SUM, MPI_COMM_WORLD, 0, ierr)
!   if (myid .eq.0) then
!     print *, v11, v12
!     print *, err
!   endif

!   ! Call Visual3DContour(f1=f_end)
!   Call Visual3DContour(f1=f_beg, f2=f_end)
!   Call Visual2DContour(f1=f_beg, f2=f_end, slice_dir=3, slice_coord=1)


!   Call MPI_FINALIZE(ierr)


! end Subroutine test3

