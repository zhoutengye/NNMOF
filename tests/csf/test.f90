#include "param.h"
# if defined TEST1
# define test_case test1
# elif defined TEST2
# define test_case test2
# elif defined TEST3
# define test_case test3
# endif
program test
  use ModGlobal
  Implicit None

  call test_case

end program test

!===============
! test1:
! (1) Gen_and_save_data
! (1) Check single processor VOF to LS
!===============
Subroutine test1
  Use ModGlobal
  Use ModTools
  Use ModVOFFunc
  Use ModNavierStokes
  Use InitSH
  Implicit None
  Integer :: i
  real(sp), allocatable :: Phi_DS(:,:,:)
  real(sp), allocatable :: LS(:,:,:)
  Character(80) :: data_name

  Call Init(inputfield=.false.)
  Call InitNavierStokes(U, V, W, Phi, P)
  call InitVolume
  Allocate(Phi_DS(0:nl(1)+1,0:nl(2)+1,0:nl(3)+1))
  Allocate(LS(0:nl(1)+1,0:nl(2)+1,0:nl(3)+1))
  Call Smooth_Phi(Phi, Phi_DS, nl(1), nl(2), nl(3))
  Phi = Phi_DS
  Call Phi_bc%SetBCS(Phi)
  Call Smooth_Phi(Phi, Phi_DS, nl(1), nl(2), nl(3))
  ! Call Visual3DContour(LS)
  Call Visual3DContour(Phi_DS)
  Call Get_LS_3D_init(Phi_DS, Ls, dl(1), dl(2), dl(3), nl(1), nl(2), nl(3))
  ! Call Visual3DContour(LS,slice_dir=3,slice_coord=25)
  Call Visual3DContour(LS)

  data_name = 'init'
  Call HDF5WriteFrame(data_name)

  Call Finalize

end Subroutine test1

!===============
! test2:
! (1) Test MPI VOF to LS
!===============
Subroutine test2
  Use ModGlobal
  Use ModTools
  Use ModVOFFunc
  Use ModNavierStokes
  Implicit None
  Integer :: i, j, k
  Real(sp) :: vol1, vol2, vol11, vol21
  real(sp), allocatable :: Phi_DS(:,:,:)
  real(sp), allocatable :: LS(:,:,:)


  Call Init(inputfield=.true.)
  Call InitNavierStokes(U, V, W, Phi, P)

  Allocate(Phi_DS(0:nl(1)+1,0:nl(2)+1,0:nl(3)+1))
  Allocate(LS(0:nl(1)+1,0:nl(2)+1,0:nl(3)+1))
  Call Smooth_Phi(Phi, Phi_DS, nl(1), nl(2), nl(3))
  Phi = Phi_DS
  Call Phi_bc%SetBCS(Phi)
  Call Smooth_Phi(Phi, Phi_DS, nl(1), nl(2), nl(3))
  ! Call Visual3DContour(LS)
  Call Visual3DContour(Phi_DS)
  Call Get_LS_3D_init(Phi_DS, Ls, dl(1), dl(2), dl(3), nl(1), nl(2), nl(3))
  ! Call Visual3DContour(LS,slice_dir=3,slice_coord=25)
  Call Visual3DContour(LS)

  ! data_name = 'init'
  ! Call HDF5WriteFrame(data_name)

  Call Finalize

End Subroutine Test2

!===============
! test3:
! (1) Test vector field and pressure
!===============
Subroutine test3
  Use ModGlobal
  Use ModTools
  Use ModNavierStokes
  Implicit None
  real(sp), allocatable :: dudt(:,:,:)
  real(sp), allocatable :: dvdt(:,:,:)
  real(sp), allocatable :: dwdt(:,:,:)
  real(sp), allocatable :: div(:,:,:)


  Call Init(inputfield=.true.)
  Call InitNavierStokes(U, V, W, Phi, P)

  CSF => CSF_VOSET
  Call TwoPhaseFlow(U, V, W, Phi, P, cx, cy, cz)

  Call Visual3DContour(P)
  Call Visual2DContour(P,slice_dir=3,slice_coord=25)
  Call Visual3DQuiver(V,U,W)

  Call Finalize

End Subroutine Test3

