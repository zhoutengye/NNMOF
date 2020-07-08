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
  real(sp), allocatable :: Phi_DS(:,:,:)
  real(sp), allocatable :: LS(:,:,:)
  Character(80) :: data_name
  Integer :: i,j,k

  Call Init(inputfield=.false.)
  Call InitNavierStokes(U, V, W, Phi, P)
  ShapeLevelSet => ShapeDropImpact
  call InitVolume
  Call Phi_bc%SetBCS(Phi)
  Allocate(Phi_DS(0:nl(1)+1,0:nl(2)+1,0:nl(3)+1))
  Allocate(LS(0:nl(1)+1,0:nl(2)+1,0:nl(3)+1))
  ! Call Smooth_Phi(Phi, Phi_DS, nl(1), nl(2), nl(3))
  ! Phi = Phi_DS
  ! Call Visual3DContour(LS)
  Call Visual3DContour(Phi)
  Call Get_LS_3D_init(Phi_DS, Ls, dl(1), dl(2), dl(3), nl(1), nl(2), nl(3))
  ! Call Visual3DContour(LS,slice_dir=3,slice_coord=25)
  ! Call Visual3DContour(LS)

  ! print *, Phi(50,50:51,25)
  ! print *, ls(50,41:60,25)
  ! print *, ls(50,31:70,25)
  Do k = 1, nl(3)
    Do j = 1, nl(2)
      Do i = 1, nl(1)
        if (k>11 .and. Phi(i,j,k) > 0.01_sp) w(i,j,k) = -1.0_sp
      End Do
    End Do
  End Do

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
  Integer :: i
  real(sp), allocatable :: dudt(:,:,:)
  real(sp), allocatable :: dvdt(:,:,:)
  real(sp), allocatable :: dwdt(:,:,:)
  real(sp), allocatable :: div(:,:,:)


  Call Init(inputfield=.true.)
  Call InitNavierStokes(U, V, W, Phi, P)

  ! U = 0.0_sp
  ! V = 0.0_sp
  ! W = 0.0_sp
  ! P = 0.0_sp
  ! Phi = 0.0_sp
  ! Phi(41:60,41:60,16:35)= 1.0_sp
  CSF => CSF_VOSET
  Do i = 1,10
    Call TwoPhaseFlow(U, V, W, Phi, P)
    print *, n_iter
  End Do


  ! Call Visual2DContour(P,slice_dir=3,slice_coord=25)
  Call Visual3DQuiver(V,U,W)
  Call Visual3Dcontour(P)

  Call Finalize

End Subroutine Test3

Subroutine test4
  Use ModGlobal
  Use ModTools
  Use ModNavierStokes
  Implicit None
  Integer :: i, j, k
  real(sp), allocatable :: dudt(:,:,:)
  real(sp), allocatable :: dvdt(:,:,:)
  real(sp), allocatable :: dwdt(:,:,:)
  real(sp), allocatable :: div(:,:,:)


  Call Init(inputfield=.true.)

  U = 0.0_sp
  V = 0.0_sp
  W = 0.0_sp
  ! P = 0.0_sp
  ! Phi = 0.0_sp
  ! Phi(:,:,0:10)= 1.0_sp
  Call InitNavierStokes(U, V, W, Phi, P)
  ! CSF => CSF_VOSET
  ! Do k = 1, nl(3)
  !   Do j = 1, nl(2)
  !     Do i = 1, nl(1)
  !       if (k>5 .and. Phi(i,j,k) > 0.01_sp) W(i,j,k) = -1.0_sp
  !     End Do
  !   End Do
  ! End Do
  Call BC_UVW(U, V, W)
  Call Visual3DContour(Phi)
  Do While (time < tend)
    Call TwoPhaseFlow(U, V, W, Phi, P)
    if (myid .eq. 0 ) print *, time, n_iter
    Call WriteFieldData
    time = time + dt
  End Do
  Call Visual3Dcontour(Phi)
  ! Call Visual3Dquiver(V,U,W)

  Call Finalize

End Subroutine Test4

