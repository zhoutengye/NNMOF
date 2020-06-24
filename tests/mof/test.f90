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
  use mod_cg3_polyhedron
  use MODSussman
  Implicit None
  Real(sp) :: f
  Real(sp) :: c(3)
  Real(sp) :: norm(3)
  Real(sp) :: init_norm(3)
  Real(sp) :: angle(2)
  type(t_polyhedron) :: poly
  Real(sp) :: ddx(3)

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
  f = 2.0/27.0_sp
  c = (/ 1.0/8.0_sp, 1.0/6.0_sp, 1.0/4.0_sp /)
  init_norm = (/ 1.0/3.0_sp, 1.0/3.0_sp, 1.0/4.0_sp /)

  ! f = 1.0/12.0_sp
  ! c = (/ 1.0/4.0_sp, 1.0/4.0_sp, 1.0/8.0_sp /)
  ! init_norm = (/ 1.0/1.0_sp, 1.0/8.0_sp, 1.0/8.0_sp /)
  ! init_norm = (/ 1.0/3.0_sp, 1.0/3.0_sp,  1.0/4.0_sp /)
  f = 1.0/12.0_sp
  c = (/ 3.0/4.0_sp, 3.0/4.0_sp, 7.0/8.0_sp /)

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
      Call MOFLemoine_GaussNewton(f,c,norm)
      ! Call MOFLemoine_GaussNewton(f,c,norm, init_norm)
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
      Call MOFZY(f,c,norm)
      ! Call MOFZY(f,c,norm,init_norm)
      nnn2 = nnn2 + mof_niter
    End Do
    Call cpu_time(ttt2)
    print *, "===ZY======================"
    print *, "Norm", norm
    print *, "CPU time", ttt2-ttt1
    print *, "Averaged iteration", dble(nnn2/dble(nnn))
    print *, ""

    ddx = 1.0_sp
    Call create_cuboid(ddx, poly)

    Call cpu_time(ttt1)
    c = c + 0.5_sp
    nnn1 = 0
    Do ii = 1, nnn
      init_norm = init_norm + 0.001
      Call MOFLemoine_BFGS(poly,f,c,norm)
      ! Call MOFLemoine_BFGS(poly,f,c,norm,init_norm)
      nnn1 = nnn1 + mof_niter
    End Do
    Call cpu_time(ttt2)
    print *, "===Lemoine BFGS================="
    print *, "Norm", norm
    print *, "CPU time", ttt2-ttt1
    print *, "Averaged iteration", dble(nnn1/dble(nnn))
    print *, ""

    Call MOFInit3d
    Call cpu_time(ttt1)
    c = c - 0.5_sp
    nnn1 = 0
    ! sussman%mof_verbose    = 1          ! 0 or 1
    Do ii = 1, nnn
      init_norm = init_norm + 0.001
      Call MOFSussmanGaussNewton(f,c,norm)
      ! Call MOFSussman(f,c,norm,init_norm)
      nnn1 = nnn1 + mof_niter
    End Do
    Call cpu_time(ttt2)
    Call Normalization1(norm)
    print *, "===Sussman GaussNewton================="
    print *, "Norm", norm
    print *, "CPU time", ttt2-ttt1
    print *, "Averaged iteration", dble(nnn1/dble(nnn))
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

subroutine create_cuboid(c, poly)
  use mod_cg3_polyhedron
  use mod_cg3_complete_polyhedron_structure
  double precision, dimension(3), intent(in) :: c
  type(t_polyhedron), intent(out) :: poly
  integer :: error_id

  poly%nb_points = 8
  poly%nb_faces = 6

  allocate(poly%point(3, poly%nb_points))

  poly%point(:,1) = [0.0d0, 0.0d0, 0.0d0]
  poly%point(:,2) = [0.0d0, 0.0d0, c(3) ]
  poly%point(:,3) = [0.0d0, c(2) , 0.0d0]
  poly%point(:,4) = [0.0d0, c(2) , c(3) ]
  poly%point(:,5) = [c(1) , 0.0d0, 0.0d0]
  poly%point(:,6) = [c(1) , 0.0d0, c(3) ]
  poly%point(:,7) = [c(1) , c(2) , 0.0d0]
  poly%point(:,8) = [c(1) , c(2) , c(3) ]

  allocate(poly%face(poly%nb_faces))

  poly%face(1)%size = 4
  allocate(poly%face(1)%id(poly%face(1)%size))
  poly%face(1)%id = [8, 4, 2, 6]

  poly%face(2)%size = 4
  allocate(poly%face(2)%id(poly%face(2)%size))
  poly%face(2)%id= [8, 6, 5, 7]

  poly%face(3)%size = 4
  allocate(poly%face(3)%id(poly%face(3)%size))
  poly%face(3)%id = [8, 7, 3, 4]

  poly%face(4)%size = 4
  allocate(poly%face(4)%id(poly%face(4)%size))
  poly%face(4)%id = [4, 3, 1, 2]

  poly%face(5)%size = 4
  allocate(poly%face(5)%id(poly%face(5)%size))
  poly%face(5)%id = [1, 3, 7, 5]

  poly%face(6)%size = 4
  allocate(poly%face(6)%id(poly%face(6)%size))
  poly%face(6)%id = [2, 1, 5, 6]

  call cg3_complete_polyhedron_structure(poly, error_id)
end subroutine create_cuboid

