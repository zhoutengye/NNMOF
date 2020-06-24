Module ModSussman
  ! from Sussman MOF
  Use probcommon_module
  Use MOF_routines_module
  Implicit None

  Type varsussman
    Integer :: caller_id
    Integer :: bfact
    Integer :: nhalf0
    Integer :: use_ls_data
    Integer :: mof_verbose
    Integer :: nmat
    Integer :: sdim
    Integer :: nmax
    Integer :: continuous_mof

    ! *** These two variables are defined in module probcommon_module
    ! Integer :: levelrz
    ! Integer :: ngeom_recon

    Integer :: order_algorithm_in     ! all components equal to nmat+1
    Integer :: MOFITERMAX_in 
    Integer :: MOF_DEBUG_RECON_in 
    Integer :: nthreads 
    Integer :: MOF_TURN_OFF_LS_in 
    Integer :: nmax_in

    Real(8), Allocatable :: LS_stencil(:,:,:,:)

    Real(8), Allocatable :: dx(:)

    Real(8), Allocatable :: xtetlist_vof(:,:,:)
    Real(8), Allocatable :: xtetlist_cen(:,:,:)
    Real(8), Allocatable :: xsten0(:,:)
    Real(8), Allocatable :: mofdata(:)
    Real(8), Allocatable :: multi_centroidA(:,:)

    Integer :: Linear_recon_FLAG !  1: succesfull reconstruction 0: invalid line information 
  End type varsussman

  Type(varsussman) :: sussman

Contains
  Subroutine MOFInit3d
    Implicit None

    sussman%caller_id      = 0
    sussman%bfact          = 1
    sussman%nhalf0         = 3
    sussman%use_ls_data    = 0          ! LS_stencil -> not used if use_ls_data=0
    sussman%nmax           = 1000       ! sufficient large
    sussman%mof_verbose    = 0          ! 0 or 1
    sussman%sdim           = 3          ! Number of imensions
    sussman%nmat           = 1          ! Numver of materials
    ngeom_recon    = 0
    sussman%continuous_mof = 0
    levelrz        = 0
    sussman%order_algorithm_in = 1      ! all components equal to nmat+1
    sussman%MOFITERMAX_in      = 10     ! MOFITERMAX_in=10 is a reasonable choice
    sussman%MOF_DEBUG_RECON_in = 10     ! Integer 
    sussman%nthreads           = 1      ! nthreads=1
    sussman%MOF_TURN_OFF_LS_in = 1      ! MOF_TURN_OFF_LS_in=1 is a reasonable choice
    sussman%nmax_in            = 500    ! nmax=500 is a reasonable choice.

    ngeom_recon = 2 * sussman%sdim + 3
    num_materials = sussman%nmat
    sussman%nmax_in = sussman%nmax

    Allocate(sussman%LS_stencil(-1:1,-1:1,-1:1,sussman%nmat))
    Allocate(sussman%dx(sussman%sdim))
    Allocate(sussman%xtetlist_vof(sussman%sdim+1,sussman%sdim,sussman%nmax))
    Allocate(sussman%xtetlist_cen(sussman%sdim+1,sussman%sdim,sussman%nmax)) 
    Allocate(sussman%xsten0(-sussman%nhalf0:sussman%nhalf0,sussman%sdim))
    Allocate(sussman%mofdata(sussman%nmat*ngeom_recon))
    Allocate(sussman%multi_centroidA(sussman%nmat,sussman%sdim))


    sussman%LS_stencil      = 0.d0
    sussman%dx              = 1.d0
    sussman%xtetlist_vof    = 0.d0
    sussman%xtetlist_cen    = 0.d0
    sussman%xsten0          = 0.d0
    sussman%mofdata         = 1.d0
    sussman%multi_centroidA = 0.d0

    sussman%xsten0(-3,1) = -1.5d0
    sussman%xsten0(-2,1) = -1.0d0
    sussman%xsten0(-1,1) = -0.5d0
    sussman%xsten0(0,1)  = 0.d0
    sussman%xsten0(1,1)  = 0.5d0
    sussman%xsten0(2,1)  = 1.0d0
    sussman%xsten0(3,1)  = 1.5d0

    sussman%xsten0(-3,2) = -1.5d0
    sussman%xsten0(-2,2) = -1.0d0
    sussman%xsten0(-1,2) = -0.5d0
    sussman%xsten0(0,2)  = 0.d0
    sussman%xsten0(1,2)  = 0.5d0
    sussman%xsten0(2,2)  = 1.0d0
    sussman%xsten0(3,2)  = 1.5d0

    sussman%xsten0(-3,3) = -1.5d0
    sussman%xsten0(-2,3) = -1.0d0
    sussman%xsten0(-1,3) = -0.5d0
    sussman%xsten0(0,3)  = 0.d0
    sussman%xsten0(1,3)  = 0.5d0
    sussman%xsten0(2,3)  = 1.0d0
    sussman%xsten0(3,3)  = 1.5d0

    Call initmof( &
        sussman%order_algorithm_in, &
        sussman%nmat,&
        sussman%MOFITERMAX_in, &
        sussman%MOF_DEBUG_RECON_in, &
        sussman%MOF_TURN_OFF_LS_in, &
        sussman%nthreads, &
        sussman%nmax_in)

  End Subroutine MOFInit3d

End Module ModSussman
