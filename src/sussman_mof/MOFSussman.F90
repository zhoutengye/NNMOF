!================================================================================
!
! Machine Learning module for MOF
!
! 
!
!
!
!
!
!=================================================================================
! 
! By Zhouteng Ye at FSU
! Last update: August 24th, 2018
!
!=================================================================================

Module SussmanMOF

    Use probcommon_module
    Implicit None
    Save

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

    Real(8), Allocatable :: xtetlist_vof(:,:,:)
    Real(8), Allocatable :: xtetlist_cen(:,:,:)
    Real(8), Allocatable :: xsten0(:,:)
    Real(8), Allocatable :: mofdata(:)
    Real(8), Allocatable :: multi_centroidA(:,:)


    Integer :: Linear_recon_FLAG !  1: succesfull reconstruction 0: invalid line information 

Contains

    Subroutine Initialize_SussmanMOF

        Implicit None

        caller_id      = 0
        bfact          = 1
        nhalf0         = 3
        use_ls_data    = 1          ! LS_stencil -> not used if use_ls_data=0
        nmax           = 1000       ! sufficient large
        mof_verbose    = 1          ! 0 or 1
        sdim           = 3          ! Number of imensions
        nmat           = 2          ! Numver of materials
        ngeom_recon    = 0
        continuous_mof = 0
        levelrz        = 0

        order_algorithm_in = 1      ! all components equal to nmat+1
        MOFITERMAX_in      = 200     ! MOFITERMAX_in=10 is a reasonable choice
        MOF_DEBUG_RECON_in = 10     ! Integer 
        nthreads           = 1      ! nthreads=1
        MOF_TURN_OFF_LS_in = 0      ! MOF_TURN_OFF_LS_in=1 is a reasonable choice
        nmax_in            = 5000    ! nmax=500 is a reasonable choice.

        ngeom_recon = 2 * sdim + 3
        num_materials = nmat
        nmax_in = nmax

        Allocate(LS_stencil(-1:1,-1:1,-1:1,nmat))
        Allocate(xtetlist_vof(sdim+1,sdim,nmax))
        Allocate(xtetlist_cen(sdim+1,sdim,nmax)) 
        Allocate(xsten0(-nhalf0:nhalf0,sdim))
        Allocate(mofdata(nmat*ngeom_recon))
        Allocate(multi_centroidA(nmat,sdim))


        LS_stencil      = 0.d0
        xtetlist_vof    = 0.d0
        xtetlist_cen    = 0.d0
        xsten0          = 0.d0
        mofdata         = 1.d0
        multi_centroidA = 0.d0

        xsten0(-3,1) = -1.5d0
        xsten0(-2,1) = -1.0d0
        xsten0(-1,1) = -0.5d0
        xsten0(0,1)  = 0.d0
        xsten0(1,1)  = 0.5d0
        xsten0(2,1)  = 1.0d0
        xsten0(3,1)  = 1.5d0

        xsten0(-3,2) = -1.5d0
        xsten0(-2,2) = -1.0d0
        xsten0(-1,2) = -0.5d0
        xsten0(0,2)  = 0.d0
        xsten0(1,2)  = 0.5d0
        xsten0(2,2)  = 1.0d0
        xsten0(3,2)  = 1.5d0

        xsten0(-3,3) = -1.5d0
        xsten0(-2,3) = -1.0d0
        xsten0(-1,3) = -0.5d0
        xsten0(0,3)  = 0.d0
        xsten0(1,3)  = 0.5d0
        xsten0(2,3)  = 1.0d0
        xsten0(3,3)  = 1.5d0

        Call initmof( &
            order_algorithm_in, &
            nmat,MOFITERMAX_in, &
            MOF_DEBUG_RECON_in, &
            MOF_TURN_OFF_LS_in, &
            nthreads, &
            nmax_in)

    End Subroutine Initialize_SussmanMOF


  End Module SussmanMOF
