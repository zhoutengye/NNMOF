!=================
! Includes:
!     (1) Derived types and procedure pointer
!     (2) Hypre interface for public
!       (2-1) Hypre Initialize
!       (2-2) Hypre Poisson
!     (3) Construct the matrix
!       (3-1) Construct Sparse matrix and RHS
!     (4) Set Solvers
!       (4-1) Set solver SMG
!       (4-2) Set solver PFMG
!       (4-3) Set solver BicgSTAB
!       (4-4) Set solver GMRES
!     (5) Set Pre-conditioner
!       (5-1) Set preconditioner SMG
!       (5-2) Set preconditioner PFMG
!     (6) Solvers
!       (6-1) SMG
!       (6-2) PFMG
!       (6-3) BicgSTAB
!       (6-4) GMRES
!-----------------
! Author: Zhouteng Ye (yzt9zju@gmail.com)
!-----------------
! Note:
!  1) The following variables are not passed from the subroutine interface,
!  but comes from the ModGlobal
!     mpi variables
!     dt: time step
!     nl: number of total grids, with rank 3
!     dl: grid size, with rank 3
!
!  2)  List of solvers:
!    In input.namelist, set hypre_solver:
!       1: SMG
!       2: PFMG
!       3: BicgSTAB
!       4: GMRES
!    List of pre_conditioner:
!       0: No preconditioner
!       1: SMG
!       2: FPMG
! Working pairs are
!     |------------|-----|------|----------|-------|
!     |            | SMG | FPMG | BicgSTAB | GMRES |
!     |------------|-----|------|----------|-------|
!     |   no pre:  | Yes | Yes  | Yes      | Yes   |
!     |   SMG:     |     |      | Yes      | Yes   |
!     |   FPMG:    |     |      | Yes      | Yes   |
!     |------------|-----|------|----------|-------|
!=================

Module HypreStruct
  Use ModGlobal
  Use ModTools
  Implicit None
  Private
  Public :: Hypre_Initialize, Hypre_Poisson

  ! (1) Derived types 
  Type :: Hypre_para
    Integer :: ier
    Integer*8 :: grid
    Integer*8 :: stencil
    Integer*8 :: A
    Integer*8 :: b
    Integer*8 :: x
    Integer*8 :: solver
    Integer*8 :: precond
    integer*8 :: precond_id
    Integer :: ilower(3)
    Integer :: iupper(3)
    Integer :: n_grid
  End Type Hypre_para

  Type(Hypre_para) :: Hypre

  ! (1) Procedure pointers
  Procedure(InterfaceSparse), Pointer :: Sparse_Matrix_Vector => Sparse_Matrix_Vector_Neumann_0
  Procedure(), Pointer :: Hypre_Solver => Hypre_Solver_GMRES
  Procedure(), Pointer :: Hypre_SetSolver => Hypre_SetSolver_GMRES
  Procedure(), Pointer :: Hypre_SetPreConditioner => Hypre_SetPreConditioner_SMG

  Integer :: nxyz
  Real(sp) :: iter_tolerance
  Integer  :: iter_max
  Integer :: niter
  Real(sp) :: residual_norm

  Interface
    Subroutine InterfaceSparse(P, Div, Rhox, Rhoy, Rhoz, flag)
      import
      Implicit None
      Real(sp), intent(in), dimension(0:,0:,0:) :: P
      Real(sp), intent(in), dimension(0:,0:,0:) :: Rhox, Rhoy, Rhoz
      Real(sp), intent(in), dimension(:,:,:) :: Div
      Integer, intent(in), dimension(:,:,:) :: flag
    End Subroutine InterfaceSparse
  End Interface

Contains

  !=======================
  ! (2-1) Initialize hypre solver
  !   Setup solver, preconditioner, maximum iteration step
  !-----------------------
  ! Inputs:
  !    tol: Relative tolerance for iteration
  !    nnmax: Maximum iteration steps
  !    set_hypre_solver: Assin the type of hypre solver
  !    set_hypre_preconditioner: Assin the type of hypre preconditioner
  !========================
  Subroutine Hypre_Initialize(tol, nmax, set_hypre_solver, set_Hypre_PreConditioner)

    Implicit None
    Real(sp), Intent(In) :: tol
    Integer , Intent(In) :: nmax
    Integer , Intent(In) :: set_hypre_solver, set_hypre_preConditioner
    External HYPRE_StructGridCreate, Hypre_structGridSetExtents, HYPRE_StructGridAssemble
    External HYPRE_StructStencilCreate, HYPRE_StructStencilSetElement

    iter_tolerance = tol
    iter_max = nmax

    ! Set pre conditioner
    hypre%precond_id = 9
    if (set_hypre_preconditioner .eq. 0) then
      hypre%precond_id = 9
    else if (set_Hypre_PreConditioner .eq. 1) then
      hypre%precond_id = 0
      Hypre_SetPreConditioner => Hypre_SetPreConditioner_SMG
    else if (set_Hypre_PreConditioner .eq. 2) then
      hypre%precond_id = 1
      Hypre_SetPreConditioner => Hypre_SetPreConditioner_PFMG
    else
      If ( myid .eq. 0 ) then
        print *, "======Fatal Error=============================="
        print *, "Incorrect preconditioner, should be:"
        print *, "0: No preconditioner"
        print *, "1: SMG preconditioner"
        print *, "2: PFMG preconditioner"
        print *, "==============================================="
      End If
      Call MPI_Finalize(ierr)
      stop
    endif

    ! Set iterator
    if (Set_Hypre_Solver .eq. 1) then
      Hypre_SetSolver => Hypre_SetSolver_SMG
      Hypre_Solver => Hypre_Solver_SMG
    else if (Set_Hypre_Solver .eq. 2) then
      Hypre_SetSolver => Hypre_SetSolver_PFMG
      Hypre_Solver => Hypre_Solver_PFMG
    else if (Set_Hypre_Solver .eq. 3) then
      Hypre_SetSolver => Hypre_SetSolver_BicgSTAB
      Hypre_Solver => Hypre_Solver_BicgSTAB
    else if (Set_Hypre_Solver .eq. 4) then
      Hypre_SetSolver => Hypre_SetSolver_GMRES
      Hypre_Solver => Hypre_Solver_GMRES
    else
      If ( myid .eq. 0 ) then
        print *, "======Fatal Error=============================="
        print *, "Incorrect hypre solver, should be:"
        print *, "1: SMG"
        print *, "2: PFMG"
        print *, "3: BicgSTAB"
        print *, "4: GMRES"
        print *, "==============================================="
      End If
      Call MPI_Finalize(ierr)
      stop
    endif

    Call Hypre_SetPreConditioner
    Call Hypre_SetSolver

    ! Set up the grid index
    Block
      Hypre%ilower(1) = coord(1) * nl(1) + 1
      Hypre%iupper(1) = ( coord(1) + 1 ) * nl(1)
      Hypre%ilower(2) = coord(2) * nl(2) + 1
      Hypre%iupper(2) = ( coord(2) + 1 ) * nl(2)
      Hypre%ilower(3) = 1
      Hypre%iupper(3) = nl(3)
      nxyz = nl(1) * nl(2) * nl(3) 
    End Block
    
    !-1. Set up a grid. 
    
    ! Create an empty 3D grid object */
    Block
      Call HYPRE_StructGridCreate(MPI_COMM_WORLD, 3, hypre%grid, Hypre%ier)
      Call HYPRE_StructGridSetExtents(hypre%grid, hypre%ilower, hypre%iupper, Hypre%ier)

      Call HYPRE_StructGridAssemble(hypre%grid, Hypre%ier)
    End Block

    !--2. Define the discretization stencil------------------
    ! Create an empty 2D, 5-pt stencil object
    Block
      Integer :: entry
      Integer :: offsets(7,3)
      Integer :: offset(3)
      Call HYPRE_StructStencilCreate(3, 7, hypre%stencil, Hypre%ier)

      ! Define the geometry of the stencil. Each represents a
      ! relative offset (in the index space)
      offsets(1,:) = [ 0, 0, 0]
      offsets(2,:) = [-1, 0, 0]
      offsets(3,:) = [ 1, 0, 0]
      offsets(4,:) = [ 0,-1, 0]
      offsets(5,:) = [ 0, 1, 0]
      offsets(6,:) = [ 0, 0,-1]
      offsets(7,:) = [ 0, 0, 1]

      ! Assign each of the 7 stencil entries
      Do entry = 1,7
        offset = offsets(entry,:)
        Call HYPRE_StructStencilSetElement(hypre%stencil, entry-1, offset, Hypre%ier);
      End Do
    End Block

  End Subroutine Hypre_Initialize

  !=======================
  ! (2-2) Solve the Poisson equatiosn
  ! Key steps:
  !   (1) Construct sparse matrix and RHS vector
  !   (2) Solve Ap = Div
  !   (3) Get the value of p
  !-----------------------
  ! Inputs:
  !    P: Pressure
  !    Rhox: x face-centered conponent of 1 / Rho
  !    Rhoy: y face-centered conponent of 1 / Rho
  !    Rhoz: z face-centered conponent of 1 / Rho
  !    Div: Divergence
  !    flag: used to indicate fluid region
  ! Outputs:
  !    P: Pressure
  !    n_iter: Number of iteration
  !    res: Relative residual norm
  !========================
  Subroutine Hypre_Poisson(P, Rhox, Rhoy, Rhoz, Div, flag, n_iter, res)
    Implicit None
    real(sp), intent(inout), dimension(0:,0:,0:) :: P
    real(sp), intent(in),    dimension(0:,0:,0:) :: Rhox, Rhoy, Rhoz
    real(sp), intent(in),    dimension(:,:,:) :: Div
    Integer, intent(in),    dimension(:,:,:) :: flag
    Integer, intent(out) :: n_iter
    Real(sp), intent(out) :: res
    Integer  :: kk
    Real(SP) :: PP(nxyz)
    Integer :: i, j, k
    External HYPRE_StructVectorGetBoxValues

    Call Sparse_Matrix_Vector(P, Div, Rhox, Rhoy, Rhoz, flag)
    Call Hypre_Solver

    call HYPRE_StructVectorGetBoxValues(hypre%x, hypre%ilower, hypre%iupper, PP, Hypre%ier)
    kk = 1
    Do k = 1, nl(3)
      Do j = 1, nl(2)
        Do i = 1, nl(1)
          P(i,j,k) = PP(kk)
          kk = kk + 1
        EndDo
      EndDo
    EndDo

    n_iter = niter
    res = residual_norm

  End Subroutine Hypre_Poisson

  !=======================
  ! (3-1) Construct Sparse matrix and RHS
  !-----------------------
  ! Inputs:
  !    P: Pressure
  !    Rhox: x face-centered conponent of 1 / Rho
  !    Rhoy: y face-centered conponent of 1 / Rho
  !    Rhoz: z face-centered conponent of 1 / Rho
  !    Div: Divergence
  !    flag: used to indicate fluid region
  ! Outputs:
  !    P: Pressure
  !    n_iter: Number of iteration
  !    res: Relative residual norm
  !========================
  Subroutine Sparse_Matrix_Vector_Neumann_0(P, Div, Rhox, Rhoy, Rhoz, flag)
    Implicit None
    real(sp), intent(in), dimension(0:,0:,0:) :: P
    real(sp), intent(in), dimension(0:,0:,0:) :: Rhox, Rhoy, Rhoz
    real(sp), intent(in), dimension(:,:,:) :: Div
    Integer , intent(in), dimension(:,:,:) :: flag
    Integer :: i, j, k
    External :: HYPRE_StructMatrixCreate, HYPRE_StructMatrixInitialize
    External :: HYPRE_StructMatrixSetBoxValues, HYPRE_StructMatrixAssemble
    External :: HYPRE_StructVectorCreate, HYPRE_StructVectorInitialize
    External :: HYPRE_StructVectorSetBoxValues

    Block
      Integer :: stencil_indices(7)
      Integer :: nentries
      Integer :: nvalues
      Real(SP), Allocatable :: values(:)
      Integer :: kk
      Integer :: i,j,k
      ! Create an empty matrix object
      Call HYPRE_StructMatrixCreate(MPI_COMM_WORLD, hypre%grid, hypre%stencil, hypre%A, Hypre%ier);

      ! Indicate that the matrix coefficients are ready to be set
      Call HYPRE_StructMatrixInitialize(hypre%A, Hypre%ier)

      ! Set the matrix coefficients.  Each processor assigns coefficients
      ! for the boxes in the grid that it owns. Note that the coefficients
      ! associated with each stencil entry may vary from grid point to grid
      ! point if desired.  Here, we first set the same stencil entries for
      ! each grid point.  Then we make modifications to grid points near
      ! the boundary.

      ! labels for the stencil entries these correspond to the offsets
      ! defined above
      stencil_indices = [0,1,2,3,4,5,6]
      nentries = 7
      nvalues = nxyz * nentries
      Allocate(values(nvalues))

      kk = 1
      Do k = 1, nl(3)
        Do j = 1, nl(2)
          Do i = 1, nl(1)

          If ( left  .eq. MPI_PROC_NULL .and. i .eq. 1 ) then
            values(kk+1) = 0.0_SP
          Else
            values(kk+1) = RHOx(i-1,j,k)   / dl(1) / dl(1)
          End If

          If ( right .eq. MPI_PROC_NULL .and. i .eq. nl(1) ) then
            values(kk+2) = 0.0_SP
          Else
            values(kk+2) = Rhox(i,j,k) / dl(1) / dl(1)
          End If

          If ( front .eq. MPI_PROC_NULL .and. j .eq. 1 ) then
            values(kk+3) = 0.0_SP
          Else
            values(kk+3) = Rhoy(i,j-1,k)   / dl(2) / dl(2)
          End If

          If ( back  .eq. MPI_PROC_NULL .and. j .eq. nl(2) ) then
            values(kk+4) = 0.0_SP
          Else
            values(kk+4) = Rhoy(i,j,k) / dl(2) / dl(2)
          End If

          If ( bottom .eq. MPI_PROC_NULL .and. k .eq. 1 ) then
            values(kk+5) = 0.0_SP
          Else
            values(kk+5) = Rhoz(i,j,k-1) / dl(3) / dl(3)
          End If

          If ( top .eq. MPI_PROC_NULL .and. k .eq. nl(3) ) then
            values(kk+6) = 0.0_SP
          Else
            values(kk+6) = Rhoz(i,j,k) / dl(3) / dl(3)
          End If

          ! If (flag_c(i,j) .le. 0) Then
            ! Values(kk+1:kk+4) = 0.0_SP
          ! End If

          kk = kk + 7

        EndDo
      EndDo
    EndDo

    kk = 1
    Do k = 1, nl(3)
      Do j = 1, nl(2)
        Do i = 1, nl(1)
          values(kk) = - ( values(kk+1) + values(kk+2) + values(kk+3) + values(kk+4) + values(kk+5) + values(kk+6) )
          kk = kk + 7
        EndDo
      EndDo
    EndDo


      Call HYPRE_StructMatrixSetBoxValues(hypre%A, hypre%ilower, hypre%iupper, nentries, &
          stencil_indices, values, Hypre%ier)
      ! This is a collective call finalizing the matrix assembly.
      ! The matrix is now ``ready to be used'' 
      Call HYPRE_StructMatrixAssemble(hypre%A, Hypre%ier)
      
    End Block

    !--4. Set up Struct Vectors for b and x.
    ! Each processor sets the vectors corresponding to its boxes.
    Block
      Integer :: nvalues
      Real(SP), Allocatable :: values(:)
      Integer :: kk

      nvalues = nxyz
      Allocate(values(nvalues))

      Call HYPRE_StructVectorCreate(MPI_COMM_WORLD, hypre%grid, hypre%b, Hypre%ier)
      Call HYPRE_StructVectorCreate(MPI_COMM_WORLD, hypre%grid, hypre%x, Hypre%ier)

      Call HYPRE_StructVectorInitialize(hypre%b, Hypre%ier)
      Call HYPRE_StructVectorInitialize(hypre%x, Hypre%ier)

      kk = 1
      Do k = 1, nl(3)
        Do j = 1, nl(2)
          Do i = 1, nl(1)
            values(kk) = Div(i,j,k) / dt
            If (flag(i,j,k) .eq. 0) values(kk) = 0.0_SP
            kk = kk + 1
          EndDo
        EndDo
      EndDo

      Call HYPRE_StructVectorSetBoxValues(hypre%b, hypre%ilower, hypre%iupper, values, Hypre%ier)

      kk = 1
      Do k = 1, nl(3)
        Do j = 1, nl(2)
          Do i = 1, nl(1)
            values(kk) = P(i,j,k)
            If (flag(i,j,k) .lt. 0) values(kk) = 0.0_SP
            kk = kk + 1
          EndDo
        EndDo
      EndDo

      Call HYPRE_StructVectorSetBoxValues(hypre%x, hypre%ilower, hypre%iupper, values, Hypre%ier)
    End Block


  end Subroutine Sparse_Matrix_Vector_Neumann_0

  !=============================
  ! (4-1) Set up SMG Solver
  !=============================
  Subroutine Hypre_SetSolver_SMG
    Implicit None
    External :: HYPRE_StructSMGCreate, HYPRE_StructSMGSetTol, HYPRE_StructSMGSetMaxIter
    Call HYPRE_StructSMGCreate(MPI_COMM_WORLD, hypre%solver, Hypre%ier)
    Call HYPRE_StructSMGSetTol(hypre%solver, iter_tolerance, Hypre%ier)
    Call HYPRE_StructSMGSetMaxIter(hypre%solver, iter_max, Hypre%ier)
  End Subroutine Hypre_SetSolver_SMG
  !=============================
  ! (4-2) Set up PFMG Solver
  !=============================
  Subroutine Hypre_SetSolver_PFMG
    Implicit None
    External :: HYPRE_StructPFMGCreate, HYPRE_StructPFMGSetTol, HYPRE_StructPFMGSetMaxIter
    Call HYPRE_StructPFMGCreate(MPI_COMM_WORLD, hypre%solver, Hypre%ier)
    Call HYPRE_StructPFMGSetTol(hypre%solver, iter_tolerance, Hypre%ier)
    Call HYPRE_StructPFMGSetMaxIter(hypre%solver, iter_max, Hypre%ier)
  End Subroutine Hypre_SetSolver_PFMG
  !=============================
  ! (4-3) Set up BiCGSTAB Solver
  !=============================
  Subroutine Hypre_SetSolver_BicgSTAB
    Implicit None
    External :: HYPRE_StructBICGSTABCreate, HYPRE_StructBICGSTABSetTol, HYPRE_StructBICGSTABSetMaxIter
    Call HYPRE_StructBICGSTABCreate(MPI_COMM_WORLD, hypre%solver, Hypre%ier)
    Call HYPRE_StructBICGSTABSetTol(hypre%solver, iter_tolerance, Hypre%ier)
    Call HYPRE_StructBICGSTABSetMaxIter(hypre%solver, iter_max, Hypre%ier)
  End Subroutine Hypre_SetSolver_BicgSTAB
  !=============================
  ! (4-4) Set up GMRES Solver
  !=============================
  Subroutine Hypre_SetSolver_GMRES
    Implicit None
    External :: HYPRE_StructGMRESCreate, HYPRE_StructGMRESSetTol, HYPRE_StructGMRESSetMaxIter
    Call HYPRE_StructGMRESCreate(MPI_COMM_WORLD, hypre%solver, Hypre%ier)
    Call HYPRE_StructGMRESSetTol(hypre%solver, iter_tolerance, Hypre%ier)
    Call HYPRE_StructGMRESSetMaxIter(hypre%solver, iter_max, Hypre%ier)
  End Subroutine Hypre_SetSolver_GMRES

  !=============================
  ! (6-1) SMG Solver
  !=============================
  Subroutine Hypre_Solver_SMG
    Implicit None
    External :: HYPRE_StructSMGSetup, HYPRE_StructSMGSolve
    External :: HYPRE_StructSMGGetNumIterations, HYPRE_Structsmggetfinalrelative
    Call HYPRE_StructSMGSetup(hypre%solver, hypre%A, hypre%b, hypre%x, Hypre%ier)
    Call HYPRE_StructSMGSolve(hypre%solver, hypre%A, hypre%b, hypre%x, Hypre%ier)
    Call HYPRE_StructSMGGetNumIterations(hypre%solver, niter, Hypre%ier)
    Call HYPRE_Structsmggetfinalrelative(hypre%solver, residual_norm, hypre%ier)
  End Subroutine Hypre_Solver_SMG
  !=============================
  ! (6-2) PFMG Solver
  !=============================
  Subroutine Hypre_Solver_PFMG
    Implicit None
    External :: HYPRE_StructPFMGSetup, HYPRE_StructPFMGSolve
    External :: HYPRE_StructPFMGGetNumIteration, HYPRE_Structpfmggetfinalrelativ
    Call HYPRE_StructPFMGSetup(hypre%solver, hypre%A, hypre%b, hypre%x, Hypre%ier)
    Call HYPRE_StructPFMGSolve(hypre%solver, hypre%A, hypre%b, hypre%x, Hypre%ier)
    Call HYPRE_StructPFMGGetNumIteration(hypre%solver, niter, Hypre%ier)
    Call hypre_structpfmggetfinalrelativ(hypre%solver, residual_norm, hypre%ier)
  End Subroutine Hypre_Solver_PFMG
  !=============================
  ! (6-3) BicgSTAB Solver
  !=============================
  Subroutine Hypre_Solver_BicgSTAB
    Implicit None
    External :: HYPRE_StructBICGSTABSetPrecond
    External :: HYPRE_StructbicgstabSetup, HYPRE_StructbicgstabSolve
    External :: HYPRE_StructbicgstabGetNumItera, HYPRE_Structbicgstabgetfinalrel
    Call HYPRE_StructBICGSTABSetPrecond(hypre%solver, hypre%precond_id, hypre%precond, Hypre%ier)
    Call HYPRE_StructBICGSTABSetup(hypre%solver, hypre%A, hypre%b, hypre%x, Hypre%ier)
    Call HYPRE_StructBICGSTABSolve(hypre%solver, hypre%A, hypre%b, hypre%x, Hypre%ier)
    Call hypre_structbicgstabgetnumitera(hypre%solver, niter, Hypre%ier)
    Call hypre_structbicgstabgetfinalrel(hypre%solver, residual_norm, hypre%ier)
  End Subroutine Hypre_Solver_BicgSTAB
  !=============================
  ! (6-4) GMRES Solver
  !=============================
  Subroutine Hypre_Solver_GMRES
    Implicit None
    External :: HYPRE_StructGMRESSetPrecond
    External :: HYPRE_StructgmresSetup, HYPRE_StructgmresSolve
    External :: HYPRE_StructgmresGetNumIteratio, HYPRE_Structgmresgetfinalrelati
    Call HYPRE_StructgmresSetPrecond(hypre%solver, hypre%precond_id, hypre%precond, Hypre%ier)
    Call HYPRE_StructGMRESSetup(hypre%solver, hypre%A, hypre%b, hypre%x, Hypre%ier)
    Call HYPRE_StructGMRESSolve(hypre%solver, hypre%A, hypre%b, hypre%x, Hypre%ier)
    Call HYPRE_StructGMRESGetNumIteratio(hypre%solver, niter, Hypre%ier)
    Call hypre_structGMRESgetfinalrelati(hypre%solver, residual_norm, hypre%ier)
  End Subroutine Hypre_Solver_GMRES

  Subroutine Hypre_SetPreConditioner_SMG
    Implicit None
    External :: HYPRE_StructSMGCreate
    External :: HYPRE_StructSMGSetMemoryUse, HYPRE_StructSMGSetMaxIter
    External :: HYPRE_StructSMGSetTol, HYPRE_StructSMGSetZeroGuess
    External :: HYPRE_StructSMGSetNumPreRelax, HYPRE_StructSMGSetNumPostRelax
    Call HYPRE_StructSMGCreate(MPI_COMM_WORLD, hypre%precond, Hypre%ier)
    Call HYPRE_StructSMGSetMemoryUse(hypre%precond, 0, Hypre%ier)
    Call HYPRE_StructSMGSetMaxIter(hypre%precond, 1, Hypre%ier)
    Call HYPRE_StructSMGSetTol(hypre%precond, 0.0, Hypre%ier)
    Call HYPRE_StructSMGSetZeroGuess(hypre%precond, Hypre%ier)
    Call HYPRE_StructSMGSetNumPreRelax(hypre%precond, 1, Hypre%ier)
    Call HYPRE_StructSMGSetNumPostRelax(hypre%precond, 1, Hypre%ier)
  End Subroutine Hypre_SetPreConditioner_SMG

  Subroutine Hypre_SetPreConditioner_PFMG
    Implicit None
    External :: HYPRE_StructPFMGCreate
    External :: HYPRE_StructPFMGSetMaxIter
    External :: HYPRE_StructPFMGSetTol, HYPRE_StructPFMGSetZeroGuess
    External :: HYPRE_StructPFMGSetNumPreRelax, HYPRE_StructPFMGSetNumPostRelax
    Call HYPRE_StructPFMGCreate(MPI_COMM_WORLD, hypre%precond, Hypre%ier)
    Call HYPRE_StructPFMGSetMaxIter(hypre%precond, 1, Hypre%ier)
    Call HYPRE_StructPFMGSetTol(hypre%precond, 0.0, Hypre%ier)
    Call HYPRE_StructPFMGSetZeroGuess(hypre%precond, Hypre%ier)
    Call HYPRE_StructPFMGSetNumPreRelax(hypre%precond, 1, Hypre%ier)
    Call HYPRE_StructPFMGSetNumPostRelax(hypre%precond, 1, Hypre%ier)
  End Subroutine Hypre_SetPreConditioner_PFMG

End Module HypreStruct
