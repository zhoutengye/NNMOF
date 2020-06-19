Module HypreStruct
  Use ModGlobal
  Implicit None

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

  Procedure(InterfaceSparse), Pointer :: Sparse_Matrix_Vector => Sparse_Matrix_Vector_Neumann_0
  Procedure(), Pointer :: Hypre_Solver => Hypre_Solver_BicgSTAB
  Procedure(), Pointer :: Hypre_SetSolver => Hypre_SetSolver_BicgSTAB
  Procedure(), Pointer :: Hypre_SetPreConditioner => Hypre_SetPreConditioner_SMG
  Procedure(), Pointer :: Hypre_PreConditioner => Hypre_PreConditioner_SMG

  Integer :: nxyz
  Real(sp) :: iter_tolerance
  Integer  :: iter_max

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

  Subroutine Hypre_Initialize(tol, nmax)

    Implicit None
    Real(sp), Intent(In) :: tol
    Integer , Intent(In) :: nmax

    iter_tolerance = tol
    iter_max = nmax

    Call HYPRE_StructBICGSTABCreate(MPI_COMM_WORLD, hypre%solver, Hypre%ier)
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

      ! Assign each of the 5 stencil entries
      Do entry = 1,7
        offset = offsets(entry,:)
        Call HYPRE_StructStencilSetElement(hypre%stencil, entry-1, offset, Hypre%ier);
      End Do
      hypre%precond_id = 9
    End Block

  End Subroutine Hypre_Initialize

  Subroutine Hypre_Poisson(P, Rhox, Rhoy, Rhoz, Div, flag)
    Implicit None
    real(sp), intent(inout), dimension(0:,0:,0:) :: P
    real(sp), intent(in),    dimension(0:,0:,0:) :: Rhox, Rhoy, Rhoz
    real(sp), intent(in),    dimension(:,:,:) :: Div
    Integer, intent(in),    dimension(:,:,:) :: flag
    Integer  :: kk
    Real(SP) :: Phi(nxyz)
    Integer :: i, j, k

    Call Sparse_Matrix_Vector(P, Div, Rhox, Rhoy, Rhoz, flag)
    Call Hypre_PreConditioner
    Call Hypre_Solver

    call HYPRE_StructVectorGetBoxValues(hypre%x, hypre%ilower, hypre%iupper, Phi, Hypre%ier)
    kk = 1
    Do k = 1, nl(3)
      Do j = 1, nl(2)
        Do i = 1, nl(1)
          P(i,j,k) = Phi(kk)
          kk = kk + 1
        EndDo
      EndDo
    EndDo


  End Subroutine Hypre_Poisson

  Subroutine Sparse_Matrix_Vector_Neumann_0(P, Div, Rhox, Rhoy, Rhoz, flag)
    Implicit None
    real(sp), intent(in), dimension(0:,0:,0:) :: P
    real(sp), intent(in), dimension(0:,0:,0:) :: Rhox, Rhoy, Rhoz
    real(sp), intent(in), dimension(:,:,:) :: Div
    Integer , intent(in), dimension(:,:,:) :: flag
    Integer :: i, j, k
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
        Do j = 1, nl(3)
          Do i = 1, nl(3)

          If ( left  .eq. MPI_PROC_NULL .and. i .eq. 1 ) then
            values(kk+1) = 0.0_SP
          Else
            values(kk+1) = RHOx(i-1,j,k)   / nl(1) / nl(1)
          End If

          If ( right .eq. MPI_PROC_NULL .and. i .eq. nl(1) ) then
            values(kk+2) = 0.0_SP
          Else
            values(kk+2) = Rhox(i,j,k) / nl(1) / nl(1)
          End If

          If ( front .eq. MPI_PROC_NULL .and. j .eq. 1 ) then
            values(kk+3) = 0.0_SP
          Else
            values(kk+3) = Rhoy(i,j-1,k)   / nl(2) / nl(2)
          End If

          If ( back  .eq. MPI_PROC_NULL .and. j .eq. nl(2) ) then
            values(kk+4) = 0.0_SP
          Else
            values(kk+4) = Rhoy(i,j,k) / nl(2) / nl(2)
          End If

          If ( bottom .eq. MPI_PROC_NULL .and. k .eq. 1 ) then
            values(kk+5) = 0.0_SP
          Else
            values(kk+5) = Rhoz(i,j,k-1) / nl(3) / nl(3)
          End If

          If ( top .eq. MPI_PROC_NULL .and. k .eq. nl(3) ) then
            values(kk+6) = 0.0_SP
          Else
            values(kk+6) = Rhoz(i,j,k) / nl(3) / nl(3)
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
          values(kk) = - (values(kk+1) + values(kk+2) + values(kk+3) + values(kk+4) )
          kk = kk + 5
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


  Subroutine Hypre_SetSolver_BicgSTAB
    Implicit None
    Call HYPRE_StructBICGSTABCreate(MPI_COMM_WORLD, hypre%solver, Hypre%ier)
    Call HYPRE_StructBICGSTABSetTol(hypre%solver, iter_tolerance, Hypre%ier)
    ! Call HYPRE_StructBiCGSTABSetPrintLev(hypre%solver, 2, Hypre%ier)
    Call HYPRE_StructBICGSTABSetMaxIter(hypre%solver, iter_max, Hypre%ier)

  End Subroutine Hypre_SetSolver_BicgSTAB

  Subroutine Hypre_Solver_BicgSTAB
    Implicit None

    Call HYPRE_StructBICGSTABSetup(hypre%solver, hypre%A, hypre%b, hypre%x, Hypre%ier)
    Call HYPRE_StructBICGSTABSolve(hypre%solver, hypre%A, hypre%b, hypre%x, Hypre%ier)

  End Subroutine Hypre_Solver_BicgSTAB

  Subroutine Hypre_SetPreConditioner_SMG
    Implicit None

    Call HYPRE_StructSMGCreate(MPI_COMM_WORLD, hypre%precond, Hypre%ier)
    Call HYPRE_StructSMGSetMemoryUse(hypre%precond, 0, Hypre%ier)
    Call HYPRE_StructSMGSetMaxIter(hypre%precond, 1, Hypre%ier)
    Call HYPRE_StructSMGSetTol(hypre%precond, 0.0, Hypre%ier)
    Call HYPRE_StructSMGSetZeroGuess(hypre%precond, Hypre%ier)
    Call HYPRE_StructSMGSetNumPreRelax(hypre%precond, 1, Hypre%ier)
    Call HYPRE_StructSMGSetNumPostRelax(hypre%precond, 1, Hypre%ier)
    Hypre%precond_id = 0

  End Subroutine Hypre_SetPreConditioner_SMG

  Subroutine Hypre_PreConditioner_SMG
    Implicit None

    Call HYPRE_StructBICGSTABSetPrecond(hypre%solver, hypre%precond_id, hypre%precond, Hypre%ier)

  End Subroutine Hypre_PreConditioner_SMG

End Module HypreStruct
