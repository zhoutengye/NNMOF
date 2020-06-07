module module_CUDA
	use module_grid
	use module_BC
	integer :: dimM(3)

	!--------------------------------------C/CU interface--------------------------------------------
	INTERFACE
		SUBROUTINE linearCUDA(A,p,dimM,Ng,maxError,beta,maxit,it,ierr,norm,tres2) BIND(C, name="linearCUDA_")
			USE, INTRINSIC :: ISO_C_BINDING, ONLY: C_INT
			IMPLICIT NONE
			INTEGER (C_INT) ::  dimM(0:2), Ng, maxit, it, ierr, norm
			real*8 :: A(0:dimM(0)-1,0:dimM(1)-1,0:dimM(2)-1,0:7), p(0:dimM(0)-1,0:dimM(1)-1,0:dimM(2)-1)
			real*8 :: maxError, beta, tres2
		END SUBROUTINE linearCUDA
	END INTERFACE

	contains

	subroutine NewSolver_CUDA(A,p,maxError,beta,maxit,it,ierr,norm,tres2)
		use module_grid
		use module_BC
		use module_IO
		use module_freesurface
		implicit none
		include 'mpif.h'
		real(8), dimension(imin:imax,jmin:jmax,kmin:kmax), intent(inout) :: p
		real(8), dimension(is:ie,js:je,ks:ke,8), intent(in) :: A
		real(8), intent(in)  :: beta, maxError
		real(8), intent(out) :: tres2
		integer, intent(in) :: maxit
		integer, intent(out) :: it, ierr
		real(8) :: res2,intvol
		integer :: i,j,k
		logical :: mask(imin:imax,jmin:jmax,kmin:kmax)
		integer, intent(in) :: norm
		integer, save :: itime=0
		dimM = (/imax-imin+1,jmax-jmin+1,kmax-kmin+1/)
	!~ 	allocate(pvec((imax-imin)*(jmax-jmin)*(kmax-kmin)))
		call linearCUDA(A,p,dimM,Ng,maxError,beta,maxit,it,ierr,norm,tres2)
	end subroutine NewSolver_CUDA




end module module_CUDA

