subroutine apply_BC_CUDA(pvec)
	use module_grid
	use module_BC
	implicit none
	include 'mpif.h'
	real(8) :: p(imin:imax,jmin:jmax,kmin:kmax)
	real(8), intent(inout) :: pvec((imax-imin+1)*(jmax-jmin+1)*(kmax-kmin+1)) ! INEFFICIENT???
	integer :: req(12),sta(MPI_STATUS_SIZE,12)
	integer :: L,ierr,it
	logical :: exist

	p = reshape(pvec,(/imax-imin+1,jmax-jmin+1,kmax-kmin+1/))
	
	call ghost_x(p,1,req( 1: 4))
	call ghost_y(p,1,req( 5: 8))
	call ghost_z(p,1,req( 9:12))
	call MPI_WAITALL(12,req,sta,ierr)


	pvec = reshape(p,(/(imax-imin+1)*(jmax-jmin+1)*(kmax-kmin+1)/))
end subroutine apply_BC_CUDA

subroutine catch_divergence_CUDA(res2,ierr,it)
		use module_grid
		use module_BC
		use module_IO
		implicit none

        real(8), intent(in) :: res2
        integer, intent(out) :: ierr
        integer, intent(in) :: it
        logical :: diverged=.false.
        logical :: extended=.true.
        integer :: l, i, j, k
		if(diverged) then	
			OPEN(UNIT=88,FILE=TRIM(out_path)//'/message-rank-'//TRIM(int2text(rank,padding))//'.txt')
			write(88,*) "ijk rank",i,j,k,rank
			write(88,*) 'A or p is NaN after',it,'iterations at rank ',rank
			close(88)
			if(rank<=30) print*,'A or p is NaN after',it,'iterations at rank ',rank
			call pariserror("A or p is NaN")
		endif
        if ((res2*npx*npy*npz)>1.d16 ) then
            if(rank<=30) print*,'Pressure solver diverged after',it,'iterations at rank ',rank
            call pariserror("newsolver CUDA error")
        else if (res2 .ne. res2) then
            if(rank<=30) print*, 'it:',it,'Pressure residual value is invalid at rank', rank
            call pariserror("newsolver CUDA error ")
        else
            ierr=0
        endif
end subroutine catch_divergence_CUDA

subroutine collect_res2_CUDA(res2,tres2,ierr,it)
		use module_grid
		use module_BC
		use module_IO
		implicit none

		include 'mpif.h'
		real(8), intent(in) :: res2
		real(8), intent(inout) :: tres2
		integer, intent(inout) :: ierr
		integer, intent (in) :: it
		310         format(I6,'  ',(e14.5))
			CHARACTER (LEN=11) :: filename
	  logical :: exist

		
		call MPI_ALLREDUCE(res2, tres2, 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_Comm_Cart, ierr)
		if(rank==0.and.mod(it,10) == 0.and.recordconvergence) write(89,310) it, tres2
end subroutine collect_res2_CUDA

