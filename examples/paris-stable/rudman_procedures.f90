
!---------------------------------------------------sub-grid to coarse grid vof computation---------------------
subroutine compute_vof_fine2coarse(c,c_flag,cf)
   use module_rudman_init
   use module_bc
   implicit none
   include 'mpif.h'
   integer :: i,j,k,isf,jsf,ksf,i1,j1,k1,iini,iend,jini,jend,kini,kend
   real(8),dimension(idx_c(1,1):idx_c(1,2),idx_c(2,1):idx_c(2,2),idx_c(3,1):idx_c(3,2)), intent(out) :: c 
   integer,dimension(idx_c(1,1):idx_c(1,2),idx_c(2,1):idx_c(2,2),idx_c(3,1):idx_c(3,2)), intent(inout) :: c_flag 
   real(8),dimension(idx_f(1,1):idx_f(1,2),idx_f(2,1):idx_f(2,2),idx_f(3,1):idx_f(3,2)), intent(in) :: cf 

   iini=idx_c(1,1) + Ng ; iend=idx_c(1,2) - Ng
   jini=idx_c(2,1) + Ng ; jend=idx_c(2,2) - Ng
   kini=idx_c(3,1) + Ng ; kend=idx_c(3,2) - Ng 

   isf = idx_f(1,1)+Ng ; jsf = idx_f(2,1)+Ng ; ksf = idx_f(3,1)+Ng ! Subgrid Indices

   c = 0.d0 ; c_flag = 0
   
   do i=iini,iend ; do j=jini,jend ; do k=kini,kend   

      do i1=0,1 ;  do j1=0,1 ;  do k1=0,1
        c(i,j,k) = c(i,j,k) + &    
                     cf( 2*(i-iini)+isf+i1 , 2*(j-jini)+jsf+j1 , 2*(k-kini)+ksf+k1 )*0.125d0 
      enddo ; enddo ; enddo                         
 
   enddo ; enddo ; enddo

   call rudman_vof_flags_and_clip(c,c_flag,idx_c) ! assign flags
   call vof_bc_rudman(c,c_flag,idx_c,coarse)  ! use flags to assign stag vof bc's
 
   call do_all_ghost(c)
   call do_all_ighost(c_flag)

end subroutine compute_vof_fine2coarse


!-----------------------------------------------------inter-proc communication procedures for sub-grid-----------------
subroutine ghost_x_fine(F,ngh,req)
  use module_rudman_init , only: idx_f
  use module_grid
  implicit none
  include 'mpif.h'
  integer, intent(in) :: ngh ! number of ghost cell layers to fill
  integer, intent(out) :: req(4)
  real(8), dimension(idx_f(1,1):idx_f(1,2),idx_f(2,1):idx_f(2,2),idx_f(3,1):idx_f(3,2)), intent(inout) :: F
  integer :: ierr !,sta(MPI_STATUS_SIZE,4)
  integer, save :: srcL_f, srcR_f, destL_f, destR_f, face_f(2)
  logical, save :: first_time_f = .true.
  integer :: iini,iend,jini,jend,kini,kend

  iini=idx_f(1,1) + ngh ; iend=idx_f(1,2) - ngh 
  jini=idx_f(2,1) + ngh ; jend=idx_f(2,2) - ngh
  kini=idx_f(3,1) + ngh ; kend=idx_f(3,2) - ngh 


  if(first_time_f) then
   first_time_f = .false.
   call init_ghost_fine(srcL_f, srcR_f, destL_f, destR_f, face_f, 1)
  endif

  call MPI_IRECV(F(iini-ngh ,   jini-ngh , kini-ngh),1,face_f(ngh),srcR_f ,0,MPI_COMM_Cart,req(1),ierr)
  call MPI_ISEND(F(iend-ngh+1 , jini-ngh , kini-ngh),1,face_f(ngh),destR_f,0,MPI_COMM_Cart,req(2),ierr)
  call MPI_IRECV(F(iend+1 ,     jini-ngh , kini-ngh),1,face_f(ngh),srcL_f ,0,MPI_COMM_Cart,req(3),ierr)
  call MPI_ISEND(F(iini ,       jini-ngh , kini-ngh),1,face_f(ngh),destL_f,0,MPI_COMM_Cart,req(4),ierr)
end subroutine ghost_x_fine


subroutine ghost_y_fine(F,ngh,req)
  use module_rudman_init , only: idx_f
  use module_grid
  implicit none
  include 'mpif.h'
  integer, intent(in) :: ngh
  integer, intent(out) :: req(4)
  real(8), dimension(idx_f(1,1):idx_f(1,2),idx_f(2,1):idx_f(2,2),idx_f(3,1):idx_f(3,2)), intent(inout) :: F
  integer :: ierr !,sta(MPI_STATUS_SIZE,4)
  integer, save :: srcL_f, srcR_f, destL_f, destR_f, face_f(2)
  logical, save :: first_time_f = .true.
  integer :: iini,iend,jini,jend,kini,kend

  iini=idx_f(1,1) + ngh ; iend=idx_f(1,2) - ngh 
  jini=idx_f(2,1) + ngh ; jend=idx_f(2,2) - ngh
  kini=idx_f(3,1) + ngh ; kend=idx_f(3,2) - ngh 

  if(first_time_f) then
   first_time_f = .false.
   call init_ghost_fine(srcL_f, srcR_f, destL_f, destR_f, face_f, 2)
  endif

  call MPI_IRECV(F(iini-ngh , jini-ngh   , kini-ngh),1,face_f(ngh),srcR_f ,0,MPI_COMM_Cart,req(1),ierr)
  call MPI_ISEND(F(iini-ngh , jend-ngh+1 , kini-ngh),1,face_f(ngh),destR_f,0,MPI_COMM_Cart,req(2),ierr)
  call MPI_IRECV(F(iini-ngh , jend+1     , kini-ngh),1,face_f(ngh),srcL_f ,0,MPI_COMM_Cart,req(3),ierr)
  call MPI_ISEND(F(iini-ngh , jini       , kini-ngh),1,face_f(ngh),destL_f,0,MPI_COMM_Cart,req(4),ierr)
end subroutine ghost_y_fine


subroutine ghost_z_fine(F,ngh,req)
  use module_rudman_init , only: idx_f
  use module_grid
  implicit none
  include 'mpif.h'
  integer, intent(in) :: ngh
  integer, intent(out) :: req(4)
  real(8), dimension(idx_f(1,1):idx_f(1,2),idx_f(2,1):idx_f(2,2),idx_f(3,1):idx_f(3,2)), intent(inout) :: F
  integer :: ierr !,sta(MPI_STATUS_SIZE,4)
  integer, save :: srcL_f, srcR_f, destL_f, destR_f, face_f(2)
  logical, save :: first_time_f = .true.
  integer :: iini,iend,jini,jend,kini,kend

  iini=idx_f(1,1) + ngh ; iend=idx_f(1,2) - ngh 
  jini=idx_f(2,1) + ngh ; jend=idx_f(2,2) - ngh
  kini=idx_f(3,1) + ngh ; kend=idx_f(3,2) - ngh 

  if(first_time_f) then
   first_time_f = .false.
   call init_ghost_fine(srcL_f, srcR_f, destL_f, destR_f, face_f, 3)
  endif

  call MPI_IRECV(F(iini-ngh , jini-ngh , kini-ngh  ),1,face_f(ngh),srcR_f ,0,MPI_COMM_Cart,req(1),ierr)
  call MPI_ISEND(F(iini-ngh , jini-ngh , kend-ngh+1),1,face_f(ngh),destR_f,0,MPI_COMM_Cart,req(2),ierr)
  call MPI_IRECV(F(iini-ngh , jini-ngh , kend+1    ),1,face_f(ngh),srcL_f ,0,MPI_COMM_Cart,req(3),ierr)
  call MPI_ISEND(F(iini-ngh , jini-ngh , kini      ),1,face_f(ngh),destL_f,0,MPI_COMM_Cart,req(4),ierr)
end subroutine ghost_z_fine


subroutine ighost_x_fine(F,ngh,req)
  use module_rudman_init , only: idx_f
  use module_grid
  implicit none
  include 'mpif.h'
  integer, intent(in) :: ngh ! number of ghost cell layers to fill
  integer, intent(out) :: req(4)
  integer, dimension(idx_f(1,1):idx_f(1,2),idx_f(2,1):idx_f(2,2),idx_f(3,1):idx_f(3,2)), intent(inout) :: F
  integer :: jlen, klen, ierr !,sta(MPI_STATUS_SIZE,4)
  integer, save :: srcL_f, srcR_f, destL_f, destR_f, face_f(2)
  logical, save :: first_time_f = .true. 
  integer :: iini,iend,jini,jend,kini,kend

  iini=idx_f(1,1) + ngh ; iend=idx_f(1,2) - ngh 
  jini=idx_f(2,1) + ngh ; jend=idx_f(2,2) - ngh
  kini=idx_f(3,1) + ngh ; kend=idx_f(3,2) - ngh 
 
  if(first_time_f) then
   first_time_f = .false.
  jlen=idx_f(2,2)-idx_f(2,1)+1; klen=idx_f(3,2)-idx_f(3,1)+1; !ilen=ngh
  call para_type_block3a(iini-ngh, iend+ngh, jini-ngh, jend+ngh, 1, jlen, klen, MPI_INTEGER, face_f(1))
  call para_type_block3a(iini-ngh, iend+ngh, jini-ngh, jend+ngh, 2, jlen, klen, MPI_INTEGER, face_f(2))
  call MPI_CART_SHIFT(MPI_COMM_CART, 0, 1, srcR_f, destR_f, ierr)
  call MPI_CART_SHIFT(MPI_COMM_CART, 0,-1, srcL_f, destL_f, ierr)
  endif

  call MPI_IRECV(F(iini-ngh   , jini-ngh , kini-ngh),1,face_f(ngh),srcR_f ,0,MPI_COMM_Cart,req(1),ierr)
  call MPI_ISEND(F(iend-ngh+1 , jini-ngh , kini-ngh),1,face_f(ngh),destR_f,0,MPI_COMM_Cart,req(2),ierr)
  call MPI_IRECV(F(iend+1     , jini-ngh , kini-ngh),1,face_f(ngh),srcL_f ,0,MPI_COMM_Cart,req(3),ierr)
  call MPI_ISEND(F(iini       , jini-ngh , kini-ngh),1,face_f(ngh),destL_f,0,MPI_COMM_Cart,req(4),ierr)
end subroutine ighost_x_fine


subroutine ighost_y_fine(F,ngh,req)
  use module_rudman_init , only: idx_f
  use module_grid
  implicit none
  include 'mpif.h'
  integer, intent(in) :: ngh
  integer, intent(out) :: req(4)
  integer, dimension(idx_f(1,1):idx_f(1,2),idx_f(2,1):idx_f(2,2),idx_f(3,1):idx_f(3,2)), intent(inout) :: F
  integer :: ilen, klen, ierr !,sta(MPI_STATUS_SIZE,4)
  integer, save :: srcL_f, srcR_f, destL_f, destR_f, face_f(2)
  logical, save :: first_time_f = .true. 
  integer :: iini,iend,jini,jend,kini,kend

  iini=idx_f(1,1) + ngh ; iend=idx_f(1,2) - ngh 
  jini=idx_f(2,1) + ngh ; jend=idx_f(2,2) - ngh
  kini=idx_f(3,1) + ngh ; kend=idx_f(3,2) - ngh 
 
  if(first_time_f) then
   first_time_f = .false.
  klen=idx_f(3,2)-idx_f(3,1)+1; ilen=idx_f(1,2)-idx_f(1,1)+1; !jlen=ngh
  call para_type_block3a(iini-ngh, iend+ngh, jini-ngh, jend+ngh, ilen, 1, klen, MPI_INTEGER, face_f(1))
  call para_type_block3a(iini-ngh, iend+ngh, jini-ngh, jend+ngh, ilen, 2, klen, MPI_INTEGER, face_f(2))
  call MPI_CART_SHIFT(MPI_COMM_CART, 1, 1, srcR_f, destR_f, ierr)
  call MPI_CART_SHIFT(MPI_COMM_CART, 1,-1, srcL_f, destL_f, ierr)
  endif

  call MPI_IRECV(F(iini-ngh , jini-ngh   , kini-ngh),1,face_f(ngh),srcR_f ,0,MPI_COMM_Cart,req(1),ierr)
  call MPI_ISEND(F(iini-ngh , jend-ngh+1 , kini-ngh),1,face_f(ngh),destR_f,0,MPI_COMM_Cart,req(2),ierr)
  call MPI_IRECV(F(iini-ngh , jend+1     , kini-ngh),1,face_f(ngh),srcL_f ,0,MPI_COMM_Cart,req(3),ierr)
  call MPI_ISEND(F(iini-ngh , jini       , kini-ngh),1,face_f(ngh),destL_f,0,MPI_COMM_Cart,req(4),ierr)
end subroutine ighost_y_fine


subroutine ighost_z_fine(F,ngh,req)
  use module_rudman_init , only: idx_f
  use module_grid
  implicit none
  include 'mpif.h'
  integer, intent(in) :: ngh
  integer, intent(out) :: req(4)
  integer, dimension(idx_f(1,1):idx_f(1,2),idx_f(2,1):idx_f(2,2),idx_f(3,1):idx_f(3,2)), intent(inout) :: F
  integer :: ilen, jlen, ierr !,sta(MPI_STATUS_SIZE,4)
  integer, save :: srcL_f, srcR_f, destL_f, destR_f, face_f(2)
  logical, save :: first_time_f = .true. 
  integer :: iini,iend,jini,jend,kini,kend

  iini=idx_f(1,1) + ngh ; iend=idx_f(1,2) - ngh 
  jini=idx_f(2,1) + ngh ; jend=idx_f(2,2) - ngh
  kini=idx_f(3,1) + ngh ; kend=idx_f(3,2) - ngh 
 
  if(first_time_f) then
   first_time_f = .false.
  ilen=idx_f(1,2)-idx_f(1,1)+1; jlen=idx_f(2,2)-idx_f(2,1)+1; !klen=ngh
  call para_type_block3a(iini-ngh, iend+ngh, jini-ngh, jend+ngh, ilen, jlen, 1, MPI_INTEGER, face_f(1))
  call para_type_block3a(iini-ngh, iend+ngh, jini-ngh, jend+ngh, ilen, jlen, 2, MPI_INTEGER, face_f(2))
  call MPI_CART_SHIFT(MPI_COMM_CART, 2, 1, srcR_f, destR_f, ierr)
  call MPI_CART_SHIFT(MPI_COMM_CART, 2,-1, srcL_f, destL_f, ierr)
  endif

  call MPI_IRECV(F(iini-ngh , jini-ngh , kini-ngh  ),1,face_f(ngh),srcR_f ,0,MPI_COMM_Cart,req(1),ierr)
  call MPI_ISEND(F(iini-ngh , jini-ngh , kend-ngh+1),1,face_f(ngh),destR_f,0,MPI_COMM_Cart,req(2),ierr)
  call MPI_IRECV(F(iini-ngh , jini-ngh , kend+1    ),1,face_f(ngh),srcL_f ,0,MPI_COMM_Cart,req(3),ierr)
  call MPI_ISEND(F(iini-ngh , jini-ngh , kini      ),1,face_f(ngh),destL_f,0,MPI_COMM_Cart,req(4),ierr)
end subroutine ighost_z_fine


subroutine init_ghost_fine(srcL_f, srcR_f, destL_f, destR_f, face_f, dir)
  use module_rudman_init , only : idx_f
  use module_grid
  implicit none
  include 'mpif.h'
  integer :: ilen, jlen, klen, ierr, dir
  integer :: srcL_f, srcR_f, destL_f, destR_f, face_f(2)
  integer :: iini,iend,jini,jend,kini,kend

  iini=idx_f(1,1) ; iend=idx_f(1,2) 
  jini=idx_f(2,1) ; jend=idx_f(2,2)
  kini=idx_f(3,1) ; kend=idx_f(3,2) 

  SELECT CASE (dir)
  CASE (1)
      jlen=jend-jini+1; klen=kend-kini+1; !ilen=ngh
      call para_type_block3a(iini, iend, jini, jend, 1, jlen, klen, MPI_DOUBLE_PRECISION, face_f(1))
      call para_type_block3a(iini, iend, jini, jend, 2, jlen, klen, MPI_DOUBLE_PRECISION, face_f(2))
      call MPI_CART_SHIFT(MPI_COMM_CART, 0, 1, srcR_f, destR_f, ierr)
      call MPI_CART_SHIFT(MPI_COMM_CART, 0,-1, srcL_f, destL_f, ierr)
  CASE (2)
        klen=kend-kini+1; ilen=iend-iini+1; !jlen=ngh
        call para_type_block3a(iini, iend, jini, jend, ilen, 1, klen, MPI_DOUBLE_PRECISION, face_f(1))
        call para_type_block3a(iini, iend, jini, jend, ilen, 2, klen, MPI_DOUBLE_PRECISION, face_f(2))
        call MPI_CART_SHIFT(MPI_COMM_CART, 1, 1, srcR_f, destR_f, ierr)
        call MPI_CART_SHIFT(MPI_COMM_CART, 1,-1, srcL_f, destL_f, ierr)
  CASE (3)
        ilen=iend-iini+1; jlen=jend-jini+1; !klen=ngh
        call para_type_block3a(iini, iend, jini, jend, ilen, jlen, 1, MPI_DOUBLE_PRECISION, face_f(1))
        call para_type_block3a(iini, iend, jini, jend, ilen, jlen, 2, MPI_DOUBLE_PRECISION, face_f(2))
        call MPI_CART_SHIFT(MPI_COMM_CART, 2, 1, srcR_f, destR_f, ierr)
        call MPI_CART_SHIFT(MPI_COMM_CART, 2,-1, srcL_f, destL_f, ierr)
  END SELECT 

end subroutine init_ghost_fine



!-------------------------------------------END GHOST OPERATIONS-----------------------------------------------

!--------------------------------------------INDEX INTERCHANGE OPERATIONS-------------------------------------

subroutine update_bounds_fine() ! updates fine grid bounds for Nx,Ny,Nz,Nxt.Nyt,Nzt
   use module_grid
   implicit none
    
    Nx = 2*Nx ; Ny = 2*Ny ; Nz = 2*Nz 
    Nxt = Nx+2*Ng ; Nyt = Ny+2*Ng ; Nzt = Nz+2*Ng 
   
end subroutine update_bounds_fine


subroutine update_index(n)  ! general setup routine for indices
   use module_grid
   use module_bc
   implicit none
   integer, intent(in) :: n(3)

      is=coords(1)*n(1)+1+Ng; imin=is-Ng
      js=coords(2)*n(2)+1+Ng; jmin=js-Ng
      ks=coords(3)*n(3)+1+Ng; kmin=ks-Ng
      ie = is + n(1) - 1
      je = js + n(2) - 1
      ke = ks + n(3) - 1
      imax = ie + Ng
      jmax = je + Ng
      kmax = ke + Ng
      ieu=ie; if(bdry_cond(1)/=1 .and. coords(1)==nPx-1) ieu=ie-1
      jev=je; if(bdry_cond(2)/=1 .and. coords(2)==nPy-1) jev=je-1
      kew=ke; if(bdry_cond(3)/=1 .and. coords(3)==nPz-1) kew=ke-1
    
end subroutine update_index


subroutine update_index_fine() ! updates indices is,ie,imin,imax etc, Nx etc are setup for fine grid
   use module_grid
   implicit none
   integer :: tmp(3)

    call update_bounds_fine()

    Mx = Nx/nPx ; My = Ny/nPy ; Mz = Nz/nPz
    tmp(1) = Mx ; tmp(2) = My ; tmp(3) = Mz
    
    call update_index(tmp) 
   
end subroutine update_index_fine


subroutine init_idx()
   use module_rudman_init
   implicit none

   call update_index_fine()
     
   idx_f(1,1) = imin ; idx_f(1,2) = imax
   idx_f(2,1) = jmin ; idx_f(2,2) = jmax
   idx_f(3,1) = kmin ; idx_f(3,2) = kmax

   call revert_index_coarse()

   idx_c(1,1) = imin ; idx_c(1,2) = imax
   idx_c(2,1) = jmin ; idx_c(2,2) = jmax
   idx_c(3,1) = kmin ; idx_c(3,2) = kmax
   

end subroutine init_idx


subroutine switchgrid(n)
   use module_rudman_init
   implicit none
   integer, intent(in) :: n  ! n=1 means switch to coarse, n=2 means switch to fine
 ! change indices to fine grid and setup fine grid coordinates

     deallocate(x, xh, dx, dxh, &
                y, yh, dy, dyh, &
                z, zh, dz, dzh  ) 
    
     if(n==2) then
       call update_index_fine() 
       subgrid  = .true. 
     elseif(n==1) then
       call revert_index_coarse()
       subgrid = .false.
     else 
       print * , "unknown grid"
     endif
    
     allocate(x(Nxt), xh(Nxt), dx(Nxt), dxh(Nxt), &
              y(Nyt), yh(Nyt), dy(Nyt), dyh(Nyt), & 
              z(Nzt), zh(Nzt), dz(Nzt), dzh(Nzt)  )
    
     call setup_coord()

end subroutine switchgrid


subroutine switchvof(flag)
   use module_rudman_init
   use module_vof
   implicit none
   integer , intent(in) :: flag 
   integer , dimension(3,2) :: n = 0 

   ! flag == 1 : to change to coarse grid containers
   ! flag == 2 : to change to fine grid containers

     if(flag==1) then
       n = idx_c
     elseif(flag==2) then
       n = idx_f
     else
      print * , "error in vof container index change"
     endif

     deallocate(cvof,vof_flag,vof_phase,tmp_flag) 
     
     allocate(cvof(n(1,1):n(1,2),n(2,1):n(2,2),n(3,1):n(3,2)), &
              vof_flag(n(1,1):n(1,2),n(2,1):n(2,2),n(3,1):n(3,2)), &
              vof_phase(n(1,1):n(1,2),n(2,1):n(2,2),n(3,1):n(3,2)), & 
              tmp_flag(n(1,1):n(1,2),n(2,1):n(2,2),n(3,1):n(3,2)) )

     cvof = 0.d0 ; vof_flag = 0 ; vof_phase = 1 ; tmp_flag = 0

end subroutine switchvof


!-------------------------------------------------------fine --> coarse ---------------------------------
subroutine revert_bounds_coarse() ! reverts to coarse grid bounds for Nx,Ny,Nz,Nxt.Nyt,Nzt
   use module_grid
   implicit none
  
    Nx = Nx/2 ; Ny = Ny/2 ; Nz = Nz/2 
    Nxt = Nx+2*Ng ; Nyt = Ny+2*Ng ; Nzt = Nz+2*Ng 
   
end subroutine revert_bounds_coarse


subroutine revert_index_coarse() ! reverts indices is,ie,imin,imax etc to coarse grid 
   use module_grid
   implicit none
   integer :: tmp(3) 
  
    call revert_bounds_coarse() ! resets Nx,Nxt etc to coarse grid 
  
    Mx = Nx/nPx ; My = Ny/nPy ; Mz = Nz/nPz
    tmp(1) = Mx ; tmp(2) = My ; tmp(3) = Mz
   
    call update_index(tmp) ! reverts indices to coarse grid
   
end subroutine revert_index_coarse
!-------------------------------------------------------END INDEX OPERATIONS-----------------------------


!-----------------------------------------------SUBGRID VOF ADVECTION ROUTINES-----------------------------
subroutine compute_vel_coarse2fine(d,uf,uc)
   use module_grid
   use module_rudman_init, only: idx_c,idx_f
   implicit none
   include 'mpif.h'
   integer, intent(in) :: d
   real(8), dimension(idx_f(1,1):idx_f(1,2),idx_f(2,1):idx_f(2,2),idx_f(3,1):idx_f(3,2)), intent(out) :: uf  ! subgrid velocity
   real(8), dimension(idx_c(1,1):idx_c(1,2),idx_c(2,1):idx_c(2,2),idx_c(3,1):idx_c(3,2)), intent(in) :: uc   ! coarse velocity
   integer :: i,j,k,isc,iec,jsc,jec,ksc,kec,isf,ief,jsf,jef,ksf,kef,i1,j1,k1

   ! COARSE grid starting indices  
    isc = idx_c(1,1) + Ng ; jsc = idx_c(2,1) + Ng ; ksc = idx_c(3,1) + Ng
    iec = idx_c(1,2) - Ng ; jec = idx_c(2,2) - Ng ; kec = idx_c(3,2) - Ng

   ! FINE grid starting indices
    isf = idx_f(1,1) + Ng ; jsf = idx_f(2,1) + Ng ; ksf = idx_f(3,1) + Ng
    ief = idx_f(1,2) - Ng ; jef = idx_f(2,2) - Ng ; kef = idx_f(3,2) - Ng

    uf = 0.d0

    select case(d)
!----------------------------------------------------X projection------------------------------------------
     case(1)   
   ! for proper indexing reference check notes

   ! ASSIGNMENT
      do i=isc-1,iec ; do j=jsc,jec ; do k=ksc,kec 
            do j1=0,1 ; do k1=0,1  
               uf( 2*(i-isc) + (isf+1) , 2*(j-jsc)+jsf+j1 , 2*(k-ksc)+ksf+k1 ) = uc(i,j,k)
            enddo ; enddo
      enddo ; enddo ; enddo


   ! AVERAGING 
      do i=isc-1,iec-1 ; do j=jsc,jec ; do k=ksc,kec 
            do j1=0,1 ; do k1=0,1  
               uf( 2*(i-isc) + (isf+2) , 2*(j-jsc)+jsf+j1 , 2*(k-ksc)+ksf+k1 ) = 0.5d0*(uc(i,j,k)+uc(i+1,j,k))
            enddo ; enddo
      enddo ; enddo ; enddo

!-----------------------------------------------------Y projection-----------------------------------------

     case(2)

  ! ASSIGNMENT
      do i=isc,iec ; do j=jsc-1,jec ; do k=ksc,kec 
            do i1=0,1 ; do k1=0,1  
               uf( 2*(i-isc)+isf+i1 , 2*(j-jsc)+(jsf+1) , 2*(k-ksc)+ksf+k1 ) = uc(i,j,k)
            enddo ; enddo
      enddo ; enddo ; enddo

   ! AVERAGING 
      do i=isc,iec ; do j=jsc-1,jec-1 ; do k=ksc,kec 
            do i1=0,1 ; do k1=0,1  
               uf( 2*(i-isc)+isf+i1 , 2*(j-jsc) + (jsf+2) , 2*(k-ksc)+ksf+k1 ) = 0.5d0*(uc(i,j,k)+uc(i,j+1,k))
            enddo ; enddo
      enddo ; enddo ; enddo

!-----------------------------------------------------Z projection-----------------------------------------

     case(3)

  ! ASSIGNMENT
      do i=isc,iec ; do j=jsc,jec ; do k=ksc-1,kec 
            do i1=0,1 ; do j1=0,1  
               uf( 2*(i-isc) +isf+i1 , 2*(j-jsc)+jsf+j1 , 2*(k-ksc)+ (ksf+1) ) = uc(i,j,k)
            enddo ; enddo
      enddo ; enddo ; enddo

   ! AVERAGING 
      do i=isc,iec ; do j=jsc,jec ; do k=ksc-1,kec-1 
            do i1=0,1 ; do j1=0,1  
               uf( 2*(i-isc)+isf+i1 , 2*(j-jsc)+jsf+j1 , 2*(k-ksc) + (ksf+2) ) = 0.5d0*(uc(i,j,k)+uc(i,j,k+1))
            enddo ; enddo
      enddo ; enddo ; enddo

     case default
       print * , "error projecting velocity to subgrid !! "   

    end select

end subroutine compute_vel_coarse2fine


subroutine advect_vof_fine_centered(us,vf,vf_flag,d,mask,dtf)
   use module_rudman_init
   use module_rudman_ghost
   implicit none
   include 'mpif.h'
   integer :: i,j,k
   integer :: i0,j0,k0
   integer :: i1,j1,k1
   integer :: iend,jend,kend
   integer :: iini,jini,kini
   integer , intent(in) :: d ! sweep direction
   real(8) , intent(inout) :: dtf ! halved time-step for subgrid
   real(8) , dimension(idx_f(1,1):idx_f(1,2),idx_f(2,1):idx_f(2,2),idx_f(3,1):idx_f(3,2)), intent(in) :: us ! subgrid velocity 
   real(8) , dimension(idx_f(1,1):idx_f(1,2),idx_f(2,1):idx_f(2,2),idx_f(3,1):idx_f(3,2)), intent(inout) :: vf !subgrid vof 
   integer , dimension(idx_f(1,1):idx_f(1,2),idx_f(2,1):idx_f(2,2),idx_f(3,1):idx_f(3,2)), intent(inout) :: vf_flag !subgrid flag 
   real(8) , dimension(idx_f(1,1):idx_f(1,2),idx_f(2,1):idx_f(2,2),idx_f(3,1):idx_f(3,2)), intent(in) :: mask !subgrid WY mask 
   real(8) , dimension(idx_f(1,1):idx_f(1,2),idx_f(2,1):idx_f(2,2),idx_f(3,1):idx_f(3,2)) :: vof1 ! vof1 fluxes fine
   real(8) , dimension(idx_f(1,1):idx_f(1,2),idx_f(2,1):idx_f(2,2),idx_f(3,1):idx_f(3,2)) :: vof3 ! vof3 fluxes fine
   real(8) :: dxyz,vof
   real(8) :: a1,a2,alpha
   real(8) :: al3d, fl3d, x0(3), deltax(3)
   real(8) :: mxyz(3),stencil3x3(-1:1,-1:1,-1:1)
   intrinsic dmax1,dmin1

   ! time-step needs to be halved as we are computing fine grid fluxes. 

   i0=0; j0=0; k0=0 ; vf_flag = 0

   iini=idx_f(1,1) + Ng ; iend=idx_f(1,2) - Ng ! limits for do loop initially set to fine grid
   jini=idx_f(2,1) + Ng ; jend=idx_f(2,2) - Ng
   kini=idx_f(3,1) + Ng ; kend=idx_f(3,2) - Ng 

   dxyz = cell_size(fine) ! SUBGRID SIZE 

   if (d == 1) then
      i0=1  
   else if (d == 2) then
      j0=1  
   else if (d == 3) then
      k0=1  
   endif

   ! Eulerian fluxes calculated based on previous time-step values of vof 

   !do i=iini-i0,iend+i0 ; do j=jini-j0,jend+j0 ; do k=kini-k0,kend+k0
   do i=iini-1,iend+1 ; do j=jini-1,jend+1 ; do k=kini-1,kend+1

      a2 = us(i,j,k)*(dtf/dxyz)
      a1 = us(i-i0,j-j0,k-k0)*(dtf/dxyz)

      vof1(i,j,k) = 0.d0
      vof3(i,j,k) = 0.d0

      if (vf(i,j,k) == 1.0d0) then   ! for pure phase 2 (bulk)

          vof1(i,j,k) = DMAX1(-a1,0.d0)
          vof3(i,j,k) = DMAX1(a2,0.d0)

      else if ((vf(i,j,k) .gt. 0.d0) .and. (vf(i,j,k) .lt. 1.d0)) then 

          do i1=-1,1; do j1=-1,1; do k1=-1,1
             stencil3x3(i1,j1,k1) = vf(i+i1,j+j1,k+k1)
          enddo;enddo;enddo

          call mycs(stencil3x3,mxyz)
          alpha = al3d(mxyz,vf(i,j,k))

          x0=0d0
          deltax=1d0
 
          if(a1<0d0) then
              deltax(d)=-a1
              vof1(i,j,k) = fl3d(mxyz,alpha,x0,deltax)
          endif

          if(a2>0d0) then
              x0(d)=1d0-a2
              deltax(d)=a2
              vof3(i,j,k) = fl3d(mxyz,alpha,x0,deltax)
          endif

      endif

   enddo ; enddo ; enddo 

   call do_all_ghost_fine(vof1)  
   call do_all_ghost_fine(vof3)

   do i=iini,iend ; do j=jini,jend ; do k=kini,kend  

     a2 = us(i,j,k)*(dtf/dxyz)
     a1 = us(i-i0,j-j0,k-k0)*(dtf/dxyz)

     vf(i,j,k) = vf(i,j,k) - (vof3(i,j,k) - vof1(i+i0,j+j0,k+k0)) + & 
                  (vof3(i-i0,j-j0,k-k0) - vof1(i,j,k)) + mask(i,j,k)*(a2-a1)

   enddo ; enddo ; enddo

   call rudman_vof_flags_and_clip(vf,vf_flag,idx_f)
   call vof_bc_rudman(vf,vf_flag,idx_f,fine)

   call do_all_ghost_fine(vf)
   call do_all_ighost_fine(vf_flag)

end subroutine advect_vof_fine_centered
!--------------------------------------------VOF advection routines end------------------------------------


!--------------------------------------------------INIT ROUTINES--------------------------------------------

subroutine rudman_init_vof_pre()
   use module_rudman_init
   use module_vof
   implicit none
  include 'mpif.h'

   call switchgrid(2) ! switching to fine grid       
   call switchvof(2) ! allocating cvof containers to fine grid
   cvof = 0.d0 ; vof_flag = 3

end subroutine rudman_init_vof_pre


subroutine rudman_init_vof_post()
   use module_rudman_init
   use module_rudman_ghost
   use module_vof
   implicit none
  include 'mpif.h'
   real(8) :: test

   call rudman_vof_flags_and_clip(cvof,vof_flag,idx_f)
   call vof_bc_rudman(cvof,vof_flag,idx_f,fine)
  
   call do_all_ghost_fine(cvof)
   call do_all_ighost_fine(vof_flag)

   call init_rudman_vof() ! initializes cvof_f and vof_flag_fine containers
   cvof_f = cvof ; vof_flag_f = vof_flag ! cvof_f and vof_flag_fine assigned fine grid values

   call switchgrid(1) ! indices and co-ordinate system reverted to coarse grid
   call switchvof(1) ! vof related containers allocated to coarse grid 

   call compute_vof_fine2coarse(cvof,vof_flag,cvof_f) ! fine to coarse averaging routine

end subroutine rudman_init_vof_post

!----------------------------------------------------------Initialization routines end----------------------------


!--------------------------------------------------------Subgrid Eulerian advection ends--------------------------

subroutine rudman_advect_vof_mom_stg(tswap,t)
   use module_rudman_init
   use module_rudman_ghost
   use module_rudman_advection
   use module_vof
   implicit none
   include 'mpif.h'
   integer , intent(in) :: tswap
   real(8) , intent(in) :: t
   integer :: i,j,k,iini,iend,jini,jend,kini,kend 

   iini=idx_c(1,1) + Ng ; iend=idx_c(1,2) - Ng
   jini=idx_c(2,1) + Ng ; jend=idx_c(2,2) - Ng
   kini=idx_c(3,1) + Ng ; kend=idx_c(3,2) - Ng 

   call init_adv(cvof) ! allocation for vel_fine --> fine grid velocity containers

   call project_coarse2fine() ! coarse grid velocities are projected to vel_fine 

   call vofandmomsweepsstaggered_fine(tswap,t) ! Subgrid mom-vof coupled Eulerian WY Advection

   call vofsweeps_fine(tswap,t) ! Subgrid vof Eulerian WY advection
   
   call compute_vof_fine2coarse(cvof,vof_flag,cvof_f) ! averaging routine for fine to coarse grid vof

end subroutine rudman_advect_vof_mom_stg
!-------------------------------------------------MOM VOF Advection ends---------------------------------


!------------------------------------------------PROJECTION FOR RUDMAN ROUTINE--------------------------------------

subroutine project_velocity_rudman()
   use module_rudman_init , only : idx_c
   use module_flow
   use module_grid
   implicit none
   integer :: i,j,k
   real(8) , dimension(idx_c(1,1):idx_c(1,2),idx_c(2,1):idx_c(2,2),idx_c(3,1):idx_c(3,2)) :: rhostg

   rhostg = 0.d0 

   call compute_vof_coarse_stg(rhostg,1)

   do k=ks,ke;  do j=js,je; do i=is,ieu    ! CORRECT THE u-velocity 
      u(i,j,k)=u(i,j,k)-dt*(umask(i,j,k)/dxh(i))*(p(i+1,j,k)-p(i,j,k))/rhostg(i,j,k)
   enddo; enddo; enddo

   call compute_vof_coarse_stg(rhostg,2)

   do k=ks,ke;  do j=js,jev; do i=is,ie    ! CORRECT THE v-velocity
      v(i,j,k)=v(i,j,k)-dt*(vmask(i,j,k)/dyh(j))*(p(i,j+1,k)-p(i,j,k))/rhostg(i,j,k)
   enddo; enddo; enddo

   call compute_vof_coarse_stg(rhostg,3)

   do k=ks,kew;  do j=js,je; do i=is,ie   ! CORRECT THE w-velocity
      w(i,j,k)=w(i,j,k)-dt*(wmask(i,j,k)/dzh(k))*(p(i,j,k+1)-p(i,j,k))/rhostg(i,j,k)
   enddo; enddo; enddo

end subroutine project_velocity_rudman










