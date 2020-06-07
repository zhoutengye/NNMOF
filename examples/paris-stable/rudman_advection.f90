!---------------------------------------coarse grid staggered vof for projection---------------------------------

subroutine compute_vof_coarse_stg(rho_stg,dir)
   use module_rudman_init
   use module_flow
   use module_bc
   implicit none
   include 'mpif.h'
   integer , intent(in) :: dir 
   real(8) , dimension(idx_c(1,1):idx_c(1,2),idx_c(2,1):idx_c(2,2),idx_c(3,1):idx_c(3,2)), intent(out) :: rho_stg
   real(8) , dimension(idx_c(1,1):idx_c(1,2),idx_c(2,1):idx_c(2,2),idx_c(3,1):idx_c(3,2)) :: vof 
   integer , dimension(idx_c(1,1):idx_c(1,2),idx_c(2,1):idx_c(2,2),idx_c(3,1):idx_c(3,2)) :: v_flag
   real(8) :: vofavg
   integer :: iend,jend,kend
   integer :: iini,jini,kini
   integer :: i,j,k,isf,jsf,ksf,i1,j1,k1
   integer :: ii,jj,kk,ip,jp,kp

  ! mom and vof are in staggered configuration , FINE VOF is used to calculate COARSE STAGGERED mom and vof
  ! vel(:,:,:,1) --> u , vel(:,:,:,2) --> v , vel(:,:,:,3) --> w 
  ! vel(:,:,:,4) --> umask , vel(:,:,:,5) --> vmask , vel(:,:,:,6) --> wmask 

   vof = 0.d0 ; v_flag = 3 
   rho_stg = 0.d0 ; vofavg = 0.d0 
   ii = 0 ; jj = 0 ; kk = 0 
   ip = 0 ; jp = 0 ; kp = 0

   ! Coarse grid indices for the do loops

   iini=idx_c(1,1) + Ng ; iend=idx_c(1,2) - Ng
   jini=idx_c(2,1) + Ng ; jend=idx_c(2,2) - Ng
   kini=idx_c(3,1) + Ng ; kend=idx_c(3,2) - Ng 

   isf = idx_f(1,1)+Ng ; jsf = idx_f(2,1)+Ng ; ksf = idx_f(3,1)+Ng ! Subgrid Indices

   if (dir==1) then
      ii = 1 
        if (bdry_cond(1)/=1 .and. coords(1)==nPx-1) ip=1
   elseif (dir==2) then
      jj = 1 
        if (bdry_cond(2)/=1 .and. coords(2)==nPy-1) jp=1
   elseif (dir==3) then
      kk = 1 
        if (bdry_cond(3)/=1 .and. coords(3)==nPz-1) kp=1
   endif

   do i=iini,iend-ip ; do j=jini,jend-jp ; do k=kini,kend-kp   
      do i1=0,1 ;  do j1=0,1 ;  do k1=0,1
        vof(i,j,k) = vof(i,j,k) + &    ! Staggered Coarse VOF
                     cvof_f( 2*(i-iini)+isf+i1+ii , 2*(j-jini)+jsf+j1+jj , 2*(k-kini)+ksf+k1+kk )*0.125d0 
      enddo ; enddo ; enddo                         
   enddo ; enddo ; enddo


   call rudman_vof_flags_and_clip(vof,v_flag,idx_c) ! assign staggered flags
   call vof_bc_rudman(vof,v_flag,idx_c,coarse)  ! use flags to assign stag vof bc's

   call do_all_ghost(vof)
   call do_all_ighost(v_flag)


   do i=iini,iend-ip ; do j=jini,jend-jp ; do k=kini,kend-kp   
      vofavg = vof(i,j,k)
      rho_stg(i,j,k) = vofavg*rho2 + (1.d0-vofavg)*rho1
   enddo ; enddo ; enddo 

   call do_all_ghost(rho_stg)

!   do i=iini,iend-ip ; do j=jini,jend-jp ; do k=kini,kend-kp   
!      if(rho_stg(i,j,k) .ne. rho1 .and. rho_stg(i,j,k) .ne. rho2) then
!          print * , rho_stg(i,j,k) , i , j ,k 
!              !  if(rho_stg(i,j,k) < 1.d-12) print * , "Function" , i,j,k
!     endif
!   enddo ; enddo ; enddo 



end subroutine compute_vof_coarse_stg


!-----------------------------------------------------------------------------------------------------------------

!------------------------------------coarse stg mom and vof calculation-------------------------------------------
subroutine compute_mom_vof_coarse_stg(mom,vof,v_flag,dir)
   use module_rudman_init
   use module_flow
   use module_bc
   implicit none
   include 'mpif.h'
   integer , intent(in) :: dir 
   real(8) , dimension(idx_c(1,1):idx_c(1,2),idx_c(2,1):idx_c(2,2),idx_c(3,1):idx_c(3,2)), intent(inout) :: mom
   real(8) , dimension(idx_c(1,1):idx_c(1,2),idx_c(2,1):idx_c(2,2),idx_c(3,1):idx_c(3,2)), intent(inout) :: vof
   integer , dimension(idx_c(1,1):idx_c(1,2),idx_c(2,1):idx_c(2,2),idx_c(3,1):idx_c(3,2)), intent(inout) :: v_flag
   real(8) , dimension(idx_c(1,1):idx_c(1,2),idx_c(2,1):idx_c(2,2),idx_c(3,1):idx_c(3,2)) :: vel 
   real(8) :: vofavg,rhoavg
   integer :: iend,jend,kend
   integer :: iini,jini,kini
   integer :: i,j,k,isf,jsf,ksf
   integer :: ii,jj,kk,i1,j1,k1,i2,j2,k2 
   integer :: ilim, jlim, klim
   real(8) :: fac

  ! mom and vof are in staggered configuration , FINE VOF is used to calculate COARSE STAGGERED mom and vof
  ! vel(:,:,:,1) --> u , vel(:,:,:,2) --> v , vel(:,:,:,3) --> w 
  ! vel(:,:,:,4) --> umask , vel(:,:,:,5) --> vmask , vel(:,:,:,6) --> wmask 

   mom = 0.d0 ; vof = 0.d0 ; v_flag = 0 ; vel =0.d0
   rhoavg = 0.d0 ; vofavg = 0.d0  
   ii = 0 ; jj = 0 ; kk = 0
   i2 = 0 ; j2 = 0 ; k2 = 0
   ilim = 0 ; jlim = 0 ; klim = 0

   ! Coarse grid indices for the do loops

   iini=idx_c(1,1) + Ng ; iend=idx_c(1,2) - Ng
   jini=idx_c(2,1) + Ng ; jend=idx_c(2,2) - Ng
   kini=idx_c(3,1) + Ng ; kend=idx_c(3,2) - Ng 

   isf = idx_f(1,1)+Ng ; jsf = idx_f(2,1)+Ng ; ksf = idx_f(3,1)+Ng ! Subgrid Indices

   if (dir==1) then
      ii = 1 ; vel = u 
   elseif (dir==2) then
      jj = 1 ; vel = v  
   elseif (dir==3) then
      kk = 1 ; vel = w 
   endif

   call do_all_ghost(vel)

   select case(dir)

    case(1)

       do i=iini-1,iend ; do j=jini,jend ; do k=kini,kend  

       if(i < iini) then
          ilim = 0 ; fac = 0.25d0 ; i2 = -1
       elseif(i == iend) then
          ilim = 2*(i-iini)+1 ; fac = 0.25d0 ; i2 = -1    
       else 
          ilim = 2*(i-iini)+1 ; fac = 0.125d0 ; i2 = 0
       endif     
               
          do i1=0,1+i2 ;  do j1=0,1 ;  do k1=0,1
            vof(i,j,k) = vof(i,j,k) + &    ! Staggered Coarse VOF
                         cvof_f( ilim+isf+i1 , 2*(j-jini)+jsf+j1, 2*(k-kini)+ksf+k1 )*fac 
          enddo ; enddo ; enddo                         

       enddo ; enddo ; enddo

    case(2)

       do i=iini,iend ; do j=jini-1,jend ; do k=kini,kend  

       if(j < jini) then
          jlim = 0 ; fac = 0.25d0 ; j2 = -1
       elseif(j == jend) then
          jlim = 2*(j-jini)+1 ; fac = 0.25d0 ; j2 = -1    
       else 
          jlim = 2*(j-jini)+1 ; fac = 0.125d0 ; j2 = 0
       endif     
               
          do i1=0,1 ;  do j1=0,1+j2 ;  do k1=0,1
            vof(i,j,k) = vof(i,j,k) + &    ! Staggered Coarse VOF
                         cvof_f( 2*(i-iini)+isf+i1 , jlim+jsf+j1, 2*(k-kini)+ksf+k1 )*fac 
          enddo ; enddo ; enddo                         

       enddo ; enddo ; enddo
    
    case(3)

       do i=iini,iend ; do j=jini,jend ; do k=kini-1,kend  

       if(k < kini) then
          klim = 0 ; fac = 0.25d0 ; k2 = -1
       elseif(k == kend) then
          klim = 2*(k-kini)+1 ; fac = 0.25d0 ; k2 = -1    
       else 
          klim = 2*(k-kini)+1 ; fac = 0.125d0 ; k2 = 0
       endif     
               
          do i1=0,1 ;  do j1=0,1 ;  do k1=0,1+k2
            vof(i,j,k) = vof(i,j,k) + &    ! Staggered Coarse VOF
                         cvof_f( 2*(i-iini)+isf+i1 , 2*(j-jini)+jsf+j1, klim+ksf+k1 )*fac 
          enddo ; enddo ; enddo                         

       enddo ; enddo ; enddo

    case default
            print * , "Momentum Initialization Error"

   end select


   call rudman_vof_flags_and_clip(vof,v_flag,idx_c) ! assign staggered flags
   call vof_bc_rudman(vof,v_flag,idx_c,coarse)  ! use flags to assign stag vof bc's

   call do_all_ghost(vof)
   call do_all_ighost(v_flag)

 
   do i=iini-1,iend+1 ; do j=jini-1,jend+1 ; do k=kini-1,kend+1   
      vofavg = vof(i,j,k)
      rhoavg = vofavg*rho2 + (1.d0-vofavg)*rho1
      mom(i,j,k) = vel(i,j,k)*rhoavg   ! Staggered Coarse MOM 
   enddo ; enddo ; enddo


   if (dir==1) then
     call mom_bc_rudman(u,vof,mom,1,umask,rho1,rho2,0.d0,idx_c)
   elseif (dir==2) then
     call mom_bc_rudman(v,vof,mom,2,vmask,rho1,rho2,0.d0,idx_c)
   elseif (dir==3) then
     call mom_bc_rudman(w,vof,mom,3,wmask,rho1,rho2,0.d0,idx_c)
   endif

   call do_all_ghost(mom)
   
end subroutine compute_mom_vof_coarse_stg
!----------------------------------------------------------------------------------------------------------------


!-----------------------------------------compute weymouth yue coefficients---------------------------------------
subroutine compute_wy_mask(vof,cc,n)
   implicit none
   include 'mpif.h'
   integer , dimension(3,2) , intent(in) :: n
   real(8) , dimension(n(1,1):n(1,2),n(2,1):n(2,2),n(3,1):n(3,2)), intent(in) :: vof 
   real(8) , dimension(n(1,1):n(1,2),n(2,1):n(2,2),n(3,1):n(3,2)), intent(out) :: cc
   integer :: i,j,k

   cc = 0.d0

   do i=n(1,1),n(1,2) ; do j=n(2,1),n(2,2) ; do k=n(3,1),n(3,2)
     if(vof(i,j,k) > 0.5d0) cc(i,j,k) = 1.d0
   enddo ; enddo ; enddo

end subroutine compute_wy_mask
!------------------------------------------------------------------------------------------------------------------


!---------------------------------------subgrid to coarse grid fluxes----------------------------------------------
subroutine compute_vof_flux_coarse_stg(us,vf,flux,d,dtf)
   use module_rudman_init
   use module_rudman_ghost
   use module_bc
   implicit none
   include 'mpif.h'
   integer :: i,j,k
   integer :: isf,jsf,ksf
   integer :: i0,j0,k0
   integer :: i1,j1,k1
   integer :: i2,j2,k2
   integer :: iend,jend,kend
   integer :: iini,jini,kini
   integer , intent(in) :: d ! sweep direction
   real(8) , intent(in) :: dtf ! half the original timestep
   real(8) , dimension(idx_f(1,1):idx_f(1,2),idx_f(2,1):idx_f(2,2),idx_f(3,1):idx_f(3,2)), intent(in) :: us 
   real(8) , dimension(idx_f(1,1):idx_f(1,2),idx_f(2,1):idx_f(2,2),idx_f(3,1):idx_f(3,2)), intent(in) :: vf !subgrid vof 
   real(8) , dimension(idx_c(1,1):idx_c(1,2),idx_c(2,1):idx_c(2,2),idx_c(3,1):idx_c(3,2),2), intent(out) :: flux
   real(8) , dimension(idx_f(1,1):idx_f(1,2),idx_f(2,1):idx_f(2,2),idx_f(3,1):idx_f(3,2)) :: vof1 ! vof1 fluxes fine
   real(8) , dimension(idx_f(1,1):idx_f(1,2),idx_f(2,1):idx_f(2,2),idx_f(3,1):idx_f(3,2)) :: vof3 ! vof3 fluxes fine
   real(8) :: dxyz,vof
   real(8) :: a1,a2,alpha
   real(8) :: al3d, fl3d, x0(3), deltax(3)
   real(8) :: mxyz(3),stencil3x3(-1:1,-1:1,-1:1)
   intrinsic dmax1,dmin1
  
   ! time-step needs to be halved as we are computing fine grid fluxes. 
 
   flux = 0.d0 
   i0=0; j0=0; k0=0
   i1=0; j1=0; k1=0
   i2=0; j2=0; k2=0
   mxyz = 0.d0 ; stencil3x3 = 0.d0

   ! i0,j0,k0 are for sweep direction
   ! i2,j2,k2 are for which momentum is being swept - rho_u,rho_v,rho_w
   ! d is for the sweeping direction, dir is for which momentum is swept

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

   do i=iini,iend ; do j=jini,jend ; do k=kini,kend

      a2 = us(i,j,k)*(dtf/dxyz)
      a1 = us(i-i0,j-j0,k-k0)*(dtf/dxyz)

      vof1(i,j,k) = 0.d0
      vof3(i,j,k) = 0.d0

      if (vf(i,j,k) == 1.0d0) then   ! for pure phase 2 (bulk)

          vof1(i,j,k) = DMAX1(-a1,0.d0)
          vof3(i,j,k) = DMAX1(a2,0.d0)
      else if ((vf(i,j,k) .gt. 0.d0) .and. (vf(i,j,k) .lt. 1.d0)) then 
! Eulerian fluxes calculated based on previous time-step values of vof 

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

   call do_all_ghost_fine(vof1)  ! do not know exactly if this is needed
   call do_all_ghost_fine(vof3)

   ! do loop limits need to be changed to coarse grid values to compute coarse grid fluxes

   iini=idx_c(1,1) + Ng ; iend=idx_c(1,2) - Ng 
   jini=idx_c(2,1) + Ng ; jend=idx_c(2,2) - Ng
   kini=idx_c(3,1) + Ng ; kend=idx_c(3,2) - Ng 
   
   isf = idx_f(1,1)+Ng ; jsf = idx_f(2,1)+Ng ; ksf = idx_f(3,1)+Ng ! starting indices for fine grid

   select case(d)

   ! flux(:,:,:,1) --> coarse vof1 flux || flux(:,:,:,2) --> coarse vof3 flux
!------------------------------------------------------------------------------------------------------
    case(1) ! addition of X direction fine fluxes

      flux = 0.d0      
      do i=iini-1,iend-1 ; do j=jini,jend ; do k=kini,kend
         do j1=0,1 ; do k1=0,1

           flux(i,j,k,2) = flux(i,j,k,2) +  &  ! vof3 coarse flux
                           0.25d0*vof3( 2*(i-iini) + (isf+2) , 2*(j-jini) + jsf+j1 , 2*(k-kini) + ksf+k1 ) 

           flux(i+1,j,k,1) = flux(i+1,j,k,1) +  &  ! vof1 coarse flux
                           0.25d0*vof1( 2*(i-iini) + (isf+3) , 2*(j-jini) + jsf+j1 , 2*(k-kini) + ksf+k1 ) 

         enddo ; enddo
      enddo ; enddo ; enddo
!------------------------------------------------------------------------------------------------------
    case(2) ! addition of Y direction fine fluxes

      flux = 0.d0
      do i=iini,iend ; do j=jini-1,jend-1 ; do k=kini,kend
         do i1=0,1 ; do k1=0,1

           flux(i,j,k,2) = flux(i,j,k,2) +  &  ! vof3 coarse flux
                           0.25d0*vof3( 2*(i-iini) + isf+i1 , 2*(j-jini) + (jsf+2) , 2*(k-kini) + ksf+k1 ) 

           flux(i,j+1,k,1) = flux(i,j+1,k,1) +  &  ! vof1 coarse flux
                           0.25d0*vof1( 2*(i-iini) + isf+i1 , 2*(j-jini) + (jsf+3) , 2*(k-kini) + ksf+k1 ) 

         enddo ; enddo
      enddo ; enddo ; enddo
    
!------------------------------------------------------------------------------------------------------
    case(3) ! addition of Z direction fine fluxes

      flux = 0.d0     
      do i=iini,iend ; do j=jini,jend ; do k=kini-1,kend-1
         do i1=0,1 ; do j1=0,1

           flux(i,j,k,2) = flux(i,j,k,2) +  &  ! vof3 coarse flux
                           0.25d0*vof3( 2*(i-iini) + isf+i1 , 2*(j-jini) + jsf+j1 , 2*(k-kini) + (ksf+2) ) 

           flux(i,j,k+1,1) = flux(i,j,k+1,1) +  &  ! vof1 coarse flux
                           0.25d0*vof1( 2*(i-iini) + isf+i1 , 2*(j-jini) + jsf+j1 , 2*(k-kini) + (ksf+3) ) 

         enddo ; enddo
      enddo ; enddo ; enddo
!------------------------------------------------------------------------------------------------------
    case default
     print * , "error in flux addition from fine to coarse !!"

   end select

   call do_all_ghost(flux(:,:,:,1))
   call do_all_ghost(flux(:,:,:,2))

end subroutine compute_vof_flux_coarse_stg
!--------------------------------------------------------------------------------------------------------------------


!--------------------------------------advecting coarse stg vof and mom along direction-------------------------------
subroutine advect_mom_vof_coarse_stg(us,usmask,mom,vof,v_flag,flux,wymask,d,dir,t,dtc)
   use module_rudman_init
   use module_flow
   use module_bc
   implicit none
   include 'mpif.h'
   real(8) , intent(in) :: t,dtc ! coarse grid time-step used
   integer , intent(in) :: d,dir
   real(8) , dimension(idx_c(1,1):idx_c(1,2),idx_c(2,1):idx_c(2,2),idx_c(3,1):idx_c(3,2)), intent(in) :: us,usmask 
   real(8) , dimension(idx_c(1,1):idx_c(1,2),idx_c(2,1):idx_c(2,2),idx_c(3,1):idx_c(3,2)), intent(inout) :: mom,vof
   integer , dimension(idx_c(1,1):idx_c(1,2),idx_c(2,1):idx_c(2,2),idx_c(3,1):idx_c(3,2)), intent(inout) :: v_flag
   real(8) , dimension(idx_c(1,1):idx_c(1,2),idx_c(2,1):idx_c(2,2),idx_c(3,1):idx_c(3,2),2), intent(in) :: flux
   real(8) , dimension(idx_c(1,1):idx_c(1,2),idx_c(2,1):idx_c(2,2),idx_c(3,1):idx_c(3,2),2), intent(in) :: wymask
   real(8) , dimension(idx_c(1,1):idx_c(1,2),idx_c(2,1):idx_c(2,2),idx_c(3,1):idx_c(3,2)) :: vof1,vof3 
   real(8) , dimension(idx_c(1,1):idx_c(1,2),idx_c(2,1):idx_c(2,2),idx_c(3,1):idx_c(3,2)) :: mom1,mom3 
   integer :: i,j,k
   integer :: i0,j0,k0
   integer :: i1,j1,k1
   integer :: i2,j2,k2
   integer :: iini,jini,kini,iend,jend,kend
   real(8) :: a1, a2, uadv1, uadv3
   real(8) :: ro1, ro2, ro3
   real(8) :: u1,u2,u3
   real(8) :: dxyz
   real(8) :: test1,test2
   integer :: a,b,c
   real(8) :: alpha
   real(8) :: al3d, fl3d, x0(3), deltax(3)
   real(8) :: mxyz(3),stencil3x3(-1:1,-1:1,-1:1)
   intrinsic dmax1,dmin1

   ! i0,j0,k0 are for sweep direction
   ! i2,j2,k2 are for which momentum is being swept - rho_u,rho_v,rho_w
   ! d is for the sweeping direction, dir is for which momentum is swept

   i0=0; j0=0; k0=0
   i2=0; j2=0; k2=0
   v_flag = 1

   mom1 = 0.d0 ; mom3= 0.d0 ; vof1 = 0.d0 ; vof3 = 0.d0

   vof1(:,:,:) = flux(:,:,:,1) 
   vof3(:,:,:) = flux(:,:,:,2)

   dxyz = cell_size(coarse)  ! COARSE CELL SIZE

   if (d == 1) then
      i0=1  
   else if (d == 2) then
      j0=1  
   else if (d == 3) then
      k0=1  
   endif

   if (dir == 1) then
      i2=1  
   else if (dir == 2) then
      j2=1  
   else if (dir == 3) then
      k2=1  
   endif

!---------------------coarse grid indices set up for loop-----------------------------

   iini=idx_c(1,1) + Ng ; iend=idx_c(1,2) - Ng 
   jini=idx_c(2,1) + Ng ; jend=idx_c(2,2) - Ng
   kini=idx_c(3,1) + Ng ; kend=idx_c(3,2) - Ng 

!-------------------------coarse mom stagg. flux computation---------------------------------------------
   
   do i=iini-i2,iend+i2 ; do j=jini-j2,jend+j2 ; do k=kini-k2,kend+k2 

     a1 = 0.5d0*(us(i-i0,j-j0,k-k0)+us(i-i0+i2,j-j0+j2,k-k0+k2))*dtc/dxyz
     a2 = 0.5d0*(us(i,j,k)+us(i+i2,j+j2,k+k2))*dtc/dxyz

     ro1 = (rho2*vof(i-i0,j-j0,k-k0) + rho1*(1.d0 - vof(i-i0,j-j0,k-k0)))
     ro2 = (rho2*vof(i,j,k) + rho1*(1.d0 - vof(i,j,k)))
     ro3 = (rho2*vof(i+i0,j+j0,k+k0) + rho1*(1.d0 - vof(i+i0,j+j0,k+k0)))

     u1 = mom(i-i0,j-j0,k-k0)/ro1
     u2 = mom(i,j,k)/ro2
     u3 = mom(i+i0,j+j0,k+k0)/ro3

     uadv1 = interpole3(u1,u2,u3,AdvectionScheme,-0.5d0*(1.d0+a1))           
     uadv3 = interpole3(u1,u2,u3,AdvectionScheme,0.5d0*(1.d0-a2))           

!     uadv1 =  interpole_quad(u1,u2,u3,-0.5d0*(1.d0+a1))
!     uadv3 =  interpole_quad(u1,u2,u3,0.5d0*(1.d0-a2)) 

  !   if(i<iini .or. i>iend .or. j<jini .or. j>jend .or. k<kini .or. k>kend) then
  ! !  if(i<iini .or. i>iend .or. j<jini .or. j>jend) then

  !      if (vof(i,j,k) == 1.0d0) then   ! for pure phase 2 (bulk)
  !          vof1(i,j,k) = DMAX1(-a1,0.d0)
  !          vof3(i,j,k) = DMAX1(a2,0.d0)
  !      else if ((vof(i,j,k) .gt. 0.d0) .and. (vof(i,j,k) .lt. 1.d0)) then 

  !          do i1=-1,1; do j1=-1,1; do k1=-1,1
  !             stencil3x3(i1,j1,k1) = vof(i+i1,j+j1,k+k1)
  !          enddo;enddo;enddo

  !          call mycs(stencil3x3,mxyz)


  !          alpha = al3d(mxyz,vof(i,j,k))
  !          x0=0d0
  !          deltax=1d0
 
  !          if(a1<0d0) then
  !              deltax(d)=-a1
  !              vof1(i,j,k) = fl3d(mxyz,alpha,x0,deltax)
  !          endif

  !          if(a2>0d0) then
  !              x0(d)=1d0-a2
  !              deltax(d)=a2
  !              vof3(i,j,k) = fl3d(mxyz,alpha,x0,deltax)
  !          endif

  !      endif
  ! !         print * , "PASS", i,j,k,vof(i,j,k) 
  !
  !   endif

     mom1(i,j,k) = (rho2*vof1(i,j,k) + rho1*(-a1 - vof1(i,j,k)))*uadv1
     mom3(i,j,k) = (rho2*vof3(i,j,k) + rho1*( a2 - vof3(i,j,k)))*uadv3


   enddo ; enddo ; enddo 

   call do_all_ghost(mom1)
   call do_all_ghost(mom3)  

   
   
!-------------------------advection of vof and momentum-------------------------------------------------


   do i=iini,iend ; do j=jini,jend ; do k=kini,kend  

     a2 = 0.5d0*(us(i,j,k)+us(i+i2,j+j2,k+k2))*dtc/dxyz
     a1 = 0.5d0*(us(i-i0,j-j0,k-k0)+us(i-i0+i2,j-j0+j2,k-k0+k2))*dtc/dxyz


     vof(i,j,k) = vof(i,j,k) - (vof3(i,j,k) - vof1(i+i0,j+j0,k+k0)) + & 
                  (vof3(i-i0,j-j0,k-k0) - vof1(i,j,k)) + wymask(i,j,k,1)*(a2-a1)


     mom(i,j,k)= mom(i,j,k) -  (mom3(i,j,k) - mom1(i+i0,j+j0,k+k0)) + & 
                 (mom3(i-i0,j-j0,k-k0) - mom1(i,j,k)) + &
                 (rho2-rho1)*wymask(i,j,k,1)*wymask(i,j,k,2)*(a2-a1)


   enddo ; enddo ; enddo


   ! Coarse grid indices passed as arguments in BC implementation

   call rudman_vof_flags_and_clip(vof,v_flag,idx_c)
   call vof_bc_rudman(vof,v_flag,idx_c,coarse)

   call mom_bc_rudman(us,vof,mom,dir,usmask,rho1,rho2,t,idx_c)

   call do_all_ghost(mom)
   call do_all_ghost(vof) 
   call do_all_ighost(v_flag)

end subroutine advect_mom_vof_coarse_stg
!--------------------------------------------------------------------------------------------------------------


!-----------------------------------------calculation of new velocity-----------------------------------------
subroutine rudman_get_velocity_coarse(mom,vof,rho1,rho2,vel,dvel,delt,dir)
   use module_rudman_init
   use module_BC
   implicit none
   include 'mpif.h'
   integer, intent(in) :: dir
   real(8), intent(in) :: delt
   real(8), dimension(idx_c(1,1):idx_c(1,2),idx_c(2,1):idx_c(2,2),idx_c(3,1):idx_c(3,2)), intent(inout) :: mom,vof
   real(8), dimension(idx_c(1,1):idx_c(1,2),idx_c(2,1):idx_c(2,2),idx_c(3,1):idx_c(3,2)), intent(inout) :: vel,dvel
   real(8), intent(in) :: rho1,rho2
   real(8) :: rhoavg
   integer :: i,j,k,ip,jp,kp  
   integer :: iend,jend,kend
   integer :: iini,jini,kini


   ip = 0 ; jp = 0 ; kp = 0

   rhoavg = 0.d0

   iini = 0 ; iend = 0 
   jini = 0 ; jend = 0 
   kini = 0 ; kend = 0 

   if (dir == 1 .and. bdry_cond(1)/=1 .and. coords(1)==nPx-1) then
      ip=1
   elseif (dir == 2 .and. bdry_cond(2)/=1 .and. coords(2)==nPy-1) then
      jp=1
   elseif (dir == 3 .and. bdry_cond(3)/=1 .and. coords(3)==nPz-1) then
      kp=1
   endif

!---------------------coarse grid indices set up for loop-----------------------------
   iini=idx_c(1,1) + Ng ; iend=idx_c(1,2) - Ng 
   jini=idx_c(2,1) + Ng ; jend=idx_c(2,2) - Ng
   kini=idx_c(3,1) + Ng ; kend=idx_c(3,2) - Ng 
!-------------------------------------------------------------------------------------

   do i=iini,iend-ip ; do j=jini,jend-ip ; do k=kini,kend-ip
   
     if ((vof(i,j,k)>0.d0) .and. (vof(i,j,k)<1.d0)) then  ! should be less than or equal to 1
        rhoavg = rho2*vof(i,j,k) + (1.d0 - vof(i,j,k))*rho1
        dvel(i,j,k) = (mom(i,j,k)/rhoavg - vel(i,j,k))/delt
     endif

   enddo ; enddo ; enddo

   call do_all_ghost(dvel)

end subroutine rudman_get_velocity_coarse

!--------------------------------------------------------------------------------------------------------


