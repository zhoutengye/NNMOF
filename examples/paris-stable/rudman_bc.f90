! Code this part later for injection conditions

!function rudman_inject(j,k,t)
!  use module_2phase
!  use module_grid
!  use module_flow
!  implicit none
!  integer :: j,k
!  real(8) :: t
!  real(8) :: uinject
!  real(8) :: ryz, low_gas_radius, NozzleThickness
!  real(8), parameter :: PI = 3.14159265359d0
!  real :: erf
!  real(8) :: temp
!  integer :: rand_seed
!  uinject=0d0
!
!  if(radius_gap_liqgas==0d0) then
!     low_gas_radius = radius_liq_inject
!  else
!     low_gas_radius = radius_gap_liqgas
!  endif
!
!  if (inject_type==1) then      ! uniform inflow
!     uinject = 1.d0
!  elseif( inject_type==2 ) then ! pulsed round jet
!     !tdelay_gas_inject = 0.01d0
!     if( (y(j) - jetcenter_yc)**2.d0 + (z(k) - jetcenter_zc)**2.d0 .lt. radius_liq_inject**2.d0 ) then 
!        uinject=uliq_inject*(1.d0+0.05d0*SIN(10.d0*2.d0*PI*t))
!     end if ! y(j)
!  else if ( inject_type == 3 ) then ! 2d coflowing jet
!     NozzleThickness = NozzleThick2Cell*dx(is) 
!     if ( y(j) <= radius_liq_inject ) then 
!        uinject = uliq_inject & 
!                 *erf( (radius_liq_inject - y(j))/blayer_gas_inject ) &
!                 *(1.d0 + erf((time-tdelay_gas_inject*0.5d0)/(tdelay_gas_inject*0.25d0)) )*0.5d0
!     else if ( y(j) > radius_liq_inject+NozzleThickness .and. y(j) <= radius_gas_inject ) then
!        uinject = ugas_inject & 
!                 *erf( (y(j) -   radius_liq_inject - NozzleThickness)/blayer_gas_inject ) & 
!                 *erf( (radius_gas_inject - y(j))/blayer_gas_inject ) & 
!                 *(1.d0 + erf((time-tdelay_gas_inject*0.5d0)/(tdelay_gas_inject*0.25d0)) )*0.5d0
!        call random_number(temp)
!        uinject = uinject*(1.d0+uinjectPertAmp*(temp-0.5d0))
!     else 
!        uinject = 0.d0 
!     end if  !
!  else if ( inject_type == 4 ) then ! 3d coaxial jet
!     !tdelay_gas_inject = 0.d-2
!     ryz = sqrt( (yh(j) - jetcenter_yc)**2.d0 + (zh(k) - jetcenter_zc)**2.d0 )
!     if ( ryz < radius_liq_inject ) then 
!        uinject = uliq_inject & 
!                 *erf( (radius_liq_inject - ryz)/blayer_gas_inject ) &
!                 *(1.d0 + erf((time-tdelay_gas_inject*0.5d0)/(tdelay_gas_inject*0.25d0)) )*0.5d0
!     else if ( ryz > low_gas_radius .and. ryz < radius_gas_inject ) then
!        uinject = ugas_inject & 
!                 *erf( (ryz - low_gas_radius)/blayer_gas_inject ) & 
!                 *erf( (radius_gas_inject - ryz)/blayer_gas_inject ) & 
!                 *(1.d0 + erf((time-tdelay_gas_inject*0.5d0)/(tdelay_gas_inject*0.25d0)) )*0.5d0
!     else 
!        uinject = 0.d0 
!     end if  !
!  else if ( inject_type == 5 ) then ! 2d coflowing with liquid on top of gas
!     if ( y(j) <= radius_gas_inject ) then 
!        uinject = ugas_inject*y(j)/radius_gas_inject
!     else if ( y(j) > radius_gas_inject .and. y(j) <= radius_liq_inject ) then
!        uinject = uliq_inject
!     else 
!        uinject = 0.d0 
!     end if  !
!  end if ! y(j), z(k)
!end function rudman_uinject


!------------------------------------GET GLAGS AND CLIP VOF-------------------------------------------------
subroutine rudman_vof_flags_and_clip(vof,vflag,n)
   implicit none
   include 'mpif.h'
   integer , dimension(3,2) , intent(in) :: n
   real(8) , dimension(n(1,1):n(1,2),n(2,1):n(2,2),n(3,1):n(3,2)), intent(inout) :: vof
   integer , dimension(n(1,1):n(1,2),n(2,1):n(2,2),n(3,1):n(3,2)), intent(inout) :: vflag
   integer :: i,j,k
   integer :: iini,iend,jini,jend,kini,kend
   real(8) , parameter :: CLIP = 1d-8


   ! setting up limmits for do loops according to fine or coarse   
   iini=n(1,1) ; iend=n(1,2)
   jini=n(2,1) ; jend=n(2,2) 
   kini=n(3,1) ; kend=n(3,2)  
  
   ! flags set to 0 for empty, 1 for full, 2 for mixed
   do i=iini,iend ; do j=jini,jend ; do k=kini,kend

      if (vof(i,j,k) .le. CLIP ) then
         vflag(i,j,k) = 0
         vof(i,j,k)=0.d0
      elseif (vof(i,j,k) .ge. (1.d0 - CLIP)) then
         vflag(i,j,k) = 1
         vof(i,j,k) = 1.d0
      else
         vflag(i,j,k) = 2
      endif

   enddo ; enddo ; enddo 

end subroutine rudman_vof_flags_and_clip 


!------------------------------------VOF BOUNDARY CONDITIONS--------------------------------------------------
subroutine vof_bc_rudman(cv,fl,nn,grid) 
    use module_flow , only: itimestep
    use module_vof
    use module_grid
    use module_BC
    use module_2phase
    use module_solid
    implicit none
    include 'mpif.h'
    integer, dimension(3,2), intent(in) :: nn 
    integer, intent(in) :: grid
    real(8), dimension(nn(1,1):nn(1,2),nn(2,1):nn(2,2),nn(3,1):nn(3,2)), intent(inout) :: cv  ! cvof
    integer, dimension(nn(1,1):nn(1,2),nn(2,1):nn(2,2),nn(3,1):nn(3,2)), intent(inout) :: fl  ! vof_flag
    integer :: fb(6),d,l,m,n,c(3),try(2:3),sign,orient,flag,flhere
    real(8) :: cvhere
    real(8) :: t1,t2
    logical :: switch = .true.

    do orient=1,6
       cond = vofbdry_cond(orient)
       if(cond=='wet') then
          fb(orient)=1
       else if(cond=='dry') then
          fb(orient)=0
       else if(cond=='periodic') then
          fb(orient) = 3  !  will do nothing
       else if(cond=='outflow') then
          fb(orient) = 4  ! will copy inflow
       else if(cond=='jet') then
          if(orient /= 1) call pariserror("jet only at x-")
          fb(orient) = 2
       else if(cond=='90deg') then 
          fb(orient) = 5
       else
          call pariserror("this vofbc not implemented")
       endif
    enddo
! 
    do orient=1,6
       ! orientation order:    
       ! x- y- z- x+ y+ z+
       d = orient
       sign=-1
       if(orient>3) then
          d = orient-3
          sign=1
       endif
       flag=fb(orient)
! sort directions so that try(1) = d and try(2),try(3) are any other two directions. 
       m=1
       n=2  
       do while (m.le.3)
          if(m.ne.d) then
             try(n) = m
             n=n+1
          endif
          m=m+1
       enddo
! end sort

       if(coords(d)==proclimit(d,sign)) then
          do l=rudman_coordstart(try(2),grid)-2,rudman_coordend(try(2),grid)+2
             do m=rudman_coordstart(try(3),grid)-2,rudman_coordend(try(3),grid)+2
                c(try(2)) = l; c(try(3)) = m
                c(d)=rudman_coordlimit(d,sign,grid) + sign
                if(flag<2) then
                   cv(c(1),c(2),c(3))=dble(flag)
                   fl(c(1),c(2),c(3))=flag
                   c(d) = c(d) + sign
                   cv(c(1),c(2),c(3))=dble(flag)
                   fl(c(1),c(2),c(3))=flag
                elseif(flag==2) then  ! jet
                   cv(c(1),c(2),c(3))=dble(inject(c(2),c(3)))
                   fl(c(1),c(2),c(3))=inject(c(2),c(3))
                elseif(flag==4) then
                   c(d)=rudman_coordlimit(d,sign,grid)
                   cvhere=cv(c(1),c(2),c(3))
                   flhere=fl(c(1),c(2),c(3))
                   c(d) = c(d) + sign
                   cv(c(1),c(2),c(3))=cvhere
                   fl(c(1),c(2),c(3))=flhere
                   c(d) = c(d) + sign
                   cv(c(1),c(2),c(3))=cvhere
                   fl(c(1),c(2),c(3))=flhere
                elseif(flag==5) then !90deg 
                   if (d==1 .and. sign==-1) then 
                      cv(c(1),c(2),c(3))=cv(c(1)+1,c(2),c(3))
                      fl(c(1),c(2),c(3))=fl(c(1)+1,c(2),c(3))
                   else if (d==1 .and. sign==1) then
                      cv(c(1),c(2),c(3))=cv(c(1)-1,c(2),c(3))
                      fl(c(1),c(2),c(3))=fl(c(1)-1,c(2),c(3))
                   else if (d==2 .and. sign==-1) then 
                      cv(c(1),c(2),c(3))=cv(c(1),c(2)+1,c(3))
                      fl(c(1),c(2),c(3))=fl(c(1),c(2)+1,c(3))
                   else if (d==2 .and. sign==1) then
                      cv(c(1),c(2),c(3))=cv(c(1),c(2)-1,c(3))
                      fl(c(1),c(2),c(3))=fl(c(1),c(2)-1,c(3))
                   else if (d==3 .and. sign==-1) then 
                      cv(c(1),c(2),c(3))=cv(c(1),c(2),c(3)+1)
                      fl(c(1),c(2),c(3))=fl(c(1),c(2),c(3)+1)
                   else if (d==3 .and. sign==1) then
                      cv(c(1),c(2),c(3))=cv(c(1),c(2),c(3)-1)
                      fl(c(1),c(2),c(3))=fl(c(1),c(2),c(3)-1)
                   end if !d
                endif
             enddo
          enddo
       endif
    enddo

        contains

             function rudman_coordstart(d,grid)
                 use module_rudman_init , only : idx_f, idx_c
                 implicit none
                 integer, intent(in) :: d, grid
                 integer :: rudman_coordstart
             
                 if      (d==1) then
                     if  (grid == 1) then
                        rudman_coordstart = idx_c(1,1) + 2
                     elseif(grid==2) then 
                        rudman_coordstart = idx_f(1,1) + 2
                     endif
                 else if (d==2) then
                     if  (grid == 1) then
                        rudman_coordstart = idx_c(2,1) + 2
                     elseif(grid==2) then 
                        rudman_coordstart = idx_f(2,1) + 2
                     endif
                 else if (d==3) then
                     if  (grid == 1) then
                        rudman_coordstart = idx_c(3,1) + 2
                     elseif(grid==2) then 
                        rudman_coordstart = idx_f(3,1) + 2
                     endif
                  else
                     rudman_coordstart = -1
                     call pariserror("rudman_coordstart: wrong sign")
                 endif
             
             end function rudman_coordstart        

             function rudman_coordend(d,grid)
               use module_rudman_init , only : idx_f, idx_c
               implicit none
               integer, intent(in) :: d, grid
               integer :: rudman_coordend

               if      (d==1) then
                   if  (grid == 1) then
                      rudman_coordend = idx_c(1,2) - 2
                   elseif(grid==2) then 
                      rudman_coordend = idx_f(1,2) - 2
                   endif
               else if (d==2) then
                   if  (grid == 1) then
                      rudman_coordend = idx_c(2,2) - 2
                   elseif(grid==2) then 
                      rudman_coordend = idx_f(2,2) - 2
                   endif
               else if (d==3) then
                   if  (grid == 1) then
                      rudman_coordend = idx_c(3,2) - 2
                   elseif(grid==2) then 
                      rudman_coordend = idx_f(3,2) - 2
                   endif
               else
                   rudman_coordend = -1
                   call pariserror("rudman_coordend: wrong sign")
               endif

             end function rudman_coordend

             function rudman_coordlimit(d,sign,grid)
               integer, intent(in) :: d,sign,grid
               integer :: rudman_coordlimit
               if(sign==1) then
                   rudman_coordlimit = rudman_coordend(d,grid)
               else if(sign==-1) then
                   rudman_coordlimit = rudman_coordstart(d,grid)   
               else 
                   rudman_coordlimit = -1
                   call pariserror("rudman_coordlimit: wrong sign")
               endif
             end function rudman_coordlimit

             function inject(j,k)
                use module_grid
                implicit none
                integer :: j,k
                integer :: inject
                inject=0
                if ( inject_type == 2 .or. inject_type == 4) then 
                   if ((y(j) - jetcenter_yc)**2 + (z(k) - jetcenter_zc)**2.lt.radius_liq_inject**2) inject=1
                else if ( inject_type == 3 .or. inject_type == 6 ) then
                   if ((y(j) - jetcenter_yc) <= radius_liq_inject ) inject = 1 
                else if ( inject_type == 5 ) then
                   if ((y(j) - jetcenter_yc) > radius_gas_inject .and. &  
                       (y(j) - jetcenter_yc) < radius_liq_inject ) inject = 1 
                end if ! inject_type
              end function inject

              function xcoord(d,i)
                implicit none
                real(8) :: xcoord
                integer, intent(in) :: d,i
                if(d==1) then
                   xcoord = x(i)
                else if(d==2) then
                   xcoord = y(i)
                else
                   xcoord = z(i)
                endif
              end function xcoord

end subroutine vof_bc_rudman


!----------------------------------VELOCITY BOUNDARY CONDITIONS---------------------------------------------------
subroutine vel_bc_rudman(vel,t,dt,AfterProjection,n)
    use module_rudman_init , only : dxh_f , dyh_f , dzh_f
    use module_bc
    use module_grid
    use module_2phase
    use module_freesurface
    implicit none
    include 'mpif.h'
    integer, dimension(3,2), intent(in) :: n 
    real(8), dimension(n(1,1):n(1,2),n(2,1):n(2,2),n(3,1):n(3,2),6), intent(inout) :: vel
    real(8), dimension(n(1,1):n(1,2),n(2,1):n(2,2),n(3,1):n(3,2)) :: u,v,w
    real(8), dimension(n(1,1):n(1,2),n(2,1):n(2,2),n(3,1):n(3,2)) :: umask,vmask,wmask
    integer, intent (in) :: AfterProjection
    real(8) :: t,dt,fluxin,tfluxin,uaverage
    real(8) :: fluxout(6),tfluxout(6),tfluxout_all,fluxratio,uinj
    integer :: i,j,k,ierr,seed
    integer :: iini,iend,jini,jend,kini,kend
    integer :: isn,ien,jsn,jen,ksn,ken

 ! as of now dxh_f,dyh_f,dzh_f are only for subgrid in this subroutine
 ! this can be easily converted to a general routine for all grids
     
 ! INITIALIZATION OF INDICES FOR COARSE GRID BC's
    iini = n(1,1)  ; iend = n(1,2) 
    jini = n(2,1)  ; jend = n(2,2)
    kini = n(3,1)  ; kend = n(3,2) 

    isn = n(1,1)+2 ; ien = n(1,2)-2
    jsn = n(2,1)+2 ; jen = n(2,2)-2
    ksn = n(3,1)+2 ; ken = n(3,2)-2

    u = 0.d0 ; v = 0.d0 ; w = 0.d0 
    umask = 0.d0 ; vmask = 0.d0 ; wmask = 0.d0

    ! definition of velocities according to arguments received
    u = vel(:,:,:,1) ; v = vel(:,:,:,2) ; w = vel(:,:,:,3)
    umask = vel(:,:,:,4) ; vmask = vel(:,:,:,5) ; wmask = vel(:,:,:,6) 

    ! solid obstacles  NO SUPPORT FOR SUBGRID VERSION YET
 !   u = u*umask
 !   v = v*vmask
 !   w = w*wmask

    ! OLD_BDRY using ifndef directives have been removed 


    ! --------------------------------------------------------------------------------------------
    ! Inflow BC
    ! --------------------------------------------------------------------------------------------
    
    ! inflow boundary condition x- with injection
    fluxin=0
    if(bdry_cond(1)==3 .and. coords(1)==0    ) then
       seed = 1317*(INT(t/1.23d-10)+1)
       call random_seed(seed)
       do j=jini,jend
          do k=kini,kend
             uinj = 1.d0 !uinject(j,k,t)
             u(isn-1,j,k)=WallVel(1,1)*uinj*umask(isn-1,j,k)
             u(isn-2,j,k)=WallVel(1,1)*uinj*umask(isn-2,j,k)
             v(isn-1,j,k)=0d0
             w(isn-1,j,k)=0d0
          enddo
       enddo
       do j=jsn,jen
          do k=ksn,ken
             fluxin = fluxin + u(isn-1,j,k)
          enddo
       enddo
    endif
    call MPI_ALLREDUCE(fluxin, tfluxin, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_Comm_Cart, ierr)
    uaverage=tfluxin/(ny*nz)
    
    ! inflow boundary condition y-
    if(bdry_cond(2)==3 .and. coords(2)==0    ) then
       do i=iini,iend
          do k=kini,kend
             v(i,jsn-1,k)=WallVel(3,2)
             v(i,jsn-2,k)=WallVel(3,2)
             u(i,jsn-1,k)=0d0
             w(i,jsn-1,k)=0d0
          enddo
       enddo
    endif
    ! inflow on z-
    if(bdry_cond(3)==3 .and. coords(3)==0   ) then
       do i=iini,iend
          do j=jini,jend
             w(i,j,ksn-1)= WallVel(5,3)
             w(i,j,ksn-2)= WallVel(5,3)
             u(i,j,ksn-1)=0d0
             v(i,j,ksn-1)=0d0
          enddo
       enddo
    endif
    ! inflow on x+
    if(bdry_cond(4)==3 .and. coords(1)==nPx-1   ) then
       do j=jini,jend
          do k=kini,kend
             u(ien,j,k)= WallVel(2,1)
             u(ien+1,j,k)= WallVel(2,1)
             v(ien+1,j,k)=0d0
             w(ien+1,j,k)=0d0
          enddo
       enddo
    endif
    ! inflow on y+
    if(bdry_cond(5)==3 .and. coords(2)==nPy-1   ) then
       do i=iini,iend
          do k=kini,kend
             v(i,jen,k)= WallVel(4,2)
             v(i,jen+1,k)= WallVel(4,2)
             u(i,jen+1,k)=0d0
             w(i,jen+1,k)=0d0
          enddo
       enddo
    endif
    ! inflow on z+
    if(bdry_cond(6)==3 .and. coords(3)==nPz-1   ) then
       do i=iini,iend
          do j=jini,jend
             w(i,j,ken)= WallVel(6,3)
             w(i,j,ken+1)= WallVel(6,3)
             u(i,j,ken+1)=0d0
             v(i,j,ken+1)=0d0
          enddo
       enddo
    endif

    ! --------------------------------------------------------------------------------------------
    ! Wall BC with velocity specified
    ! --------------------------------------------------------------------------------------------
    ! wall boundary condition x-
    if(bdry_cond(1)==0 .and. coords(1)==0    ) then
        u(isn-1,:,:)=0d0
        u(isn-2,:,:)=-u(isn,:,:)
        v(isn-1,:,:)=2*WallVel(1,2)-v(isn,:,:)
        w(isn-1,:,:)=2*WallVel(1,3)-w(isn,:,:)
    endif
    
    ! wall boundary condition x+
    if(bdry_cond(4)==0 .and. coords(1)==nPx-1) then  
        u(ien  ,:,:)=0d0
        u(ien+1,:,:)=-u(ien-1,:,:)                  ! second order extrapolation
        v(ien+1,:,:)=2*WallVel(2,2)-v(ien,:,:)      ! second order extrapolation
        w(ien+1,:,:)=2*WallVel(2,3)-w(ien,:,:)      ! second order extrapolation
    endif

    ! wall boundary condition y-
    if(bdry_cond(2)==0 .and. coords(2)==0    ) then
        v(:,jsn-1,:)=0d0
        v(:,jsn-2,:)=-v(:,jsn,:)
        u(:,jsn-1,:)=2*WallVel(3,1)-u(:,jsn,:)
        w(:,jsn-1,:)=2*WallVel(3,3)-w(:,jsn,:)
    endif
    ! wall boundary condition y+
    if(bdry_cond(5)==0 .and. coords(2)==nPy-1) then
        v(:,jen  ,:)=0d0
        v(:,jen+1,:)=-v(:,jen-1,:)
        u(:,jen+1,:)=2*WallVel(4,1)-u(:,jen,:)
        w(:,jen+1,:)=2*WallVel(4,3)-w(:,jen,:)
    endif
    ! wall boundary condition z-
    if(bdry_cond(3)==0 .and. coords(3)==0    ) then
        w(:,:,ksn-1)=0d0
        w(:,:,ksn-2)=-w(:,:,ksn)
        u(:,:,ksn-1)=2*WallVel(5,1)-u(:,:,ksn)
        v(:,:,ksn-1)=2*WallVel(5,2)-v(:,:,ksn)
    endif
    ! wall boundary condition z+
    if(bdry_cond(6)==0 .and. coords(3)==nPz-1) then
        w(:,:,ken  )=0d0
        w(:,:,ken+1)=-w(:,:,ken-1)
        u(:,:,ken+1)=2*WallVel(6,1)-u(:,:,ken)
        v(:,:,ken+1)=2*WallVel(6,2)-v(:,:,ken)
    endif
    
    ! --------------------------------------------------------------------------------------------
    ! Wall BC with shear specified (Shear is by default zero then is eqv to slip-wall BC)
    ! --------------------------------------------------------------------------------------------
    ! wall boundary condition: shear
    if(bdry_cond(1)==2 .and. coords(1)==0    ) then
        u(isn-1,:,:)=0d0
        u(isn-2,:,:)=-u(isn,:,:)
        v(isn-1,:,:) = -dxh_f(isn-1)*WallShear(1,2)+v(isn,:,:)
        w(isn-1,:,:) = -dxh_f(isn-1)*WallShear(1,3)+w(isn,:,:)
    endif
    if(bdry_cond(4)==2 .and. coords(1)==nPx-1) then
        u(ien  ,:,:)=0d0
        u(ien+1,:,:)=-u(ien-1,:,:)
        v(ien+1,:,:) = dxh_f(ien)*WallShear(2,2)+v(ien,:,:)
        w(ien+1,:,:) = dxh_f(ien)*WallShear(2,3)+w(ien,:,:)
    endif
    if(bdry_cond(2)==2 .and. coords(2)==0    ) then
        v(:,jsn-1,:)=0d0
        v(:,jsn-2,:)=-v(:,jsn,:)
        u(:,jsn-1,:) = -dyh_f(jsn-1)*WallShear(3,1)+u(:,jsn,:)
        w(:,jsn-1,:) = -dyh_f(jsn-1)*WallShear(3,3)+w(:,jsn,:)
    endif
    if(bdry_cond(5)==2 .and. coords(2)==nPy-1) then
        v(:,jen  ,:)=0d0
        v(:,jen+1,:)=-v(:,jen-1,:)
        u(:,jen+1,:) = dyh_f(jen)*WallShear(4,1)+u(:,jen,:)
        w(:,jen+1,:) = dyh_f(jen)*WallShear(4,3)+w(:,jen,:)
    endif
    if(bdry_cond(3)==2 .and. coords(3)==0    ) then
        w(:,:,ksn-1)=0d0
        w(:,:,ksn-2)=-w(:,:,ksn)
        u(:,:,ksn-1) = -dzh_f(ksn-1)*WallShear(5,1)+u(:,:,ksn)
        v(:,:,ksn-1) = -dzh_f(ksn-1)*WallShear(5,2)+v(:,:,ksn)
    endif
    if(bdry_cond(6)==2 .and. coords(3)==nPz-1) then
        w(:,:,ken  )=0d0
        w(:,:,ken+1)=-w(:,:,ken-1)
        u(:,:,ken+1) = dzh_f(ken)*WallShear(6,1)+u(:,:,ken)
        v(:,:,ken+1) = dzh_f(ken)*WallShear(6,2)+v(:,:,ken)
    endif

    ! --------------------------------------------------------------------------------------------
    ! Outflow BC with pressure specified (by default zero)
    ! Set zero normal velocity gradient 
    ! --------------------------------------------------------------------------------------------
    if (bdry_cond(1)==5 .and. coords(1)==0)then
       u(isn-1,:,:)=u(isn,:,:)
       v(isn-1,:,:)=v(isn,:,:)
       w(isn-1,:,:)=w(isn,:,:)
    endif
    
    if (bdry_cond(4)==5 .and. coords(1)==nPx-1)then
       u(ien,:,:)=u(ien-1,:,:)
       v(ien+1,:,:)=v(ien,:,:)
       w(ien+1,:,:)=w(ien,:,:)
    endif
    
    if (bdry_cond(2)==5 .and. coords(2)==0)then
       v(:,jsn-1,:)=v(:,jsn,:)
       u(:,jsn-1,:)=u(:,jsn,:)
       w(:,jsn-1,:)=w(:,jsn,:)
    endif
    
    if (bdry_cond(5)==5 .and. coords(2)==nPy-1)then
       v(:,jen,:)=v(:,jen-1,:)
       u(:,jen+1,:)=u(:,jen,:)
       w(:,jen+1,:)=w(:,jen,:)
    endif
    
    if (bdry_cond(3)==5 .and. coords(3)==0)then
       w(:,:,ksn-1)=w(:,:,ksn)
       u(:,:,ksn-1)=u(:,:,ksn)
       v(:,:,ksn-1)=v(:,:,ksn)
    endif
    
    if (bdry_cond(6)==5 .and. coords(3)==nPz-1)then
       w(:,:,ken)=w(:,:,ken-1)
       u(:,:,ken+1)=u(:,:,ken)
       v(:,:,ken+1)=v(:,:,ken)    
    endif 
    
    ! --------------------------------------------------------------------------------------------
    ! Outflow BC with velocity specified
    ! Pressure gradient is set to zero in Poisson_BC
    ! same velocity as opposing inflow. ! @generalize this !!
    ! --------------------------------------------------------------------------------------------
    if(bdry_cond(4)==4 .and. coords(1)==nPx-1) then
        if ( OutVelSpecified ) then 
           do j=jsn,jen; do k=ksn,ken
              u(ien,j,k)=uout(j)
           end do; end do
        else 
           u(ien ,:,:)=uaverage
        end if ! OutVelSpecified
        v(ien+1,:,:)=v(ien,:,:)
        w(ien+1,:,:)=w(ien,:,:)
    endif

    ! --------------------------------------------------------------------------------------------
    ! Convective BC  (du/dt+Un*du/dx =0) 
    ! --------------------------------------------------------------------------------------------
    !flux = 0.d0 
    !if(bdry_cond(4)==6 .and. coords(1)==nPx-1) then
    !   do j=js,je
    !     do k=ks,ke
    !        flux = flux + u(ie-1,j,k)
    !     enddo
    !   enddo
    !end if ! bdry_cond(4)==6

    !call MPI_ALLREDUCE(flux, tflux, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_Comm_Cart, ierr)
    !uaverage=tflux/(ny*nz)

    !if(bdry_cond(4)==6 .and. coords(1)==nPx-1) then
    !   u(ie  ,:,:)=u(ie  ,:,:)-dt*uaverage*(u(ie-1,:,:)-u(ie-2,:,:))/dx(ie-1)
    !   v(ie+1,:,:)=v(ie+1,:,:)-dt*uaverage*(v(ie  ,:,:)-v(ie-1,:,:))/dx(ie-1)
    !   w(ie+1,:,:)=w(ie+1,:,:)-dt*uaverage*(w(ie  ,:,:)-w(ie-1,:,:))/dx(ie-1)
    !end if ! bdry_cond(4)==6
     
    ! --------------------------------------------------------------------------------------------
    ! New pressure BC  (du/dn=0 & p=p0) 
    ! --------------------------------------------------------------------------------------------
    ! Note: for pressure BC, vel-BC apply after projection 
    fluxout(:) = 0.d0 
    if      (bdry_cond(1)==6 .and. coords(1)==0     & 
            .and. AfterProjection == 1 .and. .not.LateralBdry(1)) then
       do j=jsn,jen
         do k=ksn,ken
            fluxout(1) = fluxout(1) + u(isn ,j,k)
         enddo
       enddo
    else if (bdry_cond(4)==6 .and. coords(1)==nPx-1 & 
            .and. AfterProjection == 1 .and. .not.LateralBdry(4)) then
       do j=jsn,jen
         do k=ksn,ken
            fluxout(4) = fluxout(4) + u(ien-1,j,k)
         enddo
       enddo
    else if (bdry_cond(2)==6 .and. coords(2)==0     & 
            .and. AfterProjection == 1 .and. .not.LateralBdry(2)) then
       do i=isn,ien
         do k=ksn,ken
            fluxout(2) = fluxout(2) + v(i,jsn,k)
         enddo
       enddo
    else if (bdry_cond(5)==6 .and. coords(2)==nPy-1 & 
            .and. AfterProjection == 1.and. .not.LateralBdry(5)) then
       do i=isn,ien
         do k=ksn,ken
            fluxout(5) = fluxout(5) + v(i,jen-1,k)
         enddo
       enddo
    else if (bdry_cond(3)==6 .and. coords(3)==0     & 
            .and. AfterProjection == 1 .and. .not.LateralBdry(3)) then
       do i=isn,ien
         do j=jsn,jen
            fluxout(3) = fluxout(3) + w(i,j,ksn)
         enddo
       enddo
    else if (bdry_cond(6)==6 .and. coords(3)==nPz-1 & 
            .and. AfterProjection == 1 .and. .not.LateralBdry(6)) then
       do i=isn,ien
         do j=jsn,jen
            fluxout(6) = fluxout(6) + w(i,j,ken-1)
         enddo
       enddo
    endif
    call MPI_ALLREDUCE(fluxout, tfluxout, 6, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_Comm_Cart, ierr)
    tfluxout_all = sum(tfluxout)
    fluxratio = min(ABS(tfluxin/(tfluxout_all+1.d-16)),MaxFluxRatioPresBC)  ! fluxratio is capped 
    if      (bdry_cond(1)==6 .and. coords(1)==0     & 
      .and. AfterProjection == 1) then
       u(isn-1,:,:)=u(isn,:,:)*fluxratio
       v(isn-1,:,:)=v(isn,:,:)
       w(isn-1,:,:)=w(isn,:,:)
    else if (bdry_cond(4)==6 .and. coords(1)==nPx-1 .and. AfterProjection == 1) then
       u(ien,:,:) = u(ien-1,:,:)*fluxratio 
       v(ien+1,:,:)=v(ien,:,:)
       w(ien+1,:,:)=w(ien,:,:)
    else if (bdry_cond(2)==6 .and. coords(2)==0     .and. AfterProjection == 1) then
       v(:,jsn-1,:)=v(:,jsn,:)*fluxratio
       u(:,jsn-1,:)=u(:,jsn,:)
       w(:,jsn-1,:)=w(:,jsn,:)
    else if (bdry_cond(5)==6 .and. coords(2)==nPy-1 .and. AfterProjection == 1) then
       v(:,jen,:)=v(:,jen-1,:)*fluxratio
       u(:,jen+1,:)=u(:,jen,:)
       w(:,jen+1,:)=w(:,jen,:)
    else if (bdry_cond(3)==6 .and. coords(3)==0     .and. AfterProjection == 1) then
       w(:,:,ksn-1)=w(:,:,ksn)*fluxratio
       u(:,:,ksn-1)=u(:,:,ksn)
       v(:,:,ksn-1)=v(:,:,ksn)
    else if (bdry_cond(6)==6 .and. coords(3)==nPz-1 .and. AfterProjection == 1) then
       w(:,:,ken)=w(:,:,ken-1)*fluxratio
       u(:,:,ken+1)=u(:,:,ken)
       v(:,:,ken+1)=v(:,:,ken)    
    end if ! bdry_cond()

    ! Reassignment of values after bc's are applied
    vel(:,:,:,1) = u ; vel(:,:,:,2) = v ; vel(:,:,:,3) = w

end subroutine vel_bc_rudman


!-----------------------------------MOMENTUM BOUNDARY CONDITIONS-------------------------------------------------
subroutine mom_bc_rudman(u,c,mom,d,mask,rho1,rho2,t,n)
    use module_rudman_init , only : dxh_f , dyh_f , dzh_f
    use module_bc
    use module_grid
    use module_2phase
    implicit none
    include 'mpif.h'
    integer, dimension(3,2), intent(in) :: n
    real(8), dimension(n(1,1):n(1,2),n(2,1):n(2,2),n(3,1):n(3,2)), intent(inout) :: mom, u
    real(8), dimension(n(1,1):n(1,2),n(2,1):n(2,2),n(3,1):n(3,2)), intent(in) :: c
    real(8), dimension(n(1,1):n(1,2),n(2,1):n(2,2),n(3,1):n(3,2)), intent(in) :: mask
    real(8), intent(in) :: rho1,rho2
    integer, intent(in) :: d
    real(8) :: t,fluxin,tfluxin,uaverage, uinj
    integer :: i,j,k,ierr,seed
    integer :: iini,iend,jini,jend,kini,kend
    integer :: isn,ien,jsn,jen,ksn,ken


 ! dxh_f, dyh_f dzh_f are only for subgrid at this moment in this routine
 ! can be easily converted to a generalized routine for both grids 
 ! INITIALIZATION OF INDICES FOR COARSE GRID BC's
    iini = n(1,1)  ; iend = n(1,2) 
    jini = n(2,1)  ; jend = n(2,2)
    kini = n(3,1)  ; kend = n(3,2) 

    isn = n(1,1)+2 ; ien = n(1,2)-2
    jsn = n(2,1)+2 ; jen = n(2,2)-2
    ksn = n(3,1)+2 ; ken = n(3,2)-2

 ! all indices for the loops depend whether we are on subgrid or coarse grid

    ! solid obstacles : NO SUPPORT FOR SUBGRID YET 
  !  u = u*mask
  !  mom = mom*mask

    ! OLD_BDRY using ifndef directives have been removed 

    ! inflow boundary condition y-
    fluxin=0
    ! inflow boundary condition x- with injection
    if(bdry_cond(1)==3 .and. coords(1)==0    ) then
        seed = 1317*(INT(t/1.23d-10)+1)
        call random_seed(seed)
        if (d.eq.1) then
        do j=jini,jend
          do k=kini,kend
             uinj = 1.d0  !uinject(j,k,t)
             mom(isn-1,j,k)=WallVel(1,1)*uinj*(rho2*c(isn-1,j,k) + rho1*(1.d0 - c(isn-1,j,k)))
             mom(isn-2,j,k)=WallVel(1,1)*uinj*(rho2*c(isn-1,j,k) + rho1*(1.d0 - c(isn-1,j,k)))
          enddo
        enddo
        do j=jini+2,jend-2
          do k=kini+2,kend-2
             fluxin = fluxin + u(isn-1,j,k)
          enddo
        enddo
        else
            mom(isn-1,:,:) = 0.d0
        endif
    endif
    call MPI_ALLREDUCE(fluxin, tfluxin, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_Comm_Cart, ierr)
    uaverage=tfluxin/(ny*nz)

    if(bdry_cond(2)==3 .and. coords(2)==0    ) then
       if (d.eq.2) then
        do i=iini,iend
          do k=kini,kend
             mom(i,jsn-1,k)=WallVel(3,2)*(rho2*c(i,jsn-1,k) + rho1*(1.d0 - c(i,jsn-1,k)))
             mom(i,jsn-2,k)=WallVel(3,2)*(rho2*c(i,jsn-2,k) + rho1*(1.d0 - c(i,jsn-2,k)))
          enddo
        enddo
       else
           mom(:,jsn-1,:) = 0.d0
       endif
    endif
    ! inflow on z-
    if(bdry_cond(3)==3 .and. coords(3)==0   ) then
       if (d.eq.3) then
         do i=iini,iend
          do j=jini,jend
             mom(i,j,ksn-1)= WallVel(5,3)*(rho2*c(i,j,ksn-1) + rho1*(1.d0 - c(i,j,ksn-1)))
             mom(i,j,ksn-2)= WallVel(5,3)*(rho2*c(i,j,ksn-2) + rho1*(1.d0 - c(i,j,ksn-2)))
          enddo
         enddo
     else
         mom(:,:,ksn-1) = 0.d0
     endif
    endif
    ! inflow on x+
    if(bdry_cond(4)==3 .and. coords(1)==nPx-1   ) then
        if (d.eq.1) then
         do j=jini,jend
          do k=kini,kend
             mom(ien,j,k)  = WallVel(2,1)*(rho2*c(ien,j,k) + rho1*(1.d0 - c(ien,j,k)))
             mom(ien+1,j,k)= WallVel(2,1)*(rho2*c(ien+1,j,k) + rho1*(1.d0 - c(ien+1,j,k)))
          enddo
         enddo
       else
           mom(ien,:,:) = 0.d0
       endif
    endif

     ! inflow on y+
    if(bdry_cond(5)==3 .and. coords(2)==nPy-1   ) then
        if (d.eq.2) then
           do i=iini,iend
              do k=kini,kend
                mom(i,jen,k)  = WallVel(4,2)*(rho2*c(i,jen,k) + rho1*(1.d0 - c(i,jen,k)))
                mom(i,jen+1,k)= WallVel(4,2)*(rho2*c(i,jen+1,k) + rho1*(1.d0 - c(i,jen+1,k)))
              enddo
            enddo
        else
            mom(:,jen,:) = 0.d0
        endif
    endif
    ! inflow on z+
    if(bdry_cond(6)==3 .and. coords(3)==nPz-1   ) then
       if (d.eq.3) then
        do i=iini,iend
          do j=jini,jend
             mom(i,j,ken)  = WallVel(6,3)*(rho2*c(i,j,ken) + rho1*(1.d0 - c(i,j,ken)))
             mom(i,j,ken+1)= WallVel(6,3)*(rho2*c(i,j,ken+1) + rho1*(1.d0 - c(i,j,ken+1)))
          enddo
        enddo
       else
           mom(:,:,ken) = 0.d0
       endif
    endif    

    ! wall boundary condition
    if(bdry_cond(1)==0 .and. coords(1)==0    ) then
        if (d.eq.1) then
            mom(isn-1,:,:)=0d0
            mom(isn-2,:,:)=-mom(isn,:,:) ! to be corrected, it is not exact
        else !y,z
            mom(isn-1,:,:)=(2*WallVel(1,2)-u(isn,:,:))*(rho2*c(isn,:,:) + rho1*(1.d0 - c(isn,:,:))) !CHECK!!
        endif
    endif
    

    if(bdry_cond(4)==0 .and. coords(1)==nPx-1) then 
        if (d.eq.1) then
            mom(ien  ,:,:)=0d0
            mom(ien+1,:,:)=-mom(ien-1,:,:)
        else
            mom(ien+1,:,:)=(2*WallVel(2,2)-u(ien,:,:))**(rho2*c(ien,:,:) + rho1*(1.d0 - c(ien,:,:)))
        endif
    endif
    

    if(bdry_cond(2)==0 .and. coords(2)==0    ) then
        if (d.eq.2) then
            mom(:,jsn-1,:)=0d0
            mom(:,jsn-2,:)=-mom(:,jsn,:)
        else
            mom(:,jsn-1,:)=(2*WallVel(3,1)-u(:,jsn,:))*(rho2*c(:,jsn,:) + rho1*(1.d0 - c(:,jsn,:)))
        endif
    endif
    if(bdry_cond(5)==0 .and. coords(2)==nPy-1) then
        if (d.eq.2) then
            mom(:,jen  ,:)=0d0
            mom(:,jen+1,:)=-mom(:,jen-1,:)
        else
            mom(:,jen+1,:)=(2*WallVel(4,1)-u(:,jen,:))*(rho2*c(:,jen,:) + rho1*(1.d0 - c(:,jen,:)))
        endif
    endif
    if(bdry_cond(3)==0 .and. coords(3)==0    ) then
        if (d.eq.3) then
            mom(:,:,ksn-1)=0d0
            mom(:,:,ksn-2)=-mom(:,:,ksn)
        else
            mom(:,:,ksn-1)=(2*WallVel(5,1)-u(:,:,ksn))*(rho2*c(:,:,ksn) + rho1*(1.d0 - c(:,:,ksn)))
        endif
    endif
    if(bdry_cond(6)==0 .and. coords(3)==nPz-1) then
        if (d.eq.3) then
            mom(:,:,ken  )=0d0
            mom(:,:,ken+1)=-mom(:,:,ken-1)
        else
            mom(:,:,ken+1)=(2*WallVel(6,1)-u(:,:,ken))*(rho2*c(:,:,ken) + rho1*(1.d0 - c(:,:,ken)))
        endif
    endif
    ! wall boundary condition: shear
    if(bdry_cond(1)==2 .and. coords(1)==0    ) then
        if (d.eq.1) then
            mom(isn-1,:,:) = 0.d0    
            mom(isn-2,:,:) = -mom(isn,:,:)    
        elseif (d.eq.2) then
            mom(isn-1,:,:) = (-dxh(isn-1)*WallShear(1,2)+u(isn,:,:)) * &
                            (rho2*c(isn,:,:) + rho1*(1.d0 - c(isn,:,:)))
        elseif (d.eq.3) then
            mom(isn-1,:,:) = (-dxh(isn-1)*WallShear(1,3)+u(isn,:,:)) * &
                            (rho2*c(isn,:,:) + rho1*(1.d0 - c(isn,:,:)))
        endif
    endif
    if(bdry_cond(4)==2 .and. coords(1)==nPx-1) then
        if (d.eq.1) then
            mom(ien,:,:) = 0.d0
            mom(ien+1,:,:) = -mom(ien-1,:,:)    
        elseif (d.eq.2) then
            mom(ien+1,:,:) = (dxh(ien)*WallShear(2,2)+u(ien,:,:)) * &
                            (rho2*c(ien,:,:) + rho1*(1.d0 - c(ien,:,:)))
        elseif (d.eq.3) then
            mom(ien+1,:,:) = (dxh(ien)*WallShear(2,3)+u(ien,:,:)) * &
                            (rho2*c(ien,:,:) + rho1*(1.d0 - c(ien,:,:)))
        endif
    endif
    if(bdry_cond(2)==2 .and. coords(2)==0    ) then
        if (d.eq.1) then
            mom(:,jsn-1,:) = (-dyh(jsn-1)*WallShear(3,1)+u(:,jsn,:)) * &
                            (rho2*c(:,jsn,:) + rho1*(1.d0 - c(:,jsn,:)))
        elseif (d.eq.2) then
             mom(:,jsn-1,:) = 0.d0
             mom(:,jen-2,:) = -mom(:,jsn,:)   
        elseif (d.eq.3) then
            mom(:,jsn-1,:) = (-dyh(jsn-1)*WallShear(3,3)+u(:,jsn,:)) * &
                            (rho2*c(:,jsn,:) + rho1*(1.d0 - c(:,jsn,:)))
        endif
    endif
    if(bdry_cond(5)==2 .and. coords(2)==nPy-1) then
        if (d.eq.1) then
            mom(:,jen+1,:) = (dyh(jen)*WallShear(4,1)+u(:,jen,:)) * &
                            (rho2*c(:,jen,:) + rho1*(1.d0 - c(:,jen,:)))
        elseif (d.eq.2) then
            mom(:,jen,:) = 0.d0
            mom(:,jen+1,:) = -mom(:,jen-1,:)   
        elseif (d.eq.3) then
            mom(:,jen+1,:) = (dyh(jen)*WallShear(4,3)+u(:,jen,:)) * &
                            (rho2*c(:,jen,:) + rho1*(1.d0 - c(:,jen,:)))
        endif
    endif
    if(bdry_cond(3)==2 .and. coords(3)==0    ) then
        if (d.eq.1) then
            mom(:,:,ksn-1) = (-dzh(ksn-1)*WallShear(5,1)+u(:,:,ksn)) * &
                            (rho2*c(:,:,ksn) + rho1*(1.d0 - c(:,:,ksn)))
        elseif (d.eq.2) then
            mom(:,:,ksn-1) = (-dzh(ksn-1)*WallShear(5,2)+u(:,:,ksn)) * &
                            (rho2*c(:,:,ksn) + rho1*(1.d0 - c(:,:,ksn)))
        elseif (d.eq.3) then
            mom(:,:,ksn-1) = 0.d0
            mom(:,:,ksn-2) = -mom(:,:,ksn)    
        endif
    endif
    if(bdry_cond(6)==2 .and. coords(3)==nPz-1) then
        if (d.eq.1) then
            mom(:,:,ken+1) = (dzh(ken)*WallShear(6,1)+u(:,:,ken)) * &
                            (rho2*c(:,:,ken) + rho1*(1.d0 - c(:,:,ken)))
        elseif (d.eq.2) then
            mom(:,:,ken+1) = (dzh(ken)*WallShear(6,2)+u(:,:,ken)) * &
                            (rho2*c(:,:,ken) + rho1*(1.d0 - c(:,:,ken)))
        elseif (d.eq.2) then
            mom(:,:,ken) = 0.d0
            mom(:,:,ken+1) = -mom(:,:,ken-1)   
        endif
    endif
    
    ! outflow boundary condition
    if(bdry_cond(4)==4 .and. coords(1)==nPx-1) then  
        if (d.eq.1) then
            mom(ie  ,:,:)= mom(ie-1,:,:)
            mom(ie+1,:,:)=-mom(ie,:,:)
        else
            mom(ien+1,:,:)=mom(ien,:,:)
        endif
    endif

    !Set zero normal velocity gradient for pressure boundary condition
    if (bdry_cond(1)==5 .and. coords(1)==0) then
        if (d.eq.1) then
           mom(isn-1,:,:) = mom(isn,:,:)
       else
           mom(isn-1,:,:)=mom(isn,:,:)
       endif
    endif
    
    if (bdry_cond(4)==5 .and. coords(1)==nPx-1) then
       if (d.eq.1) then
           mom(ien,:,:)  = mom(ien-1,:,:)
       else
           mom(ien+1,:,:)=mom(ien,:,:)
       endif
    endif

    if (bdry_cond(2)==5 .and. coords(2)==0) then
        if (d.eq.2) then
           mom(:,jsn-1,:)= mom(:,jsn,:)
           mom(:,jsn-2,:)= mom(:,jsn,:)
        else
           mom(:,jsn-1,:)=mom(:,jsn,:)
        endif
    endif
    
    if (bdry_cond(5)==5 .and. coords(2)==nPy-1) then
       if (d.eq.2) then
           mom(:,jen,:)  = mom(:,jen-1,:)
           mom(:,jen+1,:)= mom(:,jen-1,:)
       else
           mom(:,jen+1,:)=mom(:,jen,:)
       endif
    endif
    
    if (bdry_cond(3)==5 .and. coords(3)==0) then
        if (d.eq.3) then
           mom(:,:,ksn-1)= mom(:,:,ksn)
           mom(:,:,ksn-2)= mom(:,:,ksn)
        else
           mom(:,:,ksn-1)=u(:,:,ksn)
        endif
    endif
    
    if (bdry_cond(6)==5 .and. coords(3)==nPz-1) then
        if (d.eq.3) then
            mom(:,:,ken)  = mom(:,:,ken-1)
            mom(:,:,ken+1)= mom(:,:,ken-1)
        else
            mom(:,:,ken+1)=mom(:,:,ken)
        endif
    endif



end subroutine mom_bc_rudman







